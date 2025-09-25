"""Compute evaluation metrics for transplanted ligands."""

import json
import threading
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import distance_matrix

from .utils import get_coordinates, parse_pdb, parse_sdf


class Evaluator:
    """Thread-safe helper that aggregates RMSD- and clash-based metrics."""

    def __init__(self, output_dir: Path) -> None:
        """Initialise the output directory and supporting synchronisation."""
        self.output_dir = output_dir / "Evaluation"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.output_file_json = self.output_dir / "evaluation.json"
        self.data: dict = {}
        self._io_lock = threading.Lock()

    def evaluate(
        self,
        uniprot_id: str,
        stage: str,
        alphafold_structure: Path,
        siena_structures: list[dict],
        ligand_structures: dict[str, list[dict]],
    ) -> None:
        """Compute, merge, and persist metrics for a pipeline stage.

        Args:
            uniprot_id: Identifier of the protein currently being processed.
            stage: Pipeline stage label (e.g. ``"transplant"`` or ``"jamda"``).
            alphafold_structure: Path to the AlphaFold reference structure.
            siena_structures: Siena alignment metadata including ensemble paths.
            ligand_structures: Mapping of PDB code to ligands generated at the
                current stage.

        Side Effects:
            Updates ``evaluation.json`` atomically so downstream reporting can
            consume the aggregated metrics.
        """
        # Build a quick lookup of the best (or first) alignment per PDB to get
        # ensemble paths and RMSDs
        pdb_alignment_map: dict[str, dict] = {}
        for alignment in siena_structures:
            pdb_code = alignment.get("pdb_code")
            if not pdb_code:
                continue
            # Prefer the first occurrence (list is already sorted by quality)
            if pdb_code not in pdb_alignment_map:
                ep_val = alignment.get("ensemble_path")
                ep_path = Path(ep_val) if ep_val is not None else None
                pdb_alignment_map[pdb_code] = {
                    "ensemble_path": ep_path,
                    "all_atom_rmsd": alignment.get("all_atom_rmsd"),
                    "backbone_rmsd": alignment.get("backbone_rmsd"),
                    "active_site_identity": alignment.get("active_site_identity"),
                }

        # Precompute metrics for ALL ligands available per PDB code
        computed: list[dict] = []
        for pdb_code, lig_list in ligand_structures.items():
            align_info = pdb_alignment_map.get(pdb_code)
            ensemble_path = align_info.get("ensemble_path") if align_info else None
            active_site_identity = (
                align_info.get("active_site_identity") if align_info else None
            )
            for ligand_dict in lig_list:
                ligand_id = ligand_dict.get("ligand_id")
                p = ligand_dict.get("path")
                ligand_path: Path | None = Path(p) if p else None
                # Optional reference path to donor ligand (pre-optimization)
                # If absent, fall back to current ligand path (pre-jamda stage)
                ref_p = ligand_dict.get("ref_path") or ligand_dict.get("path")
                ligand_ref_path: Path | None = Path(ref_p) if ref_p else None

                if ligand_path is not None and ensemble_path is not None:
                    local_rmsd = self._compute_local_rmsd(
                        alphafold_structure, Path(ensemble_path), ligand_path
                    )
                    tcs_results = self._compute_tcs(alphafold_structure, ligand_path)
                    # LEV: all heavy atoms of ligand + protein atoms within 6 Å
                    lev_results = None
                    try:
                        if ligand_ref_path is not None and active_site_identity == 1.0:
                            lev_results = self._compute_lev(
                                alphafold_structure,
                                Path(ensemble_path),
                                ligand_path,
                                ligand_ref_path,
                            )
                    except Exception:
                        lev_results = None
                else:
                    local_rmsd = None
                    tcs_results = None
                    lev_results = None

                computed.append(
                    {
                        "pdb_code": pdb_code,
                        "ligand_id": ligand_id,
                        "stage": stage,
                        "local_rmsd": local_rmsd,
                        "tcs": tcs_results,
                        "lev": lev_results,
                        "all_atom_rmsd": (
                            align_info.get("all_atom_rmsd") if align_info else None
                        ),
                        "backbone_rmsd": (
                            align_info.get("backbone_rmsd") if align_info else None
                        ),
                        "active_site_identity": active_site_identity,
                    }
                )

        # Merge and write under a lock
        with self._io_lock:
            try:
                with open(self.output_file_json) as f:
                    data = json.load(f)
            except (FileNotFoundError, json.JSONDecodeError):
                data = {}

            if uniprot_id not in data:
                data[uniprot_id] = {}

            for entry in computed:
                pdb_code = entry["pdb_code"]
                ligand_id = entry["ligand_id"]
                if pdb_code not in data[uniprot_id]:
                    data[uniprot_id][pdb_code] = {}

                if entry["all_atom_rmsd"] is not None:
                    data[uniprot_id][pdb_code]["all_atom_rmsd"] = entry["all_atom_rmsd"]
                if entry["backbone_rmsd"] is not None:
                    data[uniprot_id][pdb_code]["backbone_rmsd"] = entry["backbone_rmsd"]
                if entry["active_site_identity"] is not None:
                    data[uniprot_id][pdb_code]["active_site_identity"] = entry[
                        "active_site_identity"
                    ]

                if ligand_id not in data[uniprot_id][pdb_code]:
                    data[uniprot_id][pdb_code][ligand_id] = {}

                stage_key = entry["stage"]
                if stage_key not in data[uniprot_id][pdb_code][ligand_id]:
                    data[uniprot_id][pdb_code][ligand_id][stage_key] = {}

                data[uniprot_id][pdb_code][ligand_id][stage_key]["local_rmsd"] = entry[
                    "local_rmsd"
                ]
                data[uniprot_id][pdb_code][ligand_id][stage_key]["tcs"] = entry["tcs"]
                # Store LEV when available
                data[uniprot_id][pdb_code][ligand_id][stage_key]["lev"] = entry["lev"]

            with open(self.output_file_json, "w") as f:
                json.dump(data, f, indent=2)

            self.data = data

    def _compute_rmsd(self, coords_ref: np.ndarray, coords_target: np.ndarray) -> float:
        """Return the Kabsch RMSD between two aligned coordinate sets."""
        ref_centroid = coords_ref.mean(axis=0)
        tgt_centroid = coords_target.mean(axis=0)
        ref = coords_ref - ref_centroid
        tgt = coords_target - tgt_centroid
        C = np.dot(ref.T, tgt)
        V, S, Wt = np.linalg.svd(C)
        d = np.sign(np.linalg.det(np.dot(Wt.T, V.T)))
        R = np.dot(Wt.T, np.diag([1, 1, d])).dot(V.T)
        return float(np.sqrt(np.mean(np.sum((np.dot(ref, R) - tgt) ** 2, axis=1))))

    def _compute_local_rmsd(
        self,
        alphafold_protein: Path,
        experimental_protein: Path,
        ligand: Path,
        radius: float = 6.0,
        trim_fraction: float = 0.1,
    ) -> float:
        """AlphaFill-style local r.m.s.d. of protein backbone environment.

        Definition (AlphaFill): RMSD of local structural alignment using all
        protein backbone atoms (N, CA, C, O) within 6 Å of the ligand.

        Implementation:
        - Select backbone atoms in AF and donor (ensemble) within `radius` of
          ligand heavy atoms.
        - Match sets via Hungarian assignment on distances.
        - Perform Kabsch alignment; optionally trim outliers and recompute.
        """
        af_parsed = parse_pdb(alphafold_protein)
        ex_parsed = parse_pdb(experimental_protein)
        lig_coords = get_coordinates(parse_sdf(ligand))

        # Collect AF and donor backbone coordinates; include CA as backbone.
        backbone_names = {"N", "CA", "C", "O", "OXT"}
        af_bb = np.array(
            [
                [a["x"], a["y"], a["z"]]
                for a in af_parsed["atoms"]
                if (
                    a.get("atom_name") in backbone_names
                    and (a.get("element") or "?") != "H"
                )
            ]
        )
        ex_bb = np.array(
            [
                [a["x"], a["y"], a["z"]]
                for a in ex_parsed["atoms"]
                if (
                    a.get("atom_name") in backbone_names
                    and (a.get("element") or "?") != "H"
                )
            ]
        )

        # Fallback: if no backbone parsed (e.g., minimal test), fall back to CA-only
        if af_bb.size == 0 or ex_bb.size == 0:
            af_bb = get_coordinates(af_parsed, atom_types=["CA"])  # may be empty
            ex_bb = get_coordinates(ex_parsed, atom_types=["CA"])  # may be empty

        if af_bb.size == 0 or ex_bb.size == 0:
            return float("inf")

        # Select neighbors around ligand
        af_dist = distance_matrix(af_bb, lig_coords)
        ex_dist = distance_matrix(ex_bb, lig_coords)
        af_within = np.any(af_dist <= radius, axis=1)
        ex_within = np.any(ex_dist <= radius, axis=1)
        af_sel = af_bb[af_within]
        ex_sel = ex_bb[ex_within]

        if len(af_sel) == 0 or len(ex_sel) == 0:
            return float("inf")

        # One-to-one assignment using Hungarian algorithm on distances
        D = distance_matrix(af_sel, ex_sel)
        if D.size == 0:
            return float("inf")
        row_ind, col_ind = linear_sum_assignment(D)
        af_m = af_sel[row_ind]
        ex_m = ex_sel[col_ind]

        if len(af_m) < 3:
            return float("inf")

        def _kabsch_rmsd(A: np.ndarray, B: np.ndarray) -> tuple[float, np.ndarray]:
            A_c = A - A.mean(axis=0)
            B_c = B - B.mean(axis=0)
            C = np.dot(A_c.T, B_c)
            V, S, Wt = np.linalg.svd(C)
            d = np.sign(np.linalg.det(np.dot(Wt.T, V.T)))
            R = np.dot(Wt.T, np.diag([1, 1, d])).dot(V.T)
            A_rot = np.dot(A_c, R)
            diff = A_rot - B_c
            rmsd = float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))
            res = np.sqrt(np.sum(diff**2, axis=1))
            return rmsd, res

        rmsd, residuals = _kabsch_rmsd(af_m, ex_m)

        # Optional trimming of top residuals
        if 0.0 < trim_fraction < 0.5 and len(residuals) >= 10:
            k = int(round(len(residuals) * (1 - trim_fraction)))
            if k >= 3:  # need at least 3 points for stable Kabsch
                keep_idx = np.argsort(residuals)[:k]
                rmsd, _ = _kabsch_rmsd(af_m[keep_idx], ex_m[keep_idx])

        return rmsd

    def _compute_lev(
        self,
        alphafold_protein: Path,
        experimental_protein: Path,
        ligand_transplanted: Path,
        ligand_donor: Path,
        radius: float = 6.0,
    ) -> dict | None:
        """Compute LEV score and its components.

        Returns a dictionary with detailed RMSD values for the ligand, the
        binding site, and the combined local environment.
        """
        # Parse structures
        af_p = parse_pdb(alphafold_protein)
        ex_p = parse_pdb(experimental_protein)
        lig_t = parse_sdf(ligand_transplanted)
        lig_d = parse_sdf(ligand_donor)

        # Ligand heavy atoms only
        lig_t_coords = np.array(
            [
                [a["x"], a["y"], a["z"]]
                for a in lig_t["atoms"]
                if (a.get("element") or "?") != "H"
            ]
        )
        lig_d_coords = np.array(
            [
                [a["x"], a["y"], a["z"]]
                for a in lig_d["atoms"]
                if (a.get("element") or "?") != "H"
            ]
        )

        if lig_t_coords.size == 0 or lig_d_coords.size == 0:
            return None

        # Protein heavy atoms (non-hydrogen)
        def _heavy(atom: dict) -> bool:
            el = atom.get("element") or "?"
            return el != "H"

        af_atoms = np.array(
            [[a["x"], a["y"], a["z"]] for a in af_p["atoms"] if _heavy(a)]
        )
        ex_atoms = np.array(
            [[a["x"], a["y"], a["z"]] for a in ex_p["atoms"] if _heavy(a)]
        )

        # Select neighbors within radius of ligand heavy atoms
        af_dmat = distance_matrix(af_atoms, lig_t_coords)
        ex_dmat = distance_matrix(ex_atoms, lig_d_coords)
        af_sel = af_atoms[np.any(af_dmat <= radius, axis=1)]
        ex_sel = ex_atoms[np.any(ex_dmat <= radius, axis=1)]

        # Match protein neighbors via Hungarian
        if len(af_sel) == 0 or len(ex_sel) == 0:
            prot_A = np.empty((0, 3))
            prot_B = np.empty((0, 3))
        else:
            Dp = distance_matrix(af_sel, ex_sel)
            rI, cI = linear_sum_assignment(Dp)
            prot_A = af_sel[rI]
            prot_B = ex_sel[cI]

        # Map ligand atoms: prefer index-wise if same count, else Hungarian
        if len(lig_t_coords) == len(lig_d_coords) and len(lig_t_coords) >= 1:
            lig_A = lig_t_coords
            lig_B = lig_d_coords
        else:
            Dl = distance_matrix(lig_t_coords, lig_d_coords)
            lr, lc = linear_sum_assignment(Dl)
            lig_A = lig_t_coords[lr]
            lig_B = lig_d_coords[lc]

        # Combined coordinates for overall LEV score
        A_combined = np.vstack([prot_A, lig_A]) if prot_A.size else lig_A
        B_combined = np.vstack([prot_B, lig_B]) if prot_B.size else lig_B

        if len(A_combined) < 3:
            return None

        # Kabsch RMSD calculation function
        def _calculate_rmsd(coords_A, coords_B):
            if len(coords_A) < 1:
                return 0.0
            A_c = coords_A - coords_A.mean(axis=0)
            B_c = coords_B - coords_B.mean(axis=0)
            C = np.dot(A_c.T, B_c)
            V, S, Wt = np.linalg.svd(C)
            d = np.sign(np.linalg.det(np.dot(Wt.T, V.T)))
            R = np.dot(Wt.T, np.diag([1, 1, d])).dot(V.T)
            A_rot = np.dot(A_c, R)
            diff = A_rot - B_c
            return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))

        # Calculate all component RMSDs
        lev_rmsd = _calculate_rmsd(A_combined, B_combined)
        bindingsite_rmsd = _calculate_rmsd(prot_A, prot_B) if prot_A.size > 0 else 0.0
        transplant_rmsd = _calculate_rmsd(lig_A, lig_B) if lig_A.size > 0 else 0.0

        return {
            "local_environment_rmsd": lev_rmsd,
            "bindingsite_rmsd": bindingsite_rmsd,
            "transplant_rmsd": transplant_rmsd,
            "bindingsite_atom_count": len(prot_A),
            "transplant_atom_count": len(lig_A),
        }

    def _compute_tcs(
        self,
        alphafold_protein: Path,
        ligand: Path,
        distance_threshold: float = 4.0,
        exclude_monoatomic_ions: bool = True,
    ) -> dict:
        """Calculate the steric clash score (TCS) for a ligand placement.

        Args:
            alphafold_protein: AlphaFold structure containing the receptor.
            ligand: Ligand SDF file to analyse.
            distance_threshold: Maximum distance in angstroms for considering
                potential clashes.
            exclude_monoatomic_ions: Drop simple monoatomic ions before scoring
                to avoid inflated clashes.

        Returns:
            A dictionary summarising the clash score, clash locations, and atom
            counts involved in the calculation.
        """
        vdw_radii = {
            "H": 1.20,
            "C": 1.70,
            "N": 1.55,
            "O": 1.52,
            "F": 1.47,
            "P": 1.80,
            "S": 1.80,
            "Cl": 1.75,
            "Br": 1.85,
            "I": 1.98,
            "Na": 2.27,
            "Mg": 1.73,
            "K": 2.75,
            "Ca": 2.31,
            "Fe": 2.00,
            "Zn": 1.39,
            "Se": 1.90,
            "B": 1.92,
            "Si": 2.10,
            "Al": 1.84,
            "Li": 1.82,
            "Be": 1.53,
            "Ne": 1.54,
            "Ar": 1.88,
            "Kr": 2.02,
            "Xe": 2.16,
            "Rn": 2.20,
            "He": 1.40,
        }
        mono_ions = {"Na", "Mg", "K", "Ca", "Zn", "Fe", "Cl", "Br", "I", "F"}

        protein_data = parse_pdb(alphafold_protein)
        ligand_data = parse_sdf(ligand)

        protein_atoms = protein_data["atoms"]
        ligand_atoms = ligand_data["atoms"]

        if exclude_monoatomic_ions:
            ligand_atoms = [a for a in ligand_atoms if a["element"] not in mono_ions]

        if not ligand_atoms:
            return {
                "score": 0.0,
                "clash_count": 0,
                "transplant_atom_count": 0,
                "poly_atom_count": 0,
                "clashes": [],
            }

        protein_coords = np.array([[a["x"], a["y"], a["z"]] for a in protein_atoms])
        ligand_coords = np.array([[a["x"], a["y"], a["z"]] for a in ligand_atoms])

        distances = distance_matrix(protein_coords, ligand_coords)
        close_pairs = np.where(distances <= distance_threshold)

        poly_atom_indices = set(close_pairs[0])
        poly_atom_count = len(poly_atom_indices)

        if len(close_pairs[0]) == 0:
            return {
                "score": 0.0,
                "clash_count": 0,
                "transplant_atom_count": len(ligand_atoms),
                "poly_atom_count": 0,
                "clashes": [],
            }

        clashes = []
        overlaps_sq = []
        for pi, li in zip(close_pairs[0], close_pairs[1], strict=True):
            protein_atom = protein_atoms[pi]
            ligand_atom = ligand_atoms[li]

            pr = vdw_radii.get(protein_atom["element"], 1.70)
            lr = vdw_radii.get(ligand_atom["element"], 1.70)
            vdw_d = pr + lr
            d = distances[pi, li]
            overlap = vdw_d - d

            if overlap > 0:
                overlaps_sq.append(overlap**2)
                clashes.append(
                    {
                        "protein_atom": {
                            "id": protein_atom["atom_name"],
                            "seq_id": protein_atom["residue_number"],
                        },
                        "ligand_atom_id": ligand_atom["atom_name"],
                        "distance": round(d, 4),
                        "vdw_overlap": round(overlap, 4),
                    }
                )

        score = float(np.sqrt(np.mean(overlaps_sq))) if overlaps_sq else 0.0

        return {
            "score": score,
            "clash_count": len(clashes),
            "transplant_atom_count": len(ligand_atoms),
            "poly_atom_count": poly_atom_count,
            "clashes": clashes,
        }
