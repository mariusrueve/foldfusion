import json
import threading
from pathlib import Path

import numpy as np
from scipy.spatial import distance_matrix

from .utils import get_coordinates, parse_pdb, parse_sdf


class Evaluator:
    """Thread-safe evaluator for Local RMSD and TCS metrics."""

    def __init__(self, output_dir: Path) -> None:
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
                }

        # Precompute metrics for ALL ligands available per PDB code
        computed: list[dict] = []
        for pdb_code, lig_list in ligand_structures.items():
            align_info = pdb_alignment_map.get(pdb_code)
            ensemble_path = align_info.get("ensemble_path") if align_info else None
            for ligand_dict in lig_list:
                ligand_id = ligand_dict.get("ligand_id")
                p = ligand_dict.get("path")
                ligand_path: Path | None = Path(p) if p else None

                if ligand_path is not None and ensemble_path is not None:
                    local_rmsd = self._compute_local_rmsd(
                        alphafold_structure, Path(ensemble_path), ligand_path
                    )
                    tcs = self._compute_tcs(alphafold_structure, ligand_path)
                else:
                    local_rmsd = None
                    tcs = None

                computed.append(
                    {
                        "pdb_code": pdb_code,
                        "ligand_id": ligand_id,
                        "stage": stage,
                        "local_rmsd": local_rmsd,
                        "tcs": tcs,
                        "all_atom_rmsd": (
                            align_info.get("all_atom_rmsd") if align_info else None
                        ),
                        "backbone_rmsd": (
                            align_info.get("backbone_rmsd") if align_info else None
                        ),
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

                if ligand_id not in data[uniprot_id][pdb_code]:
                    data[uniprot_id][pdb_code][ligand_id] = {}

                stage_key = entry["stage"]
                if stage_key not in data[uniprot_id][pdb_code][ligand_id]:
                    data[uniprot_id][pdb_code][ligand_id][stage_key] = {}

                data[uniprot_id][pdb_code][ligand_id][stage_key]["local_rmsd"] = entry[
                    "local_rmsd"
                ]
                data[uniprot_id][pdb_code][ligand_id][stage_key]["tcs"] = entry["tcs"]

            with open(self.output_file_json, "w") as f:
                json.dump(data, f, indent=2)

            self.data = data

    def _compute_rmsd(self, coords_ref: np.ndarray, coords_target: np.ndarray) -> float:
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
    ) -> float:
        af_coords = get_coordinates(parse_pdb(alphafold_protein))
        ex_coords = get_coordinates(parse_pdb(experimental_protein))
        lig_coords = get_coordinates(parse_sdf(ligand))

        af_dist = distance_matrix(af_coords, lig_coords)
        ex_dist = distance_matrix(ex_coords, lig_coords)

        af_within = np.any(af_dist <= radius, axis=1)
        ex_within = np.any(ex_dist <= radius, axis=1)

        af_sel = af_coords[af_within]
        ex_sel = ex_coords[ex_within]

        if len(af_sel) == 0 or len(ex_sel) == 0:
            return float("inf")

        inter = distance_matrix(af_sel, ex_sel)
        if len(af_sel) <= len(ex_sel):
            closest = np.argmin(inter, axis=1)
            af_final = af_sel
            ex_final = ex_sel[closest]
        else:
            closest = np.argmin(inter, axis=0)
            af_final = af_sel[closest]
            ex_final = ex_sel

        return self._compute_rmsd(af_final, ex_final)

    def _compute_tcs(
        self,
        alphafold_protein: Path,
        ligand: Path,
        distance_threshold: float = 4.0,
        exclude_monoatomic_ions: bool = True,
    ) -> float:
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

        protein_atoms = protein_data["atoms"] + protein_data["hetero_atoms"]
        ligand_atoms = ligand_data["atoms"]

        if exclude_monoatomic_ions:
            ligand_atoms = [a for a in ligand_atoms if a["element"] not in mono_ions]

        if not ligand_atoms:
            return 0.0

        protein_coords = np.array([[a["x"], a["y"], a["z"]] for a in protein_atoms])
        ligand_coords = np.array([[a["x"], a["y"], a["z"]] for a in ligand_atoms])
        protein_elements = [a["element"] for a in protein_atoms]
        ligand_elements = [a["element"] for a in ligand_atoms]

        distances = distance_matrix(protein_coords, ligand_coords)
        close_pairs = np.where(distances <= distance_threshold)
        if len(close_pairs[0]) == 0:
            return 0.0

        overlaps_sq = []
        for pi, li in zip(close_pairs[0], close_pairs[1], strict=True):
            pr = vdw_radii.get(protein_elements[pi], 1.70)
            lr = vdw_radii.get(ligand_elements[li], 1.70)
            vdw_d = pr + lr
            d = distances[pi, li]
            overlap = vdw_d - d
            if overlap > 0:
                overlaps_sq.append(overlap**2)

        if not overlaps_sq:
            return 0.0

        return float(np.sqrt(np.mean(overlaps_sq)))
