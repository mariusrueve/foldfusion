import json
from pathlib import Path

import numpy as np
from scipy.spatial import distance_matrix

from .utils import get_coordinates, parse_pdb, parse_sdf
from .visualizer import (
    visualize_local_rmsd_region,
    visualize_protein_with_radius,
)


class Evaluator:
    """
    Evaluator for protein structure alignment and ligand transplant quality assessment.

    This class provides metrics to evaluate the quality of protein structure alignments
    and ligand transplants, specifically for AlphaFold models enhanced with experimental
    ligand binding sites.

    Metrics:
    - Local RMSD: Root Mean Square Deviation of protein atoms near ligand binding sites.
      Lower values indicate better structural alignment (good: <0.92Å, poor: >3.10Å).
    - TCS (Transplant Clash Score): Van der Waals overlap score for ligand-protein clashes.
      Lower values indicate fewer atomic clashes and better transplant quality
      (good: <0.64Å, poor: >1.27Å).

    The evaluator processes multiple protein structures, alignments, and ligands,
    storing results in a hierarchical JSON format organized by:
    UniProt ID -> PDB Code -> Ligand ID -> Stage -> Metrics
    """

    def __init__(self, output_dir: Path) -> None:
        """
        Initialize the Evaluator.

        Args:
            output_dir: Directory where evaluation results will be stored.
                       Creates an 'Evaluation' subdirectory for output files.
        """
        self.output_dir = output_dir / "Evaluation"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Quality thresholds for Local RMSD (Angstroms)
        # Lower values indicate better structural alignment
        # good: < medium threshold, medium: medium-low range, poor: > low threshold
        self.local_rmsd_thresholds = {"medium": 0.92, "low": 3.10}

        # Quality thresholds for TCS (Angstroms)
        # Lower values indicate fewer clashes and better transplant quality
        # good: < medium threshold, medium: medium-low range, poor: > low threshold
        self.tcs_thresholds = {"medium": 0.64, "low": 1.27}

        self.output_file_json = self.output_dir / "evaluation.json"
        self.data = {}

    def evaluate(
        self,
        uniprot_id: str,
        stage: str,
        alphafold_structure,
        siena_structures,
        ligand_structures,
    ):
        """
        Evaluate protein structure alignments and ligand transplant quality.

        This method computes Local RMSD and TCS metrics for all ligand-protein
        combinations and stores the results in a hierarchical JSON structure.

        Args:
            uniprot_id: UniProt identifier for the target protein
            stage: Processing stage identifier (e.g., 'initial', 'refined')
            alphafold_structure: Path to AlphaFold protein structure file
            siena_structures: List of alignment structures from SIENA, each containing:
                - pdb_code: PDB code of the template structure
                - ensemble_path: Path to the aligned ensemble structure
                - ligand_pdb_code: Ligand identifier from the template
            ligand_structures: Dictionary mapping PDB codes to lists of ligand structures,
                each containing:
                - ligand_id: Ligand identifier
                - path: Path to ligand structure file

        The method updates the internal data structure and saves results to JSON.
        Results are organized as: uniprot_id -> pdb_code -> ligand_id -> stage -> metrics
        where metrics include 'local_rmsd' and 'tcs'.
        """
        try:
            with open(self.output_file_json, "r") as f:
                self.data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            self.data = {}

        # Ensure the uniprot_id key exists
        if uniprot_id not in self.data:
            self.data[uniprot_id] = {}

        # Process each alignment structure
        for alignment in siena_structures:
            pdb_code = alignment["pdb_code"]
            ensemble_path = alignment["ensemble_path"]
            ligand_pdb_code = alignment["ligand_pdb_code"]

            if pdb_code not in self.data[uniprot_id]:
                self.data[uniprot_id][pdb_code] = {}

            if ligand_pdb_code not in self.data[uniprot_id][pdb_code]:
                self.data[uniprot_id][pdb_code][ligand_pdb_code] = {}

            if stage not in self.data[uniprot_id][pdb_code][ligand_pdb_code]:
                self.data[uniprot_id][pdb_code][ligand_pdb_code][stage] = {}

            # Find the ligand structure matching the ligand_pdb_code
            ligand_path = None
            for pdb_code_key in ligand_structures:
                for ligand_dict in ligand_structures[pdb_code_key]:
                    if ligand_dict["ligand_id"] == ligand_pdb_code:
                        ligand_path = ligand_dict["path"]
                        break

            if ligand_path is not None:
                local_rmsd = self._compute_local_rmsd(
                    alphafold_structure,
                    ensemble_path,
                    ligand_path,
                )
                tcs = self._compute_tcs(
                    alphafold_structure,
                    ligand_path,
                )
                self.data[uniprot_id][pdb_code][ligand_pdb_code][stage][
                    "local_rmsd"
                ] = local_rmsd
                self.data[uniprot_id][pdb_code][ligand_pdb_code][stage]["tcs"] = tcs
            else:
                self.data[uniprot_id][pdb_code][ligand_pdb_code][stage][
                    "local_rmsd"
                ] = None
                self.data[uniprot_id][pdb_code][ligand_pdb_code][stage]["tcs"] = None

        with open(self.output_file_json, "w") as f:
            json.dump(self.data, f, indent=2)

    def _compute_rmsd(self, coords_ref: np.ndarray, coords_target: np.ndarray) -> float:
        """
        Compute Root Mean Square Deviation (RMSD) between two sets of coordinates.

        Uses the Kabsch algorithm to find the optimal rotation matrix that minimizes
        the RMSD between reference and target coordinate sets after centering.

        Args:
            coords_ref: Reference coordinates array of shape (n_atoms, 3)
            coords_target: Target coordinates array of shape (n_atoms, 3)

        Returns:
            float: RMSD value in Angstroms. Lower values indicate better alignment.
        """
        ref_centroid = coords_ref.mean(axis=0)
        tgt_centroid = coords_target.mean(axis=0)
        ref = coords_ref - ref_centroid
        tgt = coords_target - tgt_centroid
        C = np.dot(ref.T, tgt)
        V, S, Wt = np.linalg.svd(C)
        d = np.sign(np.linalg.det(np.dot(Wt.T, V.T)))
        R = np.dot(Wt.T, np.diag([1, 1, d])).dot(V.T)
        rmsd = np.sqrt(np.mean(np.sum((np.dot(ref, R) - tgt) ** 2, axis=1)))
        return rmsd

    def _compute_local_rmsd(
        self,
        alphafold_protein: Path,
        experimental_protein: Path,
        ligand: Path,
        radius=6.0,
    ):
        """
        Compute Local RMSD for protein atoms within a specified radius of the ligand.

        This metric evaluates the structural similarity between AlphaFold and experimental
        protein structures specifically in the ligand binding region. Only protein atoms
        within the specified radius of any ligand atom are considered.

        Args:
            alphafold_protein: Path to AlphaFold protein structure (PDB format)
            experimental_protein: Path to experimental protein structure (PDB format)
            ligand: Path to ligand structure (SDF format)
            radius: Distance cutoff in Angstroms for selecting binding site atoms (default: 6.0)

        Returns:
            float: Local RMSD value in Angstroms.
                  - Lower values indicate better binding site alignment
                  - Good: < 0.92Å, Medium: 0.92-3.10Å, Poor: > 3.10Å
                  - Returns inf if no atoms are found within radius
        """
        alphafold_coords = get_coordinates(parse_pdb(alphafold_protein))
        experimental_coords = get_coordinates(parse_pdb(experimental_protein))
        ligand_coords = get_coordinates(parse_sdf(ligand))

        # Calculate distances from each protein atom to all ligand atoms
        af_distances = distance_matrix(alphafold_coords, ligand_coords)
        ex_distances = distance_matrix(experimental_coords, ligand_coords)

        # Find atoms within radius of any ligand atom
        af_within_radius = np.any(af_distances <= radius, axis=1)
        ex_within_radius = np.any(ex_distances <= radius, axis=1)

        # Get coordinates of atoms within radius for each structure
        af_selected = alphafold_coords[af_within_radius]
        ex_selected = experimental_coords[ex_within_radius]

        if len(af_selected) == 0 or len(ex_selected) == 0:
            return float("inf")

        # Find the best matching pairs of atoms between the two selections
        # Calculate distance matrix between selected atoms from both structures
        inter_distances = distance_matrix(af_selected, ex_selected)

        # Use the smaller set as reference and find closest matches
        if len(af_selected) <= len(ex_selected):
            # For each AF atom, find the closest EX atom
            closest_indices = np.argmin(inter_distances, axis=1)
            af_final = af_selected
            ex_final = ex_selected[closest_indices]
        else:
            # For each EX atom, find the closest AF atom
            closest_indices = np.argmin(inter_distances, axis=0)
            af_final = af_selected[closest_indices]
            ex_final = ex_selected

        return self._compute_rmsd(af_final, ex_final)

    def _compute_tcs(
        self,
        alphafold_protein: Path,
        ligand: Path,
        distance_threshold=4.0,
        exclude_monoatomic_ions=True,
    ):
        """
        Compute Transplant Clash Score (TCS).

        The TCS evaluates the physical plausibility of ligand transplants by measuring
        van der Waals overlaps between ligand and protein atoms. It provides an
        orthogonal quality measure to Local RMSD by focusing on atomic clashes.

        The calculation involves:
        1. Finding all protein-ligand atom pairs within distance_threshold (4Å)
        2. Computing van der Waals overlaps for these pairs
        3. Taking the root mean square of positive overlaps (clashes)

        Args:
            alphafold_protein: Path to AlphaFold protein structure (PDB format)
            ligand: Path to ligand structure (SDF format)
            distance_threshold: Distance cutoff in Angstroms for clash detection (default: 4.0)
            exclude_monoatomic_ions: Whether to exclude monoatomic ions from calculation
                                   to avoid biasing results (default: True)

        Returns:
            float: TCS value in Angstroms.
                  - Lower values indicate fewer clashes and better transplant quality
                  - Good: < 0.64Å, Medium: 0.64-1.27Å, Poor: > 1.27Å
                  - 0.0 indicates no clashes or no close contacts

        Note:
            High TCS values may indicate:
            - Local inaccuracies in the AlphaFold model
            - Suboptimal ligand transplant positioning
            - Incompatible binding site geometry
        """
        # Van der Waals radii in Angstroms for common elements
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

        # Monoatomic ions to potentially exclude
        monoatomic_ions = {"Na", "Mg", "K", "Ca", "Zn", "Fe", "Cl", "Br", "I", "F"}

        # Parse structures
        protein_data = parse_pdb(alphafold_protein)
        ligand_data = parse_sdf(ligand)

        # Get all protein atoms (including HETATM records)
        protein_atoms = protein_data["atoms"] + protein_data["hetero_atoms"]
        ligand_atoms = ligand_data["atoms"]

        # Filter out monoatomic ions from ligand if requested
        if exclude_monoatomic_ions:
            ligand_atoms = [
                atom for atom in ligand_atoms if atom["element"] not in monoatomic_ions
            ]

        if len(ligand_atoms) == 0:
            return 0.0  # No atoms to calculate TCS for

        # Extract coordinates and elements
        protein_coords = np.array(
            [[atom["x"], atom["y"], atom["z"]] for atom in protein_atoms]
        )
        ligand_coords = np.array(
            [[atom["x"], atom["y"], atom["z"]] for atom in ligand_atoms]
        )

        protein_elements = [atom["element"] for atom in protein_atoms]
        ligand_elements = [atom["element"] for atom in ligand_atoms]

        # Calculate distance matrix between protein and ligand atoms
        distances = distance_matrix(protein_coords, ligand_coords)

        # Find atom pairs within distance threshold
        close_pairs = np.where(distances <= distance_threshold)

        if len(close_pairs[0]) == 0:
            return 0.0  # No close contacts

        overlaps_squared = []

        # Calculate van der Waals overlaps for close pairs
        for protein_idx, ligand_idx in zip(close_pairs[0], close_pairs[1]):
            protein_element = protein_elements[protein_idx]
            ligand_element = ligand_elements[ligand_idx]

            # Get van der Waals radii (use default if element not found)
            protein_radius = vdw_radii.get(
                protein_element, 1.70
            )  # Default to carbon radius
            ligand_radius = vdw_radii.get(ligand_element, 1.70)

            # Calculate expected van der Waals distance and actual distance
            vdw_distance = protein_radius + ligand_radius
            actual_distance = distances[protein_idx, ligand_idx]

            # Calculate overlap (positive values indicate clash)
            overlap = vdw_distance - actual_distance

            if overlap > 0:  # Only consider positive overlaps (clashes)
                overlaps_squared.append(overlap**2)

        if len(overlaps_squared) == 0:
            return 0.0  # No clashes detected

        # Calculate root mean square of overlaps
        tcs = np.sqrt(np.mean(overlaps_squared))

        return tcs
