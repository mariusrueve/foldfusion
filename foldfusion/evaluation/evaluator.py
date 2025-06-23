import json
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from scipy.spatial import distance_matrix

from .utils import get_coordinates, parse_pdb, parse_sdf
from .visualizer import (
    visualize_local_rmsd_region,
    visualize_protein_with_radius,
)


class Evaluator:
    def __init__(self, output_dir: Path) -> None:
        """
        Initialize the Evaluator.

        Args:
            output_dir: Legacy parameter for backward compatibility. If provided, takes precedence over config.
            config: Configuration object containing output directory and other settings.
        """
        self.output_dir = output_dir / "Evaluation"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.local_rmsd_thresholds = {"medium": 0.92, "low": 3.10}
        self.tcs_thresholds = {"medium": 0.64, "low": 1.27}
        self.output_file_json = self.output_dir / "evaluation.json"

    def _compute_rmsd(self, coords_ref: np.ndarray, coords_target: np.ndarray) -> float:
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

    def calculate_local_rmsd(
        self,
        uniprot_id,
        optimization_stage: str,
        alphafold_protein: Path,
        experimental_protein: Path,
        ligand: Path,
        radius=6.0,
    ):
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

        # Calculate RMSD between matched coordinates
        rmsd = self._compute_rmsd(af_final, ex_final)

        # Save the local RMSD result to JSON file
        # Load existing data or create new structure
        try:
            with open(self.output_file_json, "r") as f:
                data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            data = {}

        # Ensure the uniprot_id key exists
        if uniprot_id not in data:
            data[uniprot_id] = {}

        # Ensure the local_rmsd key exists
        if optimization_stage not in data[uniprot_id]:
            data[uniprot_id][optimization_stage] = {}

        # Save the score
        data[uniprot_id][optimization_stage]["local_rmsd"] = rmsd

        # Write back to file
        with open(self.output_file_json, "w") as f:
            json.dump(data, f, indent=2)
        return rmsd

    def compute_tcs(self):
        pass


if __name__ == "__main__":
    af_model = Path("tests/data/alphafold/AF-Q9Y233-F1-model_v4_processed.pdb")
    ex_model = Path("tests/data/Siena/ensemble/5AXP_13.pdb")
    ligand = Path("tests/data/LigandExtractor/5AXP/4LK_A_1003.sdf")
    e = Evaluator(Path("/home/stud2022/mrueve/Downloads/output"))

    # Calculate local RMSD
    uniprot_id = "Q9Y233"
    optimization_stage = "test"
    radius = 6.0

    rmsd_value = e.calculate_local_rmsd(
        uniprot_id=uniprot_id,
        optimization_stage=optimization_stage,
        alphafold_protein=af_model,
        experimental_protein=ex_model,
        ligand=ligand,
        radius=radius,
    )

    print(f"Local RMSD calculated: {rmsd_value:.3f} Å")

    # Create visualizations
    print("Creating protein structure visualization...")
    fig1, html1 = visualize_protein_with_radius(
        alphafold_protein=af_model,
        experimental_protein=ex_model,
        ligand=ligand,
        radius=radius,
        output_file=e.output_dir / "protein_comparison.html",
    )

    print("Creating local RMSD region visualization...")
    fig2 = visualize_local_rmsd_region(
        alphafold_protein=af_model,
        experimental_protein=ex_model,
        ligand=ligand,
        radius=radius,
        rmsd_value=rmsd_value,
        output_file=e.output_dir / "local_rmsd_analysis.html",
    )

    print(f"Visualizations saved to: {e.output_dir}")
    print("Done!")
