from foldfusion.evaluation.utils import parse_pdb, parse_sdf, get_coordinates
from pathlib import Path
import numpy as np
from scipy.spatial import distance_matrix
from typing import Optional

import plotly.graph_objects as go


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

    def _calculate_local_rmsd(
        self,
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
        return rmsd

    def compute_tcs(self):
        pass


if __name__ == "__main__":
    af_model = Path("tests/data/alphafold/AF-Q9Y233-F1-model_v4_processed.pdb")
    ex_model = Path("tests/data/Siena/ensemble/5AXP_13.pdb")
    ligand = Path("tests/data/LigandExtractor/5AXP/4LK_A_1003.sdf")
    e = Evaluator(Path("/home/stud2022/mrueve/Downloads/output"))
