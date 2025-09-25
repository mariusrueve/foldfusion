"""Wrapper around the JAMDA scoring CLI used for ligand optimisation."""

import logging
from pathlib import Path

from .tool import Tool

logger = logging.getLogger(__name__)

LigandRecord = dict[str, str]
LigandMap = dict[str, list[LigandRecord]]


class JamdaScorer(Tool):
    """Optimise transplanted ligands with JAMDA and collect the results."""

    def __init__(
        self,
        executable: Path,
        alphafold_pdb_path: Path,
        ligand_structure: LigandMap,
        output_dir: Path,
    ):
        """Configure the JAMDA scorer wrapper.

        Args:
            executable: Path to the JAMDA CLI executable.
            alphafold_pdb_path: AlphaFold structure used as the receptor model.
            ligand_structure: Transplanted ligands to score, keyed by PDB code.
            output_dir: Pipeline output directory; a ``JamdaScorer/`` folder is
                created inside it.
        """
        self.executable = executable
        self.alphafold_pdb_path = alphafold_pdb_path
        self.ligand_structure: LigandMap = ligand_structure
        self.output_dir = output_dir / "JamdaScorer"

    def run(self) -> Path:
        """Optimise each ligand with JAMDA and expose the resulting mapping.

        Returns:
            Absolute path to the ``JamdaScorer/`` output directory. The method
            also publishes ``self.optimized_ligand_structure`` mirroring the
            input structure but pointing to the optimised ligands and keeping a
            ``ref_path`` to the donor molecule.
        """
        optimized_ligand_structure: LigandMap = {}
        for pdb_code, ligands in self.ligand_structure.items():
            ligand_output_dir = self.output_dir / pdb_code
            ligand_output_dir.mkdir(parents=True, exist_ok=True)
            optimized_ligand_structure[pdb_code] = []
            for ligand in ligands:
                sdf_file = Path(ligand["path"])
                ligand_id = ligand["ligand_id"]
                output_sdf = ligand_output_dir / f"{ligand_id}.sdf"
                self.command = [
                    str(self.executable),
                    "-i",
                    str(self.alphafold_pdb_path),
                    "-m",
                    str(sdf_file),
                    "-o",
                    str(output_sdf),
                    "--optimize",
                ]
                try:
                    super().run()
                except Exception as e:
                    logger.error(
                        "JAMDA scoring failed for %s %s: %s",
                        pdb_code,
                        ligand_id,
                        e,
                    )
                    continue
                optimized_ligand_structure[pdb_code].append(
                    {
                        "ligand_id": ligand_id,
                        "path": str(output_sdf.absolute()),
                        "sdf_file": f"{ligand_id}.sdf",
                        # Keep original donor ligand path to enable LEV computation
                        "ref_path": str(sdf_file.absolute()),
                    }
                )

        # Expose results and return output directory to match base class contract
        self.optimized_ligand_structure: LigandMap = optimized_ligand_structure
        return self.output_dir
