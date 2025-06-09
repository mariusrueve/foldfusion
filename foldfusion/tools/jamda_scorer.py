from .tool import Tool
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class JamdaScorer(Tool):
    def __init__(
        self,
        executable: Path,
        alphafold_pdb_path: Path,
        ligand_extractor_output_path: Path,
        output_dir: Path,
    ):
        self.executable = executable
        self.alphafold_pdb_path = alphafold_pdb_path
        self.ligand_extractor_outout_path = ligand_extractor_output_path
        self.output_dir = output_dir / "JamdaScorer"

    def run(self):
        # Get all subfolders in ligand_extractor_output_path
        subfolders = [
            f for f in self.ligand_extractor_outout_path.iterdir() if f.is_dir()
        ]
        logger.info(
            f"Found {len(subfolders)} subfolders in ligand extractor output directory:"
            + f" {subfolders}"
        )
        # Process each subfolder
        for subfolder in subfolders:
            # Get the last part of the subfolder path (the directory name)
            pdb_code = subfolder.name

            logger.info(f"Processing PDB: {pdb_code}")

            ligand_output_dir = self.output_dir / pdb_code
            ligand_output_dir.mkdir(parents=True, exist_ok=True)

            # Get all .sdf files in the current subfolder
            sdf_files = list(subfolder.glob("*.sdf"))
            logger.info(
                f"Found {len(sdf_files)} .sdf files in {pdb_code}: {[f.name for f in sdf_files]}"
            )

            if not sdf_files:
                logger.warning(f"No .sdf files found in {pdb_code}, skipping")
                continue

            # Process each .sdf file
            for sdf_file in sdf_files:
                logger.info(f"Processing PDB: {pdb_code}, SDF file: {sdf_file.name}")
                self.command = [
                    self.executable,
                    "-i",
                    self.alphafold_pdb_path,
                    "-m",
                    sdf_file,
                    "-o",
                    ligand_output_dir / Path(str(sdf_file.name)),
                    "--optimize",
                ]
                super().run()

        return self.output_dir
