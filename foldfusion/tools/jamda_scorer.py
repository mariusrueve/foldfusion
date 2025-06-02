from .tool import Tool
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class JamdaScorer(Tool):
    def __init__(
        self, config: dict, alphafold_pdb_path: Path, ligand_extractor_output_path: Path
    ):
        super().__init__(config)
        self.load_executable_config("jamda_scorer")
        self.alphafold_pdb_path = alphafold_pdb_path
        self.ligand_extractor_outout_path = ligand_extractor_output_path
        self.output_dir = self.output_dir / "jamda_scorer"

    def run_all(self):
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
                self.run()


if __name__ == "__main__":
    from foldfusion.utils.config import Config
    import logging

    # Setup basic logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler()],
    )
    logger = logging.getLogger(__name__)

    config_path = Path("config.toml")
    logger.info(f"Loading configuration from {config_path.resolve()}")
    config = Config(config_path)
    config_dict = config.dict

    alphafold_pdb_path = Path(
        "/home/stud2022/mrueve/Downloads/output/alphafold/AF-Q9Y233-F1-model_v4.pdb"
    )
    ligand_extractor_outout_path = Path(
        "/home/stud2022/mrueve/Downloads/output/ligand_extractor"
    )
    js = JamdaScorer(config_dict, alphafold_pdb_path, ligand_extractor_outout_path)
    js.run_all()
