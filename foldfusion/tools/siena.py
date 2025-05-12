from foldfusion.tools.tool import Tool
from pathlib import Path
import subprocess
import logging
import pandas as pd

logger = logging.getLogger(__name__)


class Siena(Tool):
    def __init__(self, config: dict, dogsite3_output_edf: Path):
        super().__init__(config)
        self.load_executable_config("siena")
        self.dogsite3_output = dogsite3_output_edf
        self.output_dir = self.output_dir / "siena"
        self.siena_db_path = self.initialize_siena_database()
        self.assemble_command()

    def initialize_siena_database(self):
        siena_db_config = self.config.get("siena_db", {})
        siena_db_path = Path(siena_db_config.get("siena_db", ""))
        pdb_files_path = Path(siena_db_config.get("pdb_files", ""))

        if (not siena_db_path or not siena_db_path.exists()) and (
            not pdb_files_path or not pdb_files_path.exists()
        ):
            raise ValueError(
                "Both siena_db and pdb_files paths are invalid or empty. "
                + "At least one has to be valid/defined"
            )
        if siena_db_path == "" or not siena_db_path.exists():
            logger.info("siena_db was not found, creating new one")
            executable_path = siena_db_config.get("executable")
            pdb_dir = siena_db_config.get("pdb_directory", "")
            pdb_format = siena_db_config.get("format", 1)
            command = [
                executable_path,
                "--database",
                "siena_db",
                "--directory",
                pdb_dir,
                "--format",
                str(pdb_format),
            ]

            # Execute the command
            try:
                _ = subprocess.run(command)
                siena_db_path = Path.cwd() / "siena_db"
                logger.info(f"New siena database was created at {siena_db_path}")

            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"Failed to run {self.tool_name}: {e}")
        return siena_db_path

    def assemble_command(self):
        siena_config = self.config.get("siena", {})
        executable_path = siena_config.get("executable", "")
        optional_arguments = siena_config.get("optional_arguments", "")
        self.command = [
            executable_path,
            "-e",
            self.dogsite3_output,
            "-b",
            self.siena_db_path,
            "-o",
            ".",
            optional_arguments,
        ]

    def get_best_alignments(self, n_alignments: int) -> list[list[str]]:
        csv_file = self.output_dir / "resultStatistic.csv"
        if not csv_file.exists():
            raise FileExistsError(f"{csv_file} does not exist!")
        df = pd.read_csv(str(csv_file), delimiter=";")
        # Corrected the 'by' parameter for sort_values
        df = df.sort_values(by=["Backbone RMSD", "All atom RMSD"], ascending=True)
        # Corrected column selection and conversion to list
        results = df.head(n_alignments)[["PDB code", "PDB chains"]].values.tolist()
        # Strip whitespace from both values in each result entry
        results = [
            [pdb_code.strip(), pdb_chains.strip()] for pdb_code, pdb_chains in results
        ]
        logger.info(f"Best alignments are:\n{results}")
        return results
