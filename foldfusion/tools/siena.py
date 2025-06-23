"""Module for running the Siena tool for protein structure alignment."""

import logging
from pathlib import Path

import pandas as pd

from .tool import Tool

logger = logging.getLogger(__name__)


class Siena(Tool):
    """Runs the Siena tool to align protein structures and find similar binding sites.

    Attributes:
        edf (Path): Path to the input EDF file from DoGSite3.
        siena_db (Path): Path to the Siena database.
        pdb_directory (Path): Directory containing PDB files for reference.
        output_dir (Path): Directory to save Siena output files.
    """

    def __init__(
        self,
        executable: Path,
        edf: Path,
        siena_db: Path,
        pdb_directory: Path,
        output_dir: Path,
    ):
        """Initializes the Siena tool.

        Args:
            edf (Path): Path to the EDF file from DoGSite3.
            siena_db (Path): Path to the Siena database.
            pdb_directory (Path): Directory containing PDB files for reference.

        Raises:
            FileNotFoundError: If the EDF file or Siena database does not exist.
            ValueError: If any of the required paths are invalid.
        """
        if not edf.exists():
            logger.error(f"EDF file not found: {edf}")
            raise FileNotFoundError(f"EDF file not found: {edf}")
        if not siena_db.exists():
            logger.error(f"Siena database not found: {siena_db}")
            raise FileNotFoundError(f"Siena database not found: {siena_db}")
        if not pdb_directory.exists() or not pdb_directory.is_dir():
            logger.error(f"Invalid PDB directory: {pdb_directory}")
            raise ValueError(f"Invalid PDB directory: {pdb_directory}")
        self.executable = executable
        self.edf = edf
        self.siena_db = siena_db
        self.pdb_directory = pdb_directory
        self.output_dir = output_dir / "Siena"

        # Initialize the base Tool class with command and output directory
        super().__init__(self.get_command(), self.output_dir)
        logger.debug(f"Siena initialized with EDF: {edf}, DB: {siena_db}")

    def get_command(self) -> list[str]:
        """Assembles the command and arguments for running Siena.

        Returns:
            list[str]: The command and arguments as a list of strings.

        Raises:
            ValueError: If required paths are not set.
        """
        if not self.edf or not self.siena_db:
            logger.error("EDF file or Siena database path is not set")
            raise ValueError("EDF file and Siena database paths must be set")

        # Required arguments with default values based on the help output
        command = [
            str(self.executable.resolve()),
            "--edf",
            str(self.edf.resolve()),
            "--database",
            str(self.siena_db.resolve()),
            "--output",
            ".",
        ]
        logger.debug(f"Siena command assembled: {' '.join(command)}")
        return command

    def get_best_alignments(self, n_alignments: int) -> list:
        """Retrieves the best n alignments from Siena's output.

        Args:
            n_alignments (int): Number of best alignments to return.

        Returns:
            list[list[str]]: List of [PDB code, PDB chains, ensemble_path] triplets for the best alignments.

        Raises:
            FileNotFoundError: If the results CSV file does not exist.
            ValueError: If the CSV file is malformed or empty.
        """
        csv_file = self.output_dir / "resultStatistic.csv"
        if not csv_file.exists():
            logger.error(f"Results file not found: {csv_file}")
            raise FileNotFoundError(f"Results file not found: {csv_file}")

        try:
            df = pd.read_csv(str(csv_file), delimiter=";")
            if df.empty:
                logger.error("Results file is empty")
                raise ValueError("Results file is empty")

            # Sort by RMSD values (lower is better)
            df = df.sort_values(by=["Backbone RMSD", "All atom RMSD"], ascending=True)

            # Get top n results
            top_results = df.head(n_alignments)

            results = []
            ensemble_dir = self.output_dir / "ensemble"

            for index, (_, row) in enumerate(top_results.iterrows()):
                pdb_code = row["PDB code"].strip()
                pdb_chains = row["PDB chains"].strip()

                # Find corresponding ensemble file - pattern appears to be {PDB_code}_{index}.pdb
                # We use the original dataframe index + some offset for the ensemble file naming
                ensemble_files = list(ensemble_dir.glob(f"{pdb_code}_*.pdb"))

                if ensemble_files:
                    # Sort to get consistent ordering and take the first match
                    ensemble_files.sort()
                    ensemble_path = str(ensemble_files[0])
                else:
                    # If no ensemble file found, use empty string or None
                    ensemble_path = ""
                    logger.warning(f"No ensemble file found for {pdb_code}")

                results.append([pdb_code, pdb_chains, ensemble_path])

            logger.info(f"Found {len(results)} best alignments")
            return results

        except pd.errors.EmptyDataError:
            logger.error("Results file is empty or malformed")
            raise ValueError("Results file is empty or malformed")
        except Exception as e:
            logger.error(f"Error processing results file: {e}")
            raise
