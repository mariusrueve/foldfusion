"""Module for running the Siena tool for protein structure alignment."""

import json
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
            "--identity",
            "0.85",
        ]
        logger.debug(f"Siena command assembled: {' '.join(command)}")
        return command

    def get_best_alignments(self, n_alignments: int) -> list[dict]:
        """Retrieves the best n alignments from Siena's output.

        Args:
            n_alignments (int): Number of best alignments to return.

        Returns:
            list[dict]: List of dictionaries containing alignment data with keys:
                - pdb_code: PDB structure code
                - chain: PDB chain identifier
                - ensemble_path: Path to the ensemble file
                - ligand_pdb_code: Ligand PDB code
                - backbone_rmsd: Backbone RMSD value
                - all_atom_rmsd: All atom RMSD value

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

            # Clean up column names by stripping whitespace
            df.columns = df.columns.str.strip()

            # Create pdb_file_name based on PDB code and original index
            df["ensemble_file_name"] = df.apply(
                lambda row: f"{str(row['PDB code']).strip()}_{row.name + 1}.pdb", axis=1
            )

            # Sort primarily by Active site identity (higher is better),
            # then by RMSD values (lower is better)
            if "Active site identity" in df.columns:
                # Ensure numeric sorting even if CSV stores it as text
                df["Active site identity"] = pd.to_numeric(
                    df["Active site identity"], errors="coerce"
                ).fillna(0)
                df = df.sort_values(
                    by=["Active site identity", "Backbone RMSD", "All atom RMSD"],
                    ascending=[False, True, True],
                )
            else:
                # Fallback to RMSD-only sorting if identity is unavailable
                df = df.sort_values(
                    by=["Backbone RMSD", "All atom RMSD"], ascending=[True, True]
                )

            # Get top n results
            top_results = df.head(n_alignments)

            results = []
            ensemble_dir = self.output_dir / "ensemble"

            for _, row in top_results.iterrows():
                pdb_code = str(row["PDB code"]).strip()
                pdb_chains = str(row["PDB chains"]).strip()
                ligand_pdb_code = str(row["Ligand PDB code"]).strip()
                ensemble_file_name = str(row["ensemble_file_name"]).strip()

                ensemble_path = (ensemble_dir / ensemble_file_name).resolve()
                if not ensemble_path.exists():
                    logger.warning(
                        f"Ensemble file not found: {ensemble_path}, skipping."
                    )
                    continue

                alignment_data = {
                    "pdb_code": pdb_code,
                    "chain": pdb_chains,
                    "ensemble_path": ensemble_path,
                    "ligand_pdb_code": ligand_pdb_code,
                    "backbone_rmsd": float(row["Backbone RMSD"]),
                    "all_atom_rmsd": float(row["All atom RMSD"]),
                }
                results.append(alignment_data)

            logger.info(f"Found {len(results)} best alignments")
            json_file = self.output_dir / "best_alignments.json"
            # Create a serializable copy of the results for JSON export
            serializable_results = [
                {**res, "ensemble_path": str(res["ensemble_path"])} for res in results
            ]
            with open(json_file, "w") as f:
                json.dump(serializable_results, f, indent=2)
            return results

        except pd.errors.EmptyDataError as e:
            logger.error("Results file is empty or malformed")
            raise ValueError("Results file is empty or malformed") from e
        except Exception as e:
            logger.error(f"Error processing results file: {e}")
            raise
