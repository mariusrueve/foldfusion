"""Module for running the DoGSiteScorer tool for binding site prediction."""

import logging
from pathlib import Path

from .tool import Tool

logger = logging.getLogger(__name__)


class Dogsite3(Tool):
    """Runs the DoGSiteScorer (DoGSite3) tool to predict binding pockets in protein structures.

    Attributes:
        pdb_file (Path): Path to the input PDB file for pocket prediction.
        output_dir (Path): Directory to save DoGSite3 output files.
    """

    def __init__(self, executable: Path, pdb_file: Path, output_dir: Path):
        """Initializes the Dogsite3 tool.

        Args:
            config (dict): The configuration dictionary.
            pdb_file (Path): The path to the PDB file to analyze.

        Raises:
            FileNotFoundError: If the provided PDB file does not exist.
        """
        if not pdb_file.exists():
            logger.error(f"Input PDB file not found: {pdb_file}")
            raise FileNotFoundError(f"Input PDB file not found: {pdb_file}")

        self.executable = executable
        self.pdb_file = pdb_file
        self.output_dir = output_dir / "Dogsite3"
        logger.debug(f"Dogsite3 initialized with PDB file: {self.pdb_file}")
        super().__init__(self.get_command(), self.output_dir)

    def get_command(self) -> list[str]:
        """Assembles the command and arguments for running DoGSite3."""
        if not self.pdb_file:
            logger.error("PDB file path is not set. Cannot assemble DoGSite3 command.")
            raise ValueError(
                "PDB file path must be set before assembling command arguments."
            )

        command = ["--proteinFile", str(self.pdb_file.resolve())]  # Use resolved path
        command = [str(self.executable)] + command + ["--writeSiteResiduesEDF"]
        logger.debug(f"DoGSite3 command assembled: {' '.join(command)}")
        return command

    def get_best_edf(self) -> Path:
        """Retrieves the path to the best pocket's EDF file and updates its reference PDB path.

        The EDF (Pocket Description File) with the suffix "_P_1_res.edf" is assumed
        to be the one corresponding to the best-ranked pocket.
        This method also updates the REFERENCE field in the EDF file to point to the
        correct input PDB file if it was originally "<NO-FILE>".

        Returns:
            Path: The absolute path to the (potentially modified) best EDF file.

        Raises:
            FileNotFoundError: If the expected EDF file is not found in the output directory.
        """
        edf_filename = "output_P_1_res.edf"
        edf_file_path = self.output_dir / edf_filename
        logger.info(f"Looking for best EDF file: {edf_file_path}")

        if not edf_file_path.is_file():
            logger.error(
                f"EDF file not found: {edf_file_path}. "
                f"Ensure the Dogsite3 tool has been run and produced this output file."
            )
            raise FileNotFoundError(
                f"EDF file not found: {edf_file_path}. "
                f"Ensure the Dogsite3 tool has been run and produced this output file."
            )

        try:
            lines = edf_file_path.read_text().splitlines()
        except Exception as e:
            logger.error(f"Error reading EDF file {edf_file_path}: {e}")
            raise

        new_lines = []
        modified = False
        # Ensure pdb_file is an absolute path for the REFERENCE field
        reference_pdb_path_str = str(self.pdb_file.resolve())

        for line in lines:
            # Strip whitespace from the line for robust comparison
            if line.strip() == "REFERENCE <NO-FILE>":
                new_lines.append(f"REFERENCE {reference_pdb_path_str}")
                modified = True
                logger.info(
                    f"Updated REFERENCE in EDF file to: {reference_pdb_path_str}"
                )
            else:
                new_lines.append(line)

        if modified:
            try:
                # Join lines with newline character and add a final newline
                content_to_write = "\n".join(new_lines) + "\n"
                edf_file_path.write_text(content_to_write)
                logger.info(
                    f"Successfully modified and saved EDF file: {edf_file_path}"
                )
            except Exception as e:
                logger.error(f"Error writing modified EDF file {edf_file_path}: {e}")
                raise
        else:
            logger.info(
                f"No modification needed for REFERENCE in EDF file: {edf_file_path}"
            )

        resolved_edf_path = edf_file_path.resolve()
        logger.info(f"The best EDF file path: {resolved_edf_path}")
        return resolved_edf_path
