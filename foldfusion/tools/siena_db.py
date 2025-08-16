"""Module for creating and managing the Siena database."""

import logging
from pathlib import Path

from .tool import Tool

logger = logging.getLogger(__name__)


class SienaDB(Tool):
    """Runs the Siena database creation tool to build a database of protein structures.

    This tool processes a directory of PDB files to create a database that can be used
    by the Siena tool for protein structure alignment.

    Attributes:
        executable (Path): Path to the Siena database creation executable.
        pdb_directory (Path): Directory containing PDB files to process.
        pdb_format (int): Format of PDB files (0=.ent.gz, 1=pdb).
        output_dir (Path): Directory to save the database files.
    """

    def __init__(
        self,
        executable: Path,
        database_path: Path | None,
        pdb_directory: Path,
        pdb_format: int,
        output_dir: Path,
    ):
        """Initializes the Siena database creation tool.

        Args:
            executable (Path): Path to the Siena database creation executable.
            database_path (Path | None): Path to the database file. If None/empty,
                defaults to a file named "siena_db" inside output_dir / "SienaDB".
            pdb_directory (Path): Directory containing PDB files to process.
            pdb_format (int): Format of PDB files (0=.ent.gz, 1=pdb).
            output_dir (Path): Directory to save the database files.

        Raises:
            FileNotFoundError: If the executable or PDB directory does not exist.
            ValueError: If the PDB format is invalid or paths are invalid.
        """
        if not executable.exists():
            logger.error(f"Siena database executable not found: {executable}")
            raise FileNotFoundError(
                f"Siena database executable not found: {executable}"
            )
        if not executable.is_file():
            logger.error(f"Siena database executable is not a file: {executable}")
            raise ValueError(f"Siena database executable is not a file: {executable}")
        if not pdb_directory.exists():
            logger.error(f"PDB directory not found: {pdb_directory}")
            raise FileNotFoundError(f"PDB directory not found: {pdb_directory}")
        if not pdb_directory.is_dir():
            logger.error(f"PDB directory is not a directory: {pdb_directory}")
            raise ValueError(f"PDB directory is not a directory: {pdb_directory}")
        if pdb_format not in (0, 1):
            logger.error(f"Invalid PDB format: {pdb_format}")
            raise ValueError("PDB format must be 0 (.ent.gz) or 1 (pdb)")

        self.executable = executable
        self.pdb_directory = pdb_directory
        self.pdb_format = pdb_format

        # Base output dir for default case (when no explicit DB path given)
        default_db_dir = output_dir / "SienaDB"
        self.output_dir = default_db_dir

        # Path selection logic
        if database_path and str(database_path).strip():
            self.database_path = Path(database_path)
            # Use parent directory of provided path as working directory
            self.output_dir = self.database_path.parent
            if self._is_siena_db_valid(self.database_path):
                logger.info(
                    "Using existing SIENA database at configured path: %s",
                    self.database_path,
                )
            else:
                logger.info(
                    "SIENA database will be created at configured path: %s",
                    self.database_path,
                )
        else:
            self.database_path = default_db_dir / "siena_db"
            if self._is_siena_db_valid(self.database_path):
                logger.info(
                    "Using existing default SIENA database: %s",
                    self.database_path,
                )
            else:
                logger.info(
                    "No SIENA database path provided; will create at default: %s",
                    self.database_path,
                )

        # Initialize the base Tool class with command and output directory
        super().__init__(self.get_command(), self.output_dir)
        logger.debug(
            "SienaDB initialized with PDB dir: %s, format: %s, database: %s",
            pdb_directory,
            pdb_format,
            self.database_path,
        )

    def run(self) -> Path:
        """
        Run the SienaDB creation tool.

        If a valid database already exists at the specified path, skip creation.
        Otherwise, create a new database.

        Returns:
            Path: The path to the output directory containing the database.
        """
        # If database already exists and is valid, skip execution.
        if self._is_siena_db_valid(self.database_path):
            logger.info(
                "Skipping SienaDB generation; existing database is valid: %s",
                self.database_path,
            )
            return self.database_path.resolve()

        _ = super().run()
        return self.database_path.resolve()

    def _is_siena_db_valid(self, siena_db_path: Path) -> bool:
        """
        Check if the existing SienaDB is valid.

        Args:
            siena_db_path: Path to the SienaDB file

        Returns:
            True if the SienaDB is valid, False otherwise
        """
        try:
            if not siena_db_path.exists():
                logger.debug(f"SienaDB does not exist at {siena_db_path}")
                return False

            if not siena_db_path.is_file():
                logger.warning(f"SienaDB path is not a file: {siena_db_path}")
                return False

            # Consider a minimal size threshold (arbitrary small size)
            if siena_db_path.stat().st_size < 1024:  # <1KB likely invalid
                logger.warning(
                    "SienaDB file too small (<1KB), treating as invalid: %s",
                    siena_db_path,
                )
                return False

            logger.debug(f"Valid SienaDB found at {siena_db_path}")
            return True
        except Exception as e:
            logger.warning(f"Error while validating SienaDB at {siena_db_path}: {e}")
            return False

    def get_command(self) -> list[str]:
        """Assemble the command list for Siena database creation.

        Returns:
            List[str]: Executable and arguments.
        """
        if not self.executable or not self.pdb_directory:
            logger.error("Executable or PDB directory path is not set")
            raise ValueError("Executable and PDB directory paths must be set")

        db_name = self.database_path.stem if self.database_path else "siena_db"
        command = [
            str(self.executable),
            "--database",
            db_name,
            "--directory",
            str(self.pdb_directory.resolve()),
            "--format",
            str(self.pdb_format),
        ]
        logger.debug("SienaDB command assembled: %s", " ".join(command))
        return command
