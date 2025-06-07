"""Module for creating and managing the Siena database."""

from typing import List
from .tool import Tool
from pathlib import Path
import logging

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
        database_path: Path,
        pdb_directory: Path,
        pdb_format: int,
        output_dir: Path,
    ):
        """Initializes the Siena database creation tool.

        Args:
            executable (Path): Path to the Siena database creation executable.
            database_path (Path): Path to the database file. If None/empty, defaults to "sienaDB" in output_dir.
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
        self.output_dir = output_dir / "SienaDB"

        # Determine the actual database path based on the logic:
        # 1. If database_path is given and valid, use it
        # 2. If database_path is given but invalid/empty, create new at that path
        # 3. If database_path is None/empty, create "sienaDB" at default output dir
        if database_path is None or str(database_path).strip() == "":
            # No path given, use default
            self.database_path = self.output_dir / "sienaDB"
            logger.info(
                f"No database path provided, using default: {self.database_path}"
            )
        else:
            # Path was provided
            self.database_path = Path(database_path)
            if self._is_siena_db_valid(self.database_path):
                logger.info(f"Using existing valid database at: {self.database_path}")
            else:
                logger.info(
                    f"Database path provided but invalid/empty, will create new database at: {self.database_path}"
                )

        # Initialize the base Tool class with command and output directory
        super().__init__(self.get_command(), self.output_dir)
        logger.debug(
            f"SienaDB initialized with PDB dir: {pdb_directory}, format: {pdb_format}, database: {self.database_path}"
        )

    def run(self) -> Path:
        """
        Run the SienaDB creation tool.

        If a valid database already exists at the specified path, skip creation.
        Otherwise, create a new database.

        Returns:
            Path: The path to the output directory containing the database.
        """
        # Check if we should use existing database or create new one
        if self._is_siena_db_valid(self.database_path):
            logger.info(f"Using existing valid SienaDB at: {self.database_path}")
            # Ensure output directory exists even when using existing DB
            self.output_dir.mkdir(parents=True, exist_ok=True)
            return self.output_dir.resolve()
        else:
            logger.info(f"Creating new SienaDB at: {self.database_path}")
            # Ensure the directory for the database file exists
            self.database_path.parent.mkdir(parents=True, exist_ok=True)
            # Call parent run method to execute the command
            return super().run()

    def _is_siena_db_valid(self, siena_db_path: Path) -> bool:
        """
        Check if the existing SienaDB is valid.

        Args:
            siena_db_path: Path to the SienaDB file

        Returns:
            True if the SienaDB is valid, False otherwise
        """
        if not siena_db_path.exists():
            logger.info(f"SienaDB does not exist at {siena_db_path}")
            return False

        if not siena_db_path.is_file():
            logger.warning(f"SienaDB path is not a file: {siena_db_path}")
            return False

        # Check if the file has some minimum size (empty files are invalid)
        if siena_db_path.stat().st_size == 0:
            logger.warning(f"SienaDB file is empty: {siena_db_path}")
            return False

        logger.info(f"Valid SienaDB found at {siena_db_path}")
        return True

    def get_command(self) -> List[str]:
        """Assembles the command and arguments for running the Siena database creation tool.

        Returns:
            List[str]: The command and arguments as a list of strings.

        Raises:
            ValueError: If required paths are not set.
        """
        if not self.executable or not self.pdb_directory:
            logger.error("Executable or PDB directory path is not set")
            raise ValueError("Executable and PDB directory paths must be set")

        # Use the stem (filename without extension) of the database path
        db_name = self.database_path.stem if self.database_path else "sienaDB"

        command = [
            str(self.executable),
            "--database",
            db_name,
            "--directory",
            str(self.pdb_directory.resolve()),
            "--format",
            str(self.pdb_format),
        ]
        logger.debug(f"SienaDB command assembled: {' '.join(command)}")
        return command
