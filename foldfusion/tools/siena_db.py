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
            # Accept both with and without .db extension; generator typically
            # creates file named exactly as provided to --database (we pass
            # stem currently). If user supplies a path ending in .db, we
            # remember both variants so we can adapt.
            configured_path = Path(database_path)
            alt_no_suffix = (
                configured_path.with_suffix("")
                if configured_path.suffix == ".db"
                else None
            )

            # Prefer an existing valid file among the candidates
            chosen_path: Path
            if self._is_siena_db_valid(configured_path):
                chosen_path = configured_path
            elif alt_no_suffix and self._is_siena_db_valid(alt_no_suffix):
                logger.debug(
                    (
                        "Configured Siena path %s not found but alternate "
                        "without suffix exists: %s"
                    ),
                    configured_path,
                    alt_no_suffix,
                )
                chosen_path = alt_no_suffix
            else:
                # None exist yet – choose the form we expect the generator to create.
                # The generator receives 'stem' (without extension) so use alt_no_suffix
                # if user provided a .db path; else use configured_path directly.
                chosen_path = alt_no_suffix or configured_path
                logger.info(
                    (
                        "SIENA database will be created at configured path: %s "
                        "(configured: %s)"
                    ),
                    chosen_path,
                    configured_path,
                )

            self.database_path = chosen_path
            self.output_dir = self.database_path.parent
            if self._is_siena_db_valid(self.database_path):
                logger.info(
                    "Using existing SIENA database at configured path: %s",
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

        # Post-run: if expected path (possibly with .db) not found but alternate exists,
        # switch to the existing one so downstream tools work.
        if not self.database_path.exists() and self.database_path.suffix == "":
            # If we used no-suffix form but a .db file appeared, adopt it.
            alt_with_db = self.database_path.with_suffix(".db")
            if alt_with_db.exists():
                logger.debug(
                    (
                        "Switching Siena DB path to file with .db suffix "
                        "found after generation: %s"
                    ),
                    alt_with_db,
                )
                self.database_path = alt_with_db
        elif not self.database_path.exists() and self.database_path.suffix == ".db":
            alt_no_suffix = self.database_path.with_suffix("")
            if alt_no_suffix.exists():
                logger.debug(
                    (
                        "Switching Siena DB path to file without .db suffix "
                        "found after generation: %s"
                    ),
                    alt_no_suffix,
                )
                self.database_path = alt_no_suffix

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

        # If user configured a path with suffix (.db), pass the full filename so
        # the generator creates exactly that file (avoids mismatch where a bare
        # name "siena_db" is produced while pipeline expects "siena_db.db").
        if self.database_path:
            if self.database_path.suffix == ".db":
                db_name = self.database_path.name  # keep extension
            else:
                db_name = self.database_path.name  # same behavior as before
        else:
            db_name = "siena_db"
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
