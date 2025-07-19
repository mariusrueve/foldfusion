import json
import logging
from pathlib import Path

import tomllib

from .logger import setup_logging

logger = logging.getLogger(__name__)


class Config:
    """
    Configuration manager that reads and provides access to application settings from a
    TOML file.

    This class parses a TOML configuration file and organizes settings into separate
    sections for different components of the application.
    """

    def __init__(self, config_path: Path):
        self.config_path = config_path
        self.dict = self.read_config()

        setup_logging(
            self.log_level,
            self.log_file,
        )
        logger.info(f"Loading configuration from {config_path}")
        logger.debug("Configuration loaded successfully")

    def read_config(self) -> dict:
        if not self.config_path.exists():
            logger.error(f"Config file not found at {self.config_path}")
            raise FileNotFoundError(f"Config file not found at {self.config_path}")
        try:
            logger.debug(f"Reading config file {self.config_path}")
            with open(
                self.config_path, "rb"
            ) as f:  # TOML must be opened in binary mode
                parsed_config = tomllib.load(f)

            logger.debug("Config file parsed successfully")
            return parsed_config
        except Exception as e:
            logger.exception(f"Failed to parse config file: {e}")
            raise ValueError(f"Failed to parse config file: {e}")

    @property
    def uniprot_ids(self) -> list[str]:
        """The UniProt ID from general settings."""
        uniprot_ids = self.dict["uniprot_ids"]
        if uniprot_ids is None:
            logger.error("uniprot_id is not configured")
            raise ValueError("uniprot_id is not configured")
        return uniprot_ids

    @property
    def output_dir(self) -> Path:
        """The base output directory as a Path object from general settings."""
        value = self.dict["output_dir"]
        if value is None:
            raise ValueError("output_dir is not configured")
        path_obj = Path(value)
        if not path_obj.exists():
            logger.warning(f"Output directory does not exist: {path_obj}")
            logger.info(f"Creating output directory: {path_obj}")
            path_obj.mkdir(parents=True, exist_ok=True)
        else:
            logger.debug(f"Using existing output directory: {path_obj}")
        return path_obj

    @property
    def log_level(self) -> str:
        """The log level from general settings."""
        value = self.dict["log_level"]
        if value is None:
            raise ValueError("log_level is not configured")

        valid_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        if value.upper() not in valid_levels:
            logger.warning(f"Invalid log level: {value}, defaulting to INFO")
            return "INFO"

        return value.upper()

    @property
    def log_file(self) -> Path:
        """The log file name from general settings."""
        value = self.dict["log_file"]
        if value is None:
            raise ValueError("log_file is not configured")
        # Check if log_file is an absolute path
        if Path(value).is_absolute():
            log_path = Path(value)
        else:
            # If relative path, join with output directory
            log_path = self.output_dir / value

        # Ensure the parent directory exists
        log_path.parent.mkdir(parents=True, exist_ok=True)

        logger.debug(f"Log file will be written to: {log_path}")
        return log_path

    def _get_executable_path(self, key: str) -> Path:
        """Helper to get and validate an executable path."""
        value = self.dict.get(key)
        if value is None:
            raise ValueError(f"{key} is not configured")
        path_obj = Path(value)
        if not path_obj.exists():
            logger.error(f"Executable not found at {path_obj}")
            raise FileNotFoundError(f"Executable not found at {path_obj}")
        if not path_obj.is_file():
            logger.error(f"Path is not a file: {path_obj}")
            raise ValueError(f"Path is not a file: {path_obj}")
        logger.debug(f"Using executable: {path_obj}")
        return path_obj

    # Tool executable properties
    @property
    def dogsite3_executable(self) -> Path:
        """The dogsite3 executable path."""
        return self._get_executable_path("dogsite3_executable")

    @property
    def siena_executable(self) -> Path:
        """
        Returns the path to the siena executable.
        """
        return self._get_executable_path("siena_executable")

    @property
    def siena_max_alignments(self) -> int:
        """
        Returns the maximum number of alignments for Siena.
        """
        value = self.dict.get("siena_max_alignments", 10)
        if not isinstance(value, int):
            logger.error("siena_max_alignments must be an integer.")
            raise TypeError("siena_max_alignments must be an integer.")
        return value

    @property
    def siena_db_executable(self) -> Path:
        """The siena database executable path."""
        return self._get_executable_path("siena_db_executable")

    @property
    def siena_db_database_path(self) -> Path:
        """The existing siena database path."""
        value = self.dict["siena_db_database_path"]
        if value is None:
            raise ValueError("siena_db_database_path is not configured")
        path_obj = Path(value)
        if not path_obj.exists() and not self.pdb_directory.exists():
            error = (
                f"siena_db_database_path ({path_obj}) and "
                f"pdb_directory ({self.pdb_directory}) do not exist!"
                "At least pdb_directory has to be valid"
            )
            logger.error(error)
            raise ValueError(error)
        if path_obj.exists():
            logger.debug(f"Using existing siena database: {path_obj}")
        else:
            logger.warning(
                f"Siena database not found at {path_obj}, "
                "will try to create from PDB directory"
            )
        return path_obj

    @property
    def pdb_directory(self) -> Path:
        """The PDB directory path."""
        value = self.dict["pdb_directory"]
        if value is None:
            raise ValueError("pdb_directory is not configured")
        path_obj = Path(value)
        if not path_obj.exists():
            logger.error(f"pdb_directory not found at {path_obj}")
            raise FileNotFoundError(f"pdb_directory not found at {path_obj}")
        if not path_obj.is_dir():
            logger.error(f"pdb_directory is not a directory: {path_obj}")
            raise ValueError(f"pdb_directory is not a directory: {path_obj}")
        logger.debug(f"Using PDB directory: {path_obj}")
        return path_obj

    @property
    def pdb_format(self) -> int:
        """The PDB file format (0=.ent.gz, 1=pdb)."""
        value = self.dict["pdb_format"]
        if value is None:
            raise ValueError("pdb_format is not configured")
        if not isinstance(value, int):
            logger.error(
                f"pdb_format must be an integer, got {type(value).__name__}: {value}"
            )
            raise ValueError(
                f"pdb_format must be an integer, got {type(value).__name__}: {value}"
            )
        if value not in [0, 1]:
            logger.error(f"pdb_format must be 0 (.ent.gz) or 1 (pdb), got: {value}")
            raise ValueError(f"pdb_format must be 0 (.ent.gz) or 1 (pdb), got: {value}")
        logger.debug(f"Using PDB format: {value}")
        return value

    @property
    def ligand_extractor_executable(self) -> Path:
        """
        Returns the path to the ligand_extractor executable.
        """
        return self._get_executable_path("ligand_extractor_executable")

    @property
    def jamda_scorer_executable(self) -> Path:
        """The jamda scorer executable path."""
        return self._get_executable_path("jamda_scorer_executable")

    def __repr__(self):
        """Return a string representation of the Config object."""
        return f"Config(config_path={self.config_path})"

    def __str__(self):
        """Return a readable string representation of the config."""
        return json.dumps(self.dict, indent=4, default=str)


if __name__ == "__main__":
    c = Config(Path("config.toml"))
    print(c.ligand_extractor_executable)
