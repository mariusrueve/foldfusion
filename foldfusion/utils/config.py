"""
Configuration Management Module

This module provides comprehensive configuration management for the FoldFusion pipeline.
It handles loading, validation, and access to all pipeline parameters from TOML
configuration files.

Key Features:
    - TOML configuration file parsing with error handling
    - Automatic logging system initialization
    - Property-based access to configuration values
    - Built-in validation for critical parameters
    - Path resolution and directory creation
    - Type-safe configuration access

The configuration system is designed to be robust and provide clear error messages
when configuration issues are encountered, helping users quickly identify and
resolve setup problems.
"""

import json
import logging
import tomllib
from pathlib import Path

from .logger import setup_logging

logger = logging.getLogger(__name__)


class Config:
    """
    Configuration manager for the FoldFusion pipeline.

    This class handles loading and validating configuration parameters from TOML files.
    It provides type-safe property access to all configuration values and automatically
    initializes the logging system based on configuration parameters.

    The configuration is organized into logical sections:
    - Logging parameters (log_level, log_file)
    - Input/output settings (uniprot_ids, output_dir)
    - External tool configurations (executable paths)
    - Pipeline parameters (alignment limits, file formats)

    Attributes:
        config_path: Path to the TOML configuration file
        dict: Parsed configuration dictionary from the TOML file

    Example:
        Loading and using configuration:

        >>> from pathlib import Path
        >>> config = Config(Path("config.toml"))
        >>> print(config.log_level)
        'INFO'
        >>> print(config.uniprot_ids)
        ['Q8CA95', 'Q9QYJ6']

    Raises:
        FileNotFoundError: If the configuration file doesn't exist
        ValueError: If the configuration file has syntax errors or missing values
    """

    def __init__(self, config_path: Path) -> None:
        """
        Initialize configuration from a TOML file.

        Args:
            config_path: Path to the TOML configuration file.

        Raises:
            FileNotFoundError: If the configuration file doesn't exist.
            ValueError: If the configuration file cannot be parsed or contains
                       invalid values.

        Note:
            The constructor automatically initializes the logging system based
            on the log_level and log_file parameters from the configuration.
        """
        self.config_path = config_path
        self.dict = self.read_config()

        # Initialize logging as early as possible
        try:
            setup_logging(self.log_level, self.log_file)
            logger.info(f"Configuration loaded successfully from: {config_path}")
            logger.debug("Logging system initialized from configuration")
        except Exception as e:
            # Fall back to basic logging if configuration-based setup fails
            logging.basicConfig(level=logging.INFO)
            logger.error(f"Failed to initialize logging from configuration: {e}")
            raise

    def read_config(self) -> dict:
        """
        Read and parse the TOML configuration file.

        Returns:
            Dictionary containing the parsed configuration data.

        Raises:
            FileNotFoundError: If the configuration file doesn't exist.
            ValueError: If the file cannot be parsed as valid TOML.

        Note:
            TOML files must be opened in binary mode for proper parsing.
        """
        if not self.config_path.exists():
            error_msg = f"Configuration file not found: {self.config_path}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)

        if not self.config_path.is_file():
            error_msg = f"Configuration path is not a file: {self.config_path}"
            logger.error(error_msg)
            raise ValueError(error_msg)

        try:
            logger.debug(f"Reading configuration file: {self.config_path}")
            with open(self.config_path, "rb") as f:
                parsed_config = tomllib.load(f)

            logger.debug("Configuration file parsed successfully")
            logger.debug(f"Configuration sections: {list(parsed_config.keys())}")
            return parsed_config

        except tomllib.TOMLDecodeError as e:
            error_msg = f"Invalid TOML syntax in configuration file: {e}"
            logger.error(error_msg)
            raise ValueError(error_msg) from e

        except Exception as e:
            error_msg = f"Failed to read configuration file: {e}"
            logger.error(error_msg)
            raise ValueError(error_msg) from e

    @property
    def uniprot_ids(self) -> list[str]:
        """
        List of UniProt IDs to process in the pipeline.

        This property supports two configuration modes:
        1. Direct list of UniProt IDs in the TOML file
        2. Path to a text file containing one UniProt ID per line

        Returns:
            List of UniProt identifiers as strings.

        Raises:
            ValueError: If uniprot_ids is not configured, is empty, or
                       file cannot be read.
            FileNotFoundError: If uniprot_ids_file path doesn't exist.

        Example:
            Configuration with direct list:
            >>> uniprot_ids = ["Q8CA95", "Q9QYJ6", "Q9Y233"]

            Configuration with file path:
            >>> uniprot_ids_file = "/path/to/uniprot_ids.txt"

            Usage:
            >>> config.uniprot_ids
            ['Q8CA95', 'Q9QYJ6', 'Q9Y233']
        """
        # Check if uniprot_ids is provided as a list
        uniprot_ids_list = self.dict.get("uniprot_ids")
        uniprot_ids_file = self.dict.get("uniprot_ids_file")

        # Validate that exactly one method is provided
        if uniprot_ids_list is not None and uniprot_ids_file is not None:
            error_msg = (
                "Both 'uniprot_ids' and 'uniprot_ids_file' are configured. "
                "Please use only one method."
            )
            logger.error(error_msg)
            raise ValueError(error_msg)

        if uniprot_ids_list is None and uniprot_ids_file is None:
            error_msg = (
                "Neither 'uniprot_ids' nor 'uniprot_ids_file' is configured. "
                "One must be provided."
            )
            logger.error(error_msg)
            raise ValueError(error_msg)

        # Handle direct list configuration
        if uniprot_ids_list is not None:
            if not isinstance(uniprot_ids_list, list):
                error_msg = "uniprot_ids must be a list of strings"
                logger.error(error_msg)
                raise TypeError(error_msg)

            if not uniprot_ids_list:
                error_msg = "uniprot_ids list cannot be empty"
                logger.error(error_msg)
                raise ValueError(error_msg)

            # Validate that all items are strings
            for i, uid in enumerate(uniprot_ids_list):
                if not isinstance(uid, str):
                    error_msg = (
                        f"uniprot_ids[{i}] must be a string, "
                        f"got {type(uid).__name__}: {uid}"
                    )
                    logger.error(error_msg)
                    raise TypeError(error_msg)

            logger.debug(f"Loaded {len(uniprot_ids_list)} UniProt IDs from direct list")
            return uniprot_ids_list

        # Handle file path configuration
        # At this point, uniprot_ids_file is guaranteed to be not None
        assert uniprot_ids_file is not None
        file_path = Path(uniprot_ids_file)

        if not file_path.exists():
            error_msg = f"UniProt IDs file not found: {file_path}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)

        if not file_path.is_file():
            error_msg = f"UniProt IDs path is not a file: {file_path}"
            logger.error(error_msg)
            raise ValueError(error_msg)

        try:
            logger.debug(f"Reading UniProt IDs from file: {file_path}")
            with open(file_path, encoding="utf-8") as f:
                lines = f.readlines()

            # Process lines: strip whitespace and filter out empty lines
            uniprot_ids_from_file = []
            for i, line in enumerate(lines, 1):
                stripped_line = line.strip()
                if stripped_line:  # Skip empty lines
                    # Basic validation: check for unusual lengths
                    if len(stripped_line) < 3 or len(stripped_line) > 20:
                        logger.warning(
                            f"Line {i}: UniProt ID '{stripped_line}' has unusual length"
                        )
                    uniprot_ids_from_file.append(stripped_line)

            if not uniprot_ids_from_file:
                error_msg = f"No valid UniProt IDs found in file: {file_path}"
                logger.error(error_msg)
                raise ValueError(error_msg)

            logger.debug(
                f"Loaded {len(uniprot_ids_from_file)} UniProt IDs "
                f"from file: {file_path}"
            )
            return uniprot_ids_from_file

        except Exception as e:
            error_msg = f"Failed to read UniProt IDs from file {file_path}: {e}"
            logger.error(error_msg)
            raise ValueError(error_msg) from e

    @property
    def output_dir(self) -> Path:
        """
        Base output directory for all pipeline results.

        Returns:
            Path object pointing to the output directory. Creates the directory
            if it doesn't exist.

        Raises:
            ValueError: If output_dir is not configured.

        Note:
            The directory will be created automatically if it doesn't exist,
            including any necessary parent directories.
        """
        value = self.dict["output_dir"]
        if value is None:
            raise ValueError("output_dir is not configured")

        path_obj = Path(value)
        if not path_obj.exists():
            logger.info(f"Creating output directory: {path_obj}")
            path_obj.mkdir(parents=True, exist_ok=True)
        else:
            logger.debug(f"Using existing output directory: {path_obj}")

        return path_obj

    @property
    def log_level(self) -> str:
        """
        Logging level for the pipeline.

        Returns:
            String representing the log level (DEBUG, INFO, WARNING, ERROR, CRITICAL).

        Raises:
            ValueError: If log_level is not configured.

        Note:
            Invalid log levels will be automatically corrected to INFO with a warning.
        """
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
        """
        Path to the log file for pipeline execution.

        Returns:
            Path object pointing to the log file. Creates parent directories
            if they don't exist.

        Raises:
            ValueError: If log_file is not configured.

        Note:
            Relative paths are resolved relative to the output directory.
            Absolute paths are used as-is.
        """
        value = self.dict["log_file"]
        if value is None:
            raise ValueError("log_file is not configured")

        # Handle both absolute and relative paths
        if Path(value).is_absolute():
            log_path = Path(value)
        else:
            log_path = self.output_dir / value

        # Ensure parent directory exists
        log_path.parent.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Log file configured at: {log_path}")

        return log_path

    def _get_executable_path(self, key: str) -> Path:
        """
        Validate and return an executable path from configuration.

        Args:
            key: Configuration key for the executable path.

        Returns:
            Path object pointing to the validated executable.

        Raises:
            ValueError: If the key is not configured or path is not a file.
            FileNotFoundError: If the executable file doesn't exist.

        Note:
            This helper method ensures all executables are properly validated
            before use in the pipeline.
        """
        value = self.dict.get(key)
        if value is None:
            raise ValueError(f"{key} is not configured")

        path_obj = Path(value)
        if not path_obj.exists():
            error_msg = f"Executable not found: {path_obj} (config key: {key})"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)

        if not path_obj.is_file():
            error_msg = f"Path is not a file: {path_obj} (config key: {key})"
            logger.error(error_msg)
            raise ValueError(error_msg)

        logger.debug(f"Validated executable for {key}: {path_obj}")
        return path_obj

    @property
    def dogsite3_executable(self) -> Path:
        """
        Path to the DoGSite3 executable for binding site prediction.

        Returns:
            Path object pointing to the DoGSite3 executable.

        Raises:
            ValueError: If dogsite3_executable is not configured.
            FileNotFoundError: If the executable doesn't exist.
        """
        return self._get_executable_path("dogsite3_executable")

    @property
    def siena_executable(self) -> Path:
        """
        Path to the SIENA executable for structure alignment.

        Returns:
            Path object pointing to the SIENA executable.

        Raises:
            ValueError: If siena_executable is not configured.
            FileNotFoundError: If the executable doesn't exist.
        """
        return self._get_executable_path("siena_executable")

    @property
    def siena_max_alignments(self) -> int:
        """
        Maximum number of structural alignments to retrieve from SIENA.

        Returns:
            Integer specifying the maximum number of alignments (default: 10).

        Raises:
            TypeError: If the value is not an integer.

        Note:
            This parameter controls the trade-off between processing time
            and alignment quality. Higher values provide more options but
            increase computational cost.
        """
        value = self.dict.get("siena_max_alignments", 10)
        if not isinstance(value, int):
            error_msg = "siena_max_alignments must be an integer"
            logger.error(error_msg)
            raise TypeError(error_msg)

        if value <= 0:
            error_msg = "siena_max_alignments must be positive"
            logger.error(error_msg)
            raise ValueError(error_msg)

        logger.debug(f"Using maximum {value} SIENA alignments")
        return value

    @property
    def siena_db_executable(self) -> Path:
        """
        Path to the SIENA database creation executable.

        Returns:
            Path object pointing to the SIENA database executable.

        Raises:
            ValueError: If siena_db_executable is not configured.
            FileNotFoundError: If the executable doesn't exist.
        """
        return self._get_executable_path("siena_db_executable")

    @property
    def siena_db_database_path(self) -> Path | None:
        """Optional path to the SIENA database file.

        Behaviour:
                - Path configured & exists: reuse it.
                - Path configured but missing: ensure parent dir and use it
                    (database created there later).
                - No path configured: return None; default inside output directory
                    will be used.

        Returns:
                Path | None: Existing / to-be-created database path, or None if
                not configured.

        Raises:
            ValueError: If a path is configured but neither the database nor
            the PDB directory exists (cannot create DB).
        """

        value = self.dict.get("siena_db_database_path")
        if value in (None, ""):
            logger.info(
                "No explicit SIENA database path configured; using default location."
            )
            return None

        path_obj = Path(value)

        # Ensure we can at least create it later (need valid PDB directory)
        if not path_obj.exists():
            if not self.pdb_directory.exists():  # triggers validation/logging
                error_msg = (
                    f"Cannot create SIENA database at {path_obj}: PDB directory "
                    f"missing: {self.pdb_directory}"
                )
                logger.error(error_msg)
                raise ValueError(error_msg)
            # Will be created later
            logger.info(
                (
                    "SIENA database will be created at configured path: %s "
                    "from PDB directory: %s"
                ),
                path_obj,
                self.pdb_directory,
            )
            return path_obj

        # Path exists
        if path_obj.is_file():
            logger.debug(f"Using existing SIENA database: {path_obj}")
            return path_obj

        # If path exists but is not a file (e.g. directory) treat as invalid
        error_msg = (
            f"Configured siena_db_database_path exists but is not a file: {path_obj}"
        )
        logger.error(error_msg)
        raise ValueError(error_msg)

    @property
    def pdb_directory(self) -> Path:
        """
        Directory containing PDB structure files.

        Returns:
            Path object pointing to the PDB directory.

        Raises:
            ValueError: If pdb_directory is not configured or not a directory.
            FileNotFoundError: If the directory doesn't exist.

        Note:
            This directory should contain PDB files in the format specified
            by the pdb_format configuration parameter.
        """
        value = self.dict["pdb_directory"]
        if value is None:
            raise ValueError("pdb_directory is not configured")

        path_obj = Path(value)
        if not path_obj.exists():
            error_msg = f"PDB directory not found: {path_obj}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)

        if not path_obj.is_dir():
            error_msg = f"PDB directory path is not a directory: {path_obj}"
            logger.error(error_msg)
            raise ValueError(error_msg)

        logger.debug(f"Using PDB directory: {path_obj}")
        return path_obj

    @property
    def pdb_format(self) -> int:
        """
        PDB file format specification.

        Returns:
            Integer indicating file format:
            - 0: Compressed format (.ent.gz)
            - 1: Standard PDB format (.pdb)

        Raises:
            ValueError: If pdb_format is not configured or has invalid value.
            TypeError: If the value is not an integer.

        Note:
            This parameter determines how PDB files are read from the
            configured PDB directory.
        """
        value = self.dict["pdb_format"]
        if value is None:
            raise ValueError("pdb_format is not configured")

        if not isinstance(value, int):
            error_msg = (
                f"pdb_format must be an integer, got {type(value).__name__}: {value}"
            )
            logger.error(error_msg)
            raise TypeError(error_msg)

        if value not in [0, 1]:
            error_msg = f"pdb_format must be 0 (.ent.gz) or 1 (.pdb), got: {value}"
            logger.error(error_msg)
            raise ValueError(error_msg)

        format_desc = ".ent.gz" if value == 0 else ".pdb"
        logger.debug(f"Using PDB format: {value} ({format_desc})")
        return value

    @property
    def ligand_extractor_executable(self) -> Path:
        """
        Path to the ligand extraction executable.

        Returns:
            Path object pointing to the ligand extractor executable.

        Raises:
            ValueError: If ligand_extractor_executable is not configured.
            FileNotFoundError: If the executable doesn't exist.
        """
        return self._get_executable_path("ligand_extractor_executable")

    @property
    def jamda_scorer_executable(self) -> Path:
        """
        Path to the JAMDA scoring executable for ligand optimization.

        Returns:
            Path object pointing to the JAMDA scorer executable.

        Raises:
            ValueError: If jamda_scorer_executable is not configured.
            FileNotFoundError: If the executable doesn't exist.
        """
        return self._get_executable_path("jamda_scorer_executable")

    @property
    def pipeline_concurrency(self) -> int:
        """Number of pipelines to run in parallel.

        Returns:
            Positive integer specifying the maximum number of UniProt IDs
            to process concurrently. Defaults to 1 (sequential).

        Raises:
            TypeError: If the configured value is not an integer.
            ValueError: If the value is less than 1.
        """
        value = self.dict.get("pipeline_concurrency", 1)
        if not isinstance(value, int):
            error_msg = "pipeline_concurrency must be an integer"
            logger.error(error_msg)
            raise TypeError(error_msg)
        if value < 1:
            error_msg = "pipeline_concurrency must be at least 1"
            logger.error(error_msg)
            raise ValueError(error_msg)
        logger.debug(f"Using pipeline concurrency: {value}")
        return value

    def __repr__(self) -> str:
        """Return a concise string representation of the Config object."""
        return f"Config(config_path={self.config_path})"

    def __str__(self) -> str:
        """Return a formatted JSON representation of the configuration."""
        return json.dumps(self.dict, indent=4, default=str)


if __name__ == "__main__":
    c = Config(Path("config.toml"))
    print(c.ligand_extractor_executable)
