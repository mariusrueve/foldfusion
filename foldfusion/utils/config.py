import logging
from pathlib import Path

import tomllib

from foldfusion.utils.logger import setup_logging
import json

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
        self.general_settings = self._fetch_tool_config("general_settings")
        setup_logging(
            self.log_level,  # Use property
            self.log_file,   # Use property
            self.output_dir, # Use property (returns Path)
        )
        logger.info(f"Loading configuration from {config_path}")

        self.dogsite3 = self._fetch_tool_config("dogsite3") # Corrected: was "siena"
        self.siena = self._fetch_tool_config("siena")
        self.siena_db = self._fetch_tool_config("siena_db")
        self.ligand_extractor = self._fetch_tool_config("ligand_extractor")
        self.jamda_scorer = self._fetch_tool_config("jamda_scorer")
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

    def _validate_config(self, tool_name: str) -> bool:
        match tool_name:
            case "general_settings":
                required_keys = ["uniprot_id", "output_dir", "log_level", "log_file"]
                tool_config = self.dict.get(tool_name, {})

                # Check if all required keys exist
                for key in required_keys:
                    if key not in tool_config:
                        logger.error(
                            f"Missing required setting '{key}' in {tool_name} configuration"
                        )
                        return False

                # Check if values are not empty
                for key, value in tool_config.items():
                    if key in required_keys and (value is None or value == ""):
                        logger.error(
                            f"Empty value for required setting '{key}' in {tool_name} configuration"
                        )
                        return False

                # Validate output_dir exists or can be created
                output_dir = Path(tool_config["output_dir"])
                if not output_dir.exists() and not output_dir.parent.exists():
                    logger.error(f"Invalid output_dir path: {output_dir}")
                    return False

                # Validate log_level is a valid level
                valid_log_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
                if tool_config["log_level"] not in valid_log_levels:
                    logger.error(f"Invalid log_level: {tool_config['log_level']}")
                    return False

                return True
            case "dogsite3" | "siena" | "ligand_extractor" | "jamda_scorer":
                required_keys = ["executable"]
                tool_config = self.dict.get(tool_name, {})

                for key in required_keys:
                    if key not in tool_config:
                        logger.error(
                            f"Missing required setting '{key}' in {tool_name} configuration"
                        )
                        return False
                    if tool_config[key] is None or tool_config[key] == "":
                        logger.error(
                            f"Empty value for required setting '{key}' in {tool_name} configuration"
                        )
                        return False
                # Validate executable path exists
                executable_path = Path(tool_config["executable"])
                if not executable_path.exists() or not executable_path.is_file():
                    logger.error(
                        f"Executable not found or is not a file for {tool_name} at {executable_path}"
                    )
                    return False
                return True
            case "siena_db":
                required_keys = ["executable", "siena_db", "pdb_directory", "format"]
                tool_config = self.dict.get(tool_name, {})

                for key in required_keys:
                    if key not in tool_config:
                        logger.error(
                            f"Missing required setting '{key}' in {tool_name} configuration"
                        )
                        return False
                    if key != "optional_arguments" and (
                        tool_config[key] is None or tool_config[key] == ""
                    ):
                        logger.error(
                            f"Empty value for required setting '{key}' in {tool_name} configuration"
                        )
                        return False

                # Validate executable path exists
                executable_path = Path(tool_config["executable"])
                if not executable_path.exists() or not executable_path.is_file():
                    logger.error(
                        f"Executable not found or is not a file for {tool_name} at {executable_path}"
                    )
                    return False

                # Validate siena_db path exists or can be created
                siena_db_path = Path(tool_config["siena_db"])
                if not siena_db_path.exists() and not siena_db_path.parent.exists():
                    logger.error(f"Invalid siena_db path: {siena_db_path}")
                    return False

                # Validate pdb_directory path exists
                pdb_directory_path = Path(tool_config["pdb_directory"])
                if not pdb_directory_path.exists() or not pdb_directory_path.is_dir():
                    logger.error(
                        f"PDB directory not found or is not a directory for {tool_name} at {pdb_directory_path}"
                    )
                    return False

                # Validate format is an integer (0 or 1, assuming based on typical usage)
                if not isinstance(
                    tool_config["format"], int
                ):  # or tool_config["format"] not in [0, 1]:
                    logger.error(
                        f"Invalid format value for {tool_name}: {tool_config['format']}. Must be an integer."
                    )
                    return False
                return True
            case _:
                logger.error(f"No validation rules defined for tool: {tool_name}")
                return False

    def _fetch_tool_config(self, tool_name: str) -> dict:
        tool_config_dict = self.dict.get(tool_name, None)

        if tool_config_dict is not None:
            logger.debug(f"Found configuration for {tool_name}")
            if self._validate_config(tool_name):
                return tool_config_dict
            else:
                logger.error(f"Invalid configuration for tool '{tool_name}'")
                raise ValueError(f"Invalid configuration for tool '{tool_name}'")
        else:
            logger.error(
                f"Configuration for tool '{tool_name}' not found in config file"
            )
            raise ValueError(
                f"Configuration for tool '{tool_name}' not found in config file"
            )

    @property
    def uniprot_id(self) -> str:
        """The UniProt ID from general settings."""
        return self.general_settings["uniprot_id"]

    @property
    def output_dir(self) -> Path:
        """The base output directory as a Path object from general settings."""
        return Path(self.general_settings["output_dir"])

    @property
    def log_level(self) -> str:
        """The log level from general settings."""
        return self.general_settings["log_level"]

    @property
    def log_file(self) -> str:
        """The log file name from general settings."""
        return self.general_settings["log_file"]

    def __repr__(self):
        """Return a string representation of the Config object."""
        return f"Config(config_path={self.config_path})"

    def __str__(self):
        """Return a readable string representation of the config."""
        return json.dumps(self.dict, indent=4, default=str)


if __name__ == "__main__":
    c = Config(Path("config.toml"))
    print(c.ligand_extractor)
