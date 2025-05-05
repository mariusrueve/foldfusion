from pathlib import Path
import tomllib


def read_config(config_path: Path) -> dict:
    """
    Read and parse a TOML configuration file.
    This function attempts to load and parse a TOML configuration file from the specified path.
    The file is opened in binary mode as required by the tomllib parser.
    Args:
        config_path (Path): Path to the TOML configuration file.
    Returns:
        dict: Dictionary containing the parsed configuration.
    Raises:
        FileNotFoundError: If the configuration file doesn't exist at the specified path.
        ValueError: If the configuration file cannot be parsed.
    """

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found at {config_path}")

    try:
        with open(config_path, "rb") as f:  # TOML must be opened in binary mode
            config = tomllib.load(f)
        return config
    except Exception as e:
        raise ValueError(f"Failed to parse config file: {e}")
