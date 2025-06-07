import logging
from logging.handlers import RotatingFileHandler  # For rotating log files
from pathlib import Path


def setup_logging(log_level: str, log_file: Path):
    """Configures logging based on the provided configuration."""
    numeric_level = getattr(logging, log_level, logging.INFO)

    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    # Get the root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(numeric_level)  # Set level on root logger

    # Remove any existing handlers to avoid duplication if main is called multiple times
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Console Handler
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(logging.Formatter(log_format))
    root_logger.addHandler(stream_handler)

    # File Handler (Rotating)
    file_handler = RotatingFileHandler(
        log_file, maxBytes=10 * 1024 * 1024, backupCount=5, encoding="utf-8"
    )
    file_handler.setFormatter(logging.Formatter(log_format))
    root_logger.addHandler(file_handler)

    # Get a logger for the main module itself
    logger = logging.getLogger(__name__)
    logger.info(f"Logging configured: Level={log_level}, File={log_file}")
