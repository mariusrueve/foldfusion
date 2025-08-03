"""
Logging Configuration Module

This module provides centralized logging configuration for the FoldFusion pipeline.
It sets up both console and file logging with rotation capabilities to manage
log file sizes and maintain system performance.

Features:
    - Configurable log levels (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    - Rotating file handlers to prevent excessive disk usage
    - Consistent formatting across all pipeline components
    - UTF-8 encoding support for international characters
    - Automatic log directory creation

The logging system is designed to provide comprehensive traceability of pipeline
execution while maintaining reasonable performance characteristics.
"""

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path


def setup_logging(log_level: str, log_file: Path) -> None:
    """
    Configure centralized logging for the FoldFusion pipeline.

    This function sets up both console and file logging handlers with consistent
    formatting and rotation policies. It ensures that log messages are properly
    captured throughout the pipeline execution.

    Args:
        log_level: Logging level as string (DEBUG, INFO, WARNING, ERROR, CRITICAL).
                  Invalid levels will default to INFO with a warning message.
        log_file: Path to the log file. Parent directories will be created if
                 they don't exist. The file will use rotating handlers to manage
                 disk space usage.

    Raises:
        OSError: If the log file path cannot be created or is not writable.
        ValueError: If the log_level parameter is None or empty.

    Note:
        This function removes any existing handlers to prevent duplication when
        called multiple times. Log files are rotated when they exceed 10MB,
        keeping up to 5 backup files.
    """
    if not log_level:
        raise ValueError("Log level cannot be None or empty")

    # Convert string level to numeric, defaulting to INFO for invalid levels
    numeric_level = getattr(logging, log_level.upper(), None)
    if numeric_level is None:
        print(f"Warning: Invalid log level '{log_level}', defaulting to INFO")
        numeric_level = logging.INFO

    # Standardized log format with timestamp, module name, level, and message
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(numeric_level)

    # Clear existing handlers to avoid duplication during reconfiguration
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Console handler for immediate feedback
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(logging.Formatter(log_format))
    stream_handler.setLevel(numeric_level)
    root_logger.addHandler(stream_handler)

    try:
        # Ensure log file directory exists
        log_file.parent.mkdir(parents=True, exist_ok=True)

        # Rotating file handler to manage log file sizes
        # Maximum 10MB per file, keep 5 backup files
        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=10 * 1024 * 1024,  # 10MB
            backupCount=5,
            encoding="utf-8"
        )
        file_handler.setFormatter(logging.Formatter(log_format))
        file_handler.setLevel(numeric_level)
        root_logger.addHandler(file_handler)

    except (OSError, PermissionError) as e:
        # Log to console if file logging fails
        root_logger.error(
            f"Failed to set up file logging to {log_file}: {e}. "
            "Continuing with console logging only."
        )

    # Log successful configuration
    logger = logging.getLogger(__name__)
    logger.info(
        f"Logging system initialized - Level: {log_level.upper()}, "
        f"File: {log_file}, Handlers: {len(root_logger.handlers)}"
    )
