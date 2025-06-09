"""Base class for all external tool wrappers in the FoldFusion pipeline."""

from pathlib import Path
import logging
import subprocess
from typing import List

logger = logging.getLogger(__name__)


class Tool:
    def __init__(self, command: List[str], output_dir: Path):
        self.command = command
        self.output_dir = output_dir

    def run(self) -> Path:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.debug(
            f"Ensured tool-specific output directory exists: {self.output_dir}"
        )

        # Log the command, ensuring all parts are strings for join
        command_str = " ".join(map(str, self.command))
        logger.debug(f"Running in directory {self.output_dir}:")
        logger.debug(f"COMMAND: {command_str}")

        try:
            process = subprocess.run(
                self.command,
                cwd=self.output_dir,  # Run in the tool's specific output directory
                check=True,  # Raise CalledProcessError on non-zero exit codes
                capture_output=True,  # Capture stdout and stderr
                text=True,  # Decode stdout and stderr as text
                encoding="utf-8",  # Specify encoding
            )
            logger.debug("Process completed successfully.")
            if process.stdout:
                logger.debug(f"STDOUT:\n{process.stdout.strip()}")
            if (
                process.stderr
            ):  # Log stderr even on success, as it might contain warnings
                logger.debug(f"STDERR:\n{process.stderr.strip()}")

            return self.output_dir.resolve()
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed with exit code {e.returncode}.")
            if e.stdout:
                logger.error(f"STDOUT:\n{e.stdout.strip()}")
            if e.stderr:
                logger.error(f"STDERR:\n{e.stderr.strip()}")
            raise RuntimeError(
                f"Failed to run. Exit code: {e.returncode}. "
                f"Command: '{command_str}'. Error: {e.stderr.strip() if e.stderr else 'N/A'}"
            ) from e
        except (
            FileNotFoundError
        ) as e:  # If the executable itself is not found at runtime by subprocess
            logger.error(
                f"Executable '{self.command[0]}' not found during subprocess execution. {e}"
            )
            raise RuntimeError(
                f"Executable '{self.command[0]}' not found. "
                f"Please check the path and permissions."
            ) from e
        except Exception as e:
            logger.error(
                f"An unexpected error occurred while running command '{command_str}': {e}"
            )
            raise RuntimeError(f"Unexpected error: {e}") from e
