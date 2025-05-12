"""Base class for all external tool wrappers in the FoldFusion pipeline."""

from pathlib import Path
import logging
import subprocess
from typing import Optional, List  # Corrected: Removed unused Union

logger = logging.getLogger(__name__)


class Tool:
    """A base class for wrapping external command-line tools.

    This class provides common functionality for configuring, running, and managing
    output for external bioinformatics tools.

    Attributes:
        config (dict): The main configuration dictionary for the pipeline.
        tool_name (Optional[str]): The name of the tool, used for logging and config lookup.
        executable (Optional[Path]): Path to the tool's executable.
        optional_arguments (List[str]): A list of optional arguments for the tool command.
        output_dir (Path): The base directory where the tool will write its output.
        command (Optional[List[str]]): The fully assembled command to be executed.
    """

    def __init__(self, config: dict):
        """Initializes the Tool class.

        Args:
            config (dict): The configuration dictionary.
        """
        if not isinstance(config, dict):
            logger.error("Tool initialized with invalid config type. Expected dict.")
            raise TypeError("Configuration must be a dictionary.")
        self.config: dict = config
        self.tool_name: Optional[str] = None
        self.executable: Optional[Path] = None
        # self.required_arguments: Optional[str] = None # This was unused, consider removing if not planned
        self.optional_arguments: List[str] = []

        run_settings = self.config.get("run_settings", {})
        output_dir_str = run_settings.get("output_dir", "output")  # Default to "output"
        self.output_dir = Path(output_dir_str).resolve()
        # Base output directory for the entire pipeline is created here if it doesn't exist.
        # Specific tool output directories (e.g., self.output_dir / self.tool_name)
        # should be created by the derived classes or within their run methods.
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Base output directory for tools ensured at: {self.output_dir}")

        self.command: Optional[List[str]] = None
        logger.debug(f"Tool base class initialized. Main output dir: {self.output_dir}")

    def load_executable_config(self, tool_name: str):
        """Loads tool-specific configuration, including executable path and optional arguments.

        Args:
            tool_name (str): The name of the tool (key in the configuration dictionary).

        Raises:
            ValueError: If the tool configuration is not found or executable is not specified.
            FileNotFoundError: If the specified executable path does not exist.
        """
        self.tool_name = tool_name
        tool_config = self.config.get(self.tool_name)
        if not tool_config or not isinstance(tool_config, dict):
            logger.error(
                f"Configuration for tool '{self.tool_name}' not found or is not a dictionary."
            )
            raise ValueError(
                f"Configuration for tool '{self.tool_name}' is missing or invalid."
            )

        executable_path_str = tool_config.get("executable")
        if not executable_path_str:
            logger.error(
                f"Executable path not specified for tool '{self.tool_name}' in configuration."
            )
            raise ValueError(
                f"Executable for '{self.tool_name}' must be specified in config."
            )

        self.executable = Path(executable_path_str).resolve()
        if not self.executable.is_file():  # More specific check for files
            logger.error(
                f"Executable for '{self.tool_name}' not found or is not a file at: {self.executable}"
            )
            raise FileNotFoundError(
                f"Executable for '{self.tool_name}' not found at {self.executable}"
            )

        optional_args_config = tool_config.get(
            "optional_arguments", []
        )  # Default to empty list
        if isinstance(optional_args_config, list):
            self.optional_arguments = [
                str(arg) for arg in optional_args_config
            ]  # Ensure all args are strings
        elif isinstance(optional_args_config, str):
            self.optional_arguments = optional_args_config.split()
        else:
            logger.warning(
                f"Optional arguments for '{self.tool_name}' are not a list or string. Ignoring. "
                f"Found type: {type(optional_args_config)}"
            )
            self.optional_arguments = []
        logger.debug(
            f"Loaded executable for {self.tool_name}: {self.executable}. "
            f"Optional args: {self.optional_arguments}"
        )

    def run(self) -> Path:
        """Runs the assembled command for the tool.

        The command is executed in the tool's specific output directory.
        (e.g., `pipeline_output_base/tool_name/`).

        Returns:
            Path: The path to the output directory where the tool was run.

        Raises:
            FileNotFoundError: If the executable was not set or not found prior to running.
            ValueError: If the command was not assembled before calling run.
            RuntimeError: If the command execution fails.
        """
        if (
            not self.executable or not self.executable.exists()
        ):  # Should be caught by load_executable_config
            logger.error(
                f"Executable for {self.tool_name or 'Unknown tool'} not found or not set: {self.executable}"
            )
            raise FileNotFoundError(
                f"Executable for {self.tool_name or 'Unknown tool'} not found: {self.executable}"
            )

        if not self.command:
            logger.error(
                f"Command for {self.tool_name or 'Unknown tool'} has not been assembled. Call assemble_command() first."
            )
            raise ValueError(
                f"Command not set for {self.tool_name or 'Unknown tool'}. Cannot run."
            )

        # Determine the specific output directory for this tool instance.
        # This is typically self.output_dir / self.tool_name, as set by derived classes.
        # For the base Tool class, if a derived class hasn't overridden self.output_dir,
        # it will use the main pipeline output directory. It's better if derived classes
        # set their specific subdirectories (e.g., self.output_dir = self.output_dir / self.tool_name).
        # The current implementation in derived classes already does this (e.g. self.output_dir / "alphafold").
        # So, self.output_dir here should already be the tool-specific one.
        current_tool_output_dir = (
            self.output_dir
        )  # This should be the tool-specific path
        current_tool_output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(
            f"Ensured tool-specific output directory exists: {current_tool_output_dir}"
        )

        # Log the command, ensuring all parts are strings for join
        command_str = " ".join(map(str, self.command))
        logger.info(f"Running {self.tool_name} in directory {current_tool_output_dir}:")
        logger.info(f"COMMAND: {command_str}")

        try:
            # Using subprocess.run
            # stdout=subprocess.DEVNULL and stderr=subprocess.PIPE can be used to capture stderr
            # For more detailed logging, one might capture stdout/stderr to files or logger.
            process = subprocess.run(
                self.command,
                cwd=current_tool_output_dir,  # Run in the tool's specific output directory
                check=True,  # Raise CalledProcessError on non-zero exit codes
                capture_output=True,  # Capture stdout and stderr
                text=True,  # Decode stdout and stderr as text
                encoding="utf-8",  # Specify encoding
            )
            logger.info(f"{self.tool_name} completed successfully.")
            if process.stdout:
                logger.debug(f"{self.tool_name} STDOUT:\n{process.stdout.strip()}")
            if (
                process.stderr
            ):  # Log stderr even on success, as it might contain warnings
                logger.debug(f"{self.tool_name} STDERR:\n{process.stderr.strip()}")

            return current_tool_output_dir.resolve()
        except subprocess.CalledProcessError as e:
            logger.error(
                f"{self.tool_name} command failed with exit code {e.returncode}."
            )
            if e.stdout:
                logger.error(f"STDOUT from {self.tool_name}:\n{e.stdout.strip()}")
            if e.stderr:
                logger.error(f"STDERR from {self.tool_name}:\n{e.stderr.strip()}")
            raise RuntimeError(
                f"Failed to run {self.tool_name}. Exit code: {e.returncode}. "
                f"Command: '{command_str}'. Error: {e.stderr.strip() if e.stderr else 'N/A'}"
            ) from e
        except (
            FileNotFoundError
        ) as e:  # If the executable itself is not found at runtime by subprocess
            logger.error(
                f"Executable '{self.command[0]}' for {self.tool_name} not found during subprocess execution. {e}"
            )
            raise RuntimeError(
                f"Executable '{self.command[0]}' for {self.tool_name} not found. "
                f"Please check the path and permissions."
            ) from e
        except Exception as e:
            logger.error(
                f"An unexpected error occurred while running {self.tool_name} command '{command_str}': {e}"
            )
            raise RuntimeError(f"Unexpected error running {self.tool_name}: {e}") from e
