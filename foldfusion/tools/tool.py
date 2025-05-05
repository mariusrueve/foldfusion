from pathlib import Path
import logging
import subprocess

logger = logging.getLogger(__name__)


class Tool:
    def __init__(self, config: dict):
        self.config: dict = config
        self.tool_name: str = None
        self.executable: Path = None
        self.required_arguments: str = None
        self.optional_arguments: str = None
        self.output_dir = Path(self.config.get("run_settings", "")["output_dir"])
        self.command: str = None

    def load_executable_config(self, tool_name: str):
        self.tool_name = tool_name
        tool_config = self.config.get(self.tool_name, {})
        self.executable = Path(tool_config.get("executable", ""))
        optional_arguments = tool_config.get("optional_arguments", "")
        # Parse optional arguments
        self.optional_arguments = []
        if optional_arguments:
            if isinstance(optional_arguments, list):
                self.optional_arguments = optional_arguments
            else:
                self.optional_arguments = str(optional_arguments).split()

    def run(self) -> Path:
        # Ensure executable is set
        if not self.executable or not self.executable.exists():
            raise FileNotFoundError(
                f"Executable for {self.tool_name} not found: {self.executable}"
            )

        print(self.output_dir)
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Log the command
        logger.info(f"Running command: {' '.join(str(c) for c in self.command)}")

        # Execute the command
        try:
            _ = subprocess.run(
                self.command, cwd=self.output_dir, stdout=subprocess.DEVNULL
            )
            logger.info(f"Command completed successfully")

            # Return the output directory
            return self.output_dir
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed with exit code {e.returncode}")
            logger.error(f"Error output: {e.stderr}")
            raise RuntimeError(f"Failed to run {self.tool_name}: {e}")
