"""
Base Tool Class for FoldFusion Pipeline

This module provides the base class for all external tool wrappers used in the
FoldFusion pipeline. It standardizes subprocess execution, error handling, and
logging across all bioinformatics tools.

Key Features:
    - Consistent subprocess execution with proper error handling
    - Comprehensive logging of command execution and output
    - Automatic output directory management
    - UTF-8 encoding support for international characters
    - Graceful handling of various failure modes

All tool implementations should inherit from this base class to ensure
consistent behavior and error reporting throughout the pipeline.
"""

import logging
import subprocess
from collections.abc import Callable
from typing import Any

try:
    # POSIX-only; used to enforce memory limits in child processes
    import resource  # type: ignore[attr-defined]
except Exception:  # pragma: no cover - non-POSIX environments
    resource = None  # type: ignore[assignment]
from pathlib import Path


def _stream_to_text(data: str | bytes | None) -> str | None:
    """Return decoded text for subprocess output regardless of input type."""

    if data is None:
        return None
    if isinstance(data, bytes):
        return data.decode("utf-8", errors="replace")
    return data


logger = logging.getLogger(__name__)


class Tool:
    """
    Base class for external bioinformatics tool wrappers in the FoldFusion pipeline.

    This class provides a standardized interface for executing external tools
    via subprocess, with comprehensive error handling and logging. It ensures
    consistent behavior across all tool implementations and provides robust
    error reporting for debugging purposes.

    All tool-specific classes should inherit from this base class and implement
    their specific command construction logic while leveraging the common
    execution infrastructure provided here.

    Attributes:
        command: List of command components to execute (executable + arguments)
        output_dir: Directory where tool output will be stored

    Example:
        Creating a custom tool wrapper:

        >>> class MyTool(Tool):
        ...     def __init__(
        ...         self, executable: Path, input_file: Path, output_dir: Path
        ...     ):
        ...         command = [str(executable), "--input", str(input_file)]
        ...         super().__init__(command, output_dir)
        ...
        >>> tool = MyTool("/usr/bin/mytool", "input.pdb", Path("output"))
        >>> result_dir = tool.run()
    """

    def __init__(self, command: list[Any], output_dir: Path) -> None:
        """
        Initialize the tool wrapper with command and output directory.

        Args:
            command: List containing the executable path and all arguments.
                    All elements will be converted to strings for subprocess.
            output_dir: Directory where tool output will be stored. Will be
                       created if it doesn't exist.

        Note:
            The command list should include the full path to the executable
            as the first element, followed by all required arguments.
        """
        self.command = [str(part) for part in command]
        self.output_dir = output_dir
        # Optional execution controls; may be set by callers after construction
        # Timeout in seconds for the tool execution
        self.timeout: int | None = None
        # Soft memory ceiling for the child process in megabytes (MB);
        # uses RLIMIT_AS when available (POSIX). None = no limit.
        self.memory_limit_mb: int | None = None

    def _preexec_with_limits(self) -> Callable[[], None] | None:
        """
        Build a pre-exec function that applies resource limits in the child.

        Returns a callable suitable for subprocess.run(preexec_fn=...), or None
        if no limits should be applied or platform support is unavailable.
        """
        if self.memory_limit_mb is None:
            return None

        if resource is None:
            logger.warning(
                "Resource limits not available on this platform; ignoring memory limit"
            )
            return None

        mem_bytes = int(self.memory_limit_mb) * 1024 * 1024

        def _apply_limits() -> None:
            try:
                # Address space limit: applies to virtual memory usage
                resource.setrlimit(  # type: ignore[attr-defined]
                    resource.RLIMIT_AS, (mem_bytes, mem_bytes)
                )
                # Data segment limit (best-effort)
                if hasattr(resource, "RLIMIT_DATA"):
                    resource.setrlimit(  # type: ignore[arg-type]
                        resource.RLIMIT_DATA, (mem_bytes, mem_bytes)
                    )
            except Exception:  # pragma: no cover - depends on OS
                # Log in parent after fork is not available; rely on failure behavior
                pass

        return _apply_limits

    def run(self) -> Path:
        """
        Execute the tool command and return the output directory path.

        This method handles the complete lifecycle of tool execution including:
        - Output directory creation
        - Command execution with proper subprocess configuration
        - Comprehensive logging of execution details
        - Error handling and reporting

        Returns:
            Absolute path to the output directory containing tool results.

        Raises:
            RuntimeError: If tool execution fails, executable is not found,
                         or any other execution error occurs. The exception
                         includes detailed error information for debugging.

        Note:
            All stdout and stderr output is captured and logged. The tool
            is executed in the specified output directory as the working
            directory.
        """
        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Created tool output directory: {self.output_dir}")

        # Prepare command string for logging (convert all parts to strings)
        command_str = " ".join(map(str, self.command))
        logger.info(f"Executing tool in directory: {self.output_dir}")
        logger.info(f"Command: {command_str}")

        try:
            # Execute the command with comprehensive configuration
            process = subprocess.run(
                self.command,
                cwd=self.output_dir,  # Set working directory
                check=True,  # Raise CalledProcessError on non-zero exit
                capture_output=True,  # Capture stdout and stderr
                text=True,  # Decode output as text
                encoding="utf-8",  # Specify encoding explicitly
                timeout=self.timeout if self.timeout is not None else None,
                preexec_fn=self._preexec_with_limits(),
                # Use a new session so timeouts don't leave orphan children
                start_new_session=True,
            )

            logger.debug("Tool execution completed successfully")

            # Log stdout if present
            if process.stdout:
                logger.debug(f"Tool stdout:\n{process.stdout.strip()}")

            # Log stderr even on success (may contain warnings)
            if process.stderr:
                logger.debug(f"Tool stderr:\n{process.stderr.strip()}")

            return self.output_dir.resolve()

        except subprocess.TimeoutExpired as e:
            # Kill the spawned process group to avoid orphans
            try:
                stdout_text = _stream_to_text(e.stdout)
                if stdout_text:
                    logger.error(f"Tool stdout before timeout:\n{stdout_text.strip()}")
                stderr_text = _stream_to_text(e.stderr)
                if stderr_text:
                    logger.error(f"Tool stderr before timeout:\n{stderr_text.strip()}")
            except Exception:
                pass
            error_msg = (
                f"Tool execution timed out after {self.timeout}s. "
                f"Command: '{command_str}'"
            )
            logger.error(error_msg)
            raise RuntimeError(error_msg) from e

        except subprocess.CalledProcessError as e:
            # Handle tool execution failure
            logger.error(f"Tool execution failed with exit code {e.returncode}")

            stdout_text = _stream_to_text(e.stdout)
            if stdout_text:
                logger.error(f"Tool stdout:\n{stdout_text.strip()}")
            stderr_text = _stream_to_text(e.stderr)
            if stderr_text:
                logger.error(f"Tool stderr:\n{stderr_text.strip()}")

            error_details = (
                e.stderr.strip() if e.stderr else "No error details available"
            )
            error_msg = (
                f"Tool execution failed. Exit code: {e.returncode}. "
                f"Command: '{command_str}'. Error: {error_details}"
            )
            raise RuntimeError(error_msg) from e

        except FileNotFoundError as e:
            # Handle executable not found
            logger.error(f"Executable '{self.command[0]}' not found: {e}")
            error_msg = (
                f"Executable '{self.command[0]}' not found. "
                f"Please verify the path and ensure the tool is properly installed."
            )
            raise RuntimeError(error_msg) from e

        except PermissionError as e:
            # Handle permission issues
            logger.error(f"Permission denied executing '{self.command[0]}': {e}")
            error_msg = (
                f"Permission denied for executable '{self.command[0]}'. "
                f"Please check file permissions and execution rights."
            )
            raise RuntimeError(error_msg) from e

        except Exception as e:
            # Handle any other unexpected errors
            logger.error(f"Unexpected error during tool execution: {e}")
            error_msg = f"Unexpected error executing '{command_str}': {str(e)}"
            raise RuntimeError(error_msg) from e
