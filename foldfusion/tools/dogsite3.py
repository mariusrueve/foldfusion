from foldfusion.tools.tool import Tool
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class Dogsite3(Tool):
    def __init__(self, config: dict, pdb_file: Path):
        super().__init__(config)
        self.load_executable_config("dogsite3")
        self.pdb_file = pdb_file
        self.output_dir = self.output_dir / "dogsite3"
        self.assemble_command_arguments()

    def assemble_command_arguments(self):
        command = ["--proteinFile", str(self.pdb_file)]
        self.command = [self.executable] + command + self.optional_arguments

    def get_best_edf(self) -> Path:
        edf_filename = "output_P_1_res.edf"
        edf_file_path = self.output_dir / edf_filename

        if not edf_file_path.is_file():
            raise FileNotFoundError(
                f"EDF file not found: {edf_file_path}. "
                f"Ensure the Dogsite3 tool has been run and produced this output file."
            )

        lines = edf_file_path.read_text().splitlines()
        new_lines = []
        modified = False
        # self.pdb_file is a Path object, convert to string for use in f-string
        reference_pdb_path_str = str(self.pdb_file)

        for line in lines:
            # Strip whitespace from the line for robust comparison
            if line.strip() == "REFERENCE <NO-FILE>":
                new_lines.append(f"REFERENCE {reference_pdb_path_str}")
                modified = True
            else:
                new_lines.append(line)

        if modified:
            # Join lines with newline character and add a final newline
            content_to_write = "\n".join(new_lines) + "\n"
            edf_file_path.write_text(content_to_write)
        logger.info(f"The best EDF file: {Path.cwd() / edf_file_path}")
        return Path.cwd() / edf_file_path
