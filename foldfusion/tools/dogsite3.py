from foldfusion.tools.tool import Tool
from pathlib import Path


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
