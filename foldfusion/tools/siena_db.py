from foldfusion.tools.tool import Tool
from pathlib import Path


class SienaDB(Tool):
    def __init__(self, config):
        super().__init__(config)
        self.load_executable_config("siena_db")
        self.output_dir = self.output_dir / "siena_db"
        self.assemble_command_arguments()
        

    def assemble_command_arguments(self):
        
    

