from foldfusion.tools.tool import Tool
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class LigandExtractor(Tool):
    def __init__(self, config, siena_dir: Path, pdb_code_list: list[str]):
        super().__init__(config)
        self.load_executable_config("ligand_extractor")
        self.output_dir = self.output_dir / "ligand_extractor"
        self.siena_dir = siena_dir
        self.pdb_code_list = pdb_code_list

    def _get_siena_pdb_path(self, pdb_code: str):
        ensemble_dir = Path.cwd() / self.siena_dir / "ensemble"
        pdb_files = list(ensemble_dir.glob("*.pdb"))

        # Find the one that contains pdb_code
        for pdb_file in pdb_files:
            if pdb_code.lower() in pdb_file.name.lower():
                return pdb_file.absolute()

        # If no matching file is found
        raise FileNotFoundError(f"No PDB file found for {pdb_code} in {ensemble_dir}")

    def _get_ligand_ids(self, pdb_path: Path, wanted_chain_id: str) -> list[str]:
        ligands = []

        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("HET "):
                    parts = line.split()
                    ligand_name = parts[1]
                    chain_id = parts[2]
                    residue_number = parts[3]
                    # Check if the chain ID contains any numbers
                    has_number_in_chain = any(char.isdigit() for char in chain_id)
                    if has_number_in_chain:
                        logger.warning(
                            f"Found number in chain ID '{chain_id}' for ligand "
                            + f"{ligand_name} at position {residue_number}"
                        )
                        # Split the chain ID after the first character
                        first_char = chain_id[0]
                        remaining_nums = chain_id[1:]
                        # Update residue number to include the numbers from chain ID
                        # Original line: residue_number = remaining_nums + residue_number
                        # Corrected line: Assign remaining_nums as the residue number
                        residue_number = remaining_nums
                        # Update chain ID to just the first character
                        chain_id = first_char
                        logger.info(
                            f"Split chain ID into '{chain_id}' and updated residue "
                            + f"number to {residue_number}"
                        )
                    if not chain_id == wanted_chain_id:
                        continue
                    ligands.append(f"{ligand_name}_{chain_id}_{residue_number}")
        return ligands

    def _get_commands_list(self) -> list:
        commands_list = []
        for code, chain in self.pdb_code_list:
            siena_pdb_path = self._get_siena_pdb_path(code)
            ligand_ids = self._get_ligand_ids(siena_pdb_path, chain)

            for id in ligand_ids:
                commands_list.append(
                    [
                        str(self.executable),
                        "-c",
                        self._get_siena_pdb_path(code),
                        "-l",
                        id,
                        "-o",
                        code,
                    ]
                )
        return commands_list

    def run_all(self):
        commands = self._get_commands_list()
        for command in commands:
            self.command = command
            self.run()
