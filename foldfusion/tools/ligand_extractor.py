import json
import logging
from pathlib import Path

from .tool import Tool
from foldfusion.evaluation.utils import parse_sdf

logger = logging.getLogger(__name__)


class LigandExtractor(Tool):
    """
    LigandExtractor
    ===============

    This class is responsible for extracting ligand information from PDB files and organizing the output in a structured format.

    Attributes:
    -----------
    executable : Path
        Path to the executable used for ligand extraction.
    siena_dir : Path
        Directory containing the Siena ensemble files.
    pdb_code_list : list
        List of dictionaries containing alignment data with keys:
        - pdb_code: PDB structure code
        - chain: PDB chain identifier
        - ensemble_path: Path to the ensemble file
        - ligand_pdb_code: Ligand PDB code
        - backbone_rmsd: Backbone RMSD value
        - all_atom_rmsd: All atom RMSD value
    output_dir : Path
        Directory where the extracted ligand files will be stored.
    ligand_structure : dict
        Dictionary containing the extracted ligands organized by PDB code.
        Structure: {pdb_code: [{"ligand_id": str, "path": str, "sdf_file": str}, ...]}

    Methods:
    --------
    __init__(executable: Path, siena_dir: Path, alignment_list: list, output_dir: Path)
        Initializes the LigandExtractor with the required parameters.

        Args:
            executable: Path to the ligand extraction executable
            siena_dir: Directory containing the Siena ensemble files
            alignment_list: List of alignment dictionaries from Siena.get_best_alignments()
            output_dir: Directory where extracted ligand files will be stored

    _get_siena_pdb_path(pdb_code: str) -> Path
        Retrieves the absolute path of the PDB file matching the given PDB code.

    _get_ligand_ids(pdb_path: Path, wanted_chain_id: str) -> list[str]
        Extracts ligand IDs from the PDB file for the specified chain ID.

    _get_commands_list() -> list
        Generates a list of commands to execute for ligand extraction.

    run() -> Path
        Executes the ligand extraction process and organizes the output.
    """

    def __init__(
        self,
        executable: Path,
        siena_dir: Path,
        alignment_list: list,
        output_dir: Path,
    ):
        self.executable = executable
        self.siena_dir = siena_dir
        self.alignment_list = alignment_list
        self.output_dir = output_dir / "LigandExtractor"

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
                            f"{ligand_name} at position {residue_number}. "
                            "Applying fixes. Please check for correctnes."
                        )
                        # Split the chain ID after the first character
                        first_char = chain_id[0]
                        remaining_nums = chain_id[1:]
                        residue_number = remaining_nums
                        # Update chain ID to just the first character
                        chain_id = first_char
                        logger.debug(
                            f"Split chain ID into '{chain_id}' and updated residue "
                            + f"number to {residue_number}"
                        )
                    if not chain_id == wanted_chain_id:
                        continue
                    ligands.append(f"{ligand_name}_{chain_id}_{residue_number}")
        return ligands

    def _get_commands_list(self) -> list:
        commands_list = []
        for alignment in self.alignment_list:
            if isinstance(alignment, dict):
                code = alignment["pdb_code"]
                chain = alignment["chain"]
                ensemble_path = alignment["ensemble_path"]
            else:
                raise DeprecationWarning(
                    "Tuple format for alignment is deprecated. Please use a dictionary "
                    "with keys: 'pdb_code', 'chain', and 'ensemble_path'."
                )

            # Ensure ensemble_path is a Path object
            if isinstance(ensemble_path, str):
                ensemble_path = Path(ensemble_path)

            ligand_ids = self._get_ligand_ids(ensemble_path, chain)

            for id in ligand_ids:
                commands_list.append(
                    [
                        str(self.executable),
                        "-c",
                        str(ensemble_path),  # Convert to string for command
                        "-l",
                        id,
                        "-o",
                        code,
                    ]
                )
        return commands_list

    def run(self) -> Path:
        commands = self._get_commands_list()
        ligand_structure = {}

        for command in commands:
            self.command = command
            super().run()

            # Extract ligand information from the command
            pdb_code = command[-1]  # The last element in the command is the pdb_code
            ligand_id = command[4]  # The ligand ID is the fourth element in the command
            output_path = self.output_dir / pdb_code / f"{ligand_id}.sdf"

            # Validate SDF file
            try:
                parsed = parse_sdf(output_path)
                if parsed["num_atoms"] == 0:
                    logger.warning(
                        f"Extracted ligand SDF is empty: {output_path}. Skipping."
                    )
                    continue
            except Exception as e:
                logger.warning(f"Failed to parse SDF {output_path}: {e}. Skipping.")
                continue

            # Organize the structure in a nested dictionary with better structure
            if pdb_code not in ligand_structure:
                ligand_structure[pdb_code] = []

            ligand_structure[pdb_code].append(
                {
                    "ligand_id": ligand_id,
                    "path": str(output_path.absolute()),
                    "sdf_file": f"{ligand_id}.sdf",
                }
            )

        # Store the ligand structure as an attribute for later retrieval
        self.ligand_structure = ligand_structure

        # Save ligand structure as JSON in output directory
        json_output_path = self.output_dir / "ligand_structure.json"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        with open(json_output_path, "w") as json_file:
            json.dump(ligand_structure, json_file, indent=2)

        logger.info(f"Ligand structure saved to {json_output_path}")

        return self.output_dir
