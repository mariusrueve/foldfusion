import json
import logging
from pathlib import Path

from foldfusion.evaluation.utils import parse_sdf

from .tool import Tool

logger = logging.getLogger(__name__)


class LigandExtractor(Tool):
    """Extract ligands from SIENA ensembles and organise the results."""

    def __init__(
        self,
        executable: Path,
        siena_dir: Path,
        alignment_list: list,
        output_dir: Path,
    ):
        """Persist configuration for ligand extraction runs.

        Args:
            executable: Path to the ligand extraction CLI.
            siena_dir: SIENA output directory that contains the ensemble files.
            alignment_list: Alignments produced by ``Siena.get_best_alignments``.
            output_dir: Pipeline output directory; a ``LigandExtractor/`` folder
                is created inside it.
        """
        self.executable = executable
        self.siena_dir = siena_dir
        self.alignment_list = alignment_list
        self.output_dir = output_dir / "LigandExtractor"

    def _get_siena_pdb_path(self, pdb_code: str) -> Path:
        """Return the ensemble PDB path that matches a SIENA alignment.

        Raises:
            FileNotFoundError: If the ensemble directory does not contain a
                matching PDB file.
        """
        ensemble_dir = Path.cwd() / self.siena_dir / "ensemble"
        pdb_files = list(ensemble_dir.glob("*.pdb"))

        # Find the one that contains pdb_code
        for pdb_file in pdb_files:
            if pdb_code.lower() in pdb_file.name.lower():
                return pdb_file.absolute()

        # If no matching file is found
        raise FileNotFoundError(f"No PDB file found for {pdb_code} in {ensemble_dir}")

    def _get_ligand_ids(self, pdb_path: Path, wanted_chain_id: str) -> list[str]:
        """Collect ligand identifiers for the chain associated with an alignment.

        Ligand IDs follow the ``<ligand>_<chain>_<resnum>`` pattern used by the
        extractor CLI.
        """
        ligands = []

        with open(pdb_path) as f:
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
                            "Found number in chain ID '%s' for ligand %s at "
                            "position %s. Applying fixes. Please check for "
                            "correctness.",
                            chain_id,
                            ligand_name,
                            residue_number,
                        )
                        # Split the chain ID after the first character
                        first_char = chain_id[0]
                        remaining_nums = chain_id[1:]
                        residue_number = remaining_nums
                        # Update chain ID to just the first character
                        chain_id = first_char
                        logger.debug(
                            "Split chain ID into '%s' and updated residue number to %s",
                            chain_id,
                            residue_number,
                        )
                    if not chain_id == wanted_chain_id:
                        continue
                    ligands.append(f"{ligand_name}_{chain_id}_{residue_number}")
        return ligands

    def _get_commands_list(self) -> list[list[str]]:
        """Build the command lines needed for each ligand extraction."""
        commands_list: list[list[str]] = []
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

            # Safely extract ligand IDs for the given alignment; log and skip on error
            try:
                ligand_ids = self._get_ligand_ids(ensemble_path, chain)
            except Exception as e:
                logger.warning(
                    "Failed to get ligand IDs for %s (chain %s): %s. "
                    "Skipping alignment.",
                    code,
                    chain,
                    e,
                )
                continue

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
        """Execute ligand extraction for each SIENA hit and cache the results.

        Returns:
            The absolute path to the ``LigandExtractor/`` output directory. The
            method also populates ``self.ligand_structure`` with the extracted
            ligand metadata and persists it as JSON for reuse by later stages.
        """
        # Build commands list defensively; on failure, continue with empty list
        try:
            commands = self._get_commands_list()
        except Exception as e:
            logger.error("Failed to build ligand extraction command list: %s", e)
            commands = []
        if not commands:
            logger.warning(
                "No ligand extraction commands generated. Continuing with "
                "empty results."
            )
        ligand_structure: dict[str, list[dict[str, str]]] = {}

        for command in commands:
            self.command = command
            # Execute each command independently; on failure, log and continue
            try:
                super().run()
            except Exception as e:
                logger.error(
                    "Ligand extraction command failed: %s | Error: %s",
                    " ".join(map(str, command)),
                    e,
                )
                continue

            # Extract ligand information from the command
            pdb_code = command[-1]  # The last element in the command is the pdb_code
            ligand_id = command[4]  # The ligand ID is the fourth element in the command
            output_path = self.output_dir / pdb_code / f"{ligand_id}.sdf"

            # Validate SDF file
            try:
                parsed = parse_sdf(output_path)
                if parsed["num_atoms"] == 0:
                    logger.warning(
                        "Extracted ligand SDF is empty: %s. Skipping.", output_path
                    )
                    continue
            except Exception as e:
                logger.warning(
                    "Failed to parse or locate SDF %s: %s. Skipping.",
                    output_path,
                    e,
                )
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

        logger.info("Ligand structure saved to %s", json_output_path)

        return self.output_dir
