import logging
from pathlib import Path

from .tool import Tool

logger = logging.getLogger(__name__)


class JamdaScorer(Tool):
    def __init__(
        self,
        executable: Path,
        alphafold_pdb_path: Path,
        ligand_structure: dict,
        output_dir: Path,
    ):
        self.executable = executable
        self.alphafold_pdb_path = alphafold_pdb_path
        self.ligand_structure = ligand_structure
        self.output_dir = output_dir / "JamdaScorer"

    def run(self):
        optimized_ligand_structure = {}
        for pdb_code, ligands in self.ligand_structure.items():
            ligand_output_dir = self.output_dir / pdb_code
            ligand_output_dir.mkdir(parents=True, exist_ok=True)
            optimized_ligand_structure[pdb_code] = []
            for ligand in ligands:
                sdf_file = Path(ligand["path"])
                ligand_id = ligand["ligand_id"]
                output_sdf = ligand_output_dir / f"{ligand_id}.sdf"
                self.command = [
                    self.executable,
                    "-i",
                    self.alphafold_pdb_path,
                    "-m",
                    sdf_file,
                    "-o",
                    output_sdf,
                    "--optimize",
                ]
                super().run()
                optimized_ligand_structure[pdb_code].append({
                    "ligand_id": ligand_id,
                    "path": str(output_sdf.absolute()),
                    "sdf_file": f"{ligand_id}.sdf",
                })
        self.optimized_ligand_structure = optimized_ligand_structure
        return optimized_ligand_structure
