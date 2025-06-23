from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np


def parse_pdb(pdb_file: Path) -> Dict:
    """
    Parse PDB file and extract protein structure information.

    Args:
        pdb_file: Path to PDB file

    Returns:
        Dictionary containing parsed PDB data
    """
    atoms = []
    header_info = {}
    hetero_atoms = []

    with open(pdb_file, "r") as f:
        for line in f:
            record_type = line[:6].strip()

            if record_type == "HEADER":
                header_info["classification"] = line[10:50].strip()
                header_info["deposition_date"] = line[50:59].strip()
                header_info["id_code"] = line[62:66].strip()

            elif record_type == "TITLE":
                title = header_info.get("title", "")
                header_info["title"] = title + " " + line[10:80].strip()

            elif record_type == "ATOM":
                atom_data = {
                    "atom_number": int(line[6:11]),
                    "atom_name": line[12:16].strip(),
                    "residue_name": line[17:20].strip(),
                    "chain_id": line[21].strip(),
                    "residue_number": int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "occupancy": float(line[54:60]) if line[54:60].strip() else 1.0,
                    "temp_factor": float(line[60:66]) if line[60:66].strip() else 0.0,
                    "element": line[76:78].strip(),
                }
                atoms.append(atom_data)

            elif record_type == "HETATM":
                hetatm_data = {
                    "atom_number": int(line[6:11]),
                    "atom_name": line[12:16].strip(),
                    "residue_name": line[17:20].strip(),
                    "chain_id": line[21].strip(),
                    "residue_number": int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "occupancy": float(line[54:60]) if line[54:60].strip() else 1.0,
                    "temp_factor": float(line[60:66]) if line[60:66].strip() else 0.0,
                    "element": line[76:78].strip(),
                }
                hetero_atoms.append(hetatm_data)

    return {
        "header": header_info,
        "atoms": atoms,
        "hetero_atoms": hetero_atoms,
        "num_atoms": len(atoms),
        "num_hetero_atoms": len(hetero_atoms),
    }


def parse_sdf(sdf_file: Path) -> Dict:
    """
    Parse SDF file and extract molecular structure information.

    Args:
        sdf_file: Path to SDF file

    Returns:
        Dictionary containing parsed SDF data compatible with get_coordinates
    """
    atoms = []
    header_info = {}

    with open(sdf_file, "r") as f:
        lines = f.readlines()

    # Parse header (first 3 lines)
    if len(lines) >= 3:
        header_info["title"] = lines[0].strip()
        header_info["software"] = lines[1].strip()
        header_info["comment"] = lines[2].strip()

    # Parse counts line (4th line)
    if len(lines) >= 4:
        counts_line = lines[3]
        num_atoms = int(counts_line[:3])
        num_bonds = int(counts_line[3:6])

        # Parse atom block
        for i in range(4, 4 + num_atoms):
            if i < len(lines):
                line = lines[i]
                atom_data = {
                    "atom_number": i - 3,
                    "atom_name": line[31:34].strip(),
                    "residue_name": "",
                    "chain_id": "",
                    "residue_number": 1,
                    "x": float(line[:10]),
                    "y": float(line[10:20]),
                    "z": float(line[20:30]),
                    "occupancy": 1.0,
                    "temp_factor": 0.0,
                    "element": line[31:34].strip(),
                }
                atoms.append(atom_data)

    return {
        "header": header_info,
        "atoms": atoms,
        "hetero_atoms": [],
        "num_atoms": len(atoms),
        "num_hetero_atoms": 0,
    }


def get_coordinates(
    parsed_data: Dict, atom_types: Optional[List[str]] = None
) -> np.ndarray:
    """
    Extract coordinates from parsed PDB data.

    Args:
        parsed_data: Dictionary from parse_pdb function
        atom_types: List of atom types to filter (e.g., ['CA', 'CB']). If None, returns all atoms.

    Returns:
        numpy array of shape (n_atoms, 3) containing x, y, z coordinates
    """
    atoms = parsed_data["atoms"]

    if atom_types:
        filtered_atoms = [atom for atom in atoms if atom["atom_name"] in atom_types]
    else:
        filtered_atoms = atoms

    coords = np.array([[atom["x"], atom["y"], atom["z"]] for atom in filtered_atoms])
    return coords
