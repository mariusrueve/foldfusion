"""Module for calculating quality scores for protein-ligand complexes."""

import numpy as np
from pathlib import Path
import re
import logging

logger = logging.getLogger(__name__)

# Bondi van der Waals radii (in Angstroms)
# Source: Bondi, A. (1964). J. Phys. Chem., 68(3), 441-451.
# Extended with a few common ions/metals with typical radii.
VDW_RADII = {
    "H": 1.20,
    "HE": 1.40,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "NE": 1.54,
    "SI": 2.10,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "AR": 1.88,
    "AS": 1.85,
    "SE": 1.90,
    "BR": 1.85,
    "KR": 2.02,
    "TE": 2.06,
    "I": 1.98,
    "XE": 2.16,
    "LI": 1.82,
    "NA": 2.27,
    "MG": 1.73,
    "AL": 1.84,  # AL can vary
    "K": 2.75,
    "CA": 2.31,
    "MN": 1.61,
    "FE": 1.56,
    "CO": 1.52,
    "NI": 1.49,
    "CU": 1.45,
    "ZN": 1.42,  # Metallic radii, can vary with charge/coordination
    # Default radius for unknown elements
    "UNKNOWN": 1.70,
}
DEFAULT_VDW_RADIUS = VDW_RADII["UNKNOWN"]


def _get_element_from_atom_name(atom_name):
    """Extracts element symbol from an atom name, guessing if needed."""
    # Remove digits and special characters from the start
    # e.g. " C1 " -> "C", "1HG" -> "H", "CA" -> "CA"
    cleaned_name = re.sub(r"^[^A-Z]*", "", atom_name.strip().upper())
    if not cleaned_name:
        return "UNKNOWN"
    if len(cleaned_name) >= 2 and cleaned_name[:2] in VDW_RADII:
        return cleaned_name[:2]
    if len(cleaned_name) >= 1 and cleaned_name[0] in VDW_RADII:
        return cleaned_name[0]
    return "UNKNOWN"


def parse_pdb_atoms(pdb_path: Path) -> list:
    """
    Parses ATOM and HETATM records from a PDB file.
    Returns a list of dictionaries: [{'element': str, 'coords': np.array([x,y,z])}].
    """
    atoms = []
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("ATOM  ") or line.startswith("HETATM"):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        # Element symbol is typically in columns 77-78
                        element_symbol = line[76:78].strip().upper()
                        if (
                            not element_symbol
                        ):  # If not in 77-78, try to guess from atom name (col 13-16)
                            atom_name = line[12:16].strip()
                            element_symbol = _get_element_from_atom_name(atom_name)

                        atoms.append(
                            {"element": element_symbol, "coords": np.array([x, y, z])}
                        )
                    except ValueError:
                        logger.warning(
                            f"Could not parse coordinates/element for line in {pdb_path}: {line.strip()}"
                        )
                        continue
    except FileNotFoundError:
        logger.error(f"PDB file not found for parsing: {pdb_path}")
        return []
    except Exception as e:
        logger.error(f"Error parsing PDB file {pdb_path}: {e}")
        return []
    return atoms


def parse_sdf_atoms(sdf_path: Path) -> list:
    """
    Parses atoms from an SDF file (very basic parser).
    Returns a list of dictionaries: [{'element': str, 'coords': np.array([x,y,z])}].
    WARNING: This is a simplified parser and may not work for all SDF variants.
             Using a dedicated library like RDKit is highly recommended for robustness.
    """
    atoms = []
    try:
        with open(sdf_path, "r") as f:
            lines = f.readlines()

        # Find the counts line (usually 4th line)
        # Example: "  6  5  0  0  0  0  0  0  0  0999 V2000"
        # The first number is num_atoms.
        counts_line_idx = -1
        for i in range(min(len(lines), 10)):  # Search in first few lines
            if "V2000" in lines[i] or "V3000" in lines[i]:
                # V2000: " aaabbb..." where aaa is num_atoms, bbb is num_bonds
                match = re.match(r"\s*(\d+)\s*(\d+)", lines[i])
                if match:
                    num_atoms = int(match.group(1))
                    counts_line_idx = i
                    break

        if counts_line_idx == -1 or num_atoms == 0:
            logger.warning(
                f"Could not find valid atom counts line or zero atoms in SDF: {sdf_path}"
            )
            return []

        # Atom block starts after the counts line
        # Each atom line: x y z element_symbol ...
        atom_block_start_idx = counts_line_idx + 1
        for i in range(num_atoms):
            line_idx = atom_block_start_idx + i
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                if len(parts) >= 4:
                    try:
                        x = float(parts[0])
                        y = float(parts[1])
                        z = float(parts[2])
                        element_symbol = parts[3].strip().upper()
                        atoms.append(
                            {"element": element_symbol, "coords": np.array([x, y, z])}
                        )
                    except ValueError:
                        logger.warning(
                            f"Could not parse atom line in SDF {sdf_path}: {lines[line_idx].strip()}"
                        )
                        continue
                else:
                    logger.warning(
                        f"Atom line too short in SDF {sdf_path}: {lines[line_idx].strip()}"
                    )
            else:
                logger.warning(
                    f"Expected more atom lines in SDF {sdf_path}, reached end of file."
                )
                break

    except FileNotFoundError:
        logger.error(f"SDF file not found for parsing: {sdf_path}")
        return []
    except Exception as e:
        logger.error(f"Error parsing SDF file {sdf_path}: {e}")
        return []

    if not atoms:
        logger.warning(f"No atoms parsed from SDF file: {sdf_path}")
    return atoms


def calculate_tcs(
    protein_atoms: list, ligand_atoms: list, distance_cutoff: float = 4.0
) -> float:
    """
    Calculates the Transplant Clash Score (TCS).

    Args:
        protein_atoms (list): List of protein atom dicts {'element': str, 'coords': np.array}.
        ligand_atoms (list): List of ligand atom dicts {'element': str, 'coords': np.array}.
        distance_cutoff (float): Cutoff distance in Angstroms to consider atom pairs.

    Returns:
        float: The calculated TCS score. Returns 0.0 if no atoms are within cutoff
               or no clashes are found.
    """
    sum_squared_overlaps = 0.0
    num_distances_considered = 0

    if not protein_atoms or not ligand_atoms:
        logger.warning("Protein or ligand atom list is empty for TCS calculation.")
        return 0.0

    for l_atom in ligand_atoms:
        # Use .get with a default for element if not in VDW_RADII
        r_l = VDW_RADII.get(l_atom["element"], DEFAULT_VDW_RADIUS)
        l_coords = l_atom["coords"]

        for p_atom in protein_atoms:
            r_p = VDW_RADII.get(p_atom["element"], DEFAULT_VDW_RADIUS)
            p_coords = p_atom["coords"]

            dist = np.linalg.norm(l_coords - p_coords)

            if dist < distance_cutoff:
                num_distances_considered += 1
                # Check for actual clash
                # overlap = (r_l + r_p) - dist
                # Using a slight tolerance for VdW radii sum, as exact radii can vary
                # and we want to penalize significant clashes.
                # Let's consider a clash if dist < (r_l + r_p - tolerance)
                # For TCS as per paper, overlap is (r_l + r_p) - dist
                # Only positive overlaps are squared and summed.

                combined_radii = r_l + r_p
                overlap = combined_radii - dist

                if overlap > 0:  # A clash
                    sum_squared_overlaps += overlap**2

    if num_distances_considered == 0:
        return 0.0  # No atoms close enough to be considered

    tcs = np.sqrt(sum_squared_overlaps / num_distances_considered)
    return tcs


if __name__ == "__main__":
    # Example Usage (requires dummy PDB and SDF files)
    # Create dummy protein.pdb
    dummy_pdb_content = """
ATOM      1  N   ALA A   1       8.000  -2.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       8.000  -1.000   0.500  1.00  0.00           C
ATOM      3  C   ALA A   1       7.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       7.000   1.000   0.500  1.00  0.00           O
ATOM      5  CB  ALA A   1       9.000  -1.000   1.500  1.00  0.00           C
"""
    with open("dummy_protein.pdb", "w") as f:
        f.write(dummy_pdb_content)

    # Create dummy ligand.sdf with a clash
    dummy_sdf_content = """MyLigand
  -ISIS-  05132513102D

  2  1  0  0  0  0  0  0  0  0999 V2000
    7.5000   -0.5000    0.2000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
$$$$
"""
    # Atom 1 (C) @ 7.5, -0.5, 0.2 will clash with C@7,0,0 and CA@8,-1,0.5
    with open("dummy_ligand.sdf", "w") as f:
        f.write(dummy_sdf_content)

    protein_atoms_list = parse_pdb_atoms(Path("dummy_protein.pdb"))
    ligand_atoms_list = parse_sdf_atoms(Path("dummy_ligand.sdf"))

    print(f"Parsed {len(protein_atoms_list)} protein atoms.")
    print(f"Parsed {len(ligand_atoms_list)} ligand atoms.")

    if protein_atoms_list and ligand_atoms_list:
        tcs_score = calculate_tcs(protein_atoms_list, ligand_atoms_list)
        print(f"Calculated TCS: {tcs_score:.4f}")

    # Clean up dummy files
    Path("dummy_protein.pdb").unlink()
    Path("dummy_ligand.sdf").unlink()
