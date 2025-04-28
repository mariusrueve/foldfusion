import logging
from pathlib import Path

# Import necessary components from other modules within the package
from .utils import run_command, get_ligands_from_siena_pdb
from .fetchers import AlphafoldFetcher

# Get a logger specific to this steps module
logger = logging.getLogger(__name__)


def fetch_alphafold_step(uniprot_id: str, output_dir: Path) -> Path:
    """
    Pipeline step to fetch the AlphaFold structure.

    Args:
        uniprot_id (str): The UniProt ID.
        output_dir (Path): The base output directory.

    Returns:
        Path: The path to the downloaded AlphaFold PDB file.

    Raises:
        Exception: If fetching fails (propagated from AlphafoldFetcher).
    """
    logger.info("--- Fetch AlphaFold Structure ---")
    try:
        # Initialize the fetcher with the base output directory
        fetcher = AlphafoldFetcher(uniprot_id, output_dir)
        # Call the method to get the model
        alphafold_pdb_path = fetcher.get_alphafold_model()
        logger.info(f"Successfully downloaded AlphaFold PDB: {alphafold_pdb_path}")
        return alphafold_pdb_path
    except Exception as e:
        # Log the error at this level and re-raise to be caught by main
        logger.error(f"AlphaFold fetch step failed: {e}")
        raise


def run_dogsite_step(
    pdb_path: Path, dogsite_executable: Path, dogsite_options: str, output_dir: Path
) -> Path:
    """
    Pipeline step to run DoGSite3.

    Args:
        pdb_path (Path): Path to the input PDB file (e.g., from AlphaFold).
        dogsite_executable (Path): Path to the DoGSite3 executable.
        dogsite_options (str): Space-separated string of additional DoGSite3 options.
        output_dir (Path): The base output directory for the entire pipeline.

    Returns:
        Path: The path to the generated DoGSite3 EDF file (for the first pocket).

    Raises:
        FileNotFoundError: If the expected EDF file is not found after execution.
        Exception: If DoGSite3 command fails (propagated from run_command).
    """
    logger.info("--- Run DoGSite3 ---")
    # Create a dedicated subdirectory for DoGSite outputs
    dogsite_output_dir = output_dir / "dogsite"
    dogsite_output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"DoGSite3 output directory: {dogsite_output_dir}")

    # Define a base name for dogsite outputs based on the PDB file name stem
    # Outputs will be placed in the dogsite_output_dir
    dogsite_output_base = dogsite_output_dir / pdb_path.stem

    # Construct DoGSite3 command list
    dogsite_cmd = [
        dogsite_executable,
        "-p",
        pdb_path,  # Input PDB file path
        "-o",
        dogsite_output_base,  # Base path for output files
    ]

    # Process and add user-specified options, ensuring EDF generation flag
    dogsite_options_list = dogsite_options.split()
    if "--writeSiteResiduesEDF" not in dogsite_options_list:
        logger.warning("Adding missing '--writeSiteResiduesEDF' to DoGSite3 options.")
        dogsite_options_list.append("--writeSiteResiduesEDF")
    # Extend the command list with the processed options
    dogsite_cmd.extend(dogsite_options_list)

    try:
        # Execute the DoGSite3 command using the utility function
        run_command(dogsite_cmd, "DoGSite3")

        # Define the expected path for the EDF file (pocket 1 residue list)
        # Based on DoGSite naming convention: <output_base>_P_1_res.edf
        expected_edf_filename = f"{dogsite_output_base.name}_P_1_res.edf"
        dogsite_edf_path = dogsite_output_dir / expected_edf_filename

        # Verify that the expected EDF file was created
        if not dogsite_edf_path.is_file():
            logger.error(
                f"DoGSite3 ran, but expected EDF file not found at: {dogsite_edf_path}"
            )
            logger.error(
                "Check DoGSite3 output/logs in the 'dogsite' subdirectory for errors "
                "or if pockets were detected."
            )
            # Attempt to list files in the output directory for debugging context
            try:
                files_in_dir = list(dogsite_output_dir.glob("*"))
                logger.error(
                    f"Files found in {dogsite_output_dir}: "
                    + f"{[f.name for f in files_in_dir]}"
                )
            except Exception as list_e:
                logger.error(f"Could not list files in {dogsite_output_dir}: {list_e}")

            # Raise FileNotFoundError as the crucial output is missing
            raise FileNotFoundError(
                f"Expected DoGSite3 EDF file not found: {dogsite_edf_path}"
            )

        logger.info(f"DoGSite3 generated EDF file: {dogsite_edf_path}")
        return dogsite_edf_path

    except Exception as e:
        # Log error specific to this step and re-raise
        logger.error(f"DoGSite3 step failed: {e}")
        # Add reminder about potential license issues
        logger.warning("Ensure the DoGSite3 license is activated if applicable.")
        raise


def modify_edf_step(edf_path: Path, pdb_path: Path) -> Path:
    """
    Pipeline step to modify the DoGSite EDF file to reference the correct PDB.

    Args:
        edf_path (Path): Path to the DoGSite EDF file.
        pdb_path (Path): Path to the reference PDB file (e.g., AlphaFold PDB).

    Returns:
        Path: The path to the (potentially modified) EDF file.
              Returns the original path even if modification fails or is skipped.
    """
    logger.info("--- Modify DoGSite EDF File ---")
    # Check if the input EDF file exists before attempting modification
    if not edf_path.is_file():
        logger.error(f"Cannot modify EDF file: File not found at {edf_path}")
        # Depending on requirements, you might raise FileNotFoundError here
        # For now, we log the error and return the non-existent path,
        # likely causing SIENA to fail later with a clearer error.
        return edf_path

    logger.info(f"Attempting to update {edf_path.name} to reference {pdb_path.name}...")
    try:
        # Read the original EDF file content
        edf_content = edf_path.read_text(encoding="utf-8")

        # Define the target line and the replacement line
        target_line = "REFERENCE <NO-FILE>"
        # Use the absolute path of the PDB file for the reference
        replacement_line = f"REFERENCE {pdb_path.resolve()}"

        # Check if the target line exists in the content
        if target_line in edf_content:
            # Perform the replacement
            updated_content = edf_content.replace(target_line, replacement_line)

            # Write the updated content back to the same file
            edf_path.write_text(updated_content, encoding="utf-8")
            logger.info(
                f"Successfully updated EDF file '{edf_path.name}' with reference to: "
                + f"{pdb_path.resolve()}"
            )
        else:
            # Log a warning if the target line wasn't found
            logger.warning(
                f"Could not find '{target_line}' in {edf_path.name}. "
                "File may already be modified or in an unexpected format. Skipping"
                + "modification."
            )

    except Exception as e:
        # Log any errors during file reading/writing
        logger.error(
            f"Failed to read or write EDF file '{edf_path.name}' during modification "
            + f"step: {e}",
            exc_info=True,  # Include traceback for unexpected errors
        )
        logger.warning(
            "Continuing pipeline with potentially unmodified EDF file due to error."
        )

    # Always return the path to the EDF file, whether modified or not
    return edf_path


def run_siena_step(
    pdb_path: Path,
    edf_path: Path,
    siena_executable: Path,
    siena_db: Path,
    siena_options: str,
    output_dir: Path,
    uniprot_id: str,
) -> Path:
    """
    Pipeline step to run SIENA.

    Args:
        pdb_path (Path): Path to the input PDB file.
        edf_path (Path): Path to the (potentially modified) DoGSite EDF file.
        siena_executable (Path): Path to the SIENA executable.
        siena_db (Path): Path to the SIENA database file.
        siena_options (str): Space-separated string of additional SIENA options.
        output_dir (Path): The base output directory for the entire pipeline.
        uniprot_id (str): UniProt ID used for naming the SIENA output subdirectory.

    Returns:
        Path: The path to the SIENA results directory.

    Raises:
        FileNotFoundError: If essential input files (PDB, EDF, DB, executable)
        are missing.
        Exception: If SIENA command fails (propagated from run_command).
    """
    logger.info("--- Run SIENA ---")
    # Define the base directory for SIENA outputs
    siena_output_base_dir = output_dir / "siena"
    siena_output_base_dir.mkdir(parents=True, exist_ok=True)

    # Define the specific output directory for this run's results
    siena_results_output_dir = siena_output_base_dir / f"{uniprot_id}_siena_results"
    # SIENA typically expects the output directory to exist
    siena_results_output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"SIENA results directory: {siena_results_output_dir}")

    # These checks are crucial before attempting to run the command
    if not pdb_path.is_file():
        logger.error(f"SIENA input PDB file not found: {pdb_path}")
        raise FileNotFoundError(f"SIENA input PDB file not found: {pdb_path}")
    if not edf_path.is_file():
        # This check is important after the modify step might have logged an error
        logger.error(f"SIENA input EDF file not found: {edf_path}")
        raise FileNotFoundError(f"SIENA input EDF file not found: {edf_path}")
    # Executable and DB existence are checked in main.py's setup,
    # but double-checking doesn't hurt
    if not siena_executable.is_file():
        logger.error(f"SIENA executable not found: {siena_executable}")
        raise FileNotFoundError(f"SIENA executable not found: {siena_executable}")
    if not siena_db.is_file():
        logger.error(f"SIENA database file not found: {siena_db}")
        raise FileNotFoundError(f"SIENA database file not found: {siena_db}")

    # Construct SIENA command list
    siena_cmd = [
        siena_executable,
        "-p",
        pdb_path,  # Input PDB file
        "-e",
        edf_path,  # Input EDF file (potentially modified)
        "-b",
        siena_db,  # SIENA database path
        "-o",
        siena_results_output_dir,  # Output directory path
    ]

    # Add any user-specified options
    if siena_options:
        siena_cmd.extend(siena_options.split())

    try:
        # Execute the SIENA command
        run_command(siena_cmd, "SIENA")
        logger.info(
            "SIENA execution completed. Results should be in:"
            + f"{siena_results_output_dir}"
        )
        return siena_results_output_dir
    except Exception as e:
        # Log error specific to this step and re-raise
        logger.error(f"SIENA step failed: {e}")
        # Add reminder about potential license/DB issues
        logger.warning(
            "Ensure the SIENA license is activated (if applicable) and the database "
            + "path is correct."
        )
        raise


def run_ligand_extractor_step(
    ligand_extractor_executable: Path, siena_pdb_path: Path, output_dir: Path
) -> Path:
    """
    Runs the ligand extractor tool for each ligand found in a SIENA PDB file.

    Args:
        ligand_extractor_executable (Path): Path to the ligand extractor executable.
        siena_pdb_path (Path): Path to the SIENA PDB file containing ligands.
        output_dir (Path): The specific directory where ligand files for this PDB
                           should be saved (e.g., .../ligand_extractor/pdb_name/).

    Returns:
        Path: The output directory containing the extracted ligand files.
    """
    logger.info("--- Run Ligand Extractor ---")
    logger.info(
        f"Ligand extractor output directory: {output_dir}"
    )  # Use the provided output_dir directly

    # Ensure the final output directory exists (it should have been created in main.py)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get the list of ligands from the SIENA PDB file
    # The format returned is like 'LIG_A_123'
    ligand_identifiers = get_ligands_from_siena_pdb(siena_pdb_path)
    logger.info(
        f"Found {len(ligand_identifiers)} ligands in {siena_pdb_path.name}: {ligand_identifiers}"
    )

    if not ligand_identifiers:
        logger.warning(
            f"No ligands found in {siena_pdb_path.name}, skipping extraction."
        )
        return output_dir  # Return the directory even if empty

    extracted_ligand_files = []  # Keep track of successfully extracted files

    for ligand_id_str in ligand_identifiers:
        # Construct the command to run the ligand extractor
        # The ligand_id_str (e.g., 'LIG_A_123') is passed directly to the tool
        le_cmd = [
            str(ligand_extractor_executable),
            "-c",
            str(siena_pdb_path),  # Ensure path is string
            "-l",
            ligand_id_str,  # Use the identifier string directly
            "-o",
            str(output_dir),  # Use the provided output_dir directly for output files
        ]
        logger.debug(f"Running Ligand Extractor command: {' '.join(le_cmd)}")
        try:
            run_command(le_cmd, f"Ligand Extractor ({ligand_id_str})")
            # Construct expected output filename based on ligand_id_str
            expected_sdf_filename = f"{ligand_id_str}.sdf"
            expected_sdf_path = output_dir / expected_sdf_filename
            if expected_sdf_path.is_file():
                logger.info(
                    f"Successfully extracted ligand {ligand_id_str} to {expected_sdf_path}"
                )
                extracted_ligand_files.append(expected_sdf_path)
            else:
                logger.warning(
                    f"Ligand Extractor ran for {ligand_id_str}, but output file {expected_sdf_path} not found."
                )

        except Exception as e:
            logger.error(
                f"Failed to extract ligand {ligand_id_str} from {siena_pdb_path.name}: {e}"
            )
            # Decide if you want to continue with other ligands or raise the exception
            # continue

    # Return the output directory containing the ligand extraction results
    # The caller (main.py) will need to find the .sdf files within this directory
    return output_dir


def run_jamda_scorer_step(
    protein_pdb_path: Path,
    ligand_sdf_path: Path,
    output_sdf_path: Path,
    jamda_scorer_executable: Path,
    jamda_scorer_options: str,
) -> Path:
    """
    Pipeline step to run JamdaScorer.

    Args:
        protein_pdb_path (Path): Path to the input protein PDB file (e.g., AlphaFold).
        ligand_sdf_path (Path): Path to the input ligand SDF file (from ligand extractor).
        output_sdf_path (Path): Path where the optimized ligand SDF file will be saved.
        jamda_scorer_executable (Path): Path to the JamdaScorer executable.
        jamda_scorer_options (str): Space-separated string of additional JamdaScorer options.

    Returns:
        Path: The path to the generated optimized SDF file.

    Raises:
        FileNotFoundError: If input files are not found or output is not created.
        Exception: If JamdaScorer command fails (propagated from run_command).
    """
    logger.info("--- Run JamdaScorer ---")
    logger.info(f"Input Protein: {protein_pdb_path.name}")
    logger.info(f"Input Ligand: {ligand_sdf_path.name}")
    logger.info(f"Output SDF: {output_sdf_path}")

    # Ensure input files exist
    if not protein_pdb_path.is_file():
        logger.error(f"JamdaScorer input PDB file not found: {protein_pdb_path}")
        raise FileNotFoundError(
            f"JamdaScorer input PDB file not found: {protein_pdb_path}"
        )
    if not ligand_sdf_path.is_file():
        logger.error(f"JamdaScorer input SDF file not found: {ligand_sdf_path}")
        raise FileNotFoundError(
            f"JamdaScorer input SDF file not found: {ligand_sdf_path}"
        )

    # Ensure output directory exists
    output_sdf_path.parent.mkdir(parents=True, exist_ok=True)

    # Construct JamdaScorer command list
    jamda_cmd = [
        str(jamda_scorer_executable),
        "-i",
        str(protein_pdb_path),
        "-l",
        str(ligand_sdf_path),
        "-o",
        str(output_sdf_path),
    ]

    # Add any user-specified options
    if jamda_scorer_options:
        jamda_cmd.extend(jamda_scorer_options.split())

    try:
        # Execute the JamdaScorer command
        run_command(jamda_cmd, f"JamdaScorer ({ligand_sdf_path.stem})")

        # Verify that the expected output file was created
        if not output_sdf_path.is_file():
            logger.error(
                f"JamdaScorer ran, but expected output SDF file not found at: {output_sdf_path}"
            )
            raise FileNotFoundError(
                f"Expected JamdaScorer output SDF file not found: {output_sdf_path}"
            )

        logger.info(f"JamdaScorer generated optimized SDF file: {output_sdf_path}")
        return output_sdf_path

    except Exception as e:
        # Log error specific to this step and re-raise
        logger.error(f"JamdaScorer step failed for ligand {ligand_sdf_path.name}: {e}")
        raise
