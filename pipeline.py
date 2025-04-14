# pipeline.py
import argparse
import logging
import subprocess
import sys
from pathlib import Path
import os
import requests  # Added for AlphafoldFetcher dependency

# --- AlphafoldFetcher Class (Integrated and slightly adapted) ---
# Get a logger for this module (or use the main script's logger)
fetcher_logger = logging.getLogger(__name__ + ".AlphafoldFetcher")


class AlphafoldFetcher:
    # Removed config dict, pass output_dir directly
    def __init__(self, uniprot_id: str, base_output_dir: Path):
        self.uniprot_id = uniprot_id
        self.logger = fetcher_logger
        # Create a specific subdirectory for AlphaFold output within the base dir
        self.output_dir = base_output_dir / "alphafold"
        self.logger.debug(f"AlphaFold output directory set to {self.output_dir}")
        self.logger.debug(
            f"AlphafoldFetcher initialized with uniprot_id: {self.uniprot_id}"
        )

    def get_alphafold_model(self) -> Path:
        self.logger.info(f"Fetching Alphafold structure for {self.uniprot_id}")

        # Try fetching v4 first, then fallback to v3 if v4 not found
        versions_to_try = ["v4", "v3"]
        pdb_path = None

        for version in versions_to_try:
            model_version_suffix = f"model_{version}"
            url = (
                "https://alphafold.ebi.ac.uk/files/"
                + f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            )
            output_filename = f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            output_path = self.output_dir / output_filename

            self.logger.info(f"Attempting to download from: {url}")

            # Ensure the specific output directory exists
            self.output_dir.mkdir(parents=True, exist_ok=True)

            try:
                response = requests.get(url, timeout=60)  # Added timeout
                response.raise_for_status()  # Raises HTTPError for bad responses

                if response.status_code == 200:
                    output_path.write_text(response.text)
                    self.logger.info(f"Structure ({version}) saved to {output_path}")
                    pdb_path = output_path
                    break  # Success, exit loop
                else:
                    # Should be caught by raise_for_status, but good practice
                    self.logger.warning(
                        f"Received status code {response.status_code} for {url} "
                        + f"(version {version})"
                    )

            except requests.exceptions.RequestException as e:
                self.logger.warning(
                    f"Failed to fetch structure version {version} for "
                    + f"{self.uniprot_id}. "
                    f"Error: {e}"
                )
                # If the error is a 404, it just means this version doesn't exist
                if (
                    isinstance(e, requests.exceptions.HTTPError)
                    and e.response.status_code == 404
                ):
                    self.logger.info(
                        f"Model version {version} not found for {self.uniprot_id}."
                    )
                # Otherwise, log the warning and continue to the next version

        if pdb_path:
            return pdb_path
        else:
            self.logger.error(
                f"Failed to fetch any structure version for {self.uniprot_id} "
                + f"after trying {versions_to_try}."
            )
            raise FileNotFoundError(
                f"Could not download AlphaFold structure for {self.uniprot_id}"
            )


# --- Main Pipeline Logic ---


def run_command(command, step_name):
    """Helper function to run a shell command and handle errors."""
    logger = logging.getLogger(__name__)
    logger.info(f"Running {step_name}...")
    logger.debug(f"Executing command: {' '.join(map(str, command))}")
    try:
        result = subprocess.run(
            command,
            check=True,  # Raise CalledProcessError on non-zero exit code
            capture_output=True,  # Capture stdout and stderr
            text=True,  # Decode output as text
            encoding="utf-8",  # Explicitly set encoding
        )
        logger.info(f"{step_name} completed successfully.")
        logger.debug(f"{step_name} stdout:\n{result.stdout}")
        if result.stderr:
            logger.debug(f"{step_name} stderr:\n{result.stderr}")
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"{step_name} failed with exit code {e.returncode}.")
        logger.error(f"Command: {' '.join(map(str, e.cmd))}")
        logger.error(f"Stdout:\n{e.stdout}")
        logger.error(f"Stderr:\n{e.stderr}")
        raise  # Re-raise the exception to stop the pipeline
    except FileNotFoundError:
        logger.error(f"Error: Executable not found for {step_name} at {command[0]}")
        logger.error(
            "Please ensure the path to the executable is correct and it has execute "
            + "permissions."
        )
        raise
    except Exception as e:
        logger.error(f"An unexpected error occurred during {step_name}: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Pipeline to fetch AlphaFold structure, run DoGSite3, and SIENA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--uniprot_id", required=True, help="UniProt ID of the target protein."
    )
    parser.add_argument(
        "--dogsite_executable",
        required=True,
        type=Path,
        help="Path to the dogsite3 executable.",
    )
    parser.add_argument(
        "--siena_executable",
        required=True,
        type=Path,
        help="Path to the siena executable.",
    )
    parser.add_argument(
        "--siena_db",
        required=True,
        type=Path,
        help="Path to the pre-generated SIENA database file.",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=Path("pipeline_output"),
        help="Base directory for all output files.",
    )
    parser.add_argument(
        "--dogsite_options",
        default="--writeSiteResiduesEDF",
        help="Additional options for DoGSite3 "
        + "(ensure --writeSiteResiduesEDF is included if overriding).",
    )
    parser.add_argument(
        "--siena_options", default="", help="Additional options for SIENA."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable debug logging."
    )

    args = parser.parse_args()

    # --- Setup Logging ---
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],  # Ensure logs go to stdout
    )
    logger = logging.getLogger(__name__)
    logger.info("Starting pipeline...")
    logger.info(f"Run arguments: {args}")

    # Create base output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Using base output directory: {args.output_dir.resolve()}")

    # --- Step 1: Fetch AlphaFold Structure ---
    try:
        fetcher = AlphafoldFetcher(args.uniprot_id, args.output_dir)
        alphafold_pdb_path = fetcher.get_alphafold_model()
        logger.info(f"Successfully downloaded AlphaFold PDB: {alphafold_pdb_path}")
    except Exception as e:
        logger.error(f"Failed Step 1 (AlphaFold Fetch): {e}")
        sys.exit(1)  # Exit if fetching fails

    # --- Step 2: Run DoGSite3 ---
    dogsite_output_dir = args.output_dir / "dogsite"
    dogsite_output_dir.mkdir(parents=True, exist_ok=True)
    # Define a base name for dogsite outputs based on the PDB file name stem
    dogsite_output_base = dogsite_output_dir / alphafold_pdb_path.stem

    # Construct DoGSite3 command
    dogsite_cmd = [
        args.dogsite_executable,
        "-p",
        alphafold_pdb_path,
        "-o",
        dogsite_output_base,
    ]
    # Add user-specified options, split by space
    # Ensure --writeSiteResiduesEDF is present
    dogsite_options_list = args.dogsite_options.split()
    if "--writeSiteResiduesEDF" not in dogsite_options_list:
        logger.warning("Adding missing '--writeSiteResiduesEDF' to DoGSite3 options.")
        dogsite_options_list.append("--writeSiteResiduesEDF")

    dogsite_cmd.extend(dogsite_options_list)

    try:
        run_command(dogsite_cmd, "DoGSite3")
        # Assume we want the first pocket's EDF file
        dogsite_edf_path = (
            dogsite_output_base.with_suffix(".edf").parent
            / f"{dogsite_output_base.name}_P_1_res.edf"
        )

        if not dogsite_edf_path.exists():
            logger.error(
                f"DoGSite3 ran, but expected EDF file not found at: {dogsite_edf_path}"
            )
            logger.error(
                "Check DoGSite3 output/logs for errors or if pockets were detected."
            )
            sys.exit(1)
        logger.info(f"DoGSite3 generated EDF file: {dogsite_edf_path}")

    except Exception as e:
        logger.error(f"Failed Step 2 (DoGSite3): {e}")
        # Potentially add reminder about license activation here
        logger.warning(
            "Ensure the DoGSite3 license is activated before running the pipeline."
        )
        sys.exit(1)

    # --- Step 2.5: Edit DoGSite EDF file to reference AlphaFold PDB ---
    logger.info("Updating DoGSite EDF file to reference the AlphaFold PDB file...")
    try:
        # Read the EDF file content
        edf_content = dogsite_edf_path.read_text()

        # Replace the "REFERENCE <NO-FILE>" line with the correct reference
        updated_content = edf_content.replace(
            "REFERENCE <NO-FILE>", f"REFERENCE {alphafold_pdb_path}"
        )

        # Write the updated content back to the file
        dogsite_edf_path.write_text(updated_content)
        logger.info(
            f"Successfully updated EDF file with reference to: {alphafold_pdb_path}"
        )
    except Exception as e:
        logger.error(
            f"Failed to update DoGSite EDF file with AlphaFold PDB reference: {e}"
        )
        # This isn't a fatal error, so we continue with the pipeline
        logger.warning("Continuing pipeline with unmodified EDF file")

    # --- Step 3: Run SIENA ---
    siena_output_dir = args.output_dir / "siena"
    siena_output_dir.mkdir(parents=True, exist_ok=True)
    # Define output directory for SIENA results specifically
    siena_results_output = siena_output_dir / f"{args.uniprot_id}_siena_results"
    siena_results_output.mkdir(parents=True, exist_ok=True)
    logger.info(f"SIENA results will be saved to: {siena_results_output}")

    # Construct SIENA command
    siena_cmd = [
        args.siena_executable,
        "-p",
        alphafold_pdb_path,
        "-e",
        dogsite_edf_path,  # Use the EDF file from DoGSite3
        "-b",
        args.siena_db,
        "-o",
        siena_results_output,
    ]

    logger.info(f"SIENA command: {' '.join(map(str, siena_cmd))}")
    # Check if the SIENA database file exists
    if not args.siena_db.exists():
        logger.error(f"SIENA database file not found: {args.siena_db}")
        sys.exit(1)
    # Check if the SIENA executable exists
    if not args.siena_executable.exists():
        logger.error(f"SIENA executable not found: {args.siena_executable}")
        sys.exit(1)
    # Check if the AlphaFold PDB file exists
    if not alphafold_pdb_path.exists():
        logger.error(f"AlphaFold PDB file not found: {alphafold_pdb_path}")
        sys.exit(1)
    # Check if the DoGSite EDF file exists
    if not dogsite_edf_path.exists():
        logger.error(f"DoGSite EDF file not found: {dogsite_edf_path}")
        sys.exit(1)
    # Check if the output directory exists
    if not siena_output_dir.exists():
        logger.error(f"Output directory for SIENA not found: {siena_output_dir}")
        sys.exit(1)

    # Add user-specified options, split by space
    if args.siena_options:
        siena_cmd.extend(args.siena_options.split())

    try:
        run_command(siena_cmd, "SIENA")
        logger.info(f"SIENA results should be in: {siena_results_output}")
    except Exception as e:
        logger.error(f"Failed Step 3 (SIENA): {e}")
        # Potentially add reminder about license activation here
        logger.warning(
            "Ensure the SIENA license is activated and the database path is correct."
        )
        sys.exit(1)

    logger.info("Pipeline finished successfully!")


if __name__ == "__main__":
    main()
