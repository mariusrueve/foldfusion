#!/usr/bin/env python3
import argparse
import logging
import subprocess
import sys
from pathlib import Path

try:
    import tomllib
except ImportError:
    # Fallback for Python < 3.11
    try:
        import toml as tomllib
    except ImportError:
        print("Error: 'tomllib' (Python 3.11+) or 'toml' package required.")
        print("Install 'toml' with: pip install toml")
        sys.exit(1)


import requests

# Import steps from the pipeline package
from pipeline import steps

# Configure root logger early. Can be moved to utils if preferred.
logging.basicConfig(
    level=logging.INFO,  # Default level, overridden by args
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
# Get a logger for the main script
logger = logging.getLogger(__name__)


def load_config(config_path: Path) -> dict:
    """Loads configuration from a TOML file."""
    logger.info(f"Loading configuration from: {config_path}")
    if not config_path.is_file():
        logger.error(f"Configuration file not found: {config_path}")
        sys.exit(1)
    try:
        with open(config_path, "rb") as f:
            config = tomllib.load(f)
        logger.debug(f"Configuration loaded: {config}")
        # Basic validation (can be expanded)
        if "executables" not in config or "options" not in config:
            raise ValueError(
                "Config file must contain 'executables' and 'options' sections."
            )
        if (
            "dogsite" not in config["executables"]
            or "siena" not in config["executables"]
        ):
            raise ValueError(
                "Config 'executables' section must contain 'dogsite' and 'siena' keys."
            )
        if "dogsite" not in config["options"] or "siena" not in config["options"]:
            raise ValueError(
                "Config 'options' section must contain 'dogsite' and 'siena' keys."
            )
        return config
    except tomllib.TOMLDecodeError as e:
        logger.error(f"Error decoding TOML file {config_path}: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Failed to load or validate config file {config_path}: {e}")
        sys.exit(1)


def setup_pipeline():
    """Parses command-line arguments, loads config, and sets up logging."""
    parser = argparse.ArgumentParser(
        description="Modular pipeline: AlphaFold -> DoGSite3 -> SIENA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Arguments specific to the run
    parser.add_argument(
        "--uniprot_id", required=True, help="UniProt ID of the target protein."
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
    # Configuration file path
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.toml"),
        help="Path to the configuration TOML file.",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable debug logging."
    )

    args = parser.parse_args()

    config = load_config(args.config)
    # Add config values to args namespace for convenience, converting paths
    args.dogsite_executable = Path(config["executables"]["dogsite"])
    args.siena_executable = Path(config["executables"]["siena"])
    args.dogsite_options = config["options"]["dogsite"]
    args.siena_options = config["options"]["siena"]

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.getLogger().setLevel(log_level)  # Set root logger level
    # Set level for specific external loggers if needed
    logging.getLogger("requests").setLevel(
        logging.WARNING if not args.verbose else logging.DEBUG
    )
    logging.getLogger("urllib3").setLevel(
        logging.WARNING if not args.verbose else logging.DEBUG
    )

    logger.info("Starting pipeline...")
    logger.info(f"Run arguments: {args}")  # Now includes config values

    # Create base output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Using base output directory: {args.output_dir.resolve()}")

    # Validate executable and DB paths (essential inputs)
    # Executable paths now come from config via args
    if not args.dogsite_executable.is_file():
        logger.error(
            f"DoGSite3 executable not found or not a file (from config): {args.dogsite_executable}"
        )
        sys.exit(1)
    if not args.siena_executable.is_file():
        logger.error(
            f"SIENA executable not found or not a file (from config): {args.siena_executable}"
        )
        sys.exit(1)
    if not args.siena_db.is_file():
        logger.error(f"SIENA database not found or not a file: {args.siena_db}")
        sys.exit(1)

    # Ensure dogsite options include EDF generation if needed by the step
    if "--writeSiteResiduesEDF" not in args.dogsite_options.split():
        logger.warning(
            "Configured DoGSite3 options missing '--writeSiteResiduesEDF'. "
            "The run_dogsite_step will add it if necessary, but consider "
            "updating config.toml."
        )

    return args


def main():
    """Main function to orchestrate the pipeline steps."""
    pipeline_status = 0  # 0 for success, 1 for failure
    args = None  # Initialize args to None
    try:
        # Setup now includes loading config and adding to args
        args = setup_pipeline()

        # Step 1: Fetch AlphaFold
        # Pass necessary args to the step function
        alphafold_pdb_path = steps.fetch_alphafold_step(
            uniprot_id=args.uniprot_id, output_dir=args.output_dir
        )

        # Executable path and options now come from args (loaded from config)
        dogsite_edf_path = steps.run_dogsite_step(
            pdb_path=alphafold_pdb_path,
            dogsite_executable=args.dogsite_executable,
            dogsite_options=args.dogsite_options,
            output_dir=args.output_dir,
        )

        modified_edf_path = steps.modify_edf_step(
            edf_path=dogsite_edf_path, pdb_path=alphafold_pdb_path
        )

        # Executable path and options now come from args (loaded from config)
        siena_results_dir = steps.run_siena_step(
            pdb_path=alphafold_pdb_path,
            edf_path=modified_edf_path,
            siena_executable=args.siena_executable,
            siena_db=args.siena_db,  # Still from command line
            siena_options=args.siena_options,
            output_dir=args.output_dir,
            uniprot_id=args.uniprot_id,
        )

        logger.info(f"SIENA results saved in: {siena_results_dir.resolve()}")

        logger.info(
            f"Pipeline finished successfully! Results in: {args.output_dir.resolve()}"
        )

    # Catch specific exceptions anticipated from steps or utils
    except (
        FileNotFoundError,
        subprocess.CalledProcessError,
        requests.exceptions.RequestException,
    ) as e:
        logger.error(f"Pipeline execution failed due to a known error: {e}")
        pipeline_status = 1
    # Catch unexpected exceptions
    except Exception as e:
        logger.exception(f"An unexpected error occurred during pipeline execution: {e}")
        pipeline_status = 1
    finally:
        if pipeline_status == 0:
            logger.info("Pipeline completed.")
        else:
            logger.error("Pipeline failed.")
        # Ensure args exists before trying to access output_dir in case setup failed
        # early
        if args:
            logger.info(f"Output directory: {args.output_dir.resolve()}")
        sys.exit(pipeline_status)


if __name__ == "__main__":
    main()
