"""Module for fetching protein structures from the AlphaFold Database."""

from foldfusion.tools.tool import Tool
from pathlib import Path
import requests
import logging

logger = logging.getLogger(__name__)


class AlphaFoldFetcher(Tool):
    """Fetches protein structures from the AlphaFold Database.

    Attributes:
        uniprot_id (str): The UniProt ID of the protein to fetch.
        output_dir (Path): The directory where the fetched PDB file will be saved.
    """

    def __init__(self, config: dict):
        """Initializes the AlphaFoldFetcher.

        Args:
            config (dict): The configuration dictionary.
        """
        super().__init__(config)
        self.tool_name = "alphafoldfetcher"
        run_settings = self.config.get("run_settings", {})
        if not run_settings or "uniprot_id" not in run_settings:
            logger.error(
                "UniProt ID not found in run_settings within the configuration."
            )
            raise ValueError("UniProt ID must be specified in run_settings.")
        self.uniprot_id = run_settings["uniprot_id"]
        self.output_dir = self.output_dir / "alphafold"
        logger.debug(f"AlphaFoldFetcher initialized for UniProt ID: {self.uniprot_id}")

    def get_alphafold_model(self) -> Path:
        """Downloads the AlphaFold PDB model for the specified UniProt ID.

        Tries fetching model versions v4 and then v3.

        Returns:
            Path: The path to the downloaded PDB file.

        Raises:
            FileNotFoundError: If the PDB file cannot be fetched from AlphaFold DB.
            requests.exceptions.RequestException: For network or HTTP errors.
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Ensured output directory exists: {self.output_dir}")

        versions = ["v4", "v3"]  # Prioritize v4
        pdb_path = None
        last_error = None

        for version in versions:
            model_version_suffix = f"model_{version}"
            url = (
                f"https://alphafold.ebi.ac.uk/files/"
                f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            )
            output_file_name = f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            output_file_path = self.output_dir / output_file_name

            logger.info(
                f"Attempting to download AlphaFold model version {version} from {url}"
            )
            try:
                response = requests.get(url, timeout=60)
                response.raise_for_status()  # Raises HTTPError for bad responses (4XX or 5XX)

                # No need to check response.status_code == 200 due to raise_for_status()
                output_file_path.write_text(response.text, encoding="utf-8")
                logger.info(
                    f"Successfully downloaded and saved model version {version} to {output_file_path}"
                )
                pdb_path = output_file_path
                last_error = None  # Reset last error on success
                break  # Exit loop if download is successful
            except requests.exceptions.HTTPError as e:
                last_error = e
                if e.response.status_code == 404:
                    logger.warning(
                        f"Model version {version} not found for UniProt ID "
                        f"{self.uniprot_id} (404 Error). Trying next version."
                    )
                else:
                    logger.error(
                        f"HTTP error occurred while fetching model version {version} "
                        f"for {self.uniprot_id}. Status code: {e.response.status_code}. Error: {e}"
                    )
            except requests.exceptions.RequestException as e:
                last_error = e
                logger.error(
                    f"Failed to fetch structure version {version} for "
                    f"{self.uniprot_id} due to a network or request issue. "
                    f"Error: {e}"
                )

        if pdb_path:
            resolved_path = pdb_path.resolve()
            logger.info(
                f"Returning resolved path for downloaded model: {resolved_path}"
            )
            return resolved_path
        else:
            error_message = (
                f"Failed to fetch any AlphaFold structure version for {self.uniprot_id} "
                f"after trying versions: {', '.join(versions)}."
            )
            logger.error(error_message)
            if last_error:
                logger.error(f"Last error encountered: {last_error}")
                # Raise a more specific error that also includes the original error
                raise FileNotFoundError(
                    f"{error_message} Last error: {last_error.__class__.__name__} - {last_error}"
                ) from last_error
            else:
                # This case should ideally not be reached if versions list is not empty
                # and an attempt was made.
                raise FileNotFoundError(error_message)
