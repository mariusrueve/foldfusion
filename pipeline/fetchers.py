import logging
import requests
from pathlib import Path

# Get a logger specific to this fetchers module
logger = logging.getLogger(__name__)


class AlphafoldFetcher:
    """Fetches AlphaFold PDB structures from the EBI database."""

    def __init__(self, uniprot_id: str, base_output_dir: Path):
        """
        Initializes the fetcher.

        Args:
            uniprot_id (str): The UniProt ID of the protein.
            base_output_dir (Path): The base directory where AlphaFold files
                                     will be saved (within an 'alphafold' subdir).
        """
        self.uniprot_id = uniprot_id
        # Use the logger obtained for this module
        self.logger = logger
        self.output_dir = base_output_dir / "alphafold"
        self.logger.debug(f"AlphaFold output directory set to {self.output_dir}")
        self.logger.debug(
            f"AlphafoldFetcher initialized with uniprot_id: {self.uniprot_id}"
        )

    def get_alphafold_model(self) -> Path:
        """
        Downloads the AlphaFold PDB file (tries v4 then v3).

        Returns:
            Path: The path to the downloaded PDB file.

        Raises:
            FileNotFoundError: If no structure could be downloaded after trying all
            versions.
            requests.exceptions.RequestException: If a non-404 download error occurs.
        """
        self.logger.info(f"Fetching Alphafold structure for {self.uniprot_id}")
        self.output_dir.mkdir(parents=True, exist_ok=True)  # Ensure dir exists

        versions_to_try = ["v4", "v3"]
        pdb_path = None
        last_error = None  # Keep track of the last error encountered

        for version in versions_to_try:
            model_version_suffix = f"model_{version}"
            # Construct URL using f-string
            url = (
                f"https://alphafold.ebi.ac.uk/files/"
                f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            )
            output_filename = f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            output_path = self.output_dir / output_filename

            self.logger.info(f"Attempting to download version {version} from: {url}")

            try:
                response = requests.get(url, timeout=60)
                response.raise_for_status()  # Check for HTTP errors (4xx, 5xx)

                # Explicit check for 200 OK (though raise_for_status should cover it)
                if response.status_code == 200:
                    # Write the downloaded content to the file
                    output_path.write_text(response.text, encoding="utf-8")
                    self.logger.info(f"Structure ({version}) saved to {output_path}")
                    pdb_path = output_path
                    last_error = None  # Reset error on success
                    break  # Success, exit the loop

            except requests.exceptions.HTTPError as e:
                last_error = e  # Store the error
                if e.response.status_code == 404:
                    self.logger.info(
                        f"Model version {version} not found for "
                        + f"{self.uniprot_id} (404)."
                    )
                else:
                    # Log other HTTP errors as warnings/errors
                    self.logger.warning(
                        f"HTTP error fetching structure version {version} for "
                        f"{self.uniprot_id}. Status Code: {e.response.status_code}. "
                        f"Error: {e}"
                    )
                # Continue to the next version attempt
            except requests.exceptions.RequestException as e:
                last_error = e  # Store the error
                # Handle other request errors (timeout, connection error, etc.)
                self.logger.warning(
                    f"Failed to fetch structure version {version} for "
                    f"{self.uniprot_id} due to a network or request issue. "
                    f"Error: {e}"
                )
                # Continue to the next version attempt

        # After trying all versions, check if we were successful
        if pdb_path:
            return pdb_path
        else:
            # If we failed, raise an error
            error_message = (
                f"Failed to fetch any structure version for {self.uniprot_id} "
                f"after trying {versions_to_try}."
            )
            self.logger.error(error_message)
            if last_error:
                # If we have a specific error from the last attempt, raise it
                # This provides more context than a generic FileNotFoundError
                self.logger.error(f"Last error encountered: {last_error}")
                # Raise FileNotFoundError but include details from the last error
                raise FileNotFoundError(
                    f"{error_message} Last error: {last_error}"
                ) from last_error
            else:
                # Should not happen if versions_to_try is not empty, but as a fallback
                raise FileNotFoundError(error_message)
