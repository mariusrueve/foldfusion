"""Module for fetching protein structures from the AlphaFold Database."""

from foldfusion.tools.tool import Tool  # Assuming Tool is defined elsewhere
from pathlib import Path
import requests
import logging

logger = logging.getLogger(__name__)


class AlphaFoldFetcher(Tool):
    """Fetches protein structures from the AlphaFold Database and processes them.

    Attributes:
        uniprot_id (str): The UniProt ID of the protein to fetch.
        output_dir (Path): The directory where the fetched PDB file will be saved.
                           The base output_dir is typically set by the Tool base class
                           or its configuration. This class will use a subdirectory
                           named 'alphafold' within that output_dir.
    """

    def __init__(self, config: dict):
        """Initializes the AlphaFoldFetcher.

        Args:
            config (dict): The configuration dictionary.
                           Expected to contain 'run_settings': {'uniprot_id': '...'}
                           and potentially 'output_dir' for the base Tool.
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

        # self.output_dir is inherited from Tool. We create a subdirectory for AlphaFold files.
        self.alphafold_output_dir = self.output_dir / "alphafold"
        logger.debug(f"AlphaFoldFetcher initialized for UniProt ID: {self.uniprot_id}")
        logger.debug(
            f"Output directory for AlphaFold files: {self.alphafold_output_dir}"
        )

    def _process_pdb_file(self, pdb_file_path: Path) -> Path:
        """
        Processes a PDB file to remove unreliable residues based on pLDDT scores.

        A residue is considered unreliable if any of its atoms have a pLDDT < 70.
        If the percentage of unreliable residues is <= 15%, ATOM lines for these
        residues are removed, and a new file with '_processed' suffix is saved.
        Otherwise, the original file path is returned.

        Args:
            pdb_file_path (Path): The path to the input PDB file.

        Returns:
            Path: The path to the processed PDB file, or the original path if
                  processing was skipped or failed.
        """
        logger.info(f"Starting PDB processing for: {pdb_file_path}")

        try:
            pdb_lines = pdb_file_path.read_text(encoding="utf-8").splitlines()
        except Exception as e:
            logger.error(f"Could not read PDB file {pdb_file_path}: {e}")
            return pdb_file_path

        unique_residues = set()
        unreliable_residue_ids = set()  # Stores (chain_id, res_seq_num, i_code) tuples

        for line in pdb_lines:
            if line.startswith("ATOM  "):  # Standard PDB ATOM record
                try:
                    # PDB format:
                    # Chain ID: column 22 (index 21)
                    # Residue sequence number: columns 23-26 (index 22-25)
                    # Insertion code: column 27 (index 26)
                    # B-factor (pLDDT): columns 61-66 (index 60-65)

                    chain_id = line[21:22].strip()
                    res_seq_num = line[22:26].strip()
                    i_code = line[26:27].strip()  # Insertion code

                    residue_id = (chain_id, res_seq_num, i_code)
                    unique_residues.add(residue_id)

                    plddt_str = line[60:66].strip()
                    plddt = float(plddt_str)

                    if plddt < 70.0:
                        unreliable_residue_ids.add(residue_id)
                except ValueError:
                    logger.warning(
                        f"Could not parse pLDDT or residue info from ATOM line: '{line.strip()}' in {pdb_file_path}"
                    )
                except IndexError:
                    logger.warning(
                        f"Malformed ATOM line (IndexError): '{line.strip()}' in {pdb_file_path}"
                    )

        if not unique_residues:
            logger.warning(
                f"No residues found in {pdb_file_path}. Skipping processing."
            )
            return pdb_file_path

        num_total_residues = len(unique_residues)
        num_unreliable_residues = len(unreliable_residue_ids)

        percentage_unreliable = (num_unreliable_residues / num_total_residues) * 100

        logger.info(
            f"Residue analysis for {pdb_file_path.name}: "
            f"Total unique residues: {num_total_residues}, "
            f"Unreliable residues (pLDDT < 70): {num_unreliable_residues} "
            f"({percentage_unreliable:.2f}%)"
        )

        if percentage_unreliable <= 15.00:
            processed_file_name = (
                f"{pdb_file_path.stem}_processed{pdb_file_path.suffix}"
            )
            processed_pdb_path = pdb_file_path.with_name(processed_file_name)

            logger.info(
                f"Percentage of unreliable residues ({percentage_unreliable:.2f}%) is <= 15%. "
                f"Writing processed file to: {processed_pdb_path}"
            )

            try:
                with open(processed_pdb_path, "w", encoding="utf-8") as outfile:
                    for line in pdb_lines:
                        if line.startswith("ATOM  "):
                            try:
                                chain_id = line[21:22].strip()
                                res_seq_num = line[22:26].strip()
                                i_code = line[26:27].strip()
                                current_residue_id = (chain_id, res_seq_num, i_code)
                                if current_residue_id in unreliable_residue_ids:
                                    continue  # Skip ATOM lines of unreliable residues
                            except (
                                IndexError
                            ):  # Should not happen if parsed above, but as safeguard
                                logger.warning(
                                    f"Skipping malformed ATOM line during write: {line.strip()}"
                                )
                                outfile.write(line + "\n")  # Write if unsure
                                continue

                        outfile.write(line + "\n")
                logger.info(
                    f"Successfully wrote processed PDB file: {processed_pdb_path}"
                )
                return processed_pdb_path
            except IOError as e:
                logger.error(
                    f"Failed to write processed PDB file {processed_pdb_path}: {e}"
                )
                return pdb_file_path  # Return original if writing fails
        else:
            logger.info(
                f"Percentage of unreliable residues ({percentage_unreliable:.2f}%) > 15%. "
                f"No residues removed. Using original file: {pdb_file_path}"
            )
            return pdb_file_path

    def get_alphafold_model(self) -> Path:
        """Downloads the AlphaFold PDB model and processes it.

        Tries fetching model versions v4 and then v3.
        The downloaded model is then processed to remove unreliable residues
        if they constitute 15% or less of the total residues.

        Returns:
            Path: The path to the downloaded (and potentially processed) PDB file.

        Raises:
            FileNotFoundError: If the PDB file cannot be fetched from AlphaFold DB.
            requests.exceptions.RequestException: For network or HTTP errors during download.
        """
        self.alphafold_output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(
            f"Ensured AlphaFold output directory exists: {self.alphafold_output_dir}"
        )

        versions = ["v4", "v3"]  # Prioritize v4
        downloaded_pdb_path = None
        last_error = None

        for version in versions:
            model_version_suffix = f"model_{version}"
            url = (
                f"https://alphafold.ebi.ac.uk/files/"
                f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            )
            output_file_name = f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            # Save to the specific alphafold subdirectory
            output_file_path = self.alphafold_output_dir / output_file_name

            logger.info(
                f"Attempting to download AlphaFold model version {version} from {url}"
            )
            try:
                response = requests.get(url, timeout=60)
                response.raise_for_status()  # Raises HTTPError for bad responses

                output_file_path.write_text(response.text, encoding="utf-8")
                logger.info(
                    f"Successfully downloaded and saved model version {version} to {output_file_path}"
                )
                downloaded_pdb_path = output_file_path
                last_error = None
                break
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
                    f"{self.uniprot_id} due to a network or request issue. Error: {e}"
                )
                # For critical network errors, might be best to stop trying other versions
                # break

        if downloaded_pdb_path:
            try:
                # Process the downloaded PDB file
                final_pdb_path = self._process_pdb_file(downloaded_pdb_path)
                resolved_path = final_pdb_path.resolve()
                logger.info(
                    f"Returning resolved path for model (potentially processed): {resolved_path}"
                )
                return resolved_path
            except Exception as e:
                # If processing fails, log the error and return the original downloaded path
                logger.error(
                    f"Error processing PDB file {downloaded_pdb_path}: {e}. "
                    f"Returning original downloaded file."
                )
                return downloaded_pdb_path.resolve()
        else:
            error_message = (
                f"Failed to fetch any AlphaFold structure version for {self.uniprot_id} "
                f"after trying versions: {', '.join(versions)}."
            )
            logger.error(error_message)
            if last_error:
                logger.error(f"Last error encountered: {last_error}")
                raise FileNotFoundError(
                    f"{error_message} Last error: {last_error.__class__.__name__} - {str(last_error)}"
                ) from last_error
            else:
                raise FileNotFoundError(error_message)
