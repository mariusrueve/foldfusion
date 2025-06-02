"""Module for fetching protein structures from the AlphaFold Database."""

from .tool import Tool  # Assuming Tool is defined elsewhere
from pathlib import Path
import requests
import logging

logger = logging.getLogger(__name__)


class AlphaFoldFetcher:
    def __init__(self, uniprot_id: str, output_dir: Path):
        self.uniprot_id = uniprot_id
        self.output_dir = output_dir
        self.af_fetcher_out_dir = self.output_dir / "alphafold"

    def _process_pdb_file(self, pdb_file_path: Path) -> Path:
        """
        Processes a PDB file to remove unreliable residues based on pLDDT scores.

        A residue is considered unreliable if any of its atoms have a pLDDT < 70.
        If the percentage of unreliable residues is <= 15%, ATOM lines for these
        residues are removed, and a new file with '_processed' suffix is saved.
        If more than 15% of residues are unreliable, the model is discarded by
        raising a ValueError.

        Args:
            pdb_file_path (Path): The path to the input PDB file.

        Returns:
            Path: The path to the processed PDB file.

        Raises:
            ValueError: If more than 15% of residues have pLDDT < 70.
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
            # Keep the model and process it by removing unreliable residues
            processed_file_name = (
                f"{pdb_file_path.stem}_processed{pdb_file_path.suffix}"
            )
            processed_pdb_path = pdb_file_path.with_name(processed_file_name)

            logger.info(
                f"Percentage of unreliable residues ({percentage_unreliable:.2f}%) is <= 15%. "
                f"Model is reliable. Writing processed file to: {processed_pdb_path}"
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
            # If more than 15% of residues are unreliable, discard the model
            logger.warning(
                f"Percentage of unreliable residues ({percentage_unreliable:.2f}%) > "
                f"15%. Model is considered unreliable and should be discarded."
            )
            # Raise an exception to indicate the model should be discarded
            raise ValueError(
                f"Model {pdb_file_path.name} has {percentage_unreliable:.2f}% "
                "unreliable residues (>15% threshold) and should be discarded."
            )

    def get_alphafold_model(self) -> Path:
        """Downloads the AlphaFold PDB model and processes it.

        Tries fetching model versions v4 and then v3.
        Models are processed to ensure they meet quality standards. Models with more
        than 15% unreliable residues (residues with pLDDT < 70) are discarded. For
        models with 15% or fewer unreliable residues, those residues are removed and a
        processed file is returned.

        Returns:
            Path: The path to the downloaded and processed PDB file.

        Raises:
            FileNotFoundError: If no PDB file can be fetched from AlphaFold DB.
            ValueError: If all fetched models have more than 15% unreliable residues.
            requests.exceptions.RequestException: For network or HTTP errors during
            download.
        """
        self.af_fetcher_out_dir.mkdir(parents=True, exist_ok=True)
        logger.info(
            f"Ensured AlphaFold output directory exists: {self.af_fetcher_out_dir}"
        )

        versions = ["v4", "v3"]  # Prioritize v4
        last_error = None
        model_quality_errors = []

        for version in versions:
            model_version_suffix = f"model_{version}"
            url = (
                f"https://alphafold.ebi.ac.uk/files/"
                f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            )
            output_file_name = f"AF-{self.uniprot_id}-F1-{model_version_suffix}.pdb"
            # Save to the specific alphafold subdirectory
            output_file_path = self.af_fetcher_out_dir / output_file_name

            logger.info(
                f"Attempting to download AlphaFold model version {version} from {url}"
            )

            # First try to download the model
            try:
                response = requests.get(url, timeout=60)
                response.raise_for_status()  # Raises HTTPError for bad responses
                output_file_path.write_text(response.text, encoding="utf-8")
                logger.info(
                    f"Successfully downloaded and saved model version {version} to "
                    f"{output_file_path}"
                )
            except requests.exceptions.HTTPError as e:
                last_error = e
                if e.response.status_code == 404:
                    logger.warning(
                        f"Model version {version} not found for UniProt ID "
                        f"{self.uniprot_id} (404 Error). Trying next version."
                    )
                    continue
                else:
                    logger.error(
                        f"HTTP error occurred while fetching model version {version} "
                        f"for {self.uniprot_id}. Status code: "
                        f"{e.response.status_code}. Error: {e}"
                    )
                    continue
            except requests.exceptions.RequestException as e:
                last_error = e
                logger.error(
                    f"Failed to fetch structure version {version} for "
                    f"{self.uniprot_id} due to a network or request issue. Error: {e}"
                )
                continue

            # If download successful, process the file
            try:
                # Process the downloaded PDB file - this will raise an exception if
                # the model has more than 15% unreliable residues
                processed_pdb_path = self._process_pdb_file(output_file_path)
                resolved_path = processed_pdb_path.resolve()
                logger.info(
                    f"Model {version} is reliable. Returning path for processed model:"
                    + f" {resolved_path}"
                )
                return resolved_path
            except ValueError as e:
                # This indicates the model has too many unreliable residues
                if "unreliable residues" in str(e):
                    logger.warning(f"{e}")
                    model_quality_errors.append((version, str(e)))
                    # Continue to next version
                    continue
                else:
                    # Re-raise unexpected ValueError
                    raise
            except Exception as e:
                # For other processing errors, log and try next version
                logger.error(
                    f"Error processing model {version} PDB file {output_file_path}: {e}"
                )
                model_quality_errors.append((version, f"Processing error: {str(e)}"))
                continue

        # If we're here, none of the models worked
        if model_quality_errors:
            error_details = ", ".join([f"{v}: {e}" for v, e in model_quality_errors])
            error_message = (
                f"All AlphaFold models for {self.uniprot_id} failed quality"
                + f" checks: {error_details}"
            )
            logger.error(error_message)
            raise ValueError(error_message)
        else:
            error_message = (
                f"Failed to fetch any AlphaFold structure version for {self.uniprot_id}"
                f" after trying versions: {', '.join(versions)}."
            )
            logger.error(error_message)
            if last_error:
                logger.error(f"Last error encountered: {last_error}")
                raise FileNotFoundError(
                    f"{error_message} Last error: {last_error.__class__.__name__} - "
                    + f"{str(last_error)}"
                ) from last_error
            else:
                raise FileNotFoundError(error_message)
