"""Module for fetching protein structures from the AlphaFold Database."""

import logging
from pathlib import Path

import requests

logger = logging.getLogger(__name__)


class AlphaFoldFetcher:
    def __init__(self, uniprot_id: str, output_dir: Path):
        self.uniprot_id = uniprot_id
        self.output_dir = output_dir
        self.af_fetcher_out_dir = self.output_dir / "AlphaFold"

    def _filter_unreliable_residues(self, pdb_file_path: Path) -> Path:
        """
        Processes a PDB file to remove unreliable residue segments based on
        pLDDT scores.

        This method removes segments of ≥5 consecutive residues with pLDDT < 50.
        This is more conservative than the previous approach and focuses on removing
        low-confidence segments rather than individual residues.

        Args:
            pdb_file_path (Path): The path to the input PDB file.

        Returns:
            Path: The path to the processed PDB file.
        """
        logger.info(f"Starting PDB processing for: {pdb_file_path}")

        try:
            pdb_lines = pdb_file_path.read_text(encoding="utf-8").splitlines()
        except Exception as e:
            logger.error(f"Could not read PDB file {pdb_file_path}: {e}")
            return pdb_file_path

        # Process PDB using segment-based trimming approach
        trimmed_lines = self._process_pdb_segments(pdb_lines)

        # Create processed file
        processed_file_name = (
            f"{pdb_file_path.stem}_processed{pdb_file_path.suffix}"
        )
        processed_pdb_path = pdb_file_path.with_name(processed_file_name)

        logger.info(f"Writing processed file to: {processed_pdb_path}")

        try:
            with open(processed_pdb_path, "w", encoding="utf-8") as outfile:
                for line in trimmed_lines:
                    outfile.write(line)
                    if not line.endswith('\n'):
                        outfile.write('\n')

            logger.info(
                f"Successfully wrote processed PDB file: {processed_pdb_path}"
            )
            return processed_pdb_path
        except OSError as e:
            logger.error(
                f"Failed to write processed PDB file {processed_pdb_path}: {e}"
            )
            return pdb_file_path  # Return original if writing fails

    def _process_pdb_segments(self, lines):
        """
        Processes PDB file by removing low-confidence segments.

        Removes segments with ≥5 consecutive residues with pLDDT < 50.

        Args:
            lines: List of PDB file lines

        Returns:
            List of processed PDB lines
        """
        residue_sequence = []
        trimmed = []
        bad_residues = set()
        current_segment = []

        # First pass: identify residue sequence and low-confidence segments
        for line in lines:
            if line.startswith("ATOM"):
                try:
                    # pLDDT is in columns 61-66 (0-indexed: 60-65)
                    plddt = float(line[60:66])
                    # Residue identifier (name + chain + number)
                    res_id = line[17:26].strip()
                    residue_sequence.append((res_id, plddt))
                except (ValueError, IndexError):
                    logger.warning(
                        f"Could not parse pLDDT from ATOM line: '{line.strip()}'"
                    )
                    continue

        # Find segments with ≥5 residues with pLDDT < 50
        for res_id, plddt in residue_sequence:
            if plddt < 50:
                if not current_segment or current_segment[-1] != res_id:
                    current_segment.append(res_id)
            else:
                if len(current_segment) >= 5:
                    bad_residues.update(current_segment)
                current_segment = []

        # Check last segment (end of sequence)
        if len(current_segment) >= 5:
            bad_residues.update(current_segment)

        logger.info(
            f"Found {len(bad_residues)} residues in low-confidence segments to remove"
        )

        # Second pass: remove detected segments
        for line in lines:
            if not line.startswith("ATOM"):
                trimmed.append(line)
                continue

            try:
                res_id = line[17:26].strip()
                if res_id not in bad_residues:
                    trimmed.append(line)
            except IndexError:
                logger.warning(f"Malformed ATOM line: '{line.strip()}'")
                trimmed.append(line)  # Keep line if parsing fails

        return trimmed

    def get_alphafold_model(self) -> Path:
        """Downloads the AlphaFold PDB model and processes it.

        Tries fetching model versions v4 and then v3.
        Models are processed using a segment-based trimming approach: segments of ≥5
        consecutive residues with pLDDT < 50 are removed. This is more conservative
        than removing individual low-confidence residues and focuses on removing
        entire unreliable segments.

        Returns:
            Path: The path to the downloaded and processed PDB file.

        Raises:
            FileNotFoundError: If no PDB file can be fetched from AlphaFold DB.
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
                # Process the downloaded PDB file using segment-based trimming
                processed_pdb_path = self._filter_unreliable_residues(
                    output_file_path
                )
                resolved_path = processed_pdb_path.resolve()

                logger.info(
                    f"Model {version} processed successfully. "
                    f"Returning path: {resolved_path}"
                )
                return resolved_path
            except Exception as e:
                # For processing errors, log and try next version
                logger.error(
                    f"Error processing model {version} PDB file {output_file_path}: {e}"
                )
                model_quality_errors.append((version, f"Processing error: {str(e)}"))
                continue

        # If we're here, none of the models worked
        if model_quality_errors:
            error_details = ", ".join([f"{v}: {e}" for v, e in model_quality_errors])
            error_message = (
                f"All AlphaFold models for {self.uniprot_id} failed during"
                + f" processing: {error_details}"
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
