"""
FoldFusion Pipeline Core Module

This module contains the main FoldFusion class that orchestrates the complete
protein structure enrichment pipeline. The pipeline combines multiple computational
tools to transplant ligands from experimental PDB structures into AlphaFold models.

The FoldFusion class manages the following key steps:
1. AlphaFold structure retrieval and preprocessing
2. Binding site prediction using DoGSite3
3. Structure alignment using SIENA
4. Ligand extraction and transplantation
5. Ligand optimization using JAMDA scorer
6. Comprehensive evaluation and quality assessment

Each step is implemented as a separate method to ensure modularity and
facilitate debugging and testing.
"""

import logging
from pathlib import Path

from .evaluation import Evaluator
from .tools import (
    AlphaFoldFetcher,
    Dogsite3,
    JamdaScorer,
    LigandExtractor,
    Siena,
    SienaDB,
)
from .utils import Config, setup_logging

logger = logging.getLogger(__name__)


class FoldFusion:
    """
    Main orchestrator for the FoldFusion protein structure enrichment pipeline.

    This class manages the complete workflow for transplanting ligands from
    experimental PDB structures into AlphaFold protein models. It coordinates
    multiple bioinformatics tools and ensures proper data flow between pipeline
    stages.

    The pipeline processes multiple UniProt IDs sequentially, handling errors
    gracefully and providing comprehensive logging of all operations. Results
    are stored in a structured output directory with separate subdirectories
    for each processed protein.

    Attributes:
        config_path: Path to the TOML configuration file
        config: Configuration object containing all pipeline parameters
        evaluator: Evaluation manager for quality assessment
        siena_db_database_path: Path to the SIENA database used for alignments

    Example:
        Basic usage with default configuration:

        >>> from foldfusion import FoldFusion
        >>> pipeline = FoldFusion("config.toml")
        >>> pipeline.run()

        The pipeline will process all UniProt IDs specified in the configuration
        file and store results in the configured output directory.
    """

    def __init__(self, config_path: str) -> None:
        """
        Initialize the FoldFusion pipeline with configuration.

        Args:
            config_path: Path to the TOML configuration file containing all
                        pipeline parameters including input UniProt IDs,
                        executable paths, and output directories.

        Raises:
            FileNotFoundError: If the configuration file does not exist.
            ValueError: If the configuration path is not a valid file.
            RuntimeError: If configuration loading fails due to syntax errors
                         or missing required parameters.

        Note:
            The constructor also initializes the logging system based on
            configuration parameters and sets up the evaluation framework.
        """
        self.config_path = Path(config_path)

        # Validate configuration file exists and is accessible
        if not self.config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.config_path}")

        if not self.config_path.is_file():
            raise ValueError(f"Configuration path is not a file: {self.config_path}")

        # Load and validate configuration
        try:
            self.config = Config(self.config_path)
            # Initialize logging as early as possible
            setup_logging(self.config.log_level, self.config.log_file)
            logger.info(
                f"FoldFusion pipeline initialized with configuration: "
                f"{self.config_path}"
            )
        except Exception as e:
            # Use print for initial error since logging may not be configured yet
            error_msg = f"Failed to initialize FoldFusion: {str(e)}"
            print(f"ERROR: {error_msg}")
            logger.error(f"Configuration initialization failed: {e}")
            raise RuntimeError(error_msg) from e

        # Initialize evaluation framework
        try:
            self.evaluator = Evaluator(self.config.output_dir)
            logger.info(
                f"Evaluation framework initialized with output directory: "
                f"{self.config.output_dir}"
            )
        except Exception as e:
            logger.error(f"Failed to initialize evaluator: {e}")
            raise RuntimeError(f"Evaluator initialization failed: {str(e)}") from e

        # Initialize database path placeholder
        self.siena_db_database_path: Path | None = None

        logger.debug("FoldFusion initialization completed successfully")

    def _fetch_alphafold_model(self, uniprot_id: str, output_dir: Path) -> Path:
        """
        Retrieve and preprocess AlphaFold structure for the given UniProt ID.

        This method downloads the AlphaFold model from the database and applies
        preprocessing steps including reliability filtering based on pLDDT scores.

        Args:
            uniprot_id: UniProt identifier for the target protein.
            output_dir: Directory where the AlphaFold model will be stored.

        Returns:
            Path to the processed AlphaFold structure file.

        Raises:
            RuntimeError: If the AlphaFold model cannot be retrieved or processed.

        Note:
            The method creates an 'AlphaFold' subdirectory within the output
            directory to organize downloaded structures.
        """
        logger.info(
            f"Initiating AlphaFold model retrieval for UniProt ID: {uniprot_id}"
        )

        try:
            af_fetcher = AlphaFoldFetcher(uniprot_id, output_dir)
            af_model_path = af_fetcher.get_alphafold_model()

            logger.info(
                f"AlphaFold model successfully retrieved and saved to: {af_model_path}"
            )
            logger.debug(
                f"AlphaFold structure file size: {af_model_path.stat().st_size} bytes"
            )

            return af_model_path

        except Exception as e:
            error_msg = f"Failed to fetch AlphaFold model for {uniprot_id}: {str(e)}"
            logger.error(error_msg)
            raise RuntimeError(error_msg) from e

    def _run_dogsite3(self, af_model_path: Path, output_dir: Path) -> Path:
        """
        Execute DoGSite3 binding site prediction on the AlphaFold structure.

        DoGSite3 identifies potential binding sites by analyzing geometric and
        physicochemical properties of the protein surface. The method returns
        the highest-scoring binding site for subsequent analysis.

        Args:
            af_model_path: Path to the preprocessed AlphaFold structure.
            output_dir: Directory where DoGSite3 results will be stored.

        Returns:
            Path to the best binding site EDF file identified by DoGSite3.

        Raises:
            RuntimeError: If DoGSite3 execution fails or no valid binding sites
                         are found.

        Note:
            DoGSite3 creates multiple output files; this method selects the
            binding site with the highest confidence score.
        """
        logger.info(
            f"Starting DoGSite3 binding site prediction for: {af_model_path.name}"
        )

        try:
            dogsite3 = Dogsite3(
                self.config.dogsite3_executable, af_model_path, output_dir
            )
            dogsite3_output_dir = dogsite3.run()
            best_edf = dogsite3.get_best_edf()

            logger.info("DoGSite3 analysis completed successfully")
            logger.info(f"Best binding site identified: {best_edf.name}")
            logger.debug(f"DoGSite3 output directory: {dogsite3_output_dir}")

            return best_edf

        except Exception as e:
            error_msg = f"DoGSite3 execution failed: {str(e)}"
            logger.error(error_msg)
            raise RuntimeError(error_msg) from e

    def _run_siena(self, best_edf: Path, output_dir: Path) -> tuple[Siena, list]:
        """
        Execute SIENA structure alignment to find template structures.

        SIENA performs structure-based alignment between the predicted binding site
        and known experimental structures in the PDB database. This identifies
        template structures with similar binding sites that can provide ligands
        for transplantation.

        Args:
            best_edf: Path to the EDF file defining the target binding site.
            output_dir: Directory where SIENA alignment results will be stored.

        Returns:
            A tuple containing:
                - siena: The SIENA tool instance with alignment results
                - best_alignments: List of best structural alignments found

        Raises:
            RuntimeError: If SIENA execution fails or no suitable alignments
                         are found.
            ValueError: If the SIENA database path is not properly initialized.

        Note:
            The number of returned alignments is limited by the configuration
            parameter siena_max_alignments to ensure manageable processing times.
        """
        if self.siena_db_database_path is None:
            raise ValueError("SIENA database path not initialized")

        logger.info(
            f"Starting SIENA structure alignment for binding site: {best_edf.name}"
        )

        try:
            siena = Siena(
                self.config.siena_executable,
                best_edf,
                self.siena_db_database_path,
                self.config.pdb_directory,
                output_dir,
            )
            siena_output_dir = siena.run()
            best_alignments = siena.get_best_alignments(
                self.config.siena_max_alignments
            )

            logger.info("SIENA alignment completed successfully")
            logger.info(f"Found {len(best_alignments)} high-quality alignments")
            logger.debug(f"SIENA output directory: {siena_output_dir}")
            logger.debug(f"Best alignments: {best_alignments}")

            return siena, best_alignments

        except Exception as e:
            error_msg = f"SIENA alignment failed: {str(e)}"
            logger.error(error_msg)
            raise RuntimeError(error_msg) from e

    def _run_ligand_extraction(
        self, siena: Siena, best_alignments: list, output_dir: Path
    ) -> dict:
        """
        Extract ligands from template structures based on SIENA alignments.

        This method uses the structural alignments to identify and extract ligands
        from experimental structures that correspond to the predicted binding site.
        The extracted ligands are then prepared for transplantation into the
        target AlphaFold structure.

        Args:
            siena: SIENA tool instance containing alignment information.
            best_alignments: List of structural alignments to process.
            output_dir: Directory where extracted ligands will be stored.

        Returns:
            Dictionary containing ligand structure information organized by
            PDB code and ligand identifier.

        Raises:
            RuntimeError: If ligand extraction fails or no valid ligands are found.

        Note:
            The method filters out ligands that are too small or have
            incompatible chemical properties for meaningful analysis.
        """
        logger.info(
            f"Starting ligand extraction from {len(best_alignments)} alignments"
        )

        try:
            ligand_extractor = LigandExtractor(
                self.config.ligand_extractor_executable,
                siena.output_dir,
                best_alignments,
                output_dir,
            )
            ligand_extractor_output = ligand_extractor.run()
            ligand_structure = ligand_extractor.ligand_structure

            logger.info("Ligand extraction completed successfully")
            logger.info(f"Extracted ligands saved to: {output_dir}")
            logger.info(f"Total ligands extracted: {len(ligand_structure)}")
            logger.debug(f"Ligand extractor output: {ligand_extractor_output}")
            logger.debug(f"Ligand structure details: {ligand_structure}")

            return ligand_structure

        except Exception as e:
            error_msg = f"Ligand extraction failed: {str(e)}"
            logger.error(error_msg)
            raise RuntimeError(error_msg) from e

    def _run_jamda_scorer(
        self, af_model_path: Path, ligand_structure: dict, output_dir: Path
    ) -> dict:
        """
        Optimize ligand positions using the JAMDA scoring system.

        JAMDA performs energy-based optimization of ligand poses within the
        target protein structure. This step refines the initial ligand
        placements to minimize steric clashes and optimize binding interactions.

        Args:
            af_model_path: Path to the target AlphaFold protein structure.
            ligand_structure: Dictionary containing ligand information from extraction.
            output_dir: Directory where optimized ligands will be stored.

        Returns:
            Dictionary containing optimized ligand structure information with
            updated coordinates and energy scores.

        Raises:
            RuntimeError: If JAMDA optimization fails or produces invalid results.

        Note:
            The optimization process may reject ligands that cannot be
            satisfactorily placed without severe steric conflicts.
        """
        logger.info(f"Starting JAMDA optimization for {len(ligand_structure)} ligands")

        try:
            jamda_scorer = JamdaScorer(
                self.config.jamda_scorer_executable,
                af_model_path,
                ligand_structure,
                output_dir,
            )
            optimized_ligand_structure = jamda_scorer.run()

            logger.info("JAMDA optimization completed successfully")
            logger.info(f"Optimized ligands: {len(optimized_ligand_structure)}")
            logger.debug(f"Optimization results: {optimized_ligand_structure}")

            return optimized_ligand_structure

        except Exception as e:
            error_msg = f"JAMDA optimization failed: {str(e)}"
            logger.error(error_msg)
            raise RuntimeError(error_msg) from e

    def _evaluate_and_log(
        self,
        uniprot_id: str,
        evaluation_name: str,
        af_model_path: Path,
        best_alignments: list,
        ligand_structure: dict,
    ) -> None:
        """
        Perform quality evaluation and log results for the current pipeline stage.

        This method computes various quality metrics including Local RMSD and
        Transplant Clash Score (TCS) to assess the quality of structural
        alignments and ligand placements. Results are stored for later analysis.

        Args:
            uniprot_id: UniProt identifier for the target protein.
            evaluation_name: Descriptive name for this evaluation stage.
            af_model_path: Path to the target protein structure.
            best_alignments: List of structural alignments to evaluate.
            ligand_structure: Dictionary containing ligand information.

        Raises:
            RuntimeError: If evaluation fails due to data processing errors.

        Note:
            Evaluation results are automatically saved to the configured
            output directory in JSON format for subsequent analysis.
        """
        logger.info(f"Starting evaluation stage: {evaluation_name} for {uniprot_id}")

        try:
            self.evaluator.evaluate(
                uniprot_id,
                evaluation_name,
                af_model_path,
                best_alignments,
                ligand_structure,
            )

            logger.info(f"Evaluation '{evaluation_name}' completed successfully")
            logger.debug(f"Evaluation data saved for UniProt ID: {uniprot_id}")

        except Exception as e:
            error_msg = f"Evaluation failed for {evaluation_name}: {str(e)}"
            logger.error(error_msg)
            raise RuntimeError(error_msg) from e

    def run(self) -> None:
        """
        Execute the complete FoldFusion pipeline for all configured UniProt IDs.

        This method orchestrates the entire pipeline workflow, processing each
        UniProt ID sequentially. It handles errors gracefully, allowing the
        pipeline to continue processing remaining proteins even if individual
        proteins fail.

        The pipeline includes the following stages for each protein:
        1. SIENA database preparation (once for all proteins)
        2. AlphaFold structure retrieval and preprocessing
        3. Binding site prediction with DoGSite3
        4. Structure alignment with SIENA
        5. Ligand extraction and initial evaluation
        6. Ligand optimization with JAMDA and final evaluation

        Raises:
            RuntimeError: If critical pipeline initialization fails.

        Note:
            Failed protein processing is logged but does not terminate the
            entire pipeline. A summary of skipped proteins is provided at
            completion.
        """
        logger.info("Starting FoldFusion pipeline execution")
        logger.info(f"Processing {len(self.config.uniprot_ids)} UniProt IDs")

        main_output_dir = self.config.output_dir

        # TODO: When the defined siena db from config.toml is already present with a
        # sufficient size the generate_database step can be skipped completely

        # Initialize SIENA database (shared across all proteins)
        try:
            logger.info("Initializing SIENA database")
            siena_db = SienaDB(
                self.config.siena_db_executable,
                # May be None -> SienaDB will use default path
                self.config.siena_db_database_path,
                self.config.pdb_directory,
                self.config.pdb_format,
                main_output_dir,
            )
            self.siena_db_database_path = siena_db.run()
            logger.info(f"SIENA database initialized at: {self.siena_db_database_path}")
        except Exception as e:
            error_msg = f"Failed to initialize SIENA database: {str(e)}"
            logger.error(error_msg)
            raise RuntimeError(error_msg) from e

        # Process each UniProt ID
        skipped_uniprot_data = {}  # Dictionary to store UniProt ID -> error reason
        successful_proteins = 0

        for i, uniprot_id in enumerate(self.config.uniprot_ids, 1):
            logger.info(
                f"Processing UniProt ID {i}/{len(self.config.uniprot_ids)}: "
                f"{uniprot_id}"
            )
            output_dir = main_output_dir / "Results" / uniprot_id

            try:
                self._pipeline(uniprot_id, output_dir)
                successful_proteins += 1
                logger.info(f"Successfully completed processing for {uniprot_id}")

            except Exception as e:
                error_reason = str(e)
                logger.error(
                    "Failed to process UniProt ID %s: %s", uniprot_id, error_reason
                )
                skipped_uniprot_data[uniprot_id] = error_reason
                continue

        # Pipeline completion summary
        logger.info("FoldFusion pipeline execution completed")
        logger.info(f"Successfully processed: {successful_proteins} proteins")

        if skipped_uniprot_data:
            skipped_ids = list(skipped_uniprot_data.keys())
            logger.warning(
                f"Skipped {len(skipped_ids)} UniProt IDs due to errors: {skipped_ids}"
            )
            # Log detailed error reasons for each failed UniProt ID
            logger.warning("Detailed failure reasons:")
            for uniprot_id, error_reason in skipped_uniprot_data.items():
                logger.warning(f"  - {uniprot_id}: {error_reason}")
        else:
            logger.info("All UniProt IDs processed successfully")

    def _pipeline(self, uniprot_id: str, output_dir: Path) -> None:
        """
        Execute the complete pipeline workflow for a single protein.

        This method processes a single UniProt ID through all pipeline stages,
        from AlphaFold structure retrieval to final ligand optimization and
        evaluation. It maintains consistent error handling and logging throughout.

        Args:
            uniprot_id: UniProt identifier for the target protein.
            output_dir: Directory where all results for this protein will be stored.

        Raises:
            RuntimeError: If any critical pipeline stage fails.

        Note:
            The method creates evaluation checkpoints after ligand extraction
            (pre-JAMDA) and after optimization (post-JAMDA) to enable
            comparative analysis of results.
        """
        logger.info(f"Starting individual pipeline for {uniprot_id}")

        # Stage 1: AlphaFold structure retrieval
        af_model_path = self._fetch_alphafold_model(uniprot_id, output_dir)

        # Stage 2: Binding site prediction
        best_edf = self._run_dogsite3(af_model_path, output_dir)

        # Stage 3: Structure alignment
        siena, best_alignments = self._run_siena(best_edf, output_dir)

        # Stage 4: Ligand extraction and initial evaluation
        ligand_structure = self._run_ligand_extraction(
            siena, best_alignments, output_dir
        )
        self._evaluate_and_log(
            uniprot_id,
            "pre-jamda",
            af_model_path,
            best_alignments,
            ligand_structure,
        )

        # Stage 5: Ligand optimization and final evaluation
        optimized_ligand_structure = self._run_jamda_scorer(
            af_model_path, ligand_structure, output_dir
        )
        self._evaluate_and_log(
            uniprot_id,
            "post-jamda",
            af_model_path,
            best_alignments,
            optimized_ligand_structure,
        )

        logger.info(f"Completed pipeline for {uniprot_id}")


if __name__ == "__main__":
    # Entry point for direct execution
    pipeline = FoldFusion("config.toml")
    pipeline.run()
