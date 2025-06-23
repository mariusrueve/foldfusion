import logging
from pathlib import Path

from .utils import Config
from .tools import (
    AlphaFoldFetcher,
    Dogsite3,
    SienaDB,
    Siena,
    LigandExtractor,
    JamdaScorer,
)
from .evaluation import Evaluator

logger = logging.getLogger(__name__)


class FoldFusion:
    """The heart of the foldfusion pipeline. It handles every step of the pipeline,
    gahters all meta data and evaluation data and save all into one directory.
    """

    def __init__(self, config_path: str):
        self.config_path = Path(config_path)

        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file does not exist: {self.config_path}")

        if not self.config_path.is_file():
            raise ValueError(f"Config path is not a file: {self.config_path}")

        try:
            self.config = Config(self.config_path)
        except Exception as e:
            logger.error(f"Configuration error: {e}")
            raise RuntimeError(
                f"Failed to initialize FoldFusion due to configuration error: {str(e)}"
            )
        logger.debug(f"Initialized FoldFusion with config: {self.config}")
        self.evaluator = Evaluator(self.config.output_dir)
        logger.debug(f"Evaluator initialized with path: {self.config.output_dir}")

    def run(self):
        main_output_dir = self.config.output_dir
        # Check if existing SienaDB is valid before creating a new one
        siena_db = SienaDB(
            self.config.siena_db_executable,
            self.config.siena_db_database_path,
            self.config.pdb_directory,
            self.config.pdb_format,
            main_output_dir,
        )
        self.siena_db_database_path = siena_db.run()
        skipped_uniprot_ids = []
        for uniprot_id in self.config.uniprot_ids:
            logger.info(f"Starting pipeline for UniProt ID: {uniprot_id}")
            output_dir = main_output_dir / uniprot_id
            try:
                self._pipeline(uniprot_id, output_dir)
            except Exception as e:
                logger.error(f"Error processing UniProt ID {uniprot_id}: {e}")
                skipped_uniprot_ids.append(uniprot_id)
                continue
        if skipped_uniprot_ids:
            logger.warning(f"Skipped UniProt IDs due to errors: {skipped_uniprot_ids}")
        else:
            logger.info("All UniProt IDs processed successfully")

    def _pipeline(self, uniprot_id: str, output_dir: Path):
        af_fetcher = AlphaFoldFetcher(uniprot_id, output_dir)
        af_model_path = af_fetcher.get_alphafold_model()
        logger.info(f"AlphaFold model was saved to {af_model_path}")

        dg3 = Dogsite3(self.config.dogsite3_executable, af_model_path, output_dir)
        dg3.run()
        best_edf = dg3.get_best_edf()

        # Run Siena with the determined database path
        siena = Siena(
            self.config.siena_executable,
            best_edf,
            self.siena_db_database_path,
            self.config.pdb_directory,
            output_dir,
        )
        siena.run()
        best_alignments = siena.get_best_alignments(10)
        logger.info(
            f"Pipeline completed successfully. Best alignments: {best_alignments}"
        )
        ligand_ex = LigandExtractor(
            self.config.ligand_extractor_executable,
            siena.output_dir,
            best_alignments,
            output_dir,
        )
        ligand_ex_output_path = ligand_ex.run()
        logger.info(f"Ligand extraction completed. Results saved to {output_dir}")

        jamda_scorer = JamdaScorer(
            self.config.jamda_scorer_executable,
            af_model_path,
            ligand_ex_output_path,
            output_dir,
        )
        jamda_scorer_output_path = jamda_scorer.run()
        logger.info(
            f"JAMDA scoring completed. Results saved to {jamda_scorer_output_path}"
        )


if __name__ == "__main__":
    ff = FoldFusion("config.toml")
