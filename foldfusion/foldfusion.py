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
from .utils import Config

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

    def _fetch_alphafold_model(self, uniprot_id: str, output_dir: Path) -> Path:
        """Fetches the AlphaFold model."""
        af_fetcher = AlphaFoldFetcher(uniprot_id, output_dir)
        af_model_path = af_fetcher.get_alphafold_model()
        logger.info(f"AlphaFold model was saved to {af_model_path}")
        return af_model_path

    def _run_dogsite3(self, af_model_path: Path, output_dir: Path) -> Path:
        """Runs Dogsite3 to find the best binding site."""
        dg3 = Dogsite3(self.config.dogsite3_executable, af_model_path, output_dir)
        dg3.run()
        return dg3.get_best_edf()

    def _run_siena(self, best_edf: Path, output_dir: Path) -> tuple[Siena, list]:
        """Runs Siena to find alignments."""
        siena = Siena(
            self.config.siena_executable,
            best_edf,
            self.siena_db_database_path,
            self.config.pdb_directory,
            output_dir,
        )
        siena.run()
        best_alignments = siena.get_best_alignments(self.config.siena_max_alignments)
        logger.info(f"Siena found best alignments: {best_alignments}")
        return siena, best_alignments

    def _run_ligand_extraction(
        self, siena: Siena, best_alignments: list, output_dir: Path
    ) -> dict:
        """Extracts ligands based on Siena alignments."""
        ligand_extractor = LigandExtractor(
            self.config.ligand_extractor_executable,
            siena.output_dir,
            best_alignments,
            output_dir,
        )
        ligand_extractor.run()
        logger.info(f"Ligand extraction completed. Results saved to {output_dir}")
        logger.debug(
            "The following ligands were extracted:\n"
            f" {ligand_extractor.ligand_structure}"
        )
        return ligand_extractor.ligand_structure

    def _run_jamda_scorer(
        self, af_model_path: Path, ligand_structure: dict, output_dir: Path
    ) -> dict:
        """Scores and optimizes ligands using JAMDA."""
        jamda_scorer = JamdaScorer(
            self.config.jamda_scorer_executable,
            af_model_path,
            ligand_structure,
            output_dir,
        )
        optimized_ligand_structure = jamda_scorer.run()
        logger.info(
            f"JAMDA scoring completed. Optimized ligands: {optimized_ligand_structure}"
        )
        return optimized_ligand_structure

    def _evaluate_and_log(
        self,
        uniprot_id: str,
        evaluation_name: str,
        af_model_path: Path,
        best_alignments: list,
        ligand_structure: dict,
    ):
        """Runs evaluation and logs the results."""
        self.evaluator.evaluate(
            uniprot_id,
            evaluation_name,
            af_model_path,
            best_alignments,
            ligand_structure,
        )
        logger.info(f"Evaluation '{evaluation_name}' completed.")

    def run(self):
        main_output_dir = self.config.output_dir
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
        af_model_path = self._fetch_alphafold_model(uniprot_id, output_dir)

        best_edf = self._run_dogsite3(af_model_path, output_dir)
        siena, best_alignments = self._run_siena(best_edf, output_dir)

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


if __name__ == "__main__":
    ff = FoldFusion("config.toml")
    ff.run()
