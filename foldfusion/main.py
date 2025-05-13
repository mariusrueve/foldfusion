"""Main script to run the FoldFusion pipeline.

This script orchestrates the different tools in the FoldFusion pipeline:
1. Fetches AlphaFold models.
2. Runs DogSite3 for pocket detection.
3. Runs SIENA for ensemble generation.
4. Extracts ligands from SIENA results.
"""

from pathlib import Path
from foldfusion.utils.config import Config
from foldfusion.tools.dogsite3 import Dogsite3
from foldfusion.tools.alphafoldfetcher import AlphaFoldFetcher
from foldfusion.tools.siena import Siena
from foldfusion.tools.ligand_extractor import LigandExtractor
from foldfusion.tools.jamda_scorer import JamdaScorer
import logging

# Setup basic logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)  # Use __name__ for logger


def main():
    """Runs the FoldFusion pipeline using a TOML configuration file."""
    config_path = Path("config.toml")
    logger.info(f"Loading configuration from {config_path.resolve()}")
    config = Config(config_path)
    config_dict = config.dict

    af_fetcher = AlphaFoldFetcher(config_dict)
    logger.info(f"Fetching AlphaFold structure for UniProtID: {af_fetcher.uniprot_id}")
    af_output_path = af_fetcher.get_alphafold_model()
    logger.info(f"AlphaFold model saved to: {af_output_path}")

    logger.info("Starting DogSite3 for binding site prediction...")
    dg3 = Dogsite3(config_dict, af_output_path)
    dg3.run()
    best_edf_path = dg3.get_best_edf()
    logger.info(f"Best EDF file from DogSite3: {best_edf_path}")

    logger.info("Starting SIENA for ensemble generation...")
    siena = Siena(config_dict, best_edf_path)
    siena.run()
    best_alignments_pdb_codes = siena.get_best_alignments(n_alignments=4)
    logger.info(
        f"Retrieved {len(best_alignments_pdb_codes)} best alignments from SIENA."
    )

    logger.info("Starting LigandExtractor to extract ligands...")
    le = LigandExtractor(config_dict, siena.output_dir, best_alignments_pdb_codes)
    le.run_all()
    logger.info("Ligand extraction completed.")
    logger.info("FoldFusion pipeline finished successfully.")

    logger.info("Starting JamdaScorer to optimize ligands ...")
    js = JamdaScorer(config_dict, af_output_path, le.output_dir)
    js.run_all()


if __name__ == "__main__":
    main()
