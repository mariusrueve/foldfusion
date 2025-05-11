from pathlib import Path
from foldfusion.utils.config import Config
from foldfusion.tools.dogsite3 import Dogsite3
from foldfusion.tools.alphafoldfetcher import AlphaFoldFetcher
from foldfusion.tools.siena import Siena
import logging

# Setup basic logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger("foldfusion")


def main():
    config = Config(Path("config.toml"))
    config = config.dict

    af_fetcher = AlphaFoldFetcher(config)
    logger.info(f"Fetching Alphafoldstructure for UniProtID: {af_fetcher.uniprot_id}")
    af_output_path = af_fetcher.get_alphafold_model()
    logger.info(f"AlphaFold model is saved to: {af_output_path}")

    logger.info("Starting Dogsite3 ...")
    dg3 = Dogsite3(config, af_output_path)
    dg3.run()
    best_edf_path = dg3.get_best_edf()

    logger.info("Starting Siena ...")
    siena = Siena(config, best_edf_path)
    siena.run()


if __name__ == "__main__":
    main()
