import logging
from pathlib import Path

from foldfusion.utils import Config
from foldfusion.tools import (
    AlphaFoldFetcher,
    Dogsite3,
    Siena,
    LigandExtractor,
    JamdaScorer,
)

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

    def run(self):
        uniprot_id = self.config.uniprot_id
        output_dir = self.config.output_dir  # This is now a Path object
        af_fetcher = AlphaFoldFetcher(uniprot_id, output_dir)
        af_model_path = af_fetcher.get_alphafold_model()
        logger.info(f"AlphaFold model was saved to {af_model_path}")

        # dg3 = Dogsite3()


if __name__ == "__main__":
    ff = FoldFusion("config.toml")
