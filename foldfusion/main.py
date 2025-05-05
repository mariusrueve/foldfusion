from pathlib import Path
from foldfusion.utils.read_config import read_config
from foldfusion.tools.dogsite3 import Dogsite3
from foldfusion.tools.alphafoldfetcher import AlphaFoldFetcher
import logging


def main():
    config = read_config(Path("config.toml"))
    af_fetcher = AlphaFoldFetcher(config=config)
    af_output_path = af_fetcher.get_alphafold_model()
    print(af_output_path)
    print(f"AlphaFold model is saved to: {af_output_path}")
    dg3 = Dogsite3(config=config, pdb_file=af_output_path)
    print("Pre dogsite run")
    dg3.run()


if __name__ == "__main__":
    main()
