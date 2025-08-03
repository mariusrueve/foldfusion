"""
FoldFusion Pipeline Entry Point

This module serves as the main entry point for the FoldFusion pipeline,
which enriches AlphaFold protein structures by transplanting ligands from
experimental PDB structures.

For more information about configuration and usage, please refer to the
README.md file in the project root directory.
"""

from foldfusion import FoldFusion

# Initialize and run the FoldFusion pipeline with the default configuration
ff = FoldFusion("config.toml")
ff.run()
