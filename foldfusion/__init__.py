"""
FoldFusion: A Pipeline for Enriching AlphaFold Structures with Experimental Ligands

FoldFusion is a computational pipeline designed to enhance AlphaFold protein structures
by transplanting ligands from experimental PDB structures. The pipeline combines
multiple bioinformatics tools to identify suitable binding sites, align structures,
extract ligands, and optimize their placement in the target protein structure.

Main Components:
    - AlphaFold structure retrieval and preprocessing
    - Binding site prediction using DoGSite3
    - Structure alignment using SIENA
    - Ligand extraction and transplantation
    - Ligand optimization using JAMDA
    - Comprehensive evaluation and quality assessment

Usage:
    from foldfusion import FoldFusion

    # Initialize with configuration file
    ff = FoldFusion("config.toml")

    # Run the complete pipeline
    ff.run()

For detailed configuration options and usage examples, please refer to the
project documentation and README.md file.
"""

from .foldfusion import FoldFusion
from .utils import setup_logging

__version__ = "1.0.0"
__author__ = "FoldFusion Development Team"
__all__ = ["FoldFusion", "setup_logging"]
