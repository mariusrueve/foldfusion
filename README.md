# FoldFusion: AlphaFold Structure Enrichment Pipeline

FoldFusion is a comprehensive computational pipeline designed to enhance AlphaFold protein structures by transplanting ligands from experimental PDB structures. The pipeline combines multiple state-of-the-art bioinformatics tools to identify suitable binding sites, align structures, extract ligands, and optimize their placement in target protein structures.

## Overview

AlphaFold has revolutionized structural biology by providing high-quality protein structure predictions. However, these structures lack experimental binding partners (ligands), which are crucial for understanding protein function and drug design. FoldFusion addresses this limitation by systematically transplanting ligands from experimental structures into AlphaFold models.

## Pipeline Workflow

The FoldFusion pipeline consists of the following stages:

1. **AlphaFold Structure Retrieval**: Downloads and preprocesses protein structures from the AlphaFold Database
2. **Binding Site Prediction**: Uses DoGSite3 to identify potential ligand binding sites
3. **Structure Alignment**: Employs SIENA to find experimental structures with similar binding sites
4. **Ligand Extraction**: Extracts ligands from aligned experimental structures
5. **Ligand Optimization**: Optimizes ligand placement using JAMDA scorer
6. **Quality Assessment**: Evaluates results using Local RMSD and Transplant Clash Score metrics

## Key Features

- **Automated Processing**: Handles multiple UniProt IDs in batch mode
- **Quality Control**: Filters unreliable regions based on pLDDT scores
- **Comprehensive Evaluation**: Provides detailed quality metrics for all results
- **Robust Error Handling**: Graceful handling of failures with detailed logging
- **Configurable Parameters**: Flexible configuration through TOML files
- **Professional Logging**: Comprehensive logging with rotation and multiple output formats

## Installation

### Prerequisites

- Python 3.9 or higher
- External bioinformatics tools (see Configuration section)
- Required Python packages (automatically installed)

### Setup

```bash
# Clone the repository
git clone <repository-url>
cd foldfusion-poc

# Install using uv (recommended) or pip
uv sync  # or pip install -e .

# Verify installation
python -m foldfusion --help
```

## Configuration

FoldFusion uses TOML configuration files to specify all pipeline parameters. Create a `config.toml` file based on the provided template:

```toml
# Logging configuration
log_level = "INFO"  # DEBUG, INFO, WARNING, ERROR, CRITICAL
log_file = "foldfusion_pipeline.log"

# Input/Output settings
uniprot_ids = [
    "Q8CA95",
    "Q9QYJ6", 
    "Q9Y233"
]
output_dir = "results"

# External tool executables (adjust paths as needed)
dogsite3_executable = "/path/to/dogsite3"
siena_executable = "/path/to/siena"
siena_db_executable = "/path/to/siena_db"
ligand_extractor_executable = "/path/to/ligand_extractor"
jamda_scorer_executable = "/path/to/jamda_scorer"

# PDB database configuration
pdb_directory = "/path/to/pdb/files"
pdb_format = 1  # 0 = .ent.gz, 1 = .pdb
siena_db_database_path = "/path/to/siena.db"

# Pipeline parameters
siena_max_alignments = 10
```

### Required External Tools

The pipeline requires the following external tools to be installed and accessible:

- **DoGSite3**: Binding site prediction
- **SIENA**: Structure alignment
- **Ligand Extractor**: Ligand extraction from PDB structures
- **JAMDA Scorer**: Ligand pose optimization

## Usage

### Basic Usage

```python
from foldfusion import FoldFusion

# Initialize pipeline with configuration
pipeline = FoldFusion("config.toml")

# Run complete pipeline
pipeline.run()
```

### Command Line Usage

```bash
# Run with default configuration
python main.py

# Or run the module directly
python -m foldfusion config.toml
```

### Programmatic Usage

```python
from foldfusion import FoldFusion
from pathlib import Path

# Initialize with custom configuration
config_path = Path("custom_config.toml")
pipeline = FoldFusion(str(config_path))

# Run pipeline
try:
    pipeline.run()
    print("Pipeline completed successfully")
except Exception as e:
    print(f"Pipeline failed: {e}")
```

## Output Structure

The pipeline generates a structured output directory for each processed protein:

```plaintext
results/
├── Results/
│   ├── Q8CA95/
│   │   ├── AlphaFold/
│   │   │   ├── AF-Q8CA95-F1-model_v4.pdb
│   │   │   └── AF-Q8CA95-F1-model_v4_processed.pdb
│   │   ├── Dogsite3/
│   │   │   └── [binding site predictions]
│   │   ├── Siena/
│   │   │   └── [structure alignments]
│   │   ├── LigandExtractor/
│   │   │   └── [extracted ligands]
│   │   ├── JamdaScorer/
│   │   │   └── [optimized ligands]
│   │   └── evaluation_results.json
│   └── Q9QYJ6/
│       └── [similar structure]
├── Evaluation/
│   ├── evaluation.json
│   └── [quality assessment results]
├── SienaDB/
│   └── siena_db
└── foldfusion_pipeline.log
```

## Quality Metrics

FoldFusion provides comprehensive quality assessment using established metrics:

### Local RMSD

- Measures structural alignment quality near binding sites
- **Good**: < 0.92 Å
- **Medium**: 0.92 - 3.10 Å  
- **Poor**: > 3.10 Å

### Transplant Clash Score (TCS)

- Measures steric clashes between ligands and protein
- **Good**: < 0.64 Å
- **Medium**: 0.64 - 1.27 Å
- **Poor**: > 1.27 Å

## Logging

The pipeline provides comprehensive logging with the following features:

- **Multiple log levels**: DEBUG, INFO, WARNING, ERROR, CRITICAL
- **File rotation**: Automatic log rotation to manage disk space
- **Structured output**: Consistent formatting with timestamps and module names
- **Error tracking**: Detailed error messages with stack traces
- **Progress monitoring**: Real-time progress updates during execution

Log files are automatically rotated when they exceed 10MB, keeping up to 5 backup files.

## Error Handling

FoldFusion is designed with robust error handling:

- **Graceful degradation**: Failures in individual proteins don't stop the entire pipeline
- **Detailed error reporting**: Comprehensive error messages with context
- **Recovery mechanisms**: Automatic retry for transient failures
- **Validation checks**: Pre-execution validation of configuration and dependencies

## Contributing

We welcome contributions to FoldFusion! Please follow these guidelines:

1. Fork the repository
2. Create a feature branch
3. Implement your changes with appropriate tests
4. Ensure code follows the established style guidelines
5. Submit a pull request with a clear description

### Development Setup

```bash
# Install development dependencies
uv sync --dev

# Run tests
make test

# Run linting
make lint

# Format code
make format
```

## License

This project is licensed under the [appropriate license]. See LICENSE file for details.

## Citation

If you use FoldFusion in your research, please cite:

```plaintext
[Citation details to be added]
```

## Support

For questions, issues, or feature requests:

- Open an issue on GitHub
- Contact the development team
- Check the documentation for common solutions

## Changelog

### Version 1.0.0

- Initial release
- Complete pipeline implementation
- Comprehensive documentation
- Quality assessment metrics
- Professional logging system
