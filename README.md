# FoldFusion

FoldFusion automates ligand transplantation by combining AlphaFold protein models with donor complexes from experimental structures. The pipeline orchestrates external tools for pocket detection, structural alignment, ligand extraction, refinement, and evaluation. It is designed for structural biologists who need reproducible ligand placement and scoring on AlphaFold predictions.

## Features

- End-to-end pipeline that fetches AlphaFold models, identifies binding pockets, aligns donor complexes, transplants ligands, and refines poses.
- Modular tool wrappers for DoGSite3, SIENA, LigandExtractor, and JAMDA with consistent logging and error handling.
- Resume-aware execution that skips completed UniProt IDs and records per-stage success or failure markers.
- Configurable timeouts, memory ceilings, and concurrency options for tuning to local workstations or compute clusters.
- Structured outputs, evaluation metrics, and logs to inspect pipeline progress and quality.

## Prerequisites

- Linux or macOS with Python 3.12.
- GNU Make.
- External binaries installed and accessible: DoGSite3 (DoGSiteScorer), SIENA plus its database generator, LigandExtractor, and JAMDA scorer.
- Network access to download AlphaFold models unless cached locally.
- Hardware guidance: 4 CPU cores and 16 GB RAM recommended; ensure ~50 GB free storage for models, intermediate files, and results.

## Installation

1. Clone the repository and enter the project folder:

   ```bash
   git clone https://github.com/mariusrueve/foldfusion.git
   cd foldfusion
   ```

2. Create the virtual environment and install dependencies:

   ```bash
   make venv
   source .venv/bin/activate
   make install-dev
   ```

   Use `make install` if you only need the runtime requirements.

## Quick Start

1. Review and update `config.toml` with local paths to data directories and external tool executables.
2. Run the pipeline with GNU Make or the CLI:

   ```bash
   make run
   # or
   foldfusion config.toml
   ```

3. Find outputs under the directory specified by `output_dir` (default `analysis/data/foldfusion_output`). Logs are saved to `foldfusion_pipeline.log`.
4. Optionally rebuild reports and plots:

   ```bash
   make analysis
   ```

## Configuration

- `log_level`, `log_file`: adjust verbosity and log location.
- `uniprot_ids_file`: newline-delimited UniProt accessions to process.
- `output_dir`: base directory for per-UniProt results.
- Tool executables: absolute paths for DoGSite3, SIENA, LigandExtractor, JAMDA, plus SIENA database settings.
- `pipeline_concurrency`: number of UniProt IDs processed in parallel.
- Robustness options: set `resume` and `skip_failed_ids` to control reruns; enable optional timeouts and memory caps per tool as needed.

Refer to `config.toml` for the full list of options and defaults. The LaTeX manual in `docs/documentation.tex` provides deeper operational guidance.

## Command-Line Interface

The installed script exposes `foldfusion` as an entry point:

```bash
foldfusion /path/to/config.toml
```

The CLI reads the configuration, executes each pipeline stage in order, and writes structured output under the configured directory. Use `--help` for runtime usage details.

## Project Layout

- `foldfusion/`: core pipeline package, including orchestrator, tool wrappers, evaluation utilities, and helpers.
- `main.py`: thin executable that forwards to the package entry point.
- `analysis/`: generated datasets, notebooks, and figures.
- `docs/`: user manual and supporting documentation.
- `scripts/`: helper scripts for data preparation and automation.
- `tests/`: `pytest` suites covering core functionality.
- `Makefile`: convenience targets for environment setup, running the pipeline, linting, testing, and analysis regeneration.

## Development

- Run tests and coverage: `pytest --cov=foldfusion --cov-report=term-missing`.
- Lint and type-check: `make lint`, `mypy foldfusion/`.
- Apply Ruff and Black formatting via `make fmt` if the target exists.
- The optional `viz` dependency group installs plotting extras for analysis notebooks.

## Troubleshooting

- Download failures: confirm UniProt accessions exist and that network access is available; the fetcher falls back from AlphaFold v4 to v3 automatically.
- Missing donor ligands: inspect the generated `ligand_structure.json` to confirm ligands are present for the chosen chains.
- SIENA database errors: verify the database path exists and that the generator has write permissions; rebuild the database if PDB sources have changed.
- JAMDA crashes or memory errors: increase `jamda_memory_mb` or adjust concurrency to reduce peak usage.
- Resume markers: remove `Results/<UniProt>/_SUCCESS` or `_FAILED` under the output directory to force a rerun for individual proteins.

## License

FoldFusion is released under the MIT License. See `LICENSE` for details.
