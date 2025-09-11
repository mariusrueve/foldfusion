# FoldFusion – Multi-Agent Orientation Guide

Purpose
FoldFusion enriches AlphaFold protein structure predictions with biologically relevant ligands/cofactors by: (1) fetching AlphaFold models, (2) predicting pockets (DoGSite3), (3) finding homologous pocket ensembles (SIENA), (4) extracting ligands from donor PDBs, (5) transplanting + optimizing (JAMDA), (6) scoring/evaluating placements.

Core Value
Produces context-enriched structures for exploratory analysis; positioned as a modular, reproducible alternative/complement to AlphaFill.

Key Entrypoints

- Pipeline orchestrator: [foldfusion/foldfusion.py](foldfusion/foldfusion/foldfusion.py) (class FoldFusion.run())
- Tools package: [foldfusion/tools/__init__.py](foldfusion/tools/__init__.py)
  - AlphaFold fetch: alphafoldfetcher.py
  - Pocket prediction: dogsite3.py
  - Pocket similarity & ensemble: siena.py
  - SIENA DB builder: siena_db.py
  - Ligand extraction: ligand_extractor.py
  - Pose optimization & scoring: jamda_scorer.py
- Configuration: [foldfusion/utils/config.py](foldfusion/foldfusion/utils/config.py)
- Dataset builder: [evaluation/data/build_benchmark_dataset.py](foldfusion/evaluation/data/build_benchmark_dataset.py)
- Analysis / visualizations: [scripts/evaluation_visualisations.py](foldfusion/scripts/evaluation_visualisations.py)
- Figures & thesis text (scientific context): [thesis/](foldfusion/thesis)

Thesis
The code is part of my master's thesis, which is written in latex in the directory `thesis/`. In the thesis I do not want to use the name `foldfusion` and keep it more abstract. Also the tone and form should be appropriate for a master's thesis and should have no "we" or "I".

High-Level Pipeline

1. Initialize config + logging
2. Build or reuse SIENA database (global)
3. For each UniProt ID:
   - Fetch AlphaFold PDB (trim low-confidence stretches: consecutive ≥5 residues pLDDT < 50)
   - Predict pocket + EDF (DoGSite3)
   - Search similar binding sites (SIENA) → ensemble PDB set
   - Extract ligands from ensemble donors
   - Pre-optimization evaluation (RMSD/TCS)
   - Optimize ligand poses (JAMDA --optimize)
   - Post-optimization evaluation
   - Persist artifacts + metadata
4. Aggregate summary + timing

External Native Tools (paths supplied in config.toml)

- dogsite3
- siena
- siena_db
- ligand_extractor
- jamda_scorer

Configuration Snapshot (minimal)

```toml
log_level = "INFO"
uniprot_ids = ["Q8CA95"]
output_dir = "results"

dogsite3_executable = "/abs/path/dogsite3"
siena_executable = "/abs/path/siena"
siena_db_executable = "/abs/path/siena_db"
ligand_extractor_executable = "/abs/path/ligand_extractor"
jamda_scorer_executable = "/abs/path/jamda_scorer"

pdb_directory = "/abs/path/pdb"   # Local PDB mirror (matching pdb_format)
pdb_format = 1                    # 0 = .ent.gz, 1 = .pdb
siena_max_alignments = 10
```

Build & Run

- Install (dev): `make install-dev` (uses editable + extras)
- Lint: `make lint`
- Format: `make format`
- Run pipeline: `python main.py` (wrapper around FoldFusion("config.toml"))
- Alt (uv): `uv run main.py`
- Background (nohup): `make nohup`

Primary Artifacts (per UniProt)
results/<uniprot_id>/
  AlphaFold/ (downloaded, trimmed model)
  Dogsite3/ (pockets + EDF)
  Siena/ (ensemble alignment outputs, resultStatistic.csv)
  LigandExtractor/ (ligand SDF/PDB groupings)
  JamdaScorer/ (optimized poses)
  evaluation.json (metrics pre/post optimization)
Global:
  logs/
  siena_database/ (if generated)
  evaluation/ (aggregated metrics)

Key Metrics (observed in evaluation scripts)

- Local ligand RMSD (pre vs post JAMDA)
- TCS (torsion clash or template consistency score—verify in code)
- Improvement deltas (% and absolute)
See: [scripts/evaluation_visualisations.py](foldfusion/scripts/evaluation_visualisations.py)

Notable Implementation Details

- AlphaFold fallback logic (v4→v3 if needed) and pLDDT-driven trimming in AlphaFoldFetcher.
- SIENA DB generated once; future optimization: skip rebuild if existing DB size > threshold (TODO present in FoldFusion.run()).
- EDF selection from DoGSite3 assumes pocket ranking consistency (select P_1); REFERENCE line normalization occurs before SIENA usage.
- Post-JAMDA evaluation currently relies on same metrics set; potential extension: energy components, clash counts.
