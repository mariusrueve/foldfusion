"""
Thesis-focused analysis script for the FoldFusion project.

This version produces ONLY the plots and tables that map directly to the thesis
Results chapter, exporting figures as PDF only and LaTeX tables ready to include
in the manuscript.

Sections covered:
- 3.2 Overall transplantation performance (Table 3.1)
- 3.3 Global quality indicators (summary stats + correlations)
- 3.4 Alignment quality analysis (by sequence identity)
- 3.5 Effect of JAMDA optimization (pre vs post)
- 3.6 External Validation against AlphaFill (FF vs AF correlations)
- 3.7 Cross-Metric Relationships (TCS/LEV/Local RMSD correlations)
- Figure: Top 20 transplanted ligands (used in Results)

Inputs:
- DATA_DIR/dataset.json
- DATA_DIR/foldfusion_output/Evaluation/evaluation.json
- AlphaFill JSON paths referenced inside dataset.json (existing)

Outputs go to: results/ (figures/ for SVGs, tables/ for LaTeX)

NOTE: This script is intentionally minimal and avoids exploratory outputs
not referenced in the thesis, so you can paste figures/tables directly.
"""

import json
import math
import re
import tomllib
import warnings
from collections import defaultdict
from collections.abc import Iterable
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 - required for 3D plotting
from scipy import stats

from foldfusion.evaluation.utils import parse_pdb, parse_sdf

# ----------------- CONFIG -----------------
DATA_DIR = Path(__file__).parent / "data"
OUTPUT_DIR = Path(__file__).parent / "results"
FIGURES_DIR = OUTPUT_DIR / "figures"
TABLES_DIR = OUTPUT_DIR / "tables"

REPO_ROOT = Path(__file__).resolve().parents[1]
THESIS_ROOT = REPO_ROOT / "thesis"
THESIS_FIGURES_DIR = THESIS_ROOT / "figures" / "03_results"
THESIS_TABLES_DIR = THESIS_ROOT / "tables" / "03_results"
THESIS_APPENDIX_TABLES_DIR = THESIS_ROOT / "tables" / "05_appendix"

FIGURE_OUTPUT_DIRS = [FIGURES_DIR, THESIS_FIGURES_DIR]
TABLE_OUTPUT_DIRS = [TABLES_DIR, THESIS_TABLES_DIR]
APPENDIX_TABLE_OUTPUT_DIRS = [TABLES_DIR, THESIS_APPENDIX_TABLES_DIR]

DATASET_JSON = DATA_DIR / "dataset.json"
FOLDFUSION_JSON = DATA_DIR / "foldfusion_output/Evaluation/evaluation.json"
PIPELINE_LOG = DATA_DIR / "foldfusion_output/foldfusion_pipeline.log"
CONFIG_PATH = Path("config.toml")

RE_RESULTS_TARGET = re.compile(r"Results/([A-Za-z0-9]+)")
RE_AF_TARGET = re.compile(r"AF-([A-Za-z0-9]+)-F1")
RE_FAILED_TARGET = re.compile(r"Failed to process UniProt ID ([A-Za-z0-9]+):")

BUFFER_LIGANDS = {
    "GOL",
    "EDO",
    "PEG",
    "PE3",
    "PE4",
    "PE8",
    "MPD",
}

METRIC_LABELS = {
    "ff_local_rmsd": ("Local RMSD", "Å"),
    "ff_tcs_score": ("TCS", "dimensionless"),
    "ff_lev_score": ("LEV", "Å"),
    "tcs_improvement": ("ΔTCS (pre - post)", "dimensionless"),
    "rmsd_change": ("ΔLocal RMSD (post - pre)", "Å"),
}

METAL_LIGANDS = {
    "ZN",
    "MG",
    "CA",
    "MN",
    "CU",
    "FE",
    "NA",
    "K",
    "CO",
    "NI",
}

LOCAL_RMSD_THRESHOLD = 1.0
TCS_THRESHOLD = 0.35
LEV_THRESHOLD = 3.0


def _is_buffer_ligand(name: str | None) -> bool:
    if not name:
        return False
    if name in BUFFER_LIGANDS:
        return True
    if name.startswith("PEG"):
        return True
    if len(name) == 3 and name.startswith("PE") and name[2].isdigit():
        return True
    return False


def _format_metric_label(metric_key: str) -> str:
    pretty, unit = METRIC_LABELS.get(metric_key, (metric_key, None))
    if unit is None:
        return pretty
    if unit == "dimensionless":
        return f"{pretty} (dimensionless)"
    return f"{pretty} ({unit})"


# Create output dirs
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
for path in [
    OUTPUT_DIR,
    FIGURES_DIR,
    TABLES_DIR,
    THESIS_FIGURES_DIR,
    THESIS_TABLES_DIR,
    THESIS_APPENDIX_TABLES_DIR,
]:
    path.mkdir(exist_ok=True, parents=True)

# Style
sns.set_theme(context="paper", style="whitegrid", palette="colorblind", font_scale=1.75)

warnings.filterwarnings("ignore", message="invalid value encountered")
pd.options.mode.use_inf_as_na = True

# ----------------- Helpers -----------------


def _save_figure(fig: Figure, stem: str):
    """Save a Matplotlib figure as PDF into all configured output directories."""

    try:
        fig.tight_layout(pad=0.2)
    except Exception:
        pass

    for directory in FIGURE_OUTPUT_DIRS:
        directory.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            (directory / stem).with_suffix(".pdf"),
            bbox_inches="tight",
            pad_inches=0,
        )


def _to_latex(
    df: pd.DataFrame,
    path: Path,
    caption: str,
    label: str,
    short_caption: str | None = None,
):
    """Export a DataFrame as a self-contained LaTeX table (tabular inside table).

    Also ensures LaTeX-safe and human-friendly column headers. For common
    snake_case names, map to readable labels (e.g., x_metric -> X-Axis Metric).
    """

    def _pretty_col(col: str) -> str:
        mapping = {
            "x_metric": "X-Axis Metric",
            "y_metric": "Y-Axis Metric",
            "metric_x": "Metric X",
            "metric_y": "Metric Y",
            "pearson_r": "Pearson r",
            "spearman_rho": "Spearman rho",
            "p_value": "p-value",
            "n": "N",
            "n_pairs": "N pairs",
            "Mean (95% CI)": r"Mean (95\% CI)",
            "Median (95% CI)": r"Median (95\% CI)",
            "Pearson (95% CI)": r"Pearson (95\% CI)",
            "Spearman (95% CI)": r"Spearman (95\% CI)",
            "Median TCS (dimensionless)": "Median TCS",
            "Median Local RMSD (\\si{\\angstrom})": r"Median Local RMSD (\si{\angstrom})",
            "Share vs previous (%)": r"Share vs previous (\%)",
            "Share Vs Previous (%)": r"Share vs previous (\%)",
            "Median per target (count)": "Median per target (count)",
            "Median Per Target (Count)": "Median per target (count)",
            "Median wall time (s)": "Median wall time (s)",
            "Median Wall Time (S)": "Median wall time (s)",
            "Memory ceiling (MB)": "Memory ceiling (MB)",
            "Memory Ceiling (Mb)": "Memory ceiling (MB)",
            "Share of placements (%)": r"Share of placements (\%)",
            "Share of targets (%)": r"Share of targets (\%)",
            "TCS pass rate (%)": r"TCS pass rate (\%)",
            "Joint pass rate (%)": r"Joint pass rate (\%)",
            "Min. Identity (Dimensionless)": "Min. identity",
            "Min. identity (dimensionless)": "Min. identity",
            "Share Of Total (%)": r"Share of total (\%)",
            "Share of total (%)": r"Share of total (\%)",
            "Success Rate (%)": r"Success rate (\%)",
            "Success rate (%)": r"Success rate (\%)",
        }
        if col in mapping:
            return mapping[col]
        # Generic fallback: title case and replace underscores with spaces
        return col.replace("_", " ").title()

    path.parent.mkdir(parents=True, exist_ok=True)
    df_latex = df.copy()
    # Round numeric columns to 3 decimals for consistent presentation
    num_cols = df_latex.select_dtypes(include=[np.number]).columns
    if len(num_cols) > 0:
        df_latex[num_cols] = df_latex[num_cols].round(3)
    df_latex.columns = [_pretty_col(c) for c in df_latex.columns]
    table_body = df_latex.to_latex(
        index=False,
        escape=False,
        na_rep="--",
        float_format=lambda x: f"{x:.3f}",
    ).strip()

    with open(path, "w") as f:
        f.write(
            "\n".join(
                [
                    "\\begin{table}[ht]",
                    "\\centering",
                    "\\resizebox{\\linewidth}{!}{%",
                    table_body,
                    "}",
                    f"\\caption{('[' + short_caption + ']' if short_caption else '')}{{{caption}}}",
                    f"\\label{{{label}}}",
                    "\\end{table}\n",
                ]
            )
        )


def _export_table(
    df: pd.DataFrame,
    filename: str,
    caption: str,
    label: str,
    short_caption: str | None = None,
    directories: Iterable[Path] | None = None,
):
    """Write LaTeX tables to both local analysis and thesis directories."""

    targets = list(directories) if directories is not None else TABLE_OUTPUT_DIRS
    for directory in targets:
        directory.mkdir(parents=True, exist_ok=True)
        _to_latex(df, directory / filename, caption, label, short_caption)


def _series_stats(x: pd.Series, q=(0.25, 0.5, 0.75)):
    x = pd.to_numeric(x, errors="coerce").dropna()
    res = {
        "n": int(x.size),
        "mean": x.mean(),
        "std": x.std(ddof=1) if x.size > 1 else np.nan,
        "min": x.min() if x.size else np.nan,
        "q25": x.quantile(q[0]) if x.size else np.nan,
        "median": x.quantile(q[1]) if x.size else np.nan,
        "q75": x.quantile(q[2]) if x.size else np.nan,
        "max": x.max() if x.size else np.nan,
    }
    return res


def _bootstrap_ci(
    series: pd.Series,
    func,
    n_boot: int = 2000,
    ci: float = 0.95,
    random_state: int = 42,
):
    """Return (lower, upper) bootstrap CI for a statistic."""

    data = pd.to_numeric(series, errors="coerce").dropna().to_numpy()
    if data.size == 0:
        return (np.nan, np.nan)
    rng = np.random.default_rng(random_state)
    stats_boot = []
    for _ in range(n_boot):
        sample = rng.choice(data, size=data.size, replace=True)
        stats_boot.append(func(sample))
    alpha = (1 - ci) / 2
    lower = np.quantile(stats_boot, alpha)
    upper = np.quantile(stats_boot, 1 - alpha)
    return lower, upper


def _fisher_ci(r: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Analytical 95% CI for correlation using Fisher z-transform."""

    if n <= 3 or np.isnan(r):
        return (np.nan, np.nan)
    r = max(min(r, 0.999999), -0.999999)
    z = np.arctanh(r)
    se = 1 / math.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)
    lo = np.tanh(z - z_crit * se)
    hi = np.tanh(z + z_crit * se)
    return lo, hi


def _mean_ci(series: pd.Series, alpha: float = 0.05) -> tuple[float, float]:
    arr = pd.to_numeric(series, errors="coerce").dropna().to_numpy()
    if arr.size < 2:
        return (np.nan, np.nan)
    mean = arr.mean()
    se = arr.std(ddof=1) / math.sqrt(arr.size)
    tcrit = stats.t.ppf(1 - alpha / 2, arr.size - 1)
    return mean - tcrit * se, mean + tcrit * se


def _wilson_ci(successes: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    z = stats.norm.ppf(1 - alpha / 2)
    phat = successes / n
    denominator = 1 + (z**2) / n
    center = (phat + (z**2) / (2 * n)) / denominator
    margin = z * math.sqrt((phat * (1 - phat) / n) + (z**2) / (4 * n**2)) / denominator
    return max(0.0, center - margin), min(1.0, center + margin)


def _format_ci_value(value: float, low: float, high: float, precision: int = 3) -> str:
    if any(np.isnan(v) for v in (value, low, high)):
        return "--"
    return f"{value:.{precision}f} [ {low:.{precision}f}, {high:.{precision}f} ]"


def _metric_precision(metric: str) -> int:
    metric_lower = metric.lower()
    if "lev" in metric_lower:
        return 2
    return 3


def _format_metric_scalar(value: float, metric: str) -> str:
    if pd.isna(value):
        return "--"
    precision = _metric_precision(metric)
    return f"{value:.{precision}f}"


def load_config_metadata(config_path: Path = CONFIG_PATH) -> dict:
    """Return configuration-derived metadata (IDs, memory ceilings)."""

    meta = {
        "uniprot_ids": [],
        "memory_limits": {},
    }
    if not config_path.exists():
        return meta

    try:
        with config_path.open("rb") as fh:
            cfg = tomllib.load(fh)
    except Exception:
        return meta

    if isinstance(cfg.get("uniprot_ids"), list):
        meta["uniprot_ids"] = [
            str(x).strip() for x in cfg["uniprot_ids"] if str(x).strip()
        ]
    elif cfg.get("uniprot_ids_file"):
        try:
            path = Path(cfg["uniprot_ids_file"]).expanduser()
            if path.exists():
                meta["uniprot_ids"] = [
                    line.strip()
                    for line in path.read_text().splitlines()
                    if line.strip()
                ]
        except Exception:
            pass

    for key in [
        "dogsite3_memory_mb",
        "siena_memory_mb",
        "jamda_memory_mb",
        "siena_db_memory_mb",
        "ligand_extractor_memory_mb",
    ]:
        value = cfg.get(key)
        if value is not None:
            try:
                meta["memory_limits"][key] = int(value)
            except (TypeError, ValueError):
                continue

    return meta


def _guess_target(message: str, fallback: str | None = None) -> str | None:
    for pattern in (RE_RESULTS_TARGET, RE_AF_TARGET):
        m = pattern.search(message)
        if m:
            return m.group(1)
    return fallback


def parse_pipeline_log(log_path: Path = PIPELINE_LOG) -> pd.DataFrame:
    """Parse the pipeline log to extract per-target stage durations."""

    if not log_path.exists():
        return pd.DataFrame()

    records: dict[str, dict[str, object]] = {}
    stage_starts: dict[str, dict[str, datetime]] = defaultdict(dict)
    current_target: str | None = None

    with log_path.open() as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split(" - ", 3)
            if len(parts) < 4:
                continue
            ts_str, _, _, message = parts
            try:
                ts = datetime.strptime(ts_str, "%Y-%m-%d %H:%M:%S,%f")
            except ValueError:
                continue

            if "Starting pipeline for UniProt ID:" in message:
                target = message.split(":")[-1].strip()
                current_target = target
                rec = records.setdefault(target, {"uniprot_id": target})
                rec["pipeline_start"] = ts
                stage_starts[target] = {}
                continue

            if "Completed pipeline for" in message:
                target = message.split("for")[-1].strip()
                rec = records.setdefault(target, {"uniprot_id": target})
                rec["pipeline_end"] = ts
                rec["pipeline_status"] = "completed"
                current_target = None
                continue

            m_failed = RE_FAILED_TARGET.search(message)
            if m_failed:
                target = m_failed.group(1)
                rec = records.setdefault(target, {"uniprot_id": target})
                rec["pipeline_end"] = ts
                rec["pipeline_status"] = "failed"
                current_target = None
                continue

            target = _guess_target(message, current_target)
            if target is None:
                continue
            if current_target is None:
                current_target = target

            rec = records.setdefault(target, {"uniprot_id": target})
            stages = stage_starts[target]

            if "Starting DoGSite3 binding site prediction" in message:
                stages["dogsite3"] = ts
            elif "DoGSite3 analysis completed successfully" in message:
                start = stages.pop("dogsite3", None)
                if start:
                    rec["dogsite3_seconds"] = (ts - start).total_seconds()
                rec["dogsite3_success"] = True
            elif "DoGSite3 analysis failed" in message:
                start = stages.pop("dogsite3", None)
                if start:
                    rec["dogsite3_seconds"] = (ts - start).total_seconds()
                rec["dogsite3_success"] = False

            elif "Starting SIENA structure alignment" in message:
                stages["siena"] = ts
            elif "SIENA alignment completed successfully" in message:
                start = stages.pop("siena", None)
                if start:
                    rec["siena_seconds"] = (ts - start).total_seconds()
                rec["siena_success"] = True
            elif "SIENA alignment failed" in message:
                start = stages.pop("siena", None)
                if start:
                    rec["siena_seconds"] = (ts - start).total_seconds()
                rec["siena_success"] = False

            elif "Starting ligand extraction" in message:
                stages["ligand"] = ts
            elif "Ligand extraction completed" in message:
                start = stages.pop("ligand", None)
                if start:
                    rec["ligand_seconds"] = (ts - start).total_seconds()

            elif "Starting JAMDA optimization" in message:
                stages["jamda"] = ts
                rec.pop("jamda_success", None)
            elif "JAMDA optimization completed" in message:
                start = stages.pop("jamda", None)
                if start:
                    rec["jamda_seconds"] = (ts - start).total_seconds()
                rec["jamda_success"] = True
            elif "No ligands available for JAMDA optimization" in message:
                stages.pop("jamda", None)
                rec.setdefault("jamda_seconds", 0.0)
                rec["jamda_success"] = False

            elif "Evaluation 'pre-jamda'" in message and "completed" in message:
                start = stages.pop("eval_pre", None)
                if start:
                    rec["eval_pre_seconds"] = (ts - start).total_seconds()
            elif "Starting evaluation stage: pre-jamda" in message:
                stages["eval_pre"] = ts
            elif "Evaluation 'post-jamda'" in message and "completed" in message:
                start = stages.pop("eval_post", None)
                if start:
                    rec["eval_post_seconds"] = (ts - start).total_seconds()
            elif "Starting evaluation stage: post-jamda" in message:
                stages["eval_post"] = ts

            elif "Initiating AlphaFold model retrieval" in message:
                stages["alphafold"] = ts
            elif "AlphaFold model successfully retrieved" in message:
                start = stages.pop("alphafold", None)
                if start:
                    rec["alphafold_seconds"] = (ts - start).total_seconds()

    rows = []
    for target, data in records.items():
        row = data.copy()
        start = row.get("pipeline_start")
        end = row.get("pipeline_end")
        if isinstance(start, datetime) and isinstance(end, datetime):
            row["pipeline_seconds"] = (end - start).total_seconds()
        for key in list(row.keys()):
            if isinstance(row[key], datetime):
                dt_obj = row[key]
                assert isinstance(dt_obj, datetime)
                row[key] = dt_obj.isoformat()
        rows.append(row)

    return pd.DataFrame(rows)


def _read_sdf_with_bonds(sdf_path: Path):
    try:
        lines = sdf_path.read_text().splitlines()
    except FileNotFoundError:
        return None, []
    if len(lines) < 4:
        return None, []
    try:
        num_atoms = int(lines[3][0:3])
        num_bonds = int(lines[3][3:6])
    except ValueError:
        return None, []

    bonds = []
    for i in range(4 + num_atoms, 4 + num_atoms + num_bonds):
        if i >= len(lines):
            break
        line = lines[i]
        if len(line) < 6:
            continue
        try:
            a1 = int(line[0:3]) - 1
            a2 = int(line[3:6]) - 1
        except ValueError:
            continue
        bonds.append((a1, a2))
    atoms = parse_sdf(sdf_path)["atoms"] if sdf_path.exists() else []
    return atoms, bonds


def _compute_neighbor_residues(
    pdb_path: Path, ligand_atoms: list[dict], cutoff: float = 4.0, limit: int = 8
):
    try:
        protein = parse_pdb(pdb_path)
    except FileNotFoundError:
        return []
    ligand_coords = np.array(
        [[atom["x"], atom["y"], atom["z"]] for atom in ligand_atoms], dtype=float
    )
    residues: dict[tuple[str, int, str], list[dict]] = defaultdict(list)
    for atom in protein.get("atoms", []):
        key = (
            atom.get("chain_id") or "A",
            atom.get("residue_number"),
            atom.get("residue_name"),
        )
        residues[key].append(atom)

    neighbors = []
    for (chain, resno, resname), atom_list in residues.items():
        coords = np.array([[a["x"], a["y"], a["z"]] for a in atom_list], dtype=float)
        if coords.size == 0:
            continue
        distances = np.linalg.norm(
            coords[:, None, :] - ligand_coords[None, :, :], axis=2
        )
        min_dist = float(distances.min())
        if min_dist > cutoff:
            continue
        ca_coords = [
            (a["x"], a["y"], a["z"]) for a in atom_list if a.get("atom_name") == "CA"
        ]
        if ca_coords:
            centroid = np.array(ca_coords[0], dtype=float)
        else:
            centroid = coords.mean(axis=0)
        neighbors.append(
            {
                "chain": chain or "A",
                "resno": resno,
                "resname": resname,
                "distance": min_dist,
                "centroid": centroid,
            }
        )

    neighbors.sort(key=lambda x: x["distance"])
    return neighbors[:limit]


def _set_equal_aspect(ax, points: np.ndarray):
    if points.size == 0:
        return
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    centers = (mins + maxs) / 2
    max_range = (maxs - mins).max() / 2 or 1.0
    ax.set_xlim(centers[0] - max_range, centers[0] + max_range)
    ax.set_ylim(centers[1] - max_range, centers[1] + max_range)
    ax.set_zlim(centers[2] - max_range, centers[2] + max_range)


# ----------------- Data loading -----------------


def load_foldfusion_data(filepath: Path) -> pd.DataFrame:
    if not filepath.exists():
        print(f"[WARN] FoldFusion results not found at: {filepath}")
        return pd.DataFrame()
    with open(filepath) as f:
        data = json.load(f)

    records = []
    for uniprot_id, pdb_data in data.items():
        for pdb_code, ligand_data in pdb_data.items():
            active_site_identity = ligand_data.get("active_site_identity")
            for ligand_id, stage_data in ligand_data.items():
                if not isinstance(stage_data, dict):  # skip meta keys
                    continue
                for stage, metrics in stage_data.items():
                    if not isinstance(metrics, dict):
                        continue
                    tcs_data = metrics.get("tcs")
                    lev_data = metrics.get("lev")
                    records.append(
                        {
                            "uniprot_id": uniprot_id,
                            "pdb_code": pdb_code,
                            "ligand_id": ligand_id,
                            "stage": stage,
                            "active_site_identity": active_site_identity,
                            "ff_local_rmsd": metrics.get("local_rmsd"),
                            "ff_tcs_score": tcs_data.get("score")
                            if isinstance(tcs_data, dict)
                            else None,
                            "ff_tcs_clash_count": tcs_data.get("clash_count")
                            if isinstance(tcs_data, dict)
                            else None,
                            "ff_tcs_atoms": tcs_data.get("transplant_atom_count")
                            if isinstance(tcs_data, dict)
                            else None,
                            "ff_tcs_poly_atoms": tcs_data.get("poly_atom_count")
                            if isinstance(tcs_data, dict)
                            else None,
                            "ff_lev_score": lev_data.get("local_environment_rmsd")
                            if isinstance(lev_data, dict)
                            else None,
                        }
                    )
    return pd.DataFrame(records)


def load_alphafill_data(dataset_filepath: Path) -> pd.DataFrame:
    if not dataset_filepath.exists():
        print(f"[WARN] Dataset JSON not found at: {dataset_filepath}")
        return pd.DataFrame()
    with open(dataset_filepath) as f:
        dataset = json.load(f)

    records = []
    for ligand_bin in dataset.get("ligand_bins", []):
        for row in ligand_bin.get("entries", []):
            uniprot_id = row.get("uniprot_id")
            alphafill_json_path = row.get("alphafill_json")
            if not uniprot_id or not alphafill_json_path:
                continue
            af_path = Path(alphafill_json_path)
            if not af_path.exists():
                continue
            try:
                with open(af_path) as f:
                    af_data = json.load(f)
            except Exception:
                continue
            for hit in af_data.get("hits") or []:
                for transplant in (hit or {}).get("transplants") or []:
                    records.append(
                        {
                            "uniprot_id": uniprot_id,
                            "pdb_code": (hit or {}).get("pdb_id"),
                            "ligand_id": (transplant or {}).get("compound_id"),
                            "af_local_rmsd": (transplant or {}).get("local_rmsd"),
                            "af_tcs_score": (transplant or {})
                            .get("clash", {})
                            .get("score"),
                            "af_lev_score": (
                                (transplant or {}).get("validation") or {}
                            ).get("local_environment_rmsd"),
                        }
                    )
    return pd.DataFrame(records)


def load_bins(dataset_filepath: Path) -> pd.DataFrame:
    if not dataset_filepath.exists():
        return pd.DataFrame()
    with open(dataset_filepath) as f:
        dataset = json.load(f)

    u2l, u2p = {}, {}
    for ligand_bin in dataset.get("ligand_bins", []):
        name = ligand_bin.get("ligand_bin")
        for entry in ligand_bin.get("entries", []):
            if entry.get("uniprot_id"):
                u2l[entry["uniprot_id"]] = name
    for protein_bin in dataset.get("protein_bins", []):
        cls = protein_bin.get("protein_class")
        for entry in protein_bin.get("entries", []):
            if entry.get("uniprot_id"):
                u2p[entry["uniprot_id"]] = cls
    all_u = set(u2l) | set(u2p)
    rows = [
        {
            "uniprot_id": u,
            "ligand_bin": u2l.get(u),
            "protein_class": u2p.get(u),
        }
        for u in all_u
    ]
    return pd.DataFrame(rows)


def load_expected_pairs(dataset_filepath: Path) -> pd.DataFrame:
    if not dataset_filepath.exists():
        return pd.DataFrame()

    with open(dataset_filepath) as f:
        dataset = json.load(f)

    records: list[dict[str, str]] = []
    for row in dataset.get("rows", []):
        uniprot_id = (row.get("uniprot_id") or "").strip()
        validation_pdb = (row.get("validation_pdb") or "").strip()
        if not uniprot_id or not validation_pdb:
            continue

        ligands = row.get("matched_ligands") or []
        if not ligands:
            raw = row.get("matched_ligand") or ""
            ligands = [part.strip() for part in raw.split("|") if part.strip()]

        for ligand in ligands:
            ligand_code = ligand.strip()
            if not ligand_code:
                continue
            records.append(
                {
                    "uniprot_id": uniprot_id,
                    "expected_pdb": validation_pdb.upper(),
                    "expected_ligand": ligand_code.upper(),
                    "bin_type": row.get("bin_type"),
                    "bin_name": row.get("bin_name"),
                }
            )

    if not records:
        return pd.DataFrame()

    expected_df = pd.DataFrame(records).drop_duplicates()
    return expected_df


def build_master() -> pd.DataFrame:
    print("[INFO] Loading data...")
    df_ff = load_foldfusion_data(FOLDFUSION_JSON)
    df_af = load_alphafill_data(DATASET_JSON)
    df_bins = load_bins(DATASET_JSON)

    if df_ff.empty:
        print("[ERROR] No FoldFusion data found; aborting.")
        return pd.DataFrame()

    print("[INFO] Merging...")
    df = pd.merge(df_ff, df_bins, on="uniprot_id", how="left")

    df["ligand_id_norm"] = df["ligand_id"].astype(str).str.split("_").str[0]
    if not df_af.empty:
        df_af["ligand_id_norm"] = df_af["ligand_id"].astype(str).str.strip()
        df = pd.merge(
            df,
            df_af,
            on=["uniprot_id", "pdb_code", "ligand_id_norm"],
            how="left",
            suffixes=("_ff", "_af"),
        )
        if "ligand_id_af" in df.columns:
            df = df.drop(columns=["ligand_id_af"])
        df = df.rename(columns={"ligand_id_ff": "ligand_id"})

    # Sanitize numerics
    for col in [
        "ff_local_rmsd",
        "ff_tcs_score",
        "ff_lev_score",
        "af_local_rmsd",
        "af_tcs_score",
        "af_lev_score",
        "active_site_identity",
        "ff_tcs_clash_count",
        "ff_tcs_atoms",
        "ff_tcs_poly_atoms",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df.replace([np.inf, -np.inf], np.nan, inplace=True)

    df["ligand_name"] = df["ligand_id"].astype(str).str.split("_").str[0]
    df["is_ion"] = df["ff_tcs_atoms"].fillna(0).astype(int).eq(1)
    df["ligand_class"] = np.where(df["is_ion"], "ion", "organic")

    # Deduplicate FoldFusion rows to avoid one-to-many expansion when merging with AlphaFill
    dedup_cols = [
        "uniprot_id",
        "pdb_code",
        "ligand_id",
        "stage",
        "ff_local_rmsd",
        "ff_tcs_score",
        "ff_lev_score",
        "ff_tcs_clash_count",
        "ff_tcs_atoms",
        "ff_tcs_poly_atoms",
        "active_site_identity",
    ]
    existing_cols = [c for c in dedup_cols if c in df.columns]
    if existing_cols:
        before = len(df)
        df = df.drop_duplicates(subset=existing_cols)
        after = len(df)
        if after < before:
            print(
                f"[INFO] Removed {before - after} duplicate FoldFusion rows before downstream merges."
            )

    df.to_csv(OUTPUT_DIR / "master_dataframe.csv", index=False)
    print(f"[INFO] Master dataframe saved: {OUTPUT_DIR / 'master_dataframe.csv'}")
    return df


def load_donor_hits(results_root: Path) -> set[tuple[str, str]]:
    donors: set[tuple[str, str]] = set()
    if not results_root.exists():
        return donors

    for entry in results_root.iterdir():
        if not entry.is_dir():
            continue
        siena_dir = entry / "Siena" / "best_alignments.json"
        if not siena_dir.exists():
            continue
        try:
            with siena_dir.open() as fh:
                alignments = json.load(fh)
        except Exception:
            continue

        if isinstance(alignments, dict):
            iterable = alignments.values()
        else:
            iterable = alignments or []

        for aln in iterable:
            if not isinstance(aln, dict):
                continue
            pdb_code = aln.get("pdb_code") or aln.get("pdb")
            if not pdb_code:
                continue
            donors.add((entry.name, str(pdb_code).upper()))

    return donors


# ----------------- Thesis outputs -----------------


def table_overall_summary(df: pd.DataFrame):
    """Table 3.1: Overall transplantation summary (post-JAMDA)."""
    summary = {
        "Targets": df["uniprot_id"].nunique(),
        "Donors": df["pdb_code"].nunique(),
        "Ligands considered": int((df["stage"] == "pre-jamda").sum()),
        "Post-refinement placements": int((df["stage"] == "post-jamda").sum()),
    }
    out = pd.DataFrame([summary])
    out.to_csv(OUTPUT_DIR / "table_overall_summary.csv", index=False)
    _export_table(
        out,
        "table_overall_summary.tex",
        caption="Overall transplantation summary (post-JAMDA).",
        label="tab:overall-summary",
        short_caption="Overall transplantation summary",
    )
    print("[OK] Table 3.1 exported (CSV + LaTeX).")


def table_attrition_funnel(df: pd.DataFrame, stage_df: pd.DataFrame, config_meta: dict):
    """Dataset transparency table showing attrition and runtime summaries."""

    if stage_df.empty:
        print("[WARN] Pipeline log unavailable; skipping attrition funnel table.")
        return

    targets_configured = len(config_meta.get("uniprot_ids", []))

    results_root = DATA_DIR / "foldfusion_output" / "Results"
    target_dirs = [p for p in results_root.iterdir() if p.is_dir()]
    targets_attempted = len(target_dirs)

    def _has_nonempty_best_alignments(path: Path) -> bool:
        try:
            with path.open() as fh:
                data = json.load(fh)
            return bool(data)
        except Exception:
            return False

    donor_dirs = [
        d
        for d in target_dirs
        if _has_nonempty_best_alignments(d / "Siena" / "best_alignments.json")
    ]
    donors = len(donor_dirs)
    pipeline_completed = sum(1 for d in donor_dirs if (d / "_SUCCESS").exists())

    if not targets_configured:
        targets_configured = targets_attempted

    pre = df[df["stage"] == "pre-jamda"].copy()
    post = df[df["stage"] == "post-jamda"].copy()
    pre_count = int(pre.shape[0])
    post_count = int(post.shape[0])
    pre_median_per_target = (
        float(pre.groupby("uniprot_id").size().median()) if not pre.empty else np.nan
    )
    post_median_per_target = (
        float(post.groupby("uniprot_id").size().median()) if not post.empty else np.nan
    )

    def _median_seconds(column: str) -> float:
        if column not in stage_df:
            return np.nan
        series = pd.to_numeric(stage_df[column], errors="coerce").dropna()
        return float(series.median()) if not series.empty else np.nan

    siena_median = _median_seconds("siena_seconds")
    jamda_median = _median_seconds("jamda_seconds")
    pipeline_median = _median_seconds("pipeline_seconds")
    ligand_median = _median_seconds("ligand_seconds")

    mem = config_meta.get("memory_limits", {})

    rows = []

    def _share(numerator: int, denominator: int | float | None) -> float | None:
        if denominator in (None, 0):
            return None
        return numerator / denominator

    rows.append(
        {
            "Stage": "Targets configured",
            "Unit": "targets",
            "Count": targets_configured,
            "Share vs previous": None,
            "Median per-target": None,
            "Median wall time (s)": None,
            "Memory ceiling (MB)": None,
        }
    )
    rows.append(
        {
            "Stage": "Targets attempted (pipeline started)",
            "Unit": "targets",
            "Count": targets_attempted,
            "Share vs previous": _share(targets_attempted, targets_configured),
            "Median per-target": None,
            "Median wall time (s)": pipeline_median,
            "Memory ceiling (MB)": None,
        }
    )
    rows.append(
        {
            "Stage": "Donor alignments found (SIENA)",
            "Unit": "targets",
            "Count": donors,
            "Share vs previous": _share(donors, targets_attempted),
            "Median per-target": None,
            "Median wall time (s)": siena_median,
            "Memory ceiling (MB)": mem.get("siena_memory_mb"),
        }
    )
    rows.append(
        {
            "Stage": "Pipelines completed",
            "Unit": "targets",
            "Count": pipeline_completed,
            "Share vs previous": _share(pipeline_completed, donors),
            "Median per-target": None,
            "Median wall time (s)": pipeline_median,
            "Memory ceiling (MB)": None,
        }
    )
    rows.append(
        {
            "Stage": "Placements before JAMDA",
            "Unit": "placements",
            "Count": pre_count,
            "Share vs previous": None,
            "Median per-target": pre_median_per_target,
            "Median wall time (s)": ligand_median,
            "Memory ceiling (MB)": mem.get("ligand_extractor_memory_mb"),
        }
    )
    rows.append(
        {
            "Stage": "Placements after JAMDA",
            "Unit": "placements",
            "Count": post_count,
            "Share vs previous": _share(post_count, pre_count),
            "Median per-target": post_median_per_target,
            "Median wall time (s)": jamda_median,
            "Memory ceiling (MB)": mem.get("jamda_memory_mb"),
        }
    )

    attrition_df = pd.DataFrame(rows)
    attrition_df.to_csv(OUTPUT_DIR / "table_attrition_funnel.csv", index=False)

    latex = attrition_df.copy()
    latex["Share vs previous"] = latex["Share vs previous"].map(
        lambda x: f"{x * 100:.1f}\\%" if pd.notnull(x) else "--"
    )
    latex["Median per-target"] = latex["Median per-target"].map(
        lambda x: f"{x:.1f}" if pd.notnull(x) else "--"
    )

    def _format_wall_time(x: float) -> str:
        if pd.isnull(x):
            return "--"
        if abs(x) < 1:
            return f"{x:.2f}".rstrip("0").rstrip(".")
        return f"{x:.1f}"

    latex["Median wall time (s)"] = latex["Median wall time (s)"].map(_format_wall_time)
    latex["Memory ceiling (MB)"] = latex["Memory ceiling (MB)"].map(
        lambda x: f"{int(x)}" if pd.notnull(x) else "--"
    )
    latex = latex.rename(
        columns={
            "Share vs previous": "Share vs previous (%)",
            "Median per-target": "Median per target (count)",
            "Median wall time (s)": "Median wall time (s)",
            "Memory ceiling (MB)": "Memory ceiling (MB)",
        }
    )

    _export_table(
        latex,
        "table_attrition_funnel.tex",
        caption=(
            "Pipeline inclusion funnel showing target attrition, placement counts, "
            "and median wall times per stage. Shares are relative to the preceding "
            "stage with matching units."
        ),
        label="tab:attrition-funnel",
        short_caption="Pipeline attrition summary",
    )
    print("[OK] Attrition funnel table exported.")


def table_expected_recovery(df: pd.DataFrame):
    """Compare benchmark expectations against recovered donors and ligands."""

    expected_pairs = load_expected_pairs(DATASET_JSON)
    if expected_pairs.empty:
        print(
            "[WARN] Benchmark manifest unavailable; skipping expected recovery table."
        )
        return

    expected_pairs = expected_pairs.copy()
    expected_pairs["expected_pdb"] = (
        expected_pairs["expected_pdb"].astype(str).str.upper()
    )
    expected_pairs["expected_ligand"] = (
        expected_pairs["expected_ligand"].astype(str).str.upper()
    )

    results_root = DATA_DIR / "foldfusion_output" / "Results"
    donor_hits = load_donor_hits(results_root)
    expected_pairs["pdb_found"] = expected_pairs.apply(
        lambda row: (row["uniprot_id"], row["expected_pdb"]) in donor_hits,
        axis=1,
    )

    if "ligand_name" not in df.columns:
        df = df.copy()
        df["ligand_name"] = df["ligand_id"].astype(str).str.split("_").str[0]

    placements = df.dropna(subset=["pdb_code", "ligand_name", "stage"]).copy()
    placements["pdb_code"] = placements["pdb_code"].astype(str).str.upper()
    placements["ligand_name"] = placements["ligand_name"].astype(str).str.upper()
    placements["stage"] = placements["stage"].astype(str).str.lower()

    pre_pairs = set()
    post_pairs = set()
    if not placements.empty:
        pre_pairs = {
            (str(row.uniprot_id), str(row.pdb_code), str(row.ligand_name))
            for row in placements.loc[
                placements["stage"] == "pre-jamda",
                ["uniprot_id", "pdb_code", "ligand_name"],
            ]
            .drop_duplicates()
            .itertuples(index=False)
        }
        post_pairs = {
            (str(row.uniprot_id), str(row.pdb_code), str(row.ligand_name))
            for row in placements.loc[
                placements["stage"] == "post-jamda",
                ["uniprot_id", "pdb_code", "ligand_name"],
            ]
            .drop_duplicates()
            .itertuples(index=False)
        }
    any_pairs = pre_pairs | post_pairs

    def _pair_in(pair_set: set[tuple[str, str, str]], row: pd.Series) -> bool:
        return (
            str(row["uniprot_id"]),
            str(row["expected_pdb"]),
            str(row["expected_ligand"]),
        ) in pair_set

    expected_pairs["ligand_found_pre"] = expected_pairs.apply(
        lambda row: _pair_in(pre_pairs, row), axis=1
    )
    expected_pairs["ligand_found_post"] = expected_pairs.apply(
        lambda row: _pair_in(post_pairs, row), axis=1
    )
    expected_pairs["ligand_found_any"] = expected_pairs.apply(
        lambda row: _pair_in(any_pairs, row), axis=1
    )

    expected_targets = expected_pairs["uniprot_id"].nunique()
    target_ids = set(expected_pairs["uniprot_id"])
    results_targets = set()
    if results_root.exists():
        results_targets = {
            entry.name for entry in results_root.iterdir() if entry.is_dir()
        }
    targets_with_output = len(target_ids & results_targets)

    per_target = expected_pairs.groupby("uniprot_id").agg(
        expected_pairs=("expected_ligand", "count"),
        expected_pdbs=("expected_pdb", "nunique"),
        expected_ligands=("expected_ligand", "nunique"),
        donor_found=("pdb_found", "any"),
        ligand_found_any=("ligand_found_any", "any"),
        ligand_found_post=("ligand_found_post", "any"),
        ligand_all_found=("ligand_found_any", "all"),
    )

    targets_with_expected_pdb = int(per_target["donor_found"].sum())
    targets_with_expected_ligand_any = int(per_target["ligand_found_any"].sum())
    targets_with_expected_ligand_post = int(per_target["ligand_found_post"].sum())
    targets_with_all_ligands = int(per_target["ligand_all_found"].sum())

    pairs_total = int(expected_pairs.shape[0])
    pairs_with_pdb = int(expected_pairs["pdb_found"].sum())
    pairs_with_pre = int(expected_pairs["ligand_found_pre"].sum())
    pairs_with_post = int(expected_pairs["ligand_found_post"].sum())
    pairs_with_any = int(expected_pairs["ligand_found_any"].sum())

    def _share(numerator: int, denominator: int | float | None) -> float | None:
        if denominator in (None, 0):
            return None
        return 100.0 * numerator / denominator

    summary_rows = [
        {
            "Category": "Target-level",
            "Metric": "Benchmark targets",
            "Count": expected_targets,
            "Baseline": expected_targets,
            "Share (%)": None,
        },
        {
            "Category": "Target-level",
            "Metric": "Targets with pipeline output directory",
            "Count": targets_with_output,
            "Baseline": expected_targets,
            "Share (%)": _share(targets_with_output, expected_targets),
        },
        {
            "Category": "Target-level",
            "Metric": "Targets with expected donor alignment",
            "Count": targets_with_expected_pdb,
            "Baseline": expected_targets,
            "Share (%)": _share(targets_with_expected_pdb, expected_targets),
        },
        {
            "Category": "Target-level",
            "Metric": "Targets with expected ligand recovered (any stage)",
            "Count": targets_with_expected_ligand_any,
            "Baseline": expected_targets,
            "Share (%)": _share(targets_with_expected_ligand_any, expected_targets),
        },
        {
            "Category": "Target-level",
            "Metric": "Targets with expected ligand recovered post-JAMDA",
            "Count": targets_with_expected_ligand_post,
            "Baseline": expected_targets,
            "Share (%)": _share(targets_with_expected_ligand_post, expected_targets),
        },
        {
            "Category": "Target-level",
            "Metric": "Targets with all expected ligands recovered",
            "Count": targets_with_all_ligands,
            "Baseline": expected_targets,
            "Share (%)": _share(targets_with_all_ligands, expected_targets),
        },
        {
            "Category": "Pair-level",
            "Metric": "Expected donor–ligand pairs",
            "Count": pairs_total,
            "Baseline": pairs_total,
            "Share (%)": None,
        },
        {
            "Category": "Pair-level",
            "Metric": "Pairs with expected donor alignment",
            "Count": pairs_with_pdb,
            "Baseline": pairs_total,
            "Share (%)": _share(pairs_with_pdb, pairs_total),
        },
        {
            "Category": "Pair-level",
            "Metric": "Pairs with expected ligand recovered (pre-JAMDA)",
            "Count": pairs_with_pre,
            "Baseline": pairs_total,
            "Share (%)": _share(pairs_with_pre, pairs_total),
        },
        {
            "Category": "Pair-level",
            "Metric": "Pairs with expected ligand recovered (post-JAMDA)",
            "Count": pairs_with_post,
            "Baseline": pairs_total,
            "Share (%)": _share(pairs_with_post, pairs_total),
        },
        {
            "Category": "Pair-level",
            "Metric": "Pairs with expected ligand recovered (any stage)",
            "Count": pairs_with_any,
            "Baseline": pairs_total,
            "Share (%)": _share(pairs_with_any, pairs_total),
        },
    ]

    summary_df = pd.DataFrame(summary_rows)
    summary_df["Share (%)"] = summary_df["Share (%)"].map(
        lambda x: "--" if pd.isna(x) else f"{float(x):.1f}\\%"
    )
    summary_df["Baseline"] = summary_df["Baseline"].map(
        lambda x: "--" if pd.isna(x) else int(x)
    )

    summary_df = summary_df.rename(columns={"Share (%)": r"Share (\%)"})

    summary_df.to_csv(OUTPUT_DIR / "table_expected_recovery.csv", index=False)
    _export_table(
        summary_df,
        "table_expected_recovery.tex",
        caption=(
            "Recovery of benchmark validation donors and ligands across the production run."
        ),
        label="tab:expected-recovery",
        short_caption="Benchmark recovery summary",
    )
    print("[OK] Expected donor/ligand recovery table exported.")

    details = expected_pairs[
        [
            "uniprot_id",
            "expected_pdb",
            "expected_ligand",
            "bin_type",
            "bin_name",
            "pdb_found",
            "ligand_found_pre",
            "ligand_found_post",
            "ligand_found_any",
        ]
    ].sort_values(["uniprot_id", "expected_pdb", "expected_ligand"])
    details.to_csv(OUTPUT_DIR / "expected_pair_recovery.csv", index=False)


def figure_stage_runtimes(stage_df: pd.DataFrame):
    if stage_df.empty:
        print("[WARN] Skipping stage runtime figure (no log data).")
        return

    metrics = [
        ("dogsite3_seconds", "DoGSite3 runtime (s)"),
        ("siena_seconds", "SIENA runtime (s)"),
        ("jamda_seconds", "JAMDA runtime (s)"),
        ("pipeline_seconds", "Total pipeline runtime (s)"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    runtime_summary = []

    for ax, (column, title) in zip(axes.flat, metrics):
        values = pd.to_numeric(
            stage_df.get(column, pd.Series(dtype=float)), errors="coerce"
        ).dropna()
        if values.empty:
            ax.axis("off")
            runtime_summary.append(
                {
                    "metric": column,
                    "n": 0,
                    "median": np.nan,
                    "mean": np.nan,
                    "p95": np.nan,
                }
            )
            continue
        sns.histplot(
            data=values.to_frame("values"),
            x="values",
            bins=40,
            kde=True,
            color="#4C72B0",
            ax=ax,
        )
        median = values.median()
        ax.axvline(median, color="#C44E52", linestyle="--", linewidth=1.5)
        ax.set_title(f"{title} (n={len(values)})")
        ax.set_xlabel("Seconds")
        ax.set_ylabel("Count")
        ax.text(
            0.95,
            0.85,
            f"median = {median:.1f}s",
            transform=ax.transAxes,
            ha="right",
            fontsize=10,
            bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none"},
        )
        runtime_summary.append(
            {
                "metric": column,
                "n": len(values),
                "median": median,
                "mean": values.mean(),
                "p95": values.quantile(0.95),
            }
        )

    plt.tight_layout()
    _save_figure(fig, "fig_stage_runtime_histograms")
    plt.close(fig)

    pd.DataFrame(runtime_summary).to_csv(
        OUTPUT_DIR / "stage_runtime_summary.csv", index=False
    )
    print("[OK] Stage runtime histograms exported (PDF + CSV).")


def figure_top20_ligands(df: pd.DataFrame):
    """Figure: Top 20 most frequently transplanted ligands (post-JAMDA)."""
    df_post = df[df["stage"] == "post-jamda"].copy()
    if df_post.empty:
        print("[SKIP] No post-JAMDA rows for top-20.")
        return
    counts = df_post["ligand_name"].value_counts().head(20).reset_index()
    counts.columns = ["ligand_name", "count"]

    fig, ax = plt.subplots(figsize=(11, 8))
    sns.barplot(data=counts, x="count", y="ligand_name", ax=ax)
    ax.set(
        title="Top 20 most frequently transplanted ligands (post-JAMDA)",
        xlabel="Number of transplants (count)",
        ylabel="Ligand",
    )
    _save_figure(fig, "fig_top20_transplanted_ligands")
    plt.close(fig)
    counts.to_csv(OUTPUT_DIR / "top20_ligands_counts.csv", index=False)
    print("[OK] Figure (top-20 ligands) saved (PDF) + counts CSV.")


def figure_top20_ligands_no_buffers(df: pd.DataFrame):
    df_post = df[df["stage"] == "post-jamda"].copy()
    df_post = df_post[~df_post["ligand_name"].map(_is_buffer_ligand)]
    if df_post.empty:
        print("[SKIP] No post-JAMDA rows for filtered top-20.")
        return
    counts = df_post["ligand_name"].value_counts().head(20).reset_index()
    counts.columns = ["ligand_name", "count"]

    fig, ax = plt.subplots(figsize=(11, 8))
    sns.barplot(data=counts, x="count", y="ligand_name", ax=ax, color="#55A868")
    ax.set(
        title="Top 20 ligands excluding buffer species (post-JAMDA)",
        xlabel="Number of transplants (count)",
        ylabel="Ligand",
    )
    _save_figure(fig, "fig_top20_transplanted_ligands_no_buffers")
    plt.close(fig)
    counts.to_csv(OUTPUT_DIR / "top20_ligands_counts_no_buffers.csv", index=False)
    print("[OK] Figure (top-20 no buffers) saved (PDF) + counts CSV.")


def table_global_quality(df: pd.DataFrame):
    """Section 3.3: Global quality indicators (post-JAMDA)."""
    post = df[df["stage"] == "post-jamda"].copy()
    if post.empty:
        print("[SKIP] No post-JAMDA rows for global quality.")
        return

    metric_map = [
        ("ff_local_rmsd", r"Local RMSD (\si{\angstrom})"),
        ("ff_tcs_score", "TCS (dimensionless)"),
        ("ff_lev_score", r"LEV (\si{\angstrom})"),
    ]

    subsets = {
        "All placements": post,
        "Organic (no ions)": post.loc[~post["is_ion"]],
        "Monoatomic ions": post.loc[post["is_ion"]],
    }

    summary_rows = []
    for subset_name, subset_df in subsets.items():
        for metric, pretty in metric_map:
            series = subset_df[metric]
            stats_dict = _series_stats(series)
            med_lo, med_hi = _bootstrap_ci(series, np.median)
            mean_lo, mean_hi = _mean_ci(series)
            summary_rows.append(
                {
                    "Subset": subset_name,
                    "Metric": pretty,
                    "N": stats_dict["n"],
                    "Mean": stats_dict["mean"],
                    "Mean_CI_low": mean_lo,
                    "Mean_CI_high": mean_hi,
                    "Median": stats_dict["median"],
                    "Median_CI_low": med_lo,
                    "Median_CI_high": med_hi,
                    "Q25": stats_dict["q25"],
                    "Q75": stats_dict["q75"],
                }
            )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUTPUT_DIR / "global_quality_stats.csv", index=False)

    latex = summary_df.copy()
    latex["Mean (95% CI)"] = latex.apply(
        lambda row: _format_ci_value(
            row["Mean"],
            row["Mean_CI_low"],
            row["Mean_CI_high"],
            precision=_metric_precision(row["Metric"]),
        ),
        axis=1,
    )
    latex["Median (95% CI)"] = latex.apply(
        lambda row: _format_ci_value(
            row["Median"],
            row["Median_CI_low"],
            row["Median_CI_high"],
            precision=_metric_precision(row["Metric"]),
        ),
        axis=1,
    )
    latex["Q25"] = latex.apply(
        lambda row: _format_metric_scalar(row["Q25"], row["Metric"]), axis=1
    )
    latex["Q75"] = latex.apply(
        lambda row: _format_metric_scalar(row["Q75"], row["Metric"]), axis=1
    )

    latex = latex[
        [
            "Subset",
            "Metric",
            "N",
            "Mean (95% CI)",
            "Median (95% CI)",
            "Q25",
            "Q75",
        ]
    ]

    _export_table(
        latex,
        "table_global_quality_stats.tex",
        caption=(
            r"Global quality indicators with 95\% confidence intervals for post-JAMDA "
            r"placements, reported for all ligands, organic ligands only, and "
            r"monoatomic ions."
        ),
        label="tab:global-quality-stats",
        short_caption="Global quality indicators",
    )

    # Correlations (Pearson + Spearman) among metrics
    corr_pairs = [
        ("ff_lev_score", "ff_tcs_score"),
        ("ff_lev_score", "ff_local_rmsd"),
        ("ff_local_rmsd", "ff_tcs_score"),
    ]
    corr_base = post[
        (~post["is_ion"])
        & post["ff_tcs_score"].notna()
        & post["ff_local_rmsd"].notna()
        & post["ff_lev_score"].notna()
    ].copy()
    if corr_base.empty:
        print("[WARN] No rows available for cross-metric correlation table.")
        corr_df = pd.DataFrame()
    else:
        rows = []
        for a, b in corr_pairs:
            subset = corr_base[[a, b]].dropna()
            if subset.empty:
                continue
            a_s = pd.to_numeric(subset[a], errors="coerce")
            b_s = pd.to_numeric(subset[b], errors="coerce")
            mask = a_s.notna() & b_s.notna()
            if not mask.any():
                continue
            pearson_r = a_s[mask].corr(b_s[mask])
            pearson_lo, pearson_hi = _fisher_ci(pearson_r, int(mask.sum()))
            spearman_rho, _ = stats.spearmanr(a_s[mask], b_s[mask])
            spearman_lo, spearman_hi = _fisher_ci(float(spearman_rho), int(mask.sum()))
            rows.append(
                {
                    "metric_x": a,
                    "metric_y": b,
                    "pearson_r": pearson_r,
                    "pearson_lo": pearson_lo,
                    "pearson_hi": pearson_hi,
                    "spearman_rho": spearman_rho,
                    "spearman_lo": spearman_lo,
                    "spearman_hi": spearman_hi,
                    "n": int(mask.sum()),
                }
            )
        corr_df = pd.DataFrame(rows)
    corr_df.to_csv(OUTPUT_DIR / "global_quality_correlations.csv", index=False)
    latex_corr = corr_df.copy()
    latex_corr["metric_x"] = latex_corr["metric_x"].map(_format_metric_label)
    latex_corr["metric_y"] = latex_corr["metric_y"].map(_format_metric_label)
    latex_corr["Pearson (95% CI)"] = latex_corr.apply(
        lambda row: _format_ci_value(
            row["pearson_r"], row["pearson_lo"], row["pearson_hi"]
        ),
        axis=1,
    )
    latex_corr["Spearman (95% CI)"] = latex_corr.apply(
        lambda row: _format_ci_value(
            row["spearman_rho"], row["spearman_lo"], row["spearman_hi"]
        ),
        axis=1,
    )
    latex_corr = latex_corr[
        [
            "metric_x",
            "metric_y",
            "n",
            "Pearson (95% CI)",
            "Spearman (95% CI)",
        ]
    ]
    latex_corr = latex_corr.rename(
        columns={"metric_x": "Metric X", "metric_y": "Metric Y", "n": "N"}
    )
    for col in ["Metric X", "Metric Y"]:
        latex_corr[col] = latex_corr[col].str.replace(
            "Å", r"\si{\angstrom}", regex=False
        )

    _export_table(
        latex_corr,
        "table_cross_metric_correlations.tex",
        caption=(
            r"Cross-metric correlations for post-JAMDA placements (organic ligands with "
            r"all metrics defined) with 95\% confidence intervals (Pearson $r$ and "
            r"Spearman $\rho$; $n$ denotes paired observations)."
        ),
        label="tab:cross-metric",
        short_caption="Cross-metric correlations",
    )

    # Spearman with active-site identity
    if "active_site_identity" in post.columns:
        try:
            mask_rmsd = (
                post["active_site_identity"].notna() & post["ff_local_rmsd"].notna()
            )
            rho_rmsd, p_rmsd = stats.spearmanr(
                post.loc[mask_rmsd, "active_site_identity"],
                post.loc[mask_rmsd, "ff_local_rmsd"],
            )
            mask_tcs = (
                post["active_site_identity"].notna()
                & post["ff_tcs_score"].notna()
                & (~post["is_ion"])
            )
            rho_tcs, p_tcs = stats.spearmanr(
                post.loc[mask_tcs, "active_site_identity"],
                post.loc[mask_tcs, "ff_tcs_score"],
            )
        except Exception:
            # Fallback: compute rho via pandas (no p-values)
            corr_mat = post[["active_site_identity", "ff_local_rmsd"]].corr(
                method="spearman"
            )
            rmsd_corr_val = corr_mat.iloc[0, 1]
            try:
                if pd.isna(rmsd_corr_val):
                    rho_rmsd = 0.0
                elif hasattr(rmsd_corr_val, "real"):
                    rho_rmsd = float(rmsd_corr_val.real)  # type: ignore
                else:
                    rho_rmsd = float(rmsd_corr_val)  # type: ignore
            except (TypeError, ValueError, AttributeError):
                rho_rmsd = 0.0
            corr_mat = post.loc[
                ~post["is_ion"], ["active_site_identity", "ff_tcs_score"]
            ].corr(method="spearman")
            tcs_corr_val = corr_mat.iloc[0, 1]
            try:
                if pd.isna(tcs_corr_val):
                    rho_tcs = 0.0
                elif hasattr(tcs_corr_val, "real"):
                    rho_tcs = float(tcs_corr_val.real)  # type: ignore
                else:
                    rho_tcs = float(tcs_corr_val)  # type: ignore
            except (TypeError, ValueError, AttributeError):
                rho_tcs = 0.0
            p_rmsd = np.nan
            p_tcs = np.nan
            mask_rmsd = (
                post["active_site_identity"].notna() & post["ff_local_rmsd"].notna()
            )
            mask_tcs = (
                post["active_site_identity"].notna()
                & post["ff_tcs_score"].notna()
                & (~post["is_ion"])
            )
        rmsd_ci = _fisher_ci(float(rho_rmsd), int(mask_rmsd.sum()))
        tcs_ci = _fisher_ci(float(rho_tcs), int(mask_tcs.sum()))
        spearman_df = pd.DataFrame(
            [
                {
                    "metric": r"Local RMSD (\si{\angstrom})",
                    "spearman_rho": rho_rmsd,
                    "ci_low": rmsd_ci[0],
                    "ci_high": rmsd_ci[1],
                    "p_value": p_rmsd,
                    "n": int(mask_rmsd.sum()),
                },
                {
                    "metric": "TCS (dimensionless)",
                    "spearman_rho": rho_tcs,
                    "ci_low": tcs_ci[0],
                    "ci_high": tcs_ci[1],
                    "p_value": p_tcs,
                    "n": int(mask_tcs.sum()),
                },
            ]
        )
        spearman_df.to_csv(OUTPUT_DIR / "identity_spearman.csv", index=False)
        latex_spear = spearman_df.copy()
        latex_spear["Spearman (95% CI)"] = latex_spear.apply(
            lambda row: _format_ci_value(
                row["spearman_rho"], row["ci_low"], row["ci_high"]
            ),
            axis=1,
        )
        latex_spear["p-value"] = latex_spear["p_value"].map(
            lambda x: f"{x:.2e}" if pd.notnull(x) else "--"
        )
        latex_spear = latex_spear.rename(columns={"metric": "Metric", "n": "N"})[
            ["Metric", "N", "Spearman (95% CI)", "p-value"]
        ]
        _export_table(
            latex_spear,
            "table_identity_spearman.tex",
            caption=(
                r"Spearman correlation of active-site identity with key metrics "
                r"(post-JAMDA), reported with 95\% confidence intervals and sample size."
            ),
            label="tab:identity-spearman",
            short_caption="Active-site identity correlations",
        )
    print("[OK] Global quality indicators exported.")


def figure_alignment_by_identity(df: pd.DataFrame):
    """Section 3.4: Distributions by active-site identity bins (post-JAMDA)."""
    post = df[df["stage"] == "post-jamda"].copy()
    if post.empty or "active_site_identity" not in post.columns:
        print("[SKIP] No data for identity-binned plots.")
        return
    post = post[post["active_site_identity"].notna()].copy()
    # High-identity bins within [0.80, 1.00] with finer resolution
    display_bins = [0.80, 0.85, 0.90, 0.95, 1.00]
    cut_bins = display_bins.copy()
    cut_bins[-1] = cut_bins[-1] + 1e-6  # ensure identity == 1.00 is retained

    def _fmt(v: float) -> str:
        s = f"{v:.2f}".rstrip("0").rstrip(".")
        return s

    labels = [
        f"{_fmt(display_bins[i])}-{_fmt(display_bins[i + 1])}"
        for i in range(len(display_bins) - 1)
    ]
    post["identity_bin"] = pd.cut(
        post["active_site_identity"],
        bins=cut_bins,
        labels=labels,
        right=False,
        include_lowest=True,
    )
    # Keep only rows that fall into the defined high-identity bins
    post = post[post["identity_bin"].notna()].copy()
    bin_counts = (
        post["identity_bin"].value_counts().reindex(labels).fillna(0).astype(int)
    )

    for metric in ["ff_local_rmsd", "ff_tcs_score"]:
        pretty, _ = METRIC_LABELS[metric]
        fig, ax = plt.subplots(figsize=(12, 7))
        sns.boxplot(
            data=post.dropna(subset=[metric, "identity_bin"]),
            x="identity_bin",
            y=metric,
            showfliers=False,
            order=labels,
            ax=ax,
        )
        median_per_bin = (
            post.dropna(subset=[metric, "identity_bin"])
            .groupby("identity_bin")[metric]
            .median()
            .reindex(labels)
        )
        ax.plot(
            list(range(len(labels))),
            median_per_bin.values.tolist(),
            color="#C44E52",
            marker="o",
            linestyle="--",
            label="Bin median",
        )
        ax.set(
            title=f"Distribution of {pretty} by active-site identity (post-JAMDA)",
            xlabel="Active-site identity bin",
            ylabel=_format_metric_label(metric),
        )
        ax.set_xticks(range(len(labels)))
        xtick_labels = [f"{lab}\n(n={bin_counts.get(lab, 0)})" for lab in labels]
        ax.set_xticklabels(xtick_labels)
        ax.legend(loc="upper right")
        _save_figure(fig, f"fig_distribution_{metric}_by_identity_bin")
        plt.close(fig)
    print("[OK] Alignment-by-identity figures exported (PDF).")


def figure_jamda_effect(df: pd.DataFrame):
    """Section 3.5: JAMDA pre vs post comparisons + summary table of deltas."""
    pre = df[df["stage"] == "pre-jamda"].copy()
    post = df[df["stage"] == "post-jamda"].copy()
    key = ["uniprot_id", "pdb_code", "ligand_id"]
    merged = pd.merge(pre, post, on=key, suffixes=("_pre", "_post"))
    merged = merged.dropna(
        subset=[
            "ff_tcs_score_pre",
            "ff_tcs_score_post",
            "ff_local_rmsd_pre",
            "ff_local_rmsd_post",
        ]
    )
    if merged.empty:
        print("[SKIP] No matched pre/post rows for JAMDA analysis.")
        return

    merged["tcs_improvement"] = merged["ff_tcs_score_pre"] - merged["ff_tcs_score_post"]
    merged["rmsd_change"] = merged["ff_local_rmsd_post"] - merged["ff_local_rmsd_pre"]
    merged.to_csv(OUTPUT_DIR / "jamda_effectiveness_pairs.csv", index=False)

    # Scatter: pre vs post TCS, colored by pre quartile
    q1 = merged["ff_tcs_score_pre"].quantile(0.25)
    q2 = merged["ff_tcs_score_pre"].quantile(0.50)
    q3 = merged["ff_tcs_score_pre"].quantile(0.75)

    def _grp(v):
        return (
            "Q1 lowest"
            if v <= q1
            else "Q2"
            if v <= q2
            else "Q3"
            if v <= q3
            else "Q4 highest"
        )

    merged["tcs_group"] = merged["ff_tcs_score_pre"].apply(_grp)

    fig, ax = plt.subplots(figsize=(8, 8))
    sns.scatterplot(
        data=merged,
        x="ff_tcs_score_pre",
        y="ff_tcs_score_post",
        hue="tcs_group",
        alpha=0.7,
        edgecolor="w",
        ax=ax,
    )
    ax.plot(
        [0, max(1.5, merged["ff_tcs_score_pre"].max() * 1.05)],
        [0, max(1.5, merged["ff_tcs_score_pre"].max() * 1.05)],
        linestyle="--",
    )
    ax.set(
        title="JAMDA refinement effect on TCS",
        xlabel="TCS before (dimensionless)",
        ylabel="TCS after (dimensionless)",
    )
    _save_figure(fig, "fig_jamda_tcs_pre_vs_post")
    plt.close(fig)

    # Violin/box plots for paired deltas
    delta_long = pd.melt(
        merged[["tcs_improvement", "rmsd_change"]],
        value_vars=["tcs_improvement", "rmsd_change"],
        var_name="metric",
        value_name="delta",
    )
    metric_labels = {
        "tcs_improvement": "ΔTCS (pre - post, dimensionless)",
        "rmsd_change": "ΔLocal RMSD (post - pre, \u212b)",
    }
    delta_long["metric_label"] = delta_long["metric"].map(metric_labels)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    for ax, (metric_key, label) in zip(axes, metric_labels.items()):
        data = delta_long[delta_long["metric"] == metric_key]
        if data.empty:
            ax.axis("off")
            continue
        sns.violinplot(
            data=data, x="metric_label", y="delta", inner=None, color="#4C72B0", ax=ax
        )
        sns.boxplot(
            data=data,
            x="metric_label",
            y="delta",
            width=0.2,
            color="white",
            showcaps=True,
            boxprops={"zorder": 3},
            ax=ax,
        )
        ax.axhline(0, color="#C44E52", linestyle="--", linewidth=1)
        median_val = data["delta"].median()
        ax.text(
            0.0,
            median_val,
            f"median={median_val:.3f}",
            ha="center",
            va="bottom" if median_val >= 0 else "top",
            fontsize=10,
        )
        ax.set_xlabel("")
        ax.set_ylabel(label)
        ax.set_title(f"{label} (n={len(data)})")
    plt.tight_layout()
    _save_figure(fig, "fig_jamda_delta_distributions")
    plt.close(fig)

    # Summary deltas
    tcs_med_ci = _bootstrap_ci(merged["tcs_improvement"], np.median)
    tcs_mean_ci = _mean_ci(merged["tcs_improvement"])
    rmsd_med_ci = _bootstrap_ci(merged["rmsd_change"], np.median)
    rmsd_mean_ci = _mean_ci(merged["rmsd_change"])
    delta = pd.DataFrame(
        {
            "metric": ["TCS (pre - post)", "Local RMSD (post - pre)"],
            "median_delta": [
                merged["tcs_improvement"].median(),
                merged["rmsd_change"].median(),
            ],
            "mean_delta": [
                merged["tcs_improvement"].mean(),
                merged["rmsd_change"].mean(),
            ],
            "q25_delta": [
                merged["tcs_improvement"].quantile(0.25),
                merged["rmsd_change"].quantile(0.25),
            ],
            "q75_delta": [
                merged["tcs_improvement"].quantile(0.75),
                merged["rmsd_change"].quantile(0.75),
            ],
            "mean_ci_low": [tcs_mean_ci[0], rmsd_mean_ci[0]],
            "mean_ci_high": [tcs_mean_ci[1], rmsd_mean_ci[1]],
            "median_ci_low": [tcs_med_ci[0], rmsd_med_ci[0]],
            "median_ci_high": [tcs_med_ci[1], rmsd_med_ci[1]],
            "n_pairs": [int(merged.shape[0]), int(merged.shape[0])],
        }
    )
    delta.to_csv(OUTPUT_DIR / "table_jamda_delta_summary.csv", index=False)
    latex = delta.copy()
    for col in [
        "median_delta",
        "median_ci_low",
        "median_ci_high",
        "mean_delta",
        "mean_ci_low",
        "mean_ci_high",
        "q25_delta",
        "q75_delta",
    ]:
        latex[col] = latex[col].map(lambda x: f"{x:.3f}")
    latex["Median (95% CI)"] = latex.apply(
        lambda row: f"{row['median_delta']}\\;[ {row['median_ci_low']}, {row['median_ci_high']} ]",
        axis=1,
    )
    latex["Mean (95% CI)"] = latex.apply(
        lambda row: f"{row['mean_delta']}\\;[ {row['mean_ci_low']}, {row['mean_ci_high']} ]",
        axis=1,
    )
    latex = latex[
        [
            "metric",
            "n_pairs",
            "Mean (95% CI)",
            "Median (95% CI)",
            "q25_delta",
            "q75_delta",
        ]
    ]
    latex = latex.rename(
        columns={
            "metric": "Metric",
            "n_pairs": "n",
            "q25_delta": "Q25",
            "q75_delta": "Q75",
        }
    )
    _export_table(
        latex,
        "table_jamda_delta_summary.tex",
        caption=(
            r"Summary of JAMDA refinement deltas (matched pairs) with 95\% confidence "
            r"intervals for mean and median changes."
        ),
        label="tab:jamda-deltas",
        short_caption="JAMDA refinement delta summary",
    )
    print("[OK] JAMDA figures (PDF) + delta table exported.")


def figure_ff_vs_af(df: pd.DataFrame):
    """Section 3.6: Correlations between FoldFusion and AlphaFill (post-JAMDA)."""
    post = df[df["stage"] == "post-jamda"].copy()
    sub = post.dropna(
        subset=["ff_local_rmsd", "af_local_rmsd", "ff_tcs_score", "af_tcs_score"]
    )
    if sub.empty:
        print("[SKIP] No overlapping FF/AF metrics for comparison.")
        return

    metric_display = {
        "local_rmsd": ("Local RMSD", "Å"),
        "tcs_score": ("TCS", "dimensionless"),
    }
    rows = []
    for metric, (pretty, unit) in metric_display.items():
        fig, ax = plt.subplots(figsize=(6, 6))
        sns.regplot(
            data=sub,
            x=f"af_{metric}",
            y=f"ff_{metric}",
            scatter_kws={"alpha": 0.3, "edgecolor": "w"},
            ax=ax,
        )
        valid = sub.dropna(subset=[f"af_{metric}", f"ff_{metric}"])
        r = valid[f"af_{metric}"].corr(valid[f"ff_{metric}"])
        pearson_lo, pearson_hi = _fisher_ci(r, int(valid.shape[0]))
        spearman_rho, _ = stats.spearmanr(valid[f"af_{metric}"], valid[f"ff_{metric}"])
        spearman_lo, spearman_hi = _fisher_ci(float(spearman_rho), int(valid.shape[0]))
        n_pts = int(valid.shape[0])
        unit_suffix = " (dimensionless)" if unit == "dimensionless" else f" ({unit})"
        ax.set(
            title=(
                f"This work vs AlphaFill: {pretty} (Pearson r={r:.2f}, "
                f"Spearman ρ={spearman_rho:.2f}, n={n_pts})"
            ),
            xlabel=f"AlphaFill {pretty}{unit_suffix}",
            ylabel=f"This work {pretty}{unit_suffix}",
        )
        _save_figure(fig, f"fig_correlation_{metric}_ff_vs_af")
        plt.close(fig)
        rows.append(
            {
                "metric": pretty,
                "pearson_r": r,
                "pearson_lo": pearson_lo,
                "pearson_hi": pearson_hi,
                "spearman_rho": float(spearman_rho),
                "spearman_lo": spearman_lo,
                "spearman_hi": spearman_hi,
                "n": n_pts,
            }
        )

    corr = pd.DataFrame(rows)
    corr.to_csv(OUTPUT_DIR / "table_ff_vs_af_correlations.csv", index=False)

    latex = corr.copy()
    latex["Pearson (95% CI)"] = latex.apply(
        lambda row: _format_ci_value(
            row["pearson_r"], row["pearson_lo"], row["pearson_hi"]
        ),
        axis=1,
    )
    latex["Spearman (95% CI)"] = latex.apply(
        lambda row: _format_ci_value(
            row["spearman_rho"], row["spearman_lo"], row["spearman_hi"]
        ),
        axis=1,
    )
    latex = latex.rename(
        columns={
            "metric": "Metric",
            "n": "N",
        }
    )[["Metric", "Pearson (95% CI)", "Spearman (95% CI)", "N"]]
    _export_table(
        latex,
        "table_ff_vs_af_correlations.tex",
        caption=(
            "Pearson and Spearman correlation between this pipeline and AlphaFill metrics "
            "(post-JAMDA)."
        ),
        label="tab:ff-af-corr",
        short_caption="FoldFusion vs AlphaFill correlations",
    )
    print("[OK] FF vs AF figures (PDF) + correlation table exported.")


def figure_cross_metric(df: pd.DataFrame):
    """Section 3.7: Cross-metric relationships inside FoldFusion (post-JAMDA)."""
    post = df[
        (df["stage"] == "post-jamda")
        & df["ff_tcs_score"].notna()
        & df["ff_lev_score"].notna()
        & df["ff_local_rmsd"].notna()
    ].copy()
    post = post[~post["is_ion"]].copy()
    if post.empty:
        print("[SKIP] No data for cross-metric plots.")
        return

    pairs = [
        ("ff_lev_score", "ff_tcs_score", "fig_corr_tcs_vs_lev"),
        ("ff_lev_score", "ff_local_rmsd", "fig_corr_rmsd_vs_lev"),
        ("ff_local_rmsd", "ff_tcs_score", "fig_corr_rmsd_vs_tcs"),
    ]
    rows = []
    for x, y, stem in pairs:
        xlab = _format_metric_label(x)
        ylab = _format_metric_label(y)
        fig, ax = plt.subplots(figsize=(6, 6))
        sns.scatterplot(data=post, x=x, y=y, alpha=0.5, ax=ax)
        valid = post[[x, y]].dropna()
        pearson_r = valid[x].corr(valid[y])
        pearson_lo, pearson_hi = _fisher_ci(pearson_r, int(valid.shape[0]))
        spearman_rho, _ = stats.spearmanr(valid[x], valid[y])
        spearman_lo, spearman_hi = _fisher_ci(float(spearman_rho), int(valid.shape[0]))
        n_pts = int(valid.shape[0])
        ax.set(
            title=(
                f"{xlab} vs {ylab} (Pearson r={pearson_r:.2f}, Spearman ρ={spearman_rho:.2f}, "
                f"n={n_pts})\nOrganic placements with all metrics defined"
            ),
            xlabel=xlab,
            ylabel=ylab,
        )
        _save_figure(fig, stem)
        plt.close(fig)
        rows.append(
            {
                "x_metric": xlab,
                "y_metric": ylab,
                "pearson_r": pearson_r,
                "pearson_lo": pearson_lo,
                "pearson_hi": pearson_hi,
                "spearman_rho": spearman_rho,
                "spearman_lo": spearman_lo,
                "spearman_hi": spearman_hi,
                "n": n_pts,
            }
        )

    corr = pd.DataFrame(rows)
    corr.to_csv(OUTPUT_DIR / "table_cross_metric_correlations_ff_only.csv", index=False)
    latex_corr = corr.copy()
    latex_corr["Pearson (95% CI)"] = latex_corr.apply(
        lambda row: _format_ci_value(
            row["pearson_r"], row["pearson_lo"], row["pearson_hi"]
        ),
        axis=1,
    )
    latex_corr["Spearman (95% CI)"] = latex_corr.apply(
        lambda row: _format_ci_value(
            row["spearman_rho"], row["spearman_lo"], row["spearman_hi"]
        ),
        axis=1,
    )
    latex_corr = latex_corr[
        [
            "x_metric",
            "y_metric",
            "n",
            "Pearson (95% CI)",
            "Spearman (95% CI)",
        ]
    ]
    latex_corr = latex_corr.rename(
        columns={"x_metric": "Metric X", "y_metric": "Metric Y", "n": "N"}
    )
    for col in ["Metric X", "Metric Y"]:
        latex_corr[col] = latex_corr[col].str.replace(
            "Å", r"\si{\angstrom}", regex=False
        )

    _export_table(
        latex_corr,
        "table_cross_metric_correlations_ff_only.tex",
        caption=(
            r"Pearson and Spearman correlations among this pipeline's metrics (post-JAMDA) "
            r"with 95\% confidence intervals; $n$ denotes paired observations."
        ),
        label="tab:ff-cross-metric",
        short_caption="FoldFusion cross-metric correlations",
    )
    print("[OK] Cross-metric figures (PDF) + table exported.")


def identity_threshold_ablation(df: pd.DataFrame):
    post = df[df["stage"] == "post-jamda"].copy()
    post = post[post["active_site_identity"].notna()]
    if post.empty:
        print("[SKIP] No data for identity threshold ablation.")
        return

    # The production pipeline enforces a minimum active-site identity of 0.85.
    # Explore tightening the cut further by sweeping thresholds from 0.85 up to 1.00.
    thresholds = np.round(np.arange(0.85, 1.001, 0.05), 2)
    total_placements = len(post)
    rows = []
    for thr in thresholds:
        subset = post.loc[post["active_site_identity"] >= thr]
        if subset.empty:
            rows.append(
                {
                    "threshold": thr,
                    "placements": 0,
                    "target_count": 0,
                    "success_rate": np.nan,
                    "median_rmsd": np.nan,
                    "median_tcs": np.nan,
                }
            )
            continue
        if "is_ion" in subset.columns:
            organic_mask = ~subset["is_ion"].fillna(False)
        else:
            organic_mask = ~subset.get(
                "ff_tcs_atoms",
                pd.Series(0, index=subset.index, dtype="int64"),
            ).fillna(0).astype(int).eq(1)

        success_mask = (
            organic_mask
            & (subset["ff_local_rmsd"] <= LOCAL_RMSD_THRESHOLD)
            & (subset["ff_tcs_score"] <= TCS_THRESHOLD)
        )
        success_count = int(success_mask.sum())
        rows.append(
            {
                "threshold": thr,
                "placements": int(len(subset)),
                "target_count": int(subset["uniprot_id"].nunique()),
                "success_rate": success_count / len(subset)
                if len(subset)
                else np.nan,
                "median_rmsd": subset["ff_local_rmsd"].median(),
                "median_tcs": subset["ff_tcs_score"].median(),
            }
        )

    sweep_df = pd.DataFrame(rows)
    sweep_df.to_csv(OUTPUT_DIR / "identity_threshold_sweep.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    (line1,) = axes[0].plot(
        sweep_df["threshold"], sweep_df["placements"], marker="o", label="Placements"
    )
    axes[0].set_xlabel("Minimum active-site identity")
    axes[0].set_ylabel("Placements retained (count)")
    axes[0].set_title("Yield vs. identity threshold")
    axes[0].grid(True, alpha=0.3)

    ax2 = axes[0].twinx()
    (line2,) = ax2.plot(
        sweep_df["threshold"],
        sweep_df["success_rate"],
        color="#C44E52",
        marker="s",
        label="Success rate",
    )
    ax2.set_ylabel(
        f"Success rate (fraction; TCS ≤ {TCS_THRESHOLD:.3f}, RMSD ≤ {LOCAL_RMSD_THRESHOLD:.2f}\u212b)"
    )
    ax2.set_ylim(0, 1)
    axes[0].legend([line1, line2], ["Placements", "Success rate"], loc="upper right")

    axes[1].plot(
        sweep_df["threshold"],
        sweep_df["median_rmsd"],
        marker="o",
        label="Median Local RMSD",
    )
    axes[1].plot(
        sweep_df["threshold"], sweep_df["median_tcs"], marker="s", label="Median TCS"
    )
    axes[1].set_xlabel("Minimum active-site identity")
    axes[1].set_ylabel("Median metric value (Local RMSD in \u212b, TCS dimensionless)")
    axes[1].set_title("Quality vs. identity threshold")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    _save_figure(fig, "fig_identity_threshold_sweep")
    plt.close(fig)

    latex = sweep_df.copy()
    latex["success_rate"] = latex["success_rate"].map(
        lambda x: f"{x * 100:.1f}\\%" if pd.notnull(x) else "--"
    )
    latex["placements_share"] = latex["placements"].map(
        lambda x: f"{(x / total_placements) * 100:.1f}\\%" if total_placements else "--"
    )
    latex.rename(
        columns={
            "threshold": "Min. identity (dimensionless)",
            "placements": "Placements",
            "target_count": "Targets",
            "success_rate": "Success rate (%)",
            "placements_share": "Share of total (%)",
            "median_rmsd": r"Median Local RMSD (\si{\angstrom})",
            "median_tcs": "Median TCS (dimensionless)",
        },
        inplace=True,
    )
    _export_table(
        latex[
            [
                "Min. identity (dimensionless)",
                "Placements",
                "Share of total (%)",
                "Targets",
                "Success rate (%)",
                r"Median Local RMSD (\si{\angstrom})",
                "Median TCS (dimensionless)",
            ]
        ],
        "table_identity_threshold_sweep.tex",
        caption=(
            "Effect of tightening the minimum active-site identity threshold on "
            "placement yield and quality (post-JAMDA). Success rates use the "
            "placement counts shown; the operational cut at 0.85 corresponds to "
            r"3592/6658 = 54.0\%."
        ),
        label="tab:identity-threshold-sweep",
        short_caption="Identity threshold sweep summary",
    )
    print("[OK] Identity threshold ablation exported.")


def appendix_ion_breakdown_table(df: pd.DataFrame):
    post = df[df["stage"] == "post-jamda"].copy()
    if post.empty:
        print("[SKIP] No data for appendix ion breakdown table.")
        return

    if "is_ion" not in post.columns:
        post["is_ion"] = post["ff_tcs_atoms"].fillna(0).astype(int).eq(1)

    organic_count = int(post[~post["is_ion"]].shape[0])
    ion_count = int(post[post["is_ion"]].shape[0])
    total = organic_count + ion_count
    if total == 0:
        print("[WARN] No placements available for appendix ion breakdown table.")
        return

    def _fraction(count: int) -> str:
        return f"\\SI{{{(count / total) * 100:.1f}}}{{\\percent}}"

    rows = [
        {
            "Category": "Organic / polyatomic ligands",
            "Placements": organic_count,
            "Fraction": _fraction(organic_count),
            "Notes": (
                "Valid \\gls{tcs} and \\gls{lev} metrics available when donors "
                "provide references"
            ),
        },
        {
            "Category": "Monoatomic ions",
            "Placements": ion_count,
            "Fraction": _fraction(ion_count),
            "Notes": (
                "\\gls{tcs} excluded from the main statistics; reported separately in plots"
            ),
        },
    ]

    latex = pd.DataFrame(rows)
    _export_table(
        latex,
        "table_appendix_ion_breakdown.tex",
        caption=(
            "Ligand-class distribution for the post-JAMDAscorer placements "
            "with defined \\gls{tcs} metrics. Totals exclude cases where ligand "
            "extraction failed (none in the final run)."
        ),
        label="tab:appendix-ion-breakdown",
        short_caption="Ligand class distribution",
        directories=APPENDIX_TABLE_OUTPUT_DIRS,
    )
    print("[OK] Appendix ion breakdown table exported.")


def ion_exclusion_analysis(df: pd.DataFrame):
    post = df[df["stage"] == "post-jamda"].copy()
    if post.empty:
        print("[SKIP] No data for ion exclusion analysis.")
        return

    all_tcs = pd.to_numeric(post["ff_tcs_score"], errors="coerce").dropna()
    organic_tcs = pd.to_numeric(
        post.loc[~post["is_ion"], "ff_tcs_score"], errors="coerce"
    ).dropna()
    if not all_tcs.empty and not organic_tcs.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(all_tcs, bins=50, alpha=0.3, label="All placements", density=True)
        ax.hist(organic_tcs, bins=50, alpha=0.3, label="Organic only", density=True)
        ax.axvline(TCS_THRESHOLD, color="#C44E52", linestyle="--", linewidth=1)
        ax.set(
            xlabel="TCS (dimensionless)",
            ylabel="Density",
            title="TCS distribution with and without monoatomic ions",
        )
        ax.legend()
        _save_figure(fig, "fig_tcs_distribution_ion_exclusion")
        plt.close(fig)

    rows = []
    for label, subset in [
        ("Placements with Local RMSD defined", post),
        ("Organic only", post.loc[~post["is_ion"]]),
        ("Monoatomic ions", post.loc[post["is_ion"]]),
    ]:
        subset = subset.dropna(subset=["ff_tcs_score", "ff_local_rmsd"])
        if subset.empty:
            rows.append(
                {
                    "Subset": label,
                    "Placements": 0,
                    "Median TCS": np.nan,
                    "Median Local RMSD": np.nan,
                    "TCS pass rate": np.nan,
                    "Joint pass rate": np.nan,
                }
            )
            continue
        tcs_pass = (subset["ff_tcs_score"] <= TCS_THRESHOLD).mean()
        joint_pass = (
            (subset["ff_tcs_score"] <= TCS_THRESHOLD)
            & (subset["ff_local_rmsd"] <= LOCAL_RMSD_THRESHOLD)
        ).mean()
        rows.append(
            {
                "Subset": label,
                "Placements": int(len(subset)),
                "Median TCS": subset["ff_tcs_score"].median(),
                "Median Local RMSD": subset["ff_local_rmsd"].median(),
                "TCS pass rate": tcs_pass,
                "Joint pass rate": joint_pass,
            }
        )

    ion_table = pd.DataFrame(rows)
    ion_table.to_csv(OUTPUT_DIR / "ion_exclusion_summary.csv", index=False)

    latex = ion_table.copy()
    for col in ["Median TCS", "Median Local RMSD"]:
        latex[col] = latex[col].map(lambda x: f"{x:.3f}" if pd.notnull(x) else "--")
    for col in ["TCS pass rate", "Joint pass rate"]:
        latex[col] = latex[col].map(
            lambda x: f"{x * 100:.1f}\\%" if pd.notnull(x) else "--"
        )
    latex.rename(
        columns={
            "Median TCS": "Median TCS (dimensionless)",
            "Median Local RMSD": r"Median Local RMSD (\si{\angstrom})",
            "TCS pass rate": "TCS pass rate (%)",
            "Joint pass rate": "Joint pass rate (%)",
        },
        inplace=True,
    )

    _export_table(
        latex,
        "table_ion_exclusion_summary.tex",
        caption=(
            "Effect of excluding monoatomic ions on TCS success metrics (post-JAMDA)."
        ),
        label="tab:ion-exclusion",
        short_caption="Ion exclusion TCS summary",
    )
    print("[OK] Ion exclusion analysis exported.")


def failure_mode_case_studies(df: pd.DataFrame):
    post = df[df["stage"] == "post-jamda"].copy()
    if post.empty:
        print("[SKIP] No data for failure mode case studies.")
        return

    def _select_case(condition, sort_key, ascending=False):
        subset = post.loc[condition].dropna(subset=[sort_key])
        if subset.empty:
            subset = post.dropna(subset=[sort_key])
        if subset.empty:
            return None
        return subset.sort_values(sort_key, ascending=ascending).iloc[0]

    cases = []

    collapsed = _select_case(
        (~post["is_ion"])
        & (post["ff_local_rmsd"] <= LOCAL_RMSD_THRESHOLD)
        & (post["ff_tcs_score"] >= 0.6),
        "ff_tcs_score",
        ascending=False,
    )
    if collapsed is not None:
        cases.append(("Collapsed loop (high TCS)", collapsed))

    metal_case = _select_case(
        (post["is_ion"]) & (post["ligand_name"].isin(METAL_LIGANDS)),
        "ff_tcs_score",
        ascending=False,
    )
    if metal_case is not None:
        cases.append(("Metal mis-coordination", metal_case))

    low_identity = _select_case(
        (~post["is_ion"])
        & (post["active_site_identity"] <= 0.70)
        & (post["ff_lev_score"] >= 4.0),
        "ff_lev_score",
        ascending=False,
    )
    if low_identity is not None:
        cases.append(("Low-identity donor", low_identity))

    if not cases:
        print("[WARN] No suitable case-study placements found.")
        return

    fig = plt.figure(figsize=(17, 5))
    all_points = []

    for idx, (label, row) in enumerate(cases, start=1):
        target = row["uniprot_id"]
        pdb_code = row["pdb_code"]
        ligand_id = row["ligand_id"]
        ligand_name = row["ligand_name"]
        tcs = row["ff_tcs_score"]
        rmsd = row["ff_local_rmsd"]
        lev = row.get("ff_lev_score")
        identity = row.get("active_site_identity")

        protein_pdb = (
            DATA_DIR
            / "foldfusion_output"
            / "Results"
            / target
            / "AlphaFold"
            / f"AF-{target}-F1-model_v4_processed.pdb"
        )
        if not protein_pdb.exists():
            fallback = protein_pdb.with_name(f"AF-{target}-F1-model_v3_processed.pdb")
            if fallback.exists():
                protein_pdb = fallback

        ligand_sdf = (
            DATA_DIR
            / "foldfusion_output"
            / "Results"
            / target
            / "JamdaScorer"
            / pdb_code
            / f"{ligand_id}.sdf"
        )
        if not ligand_sdf.exists():
            ligand_sdf = (
                DATA_DIR
                / "foldfusion_output"
                / "Results"
                / target
                / "LigandExtractor"
                / pdb_code
                / f"{ligand_id}.sdf"
            )

        ligand_atoms, bonds = _read_sdf_with_bonds(ligand_sdf)
        if not ligand_atoms:
            continue

        neighbors = _compute_neighbor_residues(protein_pdb, ligand_atoms, cutoff=4.0)

        ax = fig.add_subplot(1, len(cases), idx, projection="3d")
        colors = {
            "C": "#4F4F4F",
            "N": "#1f77b4",
            "O": "#d62728",
            "S": "#bcbd22",
            "P": "#ff7f0e",
            "H": "#b0b0b0",
            "ZN": "#17becf",
            "MG": "#17becf",
            "CA": "#9467bd",
            "FE": "#8c564b",
        }

        ligand_coords = np.array(
            [[atom["x"], atom["y"], atom["z"]] for atom in ligand_atoms], dtype=float
        )

        for i, j in bonds:
            if 0 <= i < len(ligand_atoms) and 0 <= j < len(ligand_atoms):
                p1 = ligand_coords[i]
                p2 = ligand_coords[j]
                ax.plot(*zip(p1, p2), color="#444444", linewidth=1.0, alpha=0.8)

        for atom, coord in zip(ligand_atoms, ligand_coords):
            element = atom.get("element", "C").upper()
            color = colors.get(element, "#4F4F4F")
            ax.scatter(*coord, color=color, s=40, depthshade=True)

        neighbor_points = []
        for neighbor in neighbors:
            centroid = neighbor["centroid"]
            neighbor_points.append(centroid)
            ax.scatter(centroid[0], centroid[1], centroid[2], color="#ff9896", s=35)
            label_text = f"{neighbor['resname']}{neighbor['resno']}"
            ax.text(centroid[0], centroid[1], centroid[2], label_text, fontsize=8)

        points = (
            np.vstack([ligand_coords, neighbor_points])
            if neighbor_points
            else ligand_coords
        )
        all_points.append(points)
        _set_equal_aspect(ax, points)
        ax.view_init(elev=18, azim=40)
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.zaxis.set_ticks([])
        ax.zaxis.set_ticklabels([])
        subtitle = (
            f"{label}\n{ligand_name} | TCS={tcs:.2f} | RMSD={rmsd:.2f}"
            + (f" | LEV={lev:.2f}" if pd.notnull(lev) else "")
            + (f" | identity={identity:.2f}" if pd.notnull(identity) else "")
        )
        ax.set_title(subtitle, fontsize=10)

    if all_points:
        combined = np.vstack(all_points)
        for ax in fig.axes:
            if hasattr(ax, "name") and ax.name == "3d":
                _set_equal_aspect(ax, combined)

    plt.tight_layout()
    _save_figure(fig, "fig_failure_modes_case_studies")
    plt.close(fig)
    print("[OK] Failure mode case-study figure exported.")


def table_internal_success(df: pd.DataFrame):
    """Summarize placement counts that meet quality thresholds."""

    post = df[df["stage"] == "post-jamda"].copy()
    if post.empty:
        print("[SKIP] No post-JAMDA rows for success summary.")
        return

    total = len(post)
    total_targets = post["uniprot_id"].nunique()

    rmsd_label = f"Local RMSD $\\leq$ {LOCAL_RMSD_THRESHOLD:.2f}~\\si{{\\angstrom}}"
    tcs_label = f"TCS $\\leq$ {TCS_THRESHOLD:.3f} (organic only)"
    lev_label = f"LEV $\\leq$ {LEV_THRESHOLD:.2f}~\\si{{\\angstrom}}"
    combined_label = (
        f"Local RMSD $\\leq$ {LOCAL_RMSD_THRESHOLD:.2f}~\\si{{\\angstrom}} "
        f"and TCS $\\leq$ {TCS_THRESHOLD:.3f} (organic only)"
    )

    organic_mask = ~post["is_ion"]
    masks = {
        "All placements": pd.Series(True, index=post.index),
        rmsd_label: post["ff_local_rmsd"] <= LOCAL_RMSD_THRESHOLD,
        tcs_label: organic_mask & (post["ff_tcs_score"] <= TCS_THRESHOLD),
        lev_label: post["ff_lev_score"] <= LEV_THRESHOLD,
        combined_label: (
            organic_mask
            & (post["ff_local_rmsd"] <= LOCAL_RMSD_THRESHOLD)
            & (post["ff_tcs_score"] <= TCS_THRESHOLD)
        ),
    }

    rows = []
    for label, mask in masks.items():
        mask = mask.fillna(False)
        placements = int(mask.sum())
        target_coverage = int(post.loc[mask, "uniprot_id"].nunique())
        rows.append(
            {
                "Criterion": label,
                "Placements": placements,
                "Share of placements": placements / total if total else np.nan,
                "Targets covered": target_coverage,
                "Share of targets": target_coverage / total_targets
                if total_targets
                else np.nan,
            }
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(OUTPUT_DIR / "table_internal_success.csv", index=False)

    latex = summary.copy()
    latex["Share of placements"] = latex["Share of placements"].map(
        lambda x: f"{x * 100:.1f}\\%" if pd.notnull(x) else "--"
    )
    latex["Share of targets"] = latex["Share of targets"].map(
        lambda x: f"{x * 100:.1f}\\%" if pd.notnull(x) else "--"
    )
    latex = latex.rename(
        columns={
            "Share of placements": "Share of placements (%)",
            "Share of targets": "Share of targets (%)",
        }
    )

    _export_table(
        latex,
        "table_internal_success.tex",
        caption=(
            "Placement counts from this pipeline that satisfy quality thresholds "
            "after JAMDA refinement."
        ),
        label="tab:internal-success",
        short_caption="Pipeline success counts",
    )
    print("[OK] Internal success summary exported.")


def table_ff_af_summary(df: pd.DataFrame):
    """Summaries of overlapping FoldFusion/AlphaFill metrics."""

    post = df[df["stage"] == "post-jamda"].copy()
    matched = post.dropna(
        subset=["ff_local_rmsd", "ff_tcs_score", "af_local_rmsd", "af_tcs_score"]
    )
    if matched.empty:
        print("[SKIP] No overlapping FF/AF metrics for summary table.")
        return

    def _metric_summary(series: pd.Series) -> dict:
        series = pd.to_numeric(series, errors="coerce").dropna()
        if series.empty:
            return {
                "N": 0,
                "Mean": np.nan,
                "Median": np.nan,
                "Q25": np.nan,
                "Q75": np.nan,
            }
        return {
            "N": int(series.size),
            "Mean": series.mean(),
            "Median": series.median(),
            "Q25": series.quantile(0.25),
            "Q75": series.quantile(0.75),
        }

    rows = []
    for pipeline, prefix in [("This work", "ff"), ("AlphaFill", "af")]:
        for metric, pretty in [
            ("local_rmsd", r"Local RMSD (\si{\angstrom})"),
            ("tcs_score", "TCS"),
        ]:
            stats = _metric_summary(matched[f"{prefix}_{metric}"])
            stats.update({"Pipeline": pipeline, "Metric": pretty})
            rows.append(stats)

    summary = pd.DataFrame(rows)[
        ["Pipeline", "Metric", "N", "Mean", "Median", "Q25", "Q75"]
    ]
    summary.to_csv(OUTPUT_DIR / "table_ff_af_summary.csv", index=False)
    _export_table(
        summary,
        "table_ff_af_summary.tex",
        caption=(
            r"Summary statistics for metrics shared between this pipeline and "
            r"AlphaFill on the \num{{{n}}} overlapping placements where both "
            r"Local \gls{{rmsd}} and \gls{{tcs}} are defined.".format(
                n=int(summary["N"].max())
            )
        ),
        label="tab:ff-af-summary",
        short_caption="Pipeline vs. AlphaFill summary statistics",
    )
    print("[OK] FF/AF summary table exported.")


def table_ff_af_success(df: pd.DataFrame):
    """Compare success rates between FoldFusion and AlphaFill on matched placements."""

    post = df[df["stage"] == "post-jamda"].copy()
    matched = post.dropna(
        subset=["ff_local_rmsd", "ff_tcs_score", "af_local_rmsd", "af_tcs_score"]
    )
    if matched.empty:
        print("[SKIP] No overlapping FF/AF metrics for success comparison.")
        return

    total = len(matched)
    total_targets = matched["uniprot_id"].nunique()

    rows = []
    for pipeline, prefix in [("This work", "ff"), ("AlphaFill", "af")]:
        mask = (matched[f"{prefix}_local_rmsd"] <= LOCAL_RMSD_THRESHOLD) & (
            matched[f"{prefix}_tcs_score"] <= TCS_THRESHOLD
        )
        mask = mask.fillna(False)
        placements = int(mask.sum())
        targets = int(matched.loc[mask, "uniprot_id"].nunique())
        share = placements / total if total else np.nan
        share_ci = _wilson_ci(placements, total) if total else (np.nan, np.nan)
        rows.append(
            {
                "Pipeline": pipeline,
                "Placements": placements,
                "Share of placements": share,
                "Share_low": share_ci[0],
                "Share_high": share_ci[1],
                "Targets covered": targets,
                "Share of targets": targets / total_targets
                if total_targets
                else np.nan,
            }
        )

    table = pd.DataFrame(rows)
    table.to_csv(OUTPUT_DIR / "table_ff_af_success.csv", index=False)

    latex = table.copy()
    latex["Share of placements"] = latex.apply(
        lambda row: (
            f"{row['Share of placements'] * 100:.1f}\\%\\;["
            f" {row['Share_low'] * 100:.1f}, {row['Share_high'] * 100:.1f} ]"
        )
        if pd.notnull(row["Share of placements"])
        else "--",
        axis=1,
    )
    latex["Share of targets"] = latex["Share of targets"].map(
        lambda x: f"{x * 100:.1f}\\%" if pd.notnull(x) else "--"
    )
    latex = latex.drop(columns=["Share_low", "Share_high"])

    _export_table(
        latex,
        "table_ff_af_success.tex",
        caption=(
            r"Proportion of overlapping placements that satisfy quality thresholds "
            r"for this pipeline and AlphaFill (Wilson 95\% CIs in brackets)."
        ),
        label="tab:ff-af-success",
        short_caption="Pipeline vs. AlphaFill success rates",
    )
    print("[OK] FF/AF success comparison exported.")


def figure_ff_af_distributions(df: pd.DataFrame):
    """Visualize matched FoldFusion vs AlphaFill metric distributions."""

    post = df[df["stage"] == "post-jamda"].copy()
    matched = post.dropna(
        subset=["ff_local_rmsd", "ff_tcs_score", "af_local_rmsd", "af_tcs_score"]
    )
    if matched.empty:
        print("[SKIP] No overlapping FF/AF metrics for distribution figure.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=False)

    rmsd_long = matched[["ff_local_rmsd", "af_local_rmsd"]].rename(
        columns={"ff_local_rmsd": "This work", "af_local_rmsd": "AlphaFill"}
    )
    rmsd_long = rmsd_long.melt(var_name="Pipeline", value_name="Local RMSD (\u212b)")
    sns.violinplot(
        data=rmsd_long,
        x="Pipeline",
        y="Local RMSD (\u212b)",
        inner="quartile",
        cut=0,
        ax=axes[0],
    )
    axes[0].set(
        title="Local RMSD distribution", xlabel="", ylabel="Local RMSD (\u212b)"
    )

    tcs_long = matched[["ff_tcs_score", "af_tcs_score"]].rename(
        columns={"ff_tcs_score": "This work", "af_tcs_score": "AlphaFill"}
    )
    tcs_long = tcs_long.melt(var_name="Pipeline", value_name="TCS (dimensionless)")
    sns.violinplot(
        data=tcs_long,
        x="Pipeline",
        y="TCS (dimensionless)",
        inner="quartile",
        cut=0,
        ax=axes[1],
    )
    axes[1].set(title="TCS distribution", xlabel="", ylabel="TCS (dimensionless)")

    fig.suptitle("This work vs AlphaFill metric distributions (matched placements)")
    _save_figure(fig, "fig_ff_af_metric_distributions")
    plt.close(fig)
    print("[OK] FF/AF distribution figure exported.")


def write_manifest():
    """Create a lightweight manifest (CSV) listing the thesis-ready assets."""
    asset_specs = [
        (
            "table",
            THESIS_TABLES_DIR / "table_overall_summary.tex",
            "3.2 Overall performance",
            "Overall transplantation summary (post-JAMDA).",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_attrition_funnel.tex",
            "3.1 Evaluation setup",
            "Pipeline inclusion funnel with attrition and runtime summaries.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_stage_runtime_histograms.pdf",
            "3.1 Evaluation setup",
            "Per-stage wall-time distributions.",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_internal_success.tex",
            "3.2 Overall performance",
            "FoldFusion placement counts that satisfy quality thresholds after JAMDA.",
        ),
        (
            "table",
            THESIS_APPENDIX_TABLES_DIR / "table_appendix_ion_breakdown.tex",
            "Appendix A Reproducibility",
            "Ligand-class breakdown for post-JAMDAscorer placements (organic vs ions).",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_top20_transplanted_ligands.pdf",
            "3.3 Ligand composition",
            "Top 20 most frequently transplanted ligands.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_top20_transplanted_ligands_no_buffers.pdf",
            "3.3 Ligand composition",
            "Top 20 ligands after excluding buffer species.",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_global_quality_stats.tex",
            "3.4 Quality indicators",
            r"Global quality indicators for post-JAMDA placements with 95\% CIs.",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_cross_metric_correlations.tex",
            "3.4 Quality indicators",
            r"Cross-metric Pearson and Spearman correlations with 95\% CIs.",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_identity_spearman.tex",
            "3.4 Alignment quality",
            "Spearman correlation of active-site identity with key metrics.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_distribution_ff_local_rmsd_by_identity_bin.pdf",
            "3.4 Alignment quality",
            "Local RMSD by active-site identity bin.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_distribution_ff_tcs_score_by_identity_bin.pdf",
            "3.4 Alignment quality",
            "TCS by active-site identity bin.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_jamda_tcs_pre_vs_post.pdf",
            "3.5 Optimization effects",
            "Pre vs post TCS (JAMDA refinement).",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_jamda_delta_distributions.pdf",
            "3.5 Optimization effects",
            "Distribution of JAMDA-induced ΔTCS and ΔRMSD.",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_jamda_delta_summary.tex",
            "3.5 Optimization effects",
            r"Summary of JAMDA refinement deltas with 95\% CIs.",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_ff_af_summary.tex",
            "3.6 External comparison",
            "FoldFusion vs AlphaFill metric summary.",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_ff_af_success.tex",
            "3.6 External comparison",
            "Quality-threshold success rates for FoldFusion and AlphaFill.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_ff_af_metric_distributions.pdf",
            "3.6 External comparison",
            "FoldFusion vs AlphaFill metric distributions.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_correlation_local_rmsd_ff_vs_af.pdf",
            "3.6 External comparison",
            "FF vs AF correlation (Local RMSD).",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_correlation_tcs_score_ff_vs_af.pdf",
            "3.6 External comparison",
            "FF vs AF correlation (TCS).",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_ff_vs_af_correlations.tex",
            "3.6 External comparison",
            "Pearson correlations between FoldFusion and AlphaFill.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_corr_tcs_vs_lev.pdf",
            "3.7 Cross-metric relationships",
            "TCS vs LEV.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_corr_rmsd_vs_lev.pdf",
            "3.7 Cross-metric relationships",
            "Local RMSD vs LEV.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_corr_rmsd_vs_tcs.pdf",
            "3.7 Cross-metric relationships",
            "Local RMSD vs TCS.",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_cross_metric_correlations_ff_only.tex",
            "3.7 Cross-metric relationships",
            "Pearson and Spearman correlations among FoldFusion metrics.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_failure_modes_case_studies.pdf",
            "3.8 Failure modes",
            "Representative failure cases (collapsed loop, metal mismatch, low-identity donor).",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_identity_threshold_sweep.tex",
            "3.9 Ablations",
            "Identity threshold sweep on placement yield and quality.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_identity_threshold_sweep.pdf",
            "3.9 Ablations",
            "Effect of active-site identity threshold on yield and quality.",
        ),
        (
            "table",
            THESIS_TABLES_DIR / "table_ion_exclusion_summary.tex",
            "3.9 Ablations",
            "Impact of ion exclusion on TCS performance.",
        ),
        (
            "figure",
            THESIS_FIGURES_DIR / "fig_tcs_distribution_ion_exclusion.pdf",
            "3.9 Ablations",
            "TCS distributions with and without monoatomic ions.",
        ),
    ]

    items = [
        {
            "type": asset_type,
            "file": str(path.relative_to(REPO_ROOT)),
            "section": section,
            "caption": caption,
        }
        for asset_type, path, section, caption in asset_specs
    ]

    man = pd.DataFrame(items)
    man.to_csv(OUTPUT_DIR / "thesis_manifest.csv", index=False)
    print("[OK] Manifest written: thesis_manifest.csv")


def main():
    df = build_master()
    if df.empty:
        return
    # Only thesis-relevant outputs below
    config_meta = load_config_metadata()
    stage_df = parse_pipeline_log()
    table_attrition_funnel(df, stage_df, config_meta)
    figure_stage_runtimes(stage_df)
    table_overall_summary(df)
    table_internal_success(df)
    appendix_ion_breakdown_table(df)
    figure_top20_ligands(df)
    figure_top20_ligands_no_buffers(df)
    table_global_quality(df)
    figure_alignment_by_identity(df)
    figure_jamda_effect(df)
    table_ff_af_summary(df)
    table_ff_af_success(df)
    figure_ff_af_distributions(df)
    figure_ff_vs_af(df)
    figure_cross_metric(df)
    identity_threshold_ablation(df)
    ion_exclusion_analysis(df)
    failure_mode_case_studies(df)
    write_manifest()
    print("\nDone. Thesis-ready assets are under:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
