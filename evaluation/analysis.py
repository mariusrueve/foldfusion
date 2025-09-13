import json
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn


def validate_file(path: Path, extension: str):
    if not path.exists():
        raise FileNotFoundError(f"File {path} does not exist.")
    if not path.is_file():
        raise ValueError(f"Path {path} is not a file.")
    if path.suffix != extension:
        raise ValueError(
            f"File {path} does not have the expected {extension} extension."
        )


def validate_folder(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Folder {path} does not exist.")
    if not path.is_dir():
        raise ValueError(f"Path {path} is not a folder.")


def load_bechmark_data(path: Path) -> pd.DataFrame:
    with open(path) as f:
        data = json.load(f)

    ligand_bins = data.get("ligand_bins")
    if ligand_bins is None:
        raise KeyError("Expected key 'ligand_bins' in benchmark_dataset.json")

    rows = data["rows"]

    return pd.DataFrame(rows)


def load_foldfusion_evaluation(path: Path) -> pd.DataFrame:
    with open(path) as f:
        eval_data = json.load(f)

    records = []
    for uniprot_id, pdb_map in eval_data.items():
        if not isinstance(pdb_map, dict):
            continue
        for pdb_code, metrics in pdb_map.items():
            if not isinstance(metrics, dict):
                continue
            all_atom = metrics.get("all_atom_rmsd")
            backbone = metrics.get("backbone_rmsd")
            ligand_keys = [
                k for k in metrics.keys() if k not in {"all_atom_rmsd", "backbone_rmsd"}
            ]
            deltas = []
            pre_vals = []
            post_vals = []
            for lk in ligand_keys:
                lig_block = metrics.get(lk, {})
                pre = lig_block.get("pre-jamda", {})
                post = lig_block.get("post-jamda", {})
                pre_lrmsd = pre.get("local_rmsd")
                post_lrmsd = post.get("local_rmsd")
                if isinstance(pre_lrmsd, int | float) and isinstance(
                    post_lrmsd, int | float
                ):
                    deltas.append(pre_lrmsd - post_lrmsd)
                    pre_vals.append(pre_lrmsd)
                    post_vals.append(post_lrmsd)

            foldfusion_score = float(pd.Series(deltas).mean()) if deltas else None
            pre_mean = float(pd.Series(pre_vals).mean()) if pre_vals else None
            post_mean = float(pd.Series(post_vals).mean()) if post_vals else None
            pre_min = float(pd.Series(pre_vals).min()) if pre_vals else None
            post_min = float(pd.Series(post_vals).min()) if post_vals else None

            records.append(
                {
                    "uniprot_id": uniprot_id,
                    "validation_pdb": pdb_code,
                    "all_atom_rmsd": all_atom,
                    "backbone_rmsd": backbone,
                    "ligand_measurements": len(deltas),
                    "foldfusion": foldfusion_score,
                    "pre_local_rmsd_mean": pre_mean,
                    "post_local_rmsd_mean": post_mean,
                    "pre_local_rmsd_min": pre_min,
                    "post_local_rmsd_min": post_min,
                }
            )

    return pd.DataFrame(records)


def extract_alphafill_metrics(alphafill_path: Path, matched_ligand: str | None) -> dict:
    try:
        with open(alphafill_path) as f:
            af = json.load(f)
    except FileNotFoundError:
        return {
            "alphafill_local_rmsd_mean": None,
            "alphafill_local_rmsd_min": None,
            "alphafill_transplants": 0,
        }
    except Exception:
        return {
            "alphafill_local_rmsd_mean": None,
            "alphafill_local_rmsd_min": None,
            "alphafill_transplants": 0,
        }

    target = (matched_ligand or "").strip().upper()
    vals = []
    hits = af.get("hits") or []
    for hit in hits:
        transplants = (hit.get("transplants") or []) if isinstance(hit, dict) else []
        for tr in transplants:
            comp = str(tr.get("compound_id", "")).strip().upper()
            if target and comp != target:
                continue
            lr = tr.get("local_rmsd")
            if isinstance(lr, int | float):
                vals.append(lr)

    if not vals:
        return {
            "alphafill_local_rmsd_mean": None,
            "alphafill_local_rmsd_min": None,
            "alphafill_transplants": 0,
        }

    s = pd.Series(vals)
    return {
        "alphafill_local_rmsd_mean": float(s.mean()),
        "alphafill_local_rmsd_min": float(s.min()),
        "alphafill_transplants": int(len(vals)),
    }


def load_per_ligand_evaluation(path: Path) -> pd.DataFrame:
    """Flatten evaluation JSON to per-ligand rows.

    Returns columns: uniprot_id, validation_pdb, comp_id, pre_local_rmsd,
    post_local_rmsd, delta_local_rmsd.
    """
    with open(path) as f:
        eval_data = json.load(f)

    recs = []
    for uniprot_id, pdb_map in eval_data.items():
        if not isinstance(pdb_map, dict):
            continue
        for pdb_code, metrics in pdb_map.items():
            if not isinstance(metrics, dict):
                continue
            for k, lig_block in metrics.items():
                if k in {"all_atom_rmsd", "backbone_rmsd"}:
                    continue
                if not isinstance(lig_block, dict):
                    continue
                comp_id = k.split("_")[0].strip().upper()
                pre = lig_block.get("pre-jamda") or {}
                post = lig_block.get("post-jamda") or {}
                pre_lrmsd = pre.get("local_rmsd")
                post_lrmsd = post.get("local_rmsd")
                if isinstance(pre_lrmsd, int | float) and isinstance(
                    post_lrmsd, int | float
                ):
                    recs.append(
                        {
                            "uniprot_id": uniprot_id,
                            "validation_pdb": pdb_code,
                            "comp_id": comp_id,
                            "pre_local_rmsd": float(pre_lrmsd),
                            "post_local_rmsd": float(post_lrmsd),
                            "delta_local_rmsd": float(pre_lrmsd - post_lrmsd),
                        }
                    )
    return pd.DataFrame(recs)


def bootstrap_ci(
    values: Iterable[float],
    n_boot: int = 2000,
    ci: float = 0.95,
    random_state: int = 42,
) -> tuple[float, float]:
    s = pd.Series(list(values)).dropna()
    if s.empty:
        return (float("nan"), float("nan"))
    rs = pd.Series(s.values)
    rng = pd.Index(range(len(rs)))
    q_low = (1 - ci) / 2
    q_high = 1 - q_low
    boots = []
    for i in range(n_boot):
        idx = rng.to_series().sample(
            n=len(rs), replace=True, random_state=random_state + i
        )
        boots.append(float(rs.iloc[idx.to_numpy()].median()))
    return (
        float(pd.Series(boots).quantile(q_low)),
        float(pd.Series(boots).quantile(q_high)),
    )


def cluster_bootstrap_median(
    df: pd.DataFrame,
    cluster_col: str,
    value_col: str,
    n_boot: int = 2000,
    ci: float = 0.95,
    random_state: int = 1337,
) -> tuple[float, float]:
    clusters = df[[cluster_col]].drop_duplicates()[cluster_col].tolist()
    if not clusters:
        return (float("nan"), float("nan"))
    q_low = (1 - ci) / 2
    q_high = 1 - q_low
    boots = []
    for i in range(n_boot):
        sample_clusters = pd.Series(clusters).sample(
            n=len(clusters), replace=True, random_state=random_state + i
        )
        samp_df = pd.concat(
            [df[df[cluster_col] == c] for c in sample_clusters], ignore_index=True
        )
        boots.append(float(samp_df[value_col].median()))
    b = pd.Series(boots)
    return (float(b.quantile(q_low)), float(b.quantile(q_high)))


def main():
    data_path = Path("evaluation/data/")
    benchmark_dataset_path = data_path / "benchmark_dataset.json"
    foldfusion_results_path = data_path / "foldfusion_output/Evaluation/evaluation.json"
    alphafill_folder_path = data_path / "benchmark_dataset_meta/alphafill"
    output_path = Path("evaluation/analysis_output")
    output_path.mkdir(parents=True, exist_ok=True)

    validate_file(benchmark_dataset_path, ".json")
    validate_file(foldfusion_results_path, ".json")
    validate_folder(alphafill_folder_path)
    validate_folder(output_path)

    plot_folder_path = output_path / "plots/"
    plot_folder_path.mkdir(parents=True, exist_ok=True)

    # Load structured data
    benchmark_df = load_bechmark_data(benchmark_dataset_path)
    evaluation_df = load_foldfusion_evaluation(foldfusion_results_path)

    # Merge and filter
    merged = pd.merge(
        benchmark_df,
        evaluation_df,
        how="left",
        on=["uniprot_id", "validation_pdb"],
        validate="m:1",
    )
    # Compute AlphaFill metrics per row using matched ligand
    af_rows = []
    for _, r in benchmark_df.iterrows():
        af_metrics = extract_alphafill_metrics(
            Path(r["alphafill_json"]), r.get("matched_ligand")
        )
        af_rows.append(
            {
                "uniprot_id": r["uniprot_id"],
                "matched_ligand": r.get("matched_ligand"),
                **af_metrics,
            }
        )
    alphafill_df = pd.DataFrame(af_rows)

    merged = pd.merge(
        merged,
        alphafill_df,
        how="left",
        on=["uniprot_id", "matched_ligand"],
    )

    # Retain merged if needed, but focus on per-ligand analyses below

    # Per-ligand evaluation joined with benchmark rows
    per_lig = load_per_ligand_evaluation(foldfusion_results_path)
    # Normalize case for matching
    per_lig["comp_id"] = per_lig["comp_id"].str.upper()
    bench_for_join = benchmark_df.copy()
    bench_for_join["matched_ligand"] = bench_for_join["matched_ligand"].str.upper()
    lig_merged = pd.merge(
        per_lig,
        bench_for_join[
            [
                "uniprot_id",
                "validation_pdb",
                "matched_ligand",
                "bin_type",
                "bin_name",
                "alphafill_json",
            ]
        ],
        left_on=["uniprot_id", "validation_pdb", "comp_id"],
        right_on=["uniprot_id", "validation_pdb", "matched_ligand"],
        how="inner",
    )
    # Attach AlphaFill metrics
    lig_merged = pd.merge(
        lig_merged,
        alphafill_df,
        how="left",
        on=["uniprot_id", "matched_ligand"],
    )

    # Save per-ligand results (primary analysis focus)
    per_lig_out = output_path / "per_ligand_results.csv"
    lig_merged.to_csv(per_lig_out, index=False)
    seaborn.set_theme(style="whitegrid")
    # 1) Paired violin: pre vs post local RMSD (per ligand)
    melted = lig_merged.melt(
        id_vars=["uniprot_id", "validation_pdb", "bin_name"],
        value_vars=["pre_local_rmsd", "post_local_rmsd"],
        var_name="phase",
        value_name="local_rmsd",
    ).dropna(subset=["local_rmsd"])
    plt.figure(figsize=(6, 5))
    seaborn.violinplot(data=melted, x="phase", y="local_rmsd", cut=0)
    plt.ylabel("Local RMSD (Å)")
    plt.tight_layout()
    plt.savefig(plot_folder_path / "paired_violin_pre_post.png", dpi=200)

    # 2) CDF of Δ = pre - post
    dvals = lig_merged["delta_local_rmsd"].dropna().to_numpy()
    dvals.sort()
    cdf = (np.arange(1, len(dvals) + 1) / len(dvals)) if len(dvals) else np.array([])
    plt.figure(figsize=(6, 5))
    if len(dvals):
        plt.plot(dvals, cdf, label="Δ CDF")
        plt.axvline(0.0, color="grey", linestyle="--", linewidth=1)
    plt.xlabel("Δ local RMSD (pre - post) [Å]")
    plt.ylabel("Cumulative fraction")
    plt.tight_layout()
    plt.savefig(plot_folder_path / "delta_cdf.png", dpi=200)

    # 3) Δ by bin
    plt.figure(figsize=(10, 6))
    seaborn.boxplot(data=lig_merged, x="bin_name", y="delta_local_rmsd")
    plt.xticks(rotation=45, ha="right")
    plt.xlabel("Ligand bin")
    plt.ylabel("Δ local RMSD (pre - post) [Å]")
    plt.tight_layout()
    plt.savefig(plot_folder_path / "delta_by_bin.png", dpi=200)

    # 4) AlphaFill (min) vs Δ correlation (overall + by bin)
    lig_corr = lig_merged.dropna(
        subset=["delta_local_rmsd", "alphafill_local_rmsd_min"]
    ).copy()
    lig_corr_pearson = (
        lig_corr["delta_local_rmsd"].corr(
            lig_corr["alphafill_local_rmsd_min"], method="pearson"
        )
        if not lig_corr.empty
        else None
    )
    lig_corr_spearman = (
        lig_corr["delta_local_rmsd"].corr(
            lig_corr["alphafill_local_rmsd_min"], method="spearman"
        )
        if not lig_corr.empty
        else None
    )
    pd.DataFrame(
        {
            "metric": ["pearson", "spearman"],
            "value": [lig_corr_pearson, lig_corr_spearman],
            "n": [len(lig_corr), len(lig_corr)],
        }
    ).to_csv(output_path / "alphafill_correlation_per_ligand.csv", index=False)

    # 4a) Overall scatter with regression
    if not lig_corr.empty:
        plt.figure(figsize=(6, 5))
        seaborn.regplot(
            data=lig_corr,
            x="alphafill_local_rmsd_min",
            y="delta_local_rmsd",
            scatter_kws={"s": 16, "alpha": 0.6},
            line_kws={"color": "red"},
        )
        plt.xlabel("AlphaFill local RMSD (min) [Å]")
        plt.ylabel("Δ local RMSD (pre - post) [Å]")
        plt.tight_layout()
        plt.savefig(plot_folder_path / "alphafill_min_vs_delta.png", dpi=200)

        # 4b) Faceted by bin_name
        g = seaborn.lmplot(
            data=lig_corr,
            x="alphafill_local_rmsd_min",
            y="delta_local_rmsd",
            col="bin_name",
            col_wrap=4,
            scatter_kws={"s": 12, "alpha": 0.5},
            line_kws={"color": "red"},
            height=3,
        )
        g.set_axis_labels("AlphaFill local RMSD (min) [Å]", "Δ local RMSD [Å]")
        plt.tight_layout()
        plt.savefig(plot_folder_path / "alphafill_min_vs_delta_by_bin.png", dpi=200)

    # Summaries for thesis
    # Overall improvement stats (per-ligand)
    if not lig_merged.empty:
        overall = {
            "n": int(len(lig_merged)),
            "median_delta": float(lig_merged["delta_local_rmsd"].median()),
            "mean_delta": float(lig_merged["delta_local_rmsd"].mean()),
            ">0_frac": float((lig_merged["delta_local_rmsd"] > 0).mean()),
            ">0.25_frac": float((lig_merged["delta_local_rmsd"] > 0.25).mean()),
            ">0.5_frac": float((lig_merged["delta_local_rmsd"] > 0.5).mean()),
        }
        lo, hi = bootstrap_ci(lig_merged["delta_local_rmsd"], n_boot=2000)
        overall["median_ci_low"] = lo
        overall["median_ci_high"] = hi
        pd.DataFrame([overall]).to_csv(
            output_path / "summary_overall_per_ligand.csv", index=False
        )

        # By bin_name
        grp = lig_merged.groupby("bin_name")["delta_local_rmsd"].agg(
            ["count", "median", "mean"]
        )
        grp = grp.rename(columns={"count": "n"})
        grp.to_csv(output_path / "summary_by_bin_per_ligand.csv")


if __name__ == "__main__":
    main()
