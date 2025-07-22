import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os

# Set up the plotting style
plt.style.use("seaborn-v0_8")
sns.set_palette("husl")

# Create output directory for visualizations
output_dir = "foldfusion_visualizations"
os.makedirs(output_dir, exist_ok=True)


def load_evaluation_data(json_file):
    """Load and parse the evaluation JSON data."""
    with open(json_file, "r") as f:
        return json.load(f)


def create_optimization_comparison(data):
    """Create grouped bar charts comparing pre-jamda vs post-jamda performance."""
    # Prepare data for plotting
    uniprot_ids = []
    ligand_ids = []
    pre_local_rmsd = []
    post_local_rmsd = []
    pre_tcs = []
    post_tcs = []

    for uniprot_id, pdb_data in data.items():
        for pdb_code, pdb_info in pdb_data.items():
            for ligand_id, ligand_data in pdb_info.items():
                if ligand_id not in ["all_atom_rmsd", "backbone_rmsd"]:
                    uniprot_ids.append(uniprot_id)
                    ligand_ids.append(f"{pdb_code}_{ligand_id}")
                    pre_local_rmsd.append(ligand_data["pre-jamda"]["local_rmsd"])
                    post_local_rmsd.append(ligand_data["post-jamda"]["local_rmsd"])
                    pre_tcs.append(ligand_data["pre-jamda"]["tcs"])
                    post_tcs.append(ligand_data["post-jamda"]["tcs"])

    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    # Local RMSD comparison
    x = np.arange(len(ligand_ids))
    width = 0.35

    ax1.bar(
        x - width / 2,
        pre_local_rmsd,
        width,
        label="Pre-JAMDA",
        alpha=0.8,
        color="lightcoral",
    )
    ax1.bar(
        x + width / 2,
        post_local_rmsd,
        width,
        label="Post-JAMDA",
        alpha=0.8,
        color="skyblue",
    )

    ax1.set_xlabel("Ligand Transplants")
    ax1.set_ylabel("Local RMSD (Å)")
    ax1.set_title("Local RMSD: Pre-JAMDA vs Post-JAMDA Optimization")
    ax1.set_xticks(x)
    ax1.set_xticklabels(ligand_ids, rotation=45, ha="right")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # TCS comparison
    ax2.bar(
        x - width / 2, pre_tcs, width, label="Pre-JAMDA", alpha=0.8, color="lightcoral"
    )
    ax2.bar(
        x + width / 2, post_tcs, width, label="Post-JAMDA", alpha=0.8, color="skyblue"
    )

    ax2.set_xlabel("Ligand Transplants")
    ax2.set_ylabel("Transplant Clash Score (TCS)")
    ax2.set_title("TCS: Pre-JAMDA vs Post-JAMDA Optimization")
    ax2.set_xticks(x)
    ax2.set_xticklabels(ligand_ids, rotation=45, ha="right")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(
        f"{output_dir}/optimization_comparison.png", dpi=300, bbox_inches="tight"
    )
    plt.close()

    return fig


def create_alignment_quality_plots(data):
    """Create scatter plots showing alignment quality vs transplant quality."""
    # Prepare data
    backbone_rmsd_values = []
    all_atom_rmsd_values = []
    post_local_rmsd_values = []
    pdb_codes = []
    uniprot_ids = []

    for uniprot_id, pdb_data in data.items():
        for pdb_code, pdb_info in pdb_data.items():
            backbone_rmsd = pdb_info["backbone_rmsd"]
            all_atom_rmsd = pdb_info["all_atom_rmsd"]

            # Get average post-jamda local_rmsd for this PDB
            local_rmsd_values = []
            for ligand_id, ligand_data in pdb_info.items():
                if ligand_id not in ["all_atom_rmsd", "backbone_rmsd"]:
                    local_rmsd_values.append(ligand_data["post-jamda"]["local_rmsd"])

            if local_rmsd_values:
                avg_local_rmsd = np.mean(local_rmsd_values)
                backbone_rmsd_values.append(backbone_rmsd)
                all_atom_rmsd_values.append(all_atom_rmsd)
                post_local_rmsd_values.append(avg_local_rmsd)
                pdb_codes.append(pdb_code)
                uniprot_ids.append(uniprot_id)

    # Create scatter plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Plot 1: Backbone RMSD vs Post-JAMDA Local RMSD
    ax1.scatter(
        backbone_rmsd_values,
        post_local_rmsd_values,
        alpha=0.7,
        s=100,
        c=range(len(backbone_rmsd_values)),
        cmap="viridis",
    )
    ax1.set_xlabel("Backbone RMSD (Å)")
    ax1.set_ylabel("Post-JAMDA Local RMSD (Å)")
    ax1.set_title("Protein Alignment Quality vs Ligand Fit\n(Backbone RMSD)")
    ax1.grid(True, alpha=0.3)

    # Add trend line
    z1 = np.polyfit(backbone_rmsd_values, post_local_rmsd_values, 1)
    p1 = np.poly1d(z1)
    ax1.plot(backbone_rmsd_values, p1(backbone_rmsd_values), "r--", alpha=0.8)

    # Calculate correlation
    corr1 = np.corrcoef(backbone_rmsd_values, post_local_rmsd_values)[0, 1]
    ax1.text(
        0.05,
        0.95,
        f"R = {corr1:.3f}",
        transform=ax1.transAxes,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    # Plot 2: All-atom RMSD vs Post-JAMDA Local RMSD
    ax2.scatter(
        all_atom_rmsd_values,
        post_local_rmsd_values,
        alpha=0.7,
        s=100,
        c=range(len(all_atom_rmsd_values)),
        cmap="viridis",
    )
    ax2.set_xlabel("All-atom RMSD (Å)")
    ax2.set_ylabel("Post-JAMDA Local RMSD (Å)")
    ax2.set_title("Protein Alignment Quality vs Ligand Fit\n(All-atom RMSD)")
    ax2.grid(True, alpha=0.3)

    # Add trend line
    z2 = np.polyfit(all_atom_rmsd_values, post_local_rmsd_values, 1)
    p2 = np.poly1d(z2)
    ax2.plot(all_atom_rmsd_values, p2(all_atom_rmsd_values), "r--", alpha=0.8)

    # Calculate correlation
    corr2 = np.corrcoef(all_atom_rmsd_values, post_local_rmsd_values)[0, 1]
    ax2.text(
        0.05,
        0.95,
        f"R = {corr2:.3f}",
        transform=ax2.transAxes,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()
    plt.savefig(
        f"{output_dir}/alignment_quality_analysis.png", dpi=300, bbox_inches="tight"
    )
    plt.close()

    return fig, (corr1, corr2)


def create_improvement_analysis(data):
    """Analyze and visualize the improvement from optimization."""
    improvements_data = []

    for uniprot_id, pdb_data in data.items():
        for pdb_code, pdb_info in pdb_data.items():
            for ligand_id, ligand_data in pdb_info.items():
                if ligand_id not in ["all_atom_rmsd", "backbone_rmsd"]:
                    pre_rmsd = ligand_data["pre-jamda"]["local_rmsd"]
                    post_rmsd = ligand_data["post-jamda"]["local_rmsd"]
                    pre_tcs = ligand_data["pre-jamda"]["tcs"]
                    post_tcs = ligand_data["post-jamda"]["tcs"]

                    rmsd_improvement = pre_rmsd - post_rmsd
                    tcs_improvement = pre_tcs - post_tcs
                    rmsd_improvement_pct = (rmsd_improvement / pre_rmsd) * 100
                    tcs_improvement_pct = (tcs_improvement / pre_tcs) * 100

                    improvements_data.append(
                        {
                            "uniprot_id": uniprot_id,
                            "pdb_ligand": f"{pdb_code}_{ligand_id}",
                            "rmsd_improvement": rmsd_improvement,
                            "tcs_improvement": tcs_improvement,
                            "rmsd_improvement_pct": rmsd_improvement_pct,
                            "tcs_improvement_pct": tcs_improvement_pct,
                        }
                    )

    df = pd.DataFrame(improvements_data)

    # Create visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # RMSD improvement histogram
    ax1.hist(
        df["rmsd_improvement"], bins=20, alpha=0.7, color="skyblue", edgecolor="black"
    )
    ax1.axvline(x=0, color="red", linestyle="--", alpha=0.8)
    ax1.set_xlabel("RMSD Improvement (Å)")
    ax1.set_ylabel("Frequency")
    ax1.set_title("Distribution of RMSD Improvements")
    ax1.grid(True, alpha=0.3)

    # TCS improvement histogram
    ax2.hist(
        df["tcs_improvement"], bins=20, alpha=0.7, color="lightcoral", edgecolor="black"
    )
    ax2.axvline(x=0, color="red", linestyle="--", alpha=0.8)
    ax2.set_xlabel("TCS Improvement")
    ax2.set_ylabel("Frequency")
    ax2.set_title("Distribution of TCS Improvements")
    ax2.grid(True, alpha=0.3)

    # RMSD vs TCS improvement scatter
    ax3.scatter(df["rmsd_improvement"], df["tcs_improvement"], alpha=0.7, s=60)
    ax3.axhline(y=0, color="red", linestyle="--", alpha=0.8)
    ax3.axvline(x=0, color="red", linestyle="--", alpha=0.8)
    ax3.set_xlabel("RMSD Improvement (Å)")
    ax3.set_ylabel("TCS Improvement")
    ax3.set_title("RMSD vs TCS Improvement Correlation")
    ax3.grid(True, alpha=0.3)

    # Percentage improvements
    ax4.scatter(
        df["rmsd_improvement_pct"],
        df["tcs_improvement_pct"],
        alpha=0.7,
        s=60,
        c="green",
    )
    ax4.axhline(y=0, color="red", linestyle="--", alpha=0.8)
    ax4.axvline(x=0, color="red", linestyle="--", alpha=0.8)
    ax4.set_xlabel("RMSD Improvement (%)")
    ax4.set_ylabel("TCS Improvement (%)")
    ax4.set_title("Percentage Improvements Correlation")
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/improvement_analysis.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Print summary statistics
    print("\n" + "=" * 60)
    print("OPTIMIZATION IMPROVEMENT ANALYSIS")
    print("=" * 60)
    print(f"Total transplants analyzed: {len(df)}")
    print(
        f"RMSD improvements > 0: {(df['rmsd_improvement'] > 0).sum()} ({(df['rmsd_improvement'] > 0).mean() * 100:.1f}%)"
    )
    print(
        f"TCS improvements > 0: {(df['tcs_improvement'] > 0).sum()} ({(df['tcs_improvement'] > 0).mean() * 100:.1f}%)"
    )
    print(
        f"Mean RMSD improvement: {df['rmsd_improvement'].mean():.3f} ± {df['rmsd_improvement'].std():.3f} Å"
    )
    print(
        f"Mean TCS improvement: {df['tcs_improvement'].mean():.3f} ± {df['tcs_improvement'].std():.3f}"
    )

    return df


def create_optimization_effectiveness_plot(data):
    """Create a detailed analysis of optimization effectiveness."""
    # Prepare data
    effectiveness_data = []

    for uniprot_id, pdb_data in data.items():
        for pdb_code, pdb_info in pdb_data.items():
            for ligand_id, ligand_data in pdb_info.items():
                if ligand_id not in ["all_atom_rmsd", "backbone_rmsd"]:
                    pre_rmsd = ligand_data["pre-jamda"]["local_rmsd"]
                    post_rmsd = ligand_data["post-jamda"]["local_rmsd"]
                    pre_tcs = ligand_data["pre-jamda"]["tcs"]
                    post_tcs = ligand_data["post-jamda"]["tcs"]

                    effectiveness_data.append(
                        {
                            "uniprot_id": uniprot_id,
                            "pdb_ligand": f"{pdb_code}_{ligand_id}",
                            "pre_rmsd": pre_rmsd,
                            "post_rmsd": post_rmsd,
                            "pre_tcs": pre_tcs,
                            "post_tcs": post_tcs,
                            "rmsd_change": post_rmsd - pre_rmsd,
                            "tcs_change": post_tcs - pre_tcs,
                            "both_improved": (post_rmsd < pre_rmsd)
                            and (post_tcs < pre_tcs),
                            "rmsd_improved": post_rmsd < pre_rmsd,
                            "tcs_improved": post_tcs < pre_tcs,
                        }
                    )

    df = pd.DataFrame(effectiveness_data)

    # Create the plot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Before vs After scatter for RMSD
    ax1.scatter(df["pre_rmsd"], df["post_rmsd"], alpha=0.7, s=60)
    min_val = min(df["pre_rmsd"].min(), df["post_rmsd"].min())
    max_val = max(df["pre_rmsd"].max(), df["post_rmsd"].max())
    ax1.plot(
        [min_val, max_val], [min_val, max_val], "r--", alpha=0.8, label="No change line"
    )
    ax1.set_xlabel("Pre-JAMDA Local RMSD (Å)")
    ax1.set_ylabel("Post-JAMDA Local RMSD (Å)")
    ax1.set_title("RMSD: Before vs After Optimization")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Before vs After scatter for TCS
    ax2.scatter(df["pre_tcs"], df["post_tcs"], alpha=0.7, s=60, color="orange")
    min_val = min(df["pre_tcs"].min(), df["post_tcs"].min())
    max_val = max(df["pre_tcs"].max(), df["post_tcs"].max())
    ax2.plot(
        [min_val, max_val], [min_val, max_val], "r--", alpha=0.8, label="No change line"
    )
    ax2.set_xlabel("Pre-JAMDA TCS")
    ax2.set_ylabel("Post-JAMDA TCS")
    ax2.set_title("TCS: Before vs After Optimization")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Plot 3: Improvement categories
    categories = ["Both Improved", "RMSD Only", "TCS Only", "Both Worsened"]
    counts = [
        df["both_improved"].sum(),
        (df["rmsd_improved"] & ~df["tcs_improved"]).sum(),
        (~df["rmsd_improved"] & df["tcs_improved"]).sum(),
        (~df["rmsd_improved"] & ~df["tcs_improved"]).sum(),
    ]

    colors = ["green", "lightblue", "orange", "red"]
    ax3.pie(counts, labels=categories, colors=colors, autopct="%1.1f%%", startangle=90)
    ax3.set_title("Optimization Outcomes Distribution")

    # Plot 4: Change magnitude distribution
    ax4.hist(
        df["rmsd_change"], bins=20, alpha=0.7, label="RMSD Change", color="skyblue"
    )
    ax4.axvline(x=0, color="red", linestyle="--", alpha=0.8)
    ax4.set_xlabel("Change in Local RMSD (Å)")
    ax4.set_ylabel("Frequency")
    ax4.set_title("Distribution of RMSD Changes")
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(
        f"{output_dir}/optimization_effectiveness.png", dpi=300, bbox_inches="tight"
    )
    plt.close()

    return df


def main():
    """Main function to generate selected visualizations."""
    # Load data
    data = load_evaluation_data("evaluation.json")

    print("FoldFusion Evaluation Analysis (Refactored)")
    print("=" * 50)
    print(f"Analyzing data for {len(data)} UniProt IDs")
    print(f"Saving visualizations to: {output_dir}/")

    # Generate selected visualizations only
    print("\n1. Creating alignment quality analysis...")
    fig, correlations = create_alignment_quality_plots(data)
    print(f"   Backbone RMSD correlation: {correlations[0]:.3f}")
    print(f"   All-atom RMSD correlation: {correlations[1]:.3f}")

    print("\n2. Creating improvement analysis...")
    create_improvement_analysis(data)

    print("\n3. Creating optimization comparison plots...")
    create_optimization_comparison(data)

    print("\n4. Creating optimization effectiveness analysis...")
    create_optimization_effectiveness_plot(data)

    print("\nAnalysis complete! Generated visualizations:")
    print(f"- {output_dir}/alignment_quality_analysis.png")
    print(f"- {output_dir}/improvement_analysis.png")
    print(f"- {output_dir}/optimization_comparison.png")
    print(f"- {output_dir}/optimization_effectiveness.png")


if __name__ == "__main__":
    main()
