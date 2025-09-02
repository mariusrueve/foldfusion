from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import plotly.offline as plotly_offline
from scipy.spatial import distance_matrix

from .utils import get_coordinates, parse_pdb, parse_sdf


def create_sphere_mesh(center: np.ndarray, radius: float, resolution: int = 20):
    """Create a sphere mesh for visualization."""

    # Create sphere coordinates
    phi = np.linspace(0, 2 * np.pi, resolution)
    theta = np.linspace(0, np.pi, resolution)
    phi, theta = np.meshgrid(phi, theta)

    x = center[0] + radius * np.sin(theta) * np.cos(phi)
    y = center[1] + radius * np.sin(theta) * np.sin(phi)
    z = center[2] + radius * np.cos(theta)

    return x, y, z


def visualize_protein_with_radius(
    alphafold_protein: Path,
    experimental_protein: Path,
    ligand: Path,
    radius: float = 6.0,
    output_file: Path | None = None,
) -> tuple[object | None, str]:
    """
    Create an interactive 3D visualization of proteins with ligand radius sphere.

    Args:
        alphafold_protein: Path to AlphaFold protein PDB file
        experimental_protein: Path to experimental protein PDB file
        ligand: Path to ligand SDF file
        radius: Radius around ligand to visualize
        output_file: Optional path to save HTML. If not provided and config is
            available, saves to config output directory with auto-generated
            filename

    Returns:
        Tuple of (plotly figure, HTML string)
    """

    # Parse structures
    alphafold_coords = get_coordinates(parse_pdb(alphafold_protein))
    experimental_coords = get_coordinates(parse_pdb(experimental_protein))
    ligand_coords = get_coordinates(parse_sdf(ligand))

    # Calculate ligand center for sphere
    ligand_center = ligand_coords.mean(axis=0)

    # Create sphere mesh
    sphere_x, sphere_y, sphere_z = create_sphere_mesh(ligand_center, radius)

    # Create the plot
    fig = go.Figure()

    # Add AlphaFold protein points
    fig.add_trace(
        go.Scatter3d(
            x=alphafold_coords[:, 0],
            y=alphafold_coords[:, 1],
            z=alphafold_coords[:, 2],
            mode="markers",
            marker={"size": 3, "color": "blue", "opacity": 0.7},
            name="AlphaFold Protein",
            hovertemplate=(
                "AlphaFold"
                "<br>X: %{x}<br>Y: %{y}<br>Z: %{z}<extra></extra>"
            ),
        )
    )

    # Add experimental protein points
    fig.add_trace(
        go.Scatter3d(
            x=experimental_coords[:, 0],
            y=experimental_coords[:, 1],
            z=experimental_coords[:, 2],
            mode="markers",
            marker={"size": 3, "color": "red", "opacity": 0.7},
            name="Experimental Protein",
            hovertemplate=(
                "Experimental"
                "<br>X: %{x}<br>Y: %{y}<br>Z: %{z}<extra></extra>"
            ),
        )
    )

    # Add ligand points
    fig.add_trace(
        go.Scatter3d(
            x=ligand_coords[:, 0],
            y=ligand_coords[:, 1],
            z=ligand_coords[:, 2],
            mode="markers",
            marker={"size": 8, "color": "green", "opacity": 1.0},
            name="Ligand",
            hovertemplate=(
                "Ligand" "<br>X: %{x}<br>Y: %{y}<br>Z: %{z}<extra></extra>"
            ),
        )
    )

    # Add transparent sphere around ligand
    fig.add_trace(
        go.Surface(
            x=sphere_x,
            y=sphere_y,
            z=sphere_z,
            opacity=0.2,
            colorscale=[[0, "yellow"], [1, "yellow"]],
            showscale=False,
            name=f"Radius Sphere ({radius}Å)",
            hovertemplate=f"Radius: {radius}Å<extra></extra>",
        )
    )

    # Update layout
    fig.update_layout(
        title=f"Protein Structure Comparison with {radius}Å Radius",
        scene={
            "xaxis_title": "X (Å)",
            "yaxis_title": "Y (Å)",
            "zaxis_title": "Z (Å)",
            "aspectmode": "cube",
        },
        width=1000,
        height=800,
        showlegend=True,
    )

    # Determine output file path
    save_file = None
    if output_file:
        save_file = Path(output_file)
    # Save to file if path is determined
    if save_file:
        if save_file.suffix != ".html":
            save_file = save_file.with_suffix(".html")
        fig.write_html(str(save_file))
        print(f"Visualization saved to: {save_file}")

    # Generate HTML string
    plot_html = plotly_offline.plot(fig, include_plotlyjs=True, output_type="div")

    return fig, plot_html


def visualize_local_rmsd_region(
    alphafold_protein: Path,
    experimental_protein: Path,
    ligand: Path,
    radius: float = 6.0,
    rmsd_value: float | None = None,
    output_file: Path | None = None,
):
    """
    Create visualization highlighting the atoms used in local RMSD calculation.

    Args:
        alphafold_protein: Path to AlphaFold protein PDB file
        experimental_protein: Path to experimental protein PDB file
        ligand: Path to ligand SDF file
        radius: Radius around ligand to analyze
        rmsd_value: Pre-calculated RMSD value to display
        output_file: Optional path to save HTML. If not provided and config is
            available, saves to config output directory with auto-generated
            filename
    """

    # Parse structures
    alphafold_coords = get_coordinates(parse_pdb(alphafold_protein))
    experimental_coords = get_coordinates(parse_pdb(experimental_protein))
    ligand_coords = get_coordinates(parse_sdf(ligand))

    # Calculate distances and find atoms within radius
    af_distances = distance_matrix(alphafold_coords, ligand_coords)
    ex_distances = distance_matrix(experimental_coords, ligand_coords)

    af_within_radius = np.any(af_distances <= radius, axis=1)
    ex_within_radius = np.any(ex_distances <= radius, axis=1)

    # Calculate ligand center for sphere
    ligand_center = ligand_coords.mean(axis=0)
    sphere_x, sphere_y, sphere_z = create_sphere_mesh(ligand_center, radius)

    fig = go.Figure()

    # Add all AlphaFold atoms (grayed out)
    fig.add_trace(
        go.Scatter3d(
            x=alphafold_coords[~af_within_radius, 0],
            y=alphafold_coords[~af_within_radius, 1],
            z=alphafold_coords[~af_within_radius, 2],
            mode="markers",
            marker={"size": 2, "color": "lightblue", "opacity": 0.3},
            name="AlphaFold (outside radius)",
            hovertemplate=(
                "AlphaFold (excluded)"
                "<br>X: %{x}<br>Y: %{y}<br>Z: %{z}<extra></extra>"
            ),
        )
    )

    # Add AlphaFold atoms within radius (highlighted)
    if np.any(af_within_radius):
        fig.add_trace(
            go.Scatter3d(
                x=alphafold_coords[af_within_radius, 0],
                y=alphafold_coords[af_within_radius, 1],
                z=alphafold_coords[af_within_radius, 2],
                mode="markers",
                marker={"size": 4, "color": "blue", "opacity": 0.8},
                name="AlphaFold (within radius)",
                hovertemplate=(
                    "AlphaFold (RMSD)"
                    "<br>X: %{x}<br>Y: %{y}<br>Z: %{z}<extra></extra>"
                ),
            )
        )

    # Add all experimental atoms (grayed out)
    fig.add_trace(
        go.Scatter3d(
            x=experimental_coords[~ex_within_radius, 0],
            y=experimental_coords[~ex_within_radius, 1],
            z=experimental_coords[~ex_within_radius, 2],
            mode="markers",
            marker={"size": 2, "color": "lightcoral", "opacity": 0.3},
            name="Experimental (outside radius)",
            hovertemplate=(
                "Experimental (excluded)"
                "<br>X: %{x}<br>Y: %{y}<br>Z: %{z}<extra></extra>"
            ),
        )
    )

    # Add experimental atoms within radius (highlighted)
    if np.any(ex_within_radius):
        fig.add_trace(
            go.Scatter3d(
                x=experimental_coords[ex_within_radius, 0],
                y=experimental_coords[ex_within_radius, 1],
                z=experimental_coords[ex_within_radius, 2],
                mode="markers",
                marker={"size": 4, "color": "red", "opacity": 0.8},
                name="Experimental (within radius)",
                hovertemplate=(
                    "Experimental (RMSD)"
                    "<br>X: %{x}<br>Y: %{y}<br>Z: %{z}<extra></extra>"
                ),
            )
        )

    # Add ligand points
    fig.add_trace(
        go.Scatter3d(
            x=ligand_coords[:, 0],
            y=ligand_coords[:, 1],
            z=ligand_coords[:, 2],
            mode="markers",
            marker={"size": 8, "color": "green", "opacity": 1.0},
            name="Ligand",
            hovertemplate=(
                "Ligand" "<br>X: %{x}<br>Y: %{y}<br>Z: %{z}<extra></extra>"
            ),
        )
    )

    # Add transparent sphere
    fig.add_trace(
        go.Surface(
            x=sphere_x,
            y=sphere_y,
            z=sphere_z,
            opacity=0.15,
            colorscale=[[0, "yellow"], [1, "yellow"]],
            showscale=False,
            name=f"Radius Sphere ({radius}Å)",
            hovertemplate=f"Radius: {radius}Å<extra></extra>",
        )
    )

    # Set title with RMSD value if provided
    title = f"Local RMSD Analysis (Radius: {radius}Å)"
    if rmsd_value is not None:
        title = f"Local RMSD Analysis (RMSD: {rmsd_value:.3f}Å, Radius: {radius}Å)"

    fig.update_layout(
        title=title,
        scene={
            "xaxis_title": "X (Å)",
            "yaxis_title": "Y (Å)",
            "zaxis_title": "Z (Å)",
            "aspectmode": "cube",
        },
        width=1000,
        height=800,
        showlegend=True,
    )

    # Determine output file path
    save_file = None
    if output_file:
        save_file = Path(output_file)

    # Save to file if path is determined
    if save_file:
        if save_file.suffix != ".html":
            save_file = save_file.with_suffix(".html")
        fig.write_html(str(save_file))
        print(f"Local RMSD visualization saved to: {save_file}")

    return fig
