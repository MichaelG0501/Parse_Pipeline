#!/usr/bin/env python3
"""Publication-quality RNA velocity visualizations.

Inputs:
    parse_outs/Auto_velocity/Auto_scvelo_nonPDO.h5ad

Outputs:
    parse_outs/Auto_velocity/figures/Auto_velocity_stream_by_sample.png
    parse_outs/Auto_velocity/figures/Auto_velocity_sample_network.png
"""

import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

COMMON_DIR = Path(__file__).resolve().parents[1] / "common"
sys.path.insert(0, str(COMMON_DIR))
from parse_pipeline_config import PARSE_SAMPLES, PROJECT_ROOT, VELOCITY_OUT

OUT = VELOCITY_OUT
TARGET_SAMPLES = PARSE_SAMPLES

def save_current(name: str, dpi: int = 300) -> None:
    path = OUT / "figures" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()

def sample_palette():
    return {
        "T0": "#0072B2",
        "T1": "#E69F00",
        "T2": "#009E73",
        "T4": "#D55E00",
        "R4": "#CC79A7",
        "eR4": "#56B4E9",
    }

def sample_centroids_and_velocities(adata: ad.AnnData):
    xy = adata.obsm["X_umap"]
    vv = adata.obsm["velocity_umap"]
    samples = adata.obs["sample"].astype(str).to_numpy()
    centers = {}
    velocities = {}
    for sample in TARGET_SAMPLES:
        mask = samples == sample
        if not np.any(mask):
            continue
        centers[sample] = xy[mask].mean(axis=0)
        velocities[sample] = vv[mask].mean(axis=0)
    return centers, velocities

def plot_sample_velocity_network(adata: ad.AnnData) -> None:
    palette = sample_palette()
    centers, velocities = sample_centroids_and_velocities(adata)
    edge_rows = []
    for source in TARGET_SAMPLES:
        if source not in centers:
            continue
        velocity = velocities[source]
        velocity_norm = float(np.linalg.norm(velocity))
        if velocity_norm == 0:
            continue
        for target in TARGET_SAMPLES:
            if target == source or target not in centers:
                continue
            target_direction = centers[target] - centers[source]
            target_norm = float(np.linalg.norm(target_direction))
            if target_norm == 0:
                continue
            score = float(np.dot(velocity, target_direction) / (velocity_norm * target_norm))
            edge_rows.append({"source": source, "target": target, "velocity_alignment": score})

    edge_table = pd.DataFrame(edge_rows).sort_values(["source", "velocity_alignment"], ascending=[True, False])
    
    allowed_edges = [("R4", "T0"), ("T4", "R4"), ("T4", "eR4"), ("T2", "T4")]
    
    def is_allowed(src, tgt):
        return (src, tgt) in allowed_edges
    
    edge_table["is_allowed"] = edge_table.apply(lambda row: is_allowed(row["source"], row["target"]), axis=1)
    edge_table["plotted"] = edge_table["is_allowed"]
    
    edge_table.to_csv(OUT / "tables" / "Auto_velocity_sample_direction_edges_pub.csv", index=False)

    fig, ax = plt.subplots(figsize=(9.8, 7.2))
    node_xy = np.vstack([centers[sample] for sample in TARGET_SAMPLES if sample in centers])
    margin = np.ptp(node_xy, axis=0).max() * 0.24
    ax.set_xlim(node_xy[:, 0].min() - margin, node_xy[:, 0].max() + margin)
    ax.set_ylim(node_xy[:, 1].min() - margin, node_xy[:, 1].max() + margin)

    index = {sample: i for i, sample in enumerate(TARGET_SAMPLES)}
    plotted_edges = edge_table[edge_table["plotted"]].copy()
    
    for _, edge in plotted_edges.iterrows():
        source = edge["source"]
        target = edge["target"]
        score = float(edge["velocity_alignment"])
        start = centers[source]
        end = centers[target]
        rad = 0.15 if index[source] < index[target] else -0.15
        
        # Alpha and linewidth scale with score
        alpha_val = 0.6 + 0.4 * score
        lw_val = 2.0 + 3.0 * score
        color_val = "black"
        
        arrow = FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=22,
            linewidth=lw_val,
            color=color_val,
            alpha=alpha_val,
            shrinkA=30,
            shrinkB=30,
            connectionstyle=f"arc3,rad={rad}",
            zorder=3,
        )
        ax.add_patch(arrow)

    for sample in TARGET_SAMPLES:
        if sample not in centers:
            continue
        center = centers[sample]
        ax.scatter(
            center[0],
            center[1],
            s=1200,
            color=palette[sample],
            edgecolor="black",
            linewidth=2.0,
            zorder=10,
        )
        text = ax.text(
            center[0],
            center[1],
            sample,
            ha="center",
            va="center",
            fontsize=16,
            fontweight="bold",
            color="black",
            zorder=11,
        )
        text.set_path_effects([pe.withStroke(linewidth=3.0, foreground="white")])

    ax.set_title("Sample-level velocity network: 1090 Parse", fontsize=18, pad=14)
    ax.text(
        0.5,
        -0.04,
        "Arrows show manually selected sample-to-sample centroid directions aligned with mean velocity.\nThickness reflects alignment score.",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=12,
    )
    ax.set_axis_off()
    fig.tight_layout()
    save_current("Auto_velocity_sample_network.png", dpi=300)

def main() -> None:
    adata = sc.read_h5ad(OUT / "Auto_scvelo_nonPDO.h5ad")
    scv.settings.verbosity = 3

    palette = sample_palette()
    if "sample" in adata.obs:
        # Assign colors directly to the AnnData object in scanpy's expected location
        categories = adata.obs["sample"].cat.categories
        adata.uns["sample_colors"] = [palette[c] for c in categories]

    # Plot 1: Stream plot
    # Updated alpha to 0.9 to match parse_umap_samples.R
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="sample",
        alpha=0.9,
        size=40,
        legend_loc="right margin",
        title="RNA velocity stream: 1090 Parse",
        frameon=False,
        show=False,
    )
    save_current("Auto_velocity_stream_by_sample.png", dpi=300)

    # Plot 2: Sample network plot
    plot_sample_velocity_network(adata)

if __name__ == "__main__":
    main()
