#!/usr/bin/env python3
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


WD = Path("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline")
OUT = WD / "parse_outs" / "Auto_velocity"
TARGET_SAMPLES = ["T0", "T1", "T2", "T4", "R4", "eR4"]
LOOM_CONFIG = {
    "NACT1090_A": {"suffix": "__s2", "dir": OUT / "looms/A"},
    "NACT1090_B": {"suffix": "__s1", "dir": OUT / "looms/B"},
}


def extract_barcode(cell_id: str) -> str:
    bc = str(cell_id).split(":")[-1]
    if bc.endswith("-1"):
        bc = bc[:-2]
    if bc.endswith("x"):
        bc = bc[:-1]
    return bc


def load_one(sublibrary: str, cfg: dict) -> ad.AnnData:
    looms = sorted(cfg["dir"].glob("*.loom"))
    if len(looms) != 1:
        raise FileNotFoundError(f"Expected exactly one loom in {cfg['dir']}, found {len(looms)}")
    adata = sc.read_loom(str(looms[0]), sparse=True)
    adata.var_names_make_unique()
    barcodes = [extract_barcode(x) for x in adata.obs_names]
    adata.obs_names = [f"{bc}{cfg['suffix']}" for bc in barcodes]
    adata.obs["bc_wells"] = barcodes
    adata.obs["sublibrary"] = sublibrary
    adata.obs_names_make_unique()
    return adata


def save_current(name: str, dpi: int = 220) -> None:
    path = OUT / "figures" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()


def sample_palette() -> dict[str, str]:
    return {
        "T0": "#2f6fbb",
        "T1": "#43a047",
        "T2": "#f2a541",
        "T4": "#d64b4b",
        "R4": "#6f4aa8",
        "eR4": "#00a6a6",
    }


def sample_centroids_and_velocities(adata: ad.AnnData) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
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


def add_centroid_label(ax: plt.Axes, center: np.ndarray, label: str, offset: tuple[float, float]) -> None:
    text = ax.text(
        center[0] + offset[0],
        center[1] + offset[1],
        label,
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
        color="black",
        zorder=25,
    )
    text.set_path_effects([pe.withStroke(linewidth=4, foreground="white")])


def plot_sample_velocity_directions(adata: ad.AnnData) -> None:
    palette = sample_palette()
    xy = adata.obsm["X_umap"]
    samples = adata.obs["sample"].astype(str).to_numpy()
    centers, velocities = sample_centroids_and_velocities(adata)
    fig, ax = plt.subplots(figsize=(9.2, 7.4))
    for sample in TARGET_SAMPLES:
        mask = samples == sample
        if not np.any(mask):
            continue
        ax.scatter(
            xy[mask, 0],
            xy[mask, 1],
            s=7,
            alpha=0.28,
            linewidths=0,
            color=palette[sample],
            label=sample,
        )

    spans = np.ptp(xy, axis=0)
    arrow_scale = 0.13 * float(np.linalg.norm(spans))
    label_offsets = {
        "T0": (-0.35, 0.34),
        "T1": (0.38, 0.34),
        "T2": (-0.45, 0.28),
        "T4": (-0.35, 0.28),
        "R4": (-0.42, 0.28),
        "eR4": (-0.42, -0.32),
    }
    for sample in TARGET_SAMPLES:
        if sample not in centers:
            continue
        center = centers[sample]
        direction = velocities[sample]
        norm = float(np.linalg.norm(direction))
        if norm == 0:
            continue
        direction = direction / norm * arrow_scale
        ax.annotate(
            "",
            xy=center + direction,
            xytext=center,
            arrowprops=dict(arrowstyle="-|>", color="black", lw=3.0, mutation_scale=20),
            zorder=10,
        )
        ax.scatter(center[0], center[1], s=110, color=palette[sample], edgecolor="black", linewidth=1.1, zorder=11)
        add_centroid_label(ax, center, sample, label_offsets.get(sample, (0.35, 0.35)))

    ax.set_title("RNA velocity direction: 1090 Parse", fontsize=18, pad=14)
    ax.set_axis_off()
    legend = ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        title="Sample",
        fontsize=15,
        title_fontsize=17,
        markerscale=2.4,
        borderaxespad=0.7,
    )
    for handle in legend.legend_handles:
        handle.set_alpha(1.0)
    fig.tight_layout()
    save_current("Auto_velocity_direction_by_sample.png", dpi=260)


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
    edge_table["plotted"] = edge_table["velocity_alignment"] >= 0.70
    edge_table.to_csv(OUT / "tables" / "Auto_velocity_sample_direction_edges.csv", index=False)

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
        rad = 0.18 if index[source] < index[target] else -0.18
        arrow = FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=26,
            linewidth=1.4 + 5.2 * score,
            color="black",
            alpha=0.22 + 0.62 * score,
            shrinkA=34,
            shrinkB=34,
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
            s=1850,
            color=palette[sample],
            edgecolor="black",
            linewidth=2.2,
            zorder=10,
        )
        text = ax.text(
            center[0],
            center[1],
            sample,
            ha="center",
            va="center",
            fontsize=18,
            fontweight="bold",
            color="black",
            zorder=11,
        )
        text.set_path_effects([pe.withStroke(linewidth=3.8, foreground="white")])

    ax.set_title("Sample-level velocity network: 1090 Parse", fontsize=19, pad=14)
    ax.text(
        0.5,
        -0.04,
        "Arrows show sample-to-sample centroid directions aligned with mean velocity (cosine >= 0.70).",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=12,
    )
    ax.set_axis_off()
    fig.tight_layout()
    save_current("Auto_velocity_sample_network.png", dpi=260)


def plot_r4_velocity_focus(adata: ad.AnnData) -> None:
    palette = sample_palette()
    scv.pl.velocity_embedding_grid(
        adata,
        basis="umap",
        color="sample",
        groups=["R4"],
        arrow_color="black",
        arrow_size=4.8,
        arrow_length=7.5,
        density=0.82,
        scale=0.3,
        autoscale=False,
        alpha=0.38,
        size=14,
        legend_loc="right margin",
        title="R4 velocity direction: 1090 Parse",
        frameon=False,
        show=False,
    )
    ax = plt.gca()
    xy = adata.obsm["X_umap"]
    samples = adata.obs["sample"].astype(str).to_numpy()
    for sample in TARGET_SAMPLES:
        mask = samples == sample
        if not np.any(mask):
            continue
        center = xy[mask].mean(axis=0)
        ax.scatter(center[0], center[1], s=105, color=palette[sample], edgecolor="black", linewidth=1.0, zorder=20)
        ax.text(
            center[0],
            center[1],
            sample,
            ha="center",
            va="center",
            fontsize=9.5,
            fontweight="bold",
            color="white",
            zorder=21,
        )
    save_current("Auto_velocity_R4_direction_focus.png", dpi=240)


def main() -> None:
    meta = pd.read_csv(OUT / "tables/Auto_velocity_nonPDO_cell_metadata.csv")
    meta = meta[meta["sample"].isin(TARGET_SAMPLES)].copy()
    meta = meta.set_index("combined_cell_id")

    adatas = [load_one(sublib, cfg) for sublib, cfg in LOOM_CONFIG.items()]
    adata = ad.concat(adatas, join="outer", merge="same", label="sublibrary_source", index_unique=None)
    adata = adata[adata.obs_names.isin(meta.index)].copy()
    meta = meta.loc[adata.obs_names]
    for col in meta.columns:
        adata.obs[col] = meta[col].astype(str).values
    adata.obs["sample"] = pd.Categorical(adata.obs["sample"], categories=TARGET_SAMPLES, ordered=True)

    scv.settings.verbosity = 3
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    if adata.n_vars > 3000:
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True, flavor="seurat")
    sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
    sc.tl.umap(adata, min_dist=0.35, spread=1.0)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata)

    sc.pl.umap(adata, color=["sample", "sublibrary"], frameon=False, wspace=0.35, show=False)
    save_current("Auto_umap_by_sample_sublibrary.png")

    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="sample",
        legend_loc="right margin",
        title="RNA velocity stream: 1090 Parse",
        frameon=False,
        show=False,
    )
    save_current("Auto_velocity_stream_by_sample.png")

    scv.pl.velocity_embedding_grid(
        adata,
        basis="umap",
        color="sample",
        arrow_color="black",
        arrow_size=4.4,
        arrow_length=7.0,
        density=0.72,
        scale=0.32,
        autoscale=False,
        alpha=0.48,
        size=16,
        legend_loc="right margin",
        title="RNA velocity grid: 1090 Parse",
        frameon=False,
        show=False,
    )
    save_current("Auto_velocity_grid_by_sample.png")

    scv.tl.velocity_confidence(adata)
    plot_sample_velocity_directions(adata)
    plot_sample_velocity_network(adata)
    plot_r4_velocity_focus(adata)
    sc.pl.umap(adata, color=["velocity_length", "velocity_confidence"], frameon=False, show=False)
    save_current("Auto_umap_velocity_quality.png")

    OUT.joinpath("tables").mkdir(parents=True, exist_ok=True)
    adata.obs.to_csv(OUT / "tables/Auto_scvelo_cell_metadata.csv")
    adata.write(OUT / "Auto_scvelo_nonPDO.h5ad", compression="gzip")


if __name__ == "__main__":
    main()
