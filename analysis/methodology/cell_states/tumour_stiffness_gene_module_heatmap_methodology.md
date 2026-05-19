# Tumour Stiffness Gene Module Heatmap Methodology

Script: `analysis/cell_states/parse_tumour_stiffness_gene_module_heatmap.R`

Status: active terminal figure-generation workflow.

## Purpose

This script generates a sample-resolved heatmap for a curated tumour-stiffness gene module. It reports the all-gene stiffness module, subgroup module scores, and individual gene expression in one grouped matrix.

## Inputs

- `parse_outs/Auto_parse_merged.rds`

## Processing

1. Load the merged Parse/PDO/SUR1090 Seurat object.
2. Use the `RNA` assay log-normalized data layer. If that layer is empty, normalize the RNA assay before scoring.
3. Match the requested gene symbols case-insensitively against the object row names.
4. Average log-normalized expression per sample for each detected gene.
5. Compute module activity as the mean log-normalized expression across detected member genes:
   - all stiffness genes together
   - core mechanosensing / membrane tension
   - integrin-focal adhesion-actin linkage
   - Rho GTPase cytoskeletal regulators
   - cell polarity / cortical organisation
   - actin cytoskeleton and contractility
6. Row-z-score each gene/module across samples for the heatmap colour scale using a vibrant, publication-grade blue-white-red palette.
7. Cluster genes within each biological group, while keeping module rows at the top of their groups and keeping sample columns in canonical Parse order.
8. Export a ComplexHeatmap PDF with horizontal row group titles (word-wrapped), top-aligned column labels, cell-count bars directly adjacent to the columns, increased text sizes, and raw score labels inside each cell (8pt size) to maximize clarity.

## Outputs

Output root:

`parse_outs/cell_states/tumour_stiffness_gene_module/`

Output tiers:

- `intermediate`: `Auto_parse_tumour_stiffness_matrices.rds`
- `tables`: gene availability and sample-level score CSV files
- `figures`: heatmap PDF and PNG
- `reports`: heatmap notes with the figure contract and missing-gene summary
- `logs`: `parse_outs/logs/run_summaries/parse_tumour_stiffness_gene_module_heatmap_*.txt`

## Downstream Use

This is a terminal presentation/manuscript-facing output. No active downstream script depends on it.
