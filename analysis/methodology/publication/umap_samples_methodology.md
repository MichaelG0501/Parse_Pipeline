# UMAP Samples Methodology

## Overview
This script creates a publication-quality UMAP visualization of the 6 Parse treatment-response timepoints (T0, T1, T2, T4, R4, eR4) along with the PDO sample. 

## Inputs
- `parse_outs/Auto_parse_merged.rds`: The final merged Seurat object containing UMAP coordinates (`@reductions$umap`) and cell metadata.

## Methodology
1. **Subsetting**: 
   The pipeline filters the merged Seurat object for cells belonging only to the specified samples: `PDO`, `T0`, `T1`, `T2`, `T4`, `R4`, and `eR4`.
2. **Visual Hierarchy (Leveled Legend)**:
   The user explicitly requested a "leveled" legend where "PDO" and "2.5D Cell Lines" act as groups. We implement this cleanly by configuring `ggplot2` to map `sample` as a factor including a transparent dummy level `2.5D Cell Lines`. 
3. **Randomization**:
   To prevent "overplotting bias" (where cells from the last-plotted sample artificially dominate the visual density), we randomly shuffle the row order of `plot_data` before passing it to `ggplot2`.
4. **Style**:
   Uses the project's consistent `theme_nature` variant for UMAPs (axis text and ticks removed, solid line axes, Arial font, proper legends). Colours are dynamically pulled from `parse_sample_colours`.

## Outputs
- `parse_outs/publication/umap_samples/figures/umap_samples.pdf`
- `parse_outs/publication/umap_samples/figures/umap_samples.png`

These are terminal figures intended for presentation or manuscripts. They do not feed any active downstream workflows.
