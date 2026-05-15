# QC Marker Heatmap Methodology

Script: `analysis/plotting/parse_qc_marker_heatmap.R`

Status: active diagnostic plotting script.

## Purpose

This script visualizes QC-stage marker expression before and after filtering. It is intended for inspection of broad lineage markers and housekeeping behavior, not final cell-state assignment.

## Inputs

- `parse_outs/by_samples/<sample>/Auto_<sample>_raw.rds`
- `parse_outs/by_samples/<sample>/Auto_<sample>_final.rds`
- QC sample manifest and summary files under `parse_outs/summary/`

## Processing

1. Define the fixed sample order and QC thresholds.
2. Define housekeeping genes and broad lineage marker panels.
3. For raw and final stages:
   - Load per-sample Seurat objects.
   - Extract count matrices.
   - Compute housekeeping expression.
   - Compute marker expression summaries.
   - Build sample-stage metadata.
4. Draw two large ComplexHeatmap-style marker heatmaps:
   - Raw/no filtering.
   - Final filtered singlets.
5. Write a stage summary table with plotted cell counts.

## Outputs

Output tiers:

- `tables`: `parse_outs/summary/Auto_parse_qc_heatmap_stage_summary.csv`
- `figures`: `parse_outs/plots/Visualisation_and_Heatmaps/QC_prefilter.png`
- `figures`: `parse_outs/plots/Visualisation_and_Heatmaps/QC_final.png`
- `reports`: `parse_outs/plots/Visualisation_and_Heatmaps/QC_heatmaps.pdf`
- `logs`: `parse_outs/logs/run_summaries/parse_qc_marker_heatmap_*.txt`

## Downstream Use

Terminal QC diagnostic only. The marker panels here are not the preferred cell-state definitions.
