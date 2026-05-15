# Filtering Summary Plots Methodology

Script: `analysis/summary/parse_filtering_summary_plots.R`

Status: active QC-summary plotting script.

## Purpose

This script summarizes global and per-sample cell retention across the Parse QC filters.

## Inputs

- `parse_outs/summary/Auto_parse_filtering_summary.csv`
- `parse_outs/summary/Auto_parse_thresholds.csv`
- `parse_outs/Auto_parse_merged.rds`

The merged Seurat object supplies final post-QC cell counts by `orig.ident`.

## Processing

1. Load QC filtering summary and final merged Seurat object.
2. Count final cells per sample from `merged_obj@meta.data$orig.ident`.
3. Join final counts onto the QC summary table.
4. Read the threshold table so axis labels contain the actual QC cutoffs.
5. Build three plots:
   - Global retained-cell bar plot across filters.
   - Per-sample before/after bar plot.
   - Per-sample filtering trajectory line plot.

## Outputs

Output tiers:

- `figures`: `parse_outs/plots/Filtering_and_Summary/global_filter_bar.png`
- `figures`: `parse_outs/plots/Filtering_and_Summary/per_sample_before_after.png`
- `figures`: `parse_outs/plots/Filtering_and_Summary/per_sample_filtering_steps.png`
- `logs`: `parse_outs/logs/run_summaries/parse_filtering_summary_plots_*.txt`

## Downstream Use

Terminal diagnostic figures only. No active downstream script consumes these plots.
