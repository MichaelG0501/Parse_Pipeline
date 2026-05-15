# Parse Internal Metaprogram Correlation Methodology

Script: `analysis/metaprograms/parse_metaprogram_internal_correlation.R`

Status: active terminal diagnostics workflow.

## Purpose

This script quantifies relationships among selected Parse metaprograms using gene-set overlap and UCell activity correlation.

## Inputs

- `parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds`
- `parse_outs/Auto_parse_metaprograms/Auto_parse_final_geneNMF.rds`
- Fallback: `parse_outs/Auto_parse_metaprograms/Auto_parse_merged_geneNMF.rds`
- `parse_outs/Auto_parse_metaprograms/Auto_parse_UCell_scores_filtered_nMP<k>.rds`

## Processing

1. Load selected Parse metaprogram definitions.
2. Filter negative-silhouette and low-coverage MPs.
3. Determine canonical MP order from the GeneNMF tree.
4. Load or compute UCell scores.
5. Compute global MP-MP Spearman correlation.
6. Compute per-sample correlation summaries.
7. Compute MP gene-set Jaccard overlap.
8. Draw a combined PDF with global correlation, per-sample correlation, Jaccard overlap, and per-cell UCell heatmap.
9. Save a compact summary table.

## Outputs

Output tiers:

- `figures`: `Auto_parse_nMP<k>_analysis_combined.pdf`
- `tables`: `Auto_parse_mp_correlation_summary.csv`
- `logs`: `parse_outs/logs/run_summaries/parse_metaprogram_internal_correlation_*.txt`

## Downstream Use

Terminal diagnostic workflow. No active downstream script consumes the combined PDF.
