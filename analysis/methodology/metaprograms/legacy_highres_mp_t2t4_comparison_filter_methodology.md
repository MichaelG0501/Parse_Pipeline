# Legacy T2/T4-High High-Resolution MP Comparison Methodology

Script: `analysis/metaprograms/legacy_parse_highres_mp_t2t4_comparison_filter.R`

Status: legacy comparison only.

## Purpose

This script retains high-resolution MPs where mean UCell activity in both T2 and T4 is higher than both T0 and eR4. It was kept as an alternative comparison to the active strict mean/median trend workflow.

## Inputs

Same input class as the active high-resolution trend workflow:

- `parse_outs/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds`
- `parse_outs/by_samples/<sample>/Auto_<sample>_final.rds`
- 3CA, cell-cycle, and developmental enrichment reference files.

## Processing

1. Build or load high-resolution nMP117 metaprograms.
2. Score MPs with UCell or reuse cached scores.
3. Summarize mean activity by sample.
4. Retain MPs where T2 and T4 means exceed both T0 and eR4.
5. Produce comparison figures, enrichment outputs, and MP gene tables in an isolated legacy folder.

## Outputs

Output root:

`parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/`

Output tiers:

- `intermediate`: T2/T4 high-resolution metaprogram, UCell, cell metadata, selected genes RDS files
- `tables`: trend/filter summaries, 3CA labels, selected gene CSV files
- `figures`: activity boxplots, trend plots, heatmaps, enrichment PNG/PDF files
- `reports`: MP gene summary XLSX
- `logs`: `parse_outs/logs/run_summaries/legacy_parse_highres_mp_t2t4_comparison_filter_*.txt`

## Downstream Use

No active downstream script should consume these outputs. Use only for comparison plots or sensitivity checks.
