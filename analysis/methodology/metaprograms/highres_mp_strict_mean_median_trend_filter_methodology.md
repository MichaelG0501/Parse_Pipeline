# High-Resolution Strict Mean/Median MP Trend Filter Methodology

Script: `analysis/metaprograms/parse_highres_mp_strict_mean_median_trend_filter.R`

Status: active high-resolution metaprogram workflow.

## Purpose

This script derives high-resolution Parse metaprograms from the GeneNMF program sweep and keeps only MPs with consistent mean and median activity trends across the treatment-response time course.

## Inputs

- `parse_outs/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds`
- `parse_outs/by_samples/<sample>/Auto_<sample>_final.rds`
- 3CA metaprograms: `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`
- Cell-cycle genes: `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv`
- Developmental enrichment files under `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/`

## Processing

The script supports `--mode=score`, `--mode=enrich`, `--mode=excel`, and `--mode=all`.

### Score Mode

1. Load the GeneNMF program object.
2. Set high-resolution nMP to half the total NMF programs.
3. Generate or load high-resolution metaprogram object.
4. Write MP gene and program-membership tables.
5. Load final Parse sample objects and find common genes.
6. Score high-resolution MPs with UCell or load cached UCell scores.
7. Summarize mean and median UCell activity by sample.
8. Classify trends across `T0`, `T1`, `T2`, `T4`, `R4`, and `eR4`.
9. Retain MPs where mean and median directions agree.
10. Classify recovery behavior into increase/decrease and consistent/reverted groups.
11. Generate selected-MP boxplots, mean/median trends, heatmaps, and selected gene tables.

### Enrich Mode

1. Load selected MP genes and trend summary.
2. Build 3CA non-cell-cycle labels.
3. Run GO, Hallmark, 3CA, and developmental enrichment.
4. Write enrichment RDS/PDF/PNG outputs with readable slide-scale fonts.

### Excel Mode

1. Load selected MP genes and trend summary.
2. Write an XLSX workbook of retained MP genes grouped by trend label.

## Caching

The following expensive intermediates are reused:

- High-resolution metaprogram RDS.
- UCell score matrix.
- Cell metadata RDS.
- Selected MP genes RDS.
- Trend summary CSV.
- Enrichment RDS.

This allows plot styling and figure sizing to be adjusted without repeating GeneNMF or UCell scoring.

## Outputs

Output root:

`parse_outs/Auto_parse_highres_metaprogram_trends/`

Output tiers:

- `intermediate`: `Auto_parse_highres_geneNMF_metaprograms_nMP117.rds`, `Auto_parse_highres_UCell_scores_nMP117.rds`, `Auto_parse_highres_cell_metadata_nMP117.rds`, `Auto_parse_highres_selected_mp_genes_nMP117.rds`
- `tables`: trend summary, retained summary, scored-gene counts, MP gene CSVs, 3CA label tables
- `figures`: boxplot PDFs, mean/median trend PDFs, heatmap PDF/PNG, enrichment PNGs
- `reports`: enrichment annotation PDF, MP gene summary XLSX
- `logs`: `parse_outs/logs/run_summaries/parse_highres_mp_strict_mean_median_trend_filter_*.txt`

## Downstream Use

Terminal high-resolution MP workflow. No active downstream script should treat the selected high-resolution MPs as canonical state labels unless explicitly requested.
