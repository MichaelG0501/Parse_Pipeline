# InferCNA Carroll Reference Methodology

Script: `analysis/cnv/parse_infercna_carroll_reference.R`

Status: active CNV workflow.

## Purpose

This script estimates copy-number alteration signal for Parse/PDO/SUR1090 epithelial cells using the Carroll_2023 scATLAS reference and produces slide-scale CNA summaries.

## Inputs

Per-sample Seurat objects:

- `parse_outs/by_samples/<sample>/Auto_<sample>_final.rds`

External references:

- Carroll_2023 reference: `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Carroll_2023_reference.rds`
- Gene order: `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt`

## Processing

1. Load the Carroll_2023 reference matrix.
2. Load all nine target sample objects.
3. Find common genes across reference and targets.
4. Convert target counts to CPM.
5. Combine targets and reference.
6. Run InferCNA using the external reference.
7. Save full and target-only CNA matrices.
8. Save combined metadata and reference summary.
9. Generate binned CNA heatmaps:
   - all nine samples
   - PDO versus T0 comparison
10. Save per-sample mean CNA profiles.
11. Generate CNA signal/correlation scatter plots and summary table.

## Caching

The script saves InferCNA matrices and heatmap input RDS files. Plot styling for heatmaps can be revised from saved heatmap input files without rerunning InferCNA.

## Outputs

Output root:

`parse_outs/cnv/`

Output tiers:

- `intermediate`: `Auto_parse_infercna_outs_Carroll_2023.rds`, `Auto_parse_infercna_target_outs_Carroll_2023.rds`, heatmap input RDS files
- `tables`: metadata CSV, reference summary TSV, scatter summary CSV
- `figures`: CNA heatmap PDFs and scatter PDF
- `logs`: `parse_outs/logs/run_summaries/parse_infercna_carroll_reference_*.txt`

## Downstream Use

Terminal CNV figure workflow. No active downstream script consumes its outputs.
