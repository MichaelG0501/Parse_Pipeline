# Parse Sample-Level Pseudotime Trajectory Methodology

Scripts:

- `analysis/trajectory/parse_pseudotime_helpers.R`
- `analysis/trajectory/parse_pseudotime_sample_trajectory.R`

Status: active trajectory workflow.

## Purpose

This workflow builds a Monocle3 trajectory using sample labels rather than cell-state labels. The trajectory is rooted on T0 and used as the cached basis for linear trajectory reports and sample distance maps.

## Inputs

- `parse_outs/Auto_parse_merged.rds`

Configuration:

- Samples default to `T0,T1,T2,T4,R4,eR4`.
- Override with `AUTO_PARSE_PSEUDOTIME_SAMPLES`.
- Root sample defaults to `T0`.
- Override with `AUTO_PARSE_ROOT_SAMPLE`.

## Processing

1. Load the merged Parse Seurat object.
2. Detect the sample metadata column from `sample`, `orig.ident`, `Sample`, or `sample_id`.
3. Keep configured Parse samples only.
4. Drop samples below the minimum cell threshold.
5. Normalize, select variable features, scale data, run PCA, and run UMAP.
6. Convert the Seurat object to a Monocle3 cell data set.
7. Cluster cells, learn a principal graph, and order cells using T0 root cells.
8. Extract pseudotime and replace infinite values with `NA`.
9. Save the CDS, pseudotime vector, metadata, summary table, and combined trajectory PDF.

## Outputs

Output tiers:

- `intermediate`: `parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime_cds.rds`
- `intermediate`: `parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime.rds`
- `tables`: `parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime_metadata.csv`
- `tables`: `parse_outs/summary/Auto_parse_sample_pseudotime_summary.csv`
- `figures`: `parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime_combined.pdf`
- `logs`: `parse_outs/logs/run_summaries/parse_pseudotime_sample_trajectory_*.txt`

## Downstream Use

Required by:

- `parse_pseudotime_linear_sample_report.R`
- `parse_pseudotime_sample_distance_maps.R`
