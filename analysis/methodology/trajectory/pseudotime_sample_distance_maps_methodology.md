# Parse Pseudotime Sample Distance Maps Methodology

Script: `analysis/trajectory/parse_pseudotime_sample_distance_maps.R`

Status: active terminal trajectory distance workflow.

## Purpose

This script computes multiple sample-to-sample distances from the cached Monocle3 trajectory and visualizes T0-to-sample relationships.

## Inputs

Cached trajectory assets:

- `parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime_cds.rds`
- `parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime.rds`
- `parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime_metadata.csv`

## Processing

1. Load cached trajectory assets.
2. Extract graph, UMAP, sample, and pseudotime data.
3. Compute distance matrices by:
   - UMAP centroid Euclidean distance.
   - Principal graph geodesic distance between sample centroids.
   - Principal graph geodesic distance between sample medoids.
   - Directed pseudotime mean distance.
   - Directed pseudotime median distance.
4. Convert matrices to long format.
5. Write method-specific CSV matrices.
6. Write T0-to-sample distance table to both pseudotime and summary locations.
7. Generate a method-comparison heatmap.
8. Generate an interconnected network-style T0 distance node plot.

## Outputs

Output root:

`parse_outs/pseudotime_samples/sample_distance_pseudotime/`

Output tiers:

- `intermediate`: `Auto_parse_sample_distance_matrices.rds`
- `tables`: method-specific matrix CSVs, `Auto_parse_sample_distance_long.csv`, `Auto_parse_sample_distance_from_T0.csv`
- `figures`: `Auto_parse_sample_distance_method_comparison_heatmap.pdf`
- `figures`: `Auto_parse_sample_distance_from_T0_nodeplot.pdf`
- `logs`: `parse_outs/logs/run_summaries/parse_pseudotime_sample_distance_maps_*.txt`

## Downstream Use

Terminal figure workflow. No active downstream script consumes these outputs.
