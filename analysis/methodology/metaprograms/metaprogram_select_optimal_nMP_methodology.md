# Parse Metaprogram Optimal nMP Selection Methodology

Script: `analysis/metaprograms/parse_metaprogram_select_optimal_nMP.R`

Status: active metaprogram selection workflow.

## Purpose

This script selects the default number of Parse metaprograms from the GeneNMF sweep using silhouette and WSS diagnostics.

## Inputs

`parse_outs/Auto_parse_metaprograms/Metaprogrammes_Results/Auto_parse_geneNMF_metaprograms_nMP_<k>.rds` for k values 4 to 35.

## Processing

1. Iterate over k = 4:35.
2. Load each available GeneNMF metaprogram object.
3. Convert program similarity to distance with `1 - programs.similarity`.
4. Cut the GeneNMF program tree into k clusters.
5. Compute average silhouette width.
6. Compute within-cluster sum of squares from the distance matrix.
7. Detect a WSS knee.
8. Select the optimal nMP from the diagnostics.
9. Save the selected nMP text file, metrics table, metrics figure, and default selected metaprogram object.

## Outputs

Output tiers:

- `tables`: `Auto_parse_optimal_nMP_metrics.csv`
- `figures`: `Auto_parse_optimal_nMP_metrics.png`
- `intermediate`: `Auto_parse_optimal_nMP.txt`
- `intermediate`: `Auto_parse_MP_outs_default.rds`
- `figures`: `Auto_parse_metaprograms_heatmap.png`
- `logs`: `parse_outs/logs/run_summaries/parse_metaprogram_select_optimal_nMP_*.txt`

## Downstream Use

Required by the second GeneNMF/UCell run, enrichment annotation, MP correlation, external comparisons, and cell-state scoring.
