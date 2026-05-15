# Parse Metaprogram GeneNMF Discovery Methodology

Script: `analysis/metaprograms/parse_metaprogram_geneNMF_discovery.R`

Status: active upstream metaprogram workflow.

## Purpose

This script discovers Parse-derived malignant epithelial metaprograms from the six Parse treatment-response samples and scores selected metaprograms using UCell after optimal nMP selection.

## Inputs

- `parse_outs/by_samples/T0/Auto_T0_final.rds`
- `parse_outs/by_samples/T1/Auto_T1_final.rds`
- `parse_outs/by_samples/T2/Auto_T2_final.rds`
- `parse_outs/by_samples/T4/Auto_T4_final.rds`
- `parse_outs/by_samples/R4/Auto_R4_final.rds`
- `parse_outs/by_samples/eR4/Auto_eR4_final.rds`

Optional caches:

- `parse_outs/Auto_parse_metaprograms/Auto_parse_list_geneNMF.rds`
- `parse_outs/Auto_parse_metaprograms/Auto_parse_merged_geneNMF.rds`
- `parse_outs/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds`

## Processing

1. Load the six Parse final Seurat objects.
2. Exclude `PDO`, `SUR1090_Treated`, and `SUR1090_Untreated` from GeneNMF discovery.
3. Cache the per-sample list and merged Seurat object.
4. Normalize the merged object for UCell scoring.
5. Run GeneNMF across the configured nMP range if outputs do not already exist.
6. Save each nMP metaprogram object under `Metaprogrammes_Results/`.
7. If `Auto_parse_optimal_nMP.txt` exists, load the selected nMP object.
8. Filter bad MPs with negative silhouette and low sample coverage.
9. Score retained MPs with UCell and cache the UCell score matrix.
10. Save the final GeneNMF Seurat object and diagnostic violin plot.

## Caching

This script is designed to be rerun after `parse_metaprogram_select_optimal_nMP.R`:

- First run: generates the GeneNMF sweep.
- Optimal-nMP script selects the default nMP.
- Second run: scores selected MPs with UCell.

Existing heavy objects are reused where possible.

## Outputs

Output tiers:

- `intermediate`: `Auto_parse_list_geneNMF.rds`, `Auto_parse_merged_geneNMF.rds`, `Auto_parse_geneNMF_outs.rds`, `Auto_parse_final_geneNMF.rds`
- `intermediate`: `Metaprogrammes_Results/Auto_parse_geneNMF_metaprograms_nMP_<k>.rds`
- `intermediate`: `Auto_parse_UCell_scores_filtered_nMP<k>.rds`
- `figures`: `Auto_parse_vln_origident.png`
- `logs`: `parse_outs/logs/run_summaries/parse_metaprogram_geneNMF_discovery_*.txt`

## Downstream Use

Required by optimal nMP selection, enrichment annotation, MP correlation, external reference comparison, and Approach B/noreg state assignment.
