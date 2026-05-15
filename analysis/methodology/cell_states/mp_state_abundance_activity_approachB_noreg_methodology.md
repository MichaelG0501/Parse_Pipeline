# MP State Abundance and Activity Methodology: Approach B, noreg

Script: `analysis/cell_states/parse_mp_state_abundance_activity_approachB_noreg.R`

Status: active preferred cell-state workflow.

## Purpose

This script assigns Parse/PDO/SUR1090 cells to metaprogram-derived epithelial state groups and generates abundance/activity plots. The preferred current state definition is PDO-pipeline state assignment using Approach B with noreg-normalized scores.

## Inputs

Per-sample Seurat objects:

- `parse_outs/by_samples/<sample>/Auto_<sample>_final.rds`

Parse metaprogram inputs:

- `parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds`
- `parse_outs/Auto_parse_metaprograms/Auto_parse_UCell_scores_filtered_nMP<k>.rds`

External reference metaprograms:

- scATLAS: `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
- PDO pipeline: `/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`

## Processing

1. Load all nine samples: `T0`, `T1`, `T2`, `T4`, `R4`, `eR4`, `PDO`, `SUR1090_Untreated`, and `SUR1090_Treated`.
2. Intersect genes across all samples and references.
3. Build combined counts and cell metadata.
4. Score scATLAS metaprograms with UCell or load cached scores.
5. Z-normalize scATLAS UCell scores by sample and study without regression.
6. Assign scATLAS top MP and state-group abundance for comparison.
7. Score PDO-pipeline metaprograms with UCell or load cached scores.
8. Z-normalize PDO-pipeline scores without regression.
9. Apply Approach B state assignment:
   - For each state group, take the maximum adjusted score among member MPs.
   - Assign the highest-scoring group.
   - Label cells `Unresolved` if the best group score is below 0.5.
   - Label cells `Hybrid` if the gap between the best and second-best group scores is below 0.3.
10. Save PDO-pipeline top MP and Approach B/noreg state assignments.
11. Plot state/MP abundance and MP activity boxplots, including and excluding PDO/SUR1090 samples.

## Caching

The script reuses these cached objects when complete:

- `Auto_parse_scATLAS_UCell_scores.rds`
- `Auto_parse_scATLAS_mp_adj_noreg.rds`
- `Auto_parse_topmp_assignments.rds`
- `Auto_parse_PDOpipeline_UCell_scores.rds`
- `Auto_parse_PDOpipeline_mp_adj_noreg.rds`
- `Auto_parse_PDOpipeline_mp_adj_all_noreg.rds`
- `Auto_parse_PDOpipeline_topmp_assignments.rds`
- Parse-derived UCell score matrices

Plot-only changes can be rerun from these cached matrices and assignments.

## Outputs

Output root:

`parse_outs/cell_states/`

Output tiers:

- `intermediate`: UCell score RDS files and adjusted-score RDS files
- `tables`: assignment CSVs, abundance summary CSVs, activity statistics CSVs
- `figures`: abundance PDFs and activity boxplot PDFs
- `logs`: `parse_outs/logs/run_summaries/parse_mp_state_abundance_activity_approachB_noreg_*.txt`

## Downstream Use

Active downstream dependencies:

- `metaprograms/parse_metaprogram_external_reference_correlation.R` consumes scATLAS/PDO UCell caches.

Current preferred state assignment:

- `parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds`
- column: `pdo_state`
- method: Approach B, noreg
