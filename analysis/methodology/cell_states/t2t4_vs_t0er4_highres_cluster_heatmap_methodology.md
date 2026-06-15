####################
# t2t4_vs_t0er4_highres_cluster_heatmap_methodology.md
#
# Methodology: High-resolution MP cluster comparison (T2T4 vs T0eR4)
####################

## Purpose
Visualise the changes in high-resolution MP cluster scores (derived from the PDO dataset, nMP156) between treated timepoints (T2, T4) and baseline/recovery timepoints (T0, eR4) in the Parse dataset.

## Approach
1. **MP Definition**: We load the nMP=156 metaprogram gene lists from the PDO pipeline (`Auto_pdo_flot_highres_geneNMF_metaprograms_nMP156.rds`).
2. **Scoring**: We score the individual Parse sample objects (`T0`, `T2`, `T4`, `eR4`) using `UCell::AddModuleScore_UCell`.
3. **State Mapping**: We merge the scores with the global cell metadata (`Auto_parse_all_meta.rds`) to assign the `pdo_state` (Approach B, noreg) to each cell. "3CA_mp_12 Protein maturation" and "3CA_mp_17 EMT III" are combined into a single `3CA_EMT_and_Protein_maturation` state.
4. **Aggregation**: Cells are grouped by State and Treatment (`T2T4` vs `T0eR4`). We calculate the mean delta (T2T4 - T0eR4) for each manually curated functional cluster.
5. **Visualisation**: We produce two heatmaps:
   - A Delta Heatmap showing the score change (T2T4 - T0eR4) for each cluster across the states.
   - An Absolute Score Heatmap showing the raw scores side-by-side (T0eR4 | T2T4).
   - An MP-level Delta Heatmap showing the delta for each individual MP within the clusters.
