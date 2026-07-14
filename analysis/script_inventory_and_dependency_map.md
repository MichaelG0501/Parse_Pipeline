# Parse Analysis Script Inventory and Dependency Map

This map is the working index for `analysis/`. Update it whenever a script is added, renamed, moved, marked legacy, or made a downstream dependency.

## Current Conventions

- Active scripts use informative names without the old `Auto_` script prefix.
- `legacy_` scripts are retained for comparison only and must not feed active downstream analysis.
- `delete_` scripts are redundant compatibility wrappers or duplicates. They are kept only for manual deletion by the user.
- Active state definition: PDO-pipeline state assignment, Approach B, noreg normalization.
- Shared constants and helpers live in `analysis/common/`.
- Methodology files live under `analysis/methodology/<matching_analysis_subfolder>/`.
- Existing output filenames under `parse_outs/` are kept stable unless a user explicitly approves migration.

## Shared Configuration

| Location | Purpose |
| :--- | :--- |
| `analysis/common/parse_pipeline_config.R` | R project paths, sample order, preferred state definition, output tiers, plotting defaults, thresholds, reference paths |
| `analysis/common/parse_pipeline_helpers.R` | Reusable R helpers for package checks, counts extraction, chunking, MP gene tables, z-normalization, slide-scale ggplot theme |
| `analysis/common/parse_pipeline_logging.R` | Lightweight run summary files in `parse_outs/logs/run_summaries/` |
| `analysis/common/parse_pipeline_config.py` | Python project paths, velocity paths, and sample constants |

## Run Order

| Step | Script | Depends on | Main outputs | Downstream consumers |
| :--- | :--- | :--- | :--- | :--- |
| 0 | `input_preparation/parse_prepare_sur1090_count_matrices.R` | External SUR1090 CSV matrices | `parse_outs/input/output_combined/<sample>/DGE_filtered/` | QC pipeline before this analysis tree |
| 1 | QC pipeline outside `analysis/` | Raw Parse/Trailmaker/SUR1090 inputs | `parse_outs/by_samples/*/Auto_*_final.rds`, `Auto_parse_merged.rds`, summary CSVs | All downstream scripts |
| 2 | `metaprograms/parse_metaprogram_geneNMF_discovery.R` | final per-sample Seurat objects | GeneNMF sweep RDS files; selected UCell after step 3 exists | nMP selection, enrichment, MP correlations, cell states |
| 3 | `metaprograms/parse_metaprogram_select_optimal_nMP.R` | GeneNMF sweep outputs | `Auto_parse_optimal_nMP.txt`, `Auto_parse_MP_outs_default.rds` | GeneNMF scoring rerun, enrichment, correlations, cell states |
| 4 | `metaprograms/parse_metaprogram_geneNMF_discovery.R` rerun | selected nMP text from step 3 | selected UCell scores and final GeneNMF Seurat object | MP correlation, cell-state scoring |
| 5 | `metaprograms/parse_metaprogram_enrichment_annotation.R` | selected MP object | enrichment RDS/PDF/PNG | terminal figures |
| 6 | `metaprograms/parse_metaprogram_internal_correlation.R` | selected MP object and Parse UCell | internal MP correlation PDF/summary | terminal figures |
| 7 | `cell_states/parse_mp_state_abundance_activity_approachB_noreg.R` | selected Parse MP object; scATLAS and PDO metaprogram references | preferred Approach B/noreg state assignments, UCell caches, abundance/activity PDFs | external MP correlation, terminal state figures |
| 8 | `metaprograms/parse_metaprogram_external_reference_correlation.R` | selected Parse UCell; cell-state scATLAS/PDO UCell caches | Parse-vs-scATLAS/PDO Jaccard and correlation PDFs | terminal figures |
| 9 | `cnv/parse_infercna_carroll_reference.R` | final per-sample Seurat objects; Carroll_2023 reference; gene order | InferCNA RDS, CNA heatmaps, CNA scatter | terminal CNV figures |
| 10 | `trajectory/parse_pseudotime_sample_trajectory.R` | `Auto_parse_merged.rds` | Monocle3 CDS, pseudotime vector, metadata, combined PDF | pseudotime linear report and distance maps |
| 11 | `trajectory/parse_pseudotime_linear_sample_report.R` | cached pseudotime assets | linear sample trajectory PDF and projections CSV | terminal trajectory figure |
| 12 | `trajectory/parse_pseudotime_sample_distance_maps.R` | cached pseudotime assets | distance matrices, heatmap, T0 nodeplot | terminal trajectory figures |
| 13 | `trajectory/parse_velocity_submit_pipeline.sh` | Trailmaker BAMs and metadata | velocity BAMs, looms, scVelo figures/tables | terminal velocity figures |

## Terminal Figure-Generation Scripts

These scripts produce final or presentation-facing outputs and generally do not feed other active workflows:

- `metaprograms/parse_metaprogram_enrichment_annotation.R`
- `metaprograms/parse_metaprogram_internal_correlation.R`
- `metaprograms/parse_metaprogram_external_reference_correlation.R`
- `metaprograms/parse_highres_mp_tcga_survival_volcano.R`
- `cell_states/parse_mp_state_abundance_activity_approachB_noreg.R`
- `cnv/parse_infercna_carroll_reference.R`
- `trajectory/parse_pseudotime_linear_sample_report.R`
- `trajectory/parse_pseudotime_sample_distance_maps.R`
- `trajectory/parse_velocity_scvelo_visualise.py`
- `plotting/parse_qc_marker_heatmap.R`
- `summary/parse_filtering_summary_plots.R`
- `publication/parse_state_abundance_timepoint.R`
- `publication/parse_highres_metaprogram_heatmap.R`

## Legacy and Delete Candidates

| Script | Status | Reason |
| :--- | :--- | :--- |
| `metaprograms/legacy_parse_highres_mp_t2t4_comparison_filter.R` | Legacy comparison | Alternative T2/T4-high high-resolution MP filter. Outputs are isolated under `Auto_T2T4_gt_T0eR4_filter/` and must not be used downstream. |
| `trajectory/delete_Auto_pseudotime_linear_plot_legacy_wrapper.R` | Delete candidate | Old compatibility wrapper for the linear pseudotime script. |
| `trajectory/delete_Auto_pseudotime_state_distance_matrix_legacy_wrapper.R` | Delete candidate | Old compatibility wrapper with misleading "state" name; active analysis is sample distance. |
| `trajectory/delete_parse_velocity_submit_pipeline_run2_duplicate.sh` | Delete candidate | Duplicate velocity submission wrapper. |

## Outdated Pointer Check

- Root-level PBS wrappers are currently untracked and were updated in-place to call the renamed scripts, but they are intentionally not staged.
- No active tracked script should call `analysis/metaprograms/Auto_parse_mp_correlation_with_pancancer.R`; that old path did not exist. The active external comparison script is `analysis/metaprograms/parse_metaprogram_external_reference_correlation.R`.
- Existing output filenames still begin with `Auto_parse_*` to preserve downstream compatibility.

## Output Tiers

New or substantially revised workflows should use this tiering inside their workflow output directory:

- `intermediate/`: reusable expensive objects such as RDS, H5AD, loom, or cached matrices
- `tables/`: CSV/TSV/XLSX summaries
- `figures/`: PNG/PDF/SVG presentation outputs
- `logs/`: run summaries, PBS logs, and package/session information
- `reports/`: multi-panel PDFs, HTML reports, or slide-facing summaries

Existing scripts may retain legacy output locations, but their headers and methodology files must classify which outputs act as intermediate caches versus terminal figures.

####################
## 2026-05-15 Tumour Stiffness Gene Module Heatmap Update

- Added active terminal figure script `cell_states/parse_tumour_stiffness_gene_module_heatmap.R`.
- Methodology: `analysis/methodology/cell_states/tumour_stiffness_gene_module_heatmap_methodology.md`.
- Depends on `parse_outs/Auto_parse_merged.rds`.
- Main outputs live under `parse_outs/cell_states/tumour_stiffness_gene_module/` with `intermediate/`, `tables/`, `figures/`, and `reports/` tiers.
- Terminal figures: `figures/Auto_parse_tumour_stiffness_module_gene_heatmap.pdf` and `figures/Auto_parse_tumour_stiffness_module_gene_heatmap_no_pdo_sur1090.pdf`.
- No active downstream consumers.
####################

####################
## 2026-05-15 T2+T4 versus T0 DEG and Pathway Response Update

- Added active terminal response script `cell_states/parse_t2t4_vs_t0_dge_pathway_response.R`.
- Methodology: `analysis/methodology/cell_states/t2t4_vs_t0_dge_pathway_response_methodology.md`.
- Depends on `parse_outs/by_samples/<sample>/Auto_<sample>_final.rds`, MSigDB Hallmark gene sets via `msigdbr`, the configured cell-cycle gene file, and the PDO-pipeline metaprogram reference for MP4/MP8 lineage metrics.
- Main outputs live under `parse_outs/t2t4_vs_t0_response/` with `intermediate/`, `tables/`, `figures/`, `logs/`, and `reports/` tiers.
- Terminal figures include the T2+T4 vs T0 volcano plot, Hallmark FGSEA dot plot, top-DEG sample heatmap, and real-value sample pathway/metric heatmap.
- No active downstream consumers.
####################

####################
## 2026-05-26 High-Resolution MP TCGA Survival Volcano Update

- Added active terminal figure script `metaprograms/parse_highres_mp_tcga_survival_volcano.R`.
- Methodology: `analysis/methodology/metaprograms/highres_mp_tcga_survival_volcano_methodology.md`.
- Depends on strict high-resolution MP selected genes/trend summary, legacy T2/T4-high selected genes, and the scRef TCGA ESCA metadata plus whole-profile TPM compatibility copies.
- Main outputs live under `parse_outs/highres_mp_tcga_survival/` with `intermediate/`, `tables/`, `figures/`, and `reports/` tiers.
- Terminal figures include the multi-page whole-TCGA volcano PDF and individual strict-increase, strict-decrease, and legacy-T2/T4-high volcano PDFs/PNGs for continuous, median, and q1q4 Cox variants.
- No active downstream consumers.
####################

####################
## 2026-06-08 Publication State Abundance Timepoint Update

- Added active terminal figure script `publication/parse_state_abundance_timepoint.R`.
- Methodology: `analysis/methodology/publication/state_abundance_timepoint_methodology.md`.
- Depends on `parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds` (produced by step 7).
- Main outputs live under `parse_outs/publication/state_abundance_timepoint/` with `tables/` and `figures/` tiers.
- Terminal figures: `figures/state_abundance_timepoint.pdf` and `.png` — stacked bar of PDO-pipeline state abundance across 6 Parse timepoints, Unresolved/Hybrid excluded.
- No active downstream consumers.
####################

####################
## 2026-06-08 Publication UMAP Samples Update

- Added active terminal figure script `publication/parse_umap_samples.R`.
- Methodology: `analysis/methodology/publication/umap_samples_methodology.md`.
- Depends on `parse_outs/Auto_parse_merged.rds`.
- Main outputs live under `parse_outs/publication/umap_samples/` with `figures/` tier.
- Terminal figures: `figures/umap_samples.pdf` and `.png` — UMAP visualization of the 6 Parse treatment-response timepoints plus PDO, colored by sample.
- No active downstream consumers.
####################

####################
## 2026-06-08 High-Resolution Metaprogram Heatmap Update

- Added active terminal figure script `publication/parse_highres_metaprogram_heatmap.R`.
- Methodology: `analysis/methodology/publication/highres_metaprogram_heatmap_methodology.md`.
- Depends on `parse_outs/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds` and `parse_outs/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_geneNMF_metaprograms_nMP<k>.rds`.
- Main outputs live under `parse_outs/publication/highres_metaprogram_heatmap/` with `figures/` tier.
- Terminal figures: `figures/highres_metaprogram_heatmap.pdf` and `.png` — ComplexHeatmap visualization of the high-resolution metaprogram splits.
- No active downstream consumers.
####################

### `parse_state_abundance_grid.R`
*   **Location**: `analysis/publication/`
*   **Purpose**: Generates a 2x2 grid of bar charts showing PDO-pipeline cell-state proportions across the 6 Parse timepoints.
*   **Run Status**: Terminal figure generation.
*   **Output Tiers**: `figures/`, `tables/`
*   **Legacy Status**: Active

### `parse_publication_mp13_mp28_trends.R`
*   **Location**: `analysis/publication/`
*   **Purpose**: Generates Nature-style trend plots for MP13 and MP28 across 6 Parse timepoints.
*   **Run Status**: Terminal figure generation.
*   **Output Tiers**: `figures/`
*   **Legacy Status**: Active

####################
## 2026-06-25 Timepoint SCENIC Regulon Workflow Update

- Added active terminal workflow script `cell_states/parse_timepoint_scenic_regulons.R`.
- Added PBS wrapper `parse_timepoint_scenic.sh`.
- Methodology: `analysis/methodology/cell_states/timepoint_scenic_regulons_methodology.md`.
- Depends on `parse_outs/Auto_parse_merged.rds` and the hg38 refseq-r80 mc9nr RcisTarget databases under `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/cistarget_databases_rcistarget_mc9nr/`.
- Main outputs live under `parse_outs/cell_states/timepoint_scenic/` with `intermediate/`, `tables/`, `figures/`, `logs/`, and `reports/` tiers.
- The workflow intentionally ignores states and metaprograms. It treats `T0`, `T1`, `T2`, `T4`, `R4`, and `eR4` as the analysis entities in one combined six-timepoint SCENIC run, yielding a shared regulon dictionary for comparable regulon AUC/RSS/differential activity across timepoints. Independent per-timepoint SCENIC runs are disabled because they create non-comparable regulon target dictionaries.
- Heatmap outputs use strict timepoint column order `T0`, `T1`, `T2`, `T4`, `R4`, `eR4`; the top-specificity companion output selects the top 20 RSS-specific regulons per timepoint and orders them by decreasing specificity within each timepoint.
- GENIE3 resumes from saved `int/1.3_GENIE3_weightMatrix_part_*.Rds` files with `resumePreviousRun = TRUE`; default `genie3_nparts = 100` keeps continuation chunks smaller after HPC walltime kills.
- No active downstream consumers.
####################

####################
## 2026-06-30 Centred GeneNMF Method-Comparison Workflow Update

- Added active comparison scripts under `metaprograms/centred/`: `parse_centred_metaprogram_geneNMF_discovery.R`, `parse_centred_highres_mp_strict_mean_median_trend_filter.R`, `parse_centred_highres_mp_t2t4_comparison_filter.R`, `parse_centred_t2t4_vs_t0er4_highres_cluster_heatmap.R`, `parse_compare_centred_vs_uncentred_highres_mps.R`, and `parse_centred_publication_highres_figures.R`.
- Added PBS wrapper `parse_centred_metaprogram_workflow.sh`.
- Methodology notes were appended to the original method files rather than maintaining separate centred methodology files: `metaprogram_geneNMF_discovery_methodology.md`, `highres_mp_strict_mean_median_trend_filter_methodology.md`, `legacy_highres_mp_t2t4_comparison_filter_methodology.md`, and `publication/highres_metaprogram_heatmap_methodology.md`.
- Main outputs live under `parse_outs/centred/` with high-resolution trend, T2/T4-high filter, publication figure, comparison table, figure, and report subdirectories.
- The workflow compares centred `multiNMF(center = TRUE)` retained MPs against existing uncentred retained MPs by gene-set Jaccard best match and automatic top 3CA non-cell-cycle annotation agreement. The high-resolution nMP is computed as half of the total NMF programmes from `T0`, `T1`, `T2`, `T4`, `R4`, and `eR4` (currently 234 / 2 = nMP117), not from the external PDO nMP156 object. It is terminal comparison work and does not replace canonical uncentred Parse MPs.
####################
