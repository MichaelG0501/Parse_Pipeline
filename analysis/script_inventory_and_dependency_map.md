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
- `cell_states/parse_mp_state_abundance_activity_approachB_noreg.R`
- `cnv/parse_infercna_carroll_reference.R`
- `trajectory/parse_pseudotime_linear_sample_report.R`
- `trajectory/parse_pseudotime_sample_distance_maps.R`
- `trajectory/parse_velocity_scvelo_visualise.py`
- `plotting/parse_qc_marker_heatmap.R`
- `summary/parse_filtering_summary_plots.R`

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
