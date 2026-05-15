# Shared Configuration and Logging Methodology

This document describes the common files in `analysis/common/`.

## Purpose

The cleanup centralized settings that were repeated across scripts:

- Project-root detection
- Common Parse sample order
- Preferred state definition
- Output-tier folder names
- Plotting defaults suitable for PowerPoint slides
- Reference file paths
- Frequently reused thresholds
- Lightweight run-summary logging

## Configuration Files

`parse_pipeline_config.R` defines:

- `parse_project_root()`: resolves the project root from `AUTO_PARSE_ROOT_DIR` or known Imperial HPC paths.
- `parse_paths()`: returns common paths such as `parse_outs`, `analysis_dir`, `logs`, `figures`, and `methodology`.
- `parse_output_tiers(base_dir)`: returns and optionally creates `intermediate/`, `tables/`, `figures/`, `logs/`, and `reports/`.
- `parse_samples`: `T0`, `T1`, `T2`, `T4`, `R4`, `eR4`.
- `parse_all_samples`: Parse samples plus `PDO`, `SUR1090_Untreated`, and `SUR1090_Treated`.
- `parse_state_definition`: active Approach B/noreg state assignment settings with `min_group_score = 0.5` and `hybrid_gap = 0.3`.
- `parse_plot_defaults`: DPI, PDF sizes, and slide-readable font/legend/point defaults.
- `parse_reference_paths`: external scATLAS, PDO, 3CA, cell-cycle, developmental, CNV, and velocity reference paths.

`parse_pipeline_config.py` mirrors project-root, sample, and velocity constants for Python scripts.

## Helper Functions

`parse_pipeline_helpers.R` provides reusable functions:

- `parse_load_or_stop()`
- `parse_get_counts()`
- `parse_chunk_vector()`
- `parse_write_mp_gene_table()`
- `parse_z_normalise()`
- `parse_slide_theme()`

These helpers should replace copy-pasted equivalents when scripts are next refactored.

## Run Summary Logging

`parse_pipeline_logging.R` provides:

- `parse_start_run()`
- `parse_finish_run()`

Scripts initialize a `script_run` object near the top, set `script_run_status <- "failed"`, register an `on.exit()` summary write, and set `script_run_status <- "success"` only after the final output is written.

Run summaries are written to:

`parse_outs/logs/run_summaries/<script>_<timestamp>.txt`

Each summary records:

- Start and end time
- Elapsed minutes
- Input files
- Output files
- Parameters
- Cache reuse notes when supplied
- `sessionInfo()` for R scripts

## Plotting Standard

Final or slide-facing scripts should use readable plot sizing:

- Prefer wide PDF pages for multi-panel summaries.
- Avoid small legends and row names; target at least 12-14 pt axis/legend text for dense slides.
- Save plot data or input RDS files next to expensive figures so visual styling can be regenerated without recomputing upstream results.
