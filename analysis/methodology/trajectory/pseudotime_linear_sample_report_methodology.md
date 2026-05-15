# Parse Pseudotime Linear Sample Report Methodology

Script: `analysis/trajectory/parse_pseudotime_linear_sample_report.R`

Status: active terminal trajectory plotting workflow.

## Purpose

This script projects cells onto the learned Monocle3 graph and produces a slide-facing report showing sample position, pseudotime, and sample density along pseudotime.

## Inputs

Cached trajectory assets:

- `parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime_cds.rds`
- `parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime.rds`
- `parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime_metadata.csv`

If these are missing, the helper can rebuild the trajectory.

## Processing

1. Load cached trajectory assets.
2. Extract Monocle3 principal graph structure.
3. Project cells to their nearest graph segments.
4. Compute projection distance and pseudotime summaries.
5. Build sample-colored trajectory projection panels.
6. Build pseudotime density ridges weighted by sample cell counts.
7. Save projected cell table and a multi-panel PDF report.

## Outputs

Output tiers:

- `tables`: `Auto_parse_sample_pseudotime_projections.csv`
- `reports`: `Auto_parse_pseudotime_sample_linear_report.pdf`
- `logs`: `parse_outs/logs/run_summaries/parse_pseudotime_linear_sample_report_*.txt`

## Downstream Use

Terminal figure/report workflow. No active downstream script consumes these outputs.
