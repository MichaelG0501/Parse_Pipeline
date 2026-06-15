# State Abundance Grid Methodology

## Purpose
This workflow generates a publication-quality 2x2 grid of bar charts, illustrating the evolution of PDO-pipeline cell-state proportions across the 6 Parse treatment-response timepoints (T0, T1, T2, T4, R4, eR4). The visualization matches a specific 2x2 panel aesthetic requested for presentation or publication purposes.

## Inputs
- `parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds`: The canonical PDO-pipeline state assignments using Approach B (noreg).

## Filtering & Subsetting
1. **Timepoints**: Restricted to the 6 core Parse treatment-response timepoints (`T0`, `T1`, `T2`, `T4`, `R4`, `eR4`).
2. **Cell States**: Excludes `Unresolved` and `Hybrid` states to focus exclusively on the four assigned biological states: Classic Proliferative, Basal to Intest. Meta, Stress-adaptive, and SMG-like Metaplasia.

## Calculation
For each sample, the cell count per state is calculated. These counts are then normalized as a proportion (percentage) of the total assigned cells (excluding unresolved/hybrid) within that sample.

## Visualization
- **Layout**: A 2x2 grid created via `patchwork`.
- **Aesthetic**: Custom deep colour palette mimicking specific presentation styles. Each state is visualised in its own dedicated subplot with a fixed y-axis (0-100%).
- **Outputs**:
  - `parse_outs/publication/state_abundance_grid/figures/state_abundance_grid.pdf`
  - `parse_outs/publication/state_abundance_grid/figures/state_abundance_grid.png`
  - Source data table saved as CSV in `tables/`.

## Downstream Status
Terminal figure script. No downstream dependencies.
