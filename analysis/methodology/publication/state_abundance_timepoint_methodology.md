# State Abundance Timepoint Publication Figure — Methodology

## Purpose

Publication-quality stacked bar chart showing the proportional distribution of PDO-pipeline cell states across the 6 Parse treatment-response timepoints (T0, T1, T2, T4, R4, eR4). Designed for manuscript inclusion following Nature figure contract standards.

## Script

`analysis/publication/parse_state_abundance_timepoint.R`

## Input

| File | Description |
| :--- | :--- |
| `parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds` | Pre-computed PDO-pipeline Approach B / noreg state assignments for all 9 samples |

## Filtering

1. **Samples**: Only the 6 Parse timepoints are retained (T0, T1, T2, T4, R4, eR4). PDO, SUR1090_Treated, and SUR1090_Untreated are excluded.
2. **States**: Unresolved and Hybrid cells are excluded, leaving 4 biologically assigned states:
   - Classic Proliferative
   - Basal to Intest. Meta
   - Stress-adaptive
   - SMG-like Metaplasia

## Proportions

Percentages are computed among the resolved cells only (i.e., after excluding Unresolved and Hybrid). This means all bars sum to 100%.

## Visual Design

- **Plot type**: Stacked bar chart with white segment borders
- **Proportion only**: Shows proportion of resolved cell states only; no secondary cell-count axis or point/line overlay is displayed
- **Colour palette**: Original colors matching the source pipeline:
  - Classic Proliferative: `#E41A1C` (red)
  - Basal to Intest. Meta: `#4DAF4A` (green)
  - Stress-adaptive: `#984EA3` (purple)
  - SMG-like Metaplasia: `#FF7F00` (orange)
- **Typography**: 6.5pt Arial base (Nature contract), bold axis titles, with enlarged bold 9pt x-axis labels
- **Export**: `cairo_pdf` for editable text embedding; 600 DPI PNG; dimensions 100 × 80 mm (approximately single-column Nature width)
- **Legend**: Single-row at bottom, ordered to match the stack order

## Outputs

| File | Tier | Description |
| :--- | :--- | :--- |
| `parse_outs/publication/state_abundance_timepoint/figures/state_abundance_timepoint.pdf` | Figures | Main terminal PDF |
| `parse_outs/publication/state_abundance_timepoint/figures/state_abundance_timepoint.png` | Figures | 600 DPI raster copy |
| `parse_outs/publication/state_abundance_timepoint/tables/state_abundance_timepoint.csv` | Tables | Source data: timepoint, state, cell count, proportion |

## Relationship to Existing Pipeline

This script is a lightweight publication derivative of page 3 of `Auto_parse_topmp_abundance_PDOpipeline_9samples.pdf` produced by `analysis/cell_states/parse_mp_state_abundance_activity_approachB_noreg.R` (step 7). It reads from the same cached state assignment RDS and does not recompute UCell scores.

## Downstream

Terminal figure — no active downstream consumers.
