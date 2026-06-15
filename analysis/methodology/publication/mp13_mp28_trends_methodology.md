# Methodology: Publication Trend Plots for MP13 and MP28

## Purpose
This workflow generates publication-quality trend plots for two specific metaprograms (MPs) identified during high-resolution analysis:
- **MP13** (Stress Response), selected from the T2/T4-high legacy filter.
- **MP28** (Stress-Associated Proliferation), selected from the strict mean/median trend filter.

## Inputs
- `parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_sample_ucell_summary_nMP117.csv` (for MP13)
- `parse_outs/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_sample_ucell_summary_nMP117.csv` (for MP28)

## Processing Steps
1. The script loads the aggregated mean and median UCell scores per sample (T0, T1, T2, T4, R4, eR4) from the existing high-resolution summary CSV files.
2. It filters out the specific MPs (MP13 and MP28).
3. The data is pivoted to a long format to plot both Mean (solid line) and Median (dashed line).
4. `ggplot2` is used with a Nature-contract theme:
   - 7pt and 8pt Arial fonts.
   - Minimalist classic theme.
   - 60x60mm dimensions.
   - Points colored by sample using the global `parse_sample_colours`.

## Outputs
- `parse_outs/publication/mp_trends/figures/parse_mp13_trend.pdf` (and `.png`)
- `parse_outs/publication/mp_trends/figures/parse_mp28_trend.pdf` (and `.png`)

## Downstream Usage
These standalone plots are intended for direct use in publication figures and do not have downstream script dependencies.
