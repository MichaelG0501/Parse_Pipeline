####################
# Timepoint SCENIC Regulon Methodology
####################

## Purpose

`analysis/cell_states/parse_timepoint_scenic_regulons.R` runs SCENIC on Parse treatment-response samples using timepoint as the only biological grouping variable. It intentionally does not load or use Parse/PDO state assignments, metaprograms, MP scores, or state labels.

## Inputs

- `parse_outs/Auto_parse_merged.rds`
- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/cistarget_databases_rcistarget_mc9nr/`
- Required timepoints in `orig.ident`: `T0`, `T1`, `T2`, `T4`, `R4`, `eR4`

The selected database directory is `cistarget_databases_rcistarget_mc9nr`, matching the PDO SCENIC workflow default and containing the hg38 refseq-r80 mc9nr RcisTarget feather databases:

- `hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather`
- `hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather`

## Method

The workflow uses a combined timepoint analysis: one SCENIC network is inferred from cells across all six Parse timepoints. Regulon AUC and regulon specificity scores are summarized by timepoint, giving directly comparable T0/T1/T2/T4/R4/eR4 activity matrices and differential activity tables.

Independent per-timepoint SCENIC runs are intentionally disabled because they infer separate regulon dictionaries. A TF can therefore have different target sets in different timepoint-specific runs, making downstream regulon activity comparison less clean for a time-course experiment.

By default, all available cells from the six Parse timepoints are used. A `cells_per_timepoint` argument can be supplied to use a reproducible balanced downsample.

The default combined regulon heatmap selects the union of the top RSS-specific regulons for each timepoint and plots row-scaled RSS values. The balanced top-specificity heatmap recomputes RSS after equal downsampling to the smallest timepoint cell count, selects the top 20 RSS regulons per timepoint, keeps duplicate regulons when selected for more than one timepoint, orders rows within each timepoint by decreasing RSS, and fixes columns as `T0`, `T1`, `T2`, `T4`, `R4`, `eR4`. The balanced top-gap heatmap selects regulons whose RSS is highest in that timepoint and ranks them by the gap to the next-highest timepoint.


## Outputs

Main output directory: `parse_outs/cell_states/timepoint_scenic/`

- `intermediate/`: combined SCENIC `int/` folder, regulon AUC RDS file, regulon target RDS file, mean-AUC matrix, and RSS matrix.
- `tables/`: selected cell table, database confirmation, mean-AUC/RSS CSV files, differential regulon activity tables, network edges, regulon target previews, top-20 specificity regulon tables, and run summaries.
- `figures/`: scaled RSS heatmaps, top-20 specificity RSS heatmaps, and timepoint-regulon network PDFs.
- `logs/`: reserved for workflow-local logs.
- `reports/`: reserved for later report assembly.

## Cache Behavior

SCENIC intermediate objects under `intermediate/scenic_combined/` are reused when present. Summary tables and figures are regenerated from the cached SCENIC outputs.

GENIE3 is explicitly run with `resumePreviousRun = TRUE` by default. Existing `int/1.3_GENIE3_weightMatrix_part_*.Rds` files are loaded and their target genes are skipped, so a continuation job starts from the saved parts instead of restarting the network inference. The default `genie3_nparts = 100` creates smaller remaining chunks, reducing lost work if an HPC walltime limit interrupts GENIE3 again. The final `int/1.4_GENIE3_linkList.Rds` is written only after all GENIE3 parts are complete.

## Downstream Status

Active terminal workflow. No active downstream script depends on these outputs yet.
