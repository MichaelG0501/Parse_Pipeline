####################
# T2+T4 versus T0 DEG and pathway response methodology
####################

## Script

`analysis/cell_states/parse_t2t4_vs_t0_dge_pathway_response.R`

## Purpose

This workflow compares Parse samples T2 and T4 against T0, then turns the gene-level result into interpretable pathway and sample-level summaries. It is intended as a terminal response-analysis figure workflow rather than a dependency for state assignment.

## Inputs

- `parse_outs/by_samples/<sample>/Auto_<sample>_final.rds`
- `parse_reference_paths$cell_cycle_genes`
- `parse_reference_paths$pdo_metaprograms`
- MSigDB Hallmark gene sets from `msigdbr`

## Outputs

- `parse_outs/t2t4_vs_t0_response/tables/`: DEG tables, Hallmark FGSEA tables, sample-level pathway/metric scores, and T2/T4-minus-T0 pathway score summaries.
- `parse_outs/t2t4_vs_t0_response/figures/`: volcano plot, FGSEA dot plot, top-DEG sample heatmap, and sample pathway/metric heatmap.
- `parse_outs/t2t4_vs_t0_response/intermediate/`: cached RDS list containing DEG, FGSEA, pseudobulk, and heatmap matrices.
- `parse_outs/t2t4_vs_t0_response/reports/`: short run-facing output summary.
- `parse_outs/logs/run_summaries/`: standard session and provenance summary from `parse_pipeline_logging.R`.

## Methods

Cells from T2 and T4 are pooled into a `T2_T4` response group and compared with T0 cells using an approximate matrix-based Wilcoxon rank-sum test over log-normalized expression. The script deterministically downsamples to at most 1,500 cells per group (`set.seed(1090)`), keeps genes detected in at least 5% of either group, and computes average log2 fold-change from mean back-transformed log-normalized expression. It loads the per-sample final RDS files and avoids loading the full merged object. No batch or replicate covariate is included, matching the current assumption that the Parse samples were sequenced together without a relevant batch effect.

The DEG ranking for pathway analysis is the average log2 fold-change, with positive values indicating higher expression in T2+T4. Hallmark enrichment is run with `fgsea` over MSigDB Hallmark gene sets from `msigdbr`.

For sample-level pathway heatmaps, the script loads each per-sample final RDS file, aggregates counts by Parse sample (`T0`, `T1`, `T2`, `T4`, `R4`, `eR4`), applies edgeR TMM normalisation, converts to logCPM, row-z-scores genes across sample pseudobulks, and scores each pathway as the mean z-scored expression of its genes. The PDO-FLOT-inspired heatmap uses real sample scores, not treatment deltas. Additional metrics include:

- `CCSIG`: top expressed consensus cell-cycle genes, falling back to E2F/G2M Hallmark genes if needed.
- `Intestinal Metaplasia`: PDO-pipeline MP4 gene signature.
- `Columnar Progenitor`: PDO-pipeline MP8 gene signature.
- Composite response metrics matching the PDO script logic: cell-cycle/proliferation, injury/checkpoint, adaptive/persistence, and proteostasis/transition.

## Interpretation Notes

The DEG test treats downsampled cells as observations and uses an asymptotic Wilcoxon approximation without full zero-tie correction. It should be read as a descriptive cell-level contrast for this single time-course sample, not as replicated patient-level inference. The pathway heatmap is sample descriptive: values are real pseudobulk signature scores per sample, and the accompanying contrast summary table reports `mean(T2, T4) - T0` for ranking.
