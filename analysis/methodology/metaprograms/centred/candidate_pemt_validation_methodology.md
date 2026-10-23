# Candidate expression and pEMT metaprogram validation planning

## Purpose and scope

The terminal exploratory workflow `analysis/metaprograms/centred/parse_centred_candidate_pemt_validation.R` describes 21 preliminary candidate genes across T0, T1, T2, T4, R4 and eR4, quantifies eight pEMT/plasticity MPs, and compares candidate expression between T4 MP-high and MP-low cells. It does not establish functional validation or population-level treatment effects.

The eight workbook-defined MPs are MP48, MP50, MP60, MP109, MP16, MP74, MP83 and MP9. MP75 and MP17 are excluded from this analysis by user choice. Candidate recurrence, ranking and theme annotations are imported directly from the workbook and preserved in output tables. Workbook headings are checked against the agreed MP list.

## Exact inputs

All paths are relative to `/rds/general/project/spatialtranscriptomics/live/Parse_Pipeline/`:

- `analysis/metaprograms/centred/T2T4hgih functional annotation.xlsx`, sheets `Candidate genes` and `pEMT-Plasticity`.
- `parse_outs/by_samples/<stage>/Auto_<stage>_final.rds` for the six stages listed above.
- `parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_UCell_scores_nMP117.rds`.
- The same upstream directory's `Auto_parse_highres_T2T4_cell_metadata_nMP117.rds` and `Auto_parse_highres_T2T4_selected_mp_genes_nMP117.rds`.
- Shared configuration and run logging under `analysis/common/`.

The six final sample objects define the population. Cached score cell IDs, original cell barcodes and `orig.ident` must agree exactly; no silent dropping, cell-type filtering or alternate input fallback is allowed. RNA counts must be available as a sparse counts layer. Input paths, sizes and modification times, and the workbook MD5 are preserved in the live source bundle.

## Activity definition

Original cached UCell scores are the primary scores. For each MP, high activity is a score strictly above its T0 95th percentile, calculated with R quantile type 7. This fixed per-MP numerical cutoff is applied to every stage. Ties at the threshold are assigned low; low means not above this operational threshold, not biologically inactive.

This choice addresses treatment-associated elevation relative to pretreatment rather than selecting a fixed fraction within T4. A universal absolute cutoff is not supported by UCell's author guidance: score meaning depends on signature and data quality. T0 is not assumed to be pEMT-negative. The 95th percentile is an operational choice, not a literature-established pEMT cutoff. Sensitivity outputs use T0 q90, q97.5 and q99 as well.

Each MP's percentage uses all cells from that stage as denominator. `Any_pEMT_MP` is the union of cells high in at least one of eight MPs, with each cell counted once. MP groups overlap. The union's baseline fraction is not fixed at 5%, and the union is not a multiple-testing-adjusted classification rule. Denominators and counts accompany percentages.

## Expression and T4 comparisons

Expression is calculated consistently as `log1p(count / total RNA counts * 10000)` from the raw RNA layer. Detection means raw count greater than zero. No imputation, smoothing, integrated expression or scaled expression is used. Stage summaries include cell counts, mean and median log-normalised expression, mean CP10k and percentage detected.

Within T4, each MP and the union define high versus all remaining low cells. Every gene has mean expression, detection, high-minus-low differences, `log2((mean_CP10k_high + 1)/(mean_CP10k_low + 1))`, and Mann–Whitney rank AUC. AUC 0.5 denotes no directional separation. Groups with fewer than three cells are explicitly marked insufficient, without an alternate method.

Two-sided cell-level Wilcoxon tests are exploratory; BH correction spans 189 comparisons (21 genes by nine group definitions) separately for each scoring method. These p-values do not treat the cells as independent biological replicates for treatment inference. Emphasis is on effect sizes, detection and robustness, not p-value ranking. No biological replicate field is inferred from a barcode; metadata inventories are exported for review. Median UMI and detected-gene counts identify potential depth imbalance between groups.

## Candidate self-contribution sensitivity

UCell ranks cap the effect of very high transcript counts but each signature gene still contributes to the score. Multiple candidates belong to these MPs. Original scores remain primary, as requested. The workflow additionally scores full and leave-one-gene-out signatures together in the same installed UCell version. Each candidate is removed separately from each signature for a leave-one-gene-out check; candidates remain in the rank universe and expression matrix. Each tested gene has its own high/low sensitivity grouping. Removing all candidates together is not used: MP83 contains 10 candidates among 14 signature genes, so only four genes would remain. The leave-one-out design preserves signature identity while removing the tested gene's direct score contribution; other candidates and correlated genes still contribute.

Rescoring uses only T0 and T4, `maxRank=1500`, four cores and 500-cell chunks, with a gene universe equal to the intersection of genes across the six RNA count matrices. Each signature is intersected with this universe. Reduced signatures must retain at least five genes. Full and gene-specific leave-one-out thresholds are independently calibrated to T0 q95. Full-signature rescoring is compared with the cached scores to expose version or implementation differences. Tables include classification agreement, prevalence, expression effects and signature overlap.

The prioritisation table flags positive mean differences in both original and leave-one-gene-out union comparisons, then orders by leave-one-gene-out AUC. This does not establish specificity to pEMT, causality, protein expression, druggability or independent validation. Signatures and candidates were selected using this experiment; selection bias remains.

## Outputs and reproducibility

The one current live output root is `parse_outs/centred/candidate_pemt_validation/`:

- `intermediate/candidate_pemt_replot_bundle.rds`: compact critical source data including per-cell candidate expression/detection, original scores, T0/T4 sensitivity scores, metadata, signatures, annotations and provenance. This file is a required persistent plotting input and must remain in live storage.
- `tables/`: expression, score summaries, baseline cutoffs, high-cell fractions, all cell classifications, T4 expression contrasts, signature overlap, QC and prioritisation.
- `figures/Auto_final/`: four selected figures for expression, activity distributions, prevalence and T4 expression contrasts.
- `figures/Auto_supporting/`: gene trajectories, leave-one-gene-out sensitivity and threshold sensitivity.
- `reports/`: current results, a compact linked figure index and validation results.
- Shared run summaries: `parse_outs/logs/run_summaries/parse_centred_candidate_pemt_validation_<timestamp>.txt`.

The exploratory figures use wide R/ggplot2 PDF exports with readable 13-point base text and PNG companions generated from the same plot. Violin/boxplot panels show continuous scores and thresholds; figures avoid significance stars that could exaggerate cell-level inference.

`--mode=analyse` reads upstream inputs, computes expression and sensitivity scores, and refuses to overwrite an existing source bundle. `--mode=replot` uses the persistent bundle alone, without reading the workbook, sample objects or upstream scores. It rewrites only this workflow's derived outputs. No upstream files are modified, moved or deleted. The script forces the live project root because the existing shared configuration otherwise chooses ephemeral storage.

## Execution registry

Use the existing `gnmf` conda environment. Resources: four CPUs, 48 GB RAM, four hours. Submission-only launcher: `job_scripts/metaprograms/centred/parse_centred_candidate_pemt_validation.sh` (ignored). Submit from the project root with `qsub` using that path. PBS stdout/stderr are streamed to the project root using `#PBS -koed` and absolute `-o`/`-e` directories. The R script runs directly on a compute node without reading the launcher.

No active downstream workflow depends on these results. The upstream centred discovery/filtering scripts are not rerun or changed. The workflow and its current output dependencies are recorded in AGENTS.md and `analysis/script_inventory_and_dependency_map.md`.

## References

- Andreatta and Carmona, UCell (2021): https://doi.org/10.1016/j.csbj.2021.06.043
- Author guidance on cutoff interpretation: https://github.com/carmonalab/UCell/issues/7
- Current implementation and normalization changes: https://github.com/carmonalab/UCell

An inherited `_R_CHECK_LIMIT_CORES_` environment variable caps BiocParallel at two workers. Inside the four-CPU PBS allocation the script unsets that package-check flag before scoring; outside PBS it stops if the flag is present. BLAS/OpenMP threads are limited to one by the launcher.

## Figure 1 trend heatmap

Figure 1 now displays genes as rows and T0, T1, T2, T4, R4 and eR4 as ordered columns. Colours encode `(stage mean - mean of six stage means) / sample SD of six stage means`, calculated independently for each gene without weighting stages by cell count. The symmetric blue-white-red scale uses the observed maximum absolute z-score without clipping. Printed values are unscaled mean log1p(CP10k) expression, including zero-expression cells, rounded to two decimals; these are not raw UMI counts. Constant or nonfinite rows stop execution. `tables/candidate_expression_heatmap_source.csv` preserves exact expression means, z-scores and scaling parameters.

`--mode=figure1` refreshes only Figure 1 PDF/PNG, its source table and the current figure index from the persistent live bundle. The mode remains available for scoped heatmap refreshes. The initial dot plot was copied to the ignored session archive before replacement; no input data or other figures are changed.

## Figure 04 continuous gene–MP association page

Figure 04 has two pages. Page 1 preserves the original T4 MP-high versus MP-low mean log1p(CP10k) differences. Page 2 shows all 168 T4 candidate–MP Spearman correlations, using the tested candidate's leave-one-gene-out score from the persistent bundle. All 2,600 T4 cells, including zero-expression cells, contribute; ties use average ranks. No binary activity cutoff enters the second page. Constant or nonfinite inputs stop execution.

Page 2 rows rank candidates by their median correlation across eight MPs; an explicitly labelled final summary column reports that median. The summary represents correlated programs rather than an independent ninth MP. A fixed symmetric blue–white–red scale spans rho -1 to +1, and each tile prints rho to two decimals. Near-zero rho means weak monotonic association, not proof of no relationship. Sparse detection, nonmonotonic responses, technical depth and other cell properties can affect this unadjusted association. Positive association does not establish causality or functional activation.

`tables/T4_candidate_continuous_mp_associations.csv` records original-score and leave-one-out-score rho, cell counts and gene detection for every pair. `tables/T4_candidate_continuous_association_summary.csv` contains median/minimum/maximum rho and positive-MP counts for each gene. `reports/figure4_continuous_association.md` explains interpretation. Validation independently compares Spearman rho with Pearson correlation of average ranks. The existing Figure 04 PNG remains the page 1 companion; `04_T4_candidate_high_vs_low_page2.png` previews page 2. There is one current two-page PDF, without duplicate standalone PDFs.

Use `--mode=figure4` to refresh only these outputs and the current figure index from the live bundle. The matching ignored launcher currently uses this mode. No upstream rescoring is needed. The old single-page PDF/PNG were copied into the ignored session archive before replacement.

## Supporting per-cell scatter plots

`figures/Auto_supporting/T4_candidate_per_cell_pemt_scatter.pdf` contains 21 pages in candidate-list order, one gene per page with all eight MPs arranged in a 4-by-2 grid. Every panel plots all 2,600 T4 cells against the candidate-specific leave-one-gene-out UCell score, matching Figure 04 page 2. Zero-expression cells are retained; no jitter or cell subsampling is used. Transparency exposes overlap, but coincident observations remain coincident. Y-axis limits are shared across the eight panels of a gene; x-axis ranges follow each MP's observed scores. Panel labels report Spearman rho.

Orange points connect mean expression and mean MP score within MP-score quantile bins, using decile boundaries with duplicate boundaries collapsed. All cells are included in exactly one bin per gene–MP pair. This is a descriptive mean trend, not a regression fit, confidence interval or causal model. The bin-mean table and page index live under `tables/`. A single page-1 PNG provides a readable preview without duplicating every page as a separate figure. The live source bundle remains sufficient to regenerate all scatter pages.

`--mode=scatter` refreshes these supporting outputs alone and adds their links to the current index. The matching local PBS launcher currently selects this mode. Full analyse/replot modes also produce the scatter PDF. The final Figure 04 retains two compact summary pages; the scatter collection is supporting evidence.
