# Parse scRNA-seq QC work report

Date: 2026-04-23

Input archive:
`/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline/8da03a32-a9ad-4938-9338-cbd1731f445d_filtered_matrices.zip`

Working directory:
`/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline`

Environment:
`/rds/general/user/sg3723/home/anaconda3/envs/dmtcp`

Samples processed:
`PDO`, `T0`, `T1`, `T2`, `T4`, `R4`, `eR4`

## Notes

The Trailmaker archive is a filtered Parse matrix export, not a raw FASTQ-level input. The QC workflow here therefore starts from the per-sample filtered count matrices and cell metadata already produced by Trailmaker.

## Workflow

1. Unpacked the Trailmaker archive into `parse_outs/input/output_combined`.
2. Read each sample from `count_matrix.mtx.gz`, `all_genes.csv.gz`, and `cell_metadata.csv.gz`.
3. Built one Seurat object per sample and computed mitochondrial percentage from `^MT-` genes.
4. Generated inspection violin plots for `nFeature_RNA`, `nCount_RNA`, and `percent.mt`.
5. Ran DoubletFinder per sample after log-normalization, HVG selection, scaling, PCA, neighbor graph, clustering, and UMAP.
6. Chose one QC threshold set for all seven samples after checking the observed QC distributions:
   - `percent.mt < 15`
   - `nFeature_RNA >= 2500`
   - `nFeature_RNA <= 13000`
   - average housekeeping expression `>= 3`
7. Saved post-doublet/post-mito prefilter objects and final filtered objects per sample.
8. Merged the final filtered objects, then reran normalization, HVG selection, scaling, PCA, neighbors, clustering, and UMAP on the combined object.
9. Generated summary plots, per-sample filtering plots, merged UMAP plots, and the QC heatmap adapted from the snSeq heatmap style for this epithelial-focused cohort.

## Doublet settings

Expected doublet rates were assigned from final loaded cell counts:

- `PDO`, `T0`, `T1`, `T4`, `R4`, `eR4`: `4%`
- `T2`: `8%`

`pK` was selected per sample from the maximum `BCmetric` in the DoubletFinder parameter sweep, then adjusted expected doublets were corrected by the homotypic proportion.

## Filtering summary

Overall:

| Step | Cells |
|---|---:|
| Raw input cells | 28075 |
| Doublet singlets | 26973 |
| After mito filter | 26932 |
| After gene-count filter | 26557 |
| Final after housekeeping filter | 25233 |

Per sample:

| Sample | Raw | Final | Retained |
|---|---:|---:|---:|
| PDO | 2868 | 2556 | 89.1% |
| T0 | 3945 | 3519 | 89.2% |
| T1 | 3463 | 3239 | 93.5% |
| T2 | 7011 | 6358 | 90.7% |
| T4 | 2841 | 2600 | 91.5% |
| R4 | 3413 | 2930 | 85.8% |
| eR4 | 4534 | 4031 | 88.9% |

Merged final object:

- genes: `62754`
- cells: `25233`
- clusters: `11`

Final cells per sample in the merged object:

| Sample | Cells |
|---|---:|
| PDO | 2556 |
| T0 | 3519 |
| T1 | 3239 |
| T2 | 6358 |
| T4 | 2600 |
| R4 | 2930 |
| eR4 | 4031 |

## Output files

Main data outputs:

- `parse_outs/Auto_parse_merged.rds`
- `parse_outs/Auto_parse_all_meta.rds`
- `parse_outs/by_samples/<sample>/Auto_<sample>_prefilter.rds`
- `parse_outs/by_samples/<sample>/Auto_<sample>_final.rds`

Summary tables:

- `parse_outs/summary/Auto_parse_sample_manifest.csv`
- `parse_outs/summary/Auto_parse_doublet_summary.csv`
- `parse_outs/summary/Auto_parse_filtering_summary.csv`
- `parse_outs/summary/Auto_parse_filtered_sample_summary.csv`
- `parse_outs/summary/Auto_parse_thresholds.csv`
- `parse_outs/summary/Auto_parse_qc_heatmap_stage_summary.csv`

Plots:

- `parse_outs/plots/Auto_parse_qc_inspection_violin.pdf`
- `parse_outs/plots/Auto_parse_doublet_filtering.pdf`
- `parse_outs/plots/Auto_parse_cells_filtering.pdf`
- `parse_outs/plots/Auto_parse_merged_umap.png`
- `parse_outs/plots/Auto_parse_marker_featureplots.png`
- `parse_outs/plots/Auto_parse_QC_prefilter.png`
- `parse_outs/plots/Auto_parse_QC_final.png`
- `parse_outs/plots/Auto_parse_QC_heatmaps.pdf`
- `parse_outs/summary/Auto_global_filter_bar.png`
- `parse_outs/summary/Auto_per_sample_before_after.png`
- `parse_outs/summary/Auto_per_sample_filtering_steps.png`

## Implementation note

The heavy QC and doublet run completed successfully in PBS. A late-stage plotting incompatibility affected only the merged UMAP metadata plot and the optional marker feature plot. That was corrected by patching the plotting code and regenerating the missing end-stage outputs from the saved merged Seurat object without rerunning the expensive QC stages.
