# AGENTS.md — Parse_Pipeline (Parse Biosciences scRNA-seq Pipeline)

Parse Biosciences single-cell RNA-seq QC and analysis pipeline for OAC (Oesophageal Adenocarcinoma) treatment-response samples. Runs on Imperial College HPC (PBS Pro scheduler). Computation uses **R** with **bash** PBS wrappers, plus **Python** for upstream Parse demultiplexing and RNA velocity.

## Repository Structure

```
data/                — Raw FASTQ deliveries and filtered count matrices (symlinked)
tools/               — External pipelines and utilities (Parse splitpipe, velocity, trailmaker)
analysis/            — Downstream analysis R scripts organized by topic
  common/            — Shared config, constants, helpers, logging
  input_preparation/ — Input conversion utilities for extra datasets
  cell_states/       — State assignment, abundance, activity reports
  cnv/               — InferCNA copy number analysis
  metaprograms/      — GeneNMF, enrichment, MP correlation, high-res trends
  methodology/       — Script-by-script methodology files
  plotting/          — QC heatmaps, plotting utilities
  summary/           — Cross-sample summary statistics
  trajectory/        — Pseudotime and RNA velocity workflows
parse_outs/          — All pipeline outputs
temp/                — PBS stdout/stderr logs
AGENTS.md            — Workspace rules and execution guidance
```

## Data Sources

| Source | Path | Description |
| :--- | :--- | :--- |
| Raw FASTQ batch 1 | `data/X204SC26033053-Z01-F001_01/` | Novogene sequencing delivery batch 1 |
| Raw FASTQ batch 2 | `data/X204SC26033053-Z01-F001_02/` | Novogene sequencing delivery batch 2 |
| Filtered matrices | `data/8da03a32-a9ad-4938-9338-cbd1731f445d_filtered_matrices.zip` | Parse-demultiplexed filtered count matrices |
| Trailmaker outputs | `data/trailmaker/` | Parse Trailmaker combinatorial-indexed outputs for 2 replicates |

## Tools

| Tool | Path | Purpose |
| :--- | :--- | :--- |
| Parse splitpipe | `tools/ParseBiosciences-Pipeline.1.7.2/` | Parse Biosciences official demultiplexing pipeline v1.7.2 |
| Parse scRNA scripts | `tools/parse_scRNA_pipeline/` | Custom Parse alignment/demux scripts (STARsolo-based) |
| RNA velocity | `tools/Auto_velocity/` | scVelo + velocyto analysis pipeline |

## Pipeline Execution Order

| Step | Shell Script | R Script | Scope | Conda Env |
|------|-------------|----------|-------|-----------| 
| 1 | `Auto_parse_QC_Pipeline.sh` | `Auto_parse_QC_Pipeline.R` | All samples | dmtcp |
| 2 | `Auto_parse_geneNMF.sh` | `analysis/metaprograms/parse_metaprogram_geneNMF_discovery.R` | Six Parse samples | gnmf |
| 3 | `Auto_parse_find_optimal_nmf.sh` | `analysis/metaprograms/parse_metaprogram_select_optimal_nMP.R` | Parse MP sweep | gnmf |
| 4 | `Auto_parse_geneNMF.sh` rerun | `analysis/metaprograms/parse_metaprogram_geneNMF_discovery.R` | Selected nMP UCell scoring | gnmf |
| 5 | `Auto_parse_enrichment_annotation.sh` | `analysis/metaprograms/parse_metaprogram_enrichment_annotation.R` | Selected Parse MPs | dmtcp |
| 6 | `Auto_run_mp_correlation.sh` | `analysis/metaprograms/parse_metaprogram_internal_correlation.R` | Selected Parse MPs | dmtcp |
| 7 | `Auto_parse_mp_abundance_activity_dual.sh` | `analysis/cell_states/parse_mp_state_abundance_activity_approachB_noreg.R` | All 9 samples | dmtcp |
| 8 | `Auto_parse_mp_correlation_with_pancancer.sh` | `analysis/metaprograms/parse_metaprogram_external_reference_correlation.R` | Parse vs scATLAS/PDO | dmtcp |
| 9 | `Auto_parse_infercna.sh` | `analysis/cnv/parse_infercna_carroll_reference.R` | All 9 samples | dmtcp |
| 10 | `Auto_parse_pseudotime_samples.sh` | `analysis/trajectory/parse_pseudotime_sample_trajectory.R` then plotting/distance scripts | Six Parse samples | dmtcp |

## Build / Run / Test Commands

There is no build system, linter, or test suite. All execution is via PBS `qsub`.

```bash
# Step 1 — QC pipeline
qsub Auto_parse_QC_Pipeline.sh

# GeneNMF metaprogram analysis
qsub Auto_parse_geneNMF.sh

# Optimal nMP selection
qsub Auto_parse_find_optimal_nmf.sh

# Rerun GeneNMF after optimal nMP selection to score selected MPs
qsub Auto_parse_geneNMF.sh

# Enrichment annotation
qsub Auto_parse_enrichment_annotation.sh

# MP correlation
qsub Auto_run_mp_correlation.sh

# Cross-data MP correlation (external datasets)
qsub Auto_parse_mp_correlation_with_pancancer.sh

# InferCNA
qsub Auto_parse_infercna.sh

# Preferred Approach B/noreg MP state abundance/activity report
qsub Auto_parse_mp_abundance_activity_dual.sh

# Submit full metaprogram workflow (geneNMF → optimal → enrichment → correlation)
bash Auto_parse_submit_metaprogram_workflow.sh

# Run R interactively
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/cell_states/parse_mp_state_abundance_activity_approachB_noreg.R

# For GeneNMF / UCell scripts, use the gnmf environment
source activate /rds/general/user/sg3723/home/anaconda3/envs/gnmf
```

## HPC & File Safety Rules

These rules are **mandatory** for any agent operating in this repo:

1. **Working directory**: All outputs go to `parse_outs/`. Never write outside project paths.
2. **Conda init**: Always run `eval "$(~/miniforge3/bin/conda shell.bash hook)"` before activating envs.
3. **No Heavy Workloads on Login Nodes (MANDATORY)**: Strictly prohibit running any computationally, memory, or IO intensive workloads on the login nodes, as it adversely affects other users. Any even slightly larger workloads MUST be submitted to the batch queue via PBS `qsub`.
4. **PBS required**: All analytical and heavy tasks → must create PBS `.sh` script with `#PBS` resource headers and submit to the queue. Interactive scripts on login nodes are only permitted for very light, trivial tasks.
5. **Live Logging**: Always use live streaming log file mode by adding `#PBS -koed` to the submission script.
6. **File naming**: Active analysis scripts should use informative names without the old `Auto_` script prefix. Use `legacy_` for comparison-only scripts and `delete_` for scripts the user should manually remove. Existing output filenames may keep `Auto_parse_*` for downstream compatibility.
7. **Modifying existing files**: New code MUST be wrapped in 20-hash comment blocks:
   ```r
   ####################
   # your new code here
   ####################
   ```
8. **No deleting/modifying** existing lines outside 20-hash blocks without permission.
9. **Test scripts**: Name `delete_<desc>.R` and delete immediately after use.
10. **Max concurrent PBS jobs**: 46 (throttled via `while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]`).
11. **Script headers**: Every persistent script must start with a 20-hash header documenting description, exact inputs, exact outputs, cache/replot behavior for expensive workflows, methodology file path, and downstream status.
12. **Methodology docs**: Every active or legacy script must have a methodology file under `analysis/methodology/<matching_subfolder>/`. Update `analysis/script_inventory_and_dependency_map.md` when scripts move, change dependencies, or become legacy/delete candidates.
13. **Output tiers**: New or substantially revised workflows should organize outputs into `intermediate/`, `tables/`, `figures/`, `logs/`, and `reports/` within their output directory. Existing legacy output locations can stay stable but must be classified in the header/methodology.
14. **Run summaries**: Long-running R scripts should use `analysis/common/parse_pipeline_logging.R` to write `parse_outs/logs/run_summaries/<script>_<timestamp>.txt` with start/end time, inputs, outputs, parameters, cache reuse, and session info.
15. **Plot readability**: Slide-facing plots must use readable font sizes, legend text, row/column labels, and point sizes. Prefer wide PDFs/PNGs and avoid tiny labels that cannot be read in PowerPoint.
16. **Storage policy — live vs ephemeral**: All scripts, final outputs (RDS data objects, figures, tables, logs, reports), and **all critical inputs required for replotting** must be read from and written to the `live` project path (`/rds/general/project/spatialtranscriptomics/live/Parse_Pipeline/`). **Exception**: exceptionally large intermediate/cache files (typically under `intermediate/` output tiers) should continue to be stored under the corresponding `ephemeral` path (`/rds/general/project/spatialtranscriptomics/ephemeral/Parse_Pipeline/`). **CRITICAL**: The `live` storage must be completely self-sufficient for final presentations; even if the `ephemeral` directory is completely deleted, you must still be able to easily reproduce all plots and critical information using only the files saved in `live/`. Scripts must `dir.create(..., recursive = TRUE, showWarnings = FALSE)` for ephemeral intermediate paths if they do not exist.

    **Decision checklist — before writing ANY `.rds`, `.csv`, or `.pdf` with `saveRDS` / `write.csv` / `ggsave`:**
    1. Is this file read by a *different* downstream script (not just the script that created it)?  → **Must be in `live`**.
    2. Is this a per-cell score matrix, state vector, gene list, or enrichment result that downstream plotting or annotation scripts depend on?  → **Must be in `live`**.
    3. Is this file needed to regenerate any figure or table without re-running the producing script?  → **Must be in `live`**.
    4. Is this file *only* a cache to speed up re-runs of the *same* script and can be fully regenerated from inputs already in `live`?  → Ephemeral is acceptable, but a live copy is still preferred.

    **When in doubt, save to BOTH `live` and `ephemeral`** (ephemeral as cache, live as the persistent copy). Use the pattern:
    ```r
    saveRDS(obj, file.path(outdir_ephemeral, "filename.rds"))  # cache
    saveRDS(obj, file.path(outdir_live, "filename.rds"))        # persistent
    ```
17. **Strictly forbid fallbacks**: Must use only the first and best option for any analysis. If not available or if it fails, the script should stop immediately.

### PBS Job Template
```bash
#!/bin/bash
#PBS -l select=1:ncpus=<N>:mem=<M>gb
#PBS -l walltime=<HH:MM:SS>
#PBS -N <jobname>
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/spatialtranscriptomics/live/Parse_Pipeline
cd $WD
Rscript <script>.R
echo $(date +%T)
```

## Shared Configuration and Script Map

| File | Purpose |
| :--- | :--- |
| `analysis/script_inventory_and_dependency_map.md` | Required analysis-folder map with run order, terminal figures, legacy/delete candidates, output tiers, and stale pointer checks |
| `analysis/common/parse_pipeline_config.R` | R project paths, sample/state constants, thresholds, reference paths, output tiers, plot defaults |
| `analysis/common/parse_pipeline_helpers.R` | Shared R helpers for package loading, count extraction, chunking, MP gene tables, z-normalization, and slide themes |
| `analysis/common/parse_pipeline_logging.R` | R run-summary logging |
| `analysis/common/parse_pipeline_config.py` | Python path/sample constants for velocity workflows |

Preferred current state-definition constants live in `parse_state_definition`:

- method: `Approach B`
- normalization: `noreg`
- unresolved threshold: `0.5`
- hybrid gap: `0.3`
- preferred assignment object: `parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds`

## Code Style Guidelines

### R Scripts

**Imports**: `library()` calls at top of file, one per line. No `require()`.

**Working directory**: Each script calls `setwd()` near the top, typically pointing to the output directory (e.g., the `parse_outs/` subdirectory under ``). Scripts should write outputs relative to `parse_outs/`.

**Variable naming**: `snake_case` for variables and functions. Short descriptive names:
- `tmdata`, `tmdata_list`, `merged_obj` — Seurat objects
- `mp.genes`, `mp_gene_lists` — metaprogram gene collections
- `geneNMF.programs`, `geneNMF.metaprograms` — GeneNMF outputs

**Data structures**: Heavy use of Seurat objects. Lists of Seurat objects keyed by sample name.

**Tidyverse style**: Pipe-heavy with `%>%` and `|>`. `dplyr` verbs, `tidyr::pivot_longer/wider`.

**Plotting**: `ggplot2` with `theme_minimal()` or `theme_classic()`. Plots saved via `ggsave()` or `pdf()`/`png()` + `dev.off()`. Composite layouts with `patchwork`. Slide-facing outputs must use legible text/legend/point sizes and dimensions appropriate for PowerPoint; use `analysis/common/parse_pipeline_config.R` plotting defaults where possible.

**Output format**: `.rds` files for data, `.png` / `.pdf` for plots, `.csv` for summary tables.

### Shell Scripts

**Shebang**: `#!/bin/bash` — always first line.

**Timestamps**: Every script prints `echo $(date +%T)` at start and end.

**Module loading**: `module purge` then `module load tools/dev` before conda.

### Error Handling

- Guard clauses with `stop()` for fatal conditions.
- `tryCatch()` for subsetting operations that might produce empty results.
- NULL assignment + skip pattern for low-cell objects.

### Key R Packages

**Core**: Seurat, dplyr, tidyr, purrr, ggplot2, Matrix, parallel
**Specialised**: infercna, GeneNMF, UCell, ComplexHeatmap, fgsea, msigdbr, monocle3
**Plotting**: patchwork, gridExtra, ggrepel, RColorBrewer, circlize, scales

### Naming Conventions for Cell Types

Lowercase with dots: `t.cell`, `b.cell`, `nk.cell`, `macrophage`, `fibroblast`, `endothelial`, `epithelial`, `plasma`, `dendritic`, `mast`. Multi-labels joined with `|`. Unknown/ambiguous cells labelled `unresolved`.

## Key Shared Data Objects

All paths are relative to `parse_outs/` unless absolute paths are specified.

| Object | Path | Description |
| :--- | :--- | :--- |
| Merged Seurat | `Auto_parse_merged.rds` | Merged Seurat object with all samples post-QC |
| All-cell metadata | `Auto_parse_all_meta.rds` | Per-cell metadata for all cells |
| Per-sample outputs | `by_samples/<sample>/` | Post-QC per-sample Seurat objects |
| GeneNMF outputs | `Auto_parse_metaprograms/` | Metaprogram results per nMP |
| CNV outputs | `cnv/` | InferCNA results |
| Cell state outputs | `cell_states/` | State assignments and abundance reports |
| High-res MP trends | `Auto_parse_highres_metaprogram_trends/` | High-resolution MP trend filtering outputs |
| Pseudotime outputs | `pseudotime_samples/` | Monocle3 pseudotime per sample |
| QC summary | `summary_outputs/` | Summary CSVs and QC plots |
| QC plots | `plots/` | UMAP and QC heatmap PDFs |

## Analysis Scripts

### `analysis/metaprograms/` — Metaprogram Analysis
| File | Purpose |
| :--- | :--- |
| `parse_metaprogram_geneNMF_discovery.R` | Run GeneNMF sweep and selected-nMP UCell scoring for Parse samples |
| `parse_metaprogram_select_optimal_nMP.R` | Determine optimal nMP using silhouette + WSS |
| `parse_metaprogram_enrichment_annotation.R` | Multi-database enrichment annotation of selected Parse MPs |
| `parse_metaprogram_internal_correlation.R` | Within-cohort MP Spearman/Jaccard diagnostics |
| `parse_metaprogram_external_reference_correlation.R` | Cross-dataset MP correlation and Jaccard overlap against scATLAS/PDO |
| `parse_highres_mp_strict_mean_median_trend_filter.R` | Active high-resolution MP strict mean/median trend filtering |
| `parse_highres_mp_tcga_survival_volcano.R` | TCGA whole-profile GSVA/Cox volcano plots for high-resolution Parse MP groups |
| `legacy_parse_highres_mp_t2t4_comparison_filter.R` | Legacy T2/T4-high high-resolution MP comparison; not downstream |

### `analysis/cell_states/` — Cell State Analysis
| File | Purpose |
| :--- | :--- |
| `parse_mp_state_abundance_activity_approachB_noreg.R` | Preferred Approach B/noreg MP state assignment, abundance, and activity reporting |

### `analysis/trajectory/` — Trajectory and Pseudotime Analysis
| File | Purpose |
| :--- | :--- |
| `parse_pseudotime_helpers.R` | Shared helpers for sample-level Monocle3 workflows |
| `parse_pseudotime_sample_trajectory.R` | Per-sample Monocle3 pseudotime rooted on T0 |
| `parse_pseudotime_linear_sample_report.R` | Linear pseudotime visualization/report |
| `parse_pseudotime_sample_distance_maps.R` | Pseudotime/graph/UMAP sample distance matrices and plots |
| `parse_velocity_prepare_inputs.py` | Extract barcodes, metadata, and references for velocity; downloads RepeatMasker if absent |
| `parse_velocity_run_velocyto_parse_bam.py` | Run velocyto on Parse BAMs |
| `parse_velocity_scvelo_visualise.py` | scVelo analysis and visualization |
| `parse_velocity_submit_pipeline.sh` | Submit full velocity workflow with PBS dependencies |
| `delete_Auto_pseudotime_linear_plot_legacy_wrapper.R` | Delete candidate compatibility wrapper |
| `delete_Auto_pseudotime_state_distance_matrix_legacy_wrapper.R` | Delete candidate compatibility wrapper |
| `delete_parse_velocity_submit_pipeline_run2_duplicate.sh` | Delete candidate duplicate velocity submitter |

### `analysis/cnv/` — Copy Number Variation
| File | Purpose |
| :--- | :--- |
| `parse_infercna_carroll_reference.R` | InferCNA on Parse/PDO/SUR1090 epithelial cells with Carroll_2023 reference |

### `analysis/plotting/` — QC and Plotting
| File | Purpose |
| :--- | :--- |
| `parse_qc_marker_heatmap.R` | Raw/final QC marker heatmaps |

### `analysis/summary/` — Summary Statistics
| File | Purpose |
| :--- | :--- |
| `parse_filtering_summary_plots.R` | Cross-sample QC filtering summary plots |

### `analysis/input_preparation/` — Input Preparation
| File | Purpose |
| :--- | :--- |
| `parse_prepare_sur1090_count_matrices.R` | Convert SUR1090 treated/untreated CSV matrices into Parse-compatible DGE_filtered folders |

### `analysis/publication/` — Publication-Quality Figures
| File | Purpose |
| :--- | :--- |
| `parse_survival_plot.R` | Kaplan-Meier curve for stress MP94 using TCGA ESCA |
| `parse_state_abundance_timepoint.R` | Publication-quality stacked bar of PDO-pipeline state abundance across 6 Parse timepoints (excluding Unresolved/Hybrid) |

## External Data Dependencies

| Resource | Path | Description |
| :--- | :--- | :--- |
| scATLAS epithelial | `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/EAC_Ref_epi.rds` | Reference epithelial Seurat for cross-data comparison |
| scATLAS MPs | `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds` | scATLAS metaprogram definitions |
| scATLAS UCell scores | `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/UCell_nMP19_filtered.rds` | Filtered UCell scores |
| PDO-pipeline MPs | `/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds` | PDO metaprogram definitions for Approach B/noreg state assignment |
| Carroll CNV reference | `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Carroll_2023_reference.rds` | External InferCNA reference |
| 3CA gene sets | `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv` | Pan-cancer 3CA metaprograms |
| Cell cycle genes | `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv` | Cell cycle gene annotations |
| Gene order | `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt` | Gene position file for InferCNA |
| Velocity GTF | `/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz` | RNA velocity gene annotation source |
| Velocity RepeatMasker | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz` | Downloaded by `parse_velocity_prepare_inputs.py` if absent |
| TCGA ESCA metadata and whole-profile TPM | `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/tcga_esca_meta.rds`; `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt` | Reference TCGA dataset used for whole-profile GSVA survival association plots |

## Sample Information

Parse Biosciences combinatorial barcoding platform, 2 replicates (A and B) per sample. Current samples are NACT1090 treatment-response time-course (T0, T1, T2, T4, R4, eR4).

## Velocity Analysis

RNA velocity pipeline uses Python (scVelo + velocyto) in the `bidcell_temp` conda env:
```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/bidcell_temp
```

Key scripts in `analysis/trajectory/`:
- `parse_velocity_prepare_inputs.py` — extract barcodes and coordinates from Seurat
- `parse_velocity_run_velocyto_parse_bam.py` — run velocyto on Parse BAMs
- `parse_velocity_scvelo_visualise.py` — scVelo analysis and visualization

## Directory Migration Notes

This workspace was previously nested under an `Auto_parse_qc_pipeline/` directory. It has now been flattened and fully reorganized:
- All `.sh` and `.pbs` submission scripts reside in the project root.
- All downstream analysis scripts (R and Python) reside in `analysis/` subdirectories.
- All pipeline outputs reside directly in `parse_outs/`.
- Job log files (`.o` and `.e`) are gathered in the project root.

## Critical Recurring Patterns

**MP Silhouette Filtering** (same as all pipelines):
```r
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
```

**Sample identity**: Use `orig.ident` from Seurat metadata.

**Preferred cell-state definition**: PDO-pipeline metaprogram state assignment with **Approach B** and **noreg** adjusted scores. The canonical assignment is `parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds`; the canonical adjusted score matrix is `parse_outs/cell_states/Auto_parse_PDOpipeline_mp_adj_noreg.rds`.

**Environment usage**:
- `gnmf` conda env: UCell scoring and GeneNMF scripts
- `dmtcp` conda env: general Seurat and analysis tasks
- `bidcell_temp` conda env: Python velocity analysis

## Uploading Files to Google Drive (rclone)

rclone is configured with a remote named `gdrive`. Always upload into the `IMPERIAL/` folder:

```bash
module load rclone
rclone copy <local_file> gdrive:IMPERIAL/ --progress
```

## AGENTS.md Living Document Rule

Future agents must update this file when they:
- Find a new analysis script and define its purpose.
- Locate new input or output file paths.
- Spot a recurring pattern or technical hurdle.
- Create or rename an analysis script.
- Add or change a methodology file.
- Change any script dependency, output tier, terminal/legacy status, or run order.

Also update `analysis/script_inventory_and_dependency_map.md` in the same change. Append new findings to the appropriate section. Don't rewrite existing documentation unless fixing an error.

####################
## 2026-05-15 Added Tumour Stiffness Gene Module Heatmap

- Active terminal figure script: `analysis/cell_states/parse_tumour_stiffness_gene_module_heatmap.R`.
- Methodology: `analysis/methodology/cell_states/tumour_stiffness_gene_module_heatmap_methodology.md`.
- Input: `parse_outs/Auto_parse_merged.rds`.
- Outputs: `parse_outs/cell_states/tumour_stiffness_gene_module/` with `intermediate/`, `tables/`, `figures/`, and `reports/` tiers.
- Main terminal figure: `parse_outs/cell_states/tumour_stiffness_gene_module/figures/Auto_parse_tumour_stiffness_module_gene_heatmap.pdf`.
- The heatmap uses sample columns, grouped gene/module rows, a vibrant blue-white-red color transition based on row z-score, horizontal and wrapped row titles for the gene groups, top-aligned columns with the cell-count barplot shifted directly adjacent to the columns, and raw mean log-normalized score labels scaled up to a highly legible 8pt size.
####################

####################
## 2026-05-15 Added T2+T4 versus T0 DEG and Pathway Response Workflow

- Active terminal response script: `analysis/cell_states/parse_t2t4_vs_t0_dge_pathway_response.R`.
- Methodology: `analysis/methodology/cell_states/t2t4_vs_t0_dge_pathway_response_methodology.md`.
- Inputs: `parse_outs/by_samples/<sample>/Auto_<sample>_final.rds`, MSigDB Hallmark gene sets via `msigdbr`, configured cell-cycle genes, and PDO-pipeline MP4/MP8 metaprogram genes.
- Outputs: `parse_outs/t2t4_vs_t0_response/` with `intermediate/`, `tables/`, `figures/`, `logs/`, and `reports/` tiers.
- Main terminal figures: T2+T4 vs T0 volcano plot, Hallmark FGSEA dot plot, top-DEG sample heatmap, and real-value sample pathway/metric heatmap.
- The DEG contrast pools T2 and T4 cells against T0 and intentionally includes no batch covariate under the sequenced-together assumption.
####################

####################
## 2026-06-08 Publication UMAP Samples Update

- Added active terminal figure script `analysis/publication/parse_umap_samples.R`.
- Methodology: `analysis/methodology/publication/umap_samples_methodology.md`.
- Depends on `parse_outs/Auto_parse_merged.rds`.
- Main outputs live under `parse_outs/publication/umap_samples/` with `figures/` tier.
- Terminal figures: `figures/umap_samples.pdf` and `.png` — UMAP visualization of the 6 Parse treatment-response timepoints plus PDO, colored by sample.
- No active downstream consumers.
####################

####################
## Nature-Figure Publication Skill (Selective Application)

A `nature-figure` skill is installed at `/rds/general/user/sg3723/home/nature-skills/nature-figure/`. It enforces Nature-journal visual standards. **Agents must exercise judgment to apply this skill primarily to scripts producing final, sharable results.**

### When to Apply (Clever Selection)

Do **not** apply this to every R script. Focus on scripts that synthesize data across samples or produce "Final" visualizations. Prioritize:
- Final cohort-level summary plots (abundance, survival, clinical associations)
- Cross-dataset comparison figures
- Any `Auto_` script producing figures explicitly intended for manuscript inclusion or presentation slides

### How to Apply (R Backend — PDF Priority)

1. **Figure contract**: Define the claim and evidence hierarchy first.
2. **Typography**: Use 6.5pt Arial (Nature standard) via `theme_nature_contract()`.
3. **Export policy (PDF Priority)**: **PDF is the preferred format.** SVG is not required unless requested. Use `grDevices::cairo_pdf()` to ensure font embedding.
   ```r
   save_pub_pdf <- function(plot, filename, width_mm = 183, height_mm = 120) {
     w <- width_mm / 25.4; h <- height_mm / 25.4
     grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
     if (inherits(plot, "Heatmap") || inherits(plot, "HeatmapList")) {
       ComplexHeatmap::draw(plot, merge_legend = TRUE)
     } else {
       print(plot)
     }
     dev.off()
   }
   ```
4. **Color & IA**: Use restrained palettes and follow the **overview → deviation → relationship** information architecture.

### Reference Files
- `~/nature-skills/nature-figure/SKILL.md` — full skill specification
- `~/nature-skills/nature-figure/references/r-workflow.md` — R-specific patterns

### Exceptions
- **Diagnostic/QC scripts**: Step 1-6 pipeline outputs, internal QC heatmaps, and debugging plots should use standard Seurat/ggplot2 defaults to save time.
- **Development/Test scripts**: `delete_*.R` scripts.
####################

####################
## 2026-05-26 Added High-Resolution MP TCGA Survival Volcano Workflow

- Active terminal figure script: `analysis/metaprograms/parse_highres_mp_tcga_survival_volcano.R`.
- Methodology: `analysis/methodology/metaprograms/highres_mp_tcga_survival_volcano_methodology.md`.
- Inputs: strict high-resolution increasing/decreasing MP gene lists from `parse_outs/Auto_parse_highres_metaprogram_trends/`, legacy T2/T4-high MP genes from `parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/`, and the scRef TCGA ESCA metadata/whole-profile TPM compatibility copies under `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/`.
- Outputs: `parse_outs/highres_mp_tcga_survival/` with `intermediate/`, `tables/`, `figures/`, and `reports/` tiers.
- Main terminal figures: `parse_outs/highres_mp_tcga_survival/figures/Auto_parse_highres_mp_tcga_survival_volcano_whole_tcga.pdf` plus individual strict-increase, strict-decrease, and legacy-T2/T4-high volcano PDFs/PNGs.
- Method: whole-TCGA GSVA reference mode followed by Cox survival models in EAC primary tumours. The unsuffixed group volcanoes are continuous-Cox plots; `_median` and `_q1q4` outputs are reference split-model companions.
####################

####################
## 2026-06-08 Added Publication State Abundance Timepoint Figure

- Active terminal figure script: `analysis/publication/parse_state_abundance_timepoint.R`.
- Methodology: `analysis/methodology/publication/state_abundance_timepoint_methodology.md`.
- Input: `parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds`.
- Outputs: `parse_outs/publication/state_abundance_timepoint/` with `tables/` and `figures/` tiers.
- Main terminal figures: `parse_outs/publication/state_abundance_timepoint/figures/state_abundance_timepoint.pdf` and standalone single-column legend `parse_outs/publication/state_abundance_timepoint/figures/state_abundance_legend.pdf`.
- The stacked bar shows PDO-pipeline state proportions across 6 Parse timepoints (T0, T1, T2, T4, R4, eR4), excluding Unresolved and Hybrid. Uses Nature figure contract style with 6.5pt Arial, cairo_pdf export, original color palette, proportion-only display (no cell counts or secondary axis), and enlarged 9pt bold x-axis labels.
####################

####################
## 2026-06-08 Added T2T4 versus T0eR4 High-Res MP Cluster Heatmap

- Active terminal figure script: `analysis/cell_states/parse_t2t4_vs_t0er4_highres_cluster_heatmap.R`.
- Methodology: `analysis/methodology/cell_states/t2t4_vs_t0er4_highres_cluster_heatmap_methodology.md`.
- Inputs: `parse_outs/by_samples/<sample>/Auto_<sample>_final.rds`, global metadata `Auto_parse_all_meta.rds`, and PDO `nMP156` metaprogram genes.
- Outputs: `parse_outs/cell_states/t2t4_vs_t0er4_highres_clusters/` with `intermediate/`, `tables/`, and `figures/` tiers.
- Main terminal figures: `parse_t2t4_vs_t0er4_cluster_delta_heatmap.pdf`, `parse_t2t4_vs_t0er4_cluster_absolute_heatmap.pdf`, `parse_t2t4_vs_t0er4_MP_delta_heatmap.pdf`.
- Scored using UCell, grouping T2+T4 against T0+eR4 to evaluate response and baseline difference for PDO-derived MP functional clusters.
####################

####################
## 2026-06-08 Added 2x2 Grid State Abundance Plot

- Active terminal figure script: `analysis/publication/parse_state_abundance_grid.R`.
- Methodology: `analysis/methodology/publication/state_abundance_grid_methodology.md`.
- Input: `parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds`.
- Outputs: `parse_outs/publication/state_abundance_grid/` with `figures/` and `tables/` tiers.
- Terminal figures: `figures/state_abundance_grid.pdf` and `.png` — 2x2 grid of bar charts showing PDO-pipeline cell-state proportions across the 6 Parse timepoints, styled with specific state colors.
- No active downstream consumers.
####################

####################
## 2026-06-09 Added Publication Trend Plots for MP13 and MP28

- Active terminal figure script: `analysis/publication/parse_publication_mp13_mp28_trends.R`.
- Methodology: `analysis/methodology/publication/mp13_mp28_trends_methodology.md`.
- Inputs: `Auto_parse_highres_T2T4_sample_ucell_summary_nMP117.csv` and `Auto_parse_highres_sample_ucell_summary_nMP117.csv`.
- Outputs: `parse_outs/publication/mp_trends/` with `figures/` tier.
- Terminal figures: `figures/parse_mp13_trend.pdf` and `figures/parse_mp28_trend.pdf` (with `.png` versions) — Nature-style trend plots of mean and median UCell scores across 6 Parse timepoints.
- No active downstream consumers.
####################

####################
## 2026-06-25 Added Timepoint SCENIC Regulon Workflow

- Active terminal workflow script: `analysis/cell_states/parse_timepoint_scenic_regulons.R`.
- PBS wrapper: `parse_timepoint_scenic.sh`.
- Methodology: `analysis/methodology/cell_states/timepoint_scenic_regulons_methodology.md`.
- Input: `parse_outs/Auto_parse_merged.rds`.
- Database confirmed from the PDO SCENIC workflow: `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/cistarget_databases_rcistarget_mc9nr/`, using the hg38 refseq-r80 mc9nr RcisTarget feather files. The sibling `cistarget_databases/` directory contains v10-clust feather files and is not the PDO workflow default.
- Outputs: `parse_outs/cell_states/timepoint_scenic/` with `intermediate/`, `tables/`, `figures/`, `logs/`, and `reports/` tiers.
- The workflow intentionally ignores cell-state and MP concepts. It treats `T0`, `T1`, `T2`, `T4`, `R4`, and `eR4` as entities in one combined six-timepoint SCENIC analysis for comparable regulon AUC/RSS/differential activity. Independent per-timepoint SCENIC analyses are disabled because they create non-comparable regulon target dictionaries.
- Heatmaps use strict timepoint column order `T0`, `T1`, `T2`, `T4`, `R4`, `eR4`. The default heatmap selects the union of top RSS-specific regulons per timepoint; `combined_timepoint_balanced2600_top20_specific_regulon_heatmap.pdf` selects the top 20 balanced-RSS regulons per timepoint, and `combined_timepoint_balanced2600_top20_gap_regulon_heatmap.pdf` selects the top regulons by positive RSS gap versus the next-highest timepoint.
- GENIE3 is run with `resumePreviousRun = TRUE` and `genie3_nparts = 100` by default. Existing `int/1.3_GENIE3_weightMatrix_part_*.Rds` files are reused after walltime interruption, and continuation jobs resume remaining target genes instead of restarting.
- No active downstream consumers.
####################

####################
## 2026-06-30 Added Centred GeneNMF Method-Comparison Workflow

- Active comparison scripts now live under `analysis/metaprograms/centred/`: `parse_centred_metaprogram_geneNMF_discovery.R`, `parse_centred_highres_mp_strict_mean_median_trend_filter.R`, `parse_centred_highres_mp_t2t4_comparison_filter.R`, `parse_centred_t2t4_vs_t0er4_highres_cluster_heatmap.R`, `parse_compare_centred_vs_uncentred_highres_mps.R`, and `parse_centred_publication_highres_figures.R`.
- PBS wrapper: `parse_centred_metaprogram_workflow.sh`.
- Methodology documentation is kept in the original method files because the centred steps are method-identical except for upstream `multiNMF(center = TRUE)`: `metaprogram_geneNMF_discovery_methodology.md`, `highres_mp_strict_mean_median_trend_filter_methodology.md`, `legacy_highres_mp_t2t4_comparison_filter_methodology.md`, and `publication/highres_metaprogram_heatmap_methodology.md`.
- Inputs: final per-sample Parse Seurat objects, centred GeneNMF caches under `parse_outs/centred/Auto_parse_metaprograms/`, canonical state assignments, and existing uncentred high-resolution outputs for comparison.
- Outputs: all centred method outputs live under `parse_outs/centred/`, including `Auto_parse_highres_metaprogram_trends/`, `Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/`, `comparison/`, and `publication/` figure tiers.
- The centred workflow runs `multiNMF(center = TRUE)`, derives high-resolution MPs at half the total NMF programme count from `T0`, `T1`, `T2`, `T4`, `R4`, and `eR4` (currently 234 / 2 = nMP117), runs the T2/T4-high filter from `parse_highres_mp_t2t4_comparison_filter.R`, annotates centred MPs with automatic top 3CA non-cell-cycle labels and best uncentred T2/T4-high MP gene-set matches, and compares centred versus uncentred retained MPs. It does not use the external PDO nMP156 object.
####################
