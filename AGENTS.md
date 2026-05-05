# AGENTS.md — Parse_Pipeline (Parse Biosciences scRNA-seq Pipeline)

Parse Biosciences single-cell RNA-seq QC and analysis pipeline for OAC (Oesophageal Adenocarcinoma) treatment-response samples. Runs on Imperial College HPC (PBS Pro scheduler). Computation uses **R** with **bash** PBS wrappers, plus **Python** for upstream Parse demultiplexing and RNA velocity.

## Repository Structure

```
data/                — Raw FASTQ deliveries and filtered count matrices (symlinked)
tools/               — External pipelines and utilities (Parse splitpipe, velocity, trailmaker)
analysis/            — Downstream analysis R scripts organized by topic
  cell_states/       — State assignment, abundance, pseudotime
  cnv/               — InferCNA copy number analysis
  metaprograms/      — GeneNMF, enrichment, MP correlation, high-res trends
  plotting/          — QC heatmaps, plotting utilities
  summary/           — Cross-sample summary statistics
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
| 2 | via metaprogram submission scripts | `analysis/metaprograms/Auto_parse_geneNMF.R` | All malignant epi | gnmf |
| 3 | via metaprogram submission scripts | `analysis/metaprograms/Auto_parse_find_optimal_nmf.R` | All | dmtcp |
| 4 | via metaprogram submission scripts | `analysis/metaprograms/Auto_parse_enrichment_annotation.R` | All | dmtcp |
| 5 | `analysis/cnv/Auto_parse_infercna.sh` | `analysis/cnv/Auto_parse_infercna.R` | All | dmtcp |
| 6 | `analysis/cell_states/Auto_parse_mp_abundance_activity_dual.sh` | `analysis/cell_states/Auto_parse_mp_abundance_activity_dual.R` | All | dmtcp |

## Build / Run / Test Commands

There is no build system, linter, or test suite. All execution is via PBS `qsub`.

```bash
# Step 1 — QC pipeline
qsub Auto_parse_QC_Pipeline.sh

# GeneNMF metaprogram analysis
qsub analysis/metaprograms/Auto_parse_geneNMF.sh

# Optimal nMP selection
qsub analysis/metaprograms/Auto_parse_find_optimal_nmf.sh

# Enrichment annotation
qsub analysis/metaprograms/Auto_parse_enrichment_annotation.sh

# MP correlation
qsub analysis/metaprograms/Auto_run_mp_correlation.sh

# Cross-data MP correlation (external datasets)
qsub analysis/metaprograms/Auto_parse_mp_correlation_with_pancancer.sh

# InferCNA
qsub analysis/cnv/Auto_parse_infercna.sh

# MP abundance/activity dual report
qsub analysis/cell_states/Auto_parse_mp_abundance_activity_dual.sh

# Submit full metaprogram workflow (geneNMF → optimal → enrichment → correlation)
bash analysis/metaprograms/Auto_parse_submit_metaprogram_workflow.sh

# Run R interactively
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/cell_states/Auto_parse_mp_abundance_activity_dual.R

# For GeneNMF / UCell scripts, use the gnmf environment
source activate /rds/general/user/sg3723/home/anaconda3/envs/gnmf
```

## HPC & File Safety Rules

These rules are **mandatory** for any agent operating in this repo:

1. **Working directory**: All outputs go to `parse_outs/`. Never write outside project paths.
2. **Conda init**: Always run `eval "$(~/miniforge3/bin/conda shell.bash hook)"` before activating envs.
3. **Interactive first**: Tasks under 8 cores / 64 GB → write only the `.R` script, no `.sh` wrapper. User runs interactively.
4. **PBS required**: Heavy tasks → must create PBS `.sh` script with `#PBS` resource headers.
5. **Live Logging**: Always use live streaming log file mode by adding `#PBS -koed` to the submission script.
6. **File naming**: New persistent files MUST be prefixed with `Auto_` (e.g., `Auto_analysis.R`).
7. **Modifying existing files**: New code MUST be wrapped in 20-hash comment blocks:
   ```r
   ####################
   # your new code here
   ####################
   ```
8. **No deleting/modifying** existing lines outside 20-hash blocks without permission.
9. **Test scripts**: Name `delete_<desc>.R` and delete immediately after use.
10. **Max concurrent PBS jobs**: 46 (throttled via `while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]`).

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
WD=/rds/general/project/spatialtranscriptomics/ephemeral/Parse_Pipeline
cd $WD
Rscript <script>.R
echo $(date +%T)
```

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

**Plotting**: `ggplot2` with `theme_minimal()` or `theme_classic()`. Plots saved via `ggsave()` or `pdf()`/`png()` + `dev.off()`. Composite layouts with `patchwork`.

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
| `Auto_parse_geneNMF.R` | Run GeneNMF on merged malignant epithelial cells |
| `Auto_parse_find_optimal_nmf.R` | Determine optimal nMP using silhouette + WSS |
| `Auto_parse_enrichment_annotation.R` | Multi-database enrichment annotation of Parse MPs |
| `Auto_parse_mp_correlation.R` | Within-cohort MP Spearman correlation |
| `Auto_parse_mp_correlation_external.R` | Cross-dataset MP correlation (scATLAS, PDO) |
| `Auto_parse_highres_mp_trend_filter.R` | High-resolution temporal MP trend filtering |

### `analysis/cell_states/` — Cell State Analysis
| File | Purpose |
| :--- | :--- |
| `Auto_parse_mp_abundance_activity_dual.R` | MP abundance and activity dual reporting (scATLAS + PDO + Parse MPs) |

### `analysis/trajectory/` — Trajectory and Pseudotime Analysis
| File | Purpose |
| :--- | :--- |
| `Auto_parse_pseudotime_samples.R` | Per-sample Monocle3 pseudotime |
| `Auto_parse_pseudotime_linear_plot.R` | Linear pseudotime visualisation |
| `Auto_parse_pseudotime_sample_distance_matrix.R` | Pseudotime-based state distance matrix |
| `Auto_prepare_velocity_inputs.py` | Extract barcodes and coordinates for velocity |
| `Auto_velocyto_parse_run.py` | Run velocyto on Parse BAMs |
| `Auto_scvelo_visualise.py` | scVelo analysis and visualization |

### `analysis/cnv/` — Copy Number Variation
| File | Purpose |
| :--- | :--- |
| `Auto_parse_infercna.R` | InferCNA on Parse epithelial cells |

### `analysis/plotting/` — QC and Plotting
| File | Purpose |
| :--- | :--- |
| `Auto_parse_qc_heatmap.R` | QC metric heatmap |

### `analysis/summary/` — Summary Statistics
| File | Purpose |
| :--- | :--- |
| `Auto_parse_summary.R` | Cross-sample summary statistics |

## External Data Dependencies

| Resource | Path | Description |
| :--- | :--- | :--- |
| scATLAS epithelial | `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/EAC_Ref_epi.rds` | Reference epithelial Seurat for cross-data comparison |
| scATLAS MPs | `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds` | scATLAS metaprogram definitions |
| scATLAS UCell scores | `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/UCell_nMP19_filtered.rds` | Filtered UCell scores |
| 3CA gene sets | `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv` | Pan-cancer 3CA metaprograms |
| Cell cycle genes | `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv` | Cell cycle gene annotations |
| Gene order | `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt` | Gene position file for InferCNA |

## Sample Information

Parse Biosciences combinatorial barcoding platform, 2 replicates (A and B) per sample. Current samples are NACT1090 treatment-response time-course (T0, T1, T2, T4, R4, eR4).

## Velocity Analysis

RNA velocity pipeline uses Python (scVelo + velocyto) in the `bidcell_temp` conda env:
```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/bidcell_temp
```

Key scripts in `analysis/trajectory/`:
- `Auto_prepare_velocity_inputs.py` — extract barcodes and coordinates from Seurat
- `Auto_velocyto_parse_run.py` — run velocyto on Parse BAMs
- `Auto_scvelo_visualise.py` — scVelo analysis and visualization

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
- Create a new `Auto_` script.

Append new findings to the appropriate section. Don't rewrite existing documentation unless fixing an error.

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
