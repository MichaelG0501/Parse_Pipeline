####################
# parse_metaprogram_external_reference_correlation.R
#
# Description:
#   Compares selected Parse metaprograms to scATLAS and PDO-pipeline
#   metaprograms using gene-set Jaccard overlap and cross-cell UCell expression
#   correlation.
#
# Inputs:
#   parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds
#   parse_outs/Auto_parse_metaprograms/Auto_parse_UCell_scores_filtered_nMP<k>.rds
#   parse_outs/cell_states/Auto_parse_scATLAS_UCell_scores.rds
#   parse_outs/cell_states/Auto_parse_PDOpipeline_UCell_scores.rds
#   scATLAS and PDO-pipeline GeneNMF reference RDS files
#
# Outputs:
#   parse_outs/Auto_parse_metaprograms/Auto_parse_mp_correlation/Auto_parse_mp_jaccard_vs_*.pdf
#   parse_outs/Auto_parse_metaprograms/Auto_parse_mp_correlation/Auto_parse_mp_correlation_vs_*.pdf
#   parse_outs/logs/run_summaries/parse_metaprogram_external_reference_correlation_*.txt
#
# Cache / replot:
#   Reads cached UCell score matrices from metaprogram and cell-state outputs.
#   Plot styling can be regenerated without rescoring signatures.
#
# Methodology:
#   analysis/methodology/metaprograms/metaprogram_external_reference_correlation_methodology.md
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_metaprogram_external_reference_correlation",
  input_files = c(
    "parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds",
    "parse_outs/Auto_parse_metaprograms/Auto_parse_UCell_scores_filtered_nMP<k>.rds",
    "parse_outs/cell_states/Auto_parse_scATLAS_UCell_scores.rds",
    "parse_outs/cell_states/Auto_parse_PDOpipeline_UCell_scores.rds",
    parse_reference_paths$scatlas_metaprograms,
    parse_reference_paths$pdo_metaprograms
  ),
  output_files = c(
    "parse_outs/Auto_parse_metaprograms/Auto_parse_mp_correlation/Auto_parse_mp_jaccard_vs_scATLAS.pdf",
    "parse_outs/Auto_parse_metaprograms/Auto_parse_mp_correlation/Auto_parse_mp_correlation_vs_PDO.pdf"
  )
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(ComplexHeatmap)
  library(circlize)
  library(pheatmap)
  library(dplyr)
  library(grid)
})

out_dir_root <- file.path(parse_project_root(), "parse_outs")
setwd(out_dir_root)

parse_mp_dir <- "Auto_parse_metaprograms"
out_dir <- file.path(parse_mp_dir, "Auto_parse_mp_correlation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cell_states_dir <- "cell_states"

# 1. Parse MPs
geneNMF.metaprograms <- readRDS(file.path(parse_mp_dir, "Auto_parse_MP_outs_default.rds"))
optimal_txt <- file.path(parse_mp_dir, "Auto_parse_optimal_nMP.txt")
optimal_nMP <- if (file.exists(optimal_txt)) readLines(optimal_txt, warn = FALSE)[1] else length(geneNMF.metaprograms$metaprograms.genes)

parse_mp_genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) parse_mp_genes <- parse_mp_genes[!names(parse_mp_genes) %in% paste0("MP", bad_mps)]
coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
parse_mp_genes <- parse_mp_genes[!names(parse_mp_genes) %in% names(coverage_tbl)[coverage_tbl < 0.25]]

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", names(parse_mp_genes)))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
parse_canon_mps <- paste0("MP", mp_tree_order)
parse_mp_genes <- parse_mp_genes[parse_canon_mps]

parse_ucell_path <- file.path(parse_mp_dir, paste0("Auto_parse_UCell_scores_filtered_nMP", optimal_nMP, ".rds"))
parse_ucell <- readRDS(parse_ucell_path)
parse_ucell <- parse_ucell[, parse_canon_mps[parse_canon_mps %in% colnames(parse_ucell)], drop = FALSE]

# 2. scATLAS MPs
scatlas_geneNMF <- readRDS(parse_reference_paths$scatlas_metaprograms)
scatlas_retained <- rownames(scatlas_geneNMF$metaprograms.metrics)[scatlas_geneNMF$metaprograms.metrics$silhouette >= 0]
scatlas_tree_order <- scatlas_geneNMF$programs.clusters[scatlas_geneNMF$programs.tree$order]
scatlas_tree_order <- paste0("MP", unique(scatlas_tree_order))
scatlas_tree_order <- scatlas_tree_order[scatlas_tree_order %in% scatlas_retained]

ordered_block <- function(mps, tree_order_names) {
  out <- intersect(tree_order_names, mps)
  c(out, setdiff(mps, out))
}
scatlas_cc <- intersect(c("MP1", "MP7", "MP9"), scatlas_retained)
scatlas_mp_order <- ordered_block(scatlas_cc, scatlas_tree_order)
scatlas_state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "Immune Infiltrating" = c("MP15")
)
for (state in names(scatlas_state_groups)) {
  scatlas_mp_order <- c(scatlas_mp_order, ordered_block(intersect(scatlas_state_groups[[state]], scatlas_retained), scatlas_tree_order))
}
scatlas_mp_order <- unique(c(scatlas_mp_order, setdiff(scatlas_retained, scatlas_mp_order)))

scatlas_desc <- c("MP1"="G2M Cell Cycle","MP9"="G1S Cell Cycle","MP2"="MYC-related Proliferation","MP17"="Basal-like Transition","MP14"="Hypoxia Adapted Epi.","MP5"="Epithelial IFN Resp.","MP10"="Columnar Diff.","MP8"="Intestinal Diff.","MP13"="Hypoxic Inflam. Epi.","MP7"="DNA Damage Repair","MP18"="Secretory Diff. (Intest.)","MP16"="Secretory Diff. (Gastric)","MP15"="Immune Infiltration","MP12"="Neuro-responsive Epi")
scatlas_mp_genes <- scatlas_geneNMF$metaprograms.genes[scatlas_mp_order]
names(scatlas_mp_genes) <- paste0(names(scatlas_mp_genes), "_", scatlas_desc[names(scatlas_mp_genes)])
names(scatlas_mp_genes)[grepl("_NA$", names(scatlas_mp_genes))] <- sub("_NA$", "", names(scatlas_mp_genes)[grepl("_NA$", names(scatlas_mp_genes))])

scatlas_ucell <- readRDS(file.path(cell_states_dir, "Auto_parse_scATLAS_UCell_scores.rds"))
scatlas_ucell <- scatlas_ucell[, scatlas_mp_order[scatlas_mp_order %in% colnames(scatlas_ucell)], drop = FALSE]
colnames(scatlas_ucell) <- paste0(colnames(scatlas_ucell), "_", scatlas_desc[colnames(scatlas_ucell)])
colnames(scatlas_ucell)[grepl("_NA$", colnames(scatlas_ucell))] <- sub("_NA$", "", colnames(scatlas_ucell)[grepl("_NA$", colnames(scatlas_ucell))])
rownames(scatlas_ucell) <- sub("^(.*?)__(.*)$", "\\1_\\2", rownames(scatlas_ucell))

# 3. PDO MPs
pdo_geneNMF <- readRDS(parse_reference_paths$pdo_metaprograms)
pdo_mp_genes_all <- pdo_geneNMF$metaprograms.genes
bad_mps <- which(pdo_geneNMF$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
coverage_tbl <- pdo_geneNMF$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
pdo_retained <- names(pdo_mp_genes_all)[!names(pdo_mp_genes_all) %in% c(bad_mp_names, low_coverage_mps)]

pdo_tree_order <- pdo_geneNMF$programs.clusters[pdo_geneNMF$programs.tree$order]
pdo_tree_order <- paste0("MP", rev(unique(pdo_tree_order)))
pdo_tree_order <- pdo_tree_order[pdo_tree_order %in% pdo_retained]
pdo_cc <- intersect(c("MP6", "MP7", "MP1", "MP3"), pdo_retained)

pdo_state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "Stress-adaptive" = c("MP10", "MP9"),
  "SMG-like Metaplasia" = c("MP8")
)
pdo_mp_order <- ordered_block(pdo_cc, pdo_tree_order)
for (state in names(pdo_state_groups)) {
  pdo_mp_order <- c(pdo_mp_order, ordered_block(intersect(pdo_state_groups[[state]], pdo_retained), pdo_tree_order))
}
pdo_mp_order <- unique(c(pdo_mp_order, setdiff(pdo_retained, pdo_mp_order)))

pdo_desc_full <- c("MP6"="MP6_G2M Cell Cycle","MP7"="MP7_DNA repair","MP5"="MP5_MYC-related Proliferation","MP1"="MP1_G2M checkpoint","MP3"="MP3_G1S Cell Cycle","MP8"="MP8_Columnar Progenitor","MP10"="MP10_Inflammatory Stress Epi.","MP9"="MP9_ECM Remodeling Epi.","MP4"="MP4_Intestinal Metaplasia")
pdo_mp_genes <- pdo_mp_genes_all[pdo_mp_order]
names(pdo_mp_genes) <- unname(pdo_desc_full[names(pdo_mp_genes)])
names(pdo_mp_genes)[is.na(names(pdo_mp_genes))] <- pdo_mp_order[is.na(names(pdo_mp_genes))]

pdo_ucell <- readRDS(file.path(cell_states_dir, "Auto_parse_PDOpipeline_UCell_scores.rds"))
pdo_ucell <- pdo_ucell[, pdo_mp_order[pdo_mp_order %in% colnames(pdo_ucell)], drop = FALSE]
colnames(pdo_ucell) <- unname(pdo_desc_full[colnames(pdo_ucell)])
colnames(pdo_ucell)[is.na(colnames(pdo_ucell))] <- pdo_mp_order[is.na(colnames(pdo_ucell))]
rownames(pdo_ucell) <- sub("^(.*?)__(.*)$", "\\1_\\2", rownames(pdo_ucell))

################################################################################
# Helper: Compute Jaccard Overlap Matrix
################################################################################
calc_jaccard_mat <- function(setsA, setsB) {
  setsA <- lapply(setsA, unique)
  setsB <- lapply(setsB, unique)
  universe <- unique(c(unlist(setsA), unlist(setsB)))
  
  jaccard_mat <- matrix(NA_real_, length(setsA), length(setsB), dimnames = list(names(setsA), names(setsB)))
  overlap_n_mat <- jaccard_mat
  pval_mat <- jaccard_mat
  
  for (i in seq_along(setsA)) {
    A <- setsA[[i]]
    for (j in seq_along(setsB)) {
      B <- setsB[[j]]
      inter <- length(intersect(A, B))
      uni <- length(union(A, B))
      overlap_n_mat[i, j] <- inter
      jaccard_mat[i, j] <- if (uni == 0) NA_real_ else inter / uni
      
      a <- inter; b <- length(setdiff(A, B)); cc <- length(setdiff(B, A)); d <- length(setdiff(universe, union(A, B)))
      pval_mat[i, j] <- if (any(c(a, b, cc, d) < 0)) NA_real_ else fisher.test(matrix(c(a, b, cc, d), nrow = 2), alternative = "greater")$p.value
    }
  }
  
  padj_mat <- matrix(p.adjust(as.vector(pval_mat), method = "BH"), nrow = nrow(pval_mat), ncol = ncol(pval_mat), dimnames = dimnames(pval_mat))
  stars_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat), dimnames = dimnames(padj_mat))
  stars_mat[padj_mat < 0.05] <- "*"
  stars_mat[padj_mat < 0.01] <- "**"
  stars_mat[padj_mat < 0.001] <- "***"
  
  display_mat <- matrix(paste0(overlap_n_mat, "\n", stars_mat), nrow = nrow(overlap_n_mat), ncol = ncol(overlap_n_mat), dimnames = dimnames(overlap_n_mat))
  
  list(jaccard = t(jaccard_mat), display = t(display_mat)) # Transpose so Parse is columns, External is rows (or vice versa depending on preference)
}

# Jaccard 1: Parse vs scATLAS
jacc_scatlas <- calc_jaccard_mat(parse_mp_genes, scatlas_mp_genes)

# Jaccard 2: Parse vs PDO
jacc_pdo <- calc_jaccard_mat(parse_mp_genes, pdo_mp_genes)

################################################################################
# Helper: Compute Spearman Correlation Matrix
################################################################################
calc_cor_mat <- function(ucellA, ucellB) {
  # ucellA: Parse (cells x MPs)
  # ucellB: External (cells x MPs)
  common_cells <- intersect(rownames(ucellA), rownames(ucellB))
  matA <- as.matrix(ucellA[common_cells, , drop = FALSE])
  matB <- as.matrix(ucellB[common_cells, , drop = FALSE])
  
  # Calculate pairwise spearman correlation between columns of matB and matA
  cor_mat <- cor(matB, matA, method = "spearman", use = "pairwise.complete.obs")
  cor_mat
}

# Cor 1: Parse vs scATLAS
cor_scatlas <- calc_cor_mat(parse_ucell, scatlas_ucell)

# Cor 2: Parse vs PDO
cor_pdo <- calc_cor_mat(parse_ucell, pdo_ucell)

################################################################################
# Generate PDFs
################################################################################
plot_colors <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)
cor_colors <- colorRampPalette(c("blue", "white", "red"))(100)

pdf(file.path(out_dir, "Auto_parse_mp_jaccard_vs_scATLAS.pdf"), width = 12, height = 10, useDingbats = FALSE)
pheatmap(jacc_scatlas$jaccard, cluster_rows = FALSE, cluster_cols = FALSE, border_color = "grey85",
         main = "Parse MPs vs scATLAS MPs (Jaccard Index)",
         angle_col = "90", display_numbers = jacc_scatlas$display, fontsize_number = 12, number_color = "black",
         fontsize_row = 14, fontsize_col = 14, fontsize = 14, color = plot_colors)
dev.off()

pdf(file.path(out_dir, "Auto_parse_mp_jaccard_vs_PDO.pdf"), width = 12, height = 10, useDingbats = FALSE)
pheatmap(jacc_pdo$jaccard, cluster_rows = FALSE, cluster_cols = FALSE, border_color = "grey85",
         main = "Parse MPs vs PDO MPs (Jaccard Index)",
         angle_col = "90", display_numbers = jacc_pdo$display, fontsize_number = 12, number_color = "black",
         fontsize_row = 14, fontsize_col = 14, fontsize = 14, color = plot_colors)
dev.off()

pdf(file.path(out_dir, "Auto_parse_mp_correlation_vs_scATLAS.pdf"), width = 12, height = 10, useDingbats = FALSE)
max_abs <- max(abs(cor_scatlas), na.rm = TRUE); if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
pheatmap(cor_scatlas, color = cor_colors, breaks = seq(-max_abs, max_abs, length.out = 101),
         display_numbers = TRUE, number_format = "%.2f", fontsize = 14, fontsize_number = 12, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Parse MPs vs scATLAS MPs (Global Spearman Correlation)", angle_col = "90")
dev.off()

pdf(file.path(out_dir, "Auto_parse_mp_correlation_vs_PDO.pdf"), width = 12, height = 10, useDingbats = FALSE)
max_abs <- max(abs(cor_pdo), na.rm = TRUE); if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
pheatmap(cor_pdo, color = cor_colors, breaks = seq(-max_abs, max_abs, length.out = 101),
         display_numbers = TRUE, number_format = "%.2f", fontsize = 14, fontsize_number = 12, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Parse MPs vs PDO MPs (Global Spearman Correlation)", angle_col = "90")
dev.off()

script_run_status <- "success"
cat("Saved 4 external comparison PDFs in:", out_dir, "\n")
