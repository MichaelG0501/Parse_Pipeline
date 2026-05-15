####################
# parse_metaprogram_internal_correlation.R
#
# Description:
#   Generates within-Parse metaprogram diagnostics: global Spearman correlation,
#   per-sample MP correlation, gene-set Jaccard overlap, and per-cell UCell
#   heatmap in one combined PDF.
#
# Inputs:
#   parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds
#   parse_outs/Auto_parse_metaprograms/Auto_parse_final_geneNMF.rds or Auto_parse_merged_geneNMF.rds
#   parse_outs/Auto_parse_metaprograms/Auto_parse_UCell_scores_filtered_nMP<k>.rds
#
# Outputs:
#   parse_outs/Auto_parse_metaprograms/Auto_parse_mp_correlation/Auto_parse_nMP<k>_analysis_combined.pdf
#   parse_outs/Auto_parse_metaprograms/Auto_parse_mp_correlation/Auto_parse_mp_correlation_summary.csv
#   parse_outs/logs/run_summaries/parse_metaprogram_internal_correlation_*.txt
#
# Cache / replot:
#   Reuses UCell scores if present; only rescoring is triggered when the selected
#   UCell matrix is missing.
#
# Methodology:
#   analysis/methodology/metaprograms/metaprogram_internal_correlation_methodology.md
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_metaprogram_internal_correlation",
  input_files = c(
    "parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds",
    "parse_outs/Auto_parse_metaprograms/Auto_parse_final_geneNMF.rds",
    "parse_outs/Auto_parse_metaprograms/Auto_parse_UCell_scores_filtered_nMP<k>.rds"
  ),
  output_files = c(
    "parse_outs/Auto_parse_metaprograms/Auto_parse_mp_correlation/Auto_parse_nMP<k>_analysis_combined.pdf",
    "parse_outs/Auto_parse_metaprograms/Auto_parse_mp_correlation/Auto_parse_mp_correlation_summary.csv"
  )
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(SeuratObject)
  library(ComplexHeatmap)
  library(circlize)
  library(pheatmap)
  library(dplyr)
  library(grid)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath(file.path(parse_project_root(), "analysis/metaprograms/parse_metaprogram_internal_correlation.R"))
}

script_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(script_dir, "..", ".."))
out_dir_root <- file.path(project_dir, "parse_outs")
setwd(out_dir_root)

parse_mp_dir <- file.path("Auto_parse_metaprograms")
out_dir <- file.path(parse_mp_dir, "Auto_parse_mp_correlation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

geneNMF.metaprograms <- readRDS(file.path(parse_mp_dir, "Auto_parse_MP_outs_default.rds"))
optimal_txt <- file.path(parse_mp_dir, "Auto_parse_optimal_nMP.txt")
optimal_nMP <- if (file.exists(optimal_txt)) readLines(optimal_txt, warn = FALSE)[1] else length(geneNMF.metaprograms$metaprograms.genes)

parse_obj_path <- file.path(parse_mp_dir, "Auto_parse_final_geneNMF.rds")
ucell_path <- file.path(parse_mp_dir, paste0("Auto_parse_UCell_scores_filtered_nMP", optimal_nMP, ".rds"))

if (!file.exists(parse_obj_path)) {
  parse_obj_path <- file.path(parse_mp_dir, "Auto_parse_merged_geneNMF.rds")
}
if (!file.exists(parse_obj_path)) {
  stop("Missing Parse merged/final GeneNMF object.")
}

parse_obj <- readRDS(parse_obj_path)

mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
mp.genes <- mp.genes[!names(mp.genes) %in% names(coverage_tbl)[coverage_tbl < 0.25]]

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", names(mp.genes)))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)

# Force UCell scoring if needed
if (!file.exists(ucell_path)) {
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    parse_obj <- SeuratObject::JoinLayers(parse_obj, assay = "RNA")
  }
  parse_obj <- AddModuleScore_UCell(parse_obj, features = mp.genes, ncores = 1, name = "")
  ucell_scores <- parse_obj@meta.data[, names(mp.genes), drop = FALSE]
  saveRDS(ucell_scores, ucell_path, compress = FALSE)
} else {
  ucell_scores <- readRDS(ucell_path)
}

common_cells <- intersect(Cells(parse_obj), rownames(ucell_scores))
parse_obj <- parse_obj[, common_cells]
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]

canon_mps_avail <- mp_tree_order_names[mp_tree_order_names %in% colnames(ucell_scores)]
raw_scores <- ucell_scores[, canon_mps_avail, drop = FALSE]
mod_mat_raw <- t(as.matrix(raw_scores))
mod_mat_scaled <- t(scale(as.matrix(raw_scores)))
rownames(mod_mat_raw) <- canon_mps_avail
rownames(mod_mat_scaled) <- canon_mps_avail

samples_vec <- parse_obj$orig.ident
samples <- unique(samples_vec)
mp_names <- rownames(mod_mat_raw)
message("MP names: ", paste(mp_names, collapse=", "))
n_mps <- length(mp_names)

# 1) Global correlation
cor_global <- cor(t(mod_mat_raw), method = "spearman", use = "pairwise.complete.obs")

# 2) Per-sample correlation
cor_array <- array(NA_real_, dim = c(n_mps, n_mps, length(samples)), dimnames = list(mp_names, mp_names, samples))
for (samp in samples) {
  cells_in_sample <- colnames(mod_mat_raw)[samples_vec == samp]
  if (length(cells_in_sample) < 10) next
  sub_mat <- mod_mat_raw[, cells_in_sample, drop = FALSE]
  cor_array[, , samp] <- cor(t(sub_mat), method = "spearman", use = "pairwise.complete.obs")
}
z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))
mean_rho <- matrix(NA_real_, n_mps, n_mps, dimnames = list(mp_names, mp_names))
p_vals <- matrix(NA_real_, n_mps, n_mps, dimnames = list(mp_names, mp_names))
for (i in seq_len(n_mps)) {
  for (j in seq_len(n_mps)) {
    if (i == j) { mean_rho[i, j] <- 1; p_vals[i, j] <- 0 } else {
      z_scores <- z_array[i, j, ]; z_scores <- z_scores[is.finite(z_scores)]
      if (length(z_scores) >= 3) {
        mean_rho[i, j] <- tanh(mean(z_scores))
        tt <- tryCatch(t.test(z_scores), error = function(e) NULL)
        p_vals[i, j] <- if (is.null(tt)) NA_real_ else tt$p.value
      }
    }
  }
}

# 3) Jaccard
canon_gene_sets <- mp.genes[mp_tree_order_names[mp_tree_order_names %in% names(mp.genes)]]
all_sets <- lapply(canon_gene_sets, unique)
universe <- unique(unlist(all_sets))
jaccard_mat <- matrix(NA_real_, length(all_sets), length(all_sets), dimnames = list(names(all_sets), names(all_sets)))
overlap_n_mat <- jaccard_mat; pval_mat <- jaccard_mat
for (i in seq_along(all_sets)) {
  A <- all_sets[[i]]
  for (j in seq_along(all_sets)) {
    B <- all_sets[[j]]
    inter <- length(intersect(A, B)); uni <- length(union(A, B))
    overlap_n_mat[i, j] <- inter; jaccard_mat[i, j] <- if (uni == 0) NA_real_ else inter / uni
    a <- inter; b <- length(setdiff(A, B)); cc <- length(setdiff(B, A)); d <- length(setdiff(universe, union(A, B)))
    pval_mat[i, j] <- if (any(c(a, b, cc, d) < 0)) NA_real_ else fisher.test(matrix(c(a, b, cc, d), nrow = 2), alternative = "greater")$p.value
  }
}
padj_mat <- matrix(p.adjust(as.vector(pval_mat), method = "BH"), nrow = nrow(pval_mat), ncol = ncol(pval_mat))
stars_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat))
stars_mat[padj_mat < 0.05] <- "*"; stars_mat[padj_mat < 0.01] <- "**"; stars_mat[padj_mat < 0.001] <- "***"
display_mat_jaccard <- matrix(paste0(overlap_n_mat, "\n", stars_mat), nrow = nrow(overlap_n_mat), ncol = ncol(overlap_n_mat))

# 4) Per-cell heatmap setup
set.seed(42); MAX_CELLS <- 10000
cells_to_plot <- if (ncol(mod_mat_scaled) > MAX_CELLS) sample(colnames(mod_mat_scaled), MAX_CELLS) else colnames(mod_mat_scaled)
sub_scores <- mod_mat_scaled[mp_names, cells_to_plot, drop = FALSE]
lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

# Sample annotation
sample_order <- c("T0", "T1", "T2", "T4", "R4", "eR4")
sample_labels_raw <- parse_obj$orig.ident[match(cells_to_plot, Cells(parse_obj))]
sample_labels_raw <- factor(sample_labels_raw, levels = intersect(sample_order, unique(sample_labels_raw)))

# Pre-order cells and scores by sample order
order_idx <- order(sample_labels_raw)
cells_to_plot <- cells_to_plot[order_idx]
sample_labels <- sample_labels_raw[order_idx]
sub_scores <- mod_mat_scaled[mp_names, cells_to_plot, drop = FALSE]

ha_column <- HeatmapAnnotation(
  Sample = sample_labels,
  col = list(Sample = setNames(hcl.colors(length(levels(sample_labels)), palette = "Dark 3"), levels(sample_labels))),
  show_annotation_name = FALSE
)

# 5) PDF Output
pdf(file.path(out_dir, paste0("Auto_parse_nMP", optimal_nMP, "_analysis_combined.pdf")), width = 14, height = 12, useDingbats = FALSE)
# A) Global correlation
max_abs <- max(abs(cor_global), na.rm = TRUE); if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
pheatmap(cor_global, color = colorRampPalette(c("blue", "white", "red"))(100), breaks = seq(-max_abs, max_abs, length.out = 101),
         display_numbers = TRUE, number_format = "%.2f", fontsize = 14, fontsize_number = 12, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Global MP-MP Spearman correlation (all cells)", angle_col = "90")
# B) Per-sample correlation
ht_cor <- Heatmap(mean_rho, name = "Mean Rho", col = colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red")),
                  rect_gp = gpar(col = "white", lwd = 1), cluster_rows = FALSE, cluster_columns = FALSE,
                  row_names_side = "right", column_names_side = "bottom", column_title_rot = 20, column_title_side = "top",
                  column_title_gp = gpar(fontsize = 12, fontface = "bold"), row_names_gp = gpar(fontsize = 14, fontface = "bold"), column_names_gp = gpar(fontsize = 14, fontface = "bold"),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    p <- p_vals[i, j]; rho <- mean_rho[i, j]
                    if (!is.na(p) && !is.na(rho)) {
                      stars <- if (p < 0.001) "\n***" else if (p < 0.01) "\n**" else if (p < 0.05) "\n*" else ""
                      grid.text(paste0(round(rho, 2), stars), x, y, gp = gpar(fontsize = 14))
                    }
                  })
draw(ht_cor)
# C) Jaccard
pheatmap(jaccard_mat, cluster_rows = FALSE, cluster_cols = FALSE, border_color = "grey85", main = "Parse MP Gene Set Overlap (Jaccard Index)",
         angle_col = "90", display_numbers = display_mat_jaccard, fontsize_number = 12, number_color = "black",
         fontsize_row = 14, fontsize_col = 14, fontsize = 14, color = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100))
# D) Per-cell heatmap
ht_cells <- Heatmap(sub_scores, name = "Scaled UCell", col = col_fun_sc, 
                    cluster_rows = FALSE, cluster_columns = TRUE, 
                    column_split = sample_labels, cluster_column_slices = FALSE,
                    top_annotation = ha_column,
                    show_row_dend = FALSE, show_column_dend = FALSE, row_names_side = "left", 
                    row_names_gp = gpar(fontsize = 12, fontface = "bold"),
                    show_column_names = FALSE, use_raster = TRUE, raster_quality = 5, border = TRUE,
                    column_title = "Per-cell MP Expression (Clustered by Sample)", 
                    column_title_gp = gpar(fontsize = 14, fontface = "bold"))
draw(ht_cells)
dev.off()

write.csv(data.frame(n_samples = length(samples), n_parse_mps = length(all_sets), correlation_heatmap_rows = nrow(mod_mat_raw)),
          file.path(out_dir, "Auto_parse_mp_correlation_summary.csv"), row.names = FALSE)
script_run_status <- "success"
cat("Saved combined PDF in:", out_dir, "\n")
