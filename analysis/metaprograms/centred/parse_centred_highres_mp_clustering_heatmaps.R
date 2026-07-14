####################
# parse_centred_highres_mp_clustering_heatmaps.R
#
# Description:
#   Clusters the retained T2/T4 high MPs and decrease+consistent MPs based on
#   their gene list similarity (Jaccard) and MP expression correlation (Fisher-Z
#   transformed Spearman rho). Replicates the visual style of the scRef pipeline
#   centred MP heatmaps. Enforces distinct grouping: T2/T4 High vs Decrease MPs.
#   Clustering is performed based on correlation within each group, and this
#   exact order is applied to the Jaccard heatmap.
#
# Outputs:
#   - live/.../clustering/centred_highres_mp_jaccard_heatmap.pdf
#   - live/.../clustering/centred_highres_mp_correlation_heatmap.pdf
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
})

source("analysis/common/parse_pipeline_config.R")

project_dir <- parse_project_root()
qc_dir <- file.path(project_dir, "parse_outs")
base_highres_dir <- file.path(qc_dir, "centred", "Auto_parse_highres_metaprogram_trends")

live_project_dir <- "/rds/general/project/spatialtranscriptomics/live/Parse_Pipeline"
live_outdir <- file.path(live_project_dir, "parse_outs", "centred", "Auto_parse_highres_metaprogram_trends", "clustering")
dir.create(live_outdir, recursive = TRUE, showWarnings = FALSE)

# 1. Identify Retained MPs
t2t4_file <- file.path(base_highres_dir, "Auto_T2T4_gt_T0eR4_filter", "Auto_parse_highres_T2T4_filter_summary_nMP117.csv")
strict_file <- file.path(base_highres_dir, "Auto_parse_highres_trend_summary_nMP117.csv")

t2t4 <- read.csv(t2t4_file, stringsAsFactors = FALSE)
strict <- read.csv(strict_file, stringsAsFactors = FALSE)

t2t4_mps <- t2t4$MP[t2t4$retained == TRUE]
strict_mps <- strict$MP[strict$retained == TRUE & strict$trend_type_label == "Decrease, Consistent"]

# Guarantee strict distinction
strict_mps <- setdiff(strict_mps, t2t4_mps)

retained_mps <- c(t2t4_mps, strict_mps)
cat("Found", length(retained_mps), "unique MPs to cluster.\n")

# Load MP Enrichment Annotations exactly as done in the pipeline
anno_strict_path <- file.path(base_highres_dir, "Auto_parse_highres_top_3CA_noncellcycle_nMP117.csv")
anno_t2t4_path <- file.path(base_highres_dir, "Auto_T2T4_gt_T0eR4_filter", "Auto_parse_highres_T2T4_top_3CA_noncellcycle_nMP117.csv")

if (file.exists(anno_strict_path)) {
  label_strict <- read.csv(anno_strict_path, stringsAsFactors = FALSE)
  strict <- left_join(strict, label_strict, by = "MP")
} else {
  strict$top_3ca_noncc <- NA_character_
}
strict <- strict %>%
  mutate(
    top_3ca_noncc = ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc),
    display_label = paste0(MP, " - ", top_3ca_noncc)
  )

if (file.exists(anno_t2t4_path)) {
  label_t2t4 <- read.csv(anno_t2t4_path, stringsAsFactors = FALSE)
  t2t4 <- left_join(t2t4, label_t2t4, by = "MP")
} else {
  t2t4$top_3ca_noncc <- NA_character_
}
t2t4 <- t2t4 %>%
  mutate(
    top_3ca_noncc = ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc),
    display_label = paste0(MP, " - ", top_3ca_noncc)
  )

label_map <- setNames(c(strict$display_label, t2t4$display_label), c(strict$MP, t2t4$MP))
mp_labels <- label_map[retained_mps]

# Source Group vector
source_group <- factor(
  ifelse(retained_mps %in% t2t4_mps, "T2/T4 High", "Decrease/Consistent"),
  levels = c("T2/T4 High", "Decrease/Consistent")
)
names(source_group) <- retained_mps

# Heatmap Annotations
ha <- HeatmapAnnotation(
  Group = source_group,
  col = list(Group = c("T2/T4 High" = "#E41A1C", "Decrease/Consistent" = "#377EB8")),
  show_legend = TRUE,
  annotation_legend_param = list(Group = list(title = "Trend Origin", title_gp = gpar(fontface = "bold")))
)

row_ha <- rowAnnotation(
  Group = source_group,
  col = list(Group = c("T2/T4 High" = "#E41A1C", "Decrease/Consistent" = "#377EB8")),
  show_legend = FALSE
)

# Load geneNMF to get all MP gene lists
geneNMF_path <- file.path(base_highres_dir, "Auto_parse_highres_geneNMF_metaprograms_nMP117.rds")
geneNMF.metaprograms <- readRDS(geneNMF_path)
mp_genes <- geneNMF.metaprograms$metaprograms.genes
retained_genes <- mp_genes[retained_mps]

# 2. Compute Jaccard Similarity Matrix
n_mps <- length(retained_mps)
jaccard_mat <- matrix(0, n_mps, n_mps, dimnames = list(retained_mps, retained_mps))
for (i in seq_len(n_mps)) {
  for (j in seq_len(n_mps)) {
    g1 <- retained_genes[[i]]
    g2 <- retained_genes[[j]]
    jaccard_mat[i, j] <- length(intersect(g1, g2)) / length(union(g1, g2))
  }
}
rownames(jaccard_mat) <- mp_labels[rownames(jaccard_mat)]
colnames(jaccard_mat) <- mp_labels[colnames(jaccard_mat)]

# 3. Load UCell scores and calculate Fisher-Z Correlation Matrix
ucell_path <- file.path(base_highres_dir, "Auto_parse_highres_UCell_scores_nMP117.rds")
meta_path <- file.path(base_highres_dir, "Auto_parse_highres_cell_metadata_nMP117.rds")

ucell_scores <- readRDS(ucell_path)
cell_meta <- readRDS(meta_path)
score_mat <- scale(as.matrix(ucell_scores[, retained_mps, drop = FALSE]))

cell_ids <- rownames(score_mat)
sample_vec <- cell_meta$sample[match(cell_ids, cell_meta$cell)]
sample_vec[is.na(sample_vec)] <- cell_ids[is.na(sample_vec)]
samples <- unique(sample_vec)

cor_array <- array(
  NA_real_,
  dim = c(n_mps, n_mps, length(samples)),
  dimnames = list(retained_mps, retained_mps, samples)
)

cat("Computing per-sample UCell correlations...\n")
for (samp in samples) {
  cells_in_sample <- which(sample_vec == samp)
  if (length(cells_in_sample) < 10) next
  sub_mat <- score_mat[cells_in_sample, , drop = FALSE]
  cor_array[, , samp] <- cor(sub_mat, method = "spearman")
}

z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))
mean_rho <- matrix(NA_real_, n_mps, n_mps, dimnames = list(retained_mps, retained_mps))
p_vals <- matrix(NA_real_, n_mps, n_mps, dimnames = list(retained_mps, retained_mps))

for (i in seq_len(n_mps)) {
  for (j in seq_len(n_mps)) {
    if (i == j) {
      mean_rho[i, j] <- 1
      p_vals[i, j] <- 0
      next
    }
    z_scores <- z_array[i, j, ]
    z_scores <- z_scores[is.finite(z_scores)]
    if (length(z_scores) < 3) next
    mean_rho[i, j] <- tanh(mean(z_scores))
    test_res <- tryCatch(t.test(z_scores), error = function(e) NULL)
    p_vals[i, j] <- if (!is.null(test_res)) test_res$p.value else NA_real_
  }
}
rownames(mean_rho) <- mp_labels[rownames(mean_rho)]
colnames(mean_rho) <- mp_labels[colnames(mean_rho)]

# 4. Generate Expression Correlation Heatmap (which defines the order)
cat("Generating Correlation Heatmap...\n")
col_cor <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))

ht_cor <- Heatmap(
  mean_rho,
  name = paste0("Mean Rho\n(", sum(apply(cor_array, 3, function(x) any(is.finite(x)))), " Samples)"),
  col = col_cor,
  rect_gp = gpar(col = "white", lwd = 1),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_split = source_group,
  column_split = source_group,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  top_annotation = ha,
  left_annotation = row_ha,
  column_title = paste0("Centred MP Expression Correlation (n = ", n_mps, ")"),
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  row_names_side = "right",
  column_names_side = "bottom",
  column_names_rot = 30,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  width = unit(10.5, "inch"),
  height = unit(10.5, "inch"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    p <- p_vals[i, j]
    rho <- mean_rho[i, j]
    if (is.na(p) || is.na(rho)) {
      grid.text("NA", x, y, gp = gpar(fontsize = 6, col = "grey50"))
    } else if (p < 0.001) {
      grid.text(paste0(round(rho, 2), "\n***"), x, y, gp = gpar(fontsize = 6))
    } else if (p < 0.01) {
      grid.text(paste0(round(rho, 2), "\n**"), x, y, gp = gpar(fontsize = 6))
    } else if (p < 0.05) {
      grid.text(paste0(round(rho, 2), "\n*"), x, y, gp = gpar(fontsize = 6))
    } else {
      grid.text(round(rho, 2), x, y, gp = gpar(fontsize = 6))
    }
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14)
  )
)

cor_pdf_path <- file.path(live_outdir, "centred_highres_mp_correlation_heatmap.pdf")
cor_png_path <- file.path(live_outdir, "centred_highres_mp_correlation_heatmap.png")

cairo_pdf(cor_pdf_path, width = 23.5, height = 18)
ht_cor_drawn <- draw(ht_cor, heatmap_legend_side = "left", padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()

png(cor_png_path, width = 23.5, height = 18, units = "in", res = 350, bg = "white")
draw(ht_cor_drawn, heatmap_legend_side = "left", padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()

# Extract clustering order from the correlation heatmap to enforce on Jaccard heatmap
ro <- row_order(ht_cor_drawn)
ordered_indices <- unlist(ro)

# 4b. Generate Clean Correlation Heatmap (no text, no labels, jaccard size)
ht_cor_clean <- Heatmap(
  mean_rho,
  name = paste0("Mean Rho\n(", sum(apply(cor_array, 3, function(x) any(is.finite(x)))), " Samples)"),
  col = col_cor,
  rect_gp = gpar(col = "white", lwd = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = ordered_indices,
  column_order = ordered_indices,
  row_split = source_group,
  column_split = source_group,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  top_annotation = ha,
  left_annotation = row_ha,
  column_title = paste0("Centred MP Expression Correlation (n = ", n_mps, ")"),
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  show_row_names = FALSE,
  show_column_names = FALSE,
  width = unit(205, "mm"),
  height = unit(205, "mm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    direction = "horizontal",
    legend_width = unit(8, "cm")
  )
)

cor_clean_pdf_path <- file.path(live_outdir, "centred_highres_mp_correlation_heatmap_clean.pdf")
cor_clean_png_path <- file.path(live_outdir, "centred_highres_mp_correlation_heatmap_clean.png")

cairo_pdf(cor_clean_pdf_path, width = 17.5, height = 12.8)
draw(ht_cor_clean, heatmap_legend_side = "bottom", padding = unit(c(12, 10, 14, 10), "mm"))
dev.off()

png(cor_clean_png_path, width = 17.5, height = 12.8, units = "in", res = 450, bg = "white")
draw(ht_cor_clean, heatmap_legend_side = "bottom", padding = unit(c(12, 10, 14, 10), "mm"))
dev.off()

# 5. Generate Jaccard Heatmap
cat("Generating Jaccard Heatmap...\n")
col_fun_nmf <- colorRamp2(c(0, max(jaccard_mat[upper.tri(jaccard_mat)])/2, max(jaccard_mat[upper.tri(jaccard_mat)])), c("#FFF7F3", "#FB6A4A", "#67000D"))

ht_jaccard <- Heatmap(
  jaccard_mat,
  name = "Jaccard",
  col = col_fun_nmf,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = source_group,
  column_split = source_group,
  row_order = ordered_indices,
  column_order = ordered_indices,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  top_annotation = ha,
  left_annotation = row_ha,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 30,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title = paste0("Centred MP Gene List Similarity (n = ", n_mps, ")"),
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  use_raster = TRUE,
  raster_quality = 4,
  border = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    direction = "horizontal",
    legend_width = unit(8, "cm")
  ),
  width = unit(205, "mm"),
  height = unit(205, "mm")
)

pdf_path_jaccard <- file.path(live_outdir, "centred_highres_mp_jaccard_heatmap.pdf")
png_path_jaccard <- file.path(live_outdir, "centred_highres_mp_jaccard_heatmap.png")

cairo_pdf(pdf_path_jaccard, width = 17.5, height = 12.8)
draw(ht_jaccard, heatmap_legend_side = "bottom", padding = unit(c(12, 10, 14, 10), "mm"))
dev.off()

png(png_path_jaccard, width = 17.5, height = 12.8, units = "in", res = 450, bg = "white")
draw(ht_jaccard, heatmap_legend_side = "bottom", padding = unit(c(12, 10, 14, 10), "mm"))
dev.off()

cat("Finished generating heatmaps.\n")
