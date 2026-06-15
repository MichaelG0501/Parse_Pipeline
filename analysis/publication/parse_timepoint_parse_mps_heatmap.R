####################
# parse_timepoint_parse_mps_heatmap.R
#
# Description:
#   Publication-quality heatmap showing the high-resolution NMF MP cluster
#   mean UCell scores across timepoints (T0, T1, T2, T4, R4, eR4), using the 
#   updated Parse MP list. Nature-contract figure style.
#
# Inputs:
#   parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_UCell_scores_nMP117.rds
#   parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds
#
# Outputs:
#   parse_outs/publication/timepoint_mps_heatmap/figures/
#     timepoint_parse_mps_heatmap.pdf
#     timepoint_parse_mps_heatmap.png
#
# Methodology:
#   analysis/methodology/publication/timepoint_mps_heatmap_methodology.md
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source("analysis/common/parse_pipeline_config.R")

# ── Paths ────────────────────────────────────────────────────────────────────
project_dir <- parse_project_root()
paths       <- parse_paths(project_dir)
out_base    <- file.path(paths$parse_outs, "publication", "timepoint_mps_heatmap")
fig_dir     <- file.path(out_base, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

setwd(paths$parse_outs)

# ── Definitions ──────────────────────────────────────────────────────────────
timepoint_levels <- c("T0", "T1", "T2", "T4", "R4", "eR4")

# parse_sample_colours is sourced from parse_pipeline_config.R
timepoint_cols <- parse_sample_colours[timepoint_levels]

increase_clusters <- list(
  "EMT/Remodeling" = c("MP23", "MP46", "MP81"),
  "Stress response" = c("MP13", "MP25"),
  "Stem-like plasticity" = c("MP16"),
  "G2/M" = c("MP36", "MP40", "MP78", "MP80", "MP94"),
  "Epithelial differentiation" = c("MP72", "MP77"),
  "Inflammatory-associated\n  mitotic program" = c("MP53", "MP89"),
  "Proliferative with angiogenic\nremodeling/ plasticity" = c("MP92", "MP42", "MP98", "MP55"),
  "Angiogenic / Inflammatory remodeling" = c("MP35")
)

decrease_clusters <- list(
  "Proliferation / DNA Replication​" = c("MP3", "MP27", "MP29"),
  "Epithelial Differentiation​" = c("MP47")
)


all_clusters <- c(increase_clusters, decrease_clusters)
cluster_direction <- c(
  rep("Enriched in Persisters", length(increase_clusters)),
  rep("Depleted in Persisters", length(decrease_clusters))
)
names(cluster_direction) <- names(all_clusters)

# ── Load Data ────────────────────────────────────────────────────────────────
message("Loading UCell scores...")
ucell_cache <- "Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_UCell_scores_nMP117.rds"
if (!file.exists(ucell_cache)) stop("Missing cached UCell scores: ", ucell_cache)
ucell <- readRDS(ucell_cache)
# Fix any prepended names from old cache or previous runs
rownames(ucell) <- sub("^([A-Za-z0-9]+)_", "\\1__", rownames(ucell))

message("Loading cell states (for metadata)...")
states <- readRDS("cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds")

common_cells <- intersect(rownames(ucell), states$cell)

cell_df <- states %>%
  filter(cell %in% common_cells) %>%
  filter(sample %in% timepoint_levels) %>%
  mutate(sample = factor(sample, levels = timepoint_levels))

# ── Data Processing ──────────────────────────────────────────────────────────
mp_to_cluster <- data.frame(
  MP = unname(unlist(all_clusters)),
  Cluster = rep(names(all_clusters), lengths(all_clusters)),
  stringsAsFactors = FALSE
)
mp_to_cluster$Direction <- cluster_direction[mp_to_cluster$Cluster]

available_mps <- intersect(mp_to_cluster$MP, colnames(ucell))
mp_to_cluster <- mp_to_cluster[mp_to_cluster$MP %in% available_mps, ]
unique_mps <- unique(mp_to_cluster$MP)

mp_long <- as.data.frame(ucell[cell_df$cell, unique_mps, drop=FALSE]) %>%
  mutate(cell = cell_df$cell, sample = cell_df$sample) %>%
  pivot_longer(cols = all_of(unique_mps), names_to = "MP", values_to = "score")

mp_agg <- mp_long %>%
  group_by(sample, MP) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

mp_mean_matrix <- mp_agg %>%
  pivot_wider(names_from = sample, values_from = mean_score) %>%
  left_join(mp_to_cluster, by = "MP", relationship = "many-to-many")

mp_inc <- mp_mean_matrix %>% filter(Direction == "Enriched in Persisters")
mp_dec <- mp_mean_matrix %>% filter(Direction == "Depleted in Persisters")

# ── Nature-contract Drawing ──────────────────────────────────────────────────
create_mp_ht <- function(df, title) {
  if (nrow(df) == 0) return(NULL)
  
  mat_raw <- df %>% select(MP, all_of(timepoint_levels)) %>% distinct() %>% column_to_rownames("MP") %>% as.matrix()
  
  # Ensure column order
  mat_raw <- mat_raw[, timepoint_levels[timepoint_levels %in% colnames(mat_raw)], drop = FALSE]
  
  # Row z-score
  mat_z <- t(scale(t(mat_raw)))
  colnames(mat_z) <- colnames(mat_raw)
  
  df_unique <- df %>% select(MP, Cluster) %>% distinct()
  cluster_lvls <- unique(names(all_clusters))
  df_unique <- df_unique %>% mutate(Cluster = factor(Cluster, levels = cluster_lvls)) %>% arrange(Cluster, MP)
  
  mat_raw <- mat_raw[df_unique$MP, , drop = FALSE]
  mat_z <- mat_z[df_unique$MP, , drop = FALSE]
  rownames(mat_raw) <- df_unique$MP
  rownames(mat_z) <- df_unique$MP
  
  row_split_fac <- factor(df_unique$Cluster, levels = cluster_lvls[cluster_lvls %in% df_unique$Cluster])
  
  clip <- 2
  col_f <- colorRamp2(c(-clip, 0, clip), c("#245F7B", "white", "#B63E2F"))
  
  cluster_colors <- rainbow(length(cluster_lvls))
  names(cluster_colors) <- cluster_lvls
  
  ha_row <- rowAnnotation(
    Cluster = df_unique$Cluster,
    col = list(Cluster = cluster_colors),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  
  cell_f <- function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.3f", mat_raw[i, j]), x, y, gp = gpar(fontsize = 7, fontfamily = "Arial"))
  }
  
  Heatmap(
    mat_z,
    name = "Z-Score",
    col = col_f,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
    column_names_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "Arial"),
    column_names_rot = 0,
    column_names_centered = TRUE,
    row_split = row_split_fac,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
    row_gap = unit(3, "mm"),
    cell_fun = cell_f,
    left_annotation = ha_row,
    column_names_side = "top",
    show_heatmap_legend = (title == "Enriched in Persisters"),
    heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"),
    width = unit(ncol(mat_z) * 1.5, "cm"),
    height = unit(nrow(mat_z) * 0.8, "cm"),
    column_title = paste("MPs", title),
    column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial")
  )
}

message("Generating MP-level timepoint heatmap ...")
ht_inc <- create_mp_ht(mp_inc, "Enriched in Persisters")
ht_dec <- create_mp_ht(mp_dec, "Depleted in Persisters")

pdf_path <- file.path(fig_dir, "timepoint_parse_mps_heatmap.pdf")
png_path <- file.path(fig_dir, "timepoint_parse_mps_heatmap.png")

save_pub_pdf <- function(path, width_mm = 240, height_mm = 140) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  grDevices::cairo_pdf(path, width = w, height = h, family = "Arial")
  pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  if (!is.null(ht_inc)) draw(ht_inc, newpage = FALSE, heatmap_legend_side = "bottom")
  popViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  if (!is.null(ht_dec)) draw(ht_dec, newpage = FALSE, heatmap_legend_side = "bottom")
  popViewport()
  dev.off()
}

message("Saving PDF: ", pdf_path)
save_pub_pdf(pdf_path, width_mm = 450, height_mm = 200)

message("Saving PNG: ", png_path)
png(png_path, width = 450, height = 200, units = "mm", res = 600, type = "cairo")
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
if (!is.null(ht_inc)) draw(ht_inc, newpage = FALSE, heatmap_legend_side = "bottom")
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
if (!is.null(ht_dec)) draw(ht_dec, newpage = FALSE, heatmap_legend_side = "bottom")
popViewport()
dev.off()

message("Done.")
