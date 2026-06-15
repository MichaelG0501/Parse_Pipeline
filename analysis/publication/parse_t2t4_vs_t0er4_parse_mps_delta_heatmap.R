####################
# parse_t2t4_vs_t0er4_parse_mps_delta_heatmap.R
#
# Description:
#   Publication-quality heatmap showing the high-resolution NMF MP cluster
#   delta score changes (T2T4 vs T0eR4 per state), using the updated Parse MP list.
#   Nature-contract figure style.
#
# Inputs:
#   parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_UCell_scores_nMP117.rds
#   parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds
#
# Outputs:
#   parse_outs/publication/t2t4_vs_t0er4_mp_delta/figures/
#     t2t4_vs_t0er4_parse_mps_delta_heatmap.pdf
#     t2t4_vs_t0er4_parse_mps_delta_heatmap.png
#
# Methodology:
#   analysis/methodology/publication/t2t4_vs_t0er4_mp_delta_heatmap_methodology.md
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
out_base    <- file.path(paths$parse_outs, "publication", "t2t4_vs_t0er4_mp_delta")
fig_dir     <- file.path(out_base, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

setwd(paths$parse_outs)

# ── Definitions ──────────────────────────────────────────────────────────────
state_levels <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation"
)

state_labels_split <- c(
  "Classic\nProliferative",
  "Basal to\nIntest. Meta",
  "SMG-like\nMetaplasia",
  "Stress-\nadaptive",
  "3CA EMT &\nProt. mat."
)
names(state_labels_split) <- state_levels

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Stress-adaptive"       = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8"
)

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
  rep("Increased in T2/T4", length(increase_clusters)),
  rep("Decreased in T2/T4", length(decrease_clusters))
)
names(cluster_direction) <- names(all_clusters)

# ── Load Data ────────────────────────────────────────────────────────────────
message("Loading UCell scores...")
ucell_cache <- "Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_UCell_scores_nMP117.rds"
if (!file.exists(ucell_cache)) stop("Missing cached UCell scores: ", ucell_cache)
ucell <- readRDS(ucell_cache)
# Fix any prepended names from old cache or previous runs
rownames(ucell) <- sub("^([A-Za-z0-9]+)_", "\\1__", rownames(ucell))

message("Loading cell states...")
states <- readRDS("cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds")
states <- states %>% mutate(state = pdo_state)
emt_prot_states <- c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III")
states$state[states$state %in% emt_prot_states] <- "3CA_EMT_and_Protein_maturation"

common_cells <- intersect(rownames(ucell), states$cell)

cell_df <- states %>%
  filter(cell %in% common_cells) %>%
  mutate(
    Treatment = ifelse(sample %in% c("T2", "T4"), "T2T4", 
                       ifelse(sample %in% c("T0", "eR4"), "T0eR4", NA))
  ) %>%
  filter(!is.na(Treatment)) %>%
  filter(state %in% state_levels)

state_levels <- intersect(state_levels, unique(cell_df$state))

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
  mutate(cell = cell_df$cell, state = cell_df$state, Treatment = cell_df$Treatment) %>%
  pivot_longer(cols = all_of(unique_mps), names_to = "MP", values_to = "score")

mp_agg <- mp_long %>%
  group_by(state, Treatment, MP) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

mp_mean_delta <- mp_agg %>%
  pivot_wider(names_from = Treatment, values_from = mean_score) %>%
  filter(!is.na(T0eR4), !is.na(T2T4)) %>%
  mutate(mean_delta = T2T4 - T0eR4) %>%
  left_join(
    mp_long %>%
      group_by(state, MP) %>%
      summarise(
        p_value = tryCatch(t.test(score ~ Treatment)$p.value, error = function(e) NA_real_),
        .groups = "drop"
      ),
    by = c("state", "MP")
  ) %>%
  mutate(
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  left_join(mp_to_cluster, by = "MP", relationship = "many-to-many")

mp_inc <- mp_mean_delta %>% filter(Direction == "Increased in T2/T4")
mp_dec <- mp_mean_delta %>% filter(Direction == "Decreased in T2/T4")

# ── Nature-contract Drawing ──────────────────────────────────────────────────
# Adjusted for publication standard: Helvetica/Arial font, appropriate sizes.
create_mp_ht <- function(df, title) {
  if (nrow(df) == 0) return(NULL)
  
  mat <- df %>% select(MP, state, mean_delta) %>% distinct() %>% pivot_wider(names_from = state, values_from = mean_delta) %>% column_to_rownames("MP") %>% as.matrix()
  mat <- mat[, state_levels[state_levels %in% colnames(mat)], drop = FALSE]
  
  sig <- df %>% select(MP, state, sig_label) %>% distinct() %>% pivot_wider(names_from = state, values_from = sig_label) %>% column_to_rownames("MP") %>% as.matrix()
  sig <- sig[rownames(mat), colnames(mat), drop = FALSE]
  sig[is.na(sig)] <- ""
  
  df_unique <- df %>% select(MP, Cluster) %>% distinct()
  cluster_lvls <- unique(names(all_clusters))
  df_unique <- df_unique %>% mutate(Cluster = factor(Cluster, levels = cluster_lvls)) %>% arrange(Cluster, MP)
  
  mat <- mat[df_unique$MP, , drop = FALSE]
  sig <- sig[df_unique$MP, , drop = FALSE]
  rownames(mat) <- df_unique$MP
  rownames(sig) <- df_unique$MP
  
  row_split_fac <- factor(df_unique$Cluster, levels = cluster_lvls[cluster_lvls %in% df_unique$Cluster])
  
  clip <- max(0.002, quantile(abs(mat), 0.95, na.rm = TRUE))
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
    grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 7, fontfamily = "Arial"))
  }
  
  Heatmap(
    mat,
    name = "\u0394 Score",
    col = col_f,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
    column_names_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
    column_labels = paste0("\n", state_labels_split[colnames(mat)]),
    column_names_rot = 0,
    column_names_centered = TRUE,
    row_split = row_split_fac,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
    row_gap = unit(3, "mm"),
    cell_fun = cell_f,
    left_annotation = ha_row,
    column_names_side = "top",
    show_heatmap_legend = (title == "Increased in T2/T4"),
    heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"),
    width = unit(ncol(mat) * 2.5, "cm"),
    height = unit(nrow(mat) * 0.8, "cm"),
    column_title = paste("MPs", title),
    column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial")
  )
}

message("Generating MP-level delta heatmap ...")
ht_inc <- create_mp_ht(mp_inc, "Increased in T2/T4")
ht_dec <- create_mp_ht(mp_dec, "Decreased in T2/T4")

pdf_path <- file.path(fig_dir, "t2t4_vs_t0er4_parse_mps_delta_heatmap.pdf")
png_path <- file.path(fig_dir, "t2t4_vs_t0er4_parse_mps_delta_heatmap.png")

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
