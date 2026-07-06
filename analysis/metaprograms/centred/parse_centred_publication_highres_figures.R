####################
# parse_centred_publication_highres_figures.R
#
# Description:
#   Generates centred high-resolution publication-style figures using the same
#   visual contracts as the uncentred publication scripts for high-resolution
#   metaprogram similarity, T2T4-vs-T0eR4 MP delta heatmaps, and timepoint MP
#   heatmaps. Centred T2/T4-high MPs are annotated by best uncentred MP gene-set
#   match and automatic top 3CA non-cell-cycle label.
#
# Inputs:
#   parse_outs/centred/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_geneNMF_metaprograms_nMP<half_total>.rds
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_UCell_scores_nMP<half_total>.rds
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_selected_mp_genes_nMP<half_total>.rds
#   parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds
#
# Outputs:
#   parse_outs/centred/publication/highres_metaprogram_heatmap/figures/highres_metaprogram_heatmap.{pdf,png}
#   parse_outs/centred/publication/t2t4_vs_t0er4_mp_delta/figures/t2t4_vs_t0er4_parse_mps_delta_heatmap.{pdf,png}
#   parse_outs/centred/publication/timepoint_mps_heatmap/figures/timepoint_parse_mps_heatmap.{pdf,png}
#   parse_outs/centred/publication/t2t4_best_match/tables/centred_T2T4_best_match_annotations.csv
#   parse_outs/logs/run_summaries/parse_centred_publication_highres_figures_*.txt
#
# Cache / replot:
#   Reuses centred GeneNMF, centred T2/T4-high UCell scores, selected genes,
#   and state assignments. Re-running overwrites only centred publication
#   figures and best-match annotation tables.
#
# Methodology:
#   analysis/methodology/publication/highres_metaprogram_heatmap_methodology.md
#   analysis/methodology/metaprograms/legacy_highres_mp_t2t4_comparison_filter_methodology.md
#
# Downstream status:
#   Terminal centred publication-style figure workflow.
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_centred_publication_highres_figures",
  parameters = list(method = "centred", nMP_rule = "total_nmf_programmes / 2"),
  input_files = c(
    "parse_outs/centred/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds",
    "parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_geneNMF_metaprograms_nMP<half_total>.rds",
    "parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/*",
    "parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds"
  ),
  output_files = c(
    "parse_outs/centred/publication/highres_metaprogram_heatmap/figures/highres_metaprogram_heatmap.pdf",
    "parse_outs/centred/publication/t2t4_vs_t0er4_mp_delta/figures/t2t4_vs_t0er4_parse_mps_delta_heatmap.pdf",
    "parse_outs/centred/publication/timepoint_mps_heatmap/figures/timepoint_parse_mps_heatmap.pdf"
  )
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(viridis)
  library(RColorBrewer)
})

project_dir <- parse_project_root()
paths <- parse_paths(project_dir)
setwd(paths$parse_outs)

centred_mp_dir <- file.path(paths$parse_outs, "centred", "Auto_parse_metaprograms")
centred_highres_dir <- file.path(paths$parse_outs, "centred", "Auto_parse_highres_metaprogram_trends")
centred_t2t4_dir <- file.path(centred_highres_dir, "Auto_T2T4_gt_T0eR4_filter")
best_match_dir <- file.path(paths$parse_outs, "centred", "publication", "t2t4_best_match")
best_match_table_dir <- file.path(best_match_dir, "tables")
dir.create(best_match_table_dir, recursive = TRUE, showWarnings = FALSE)

find_existing_dir <- function(candidates) {
  found <- candidates[dir.exists(candidates)]
  if (length(found) == 0) {
    stop("Could not find any of: ", paste(candidates, collapse = ", "))
  }
  found[1]
}

uncentred_highres_dir <- find_existing_dir(c(
  file.path(paths$parse_outs, "Auto_parse_highres_metaprogram_trends"),
  file.path(paths$parse_outs, "parse_outs", "Auto_parse_highres_metaprogram_trends")
))
uncentred_t2t4_dir <- file.path(uncentred_highres_dir, "Auto_T2T4_gt_T0eR4_filter")

get_nmp_from_geneNMF <- function(geneNMF_out_path) {
  if (!file.exists(geneNMF_out_path)) stop("Missing GeneNMF programmes: ", geneNMF_out_path)
  geneNMF.programs <- readRDS(geneNMF_out_path)
  total_nmf_programs <- sum(vapply(geneNMF.programs, function(x) ncol(x$w), numeric(1)))
  nMP <- as.integer(total_nmf_programs / 2)
  if (total_nmf_programs / 2 != nMP) stop("Total NMF programme count is odd: ", total_nmf_programs)
  nMP
}

nMP <- get_nmp_from_geneNMF(file.path(centred_mp_dir, "Auto_parse_geneNMF_outs.rds"))

jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  union_n <- length(union(a, b))
  if (union_n == 0) return(NA_real_)
  length(intersect(a, b)) / union_n
}

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
  "Proliferation / DNA Replication" = c("MP3", "MP27", "MP29"),
  "Epithelial Differentiation" = c("MP47")
)
all_clusters <- c(increase_clusters, decrease_clusters)
cluster_direction <- c(
  rep("Enriched in Persisters", length(increase_clusters)),
  rep("Depleted in Persisters", length(decrease_clusters))
)
names(cluster_direction) <- names(all_clusters)

build_uncentred_cluster_map <- function() {
  data.frame(
    uncentred_MP = unname(unlist(all_clusters)),
    Cluster = rep(names(all_clusters), lengths(all_clusters)),
    Direction = rep(cluster_direction[names(all_clusters)], lengths(all_clusters)),
    stringsAsFactors = FALSE
  ) |>
    dplyr::distinct(uncentred_MP, .keep_all = TRUE)
}

make_best_match_annotations <- function() {
  centred_genes_path <- file.path(centred_t2t4_dir, paste0("Auto_parse_highres_T2T4_selected_mp_genes_nMP", nMP, ".rds"))
  uncentred_genes_path <- file.path(uncentred_t2t4_dir, paste0("Auto_parse_highres_T2T4_selected_mp_genes_nMP", nMP, ".rds"))
  centred_label_path <- file.path(centred_t2t4_dir, paste0("Auto_parse_highres_T2T4_top_3CA_noncellcycle_nMP", nMP, ".csv"))
  centred_trend_path <- file.path(centred_t2t4_dir, paste0("Auto_parse_highres_T2T4_filter_summary_nMP", nMP, ".csv"))
  missing <- c(centred_genes_path, uncentred_genes_path, centred_label_path, centred_trend_path)[
    !file.exists(c(centred_genes_path, uncentred_genes_path, centred_label_path, centred_trend_path))
  ]
  if (length(missing) > 0) stop("Missing T2/T4 best-match inputs: ", paste(missing, collapse = ", "))

  centred_genes <- readRDS(centred_genes_path)
  uncentred_genes <- readRDS(uncentred_genes_path)
  label_table <- read.csv(centred_label_path, check.names = FALSE, stringsAsFactors = FALSE)
  trend_summary <- read.csv(centred_trend_path, check.names = FALSE, stringsAsFactors = FALSE)
  retained <- trend_summary$MP[trend_summary$retained == TRUE | trend_summary$retained == "TRUE"]
  centred_genes <- centred_genes[intersect(retained, names(centred_genes))]

  pairwise <- expand.grid(
    centred_MP = names(centred_genes),
    uncentred_MP = names(uncentred_genes),
    stringsAsFactors = FALSE
  ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      overlap_genes = length(intersect(centred_genes[[centred_MP]], uncentred_genes[[uncentred_MP]])),
      centred_gene_count = length(centred_genes[[centred_MP]]),
      uncentred_gene_count = length(uncentred_genes[[uncentred_MP]]),
      jaccard = jaccard(centred_genes[[centred_MP]], uncentred_genes[[uncentred_MP]])
    ) |>
    dplyr::ungroup()

  cluster_map <- build_uncentred_cluster_map()
  best <- pairwise |>
    dplyr::arrange(centred_MP, dplyr::desc(jaccard), dplyr::desc(overlap_genes), uncentred_MP) |>
    dplyr::group_by(centred_MP) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::left_join(cluster_map, by = "uncentred_MP") |>
    dplyr::left_join(label_table[, c("MP", "top_3ca_noncc", "top_3ca_noncc_p_adj")], by = c("centred_MP" = "MP")) |>
    dplyr::mutate(
      top_3ca_noncc = ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc),
      Cluster = ifelse(is.na(Cluster) | Cluster == "", paste0("Best 3CA: ", top_3ca_noncc), Cluster),
      Direction = ifelse(is.na(Direction) | Direction == "", "Enriched in Persisters", Direction),
      display_label = paste0(centred_MP, "\nbest:", uncentred_MP),
      annotation_label = paste0(centred_MP, " best:", uncentred_MP, " ", top_3ca_noncc)
    )
  write.csv(best, file.path(best_match_table_dir, "centred_T2T4_best_match_annotations.csv"), row.names = FALSE)
  best
}

best_match <- make_best_match_annotations()

plot_similarity_heatmap <- function() {
  out_base <- file.path(paths$parse_outs, "centred", "publication", "highres_metaprogram_heatmap")
  fig_dir <- file.path(out_base, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  mp_path <- file.path(centred_highres_dir, paste0("Auto_parse_highres_geneNMF_metaprograms_nMP", nMP, ".rds"))
  if (!file.exists(mp_path)) stop("Missing centred high-resolution metaprograms: ", mp_path)
  mp.res <- readRDS(mp_path)

  J <- mp.res[["programs.similarity"]]
  tree <- mp.res[["programs.tree"]]
  cl_members <- mp.res[["programs.clusters"]]
  J[J < 0] <- 0
  J[J > 1] <- 1

  labs.order <- labels(as.dendrogram(tree))
  cl_members_labels <- paste0("MP", cl_members)
  names(cl_members_labels) <- names(cl_members)
  row_split_factor <- factor(cl_members_labels, levels = unique(cl_members_labels[labs.order]))

  col_fun <- colorRamp2(seq(0, 1, length = 100), viridis(100, option = "A", direction = -1))
  anno_colors <- brewer.pal(n = max(min(nMP, 12), 3), name = "Paired")
  anno_colors_full <- rep(anno_colors, length.out = nMP)
  names(anno_colors_full) <- levels(row_split_factor)

  row_anno <- rowAnnotation(
    MP = row_split_factor,
    col = list(MP = anno_colors_full),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    width = unit(2, "mm")
  )
  col_anno <- HeatmapAnnotation(
    MP = row_split_factor,
    col = list(MP = anno_colors_full),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    height = unit(2, "mm")
  )

  ht <- Heatmap(
    J,
    name = "Cosine Similarity",
    col = col_fun,
    cluster_rows = tree,
    cluster_columns = tree,
    row_split = nMP,
    column_split = nMP,
    row_gap = unit(0.05, "mm"),
    column_gap = unit(0.05, "mm"),
    border = gpar(col = "black", lwd = 0.2),
    left_annotation = row_anno,
    top_annotation = col_anno,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = NULL,
    use_raster = TRUE,
    raster_quality = 5,
    column_title = paste0("Centred High-Resolution Metaprograms Similarity (nMP = ", nMP, ")"),
    column_title_gp = gpar(fontsize = 16, fontface = "bold")
  )

  pdf_path <- file.path(fig_dir, "highres_metaprogram_heatmap.pdf")
  png_path <- file.path(fig_dir, "highres_metaprogram_heatmap.png")
  pdf(pdf_path, width = 12, height = 11, useDingbats = FALSE)
  draw(ht)
  dev.off()
  png(png_path, width = 3600, height = 3300, res = 300)
  draw(ht)
  dev.off()
}

load_ucell_and_states <- function() {
  ucell_path <- file.path(centred_t2t4_dir, paste0("Auto_parse_highres_T2T4_UCell_scores_nMP", nMP, ".rds"))
  states_path <- file.path(paths$parse_outs, "cell_states", "Auto_parse_PDOpipeline_topmp_assignments.rds")
  if (!file.exists(ucell_path)) stop("Missing centred T2/T4 UCell scores: ", ucell_path)
  if (!file.exists(states_path)) stop("Missing state assignment file: ", states_path)
  ucell <- readRDS(ucell_path)
  rownames(ucell) <- sub("^([A-Za-z0-9]+)_", "\\1__", rownames(ucell))
  states <- readRDS(states_path)
  list(ucell = ucell, states = states)
}

state_levels <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation"
)
state_labels_split <- c(
  "Classic Proliferative" = "Classic\nProliferative",
  "Basal to Intest. Meta" = "Basal to\nIntest. Meta",
  "SMG-like Metaplasia" = "SMG-like\nMetaplasia",
  "Stress-adaptive" = "Stress-\nadaptive",
  "3CA_EMT_and_Protein_maturation" = "3CA EMT &\nProt. mat."
)
state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8"
)

plot_delta_heatmap <- function() {
  out_base <- file.path(paths$parse_outs, "centred", "publication", "t2t4_vs_t0er4_mp_delta")
  fig_dir <- file.path(out_base, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  data <- load_ucell_and_states()
  ucell <- data$ucell
  states <- data$states |>
    dplyr::mutate(state = pdo_state)
  emt_prot_states <- c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III")
  states$state[states$state %in% emt_prot_states] <- "3CA_EMT_and_Protein_maturation"

  common_cells <- intersect(rownames(ucell), states$cell)
  cell_df <- states |>
    dplyr::filter(cell %in% common_cells) |>
    dplyr::mutate(Treatment = ifelse(sample %in% c("T2", "T4"), "T2T4", ifelse(sample %in% c("T0", "eR4"), "T0eR4", NA))) |>
    dplyr::filter(!is.na(Treatment), state %in% state_levels)
  active_state_levels <- intersect(state_levels, unique(cell_df$state))
  mps <- intersect(best_match$centred_MP, colnames(ucell))

  mp_long <- as.data.frame(ucell[cell_df$cell, mps, drop = FALSE]) |>
    dplyr::mutate(cell = cell_df$cell, state = cell_df$state, Treatment = cell_df$Treatment) |>
    tidyr::pivot_longer(cols = dplyr::all_of(mps), names_to = "centred_MP", values_to = "score")
  mp_mean_delta <- mp_long |>
    dplyr::group_by(state, Treatment, centred_MP) |>
    dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Treatment, values_from = mean_score) |>
    dplyr::filter(!is.na(T0eR4), !is.na(T2T4)) |>
    dplyr::mutate(mean_delta = T2T4 - T0eR4) |>
    dplyr::left_join(
      mp_long |>
        dplyr::group_by(state, centred_MP) |>
        dplyr::summarise(p_value = tryCatch(t.test(score ~ Treatment)$p.value, error = function(e) NA_real_), .groups = "drop"),
      by = c("state", "centred_MP")
    ) |>
    dplyr::mutate(sig_label = dplyr::case_when(is.na(p_value) ~ "", p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value < 0.05 ~ "*", TRUE ~ "")) |>
    dplyr::left_join(best_match, by = "centred_MP")

  create_mp_ht <- function(df, title) {
    if (nrow(df) == 0) return(NULL)
    mat <- df |>
      dplyr::select(centred_MP, state, mean_delta) |>
      dplyr::distinct() |>
      tidyr::pivot_wider(names_from = state, values_from = mean_delta) |>
      tibble::column_to_rownames("centred_MP") |>
      as.matrix()
    mat <- mat[, active_state_levels[active_state_levels %in% colnames(mat)], drop = FALSE]
    sig <- df |>
      dplyr::select(centred_MP, state, sig_label) |>
      dplyr::distinct() |>
      tidyr::pivot_wider(names_from = state, values_from = sig_label) |>
      tibble::column_to_rownames("centred_MP") |>
      as.matrix()
    sig <- sig[rownames(mat), colnames(mat), drop = FALSE]
    sig[is.na(sig)] <- ""

    df_unique <- df |>
      dplyr::select(centred_MP, Cluster, display_label) |>
      dplyr::distinct()
    cluster_lvls <- unique(names(all_clusters))
    fallback_clusters <- setdiff(unique(df_unique$Cluster), cluster_lvls)
    cluster_lvls <- c(cluster_lvls, sort(fallback_clusters))
    df_unique <- df_unique |>
      dplyr::mutate(Cluster = factor(Cluster, levels = cluster_lvls)) |>
      dplyr::arrange(Cluster, centred_MP)
    mat <- mat[df_unique$centred_MP, , drop = FALSE]
    sig <- sig[df_unique$centred_MP, , drop = FALSE]
    row_split_fac <- factor(df_unique$Cluster, levels = cluster_lvls[cluster_lvls %in% df_unique$Cluster])
    row_labels <- setNames(df_unique$display_label, df_unique$centred_MP)
    clip <- max(0.002, quantile(abs(mat), 0.95, na.rm = TRUE))
    col_f <- colorRamp2(c(-clip, 0, clip), c("#245F7B", "white", "#B63E2F"))
    cluster_colors <- rainbow(length(cluster_lvls))
    names(cluster_colors) <- cluster_lvls
    ha_row <- rowAnnotation(Cluster = df_unique$Cluster, col = list(Cluster = cluster_colors), show_annotation_name = FALSE, show_legend = FALSE)
    cell_f <- function(j, i, x, y, w, h, fill) {
      grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 7, fontfamily = "Arial"))
    }
    Heatmap(
      mat,
      name = "Delta Score",
      col = col_f,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
      row_labels = row_labels[rownames(mat)],
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
      column_title = paste("Centred MPs", title),
      column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial")
    )
  }

  ht_inc <- create_mp_ht(mp_mean_delta |> dplyr::filter(Direction == "Enriched in Persisters"), "Increased in T2/T4")
  ht_dec <- create_mp_ht(mp_mean_delta |> dplyr::filter(Direction == "Depleted in Persisters"), "Decreased in T2/T4")
  pdf_path <- file.path(fig_dir, "t2t4_vs_t0er4_parse_mps_delta_heatmap.pdf")
  png_path <- file.path(fig_dir, "t2t4_vs_t0er4_parse_mps_delta_heatmap.png")
  grDevices::cairo_pdf(pdf_path, width = 450 / 25.4, height = 200 / 25.4, family = "Arial")
  pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  if (!is.null(ht_inc)) draw(ht_inc, newpage = FALSE, heatmap_legend_side = "bottom")
  popViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  if (!is.null(ht_dec)) draw(ht_dec, newpage = FALSE, heatmap_legend_side = "bottom")
  popViewport()
  dev.off()
  png(png_path, width = 450, height = 200, units = "mm", res = 600, type = "cairo")
  pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  if (!is.null(ht_inc)) draw(ht_inc, newpage = FALSE, heatmap_legend_side = "bottom")
  popViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  if (!is.null(ht_dec)) draw(ht_dec, newpage = FALSE, heatmap_legend_side = "bottom")
  popViewport()
  dev.off()
}

plot_timepoint_heatmap <- function() {
  out_base <- file.path(paths$parse_outs, "centred", "publication", "timepoint_mps_heatmap")
  fig_dir <- file.path(out_base, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  data <- load_ucell_and_states()
  ucell <- data$ucell
  states <- data$states
  timepoint_levels <- c("T0", "T1", "T2", "T4", "R4", "eR4")
  common_cells <- intersect(rownames(ucell), states$cell)
  cell_df <- states |>
    dplyr::filter(cell %in% common_cells, sample %in% timepoint_levels) |>
    dplyr::mutate(sample = factor(sample, levels = timepoint_levels))
  mps <- intersect(best_match$centred_MP, colnames(ucell))
  mp_long <- as.data.frame(ucell[cell_df$cell, mps, drop = FALSE]) |>
    dplyr::mutate(cell = cell_df$cell, sample = cell_df$sample) |>
    tidyr::pivot_longer(cols = dplyr::all_of(mps), names_to = "centred_MP", values_to = "score")
  mp_mean_matrix <- mp_long |>
    dplyr::group_by(sample, centred_MP) |>
    dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = sample, values_from = mean_score) |>
    dplyr::left_join(best_match, by = "centred_MP")

  create_mp_ht <- function(df, title) {
    if (nrow(df) == 0) return(NULL)
    mat_raw <- df |>
      dplyr::select(centred_MP, dplyr::all_of(timepoint_levels)) |>
      dplyr::distinct() |>
      tibble::column_to_rownames("centred_MP") |>
      as.matrix()
    mat_raw <- mat_raw[, timepoint_levels[timepoint_levels %in% colnames(mat_raw)], drop = FALSE]
    mat_z <- t(scale(t(mat_raw)))
    mat_z[!is.finite(mat_z)] <- 0
    df_unique <- df |>
      dplyr::select(centred_MP, Cluster, display_label) |>
      dplyr::distinct()
    cluster_lvls <- unique(names(all_clusters))
    fallback_clusters <- setdiff(unique(df_unique$Cluster), cluster_lvls)
    cluster_lvls <- c(cluster_lvls, sort(fallback_clusters))
    df_unique <- df_unique |>
      dplyr::mutate(Cluster = factor(Cluster, levels = cluster_lvls)) |>
      dplyr::arrange(Cluster, centred_MP)
    mat_raw <- mat_raw[df_unique$centred_MP, , drop = FALSE]
    mat_z <- mat_z[df_unique$centred_MP, , drop = FALSE]
    row_split_fac <- factor(df_unique$Cluster, levels = cluster_lvls[cluster_lvls %in% df_unique$Cluster])
    row_labels <- setNames(df_unique$display_label, df_unique$centred_MP)
    col_f <- colorRamp2(c(-2, 0, 2), c("#245F7B", "white", "#B63E2F"))
    cluster_colors <- rainbow(length(cluster_lvls))
    names(cluster_colors) <- cluster_lvls
    ha_row <- rowAnnotation(Cluster = df_unique$Cluster, col = list(Cluster = cluster_colors), show_annotation_name = FALSE, show_legend = FALSE)
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
      row_labels = row_labels[rownames(mat_z)],
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
      column_title = paste("Centred MPs", title),
      column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial")
    )
  }
  ht_inc <- create_mp_ht(mp_mean_matrix |> dplyr::filter(Direction == "Enriched in Persisters"), "Enriched in Persisters")
  ht_dec <- create_mp_ht(mp_mean_matrix |> dplyr::filter(Direction == "Depleted in Persisters"), "Depleted in Persisters")
  pdf_path <- file.path(fig_dir, "timepoint_parse_mps_heatmap.pdf")
  png_path <- file.path(fig_dir, "timepoint_parse_mps_heatmap.png")
  grDevices::cairo_pdf(pdf_path, width = 450 / 25.4, height = 200 / 25.4, family = "Arial")
  pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  if (!is.null(ht_inc)) draw(ht_inc, newpage = FALSE, heatmap_legend_side = "bottom")
  popViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  if (!is.null(ht_dec)) draw(ht_dec, newpage = FALSE, heatmap_legend_side = "bottom")
  popViewport()
  dev.off()
  png(png_path, width = 450, height = 200, units = "mm", res = 600, type = "cairo")
  pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  if (!is.null(ht_inc)) draw(ht_inc, newpage = FALSE, heatmap_legend_side = "bottom")
  popViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  if (!is.null(ht_dec)) draw(ht_dec, newpage = FALSE, heatmap_legend_side = "bottom")
  popViewport()
  dev.off()
}

plot_similarity_heatmap()
plot_delta_heatmap()
plot_timepoint_heatmap()

script_run_status <- "success"
message("parse_centred_publication_highres_figures.R completed successfully.")
