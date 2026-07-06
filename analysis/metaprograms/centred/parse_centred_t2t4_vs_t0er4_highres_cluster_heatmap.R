####################
# parse_centred_t2t4_vs_t0er4_highres_cluster_heatmap.R
#
# Description:
#   Visualises centred high-resolution Parse MP score changes per PDO-pipeline
#   cell state, comparing T2+T4 vs T0+eR4. Retained centred MPs use nMP equal
#   to total six-sample Parse NMF programmes divided by 2, and are grouped by
#   automatic top 3CA non-cell-cycle annotation.
#
# Inputs:
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_UCell_scores_nMP<half_total>.rds
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_trend_summary_nMP<half_total>.csv
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_top_3CA_noncellcycle_nMP<half_total>.csv
#   parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds
#
# Outputs:
#   parse_outs/centred/cell_states/t2t4_vs_t0er4_highres_clusters/tables/parse_centred_t2t4_vs_t0er4_mp_delta_scores.csv
#   parse_outs/centred/cell_states/t2t4_vs_t0er4_highres_clusters/tables/parse_centred_t2t4_vs_t0er4_annotation_cluster_delta_scores.csv
#   parse_outs/centred/cell_states/t2t4_vs_t0er4_highres_clusters/figures/parse_centred_t2t4_vs_t0er4_MP_delta_heatmap.pdf
#   parse_outs/centred/cell_states/t2t4_vs_t0er4_highres_clusters/figures/parse_centred_t2t4_vs_t0er4_annotation_cluster_delta_heatmap.pdf
#   parse_outs/logs/run_summaries/parse_centred_t2t4_vs_t0er4_highres_cluster_heatmap_*.txt
#
# Cache / replot:
#   Reuses centred high-resolution UCell scores and annotation tables produced
#   by parse_centred_highres_mp_strict_mean_median_trend_filter.R.
#
# Methodology:
#   analysis/methodology/metaprograms/legacy_highres_mp_t2t4_comparison_filter_methodology.md
#
# Downstream status:
#   Terminal centred comparison figure and table workflow.
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_centred_t2t4_vs_t0er4_highres_cluster_heatmap",
  parameters = list(contrast = "T2+T4 vs T0+eR4", grouping = "top_3ca_noncc"),
  input_files = c(
    "parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_UCell_scores_nMP<half_total>.rds",
    "parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_trend_summary_nMP<half_total>.csv",
    "parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_top_3CA_noncellcycle_nMP<half_total>.csv",
    "parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds"
  ),
  output_files = c(
    "parse_outs/centred/cell_states/t2t4_vs_t0er4_highres_clusters/tables/parse_centred_t2t4_vs_t0er4_mp_delta_scores.csv",
    "parse_outs/centred/cell_states/t2t4_vs_t0er4_highres_clusters/figures/parse_centred_t2t4_vs_t0er4_MP_delta_heatmap.pdf"
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
})

project_dir <- parse_project_root()
paths <- parse_paths(project_dir)
highres_dir <- file.path(paths$parse_outs, "centred", "Auto_parse_highres_metaprogram_trends")
out_dir <- file.path(paths$parse_outs, "centred", "cell_states", "t2t4_vs_t0er4_highres_clusters")
fig_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

config_files <- list.files(highres_dir, pattern = "^Auto_parse_highres_nMP[0-9]+_config\\.csv$", full.names = TRUE)
if (length(config_files) == 0) {
  stop("No centred high-resolution config found in ", highres_dir)
}
config <- read.csv(config_files[order(file.info(config_files)$mtime, decreasing = TRUE)[1]], check.names = FALSE)
nMP <- as.integer(config$nMP[1])

ucell_path <- file.path(highres_dir, paste0("Auto_parse_highres_UCell_scores_nMP", nMP, ".rds"))
trend_path <- file.path(highres_dir, paste0("Auto_parse_highres_trend_summary_nMP", nMP, ".csv"))
label_path <- file.path(highres_dir, paste0("Auto_parse_highres_top_3CA_noncellcycle_nMP", nMP, ".csv"))
states_path <- file.path(paths$parse_outs, "cell_states", "Auto_parse_PDOpipeline_topmp_assignments.rds")
missing_inputs <- c(ucell_path, trend_path, label_path, states_path)[!file.exists(c(ucell_path, trend_path, label_path, states_path))]
if (length(missing_inputs) > 0) {
  stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "))
}

ucell <- readRDS(ucell_path)
trend_summary <- read.csv(trend_path, check.names = FALSE, stringsAsFactors = FALSE)
label_table <- read.csv(label_path, check.names = FALSE, stringsAsFactors = FALSE)
states <- readRDS(states_path)

retained_mps <- trend_summary$MP[trend_summary$retained == TRUE | trend_summary$retained == "TRUE"]
retained_mps <- intersect(retained_mps, colnames(ucell))
if (length(retained_mps) == 0) {
  stop("No retained centred MPs are present in the UCell matrix.")
}

label_table <- label_table |>
  dplyr::mutate(
    top_3ca_noncc = ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc)
  )
mp_annotation <- trend_summary |>
  dplyr::filter(MP %in% retained_mps) |>
  dplyr::left_join(label_table[, c("MP", "top_3ca_noncc")], by = "MP") |>
  dplyr::mutate(
    top_3ca_noncc = ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc),
    trend_type_label = ifelse(is.na(trend_type_label) | trend_type_label == "", treatment_direction, trend_type_label),
    display_label = paste0(MP, "\n", top_3ca_noncc)
  )

normalise_cell_names <- function(x, state_cells) {
  common <- intersect(rownames(x), state_cells)
  if (length(common) > 0) return(x)
  rownames(x) <- sub("^([A-Za-z0-9]+)_", "\\1__", rownames(x))
  x
}

states <- states |>
  dplyr::mutate(state = pdo_state)
emt_prot_states <- c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III")
states$state[states$state %in% emt_prot_states] <- "3CA_EMT_and_Protein_maturation"
ucell <- normalise_cell_names(ucell, states$cell)

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

cell_df <- states |>
  dplyr::filter(cell %in% rownames(ucell)) |>
  dplyr::mutate(
    Treatment = dplyr::case_when(
      sample %in% c("T2", "T4") ~ "T2T4",
      sample %in% c("T0", "eR4") ~ "T0eR4",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(Treatment), state %in% state_levels)
state_levels <- intersect(state_levels, unique(cell_df$state))
if (nrow(cell_df) == 0 || length(state_levels) == 0) {
  stop("No cells overlap between centred UCell scores and PDO-pipeline state assignments.")
}

mp_long <- as.data.frame(ucell[cell_df$cell, retained_mps, drop = FALSE]) |>
  tibble::rownames_to_column("cell") |>
  dplyr::left_join(cell_df[, c("cell", "state", "Treatment")], by = "cell") |>
  tidyr::pivot_longer(cols = dplyr::all_of(retained_mps), names_to = "MP", values_to = "score")

mp_delta <- mp_long |>
  dplyr::group_by(state, Treatment, MP) |>
  dplyr::summarise(mean_score = mean(score, na.rm = TRUE), n_cells = dplyr::n(), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Treatment, values_from = c(mean_score, n_cells)) |>
  dplyr::filter(!is.na(mean_score_T0eR4), !is.na(mean_score_T2T4)) |>
  dplyr::mutate(mean_delta = mean_score_T2T4 - mean_score_T0eR4) |>
  dplyr::left_join(
    mp_long |>
      dplyr::group_by(state, MP) |>
      dplyr::summarise(p_value = tryCatch(stats::t.test(score ~ Treatment)$p.value, error = function(e) NA_real_), .groups = "drop"),
    by = c("state", "MP")
  ) |>
  dplyr::mutate(
    sig_label = dplyr::case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) |>
  dplyr::left_join(mp_annotation[, c("MP", "top_3ca_noncc", "trend_type_label", "display_label")], by = "MP")
write.csv(mp_delta, file.path(table_dir, "parse_centred_t2t4_vs_t0er4_mp_delta_scores.csv"), row.names = FALSE)

cluster_long <- mp_long |>
  dplyr::left_join(mp_annotation[, c("MP", "top_3ca_noncc", "trend_type_label")], by = "MP") |>
  dplyr::group_by(cell, state, Treatment, top_3ca_noncc, trend_type_label) |>
  dplyr::summarise(score = mean(score, na.rm = TRUE), .groups = "drop")

cluster_delta <- cluster_long |>
  dplyr::group_by(state, Treatment, top_3ca_noncc, trend_type_label) |>
  dplyr::summarise(mean_score = mean(score, na.rm = TRUE), n_cells = dplyr::n(), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Treatment, values_from = c(mean_score, n_cells)) |>
  dplyr::filter(!is.na(mean_score_T0eR4), !is.na(mean_score_T2T4)) |>
  dplyr::mutate(
    mean_delta = mean_score_T2T4 - mean_score_T0eR4,
    annotation_trend_id = paste(top_3ca_noncc, trend_type_label, sep = " | ")
  )
write.csv(cluster_delta, file.path(table_dir, "parse_centred_t2t4_vs_t0er4_annotation_cluster_delta_scores.csv"), row.names = FALSE)

make_delta_heatmap <- function(delta_df, row_col, label_col, output_prefix, row_height_cm = 0.75) {
  row_order <- delta_df |>
    dplyr::distinct(.data[[row_col]], .data[[label_col]]) |>
    dplyr::arrange(.data[[label_col]], .data[[row_col]])
  mat <- delta_df |>
    dplyr::select(row_id = dplyr::all_of(row_col), state, mean_delta) |>
    dplyr::distinct() |>
    tidyr::pivot_wider(names_from = state, values_from = mean_delta) |>
    tibble::column_to_rownames("row_id") |>
    as.matrix()
  mat <- mat[row_order[[row_col]], state_levels[state_levels %in% colnames(mat)], drop = FALSE]
  sig <- if ("sig_label" %in% colnames(delta_df)) {
    delta_df |>
      dplyr::select(row_id = dplyr::all_of(row_col), state, sig_label) |>
      dplyr::distinct() |>
      tidyr::pivot_wider(names_from = state, values_from = sig_label) |>
      tibble::column_to_rownames("row_id") |>
      as.matrix()
  } else {
    matrix("", nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  }
  sig <- sig[rownames(mat), colnames(mat), drop = FALSE]
  sig[is.na(sig)] <- ""

  split_values <- factor(row_order[[label_col]], levels = unique(row_order[[label_col]]))
  names(split_values) <- row_order[[row_col]]
  split_values <- split_values[rownames(mat)]

  clip <- max(0.002, stats::quantile(abs(mat), 0.95, na.rm = TRUE))
  col_fun <- circlize::colorRamp2(c(-clip, 0, clip), c("#245F7B", "white", "#B63E2F"))
  cell_fun <- function(j, i, x, y, w, h, fill) {
    grid::grid.text(sprintf("%.4f", mat[i, j]), x, y - grid::unit(1, "mm"), gp = grid::gpar(fontsize = 6.5))
    if (!is.na(sig[i, j]) && sig[i, j] != "") {
      grid::grid.text(sig[i, j], x, y + grid::unit(2.5, "mm"), gp = grid::gpar(fontsize = 9, fontface = "bold"))
    }
  }
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "Mean delta\nT2T4 - T0eR4",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = split_values,
    row_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    column_labels = state_labels_split[colnames(mat)],
    column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    column_names_rot = 0,
    column_names_centered = TRUE,
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      State = colnames(mat),
      col = list(State = state_cols[colnames(mat)]),
      show_annotation_name = FALSE
    ),
    cell_fun = cell_fun,
    width = grid::unit(ncol(mat) * 2.4, "cm"),
    height = grid::unit(max(1, nrow(mat)) * row_height_cm, "cm")
  )
  pdf_path <- file.path(fig_dir, paste0(output_prefix, ".pdf"))
  png_path <- file.path(fig_dir, paste0(output_prefix, ".png"))
  grDevices::pdf(pdf_path, width = 12, height = max(6, 0.35 * nrow(mat) + 3), useDingbats = FALSE)
  ComplexHeatmap::draw(ht, merge_legend = TRUE)
  grDevices::dev.off()
  grDevices::png(png_path, width = 12, height = max(6, 0.35 * nrow(mat) + 3), units = "in", res = 300)
  ComplexHeatmap::draw(ht, merge_legend = TRUE)
  grDevices::dev.off()
}

make_delta_heatmap(
  mp_delta,
  row_col = "MP",
  label_col = "top_3ca_noncc",
  output_prefix = "parse_centred_t2t4_vs_t0er4_MP_delta_heatmap",
  row_height_cm = 0.72
)
make_delta_heatmap(
  cluster_delta,
  row_col = "annotation_trend_id",
  label_col = "trend_type_label",
  output_prefix = "parse_centred_t2t4_vs_t0er4_annotation_cluster_delta_heatmap",
  row_height_cm = 0.9
)

script_run_status <- "success"
message("parse_centred_t2t4_vs_t0er4_highres_cluster_heatmap.R completed successfully.")
