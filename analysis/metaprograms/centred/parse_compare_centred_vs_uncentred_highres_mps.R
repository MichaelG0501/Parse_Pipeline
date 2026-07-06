####################
# parse_compare_centred_vs_uncentred_highres_mps.R
#
# Description:
#   Compares final retained centred high-resolution MPs against the existing
#   uncentred high-resolution MPs by gene-set Jaccard overlap and automatic top
#   3CA non-cell-cycle annotation.
#
# Inputs:
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_selected_mp_genes_nMP117.rds
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_trend_summary_nMP117.csv
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_top_3CA_noncellcycle_nMP117.csv
#   parse_outs/Auto_parse_highres_metaprogram_trends/ or parse_outs/parse_outs/Auto_parse_highres_metaprogram_trends/ matching uncentred files
#
# Outputs:
#   parse_outs/centred/comparison/tables/centred_vs_uncentred_highres_best_matches.csv
#   parse_outs/centred/comparison/tables/uncentred_vs_centred_highres_best_matches.csv
#   parse_outs/centred/comparison/tables/centred_vs_uncentred_annotation_summary.csv
#   parse_outs/centred/comparison/figures/centred_vs_uncentred_jaccard_heatmap.pdf
#   parse_outs/centred/comparison/reports/centred_vs_uncentred_highres_summary.txt
#   parse_outs/logs/run_summaries/parse_compare_centred_vs_uncentred_highres_mps_*.txt
#
# Cache / replot:
#   Reads existing centred and uncentred high-resolution outputs and overwrites
#   comparison tables/figures only.
#
# Methodology:
#   analysis/methodology/metaprograms/highres_mp_strict_mean_median_trend_filter_methodology.md
#
# Downstream status:
#   Terminal comparison report. Does not feed canonical Parse MP definitions.
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_compare_centred_vs_uncentred_highres_mps",
  parameters = list(match_metric = "Jaccard gene overlap", annotation = "top_3ca_noncc"),
  input_files = c(
    "parse_outs/centred/Auto_parse_highres_metaprogram_trends/*nMP117*",
    "parse_outs/Auto_parse_highres_metaprogram_trends/*nMP117* or parse_outs/parse_outs/Auto_parse_highres_metaprogram_trends/*nMP117*"
  ),
  output_files = c(
    "parse_outs/centred/comparison/tables/centred_vs_uncentred_highres_best_matches.csv",
    "parse_outs/centred/comparison/tables/centred_vs_uncentred_annotation_summary.csv",
    "parse_outs/centred/comparison/figures/centred_vs_uncentred_jaccard_heatmap.pdf"
  )
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

project_dir <- parse_project_root()
paths <- parse_paths(project_dir)
centred_dir <- file.path(paths$parse_outs, "centred", "Auto_parse_highres_metaprogram_trends")
uncentred_candidates <- c(
  file.path(paths$parse_outs, "Auto_parse_highres_metaprogram_trends"),
  file.path(paths$parse_outs, "parse_outs", "Auto_parse_highres_metaprogram_trends")
)
uncentred_dir <- uncentred_candidates[dir.exists(uncentred_candidates)][1]
if (is.na(uncentred_dir)) {
  stop("Could not locate uncentred high-resolution output directory. Checked: ", paste(uncentred_candidates, collapse = ", "))
}

out_dir <- file.path(paths$parse_outs, "centred", "comparison")
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
report_dir <- file.path(out_dir, "reports")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

read_nmp <- function(base_dir) {
  config_files <- list.files(base_dir, pattern = "^Auto_parse_highres_nMP[0-9]+_config\\.csv$", full.names = TRUE)
  if (length(config_files) == 0) {
    stop("Missing high-resolution config in ", base_dir)
  }
  config <- read.csv(config_files[order(file.info(config_files)$mtime, decreasing = TRUE)[1]], check.names = FALSE)
  as.integer(config$nMP[1])
}

load_method <- function(base_dir, method_name) {
  nMP <- read_nmp(base_dir)
  genes_path <- file.path(base_dir, paste0("Auto_parse_highres_selected_mp_genes_nMP", nMP, ".rds"))
  trend_path <- file.path(base_dir, paste0("Auto_parse_highres_trend_summary_nMP", nMP, ".csv"))
  label_path <- file.path(base_dir, paste0("Auto_parse_highres_top_3CA_noncellcycle_nMP", nMP, ".csv"))
  missing <- c(genes_path, trend_path, label_path)[!file.exists(c(genes_path, trend_path, label_path))]
  if (length(missing) > 0) {
    stop("Missing ", method_name, " comparison inputs: ", paste(missing, collapse = ", "))
  }
  genes <- readRDS(genes_path)
  trend <- read.csv(trend_path, check.names = FALSE, stringsAsFactors = FALSE)
  labels <- read.csv(label_path, check.names = FALSE, stringsAsFactors = FALSE)
  labels <- labels |>
    dplyr::mutate(top_3ca_noncc = ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc))
  retained <- trend$MP[trend$retained == TRUE | trend$retained == "TRUE"]
  genes <- genes[intersect(retained, names(genes))]
  meta <- trend |>
    dplyr::filter(MP %in% names(genes)) |>
    dplyr::left_join(labels[, c("MP", "top_3ca_noncc")], by = "MP") |>
    dplyr::mutate(method = method_name)
  list(nMP = nMP, genes = genes, meta = meta, dir = base_dir)
}

centred <- load_method(centred_dir, "centred")
uncentred <- load_method(uncentred_dir, "uncentred")

jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  union_n <- length(union(a, b))
  if (union_n == 0) return(NA_real_)
  length(intersect(a, b)) / union_n
}

pairwise <- expand.grid(
  centred_MP = names(centred$genes),
  uncentred_MP = names(uncentred$genes),
  stringsAsFactors = FALSE
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    overlap_genes = length(intersect(centred$genes[[centred_MP]], uncentred$genes[[uncentred_MP]])),
    centred_gene_count = length(centred$genes[[centred_MP]]),
    uncentred_gene_count = length(uncentred$genes[[uncentred_MP]]),
    union_genes = length(union(centred$genes[[centred_MP]], uncentred$genes[[uncentred_MP]])),
    jaccard = jaccard(centred$genes[[centred_MP]], uncentred$genes[[uncentred_MP]])
  ) |>
  dplyr::ungroup()

centred_meta <- centred$meta |>
  dplyr::select(
    centred_MP = MP,
    centred_trend_type = trend_type,
    centred_trend_type_label = trend_type_label,
    centred_direction = treatment_direction,
    centred_top_3ca_noncc = top_3ca_noncc,
    centred_numberPrograms = numberPrograms,
    centred_silhouette = silhouette
  )
uncentred_meta <- uncentred$meta |>
  dplyr::select(
    uncentred_MP = MP,
    uncentred_trend_type = trend_type,
    uncentred_trend_type_label = trend_type_label,
    uncentred_direction = treatment_direction,
    uncentred_top_3ca_noncc = top_3ca_noncc,
    uncentred_numberPrograms = numberPrograms,
    uncentred_silhouette = silhouette
  )

best_centred <- pairwise |>
  dplyr::arrange(centred_MP, dplyr::desc(jaccard), dplyr::desc(overlap_genes), uncentred_MP) |>
  dplyr::group_by(centred_MP) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup() |>
  dplyr::left_join(centred_meta, by = "centred_MP") |>
  dplyr::left_join(uncentred_meta, by = "uncentred_MP") |>
  dplyr::mutate(
    same_top_3ca_noncc = centred_top_3ca_noncc == uncentred_top_3ca_noncc,
    same_trend_type = centred_trend_type == uncentred_trend_type
  )
write.csv(best_centred, file.path(table_dir, "centred_vs_uncentred_highres_best_matches.csv"), row.names = FALSE)

best_uncentred <- pairwise |>
  dplyr::arrange(uncentred_MP, dplyr::desc(jaccard), dplyr::desc(overlap_genes), centred_MP) |>
  dplyr::group_by(uncentred_MP) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup() |>
  dplyr::left_join(centred_meta, by = "centred_MP") |>
  dplyr::left_join(uncentred_meta, by = "uncentred_MP") |>
  dplyr::mutate(
    same_top_3ca_noncc = centred_top_3ca_noncc == uncentred_top_3ca_noncc,
    same_trend_type = centred_trend_type == uncentred_trend_type
  )
write.csv(best_uncentred, file.path(table_dir, "uncentred_vs_centred_highres_best_matches.csv"), row.names = FALSE)

annotation_summary <- dplyr::bind_rows(
  centred$meta |> dplyr::transmute(method = "centred", MP, top_3ca_noncc, trend_type, trend_type_label),
  uncentred$meta |> dplyr::transmute(method = "uncentred", MP, top_3ca_noncc, trend_type, trend_type_label)
) |>
  dplyr::count(method, top_3ca_noncc, trend_type_label, name = "n_mps") |>
  dplyr::arrange(top_3ca_noncc, method, trend_type_label)
write.csv(annotation_summary, file.path(table_dir, "centred_vs_uncentred_annotation_summary.csv"), row.names = FALSE)

method_annotation_sets <- annotation_summary |>
  dplyr::group_by(method) |>
  dplyr::summarise(labels = list(unique(top_3ca_noncc)), .groups = "drop")
centred_labels <- method_annotation_sets$labels[[which(method_annotation_sets$method == "centred")]]
uncentred_labels <- method_annotation_sets$labels[[which(method_annotation_sets$method == "uncentred")]]
label_overlap <- data.frame(
  common_labels = length(intersect(centred_labels, uncentred_labels)),
  centred_only_labels = length(setdiff(centred_labels, uncentred_labels)),
  uncentred_only_labels = length(setdiff(uncentred_labels, centred_labels)),
  stringsAsFactors = FALSE
)
write.csv(label_overlap, file.path(table_dir, "centred_vs_uncentred_label_overlap_summary.csv"), row.names = FALSE)

all_pairwise_path <- file.path(table_dir, "centred_vs_uncentred_all_pairwise_jaccard.csv")
write.csv(pairwise, all_pairwise_path, row.names = FALSE)

####################
make_mp_annotation_label <- function(mp, top_3ca_noncc) {
  top_3ca_noncc <- ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc)
  paste0(mp, " | ", top_3ca_noncc)
}

centred_axis <- centred_meta |>
  dplyr::mutate(axis_label = make_mp_annotation_label(centred_MP, centred_top_3ca_noncc))
uncentred_axis <- uncentred_meta |>
  dplyr::mutate(axis_label = make_mp_annotation_label(uncentred_MP, uncentred_top_3ca_noncc))

centred_order <- unique(best_centred$centred_MP[order(best_centred$uncentred_MP, best_centred$centred_MP)])
uncentred_order <- unique(best_uncentred$uncentred_MP[order(best_uncentred$uncentred_MP)])
centred_label_map <- stats::setNames(centred_axis$axis_label, centred_axis$centred_MP)
uncentred_label_map <- stats::setNames(uncentred_axis$axis_label, uncentred_axis$uncentred_MP)

plot_df <- pairwise |>
  dplyr::mutate(
    centred_MP = factor(centred_MP, levels = centred_order, labels = centred_label_map[centred_order]),
    uncentred_MP = factor(uncentred_MP, levels = rev(uncentred_order), labels = uncentred_label_map[rev(uncentred_order)])
  )
heatmap_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = centred_MP, y = uncentred_MP, fill = jaccard)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.08) +
  ggplot2::scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, max(plot_df$jaccard, na.rm = TRUE))) +
  ggplot2::labs(x = "Centred retained MP", y = "Uncentred retained MP", fill = "Jaccard") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1, size = 7.4),
    axis.text.y = ggplot2::element_text(size = 7.4),
    axis.title = ggplot2::element_text(size = 13),
    legend.title = ggplot2::element_text(size = 11),
    legend.text = ggplot2::element_text(size = 10),
    axis.line = ggplot2::element_blank()
  )
ggplot2::ggsave(file.path(fig_dir, "centred_vs_uncentred_jaccard_heatmap.pdf"), heatmap_plot, width = 12, height = 10, useDingbats = FALSE)
ggplot2::ggsave(file.path(fig_dir, "centred_vs_uncentred_jaccard_heatmap.png"), heatmap_plot, width = 12, height = 10, dpi = 300)
####################

summary_lines <- c(
  "Centred vs uncentred high-resolution MP comparison",
  paste0("centred_dir: ", centred$dir),
  paste0("uncentred_dir: ", uncentred$dir),
  paste0("centred_nMP: ", centred$nMP),
  paste0("uncentred_nMP: ", uncentred$nMP),
  paste0("centred_retained_mps: ", length(centred$genes)),
  paste0("uncentred_retained_mps: ", length(uncentred$genes)),
  paste0("centred_best_matches_with_jaccard_ge_0.25: ", sum(best_centred$jaccard >= 0.25, na.rm = TRUE)),
  paste0("centred_best_matches_with_same_top_3ca_noncc: ", sum(best_centred$same_top_3ca_noncc, na.rm = TRUE)),
  paste0("common_top_3ca_noncc_labels: ", label_overlap$common_labels[1]),
  paste0("centred_only_top_3ca_noncc_labels: ", label_overlap$centred_only_labels[1]),
  paste0("uncentred_only_top_3ca_noncc_labels: ", label_overlap$uncentred_only_labels[1])
)
writeLines(summary_lines, file.path(report_dir, "centred_vs_uncentred_highres_summary.txt"))

script_run_status <- "success"
message("parse_compare_centred_vs_uncentred_highres_mps.R completed successfully.")
