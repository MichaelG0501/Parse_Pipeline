####################
# legacy_parse_highres_mp_t2t4_comparison_filter.R
#
# LEGACY COMPARISON SCRIPT.
#
# Description:
#   Alternative high-resolution Parse MP filter retained for comparison only.
#   It uses the same nMP117 metaprograms/UCell score logic as the active strict
#   mean+median trend script, but retains MPs where mean UCell activity in both
#   T2 and T4 is higher than both T0 and eR4.
#
# Inputs:
#   parse_outs/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds
#   parse_outs/by_samples/<sample>/Auto_<sample>_final.rds
#   3CA, cell-cycle, and developmental enrichment references.
#
# Outputs:
#   parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/*
#   parse_outs/logs/run_summaries/legacy_parse_highres_mp_t2t4_comparison_filter_*.txt
#
# Downstream status:
#   No active downstream script should consume outputs from this legacy folder.
#   Use parse_highres_mp_strict_mean_median_trend_filter.R for current results.
#
# Cache / replot:
#   Supports --mode=score, --mode=enrich, --mode=excel, and --mode=all. Cached
#   nMP117 metaprograms, UCell scores, trend tables, and selected genes are
#   reused if present.
#
# Methodology:
#   analysis/methodology/metaprograms/legacy_highres_mp_t2t4_comparison_filter_methodology.md
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

args <- commandArgs(trailingOnly = TRUE)
mode_arg <- grep("^--mode=", args, value = TRUE)
mode <- if (length(mode_arg) > 0) sub("^--mode=", "", mode_arg[1]) else "all"
mode <- match.arg(mode, c("all", "score", "enrich", "excel"))

project_dir <- parse_project_root()
qc_dir <- file.path(project_dir, "parse_outs")
parse_mp_dir <- file.path(qc_dir, "Auto_parse_metaprograms")
base_highres_dir <- file.path(qc_dir, "Auto_parse_highres_metaprogram_trends")
out_dir <- file.path(base_highres_dir, "Auto_T2T4_gt_T0eR4_filter")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

parse_samples <- parse_samples
sample_cols <- setNames(hcl.colors(length(parse_samples), palette = "Dark 3"), parse_samples)
score_ncores <- 2
max_mps_per_boxplot_page <- 10
max_mps_per_trend_page <- 12

script_run <- parse_start_run(
  "legacy_parse_highres_mp_t2t4_comparison_filter",
  parameters = list(mode = mode, legacy = TRUE, samples = paste(parse_samples, collapse = ",")),
  input_files = c(
    "parse_outs/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds",
    "parse_outs/by_samples/<sample>/Auto_<sample>_final.rds"
  ),
  output_files = "parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/*"
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

has_packages <- function(pkgs) {
  vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
}

load_or_stop <- function(pkgs) {
  available <- has_packages(pkgs)
  if (!all(available)) {
    stop("Missing required R packages: ", paste(names(available)[!available], collapse = ", "))
  }
  invisible(lapply(pkgs, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))))
}

get_counts <- function(obj) {
  suppressWarnings({
    tryCatch(
      SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts"),
      error = function(e) SeuratObject::GetAssayData(obj, assay = "RNA", slot = "counts")
    )
  })
}

chunk_vector <- function(x, n) {
  if (length(x) == 0) return(list())
  split(x, ceiling(seq_along(x) / n))
}

write_mp_gene_table <- function(mp_genes, path) {
  if (length(mp_genes) == 0) {
    gene_table <- data.frame(MP = character(), rank = integer(), gene = character(), stringsAsFactors = FALSE)
  } else {
    gene_table <- do.call(rbind, lapply(names(mp_genes), function(mp) {
      data.frame(MP = mp, rank = seq_along(mp_genes[[mp]]), gene = mp_genes[[mp]], stringsAsFactors = FALSE)
    }))
  }
  write.csv(gene_table, path, row.names = FALSE)
}

read_cell_cycle_genes <- function() {
  cc_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv"
  if (!file.exists(cc_path)) return(character())
  cc <- read.csv(cc_path, check.names = FALSE, stringsAsFactors = FALSE)
  unique(stats::na.omit(as.character(unlist(cc, use.names = FALSE))))
}

load_3ca_gene_sets <- function() {
  mp_csv <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
  if (!file.exists(mp_csv)) return(list())
  mp_list <- read.csv(mp_csv, check.names = FALSE, stringsAsFactors = FALSE)
  mp_list <- as.list(mp_list)
  mp_list <- lapply(mp_list, function(x) unique(x[x != "" & !is.na(x)]))
  names(mp_list) <- sub("^MP", "3CA_mp", names(mp_list))
  mp_list[lengths(mp_list) > 0]
}

make_3ca_label_table <- function(mp_genes, nMP) {
  sets_3ca <- load_3ca_gene_sets()
  if (length(sets_3ca) == 0 || length(mp_genes) == 0) {
    empty <- data.frame(MP = names(mp_genes), top_3ca_noncc = NA_character_, top_3ca_noncc_p_adj = NA_real_, stringsAsFactors = FALSE)
    write.csv(empty, file.path(out_dir, paste0("Auto_parse_highres_T2T4_top_3CA_noncellcycle_nMP", nMP, ".csv")), row.names = FALSE)
    return(empty)
  }

  cc_genes <- read_cell_cycle_genes()
  universe <- unique(c(unlist(mp_genes, use.names = FALSE), unlist(sets_3ca, use.names = FALSE)))
  cc_terms <- vapply(sets_3ca, function(genes) {
    cc_overlap <- length(intersect(genes, cc_genes))
    cc_overlap >= 5 && cc_overlap / length(genes) >= 0.05
  }, logical(1))

  enrich_rows <- do.call(rbind, lapply(names(mp_genes), function(mp) {
    genes <- unique(mp_genes[[mp]])
    do.call(rbind, lapply(names(sets_3ca), function(term) {
      term_genes <- unique(sets_3ca[[term]])
      overlap <- length(intersect(genes, term_genes))
      mat <- matrix(
        c(
          overlap,
          length(genes) - overlap,
          length(term_genes) - overlap,
          length(universe) - length(genes) - length(term_genes) + overlap
        ),
        nrow = 2
      )
      data.frame(
        MP = mp,
        term = term,
        overlap = overlap,
        p_value = tryCatch(stats::fisher.test(mat, alternative = "greater")$p.value, error = function(e) NA_real_),
        is_cell_cycle_term = isTRUE(cc_terms[[term]]),
        stringsAsFactors = FALSE
      )
    }))
  }))
  enrich_rows$p_adj <- stats::p.adjust(enrich_rows$p_value, method = "BH")
  write.csv(enrich_rows, file.path(out_dir, paste0("Auto_parse_highres_T2T4_3CA_base_enrichment_nMP", nMP, ".csv")), row.names = FALSE)

  label_table <- enrich_rows |>
    dplyr::filter(!is_cell_cycle_term, overlap > 0, is.finite(p_adj)) |>
    dplyr::arrange(MP, p_adj, dplyr::desc(overlap)) |>
    dplyr::group_by(MP) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::transmute(MP, top_3ca_noncc = term, top_3ca_noncc_p_adj = p_adj)

  label_table <- data.frame(MP = names(mp_genes), stringsAsFactors = FALSE) |>
    dplyr::left_join(label_table, by = "MP")
  write.csv(label_table, file.path(out_dir, paste0("Auto_parse_highres_T2T4_top_3CA_noncellcycle_nMP", nMP, ".csv")), row.names = FALSE)
  label_table
}

trend_similarity_order <- function(trend_summary) {
  retained <- trend_summary[trend_summary$retained, , drop = FALSE]
  if (nrow(retained) == 0) return(character())
  mp_vec <- retained$MP
  if (length(mp_vec) <= 2) return(mp_vec)
  mat <- as.matrix(retained[, paste0("mean_", parse_samples), drop = FALSE])
  rownames(mat) <- retained$MP
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[!is.finite(mat_scaled)] <- 0
  hc <- hclust(dist(mat_scaled), method = "ward.D2")
  rownames(mat_scaled)[hc$order]
}

make_display_labels <- function(trend_summary, nMP) {
  label_path <- file.path(out_dir, paste0("Auto_parse_highres_T2T4_top_3CA_noncellcycle_nMP", nMP, ".csv"))
  label_table <- if (file.exists(label_path)) {
    read.csv(label_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    data.frame(MP = trend_summary$MP, top_3ca_noncc = NA_character_, stringsAsFactors = FALSE)
  }
  label_df <- trend_summary |>
    dplyr::left_join(label_table, by = "MP") |>
    dplyr::mutate(
      top_3ca_noncc = ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc),
      single_programme_label = ifelse(numberPrograms == 1, "\n[1 programme]", ""),
      display_label = paste0(MP, "\n", top_3ca_noncc, single_programme_label)
    )
  setNames(label_df$display_label, label_df$MP)
}

make_activity_plot <- function(ucell_scores, cell_meta, mp_order, title_text, label_map) {
  activity_long <- as.data.frame(ucell_scores[, mp_order, drop = FALSE]) |>
    tibble::rownames_to_column("cell") |>
    dplyr::left_join(cell_meta, by = "cell") |>
    dplyr::filter(sample %in% parse_samples) |>
    tidyr::pivot_longer(cols = dplyr::all_of(mp_order), names_to = "MP", values_to = "score") |>
    dplyr::mutate(sample = factor(sample, levels = parse_samples), MP = factor(MP, levels = mp_order))

  activity_stats <- activity_long |>
    dplyr::group_by(MP) |>
    dplyr::summarise(p_value = tryCatch(stats::kruskal.test(score ~ sample)$p.value, error = function(e) NA_real_), .groups = "drop") |>
    dplyr::mutate(
      p_adj = stats::p.adjust(p_value, method = "BH"),
      significance = dplyr::case_when(is.na(p_adj) ~ "", p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", p_adj < 0.05 ~ "*", TRUE ~ "ns")
    )

  annot_df <- activity_long |>
    dplyr::group_by(MP) |>
    dplyr::summarise(y_pos = max(score, na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(activity_stats, by = "MP") |>
    dplyr::mutate(y_pos = y_pos + 0.025, label = ifelse(!is.na(p_adj) & p_adj < 0.05, significance, ""))

  p <- ggplot2::ggplot(activity_long, ggplot2::aes(x = MP, y = score, fill = sample, color = sample)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.82), width = 0.62, outlier.shape = NA, alpha = 0.78, linewidth = 0.28, color = "black") +
    ggplot2::geom_text(data = dplyr::filter(annot_df, label != ""), ggplot2::aes(x = MP, y = y_pos, label = label), inherit.aes = FALSE, size = 3.8, fontface = "bold") +
    ggplot2::scale_fill_manual(values = sample_cols, drop = FALSE, name = "Sample") +
    ggplot2::scale_color_manual(values = sample_cols, guide = "none", drop = FALSE) +
    ggplot2::scale_x_discrete(labels = label_map[mp_order]) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.14))) +
    ggplot2::labs(title = title_text, subtitle = "Per-cell UCell scores; stars mark BH-adjusted Kruskal-Wallis p < 0.05 across samples.", x = NULL, y = "UCell score") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 24),
      plot.subtitle = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 10, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black"),
      axis.line.x = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold", size = 20),
      legend.text = ggplot2::element_text(size = 20),
      plot.margin = ggplot2::margin(12, 18, 12, 12)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  list(plot = p, stats = activity_stats)
}

plot_retained_outputs <- function(ucell_scores, cell_meta, sample_summary, trend_summary, nMP) {
  retained_mps <- trend_similarity_order(trend_summary)
  if (length(retained_mps) == 0) {
    writeLines("No metaprograms passed the T2/T4 > T0/eR4 mean UCell filter.", file.path(out_dir, paste0("Auto_parse_highres_T2T4_no_retained_nMP", nMP, ".txt")))
    return(invisible(NULL))
  }

  label_map <- make_display_labels(trend_summary, nMP)
  boxplot_stats <- list()
  pdf(file.path(out_dir, paste0("Auto_parse_highres_T2T4_activity_boxplots_nMP", nMP, "_selected.pdf")), width = 20, height = 9, useDingbats = FALSE)
  chunks <- chunk_vector(retained_mps, max_mps_per_boxplot_page)
  for (i in seq_along(chunks)) {
    activity <- make_activity_plot(ucell_scores, cell_meta, chunks[[i]], paste0("T2/T4-high Parse MP activity - trend similarity page ", i), label_map)
    print(activity$plot)
    boxplot_stats[[i]] <- activity$stats
  }
  dev.off()
  write.csv(dplyr::bind_rows(boxplot_stats), file.path(out_dir, paste0("Auto_parse_highres_T2T4_activity_boxplots_nMP", nMP, "_selected_stats.csv")), row.names = FALSE)

  retained_summary <- sample_summary |>
    dplyr::filter(MP %in% retained_mps) |>
    dplyr::mutate(sample = factor(sample, levels = parse_samples), MP = factor(MP, levels = retained_mps))

  trend_long <- retained_summary |>
    dplyr::select(MP, sample, mean_score, median_score) |>
    tidyr::pivot_longer(cols = c(mean_score, median_score), names_to = "summary_stat", values_to = "score") |>
    dplyr::mutate(summary_stat = dplyr::recode(summary_stat, mean_score = "Mean", median_score = "Median"), display_label = label_map[as.character(MP)])
  trend_y_limits <- range(trend_long$score, na.rm = TRUE)
  trend_y_pad <- diff(trend_y_limits) * 0.05
  if (!is.finite(trend_y_pad) || trend_y_pad == 0) trend_y_pad <- 0.01
  trend_y_limits <- trend_y_limits + c(-trend_y_pad, trend_y_pad)

  pdf(file.path(out_dir, paste0("Auto_parse_highres_T2T4_mean_median_trends_nMP", nMP, "_selected.pdf")), width = 14, height = 10, useDingbats = FALSE)
  chunks <- chunk_vector(retained_mps, max_mps_per_trend_page)
  for (i in seq_along(chunks)) {
    chunk_labels <- label_map[chunks[[i]]]
    p <- trend_long |>
      dplyr::filter(MP %in% chunks[[i]]) |>
      dplyr::mutate(display_label = factor(display_label, levels = chunk_labels)) |>
      ggplot2::ggplot(ggplot2::aes(x = sample, y = score, group = summary_stat, linetype = summary_stat)) +
      ggplot2::geom_line(color = "grey25", linewidth = 0.5) +
      ggplot2::geom_point(ggplot2::aes(fill = sample), shape = 21, size = 2.7, color = "black") +
      ggplot2::scale_fill_manual(values = sample_cols, guide = "none") +
      ggplot2::scale_linetype_manual(values = c("Mean" = "solid", "Median" = "dashed"), name = NULL) +
      ggplot2::facet_wrap(~display_label, scales = "fixed", ncol = 4) +
      ggplot2::coord_cartesian(ylim = trend_y_limits) +
      ggplot2::labs(title = paste0("T2/T4-high mean and median UCell trend - page ", i), x = NULL, y = "UCell score") +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"), strip.text = ggplot2::element_text(face = "bold", size = 8.5), legend.position = "top")
    print(p)
  }
  dev.off()

  mat_df <- retained_summary |>
    dplyr::select(MP, sample, mean_score) |>
    tidyr::pivot_wider(names_from = sample, values_from = mean_score) |>
    as.data.frame()
  row.names(mat_df) <- mat_df$MP
  mat <- as.matrix(mat_df[retained_mps, parse_samples, drop = FALSE])
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[!is.finite(mat_scaled)] <- 0
  row_cluster <- if (nrow(mat_scaled) > 1) hclust(dist(mat_scaled), method = "ward.D2") else FALSE
  heat_cols <- grDevices::colorRampPalette(c("white", "#fee0d2", "#fc9272", "#de2d26", "#67000d"))(100)
  heat_breaks <- seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = length(heat_cols) + 1)
  if (length(unique(as.numeric(mat))) == 1) heat_breaks <- seq(0, max(mat, na.rm = TRUE) + 1e-6, length.out = length(heat_cols) + 1)
  row_anno <- trend_summary |>
    dplyr::filter(MP %in% rownames(mat)) |>
    dplyr::select(MP, mean_filter_pass, median_filter_pass, numberPrograms, silhouette) |>
    dplyr::mutate(
      mean_filter_pass = ifelse(mean_filter_pass, "mean_pass", "mean_fail"),
      median_filter_pass = ifelse(median_filter_pass, "median_pass", "median_fail"),
      numberPrograms = paste0("nProg_", numberPrograms),
      silhouette = sprintf("%.3f", silhouette)
    ) |>
    as.data.frame()
  row.names(row_anno) <- row_anno$MP
  row_anno <- row_anno[rownames(mat), c("mean_filter_pass", "median_filter_pass", "numberPrograms", "silhouette"), drop = FALSE]

  pdf(file.path(out_dir, paste0("Auto_parse_highres_T2T4_selected_mean_activity_heatmap_nMP", nMP, ".pdf")), width = 9.5, height = max(6, 0.32 * nrow(mat) + 2.5), useDingbats = FALSE)
  pheatmap::pheatmap(mat, color = heat_cols, breaks = heat_breaks, cluster_rows = row_cluster, cluster_cols = FALSE, annotation_row = row_anno, labels_row = label_map[rownames(mat)], border_color = "white", fontsize = 18, fontsize_row = 10, fontsize_col = 11, main = "T2/T4-high MP mean UCell activity")
  dev.off()
  png(file.path(out_dir, paste0("Auto_parse_highres_T2T4_selected_mean_activity_heatmap_nMP", nMP, ".png")), width = 2850, height = max(1800, 96 * nrow(mat) + 750), res = 300)
  pheatmap::pheatmap(mat, color = heat_cols, breaks = heat_breaks, cluster_rows = row_cluster, cluster_cols = FALSE, annotation_row = row_anno, labels_row = label_map[rownames(mat)], border_color = "white", fontsize = 16, fontsize_row = 10, fontsize_col = 11, main = "T2/T4-high MP mean UCell activity")
  dev.off()
}

get_or_make_metaprograms <- function(nMP, geneNMF.programs) {
  base_path <- file.path(base_highres_dir, paste0("Auto_parse_highres_geneNMF_metaprograms_nMP", nMP, ".rds"))
  out_path <- file.path(out_dir, paste0("Auto_parse_highres_T2T4_geneNMF_metaprograms_nMP", nMP, ".rds"))
  if (file.exists(base_path)) {
    geneNMF.metaprograms <- readRDS(base_path)
  } else if (file.exists(out_path)) {
    geneNMF.metaprograms <- readRDS(out_path)
  } else {
    geneNMF.metaprograms <- GeneNMF::getMetaPrograms(geneNMF.programs, metric = "cosine", specificity.weight = 5, weight.explained = 0.5, nMP = nMP, min.confidence = 0.5)
  }
  saveRDS(geneNMF.metaprograms, out_path, compress = FALSE)
  geneNMF.metaprograms
}

run_score <- function() {
  load_or_stop(c("GeneNMF", "UCell", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble", "ggplot2", "pheatmap"))

  geneNMF_out_path <- file.path(parse_mp_dir, "Auto_parse_geneNMF_outs.rds")
  if (!file.exists(geneNMF_out_path)) stop("Missing GeneNMF programmes: ", geneNMF_out_path)
  message("Loading GeneNMF programmes: ", geneNMF_out_path)
  geneNMF.programs <- readRDS(geneNMF_out_path)
  total_nmf_programs <- sum(vapply(geneNMF.programs, function(x) ncol(x$w), numeric(1)))
  nMP <- as.integer(total_nmf_programs / 2)
  if (total_nmf_programs / 2 != nMP) stop("Total NMF programme count is odd: ", total_nmf_programs)
  write.csv(data.frame(total_nmf_programs = total_nmf_programs, nMP = nMP, programmes_per_MP_target = total_nmf_programs / nMP, filter = "mean_T2_and_T4_gt_mean_T0_and_eR4", samples = paste(parse_samples, collapse = ",")), file.path(out_dir, paste0("Auto_parse_highres_T2T4_nMP", nMP, "_config.csv")), row.names = FALSE)

  geneNMF.metaprograms <- get_or_make_metaprograms(nMP, geneNMF.programs)
  mp_genes <- geneNMF.metaprograms$metaprograms.genes
  write_mp_gene_table(mp_genes, file.path(out_dir, paste0("Auto_parse_highres_T2T4_mp_genes_nMP", nMP, ".csv")))
  make_3ca_label_table(mp_genes, nMP)

  program_membership <- data.frame(
    nmf_programme = names(geneNMF.metaprograms$programs.clusters),
    MP = paste0("MP", as.integer(geneNMF.metaprograms$programs.clusters)),
    sample = sub("\\..*$", "", names(geneNMF.metaprograms$programs.clusters)),
    k = as.integer(sub("^[^.]+\\.k([0-9]+)\\..*$", "\\1", names(geneNMF.metaprograms$programs.clusters))),
    stringsAsFactors = FALSE
  )
  write.csv(program_membership, file.path(out_dir, paste0("Auto_parse_highres_T2T4_program_membership_nMP", nMP, ".csv")), row.names = FALSE)

  ucell_path <- file.path(out_dir, paste0("Auto_parse_highres_T2T4_UCell_scores_nMP", nMP, ".rds"))
  cell_meta_path <- file.path(out_dir, paste0("Auto_parse_highres_T2T4_cell_metadata_nMP", nMP, ".rds"))
  base_ucell_path <- file.path(base_highres_dir, paste0("Auto_parse_highres_UCell_scores_nMP", nMP, ".rds"))
  base_cell_meta_path <- file.path(base_highres_dir, paste0("Auto_parse_highres_cell_metadata_nMP", nMP, ".rds"))

  if (file.exists(base_ucell_path) && file.exists(base_cell_meta_path)) {
    message("Loading cached high-resolution UCell scores from: ", base_ucell_path)
    ucell_scores <- readRDS(base_ucell_path)
    cell_meta <- readRDS(base_cell_meta_path)
    saveRDS(ucell_scores, ucell_path, compress = FALSE)
    saveRDS(cell_meta, cell_meta_path, compress = FALSE)
  } else if (file.exists(ucell_path) && file.exists(cell_meta_path)) {
    message("Loading existing T2/T4 UCell scores: ", ucell_path)
    ucell_scores <- readRDS(ucell_path)
    cell_meta <- readRDS(cell_meta_path)
  } else {
    sample_files <- setNames(file.path(qc_dir, "by_samples", parse_samples, paste0("Auto_", parse_samples, "_final.rds")), parse_samples)
    missing_files <- names(sample_files)[!file.exists(sample_files)]
    if (length(missing_files) > 0) stop("Missing final sample RDS files: ", paste(missing_files, collapse = ", "))
    sample_genes <- lapply(sample_files, function(path) rownames(readRDS(path)))
    common_genes <- Reduce(intersect, sample_genes)
    counts_list <- list()
    meta_list <- list()
    for (sample in parse_samples) {
      message("Loading counts for ", sample)
      obj <- readRDS(sample_files[[sample]])
      old_cells <- colnames(obj)
      new_cells <- paste(sample, old_cells, sep = "_")
      counts <- get_counts(obj)[common_genes, , drop = FALSE]
      colnames(counts) <- new_cells
      counts_list[[sample]] <- counts
      meta_list[[sample]] <- data.frame(cell = new_cells, original_cell = old_cells, sample = sample, stringsAsFactors = FALSE)
      rm(obj, counts)
      gc()
    }
    cell_meta <- dplyr::bind_rows(meta_list)
    rownames(cell_meta) <- cell_meta$cell
    counts_all <- do.call(cbind, counts_list)
    score_features <- lapply(mp_genes, intersect, rownames(counts_all))
    score_features <- score_features[lengths(score_features) > 0]
    ucell_scores <- UCell::ScoreSignatures_UCell(matrix = counts_all, features = score_features, maxRank = 1500, chunk.size = 1000, ncores = score_ncores, force.gc = TRUE)
    ucell_scores <- as.data.frame(ucell_scores)
    colnames(ucell_scores) <- sub("_UCell$", "", colnames(ucell_scores))
    ucell_scores <- as.matrix(ucell_scores)
    saveRDS(ucell_scores, ucell_path, compress = FALSE)
    saveRDS(cell_meta, cell_meta_path, compress = FALSE)
  }

  mp_order <- intersect(names(mp_genes), colnames(ucell_scores))
  score_long <- as.data.frame(ucell_scores[, mp_order, drop = FALSE]) |>
    tibble::rownames_to_column("cell") |>
    dplyr::left_join(cell_meta[, c("cell", "sample")], by = "cell") |>
    dplyr::filter(sample %in% parse_samples) |>
    tidyr::pivot_longer(cols = dplyr::all_of(mp_order), names_to = "MP", values_to = "score")

  sample_summary <- score_long |>
    dplyr::group_by(MP, sample) |>
    dplyr::summarise(n_cells = dplyr::n(), mean_score = mean(score, na.rm = TRUE), median_score = stats::median(score, na.rm = TRUE), q1 = stats::quantile(score, 0.25, na.rm = TRUE), q3 = stats::quantile(score, 0.75, na.rm = TRUE), .groups = "drop")
  write.csv(sample_summary, file.path(out_dir, paste0("Auto_parse_highres_T2T4_sample_ucell_summary_nMP", nMP, ".csv")), row.names = FALSE)

  mean_wide <- sample_summary |>
    dplyr::select(MP, sample, mean_score) |>
    tidyr::pivot_wider(names_from = sample, values_from = mean_score, names_prefix = "mean_")
  median_wide <- sample_summary |>
    dplyr::select(MP, sample, median_score) |>
    tidyr::pivot_wider(names_from = sample, values_from = median_score, names_prefix = "median_")
  trend_summary <- dplyr::left_join(mean_wide, median_wide, by = "MP")
  trend_summary$mean_filter_pass <- with(trend_summary, mean_T2 > mean_T0 & mean_T2 > mean_eR4 & mean_T4 > mean_T0 & mean_T4 > mean_eR4)
  trend_summary$median_filter_pass <- with(trend_summary, median_T2 > median_T0 & median_T2 > median_eR4 & median_T4 > median_T0 & median_T4 > median_eR4)
  trend_summary$retained <- trend_summary$mean_filter_pass
  trend_summary$T2_minus_T0 <- trend_summary$mean_T2 - trend_summary$mean_T0
  trend_summary$T2_minus_eR4 <- trend_summary$mean_T2 - trend_summary$mean_eR4
  trend_summary$T4_minus_T0 <- trend_summary$mean_T4 - trend_summary$mean_T0
  trend_summary$T4_minus_eR4 <- trend_summary$mean_T4 - trend_summary$mean_eR4
  trend_summary$T4_minus_T2 <- trend_summary$mean_T4 - trend_summary$mean_T2

  metrics <- geneNMF.metaprograms$metaprograms.metrics
  metrics$MP <- rownames(metrics)
  metrics$numberPrograms <- suppressWarnings(as.integer(metrics$numberPrograms))
  trend_summary <- trend_summary |>
    dplyr::left_join(metrics, by = "MP")
  trend_order <- trend_similarity_order(trend_summary)
  trend_summary$trend_similarity_order <- match(trend_summary$MP, trend_order)
  trend_summary <- trend_summary |>
    dplyr::arrange(!retained, trend_similarity_order, dplyr::desc(numberPrograms), dplyr::desc(silhouette), MP)
  write.csv(trend_summary, file.path(out_dir, paste0("Auto_parse_highres_T2T4_filter_summary_nMP", nMP, ".csv")), row.names = FALSE)

  retained_mps <- trend_similarity_order(trend_summary)
  retained_genes <- mp_genes[retained_mps]
  saveRDS(retained_genes, file.path(out_dir, paste0("Auto_parse_highres_T2T4_selected_mp_genes_nMP", nMP, ".rds")), compress = FALSE)
  write_mp_gene_table(retained_genes, file.path(out_dir, paste0("Auto_parse_highres_T2T4_selected_mp_genes_nMP", nMP, ".csv")))
  write.csv(dplyr::filter(trend_summary, retained), file.path(out_dir, paste0("Auto_parse_highres_T2T4_filter_retained_nMP", nMP, ".csv")), row.names = FALSE)

  message("Retained ", length(retained_mps), " MPs with mean T2/T4 > mean T0/eR4.")
  plot_retained_outputs(ucell_scores, cell_meta, sample_summary, trend_summary, nMP)
  invisible(list(nMP = nMP, retained_mps = retained_mps))
}

load_enrichment_references <- function() {
  hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
  hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]
  sets_3ca <- load_3ca_gene_sets()
  mp_term2gene <- data.frame(term = rep(names(sets_3ca), lengths(sets_3ca)), gene = unlist(sets_3ca), row.names = NULL)
  mp_term2name <- data.frame(term = unique(mp_term2gene$term), name = unique(mp_term2gene$term))
  individual_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/"
  custom_files <- if (dir.exists(individual_dir)) list.files(individual_dir, pattern = "\\.rds$", full.names = TRUE) else character()
  custom_refs <- lapply(custom_files, readRDS)
  names(custom_refs) <- sub("\\.rds$", "", sub(".*enrich_dev_", "", basename(custom_files)))
  list(hallmark_term2gene = hallmark_term2gene, hallmark_term2name = hallmark_term2name, mp_term2gene = mp_term2gene, mp_term2name = mp_term2name, custom_refs = custom_refs)
}

run_enrichment_heatmap <- function(cluster_enrich, selected_mp_genes, refs, element, output_prefix, trend_summary, label_map, top_per_program = 8, top_n = 80, cap = 7, cols = grDevices::colorRampPalette(c("white", "#fee0d2", "#fc9272", "#de2d26", "#67000d"))(100), fontsize = 18, fontsize_row = 10, fontsize_col = 11) {
  is_custom <- !element %in% c("GO", "Hallmark", "MPs_3CA")
  df_list <- lapply(names(cluster_enrich), function(prog) {
    er <- cluster_enrich[[prog]][[element]]
    if (is.null(er)) return(NULL)
    r <- tryCatch(er@result, error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0) return(NULL)
    r_sig <- r[which(r$p.adjust < 0.05 & r$p.adjust > 0), ]
    data_source <- if (is_custom) r else r_sig
    if (nrow(data_source) == 0 && !is_custom) return(NULL)
    term <- if ("Description" %in% colnames(data_source)) data_source$Description else data_source$ID
    data.frame(Program = prog, Term = term, padj = data_source$p.adjust, Overlap = data_source$GeneRatio, stringsAsFactors = FALSE)
  })
  df <- dplyr::bind_rows(df_list)
  if (nrow(df) == 0) {
    message("  -> Skipping enrichment heatmap for ", element, ": No significant terms.")
    return(FALSE)
  }

  if (is_custom) {
    if (!element %in% names(refs$custom_refs)) return(FALSE)
    terms_use <- as.character(refs$custom_refs[[element]]$TERM2NAME$term)
  } else {
    terms_use <- df |>
      dplyr::filter(padj < 0.05) |>
      dplyr::arrange(Program, padj) |>
      dplyr::group_by(Program) |>
      dplyr::slice_head(n = top_per_program) |>
      dplyr::ungroup() |>
      dplyr::distinct(Term) |>
      dplyr::pull(Term)
    if (length(terms_use) == 0) {
      message("  -> Skipping enrichment heatmap for ", element, ": No terms after filtering.")
      return(FALSE)
    }
    if (length(terms_use) > top_n) {
      terms_use <- df |>
        dplyr::filter(Term %in% terms_use) |>
        dplyr::group_by(Term) |>
        dplyr::summarise(min_p = min(padj), .groups = "drop") |>
        dplyr::arrange(min_p) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::pull(Term)
    }
  }

  ordered_mps <- trend_similarity_order(trend_summary)
  ordered_mps <- ordered_mps[ordered_mps %in% names(selected_mp_genes)]
  full_grid <- expand.grid(Term = terms_use, Program = ordered_mps, stringsAsFactors = FALSE)
  final_df <- full_grid |>
    dplyr::left_join(df, by = c("Term", "Program")) |>
    dplyr::mutate(score = tidyr::replace_na(pmin(-log10(padj), cap), 0), display_text = tidyr::replace_na(Overlap, ""))

  mat <- final_df |>
    dplyr::select(Term, Program, score) |>
    tidyr::pivot_wider(names_from = Program, values_from = score) |>
    as.data.frame()
  row.names(mat) <- mat$Term
  mat <- as.matrix(dplyr::select(mat, -Term))
  text_mat <- final_df |>
    dplyr::select(Term, Program, display_text) |>
    tidyr::pivot_wider(names_from = Program, values_from = display_text) |>
    as.data.frame()
  row.names(text_mat) <- text_mat$Term
  text_mat <- as.matrix(dplyr::select(text_mat, -Term))
  mat <- mat[terms_use, ordered_mps[ordered_mps %in% colnames(mat)], drop = FALSE]
  text_mat <- text_mat[terms_use, colnames(mat), drop = FALSE]
  if (nrow(mat) == 0 || ncol(mat) == 0) return(FALSE)
  mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  if (!is_custom) {
    best_mp <- colnames(mat)[max.col(mat, ties.method = "first")]
    row_order <- order(match(best_mp, colnames(mat)), -rowSums(mat))
    mat <- mat[row_order, , drop = FALSE]
    text_mat <- text_mat[row_order, , drop = FALSE]
  }
  row_gaps <- if (!is_custom && nrow(mat) > 1) {
    groups <- colnames(mat)[max.col(mat, ties.method = "first")]
    which(groups[-length(groups)] != groups[-1])
  } else {
    NULL
  }
  # col_anno <- trend_summary |>
  #   dplyr::filter(MP %in% colnames(mat)) |>
  #   dplyr::select(MP, mean_filter_pass, median_filter_pass, numberPrograms) |>
  #   dplyr::mutate(
  #     mean_filter_pass = ifelse(mean_filter_pass, "mean_pass", "mean_fail"),
  #     median_filter_pass = ifelse(median_filter_pass, "median_pass", "median_fail"),
  #     numberPrograms = paste0("nProg_", numberPrograms)
  #   ) |>
  #   as.data.frame()
  # row.names(col_anno) <- col_anno$MP
  # col_anno <- col_anno[colnames(mat), c("mean_filter_pass", "median_filter_pass", "numberPrograms"), drop = FALSE]
  breaks <- seq(0, cap, length.out = length(cols) + 1)
  # EXPLICIT GRID DRAW FOR MULTI-PAGE PDF RELIABILITY
  grid::grid.newpage()
  p_obj <- pheatmap::pheatmap(
    mat, 
    display_numbers = text_mat, 
    number_color = "black", 
    fontsize_number = fontsize_row * 0.95, 
    labels_col = label_map[colnames(mat)], 
    color = cols, 
    breaks = breaks, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE, 
    gaps_row = row_gaps, 
    # annotation_col = col_anno, 
    border_color = NA, 
    show_colnames = TRUE, 
    angle_col = 45, 
    fontsize = fontsize,
    fontsize_row = fontsize_row, 
    fontsize_col = fontsize_col, 
    main = paste0(element, " Enrichment (-log10 padj)"),
    silent = TRUE
  )
  grid::grid.draw(p_obj$gtable)
  message("  -> Successfully added ", element, " enrichment page to PDF (", length(terms_use), " terms).")
  TRUE
}

run_enrich <- function() {
  load_or_stop(c("clusterProfiler", "org.Hs.eg.db", "msigdbr", "dplyr", "tidyr", "pheatmap", "ggplot2"))
  config_files <- list.files(out_dir, pattern = "^Auto_parse_highres_T2T4_nMP[0-9]+_config\\.csv$", full.names = TRUE)
  if (length(config_files) == 0) stop("No T2/T4 config found in ", out_dir, ". Run --mode=score first.")
  config <- read.csv(config_files[order(file.info(config_files)$mtime, decreasing = TRUE)[1]])
  nMP <- as.integer(config$nMP[1])
  selected_genes_path <- file.path(out_dir, paste0("Auto_parse_highres_T2T4_selected_mp_genes_nMP", nMP, ".rds"))
  trend_path <- file.path(out_dir, paste0("Auto_parse_highres_T2T4_filter_summary_nMP", nMP, ".csv"))
  if (!file.exists(selected_genes_path) || !file.exists(trend_path)) stop("Missing T2/T4 selected MP outputs. Run --mode=score first.")
  selected_mp_genes <- readRDS(selected_genes_path)
  trend_summary <- read.csv(trend_path, check.names = FALSE, stringsAsFactors = FALSE)
  selected_mp_genes <- selected_mp_genes[trend_similarity_order(trend_summary)]
  if (length(selected_mp_genes) == 0) {
    writeLines("No retained MPs; enrichment annotation was not run.", file.path(out_dir, paste0("Auto_parse_highres_T2T4_enrichment_skipped_nMP", nMP, ".txt")))
    message("No retained MPs; skipping enrichment.")
    return(invisible(NULL))
  }
  label_map <- make_display_labels(trend_summary, nMP)
  refs <- load_enrichment_references()
  message("Running enrichment for ", length(selected_mp_genes), " T2/T4-high retained MPs.")
  cluster_enrich <- lapply(names(selected_mp_genes), function(mp_name) {
    genes <- selected_mp_genes[[mp_name]]
    message("Processing enrichment for ", mp_name)
    res_GO <- clusterProfiler::enrichGO(gene = genes, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", qvalueCutoff = 0.05, readable = TRUE)
    res_H <- clusterProfiler::enricher(gene = genes, TERM2GENE = refs$hallmark_term2gene, TERM2NAME = refs$hallmark_term2name, qvalueCutoff = 0.05)
    res_M <- clusterProfiler::enricher(gene = genes, TERM2GENE = refs$mp_term2gene, TERM2NAME = refs$mp_term2name, qvalueCutoff = 0.05)
    res_custom_list <- lapply(names(refs$custom_refs), function(ref_name) {
      message("  -> Running custom enrichment: ", ref_name)
      clusterProfiler::enricher(gene = genes, TERM2GENE = refs$custom_refs[[ref_name]]$TERM2GENE, TERM2NAME = refs$custom_refs[[ref_name]]$TERM2NAME, pAdjustMethod = "BH", qvalueCutoff = 0.05)
    })
    names(res_custom_list) <- names(refs$custom_refs)
    c(list(rep_prog = mp_name, genes = genes, GO = res_GO, Hallmark = res_H, MPs_3CA = res_M), res_custom_list)
  })
  names(cluster_enrich) <- names(selected_mp_genes)
  saveRDS(cluster_enrich, file.path(out_dir, paste0("Auto_parse_highres_T2T4_cluster_enrich_nMP", nMP, ".rds")), compress = FALSE)
  output_prefix <- paste0("Auto_parse_highres_T2T4_enrich_nMP", nMP)
  # Extreme font sizes for meeting presentation clarity
  pdf_path <- file.path(out_dir, paste0("Auto_parse_highres_T2T4_enrichment_annotation_nMP", nMP, ".pdf"))
  message("Generating multi-page PDF: ", pdf_path)
  pdf(pdf_path, width = 20, height = 14, useDingbats = FALSE)
  # Ordered by user request: 3CA, Hallmark, GO, Early Embryogenesis, Organogenesis, etc.
  plotted <- c()
  plotted <- c(plotted, run_enrichment_heatmap(cluster_enrich, selected_mp_genes, refs, "MPs_3CA", output_prefix, trend_summary, label_map, top_per_program = 8, top_n = 80))
  plotted <- c(plotted, run_enrichment_heatmap(cluster_enrich, selected_mp_genes, refs, "Hallmark", output_prefix, trend_summary, label_map, top_per_program = 8, top_n = 80))
  plotted <- c(plotted, run_enrichment_heatmap(cluster_enrich, selected_mp_genes, refs, "GO", output_prefix, trend_summary, label_map, top_per_program = 6, top_n = 60))
  
  custom_order <- c(
    "Early_Embryogenesis", 
    "Organogenesis_major", 
    "Organogenesis_sub", 
    "Adult_Epithelium", 
    "Barretts_Oesophagus",
    "Normal_Development_long", 
    "Normal_Development_short"
  )
  
  for (element in custom_order) {
    if (element %in% names(refs$custom_refs)) {
      plotted <- c(plotted, run_enrichment_heatmap(cluster_enrich, selected_mp_genes, refs, element, output_prefix, trend_summary, label_map, top_per_program = 8, top_n = 80))
    }
  }
  dev.off()
  if (!any(plotted)) writeLines("Enrichment ran, but no heatmaps had plottable content.", file.path(out_dir, paste0("Auto_parse_highres_T2T4_enrichment_no_plottable_heatmaps_nMP", nMP, ".txt")))
  message("T2/T4 enrichment outputs written to: ", out_dir)
}

run_excel <- function() {
  load_or_stop(c("openxlsx", "dplyr"))
  
  config_files <- list.files(out_dir, pattern = "^Auto_parse_highres_T2T4_nMP[0-9]+_config\\.csv$", full.names = TRUE)
  if (length(config_files) == 0) stop("No T2/T4 config found in ", out_dir, ". Run --mode=score first.")
  config <- read.csv(config_files[order(file.info(config_files)$mtime, decreasing = TRUE)[1]])
  nMP <- as.integer(config$nMP[1])
  
  trend_path <- file.path(out_dir, paste0("Auto_parse_highres_T2T4_filter_summary_nMP", nMP, ".csv"))
  genes_path <- file.path(out_dir, paste0("Auto_parse_highres_T2T4_selected_mp_genes_nMP", nMP, ".rds"))
  
  if (!file.exists(trend_path) || !file.exists(genes_path)) stop("Missing T2/T4 selected MP outputs. Run --mode=score first.")
  
  trend_summary <- read.csv(trend_path, check.names = FALSE, stringsAsFactors = FALSE)
  mp_genes <- readRDS(genes_path)
  
  retained_df <- trend_summary[trend_summary$retained == TRUE | trend_summary$retained == "TRUE", ]
  
  if (nrow(retained_df) == 0) {
    message("No retained categories found for excel generation.")
    return(invisible(NULL))
  }
  
  build_mp_matrix <- function(mp_names_vec) {
    if (length(mp_names_vec) == 0) return(NULL)
    max_g <- max(sapply(mp_names_vec, function(x) length(mp_genes[[x]])))
    n_mp <- length(mp_names_vec)
    
    n_rows <- max_g + 2
    mat <- matrix(NA_character_, nrow = n_rows, ncol = n_mp)
    for (i in seq_along(mp_names_vec)) {
      mp <- mp_names_vec[i]
      n_prog <- retained_df$numberPrograms[retained_df$MP == mp]
      if (length(n_prog) > 0 && !is.na(n_prog[1])) {
        prog_label <- ifelse(n_prog[1] == 1, "1 programme", paste0(n_prog[1], " programmes"))
        mat[1, i] <- paste0(mp, " (", prog_label, ")")
      } else {
        mat[1, i] <- mp
      }
      mat[2, i] <- "" 
      genes <- mp_genes[[mp]]
      if (length(genes) > 0) {
        mat[3:(length(genes)+2), i] <- genes
      }
    }
    return(as.data.frame(mat, stringsAsFactors = FALSE))
  }
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "MP_Genes")
  
  current_col <- 1
  section_style <- openxlsx::createStyle(fontSize = 14, textDecoration = "bold", fgFill = "#FFC000")
  mp_name_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
  desc_style <- openxlsx::createStyle(fgFill = "#F2F2F2")
  
  mps_in_group <- retained_df$MP
  mps_in_group <- mps_in_group[mps_in_group %in% names(mp_genes)]
  mps_in_group <- mps_in_group[order(retained_df$trend_similarity_order[match(mps_in_group, retained_df$MP)])]
  
  if (length(mps_in_group) > 0) {
    df_group <- build_mp_matrix(mps_in_group)
    
    openxlsx::writeData(wb, sheet = 1, x = "RETAINED MPs", startCol = current_col, startRow = 1)
    openxlsx::writeData(wb, sheet = 1, x = df_group, startCol = current_col, startRow = 2, colNames = FALSE)
    
    openxlsx::addStyle(wb, sheet = 1, section_style, rows = 1, cols = current_col, gridExpand = TRUE)
    openxlsx::addStyle(wb, sheet = 1, mp_name_style, rows = 2, cols = current_col:(current_col + ncol(df_group) - 1), gridExpand = TRUE)
    openxlsx::addStyle(wb, sheet = 1, desc_style, rows = 3, cols = current_col:(current_col + ncol(df_group) - 1), gridExpand = TRUE)
    
    for (i in current_col:(current_col + ncol(df_group) - 1)) {
      openxlsx::setColWidths(wb, 1, cols = i, widths = 25)
    }
  }
  
  output_path <- file.path(out_dir, paste0("Auto_parse_highres_T2T4_mp_genes_summary_nMP", nMP, ".xlsx"))
  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
  message("Excel summary written to: ", output_path)
}

if (mode %in% c("all", "score")) {
  run_score()
}

if (mode %in% c("all", "enrich")) {
  run_enrich()
}

if (mode %in% c("all", "excel")) {
  excel_pkgs <- c("openxlsx", "dplyr")
  if (mode == "all" && !all(has_packages(excel_pkgs))) {
    message("Skipping excel generation in --mode=all because openxlsx/dplyr is missing.")
  } else {
    run_excel()
  }
}

script_run_status <- "success"
message("Done. Outputs written to: ", out_dir)
