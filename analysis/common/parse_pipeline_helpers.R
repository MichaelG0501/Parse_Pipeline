####################
# parse_pipeline_helpers.R
#
# Reusable R helpers for Parse_Pipeline workflows.
#
# Inputs:
#   Source after analysis/common/parse_pipeline_config.R, or source directly.
#
# Outputs:
#   No files are written except by parse_write_mp_gene_table().
#
# Methodology:
#   analysis/methodology/common/shared_configuration_and_logging_methodology.md
####################

if (!exists("parse_project_root")) {
  source("analysis/common/parse_pipeline_config.R")
}

parse_has_packages <- function(pkgs) {
  stats::setNames(
    vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
    pkgs
  )
}

parse_load_or_stop <- function(pkgs) {
  available <- parse_has_packages(pkgs)
  if (!all(available)) {
    stop("Missing required R packages: ", paste(names(available)[!available], collapse = ", "))
  }
  invisible(lapply(pkgs, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))))
}

parse_get_counts <- function(obj, assay = "RNA") {
  suppressWarnings({
    tryCatch(
      SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
      error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
    )
  })
}

parse_chunk_vector <- function(x, n) {
  if (length(x) == 0) return(list())
  split(x, ceiling(seq_along(x) / n))
}

parse_write_mp_gene_table <- function(mp_genes, path) {
  if (length(mp_genes) == 0) {
    gene_table <- data.frame(MP = character(), rank = integer(), gene = character(), stringsAsFactors = FALSE)
  } else {
    gene_table <- do.call(rbind, lapply(names(mp_genes), function(mp) {
      data.frame(MP = mp, rank = seq_along(mp_genes[[mp]]), gene = mp_genes[[mp]], stringsAsFactors = FALSE)
    }))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(gene_table, path, row.names = FALSE)
}

parse_z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]
  study_sd <- clust_df |>
    dplyr::group_by(.study) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(colnames(mat)), ~ stats::sd(.x, na.rm = TRUE)), .groups = "drop") |>
    tibble::column_to_rownames(".study") |>
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1
  clust_centered <- clust_df |>
    dplyr::group_by(.sample) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) |>
    dplyr::ungroup()
  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

parse_slide_theme <- function(base_size = parse_plot_defaults$slide_base_size) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = "black", size = parse_plot_defaults$slide_axis_text_size),
      legend.text = ggplot2::element_text(size = parse_plot_defaults$slide_legend_text_size),
      legend.title = ggplot2::element_text(size = parse_plot_defaults$slide_legend_title_size, face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}
