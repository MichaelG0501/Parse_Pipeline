####################
# parse_tumour_stiffness_gene_module_heatmap.R
#
# Description:
#   Sample-resolved tumour stiffness gene-expression and module-activity
#   heatmap. The figure mirrors the grouped ComplexHeatmap layout used by the
#   PDO matched-FLOT response heatmap, with sample columns replacing state
#   columns and row blocks showing the requested mechanosensing, adhesion,
#   Rho/cytoskeletal, polarity, and contractility genes.
#
# Inputs:
#   parse_outs/Auto_parse_merged.rds
#
# Outputs:
#   parse_outs/cell_states/tumour_stiffness_gene_module/figures/Auto_parse_tumour_stiffness_module_gene_heatmap.pdf
#   parse_outs/cell_states/tumour_stiffness_gene_module/figures/Auto_parse_tumour_stiffness_module_gene_heatmap_no_pdo_sur1090.pdf
#   parse_outs/cell_states/tumour_stiffness_gene_module/figures/Auto_parse_tumour_stiffness_gene_boxplot_T0_T1_T4_eR4.pdf
#
# Cache / replot:
#   The script recomputes the small requested gene/module matrices from
#   parse_outs/Auto_parse_merged.rds on each run, then writes them to
#   intermediate/ for inspection. No expensive UCell or model fitting is used.
#
# Methodology:
#   analysis/methodology/cell_states/tumour_stiffness_gene_module_heatmap_methodology.md
#
# Downstream:
#   Terminal presentation/manuscript-facing heatmap. No active downstream
#   analysis currently depends on this output.
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_helpers.R")
source("analysis/common/parse_pipeline_logging.R")

script_name <- "parse_tumour_stiffness_gene_module_heatmap"
root_dir <- parse_project_root()
paths <- parse_paths(root_dir)
out_root <- file.path(paths$parse_outs, "cell_states", "tumour_stiffness_gene_module")
out_tiers <- parse_output_tiers(out_root, create = TRUE)

script_run <- parse_start_run(
  script_name,
  parameters = list(
    assay = "RNA",
    sample_column = parse_metadata_columns$sample,
    score_type = "mean log-normalized expression; heatmap values are row z-scores across samples",
    figure_contract = "Sample-resolved stiffness module activity and component gene expression"
  ),
  input_files = c("parse_outs/Auto_parse_merged.rds"),
  output_files = c(
    file.path(out_tiers$figures, "Auto_parse_tumour_stiffness_module_gene_heatmap.pdf"),
    file.path(out_tiers$figures, "Auto_parse_tumour_stiffness_module_gene_heatmap_no_pdo_sur1090.pdf"),
    file.path(out_tiers$figures, "Auto_parse_tumour_stiffness_gene_boxplot_T0_T1_T4_eR4.pdf")
  )
)
script_run_status <- "failed"
script_run_finished <- FALSE
on.exit({
  if (!isTRUE(script_run_finished)) {
    parse_finish_run(script_run, status = script_run_status)
  }
}, add = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

setwd(paths$parse_outs)

####################
# inputs and gene modules
####################

matrices_path <- file.path(out_tiers$intermediate, "Auto_parse_tumour_stiffness_matrices.rds")

if (file.exists(matrices_path)) {
  message("Loading cached intermediate matrices: ", matrices_path)
  cached_data <- readRDS(matrices_path)
  raw_score_mat <- cached_data$raw_score_mat
  z_score_mat <- cached_data$z_score_mat
  detection_mat <- cached_data$detection_mat
  feature_info <- cached_data$feature_info
  gene_availability <- cached_data$gene_availability
  sample_order <- cached_data$sample_order
  cell_counts <- cached_data$cell_counts
  
  missing_genes <- gene_availability$requested_gene[!gene_availability$detected]
} else {
  merged_path <- file.path(paths$parse_outs, "Auto_parse_merged.rds")
  if (!file.exists(merged_path)) {
    stop("Missing input Seurat object: ", merged_path)
  }

  stiffness_gene_groups <- list(
    "Overall tumour stiffness module" = character(),
    "Core mechanosensing / membrane tension" = c("PIEZO1", "MICAL2", "MARCKS"),
    "Integrin-focal adhesion-actin linkage" = c("ITGB1", "TLN1", "EZR", "RDX", "MSN"),
    "Rho GTPase cytoskeletal regulators" = c("RHOA", "RAC1", "CDC42"),
    "Cell polarity / cortical organisation" = c("LLGL1"),
    "Actin cytoskeleton and contractility" = c(
      "ACTB", "ACTG1", "ACTN1", "ACTN4", "VCL", "PXN", "FLNA",
      "MYH9", "MYH10", "MYL9", "MYL12A", "MYL12B", "ROCK1", "ROCK2",
      "DIAPH1", "CFL1", "CAPZA1", "CAPZB", "ARPC2", "ARPC3", "WASL",
      "PFN1", "VIM"
    )
  )
  all_requested_genes <- unique(unlist(stiffness_gene_groups[names(stiffness_gene_groups) != "Overall tumour stiffness module"]))
  stiffness_gene_groups[["Overall tumour stiffness module"]] <- all_requested_genes

  message("Loading merged Seurat object: ", merged_path)
  merged_obj <- readRDS(merged_path)
  if (!inherits(merged_obj, "Seurat")) {
    stop("Input is not a Seurat object: ", merged_path)
  }
  if (!"RNA" %in% Assays(merged_obj)) {
    stop("RNA assay not found in merged Seurat object.")
  }
  DefaultAssay(merged_obj) <- "RNA"

  sample_col <- parse_metadata_columns$sample
  if (!sample_col %in% colnames(merged_obj@meta.data)) {
    sample_col <- parse_metadata_columns$sample_fallbacks[
      parse_metadata_columns$sample_fallbacks %in% colnames(merged_obj@meta.data)
    ][1]
  }
  if (is.na(sample_col) || !nzchar(sample_col)) {
    stop("No sample metadata column found. Checked: ", paste(parse_metadata_columns$sample_fallbacks, collapse = ", "))
  }

  sample_vec <- as.character(merged_obj@meta.data[[sample_col]])
  sample_order <- intersect(parse_all_samples, unique(sample_vec))
  sample_order <- c(sample_order, setdiff(sort(unique(sample_vec)), sample_order))
  sample_factor <- factor(sample_vec, levels = sample_order)
  cell_counts <- table(sample_factor)

  get_log_data <- function(obj, features) {
    assay_obj <- obj[["RNA"]]
    assay_layers <- SeuratObject::Layers(assay_obj)
    data_layers <- assay_layers[grepl("^data(\\.|$)", assay_layers)]
    if (length(data_layers) > 1) {
      message("Extracting requested genes from ", length(data_layers), " Seurat v5 RNA data layers.")
      layer_mats <- lapply(data_layers, function(layer_name) {
        layer_mat <- SeuratObject::LayerData(obj, assay = "RNA", layer = layer_name, features = features)
        layer_mat[features, , drop = FALSE]
      })
      dat <- do.call(cbind, layer_mats)
      missing_cells <- setdiff(colnames(obj), colnames(dat))
      if (length(missing_cells) > 0) {
        stop("Data layers did not cover all cells. Missing cells: ", length(missing_cells))
      }
      return(dat[features, colnames(obj), drop = FALSE])
    }
    if (length(data_layers) == 1) {
      dat <- SeuratObject::LayerData(obj, assay = "RNA", layer = data_layers[1], features = features)
      return(dat[features, colnames(obj), drop = FALSE])
    }
    dat <- suppressWarnings({
      tryCatch(
        GetAssayData(obj, assay = "RNA", layer = "data"),
        error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
      )
    })
    if (nrow(dat) == 0 || ncol(dat) == 0) {
      message("RNA data layer is empty; running NormalizeData before scoring.")
      obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
      dat <- suppressWarnings({
        tryCatch(
          GetAssayData(obj, assay = "RNA", layer = "data"),
          error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
        )
      })
    }
    dat[features, colnames(obj), drop = FALSE]
  }

  match_gene_symbols <- function(requested_genes, available_genes) {
    available_upper <- toupper(available_genes)
    matched <- vapply(requested_genes, function(gene) {
      hit <- which(available_upper == toupper(gene))
      if (length(hit) == 0) return(NA_character_)
      available_genes[hit[1]]
    }, FUN.VALUE = character(1))
    matched
  }

  zscore_rows <- function(mat) {
    mat_z <- t(scale(t(mat)))
    mat_z[!is.finite(mat_z)] <- 0
    mat_z
  }

  cluster_order_within_group <- function(mat, row_info) {
    ordered_ids <- character()
    for (block in unique(row_info$row_block)) {
      block_ids <- row_info$feature_id[row_info$row_block == block]
      module_ids <- row_info$feature_id[row_info$row_block == block & row_info$feature_type == "module"]
      gene_ids <- setdiff(block_ids, module_ids)
      if (length(gene_ids) > 2) {
        gene_mat <- mat[gene_ids, , drop = FALSE]
        gene_ids <- gene_ids[stats::hclust(stats::dist(gene_mat))$order]
      }
      ordered_ids <- c(ordered_ids, module_ids, gene_ids)
    }
    ordered_ids
  }

  gene_matches <- match_gene_symbols(all_requested_genes, rownames(merged_obj))
  gene_availability <- data.frame(
    requested_gene = names(gene_matches),
    matched_gene = unname(gene_matches),
    detected = !is.na(gene_matches),
    group = vapply(names(gene_matches), function(gene) {
      group_hits <- names(stiffness_gene_groups)[vapply(stiffness_gene_groups, function(x) gene %in% x, FUN.VALUE = logical(1))]
      paste(setdiff(group_hits, "Overall tumour stiffness module"), collapse = "; ")
    }, FUN.VALUE = character(1)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    gene_availability,
    file.path(out_tiers$tables, "Auto_parse_tumour_stiffness_gene_availability.csv"),
    row.names = FALSE
  )

  matched_genes <- stats::na.omit(unname(gene_matches))
  if (length(matched_genes) == 0) {
    stop("None of the requested stiffness genes were found in the RNA assay.")
  }
  missing_genes <- gene_availability$requested_gene[!gene_availability$detected]

  log_data <- get_log_data(merged_obj, matched_genes)
  expr_subset <- log_data[matched_genes, , drop = FALSE]
  group_mm <- Matrix::sparse.model.matrix(~ 0 + sample_factor)
  colnames(group_mm) <- levels(sample_factor)
  avg_gene_expr <- as.matrix(expr_subset %*% group_mm)
  avg_gene_expr <- sweep(avg_gene_expr, 2, as.numeric(cell_counts[colnames(avg_gene_expr)]), "/")
  avg_gene_expr <- avg_gene_expr[, sample_order, drop = FALSE]

  detection_gene <- as.matrix((expr_subset > 0) %*% group_mm)
  detection_gene <- sweep(detection_gene, 2, as.numeric(cell_counts[colnames(detection_gene)]), "/")
  detection_gene <- detection_gene[, sample_order, drop = FALSE]

  feature_info <- data.frame(
    feature_id = "module__overall",
    display_name = "Module: all stiffness genes",
    row_block = "Overall tumour stiffness module",
    feature_type = "module",
    source_genes = paste(matched_genes, collapse = ";"),
    stringsAsFactors = FALSE
  )

  raw_rows <- list("module__overall" = colMeans(avg_gene_expr[matched_genes, , drop = FALSE], na.rm = TRUE))
  detection_rows <- list("module__overall" = colMeans(detection_gene[matched_genes, , drop = FALSE], na.rm = TRUE))

  for (block in setdiff(names(stiffness_gene_groups), "Overall tumour stiffness module")) {
    requested_block_genes <- stiffness_gene_groups[[block]]
    matched_block_genes <- stats::na.omit(unname(gene_matches[requested_block_genes]))
    module_id <- paste0("module__", gsub("[^A-Za-z0-9]+", "_", tolower(block)))
    if (length(matched_block_genes) > 0) {
      raw_rows[[module_id]] <- colMeans(avg_gene_expr[matched_block_genes, , drop = FALSE], na.rm = TRUE)
      detection_rows[[module_id]] <- colMeans(detection_gene[matched_block_genes, , drop = FALSE], na.rm = TRUE)
    } else {
      raw_rows[[module_id]] <- rep(NA_real_, length(sample_order))
      detection_rows[[module_id]] <- rep(NA_real_, length(sample_order))
    }
    names(raw_rows[[module_id]]) <- sample_order
    names(detection_rows[[module_id]]) <- sample_order
    feature_info <- bind_rows(
      feature_info,
      data.frame(
        feature_id = module_id,
        display_name = "Module: group mean",
        row_block = block,
        feature_type = "module",
        source_genes = paste(matched_block_genes, collapse = ";"),
        stringsAsFactors = FALSE
      )
    )
    for (gene in requested_block_genes) {
      matched_gene <- unname(gene_matches[gene])
      if (is.na(matched_gene)) next
      gene_id <- paste0("gene__", matched_gene)
      raw_rows[[gene_id]] <- avg_gene_expr[matched_gene, sample_order]
      detection_rows[[gene_id]] <- detection_gene[matched_gene, sample_order]
      feature_info <- bind_rows(
        feature_info,
        data.frame(
          feature_id = gene_id,
          display_name = matched_gene,
          row_block = block,
          feature_type = "gene",
          source_genes = matched_gene,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  raw_score_mat <- do.call(rbind, raw_rows)
  raw_score_mat <- raw_score_mat[feature_info$feature_id, sample_order, drop = FALSE]
  detection_mat <- do.call(rbind, detection_rows)
  detection_mat <- detection_mat[feature_info$feature_id, sample_order, drop = FALSE]
  z_score_mat <- zscore_rows(raw_score_mat)
  row_order <- cluster_order_within_group(z_score_mat, feature_info)

  feature_info <- feature_info[match(row_order, feature_info$feature_id), , drop = FALSE]
  raw_score_mat <- raw_score_mat[row_order, , drop = FALSE]
  detection_mat <- detection_mat[row_order, , drop = FALSE]
  z_score_mat <- z_score_mat[row_order, , drop = FALSE]

  score_long <- as.data.frame(raw_score_mat) |>
    rownames_to_column("feature_id") |>
    pivot_longer(cols = all_of(sample_order), names_to = "sample", values_to = "mean_log_normalized_score") |>
    left_join(feature_info, by = "feature_id") |>
    left_join(
      as.data.frame(z_score_mat) |>
        rownames_to_column("feature_id") |>
        pivot_longer(cols = all_of(sample_order), names_to = "sample", values_to = "z_score"),
      by = c("feature_id", "sample")
    ) |>
    left_join(
      as.data.frame(detection_mat) |>
        rownames_to_column("feature_id") |>
        pivot_longer(cols = all_of(sample_order), names_to = "sample", values_to = "mean_detection_fraction"),
      by = c("feature_id", "sample")
    ) |>
    mutate(
      sample = factor(sample, levels = sample_order),
      cell_n = as.integer(cell_counts[as.character(sample)])
    ) |>
    arrange(match(feature_id, row_order), sample)

  utils::write.csv(
    score_long,
    file.path(out_tiers$tables, "Auto_parse_tumour_stiffness_sample_scores.csv"),
    row.names = FALSE
  )

  saveRDS(
    list(
      raw_score_mat = raw_score_mat,
      z_score_mat = z_score_mat,
      detection_mat = detection_mat,
      feature_info = feature_info,
      gene_availability = gene_availability,
      sample_order = sample_order,
      cell_counts = cell_counts
    ),
    matrices_path
  )
}

####################
# heatmap
####################

group_cols <- c(
  "Overall tumour stiffness module" = "#3B3B3B",
  "Core mechanosensing / membrane tension" = "#3182BD",
  "Integrin-focal adhesion-actin linkage" = "#33B5A5",
  "Rho GTPase cytoskeletal regulators" = "#E28E2C",
  "Cell polarity / cortical organisation" = "#7A5195",
  "Actin cytoskeleton and contractility" = "#D24B40"
)

####################
hm_clip <- max(1, stats::quantile(abs(z_score_mat), 0.98, na.rm = TRUE))
hm_col <- circlize::colorRamp2(c(-hm_clip, 0, hm_clip), c("#2166AC", "#F7F7F7", "#B2182B"))

row_labels <- feature_info$display_name
names(row_labels) <- feature_info$feature_id
row_blocks <- factor(feature_info$row_block, levels = names(group_cols))

build_stiffness_heatmap <- function(samples_use, column_fontsize = 10) {
  z_mat_use <- z_score_mat[, samples_use, drop = FALSE]
  raw_mat_use <- raw_score_mat[, samples_use, drop = FALSE]

  column_rot <- if (any(nchar(samples_use) > 6)) 45 else 0

  row_anno <- rowAnnotation(
    Group = feature_info$row_block,
    col = list(Group = group_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = unit(4, "mm")
  )

  # Column names are handled at both the top (via top_anno Samples) and additionally at the bottom
  # (via show_column_names = TRUE and column_names_side = "bottom").
  top_anno <- HeatmapAnnotation(
    Cells = anno_barplot(
      as.numeric(cell_counts[samples_use]),
      gp = gpar(fill = "#787878", col = NA),
      height = unit(13, "mm"),
      border = FALSE,
      axis_param = list(gp = gpar(fontsize = 8))
    ),
    Samples = anno_text(
      samples_use,
      rot = column_rot,
      just = "center",
      gp = gpar(fontsize = column_fontsize, fontface = "bold")
    ),
    annotation_height = unit(c(13, 6), "mm"),
    gap = unit(2, "mm"),
    show_annotation_name = FALSE
  )

  cell_fun_scores <- function(j, i, x, y, w, h, fill) {
    grid.text(
      sprintf("%.2f", raw_mat_use[i, j]),
      x,
      y,
      gp = gpar(fontsize = 8.0, col = "black")
    )
    if (feature_info$feature_type[i] == "module") {
      grid.lines(
        c(x - 0.5 * w, x + 0.5 * w),
        c(y + 0.5 * h, y + 0.5 * h),
        gp = gpar(col = "black", lwd = 1.2)
      )
      grid.lines(
        c(x - 0.5 * w, x + 0.5 * w),
        c(y - 0.5 * h, y - 0.5 * h),
        gp = gpar(col = "black", lwd = 1.2)
      )
    }
  }

  wrapped_row_titles <- c(
    "Overall tumour stiffness module" = "Overall tumour\nstiffness module",
    "Core mechanosensing / membrane tension" = "Core mechanosensing\nmembrane tension",
    "Integrin-focal adhesion-actin linkage" = "Integrin-focal adhesion-actin\nlinkage",
    "Rho GTPase cytoskeletal regulators" = "Rho GTPase\ncytoskeletal regulators",
    "Cell polarity / cortical organisation" = "Cell polarity\ncortical organisation",
    "Actin cytoskeleton and contractility" = "Actin cytoskeleton\nand contractility"
  )

  ####################
  Heatmap(
    z_mat_use,
    name = "Row z-score",
    col = hm_col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_side = "bottom",
    column_names_rot = column_rot,
    column_names_gp = gpar(fontsize = column_fontsize, fontface = "bold"),
    row_labels = row_labels,
    row_names_gp = gpar(
      fontsize = 10,
      fontface = ifelse(feature_info$feature_type == "module", "bold", "plain")
    ),
    row_names_rot = 0,
    row_split = row_blocks,
    row_title = wrapped_row_titles[levels(row_blocks)],
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 9.5, fontface = "bold"),
    row_gap = unit(2.5, "mm"),
    left_annotation = row_anno,
    top_annotation = top_anno,
    cell_fun = cell_fun_scores,
    heatmap_legend_param = list(
      title = "Expression\nz-score",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )
  ####################
}

pdf_path_full <- file.path(out_tiers$figures, "Auto_parse_tumour_stiffness_module_gene_heatmap.pdf")
message("Writing heatmap PDF: ", pdf_path_full)
grDevices::cairo_pdf(pdf_path_full, width = 13.5, height = 11, family = "Arial")
draw(
  build_stiffness_heatmap(sample_order, column_fontsize = 10),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE,
  merge_legend = TRUE
)
dev.off()

sample_order_sub <- sample_order[!sample_order %in% c("PDO", "SUR1090_Untreated", "SUR1090_Treated")]
pdf_path_sub <- file.path(out_tiers$figures, "Auto_parse_tumour_stiffness_module_gene_heatmap_no_pdo_sur1090.pdf")
message("Writing heatmap PDF: ", pdf_path_sub)
grDevices::cairo_pdf(pdf_path_sub, width = 11.5, height = 11, family = "Arial")
draw(
  build_stiffness_heatmap(sample_order_sub, column_fontsize = 10),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE,
  merge_legend = TRUE
)
dev.off()
####################

####################
# Nature-style boxplot: stiffness gene scores across T0, T1, T4, eR4
####################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
})

# ---- Nature contract ---------------------------------------------------
# Claim: Tumour-stiffness gene expression shifts across treatment timepoints.
# Evidence: Per-gene mean log-normalised expression as individual data points,
#           with pairwise Wilcoxon rank-sum tests across T0, T1, T4, eR4.
# Archetype: quantitative grid (single panel).
# Export: PDF via cairo_pdf, Arial, Nature 6.5 pt base.
# -------------------------------------------------------------------------

theme_nature_contract <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line        = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.35, colour = "black"),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 0.5),
      legend.title     = element_text(size = base_size - 0.3),
      legend.text      = element_text(size = base_size - 0.7),
      strip.text       = element_text(size = base_size - 0.3, face = "bold"),
      plot.title       = element_text(size = base_size + 0.5, face = "bold"),
      panel.grid       = element_blank()
    )
}

# Subset to the four requested timepoints
boxplot_samples <- c("T0", "T1", "T4", "eR4")
boxplot_samples <- intersect(boxplot_samples, sample_order)

if (length(boxplot_samples) >= 2) {

  # Build long data: one row per gene × sample (gene-level features only)
  gene_features <- feature_info$feature_id[feature_info$feature_type == "gene"]
  gene_display  <- setNames(feature_info$display_name[feature_info$feature_type == "gene"],
                            gene_features)

  box_mat <- raw_score_mat[gene_features, boxplot_samples, drop = FALSE]

  box_df <- as.data.frame(box_mat) |>
    tibble::rownames_to_column("feature_id") |>
    tidyr::pivot_longer(
      cols      = all_of(boxplot_samples),
      names_to  = "sample",
      values_to = "mean_log_expr"
    ) |>
    dplyr::mutate(
      gene   = gene_display[feature_id],
      sample = factor(sample, levels = boxplot_samples)
    )

  # All pairwise comparisons
  pair_list <- combn(boxplot_samples, 2, simplify = FALSE)

  # Timepoint colours from config
  tp_cols <- parse_sample_colours[boxplot_samples]

  p_box <- ggplot(box_df, aes(x = sample, y = mean_log_expr, fill = sample)) +
    geom_boxplot(
      width       = 0.55,
      outlier.shape = NA,
      linewidth   = 0.3,
      alpha       = 0.65,
      colour      = "black"
    ) +
    geom_jitter(
      width  = 0.18,
      size   = 0.7,
      alpha  = 0.75,
      colour = "grey25",
      stroke = 0
    ) +
    stat_compare_means(
      comparisons   = pair_list,
      method        = "wilcox.test",
      label         = "p.signif",
      size          = 2.2,
      bracket.size  = 0.3,
      step.increase = 0.07,
      tip.length    = 0.015
    ) +
    scale_fill_manual(values = tp_cols) +
    labs(
      x     = NULL,
      y     = "Mean log-normalised expression",
      title = "Tumour stiffness gene expression across treatment timepoints"
    ) +
    theme_nature_contract() +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 7, face = "bold")
    )

  # Export — Nature single-column width ≈ 89 mm; double-column ≈ 183 mm
  boxplot_pdf <- file.path(out_tiers$figures,
                           "Auto_parse_tumour_stiffness_gene_boxplot_T0_T1_T4_eR4.pdf")
  message("Writing stiffness boxplot PDF: ", boxplot_pdf)
  grDevices::cairo_pdf(boxplot_pdf,
                       width  = 89 / 25.4,
                       height = 100 / 25.4,
                       family = "Arial")
  print(p_box)
  dev.off()

  message("Boxplot saved: ", boxplot_pdf)
} else {
  message("Fewer than 2 of the requested boxplot samples (T0, T1, T4, eR4) found — skipping boxplot.")
}
####################

script_run_status <- "success"
parse_finish_run(script_run, status = script_run_status)
script_run_finished <- TRUE
