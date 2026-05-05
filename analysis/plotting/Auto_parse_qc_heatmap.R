suppressPackageStartupMessages({
  library("Seurat")
  library("dplyr")
  library("ComplexHeatmap")
  library("circlize")
  library("gridExtra")
  library("grid")
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline/analysis/plotting/Auto_parse_qc_heatmap.R")
}

project_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))
out_dir <- file.path(project_dir, "parse_outs")
setwd(out_dir)

####################
# Constants and Helper Functions (Merged from Auto_parse_helpers.R)
sample_order <- c("PDO", "T0", "T1", "T2", "T4", "R4", "eR4", "SUR1090_Treated", "SUR1090_Untreated")

max_mt <- 15
min_ngenes <- 2500
max_ngenes <- 13000
min_hk_expr <- 3

housekeeping_genes <- c(
  "ACTB", "GAPDH", "RPS11", "RPS13", "RPS14", "RPS15", "RPS16", "RPS18",
  "RPS19", "RPS20", "RPL10", "RPL13", "RPL15", "RPL18"
)

marker_panel <- list(
  b.cell = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
  dendritic = c("CLEC10A", "CCR7", "CD86"),
  endothelial = c("ENG", "CLEC14A", "CLDN5", "VWF", "CDH5"),
  epithelial = c("KRT7", "MUC1", "KRT19", "EPCAM"),
  fibroblast = c("COL3A1", "COL1A2", "LUM", "COL1A1", "COL6A3", "DCN"),
  macrophage = c("CSF1R", "TYROBP", "CD14", "CD163", "AIF1", "CD68"),
  mast = c("MS4A2", "CPA3", "TPSB2", "TPSAB1"),
  nk.cell = c("GNLY", "NKG7", "PRF1", "GZMB", "KLRB1"),
  plasma = c("MZB1", "JCHAIN", "DERL3"),
  t.cell = c("CD3E", "CD3D", "CD2", "CD3G"),
  housekeeping = housekeeping_genes
)

marker_genes <- unique(unlist(marker_panel, use.names = FALSE))

compute_hk_mean <- function(counts) {
  totals <- Matrix::colSums(counts)
  hk_present <- intersect(housekeeping_genes, rownames(counts))
  if (length(hk_present) == 0) {
    return(rep(NA_real_, ncol(counts)))
  }
  hk_counts <- counts[hk_present, , drop = FALSE]
  hk_cpm <- t(t(hk_counts) / totals) * 1e6
  hk_cpm[is.na(hk_cpm)] <- 0
  Matrix::colMeans(log2(as.matrix(hk_cpm) / 10 + 1))
}

get_marker_expr <- function(counts) {
  keep_genes <- intersect(marker_genes, rownames(counts))
  expr <- matrix(
    0,
    nrow = length(marker_genes),
    ncol = ncol(counts),
    dimnames = list(marker_genes, colnames(counts))
  )
  if (length(keep_genes) == 0) {
    return(expr)
  }
  totals <- Matrix::colSums(counts)
  marker_counts <- counts[keep_genes, , drop = FALSE]
  marker_cpm <- t(t(marker_counts) / totals) * 1e6
  marker_cpm[is.na(marker_cpm)] <- 0
  expr[keep_genes, ] <- log2(as.matrix(marker_cpm) / 10 + 1)
  expr
}
####################



sample_manifest <- read.csv(file.path("summary", "Auto_parse_sample_manifest.csv"), stringsAsFactors = FALSE)
sample_manifest <- sample_manifest %>% filter(sample %in% sample_order)

build_stage <- function(kind = c("raw", "prefilter", "final")) {
  kind <- match.arg(kind)
  expr_parts <- list()
  meta_parts <- list()

  for (sample_name in sample_order) {
    file_name <- if (kind == "raw") {
      file.path("by_samples", sample_name, paste0("Auto_", sample_name, "_raw.rds"))
    } else if (kind == "prefilter") {
      file.path("by_samples", sample_name, paste0("Auto_", sample_name, "_prefilter.rds"))
    } else {
      file.path("by_samples", sample_name, paste0("Auto_", sample_name, "_final.rds"))
    }
    obj <- readRDS(file_name)
    counts <- GetAssayData(obj, layer = "counts")
    expr <- get_marker_expr(counts)
    plot_cell_id <- paste(sample_name, colnames(counts), sep = "___")
    colnames(expr) <- plot_cell_id
    expr_parts[[sample_name]] <- expr

    meta_parts[[sample_name]] <- data.frame(
      plot_cell_id = plot_cell_id,
      sample = sample_name,
      nFeature_RNA = obj$nFeature_RNA,
      percent.mt = obj$percent.mt,
      hk_mean = compute_hk_mean(counts),
      stringsAsFactors = FALSE
    )
  }

  meta_df <- bind_rows(meta_parts)
  rownames(meta_df) <- meta_df$plot_cell_id
  expr_mat <- do.call(cbind, expr_parts)
  expr_mat <- expr_mat[, rownames(meta_df), drop = FALSE]
  list(expr = expr_mat, meta = meta_df)
}

plot_heatmap <- function(stage_obj, sampleid) {
  expr_data <- as.data.frame(t(stage_obj$expr))
  expr_data$sample <- factor(stage_obj$meta$sample, levels = sample_order)
  expr_data <- expr_data[order(expr_data$sample), , drop = FALSE]

  hk_avg <- matrix(stage_obj$meta$hk_mean[match(rownames(expr_data), stage_obj$meta$plot_cell_id)], nrow = 1)
  colnames(hk_avg) <- rownames(expr_data)
  rownames(hk_avg) <- "avg_hk"

  ncounts <- matrix(stage_obj$meta$nFeature_RNA[match(rownames(expr_data), stage_obj$meta$plot_cell_id)], nrow = 1)
  colnames(ncounts) <- rownames(expr_data)
  rownames(ncounts) <- "nGenes"

  heatplot <- t(as.matrix(expr_data[, marker_genes, drop = FALSE]))
  expr_max <- max(heatplot, na.rm = TRUE)
  expr_colors <- colorRamp2(c(0, round(0.6 * expr_max, 1), ceiling(expr_max)), c("#D0D0D0", "red4", "red4"))

  present_groups <- intersect(sample_order, unique(as.character(expr_data$sample)))
  len <- length(present_groups)
  heatmap_grobs <- list()
  stats_grobs <- list()

  for (i in seq_along(marker_panel)) {
    marker <- marker_panel[[i]]
    marker <- marker[marker %in% rownames(heatplot)]
    if (length(marker) == 0) {
      next
    }

    temp <- list()
    for (j in seq_along(present_groups)) {
      group_name <- present_groups[j]
      cells <- rownames(expr_data)[expr_data$sample == group_name]
      ht <- Heatmap(
        heatplot[marker, cells, drop = FALSE],
        col = expr_colors,
        show_column_names = FALSE,
        show_row_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_heatmap_legend = FALSE,
        use_raster = TRUE
      )
      temp[[j]] <- grid.grabExpr(draw(ht, newpage = FALSE, padding = unit(c(2, 1, 2, 1), "mm")))
    }

    gene_label_col <- arrangeGrob(
      grobs = lapply(marker, function(name) textGrob(name, x = unit(1, "npc"), just = "right", gp = gpar(fontsize = 40))),
      ncol = 1
    )

    temp <- c(list(gene_label_col), temp)
    text_grob <- textGrob(names(marker_panel)[i], gp = gpar(fontsize = 50, fontface = "bold"))
    rect_grob <- rectGrob(gp = gpar(fill = "grey", col = NA))
    temp[[length(present_groups) + 2]] <- gTree(children = gList(rect_grob, text_grob))

    heatmap_grobs[[length(heatmap_grobs) + 1]] <- do.call(
      arrangeGrob,
      c(temp, list(ncol = len + 2, widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
    )
  }

  for (i in seq_along(list(ncounts, hk_avg))) {
    data_mat <- list(ncounts, hk_avg)[[i]]
    temp <- list()

    for (j in seq_along(present_groups)) {
      group_name <- present_groups[j]
      cells <- rownames(expr_data)[expr_data$sample == group_name]
      if (i == 1) {
        ht <- Heatmap(
          data_mat[, cells, drop = FALSE],
          col = colorRamp2(c(min_ngenes, round(0.7 * max(ncounts), -3), ceiling(max(ncounts) / 1000) * 1000), c("#D0D0D0", "blue3", "blue3")),
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          show_heatmap_legend = FALSE,
          use_raster = TRUE
        )
      } else {
        ht <- Heatmap(
          data_mat[, cells, drop = FALSE],
          col = colorRamp2(c(0, min_hk_expr, ceiling(max(hk_avg) * 10) / 10), c("#D0D0D0", "#A0B8E6", "blue3")),
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          show_heatmap_legend = FALSE,
          use_raster = TRUE
        )
      }
      temp[[j]] <- grid.grabExpr(draw(ht, newpage = FALSE, padding = unit(c(6, 1.5, 6, 1.5), "mm")))
    }

    gene_label_col <- arrangeGrob(
      grobs = lapply(rownames(data_mat), function(name) textGrob(name, just = "center", gp = gpar(fontsize = 40))),
      ncol = 1
    )
    temp <- c(list(gene_label_col), temp)
    label <- c("Number\nof genes", "Average\nhousekeeping\nexpression")[i]
    temp[[length(present_groups) + 2]] <- gTree(children = gList(rectGrob(gp = gpar(fill = "white", col = NA)), textGrob(label, gp = gpar(fontsize = 40, fontface = "bold"))))
    stats_grobs[[length(stats_grobs) + 1]] <- do.call(
      arrangeGrob,
      c(temp, list(ncol = len + 2, widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
    )
  }

  label_row <- do.call(
    arrangeGrob,
    c(lapply(c("", present_groups, ""), function(x) textGrob(x, gp = gpar(fontsize = 50, fontface = "bold"))),
      list(ncol = len + 2, widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
  )

  count_row <- do.call(
    arrangeGrob,
    c(lapply(c("", as.numeric(table(expr_data$sample)[present_groups]), "Markers"), function(x) textGrob(x, gp = gpar(fontsize = 40, fontface = "bold"))),
      list(ncol = len + 2, widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
  )

  title_grob <- textGrob(
    paste0(
      "Mitochondria DNA < ", max_mt,
      "     &&     Minimum number of genes > ", min_ngenes,
      "     &&     Maximum number of genes < ", max_ngenes,
      "     &&     Minimum HK expression > ", min_hk_expr
    ),
    gp = gpar(fontsize = 40)
  )

  expr_legend <- Legend(
    title = "\nExpression\nlevels (E)\n",
    at = pretty(c(0, ceiling(expr_max)), n = 5),
    labels_gp = gpar(fontsize = 50),
    title_gp = gpar(fontsize = 50, fontface = "bold"),
    grid_height = unit(220, "mm"),
    legend_width = unit(40, "mm"),
    legend_height = unit(220, "mm"),
    grid_width = unit(40, "mm"),
    title_position = "topcenter",
    col_fun = expr_colors
  )

  gene_high <- max(ceiling(max(ncounts) / 1000) * 1000, min_ngenes + 1000)
  hk_high <- max(ceiling(max(hk_avg) * 10) / 10, min_hk_expr)

  gene_legend <- Legend(
    title = "\nNumber\nof genes\n",
    at = round(seq(min_ngenes, gene_high, length.out = 4), -3),
    labels_gp = gpar(fontsize = 50),
    title_gp = gpar(fontsize = 50, fontface = "bold"),
    grid_height = unit(220, "mm"),
    legend_width = unit(40, "mm"),
    legend_height = unit(220, "mm"),
    grid_width = unit(40, "mm"),
    title_position = "topcenter",
    col_fun = colorRamp2(c(min_ngenes, round(0.7 * max(ncounts), -3), gene_high), c("#D0D0D0", "blue3", "blue3"))
  )

  hk_legend <- Legend(
    title = "\nAverage\nhousekeeping\nexpression\n",
    at = seq(floor(min_hk_expr), hk_high, length.out = 4) %>% round(1),
    labels_gp = gpar(fontsize = 50),
    title_gp = gpar(fontsize = 50, fontface = "bold"),
    grid_height = unit(220, "mm"),
    legend_width = unit(40, "mm"),
    legend_height = unit(220, "mm"),
    grid_width = unit(40, "mm"),
    title_position = "topcenter",
    col_fun = colorRamp2(c(0, min_hk_expr, hk_high), c("#D0D0D0", "#A0B8E6", "blue3"))
  )

  main_content <- arrangeGrob(
    grobs = c(list(title_grob), list(label_row), list(count_row), heatmap_grobs, stats_grobs),
    ncol = 1,
    heights = c(2.5, 2, 1.5, rep(4, length(heatmap_grobs)), 4, 4)
  )

  legend_column <- arrangeGrob(
    textGrob(sampleid, gp = gpar(fontsize = 35, fontface = "bold")),
    textGrob(paste0("Total cell count: ", nrow(expr_data)), gp = gpar(fontsize = 35)),
    grid.grabExpr(draw(expr_legend)),
    grid.grabExpr(draw(gene_legend)),
    grid.grabExpr(draw(hk_legend)),
    ncol = 1,
    heights = c(0.3, 0.2, 2, 2, 2)
  )

  grid.arrange(main_content, legend_column, ncol = 2, widths = c(9.3, 1))
}

raw_stage <- build_stage("raw")
final_stage <- build_stage("final")

summary_tbl <- bind_rows(
  raw_stage$meta %>% count(stage_name = "raw_no_filtering", sample, name = "n_plotted") %>% mutate(total_cells_stage = nrow(raw_stage$meta)),
  final_stage$meta %>% count(stage_name = "final_filtered_singlets", sample, name = "n_plotted") %>% mutate(total_cells_stage = nrow(final_stage$meta))
)

write.csv(summary_tbl, file.path("summary", "Auto_parse_qc_heatmap_stage_summary.csv"), row.names = FALSE)

png(file.path("plots/Visualisation_and_Heatmaps", "QC_prefilter.png"), width = 80, height = 50, units = "in", res = 300)
plot_heatmap(raw_stage, "Raw (No Filtering)")
dev.off()

png(file.path("plots/Visualisation_and_Heatmaps", "QC_final.png"), width = 80, height = 50, units = "in", res = 300)
plot_heatmap(final_stage, "After Doublet + MITO")
dev.off()

pdf(file.path("plots/Visualisation_and_Heatmaps", "QC_heatmaps.pdf"), width = 80, height = 50, onefile = TRUE, useDingbats = FALSE)
plot_heatmap(raw_stage, "Raw (No Filtering)")
grid.newpage()
plot_heatmap(final_stage, "After Doublet + MITO")
dev.off()
