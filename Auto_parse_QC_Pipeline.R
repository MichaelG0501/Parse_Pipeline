suppressPackageStartupMessages({
  library("Seurat")
  library("DoubletFinder")
  library("patchwork")
  library("ggplot2")
  library("gridExtra")
  library("grid")
  library("dplyr")
  library("tibble")
  library("Matrix")
  library("data.table")
  library("parallel")
})

####################
# Constants and Helper Functions (Merged from Auto_parse_helpers.R)
sample_order <- c("PDO", "T0", "T1", "T2", "T4", "R4", "eR4", "SUR1090_Treated", "SUR1090_Untreated")

max_mt <- 15
min_ngenes <- 2500
max_ngenes <- 13000
min_hk_expr <- 3
n_pcs <- 30

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

read_parse_sample <- function(sample_dir, sample_name) {
  meta <- fread(file.path(sample_dir, "cell_metadata.csv.gz"))
  genes <- fread(file.path(sample_dir, "all_genes.csv.gz"))
  counts <- t(readMM(gzfile(file.path(sample_dir, "count_matrix.mtx.gz"))))
  rownames(counts) <- make.unique(genes$gene_name)
  colnames(counts) <- make.unique(meta$bc_wells)

  meta <- as.data.frame(meta)
  rownames(meta) <- colnames(counts)
  meta$orig.ident <- sample_name
  meta$sample <- sample_name

  list(
    counts = counts,
    meta = meta,
    genes = genes
  )
}

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

get_doublet_rate <- function(ncells) {
  if (ncells <= 1000) {
    return(0.008)
  }
  if (ncells <= 5000) {
    return(0.04)
  }
  if (ncells <= 10000) {
    return(0.08)
  }
  if (ncells <= 20000) {
    return(0.16)
  }
  if (ncells <= 30000) {
    return(0.24)
  }
  0.24
}
####################

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline/Auto_parse_QC_Pipeline.R")
}

project_dir <- dirname(script_path)
out_dir <- file.path(project_dir, "parse_outs")
input_dir <- file.path(out_dir, "input", "output_combined")
setwd(out_dir)



if (!dir.exists(input_dir)) {
  stop("Missing extracted input under parse_outs/input/output_combined. Run the PBS wrapper first.")
}

sample_dirs <- setNames(file.path(input_dir, sample_order, "DGE_filtered"), sample_order)
missing_samples <- names(sample_dirs)[!dir.exists(sample_dirs)]
if (length(missing_samples) > 0) {
  stop("Missing sample directories: ", paste(missing_samples, collapse = ", "))
}

dir.create("plots/Filtering_and_Summary", showWarnings = FALSE, recursive = TRUE)
dir.create("plots/Visualisation_and_Heatmaps", showWarnings = FALSE, recursive = TRUE)
dir.create("summary", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

inspect_samples <- function(tmdata_list) {
  x_features_plot <- list()
  x_count_plot <- list()
  x_mito_plot <- list()

  for (name in names(tmdata_list)) {
    obj <- tmdata_list[[name]]
    meta_df <- obj@meta.data
    meta_df$orig.ident <- factor(meta_df$orig.ident, levels = name)

    mean_nfeature <- mean(meta_df$nFeature_RNA, na.rm = TRUE)
    median_nfeature <- median(meta_df$nFeature_RNA, na.rm = TRUE)
    mean_ncount <- mean(meta_df$nCount_RNA, na.rm = TRUE)
    median_ncount <- median(meta_df$nCount_RNA, na.rm = TRUE)
    mean_percent_mt <- mean(meta_df$percent.mt, na.rm = TRUE)
    median_percent_mt <- median(meta_df$percent.mt, na.rm = TRUE)

    base_theme <- theme(
      text = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 6),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 6),
      legend.position = "none"
    )

    x_features_plot[[name]] <- ggplot(meta_df, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_hline(yintercept = mean_nfeature, linetype = "dashed", color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = median_nfeature, linetype = "solid", color = "red", linewidth = 0.5) +
      annotate("text", x = 1, y = mean_nfeature, label = paste("Mean:", round(mean_nfeature, 1)),
               hjust = 0.5, vjust = -1, size = 3, color = "blue") +
      annotate("text", x = 1, y = median_nfeature, label = paste("Median:", round(median_nfeature, 1)),
               hjust = 0.5, vjust = 1.5, size = 3, color = "red") +
      annotate("text", x = Inf, y = Inf, label = paste("NCells:", ncol(obj)),
               hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
      labs(title = name, y = "nFeature_RNA") +
      base_theme

    x_count_plot[[name]] <- ggplot(meta_df, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_hline(yintercept = mean_ncount, linetype = "dashed", color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = median_ncount, linetype = "solid", color = "red", linewidth = 0.5) +
      annotate("text", x = 1, y = mean_ncount, label = paste("Mean:", round(mean_ncount, 1)),
               hjust = 0.5, vjust = -1, size = 3, color = "blue") +
      annotate("text", x = 1, y = median_ncount, label = paste("Median:", round(median_ncount, 1)),
               hjust = 0.5, vjust = 1.5, size = 3, color = "red") +
      annotate("text", x = Inf, y = Inf, label = paste("NCells:", ncol(obj)),
               hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
      labs(title = name, y = "nCount_RNA") +
      base_theme

    x_mito_plot[[name]] <- ggplot(meta_df, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_hline(yintercept = mean_percent_mt, linetype = "dashed", color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = median_percent_mt, linetype = "solid", color = "red", linewidth = 0.5) +
      annotate("text", x = 1, y = mean_percent_mt, label = paste("Mean:", round(mean_percent_mt, 1)),
               hjust = 0.5, vjust = -1, size = 3, color = "blue") +
      annotate("text", x = 1, y = median_percent_mt, label = paste("Median:", round(median_percent_mt, 1)),
               hjust = 0.5, vjust = 1.5, size = 3, color = "red") +
      annotate("text", x = Inf, y = Inf, label = paste("NCells:", ncol(obj)),
               hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
      labs(title = name, y = "percent.mt") +
      base_theme
  }

  pdf(file.path("plots/Filtering_and_Summary", "qc_inspection_violin.pdf"), width = 16, height = 9)
  for (plot_list in list(x_features_plot, x_count_plot, x_mito_plot)) {
    # Layout 5 columns x 2 rows = 10 slots (fits all 9 samples)
    chunks <- split(plot_list, ceiling(seq_along(plot_list) / 10))
    for (chunk in chunks) {
      grid.arrange(grobs = chunk, ncol = 5, nrow = 2)
    }
  }
  dev.off()
}

run_doublet_filter <- function(obj, sample_name) {
  options(future.globals.maxSize = 8 * 1024^3)
  raw_n <- ncol(obj)

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:n_pcs, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:n_pcs, verbose = FALSE)

  sweep_obj <- paramSweep(obj, PCs = 1:n_pcs, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_obj)
  bcmvn <- find.pK(sweep_stats)
  pk_value <- suppressWarnings(as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])))
  if (!is.finite(pk_value)) {
    pk_value <- 0.01
  }

  dr <- get_doublet_rate(ncol(obj))
  homotypic <- modelHomotypic(obj$seurat_clusters)
  n_exp <- round(dr * ncol(obj))
  n_exp_adj <- max(1, round(n_exp * (1 - homotypic)))

  obj <- doubletFinder(
    obj,
    PCs = 1:n_pcs,
    pN = 0.25,
    pK = pk_value,
    nExp = n_exp_adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  class_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)
  if (length(class_col) == 0) {
    stop("DoubletFinder classification column missing for ", sample_name)
  }
  class_col <- class_col[1]

  plot_obj <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    DimPlot(obj, reduction = "umap", group.by = class_col) +
    plot_annotation(title = sample_name)

  singlets <- rownames(obj@meta.data)[obj@meta.data[[class_col]] == "Singlet"]
  obj <- subset(obj, cells = singlets)

  summary_row <- tibble(
    sample = sample_name,
    raw = raw_n,
    doublet_singlets = length(singlets),
    pK = pk_value,
    doublet_rate = dr,
    expected_doublets = n_exp,
    expected_doublets_adj = n_exp_adj
  )

  list(obj = obj, plot = plot_obj, summary = summary_row)
}

plot_umap_meta <- function(obj, column, title) {
  umap_df <- as.data.frame(Embeddings(obj, "umap"))
  colnames(umap_df)[seq_len(min(2, ncol(umap_df)))] <- c("UMAP_1", "UMAP_2")[seq_len(min(2, ncol(umap_df)))]
  umap_df[[column]] <- obj@meta.data[[column]]

  ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[column]])) +
    geom_point(size = 0.15) +
    scale_color_viridis_c() +
    labs(title = title, color = column) +
    theme_classic(base_size = 10) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

plot_umap_gene <- function(obj, gene) {
  umap_df <- as.data.frame(Embeddings(obj, "umap"))
  colnames(umap_df)[seq_len(min(2, ncol(umap_df)))] <- c("UMAP_1", "UMAP_2")[seq_len(min(2, ncol(umap_df)))]
  umap_df$expr <- FetchData(obj, vars = gene)[, 1]

  ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = expr)) +
    geom_point(size = 0.15) +
    scale_color_viridis_c() +
    labs(title = gene, color = "expr") +
    theme_classic(base_size = 10) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

main <- function() {
  raw_list <- list()
  manifest_rows <- list()

  for (sample_name in sample_order) {
    parsed <- read_parse_sample(sample_dirs[[sample_name]], sample_name)
    obj <- CreateSeuratObject(counts = parsed$counts, meta.data = parsed$meta)
    obj$orig.ident <- sample_name
    obj$sample <- sample_name
    obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
    raw_list[[sample_name]] <- obj

    manifest_rows[[sample_name]] <- tibble(
      sample = sample_name,
      sample_dir = sample_dirs[[sample_name]],
      raw_cells = ncol(obj)
    )
    
    # Save truly raw for heatmap downstream
    sample_dir <- file.path("by_samples", sample_name)
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(obj, file.path(sample_dir, paste0("Auto_", sample_name, "_raw.rds")), compress = FALSE)
  }

  manifest_df <- bind_rows(manifest_rows)
  write.csv(manifest_df, file.path("summary", "Auto_parse_sample_manifest.csv"), row.names = FALSE)

  inspect_samples(raw_list)

  doublet_plots <- list()
  doublet_rows <- list()
  singlet_list <- list()

  for (sample_name in sample_order) {
    res <- run_doublet_filter(raw_list[[sample_name]], sample_name)
    singlet_list[[sample_name]] <- res$obj
    doublet_plots[[sample_name]] <- res$plot
    doublet_rows[[sample_name]] <- res$summary
    gc()
  }

  pdf(file.path("plots/Filtering_and_Summary", "doublet_filtering.pdf"), width = 12, height = 8)
  for (sample_name in sample_order) {
    print(doublet_plots[[sample_name]])
  }
  dev.off()

  doublet_df <- bind_rows(doublet_rows)
  write.csv(doublet_df, file.path("summary", "Auto_parse_doublet_summary.csv"), row.names = FALSE)

  final_list <- list()
  filter_plots <- list()
  summary_rows <- list()

  for (sample_name in sample_order) {
    obj <- singlet_list[[sample_name]]
    sample_dir <- file.path("by_samples", sample_name)
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
    raw_n <- manifest_df$raw_cells[manifest_df$sample == sample_name]
    singlet_n <- ncol(obj)

    obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj <- subset(obj, subset = percent.mt < max_mt)
    mito_n <- ncol(obj)

    saveRDS(obj, file.path(sample_dir, paste0("Auto_", sample_name, "_prefilter.rds")), compress = FALSE)

    hk_mean <- compute_hk_mean(GetAssayData(obj, layer = "counts"))
    obj$hk_mean <- hk_mean

    keep_gene <- obj$nFeature_RNA >= min_ngenes & obj$nFeature_RNA <= max_ngenes
    gene_n <- sum(keep_gene)
    keep_hk <- keep_gene & obj$hk_mean >= min_hk_expr
    hk_n <- sum(keep_hk)

    plot_df <- data.frame(
      n_genes = obj$nFeature_RNA,
      hk_mean = obj$hk_mean,
      keep = keep_hk
    )

    filter_plots[[sample_name]] <- ggplot(plot_df, aes(x = n_genes, y = hk_mean, color = keep)) +
      geom_point(size = 0.5) +
      scale_x_continuous(trans = "log10", labels = scales::comma) +
      scale_y_continuous(trans = "log10", labels = scales::comma) +
      scale_color_manual(values = c("FALSE" = "lightgrey", "TRUE" = "black")) +
      geom_vline(xintercept = min_ngenes, linetype = "dashed", color = "red") +
      geom_vline(xintercept = max_ngenes, linetype = "dashed", color = "red") +
      geom_hline(yintercept = min_hk_expr, linetype = "dashed", color = "red") +
      annotate("text", x = Inf, y = Inf, label = paste0("NCells passed: ", hk_n),
               hjust = 1.1, vjust = 1.1, size = 3) +
      labs(title = sample_name, x = "Number of genes", y = "HK mean", color = "keep") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none")

    if (hk_n == 0) {
      stop("No cells passed QC for ", sample_name)
    }

    obj <- subset(obj, cells = colnames(obj)[keep_hk])
    obj <- NormalizeData(obj, verbose = FALSE)
    final_list[[sample_name]] <- obj
    saveRDS(obj, file.path(sample_dir, paste0("Auto_", sample_name, "_final.rds")), compress = FALSE)

    summary_rows[[sample_name]] <- tibble(
      sample = sample_name,
      raw = raw_n,
      doublet_singlets = singlet_n,
      `mito_DNA\npercentage < 15` = mito_n,
      `number of\ngenes` = gene_n,
      `housekeeping\nexpression > 3` = hk_n
    )
  }

  pdf(file.path("plots/Filtering_and_Summary", "cells_filtering.pdf"), width = 16, height = 9)
  # Layout 5 columns x 2 rows = 10 slots (fits all 9 samples)
  chunks <- split(filter_plots, ceiling(seq_along(filter_plots) / 10))
  for (chunk in chunks) {
    grid.arrange(grobs = chunk, ncol = 5, nrow = 2)
  }
  dev.off()

  summary_df <- bind_rows(summary_rows)
  write.csv(summary_df, file.path("summary", "Auto_parse_filtering_summary.csv"), row.names = FALSE)
  write.csv(
    data.frame(
      max_mt = max_mt,
      min_ngenes = min_ngenes,
      max_ngenes = max_ngenes,
      min_hk_expr = min_hk_expr,
      stringsAsFactors = FALSE
    ),
    file.path("summary", "Auto_parse_thresholds.csv"),
    row.names = FALSE
  )

  merged_obj <- merge(final_list[[1]], y = final_list[-1], add.cell.ids = sample_order, project = "parse_qc")
  merged_obj <- NormalizeData(merged_obj, verbose = FALSE)
  merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  merged_obj <- ScaleData(merged_obj, verbose = FALSE)
  merged_obj <- RunPCA(merged_obj, npcs = n_pcs, verbose = FALSE)
  merged_obj <- FindNeighbors(merged_obj, dims = 1:n_pcs, verbose = FALSE)
  merged_obj <- FindClusters(merged_obj, resolution = 0.8, verbose = FALSE)
  merged_obj <- RunUMAP(merged_obj, dims = 1:n_pcs, verbose = FALSE)

  saveRDS(merged_obj, "Auto_parse_merged.rds", compress = FALSE)
  saveRDS(merged_obj@meta.data, "Auto_parse_all_meta.rds", compress = FALSE)

  p1 <- DimPlot(merged_obj, group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
  p2 <- DimPlot(merged_obj, group.by = "orig.ident") + ggtitle("Sample")
  p3 <- plot_umap_meta(merged_obj, "nFeature_RNA", "nFeature_RNA")
  p4 <- plot_umap_meta(merged_obj, "percent.mt", "percent.mt")
  ggsave(file.path("plots/Visualisation_and_Heatmaps", "Auto_parse_merged_umap.png"), (p1 | p2) / (p3 | p4), width = 12, height = 8, dpi = 300)

  feature_genes <- intersect(c("EPCAM", "KRT19", "KRT7", "COL1A1", "PTPRC", "VWF"), rownames(merged_obj))
  if (length(feature_genes) > 0) {
    p_gene <- wrap_plots(lapply(feature_genes, function(gene) plot_umap_gene(merged_obj, gene)), ncol = 3)
    ggsave(file.path("plots/Visualisation_and_Heatmaps", "Auto_parse_marker_featureplots.png"), p_gene, width = 12, height = 8, dpi = 300)
  }

  write.csv(
    merged_obj@meta.data %>%
      count(orig.ident, name = "final_cells") %>%
      rename(sample = orig.ident) %>%
      mutate(sample = factor(sample, levels = sample_order)) %>%
      arrange(sample) %>%
      mutate(sample = as.character(sample)),
    file.path("summary", "Auto_parse_filtered_sample_summary.csv"),
    row.names = FALSE
  )
}

main()
