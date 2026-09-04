suppressPackageStartupMessages({
  library("Seurat")
  library("patchwork")
  library("ggplot2")
  library("dplyr")
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/novogene/Auto_parse_qc_pipeline/Auto_parse_finalize_outputs.R")
}

project_dir <- dirname(script_path)
out_dir <- file.path(project_dir, "parse_qc_outs")
setwd(out_dir)

source(file.path(project_dir, "Auto_parse_helpers.R"))

dir.create("plots", showWarnings = FALSE, recursive = TRUE)
dir.create("summary", showWarnings = FALSE, recursive = TRUE)

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

merged_obj <- readRDS("merged_obj.rds")

p1 <- DimPlot(merged_obj, group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
p2 <- DimPlot(merged_obj, group.by = "orig.ident") + ggtitle("Sample")
p3 <- plot_umap_meta(merged_obj, "nFeature_RNA", "nFeature_RNA")
p4 <- plot_umap_meta(merged_obj, "percent.mt", "percent.mt")
ggsave(file.path("plots", "Auto_parse_merged_umap.png"), (p1 | p2) / (p3 | p4), width = 12, height = 8, dpi = 300)

feature_genes <- intersect(c("EPCAM", "KRT19", "KRT7", "COL1A1", "PTPRC", "VWF"), rownames(merged_obj))
if (length(feature_genes) > 0) {
  p_gene <- wrap_plots(lapply(feature_genes, function(gene) plot_umap_gene(merged_obj, gene)), ncol = 3)
  ggsave(file.path("plots", "Auto_parse_marker_featureplots.png"), p_gene, width = 12, height = 8, dpi = 300)
}

write.csv(
  merged_obj@meta.data %>%
    count(orig.ident, name = "final_cells") %>%
    rename(sample = orig.ident) %>%
    mutate(sample = factor(sample, levels = sample_order)) %>%
    arrange(sample) %>%
    mutate(sample = as.character(sample)),
  file.path("summary", "filtered_sample_summary.csv"),
  row.names = FALSE
)
