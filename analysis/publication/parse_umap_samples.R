####################
# parse_umap_samples.R
#
# Description:
#   Publication-quality UMAP visualization of the 6 Parse
#   treatment-response timepoints (T0, T1, T2, T4, R4, eR4) plus PDO.
#   Cells are colored by sample with direct text labels on cluster
#   centroids (no legend). Nature-contract figure style.
#
# Inputs:
#   parse_outs/Auto_parse_merged.rds
#
# Outputs:
#   parse_outs/publication/umap_samples/figures/
#     umap_samples.pdf
#     umap_samples.png
#
# Methodology:
#   analysis/methodology/publication/umap_samples_methodology.md
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

source("analysis/common/parse_pipeline_config.R")

# ── Paths ────────────────────────────────────────────────────────────────────
project_dir <- parse_project_root()
paths       <- parse_paths(project_dir)
out_base    <- file.path(paths$parse_outs, "publication", "umap_samples")
fig_dir     <- file.path(out_base, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading merged object...\n")
# ── Load data ────────────────────────────────────────────────────────────────
merged_path <- file.path(paths$parse_outs, "Auto_parse_merged.rds")
if (!file.exists(merged_path)) stop("Missing: ", merged_path)
merged_obj <- readRDS(merged_path)

cat("Subsetting...\n")
# ── Subset & Prepare data ────────────────────────────────────────────────────
target_samples <- c("PDO", "T0", "T1", "T2", "T4", "R4", "eR4")
sub_obj <- subset(merged_obj, subset = sample %in% target_samples)

cat("Recalculating UMAP...\n")
sub_obj <- RunUMAP(sub_obj, dims = 1:30, verbose = FALSE)

cat("Extracting UMAP...\n")
# Extract UMAP coordinates
umap_coords <- Embeddings(sub_obj, "umap")
plot_data <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  sample = as.character(sub_obj$sample)
)

# Randomize plot order to prevent overplotting bias
set.seed(42)
plot_data <- plot_data[sample(nrow(plot_data)), ]
plot_data$sample <- factor(plot_data$sample, levels = target_samples)

cat("Computing cluster centroids...\n")
# ── Centroid labels ──────────────────────────────────────────────────────────
centroids <- plot_data %>%
  group_by(sample) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

cat("Setting colours...\n")
# ── Colours ──────────────────────────────────────────────────────────────────
cols <- c(
  "PDO" = unname(parse_sample_colours["PDO"]),
  "T0"  = unname(parse_sample_colours["T0"]),
  "T1"  = unname(parse_sample_colours["T1"]),
  "T2"  = unname(parse_sample_colours["T2"]),
  "T4"  = unname(parse_sample_colours["T4"]),
  "R4"  = unname(parse_sample_colours["R4"]),
  "eR4" = unname(parse_sample_colours["eR4"])
)

cat("Building plot...\n")
# ── Nature-contract theme ────────────────────────────────────────────────────
theme_nature_umap <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line         = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks        = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = unit(0.8, "mm"),
      axis.text         = element_blank(),
      axis.title        = element_text(size = base_size, face = "bold"),
      legend.position   = "none",
      plot.title        = element_blank(),
      plot.margin       = margin(4, 6, 4, 4),
      panel.grid        = element_blank()
    )
}

# ── Build plot ───────────────────────────────────────────────────────────────
p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, colour = sample)) +
  # Background shadow layer for depth
  geom_point(size = 0.6, alpha = 0.15, colour = "grey30", stroke = 0) +
  # Main data points
  geom_point(size = 0.45, alpha = 0.9, stroke = 0) +
  # Direct labels at centroids
  geom_label_repel(
    data          = centroids,
    aes(label = sample, colour = sample),
    size          = 2.8,
    fontface      = "bold",
    fill          = alpha("white", 0.85),
    label.size    = 0.2,
    label.r       = unit(1.2, "mm"),
    label.padding = unit(1.5, "mm"),
    segment.size  = 0.3,
    segment.color = "grey50",
    box.padding   = unit(4, "mm"),
    point.padding = unit(2, "mm"),
    min.segment.length = 0,
    seed          = 42,
    max.overlaps  = Inf
  ) +
  scale_colour_manual(values = cols) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  coord_fixed(ratio = 1) +
  theme_nature_umap()

cat("Saving plot...\n")

# ── Export ───────────────────────────────────────────────────────────────────
pdf_path <- file.path(fig_dir, "umap_samples.pdf")
png_path <- file.path(fig_dir, "umap_samples.png")

# Nature-contract export: PDF primary, PNG secondary
save_pub_pdf <- function(plot, path, width_mm = 120, height_mm = 100) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  grDevices::cairo_pdf(path, width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
}

save_pub_pdf(p, pdf_path, width_mm = 120, height_mm = 100)
# Use ragg for superior PNG anti-aliasing and clarity
ggsave(png_path, plot = p, width = 120, height = 100, units = "mm", dpi = 600, device = ragg::agg_png)

cat("Saved PDF : ", pdf_path, "\n")
cat("Saved PNG : ", png_path, "\n")
cat("Done.\n")
