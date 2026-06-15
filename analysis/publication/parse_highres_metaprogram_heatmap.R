####################
# parse_highres_metaprogram_heatmap.R
#
# Description:
#   Publication-quality heatmap of the high-resolution metaprograms similarity
#   matrix. It uses the "many nMP" splitting strategy (e.g., nMP=117) from the
#   legacy comparison filter. To prevent the massive number of splits from
#   turning the heatmap into a solid white block via pheatmap gaps, it uses
#   ComplexHeatmap with precisely controlled 0.2mm row/col gaps and thin borders.
#
# Inputs:
#   parse_outs/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds
#   parse_outs/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_geneNMF_metaprograms_nMP<k>.rds
#
# Outputs:
#   parse_outs/publication/highres_metaprogram_heatmap/figures/
#     highres_metaprogram_heatmap.pdf
#     highres_metaprogram_heatmap.png
#   parse_outs/logs/run_summaries/parse_highres_metaprogram_heatmap_*.txt
#
# Methodology:
#   analysis/methodology/publication/highres_metaprogram_heatmap_methodology.md
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_highres_metaprogram_heatmap",
  parameters = list(
    splitting_strategy = "total_nmf_programs / 2",
    visualization = "ComplexHeatmap"
  ),
  input_files = c(
    "parse_outs/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds",
    "parse_outs/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_geneNMF_metaprograms_nMP<k>.rds"
  ),
  output_files = c(
    "parse_outs/publication/highres_metaprogram_heatmap/figures/highres_metaprogram_heatmap.pdf",
    "parse_outs/publication/highres_metaprogram_heatmap/figures/highres_metaprogram_heatmap.png"
  )
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
  library(RColorBrewer)
})

# ── Paths ────────────────────────────────────────────────────────────────────
project_dir <- parse_project_root()
paths       <- parse_paths(project_dir)
out_base    <- file.path(paths$parse_outs, "publication", "highres_metaprogram_heatmap")
fig_dir     <- file.path(out_base, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Determine nMP ────────────────────────────────────────────────────────────
geneNMF_out_path <- file.path(paths$parse_outs, "Auto_parse_metaprograms", "Auto_parse_geneNMF_outs.rds")
if (!file.exists(geneNMF_out_path)) stop("Missing GeneNMF programmes: ", geneNMF_out_path)

message("Loading GeneNMF programmes to determine nMP...")
geneNMF.programs <- readRDS(geneNMF_out_path)
total_nmf_programs <- sum(vapply(geneNMF.programs, function(x) ncol(x$w), numeric(1)))
nMP <- as.integer(total_nmf_programs / 2)

if (total_nmf_programs / 2 != nMP) {
  stop("Total NMF programme count is odd: ", total_nmf_programs)
}

message("Calculated target nMP: ", nMP)

# ── Load High-Res Metaprograms ───────────────────────────────────────────────
highres_dir <- file.path(paths$parse_outs, "Auto_parse_highres_metaprogram_trends")
mp_path <- file.path(highres_dir, paste0("Auto_parse_highres_geneNMF_metaprograms_nMP", nMP, ".rds"))

if (!file.exists(mp_path)) {
  stop(
    "Missing high-resolution metaprograms file: ", mp_path, "\n",
    "Please run legacy_parse_highres_mp_t2t4_comparison_filter.R (or similar) ",
    "to generate the geneNMF.metaprograms for nMP=", nMP, " before plotting."
  )
}

message("Loading cached high-resolution metaprograms: ", mp_path)
mp.res <- readRDS(mp_path)

# ── Extract data for heatmap ─────────────────────────────────────────────────
J <- mp.res[["programs.similarity"]]
tree <- mp.res[["programs.tree"]]
cl_members <- mp.res[["programs.clusters"]]

# Limit similarity to [0, 1] as in GeneNMF
J[J < 0] <- 0
J[J > 1] <- 1

# Ensure cl_members follows the dendrogram order
labs.order <- labels(as.dendrogram(tree))
cl_names <- names(cl_members)

# Prepend MP prefix
cl_members_labels <- paste0("MP", cl_members)
names(cl_members_labels) <- cl_names

# The cluster assignment for each row/col
row_split_factor <- factor(cl_members_labels, levels = unique(cl_members_labels[labs.order]))

# ── Define colours ───────────────────────────────────────────────────────────
col_fun <- colorRamp2(seq(0, 1, length = 100), viridis(100, option = "A", direction = -1))

# Recreate the annotation colours used in GeneNMF
n_colors <- min(nMP, 12)
anno_colors <- brewer.pal(n = max(n_colors, 3), name = "Paired")
anno_colors_full <- rep(anno_colors, length.out = nMP)
names(anno_colors_full) <- levels(row_split_factor)

# Simple annotation tracks to show the MP groupings
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

# ── Plot Heatmap ─────────────────────────────────────────────────────────────
message("Generating ComplexHeatmap...")
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
  column_title = paste0("High-Resolution Metaprograms Similarity (nMP = ", nMP, ")"),
  column_title_gp = gpar(fontsize = 16, fontface = "bold")
)

# ── Export ───────────────────────────────────────────────────────────────────
pdf_path <- file.path(fig_dir, "highres_metaprogram_heatmap.pdf")
png_path <- file.path(fig_dir, "highres_metaprogram_heatmap.png")

message("Exporting PDF...")
pdf(pdf_path, width = 12, height = 11, useDingbats = FALSE)
draw(ht)
dev.off()

message("Exporting PNG...")
png(png_path, width = 3600, height = 3300, res = 300)
draw(ht)
dev.off()

script_run_status <- "success"
message("Saved PDF: ", pdf_path)
message("Saved PNG: ", png_path)
message("Done.")
