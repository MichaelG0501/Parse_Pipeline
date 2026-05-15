####################
# parse_metaprogram_select_optimal_nMP.R
#
# Description:
#   Selects the default Parse nMP by combining average silhouette width and WSS
#   elbow diagnostics from the GeneNMF metaprogram sweep.
#
# Inputs:
#   parse_outs/Auto_parse_metaprograms/Metaprogrammes_Results/Auto_parse_geneNMF_metaprograms_nMP_<k>.rds
#
# Outputs:
#   parse_outs/Auto_parse_metaprograms/Auto_parse_optimal_nMP.txt
#   parse_outs/Auto_parse_metaprograms/Auto_parse_optimal_nMP_metrics.csv
#   parse_outs/Auto_parse_metaprograms/Auto_parse_optimal_nMP_metrics.png
#   parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds
#   parse_outs/logs/run_summaries/parse_metaprogram_select_optimal_nMP_*.txt
#
# Cache / replot:
#   Reads existing GeneNMF sweep outputs only; rerun is quick and can regenerate
#   the metrics plot after styling changes without repeating GeneNMF.
#
# Methodology:
#   analysis/methodology/metaprograms/metaprogram_select_optimal_nMP_methodology.md
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_metaprogram_select_optimal_nMP",
  parameters = list(k_range = "4:35"),
  input_files = "parse_outs/Auto_parse_metaprograms/Metaprogrammes_Results/Auto_parse_geneNMF_metaprograms_nMP_<k>.rds",
  output_files = c(
    "parse_outs/Auto_parse_metaprograms/Auto_parse_optimal_nMP.txt",
    "parse_outs/Auto_parse_metaprograms/Auto_parse_optimal_nMP_metrics.csv",
    "parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds"
  )
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

suppressPackageStartupMessages({
  library(cluster)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(GeneNMF)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath(file.path(parse_project_root(), "analysis/metaprograms/parse_metaprogram_select_optimal_nMP.R"))
}

script_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(script_dir, "..", ".."))
out_dir <- file.path(project_dir, "parse_outs")
setwd(out_dir)

parse_mp_dir <- file.path("Auto_parse_metaprograms")
results_dir <- file.path(parse_mp_dir, "Metaprogrammes_Results")
if (!dir.exists(results_dir)) {
  stop("Missing metaprogram results dir: ", results_dir)
}

k_vals <- 4:35
avg_sil_widths <- numeric(length(k_vals))
wss_vals <- numeric(length(k_vals))

for (i in seq_along(k_vals)) {
  k <- k_vals[i]
  rds_path <- file.path(results_dir, paste0("Auto_parse_geneNMF_metaprograms_nMP_", k, ".rds"))

  if (file.exists(rds_path)) {
    mp_res <- readRDS(rds_path)

    dist_mat <- as.dist(1 - mp_res$programs.similarity)
    cluster_assignments <- cutree(mp_res$programs.tree, k = k)

    sil <- silhouette(cluster_assignments, dist = dist_mat)
    avg_sil_widths[i] <- summary(sil)$avg.width

    wss_k <- 0
    dist_m <- as.matrix(dist_mat)
    for (clust_id in unique(cluster_assignments)) {
      idx <- which(cluster_assignments == clust_id)
      if (length(idx) > 1) {
        cluster_dist <- dist_m[idx, idx]
        wss_k <- wss_k + sum(cluster_dist^2) / (2 * length(idx))
      }
    }
    wss_vals[i] <- wss_k
  } else {
    avg_sil_widths[i] <- NA
    wss_vals[i] <- NA
  }
}

df_metrics <- data.frame(nMP = k_vals, Silhouette = avg_sil_widths, WSS = wss_vals)
print(df_metrics)

find_knee <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3) {
    stop("Need at least 3 finite points to detect knee.")
  }
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  x1 <- x_norm[1]
  y1 <- y_norm[1]
  x2 <- x_norm[length(x_norm)]
  y2 <- y_norm[length(y_norm)]
  dists <- abs((y2 - y1) * x_norm - (x2 - x1) * y_norm + x2 * y1 - y2 * x1) /
    sqrt((y2 - y1)^2 + (x2 - x1)^2)
  x[which.max(dists)]
}

sil_knee <- find_knee(df_metrics$nMP, df_metrics$Silhouette)
wss_knee <- find_knee(df_metrics$nMP, df_metrics$WSS)
optimal_nMP <- sil_knee

message(paste0("Silhouette inflection point: nMP = ", sil_knee))
message(paste0("WSS elbow point: nMP = ", wss_knee))
message(paste0("Selected optimal nMP: ", optimal_nMP))

p1 <- ggplot(df_metrics, aes(x = nMP, y = Silhouette)) +
  geom_line(color = "steelblue", linewidth = 1, na.rm = TRUE) +
  geom_point(color = "steelblue", size = 3, na.rm = TRUE) +
  geom_vline(xintercept = sil_knee, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate(
    "text",
    x = sil_knee + 0.5,
    y = max(df_metrics$Silhouette, na.rm = TRUE),
    label = paste0("Inflection: ", sil_knee),
    hjust = 0,
    color = "red",
    size = 3.5
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(k_vals), max(k_vals), by = 2)) +
  labs(title = "Silhouette Analysis", x = "Number of MetaPrograms (nMP)", y = "Average Silhouette Width") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- ggplot(df_metrics, aes(x = nMP, y = WSS)) +
  geom_line(color = "darkred", linewidth = 1, na.rm = TRUE) +
  geom_point(color = "darkred", size = 3, na.rm = TRUE) +
  geom_vline(xintercept = wss_knee, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate(
    "text",
    x = wss_knee + 0.5,
    y = max(df_metrics$WSS, na.rm = TRUE) * 0.95,
    label = paste0("Elbow: ", wss_knee),
    hjust = 0,
    color = "red",
    size = 3.5
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(k_vals), max(k_vals), by = 2)) +
  labs(title = "Elbow Method (WSS)", x = "Number of MetaPrograms (nMP)", y = "Total Within Sum of Squares") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

combined <- p1 + p2
ggsave(file.path(parse_mp_dir, "Auto_parse_optimal_nMP_metrics.png"), combined, width = 16, height = 6, dpi = 300)
write.csv(df_metrics, file.path(parse_mp_dir, "Auto_parse_optimal_nMP_metrics.csv"), row.names = FALSE)
writeLines(as.character(optimal_nMP), file.path(parse_mp_dir, "Auto_parse_optimal_nMP.txt"))

selected_rds <- file.path(results_dir, paste0("Auto_parse_geneNMF_metaprograms_nMP_", optimal_nMP, ".rds"))
if (!file.exists(selected_rds)) {
  stop("Selected nMP file missing: ", selected_rds)
}

geneNMF.metaprograms <- readRDS(selected_rds)
saveRDS(geneNMF.metaprograms, file.path(parse_mp_dir, "Auto_parse_MP_outs_default.rds"), compress = FALSE)

anno_colors <- brewer.pal(n = min(optimal_nMP, 12), name = "Paired")
anno_colors <- anno_colors[seq_len(length(geneNMF.metaprograms$metaprograms.genes))]
names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

png(file.path(parse_mp_dir, "Auto_parse_metaprograms_heatmap.png"), width = 3000, height = 2500, res = 300)
plotMetaPrograms(
  geneNMF.metaprograms,
  annotation_colors = anno_colors,
  similarity.cutoff = c(0, 1)
)
dev.off()

script_run_status <- "success"
message("Saved optimal nMP outputs to ", parse_mp_dir)
