####################
# parse_publication_mp13_mp28_trends.R
#
# Description:
#   Publication-quality nature-style trend plots for MP13 and MP28 UCell activity.
#   MP13 is sourced from the legacy T2/T4 filter summary.
#   MP28 is sourced from the active strict mean/median filter summary.
#
# Inputs:
#   parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_sample_ucell_summary_nMP117.csv
#   parse_outs/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_sample_ucell_summary_nMP117.csv
#
# Outputs:
#   parse_outs/publication/mp_trends/figures/parse_mp13_trend.pdf
#   parse_outs/publication/mp_trends/figures/parse_mp13_trend.png
#   parse_outs/publication/mp_trends/figures/parse_mp28_trend.pdf
#   parse_outs/publication/mp_trends/figures/parse_mp28_trend.png
#
# Methodology:
#   analysis/methodology/publication/mp13_mp28_trends_methodology.md
####################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("analysis/common/parse_pipeline_config.R")

# ── Paths ────────────────────────────────────────────────────────────────────
project_dir <- parse_project_root()
paths <- parse_paths(project_dir)

mp13_summary_path <- file.path(paths$parse_outs, "Auto_parse_highres_metaprogram_trends", "Auto_T2T4_gt_T0eR4_filter", "Auto_parse_highres_T2T4_sample_ucell_summary_nMP117.csv")
mp28_summary_path <- file.path(paths$parse_outs, "Auto_parse_highres_metaprogram_trends", "Auto_parse_highres_sample_ucell_summary_nMP117.csv")

out_base <- file.path(paths$parse_outs, "publication", "mp_trends")
fig_dir <- file.path(out_base, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load Data ────────────────────────────────────────────────────────────────
if (!file.exists(mp13_summary_path)) stop("Missing MP13 summary: ", mp13_summary_path)
if (!file.exists(mp28_summary_path)) stop("Missing MP28 summary: ", mp28_summary_path)

summary_mp13 <- read.csv(mp13_summary_path, stringsAsFactors = FALSE) %>% filter(MP == "MP13")
summary_mp28 <- read.csv(mp28_summary_path, stringsAsFactors = FALSE) %>% filter(MP == "MP28")

if (nrow(summary_mp13) == 0) stop("MP13 not found in ", mp13_summary_path)
if (nrow(summary_mp28) == 0) stop("MP28 not found in ", mp28_summary_path)

# ── Plotting Function ────────────────────────────────────────────────────────
create_trend_plot <- function(summary_df, mp_name, title_text) {
  trend_long <- summary_df %>%
    select(MP, sample, mean_score, median_score) %>%
    pivot_longer(cols = c(mean_score, median_score), names_to = "summary_stat", values_to = "score") %>%
    mutate(
      summary_stat = recode(summary_stat, mean_score = "Mean", median_score = "Median"),
      sample = factor(sample, levels = parse_samples)
    )
  
  trend_y_limits <- range(trend_long$score, na.rm = TRUE)
  trend_y_pad <- diff(trend_y_limits) * 0.05
  if (!is.finite(trend_y_pad) || trend_y_pad == 0) trend_y_pad <- 0.01
  trend_y_limits <- trend_y_limits + c(-trend_y_pad, trend_y_pad)
  
  p <- ggplot(trend_long, aes(x = sample, y = score, group = summary_stat, linetype = summary_stat)) +
    geom_line(color = "grey25", linewidth = 0.4) +
    geom_point(aes(fill = sample), shape = 21, size = 2.0, color = "black", stroke = 0.3) +
    scale_fill_manual(values = parse_sample_colours, guide = "none") +
    scale_linetype_manual(values = c("Mean" = "solid", "Median" = "dashed"), name = NULL) +
    coord_cartesian(ylim = trend_y_limits) +
    labs(title = title_text, x = "Timepoint", y = "UCell Score") +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 7),
      axis.text.x = element_text(angle = 0, colour = "black", size = 7, face = "bold"),
      axis.text.y = element_text(colour = "black", size = 7),
      legend.position = "top",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 7),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.margin = margin(0,0,0,0)
    )
  
  p
}

# ── Generate Plots ───────────────────────────────────────────────────────────
p13 <- create_trend_plot(summary_mp13, "MP13", "MP13\n(Stress Response)")
p28 <- create_trend_plot(summary_mp28, "MP28", "MP28\n(Classical Lineage)")

save_pub_plot <- function(plot_obj, mp_name) {
  pdf_path <- file.path(fig_dir, paste0("parse_", mp_name, "_trend.pdf"))
  png_path <- file.path(fig_dir, paste0("parse_", mp_name, "_trend.png"))
  
  # Nature standard width for single panel ~ 60mm. Height ~ 60mm.
  w_mm <- 60
  h_mm <- 60
  
  ggsave(pdf_path, plot = plot_obj, width = w_mm, height = h_mm, units = "mm", device = grDevices::cairo_pdf, fallback_resolution = 600)
  ggsave(png_path, plot = plot_obj, width = w_mm, height = h_mm, units = "mm", dpi = 600, type = "cairo")
  
  message("Saved ", mp_name, " to ", pdf_path)
}

save_pub_plot(p13, "mp13")
save_pub_plot(p28, "mp28")

message("Done.")
