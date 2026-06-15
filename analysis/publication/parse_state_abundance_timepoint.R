####################
# parse_state_abundance_timepoint.R
#
# Description:
#   Publication-quality cell-state abundance visualizations for the
#   6 Parse treatment-response timepoints (T0, T1, T2, T4, R4, eR4).
#
#   1) Stacked bar chart: Proportions among assigned states (excluding
#      Unresolved/Hybrid). Classic state plotted on top, SMG on bottom.
#   2) 2x2 Grid plot: Evolution of the 4 biological states individually,
#      with proportions calculated out of all cells (including Unresolved/Hybrid).
#      Y-axis capped at 50%.
#
# Inputs:
#   parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds
#
# Outputs:
#   parse_outs/publication/state_abundance_timepoint/figures/
#     state_abundance_timepoint.pdf
#     state_abundance_timepoint.png
#   parse_outs/publication/state_abundance_grid/figures/
#     state_abundance_grid.pdf
#     state_abundance_grid.png
#
# Methodology:
#   analysis/methodology/publication/state_abundance_timepoint_methodology.md
#   analysis/methodology/publication/state_abundance_grid_methodology.md
####################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

source("analysis/common/parse_pipeline_config.R")

# ── Paths ────────────────────────────────────────────────────────────────────
project_dir <- parse_project_root()
paths       <- parse_paths(project_dir)

out_base_stack <- file.path(paths$parse_outs, "publication", "state_abundance_timepoint")
dir.create(file.path(out_base_stack, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base_stack, "tables"), recursive = TRUE, showWarnings = FALSE)

out_base_grid <- file.path(paths$parse_outs, "publication", "state_abundance_grid")
dir.create(file.path(out_base_grid, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base_grid, "tables"), recursive = TRUE, showWarnings = FALSE)

# ── Load cached state assignments ────────────────────────────────────────────
assign_path <- file.path(
  paths$parse_outs, "cell_states",
  "Auto_parse_PDOpipeline_topmp_assignments.rds"
)
if (!file.exists(assign_path)) stop("Missing: ", assign_path)
pdo_assign <- readRDS(assign_path)

# ── Base Data & Order ────────────────────────────────────────────────────────
timepoint_order <- c("T0", "T1", "T2", "T4", "R4", "eR4")

state_order <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "Stress-adaptive",
  "SMG-like Metaplasia"
)

# ─────────────────────────────────────────────────────────────────────────────
# 1. STACKED BAR CHART (Resolved cells only)
# ─────────────────────────────────────────────────────────────────────────────
df_resolved <- pdo_assign |>
  filter(
    sample %in% timepoint_order,
    !pdo_state %in% c("Unresolved", "Hybrid")
  )

prop_resolved <- df_resolved |>
  count(sample, pdo_state, name = "n") |>
  right_join(
    expand_grid(sample = factor(timepoint_order, levels = timepoint_order),
                pdo_state = factor(state_order, levels = state_order)),
    by = c("sample", "pdo_state")
  ) |>
  mutate(n = replace_na(n, 0L)) |>
  group_by(sample) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup()

# Set factor levels. state_order puts Classic at top of legend and stack.
prop_resolved <- prop_resolved |>
  mutate(pdo_state = factor(pdo_state, levels = state_order))

state_colours <- c(
  "Classic Proliferative"  = "#E41A1C",
  "Basal to Intest. Meta"  = "#4DAF4A",
  "Stress-adaptive"        = "#984EA3",
  "SMG-like Metaplasia"    = "#FF7F00"
)

theme_nature <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line        = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.35, colour = "black"),
      axis.title       = element_text(size = base_size, face = "bold"),
      axis.text        = element_text(size = base_size - 0.5, colour = "black"),
      axis.text.x      = element_text(face = "bold", size = 9),
      legend.title     = element_text(size = base_size, face = "bold"),
      legend.text      = element_text(size = base_size - 0.5),
      legend.key.size  = unit(3, "mm"),
      legend.key.width = unit(4, "mm"),
      legend.margin    = margin(0, 0, 0, 0),
      legend.box.margin = margin(-2, 0, -2, 0),
      legend.position  = "bottom",
      plot.title       = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = base_size - 0.3, colour = "grey40", hjust = 0),
      plot.margin      = margin(6, 8, 4, 4),
      panel.grid       = element_blank()
    )
}

p_stack <- ggplot(prop_resolved, aes(x = sample, y = pct, fill = pdo_state)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = state_colours,
    breaks = state_order,
    name   = NULL,
    drop   = FALSE
  ) +
  scale_y_continuous(
    name   = "Proportion (%)",
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  labs(
    x        = NULL,
    title    = "Cell-state abundance across treatment timepoints",
    subtitle = "PDO-pipeline states (Approach B / noreg) — Unresolved & Hybrid excluded"
  ) +
  theme_nature() +
  guides(fill = guide_legend(nrow = 1))

# Export Stacked
pdf_path_stack <- file.path(out_base_stack, "figures", "state_abundance_timepoint.pdf")
png_path_stack <- file.path(out_base_stack, "figures", "state_abundance_timepoint.png")

save_pub_pdf <- function(plot, path, width_mm = 100, height_mm = 80) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  grDevices::cairo_pdf(path, width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
}

save_pub_pdf(p_stack, pdf_path_stack, width_mm = 100, height_mm = 80)
ggsave(png_path_stack, plot = p_stack, width = 100, height = 80, units = "mm", dpi = 600)
write.csv(prop_resolved, file.path(out_base_stack, "tables", "state_abundance_timepoint.csv"), row.names = FALSE)


# ─────────────────────────────────────────────────────────────────────────────
# 2. 2x2 GRID PLOT (All cells)
# ─────────────────────────────────────────────────────────────────────────────
df_all <- pdo_assign |> filter(sample %in% timepoint_order)

prop_all <- df_all |>
  count(sample, pdo_state, name = "n") |>
  right_join(
    expand_grid(sample = factor(timepoint_order, levels = timepoint_order),
                pdo_state = unique(df_all$pdo_state)),
    by = c("sample", "pdo_state")
  ) |>
  mutate(n = replace_na(n, 0L)) |>
  group_by(sample) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup() |>
  filter(pdo_state %in% state_order) |>
  mutate(
    sample = factor(sample, levels = timepoint_order),
    pdo_state = factor(pdo_state, levels = state_order)
  )

state_colours_grid <- c(
  "Classic Proliferative"  = "#B2182B",
  "Basal to Intest. Meta"  = "#1B7837",
  "Stress-adaptive"        = "#762A83",
  "SMG-like Metaplasia"    = "#FF7F00"
)

format_title <- function(x) {
  # Match exact capitalization and naming from the legend
  if (x == "Classic Proliferative") return("Classic Proliferative")
  if (x == "Basal to Intest. Meta") return("Basal to Intestinal Metaplasia")
  if (x == "SMG-like Metaplasia") return("SMG-like Metaplasia")
  if (x == "Stress-adaptive") return("Stress-adaptive")
  return(x)
}

create_panel <- function(df, state_name, state_color, show_y_title = FALSE) {
  df_sub <- df |> filter(pdo_state == state_name)
  panel_title <- format_title(state_name)
  
  p <- ggplot(df_sub, aes(x = sample, y = pct)) +
    geom_col(fill = state_color, width = 0.7) +
    scale_y_continuous(
      limits = c(0, 50),
      breaks = seq(0, 50, 10),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = if (show_y_title) "Proportion (%)" else NULL,
      title = panel_title
    ) +
    theme_bw(base_size = 12, base_family = "Arial") +
    theme(
      panel.grid.major = element_line(colour = "grey80", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title.y     = element_text(face = "bold", size = 14),
      axis.text        = element_text(colour = "black", size = 12),
      axis.text.x      = element_text(face = "bold", size = 11),
      axis.text.y      = element_text(size = 11),
      plot.margin      = margin(10, 10, 10, 10)
    )
  return(p)
}

p1 <- create_panel(prop_all, "Classic Proliferative", state_colours_grid["Classic Proliferative"], show_y_title = TRUE)
p2 <- create_panel(prop_all, "Basal to Intest. Meta", state_colours_grid["Basal to Intest. Meta"], show_y_title = FALSE)
p3 <- create_panel(prop_all, "Stress-adaptive", state_colours_grid["Stress-adaptive"], show_y_title = TRUE)
p4 <- create_panel(prop_all, "SMG-like Metaplasia", state_colours_grid["SMG-like Metaplasia"], show_y_title = FALSE)

combined_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "EVOLUTION OF STATE DISTRIBUTIONS OVER EXPERIMENTAL TIME POINTS",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, family = "Arial"),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

pdf_path_grid <- file.path(out_base_grid, "figures", "state_abundance_grid.pdf")
png_path_grid <- file.path(out_base_grid, "figures", "state_abundance_grid.png")

save_pub_pdf_grid <- function(plot, path, width_in = 8, height_in = 6.5) {
  grDevices::cairo_pdf(path, width = width_in, height = height_in, family = "Arial")
  print(plot)
  dev.off()
}

save_pub_pdf_grid(combined_plot, pdf_path_grid, width_in = 8, height_in = 6.5)
ggsave(png_path_grid, plot = combined_plot, width = 8, height = 6.5, dpi = 600)

write.csv(prop_all, file.path(out_base_grid, "tables", "state_abundance_grid.csv"), row.names = FALSE)

####################
# 3. STANDALONE LEGEND PLOT (Single Column)
# Create a dummy plot with all 7 states to extract a consistent, standalone single-column legend
legend_states <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "Immune Infiltrating",
  "Unresolved",
  "Hybrid"
)

legend_colours <- c(
  "Classic Proliferative"          = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia"            = "#FF7F00",
  "Stress-adaptive"                = "#984EA3",
  "Immune Infiltrating"            = "#377EB8",
  "Unresolved"                     = "#CCCCCC",
  "Hybrid"                         = "#000000"
)

dummy_df <- data.frame(
  x = 1:7,
  y = 1:7,
  State = factor(legend_states, levels = legend_states)
)

p_dummy <- ggplot(dummy_df, aes(x = x, y = y, fill = State)) +
  geom_col() +
  scale_fill_manual(
    values = legend_colours,
    name = "State"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10, family = "Arial"),
    legend.text = element_text(size = 9, family = "Arial"),
    legend.key.size = unit(4, "mm")
  ) +
  guides(fill = guide_legend(ncol = 1, title.position = "top"))

# Extract legend function
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) == 0) return(NULL)
  return(tmp$grobs[[leg]])
}

leg <- get_legend(p_dummy)

if (!is.null(leg)) {
  pdf_path_legend <- file.path(out_base_stack, "figures", "state_abundance_legend.pdf")
  png_path_legend <- file.path(out_base_stack, "figures", "state_abundance_legend.png")
  
  # Save PDF
  grDevices::cairo_pdf(pdf_path_legend, width = 3.2, height = 3.5, family = "Arial")
  grid::grid.newpage()
  grid::grid.draw(leg)
  dev.off()
  
  # Save PNG
  png(png_path_legend, width = 3.2, height = 3.5, units = "in", res = 600)
  grid::grid.newpage()
  grid::grid.draw(leg)
  dev.off()
}
####################

cat("Success.\n")

