####################
# Auto_parse_pseudotime_samples.R
#
# Monocle3 pseudotime analysis for the Parse data using samples in place of
# cell-state labels. The trajectory is built on all Parse cells and rooted on
# the T0 sample.
#
# Output:
#   parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime_combined.pdf
#   parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime.rds
#   parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime_cds.rds
#   parse_outs/pseudotime_samples/Auto_parse_sample_pseudotime_metadata.csv
#   parse_outs/summary/Auto_parse_sample_pseudotime_summary.csv
####################

helper_candidates <- c(
  "analysis/trajectory/Auto_parse_pseudotime_helpers.R",
  "analysis/cell_states/Auto_parse_pseudotime_helpers.R",
  "Auto_parse_pseudotime_helpers.R"
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) stop("Could not find Auto_parse_pseudotime_helpers.R")
source(helper_path)

message("=== Parse sample-level pseudotime ===")

result <- build_parse_trajectory()
cds <- result$cds
pt <- result$pseudotime
sample_levels <- result$sample_levels
sample_cols <- get_sample_colours(sample_levels)

sample_vec <- as.character(colData(cds)$sample_label)
names(sample_vec) <- colnames(cds)

sample_counts <- table(factor(sample_vec, levels = sample_levels))
legend_labels <- setNames(
  paste0(sample_levels, " (", as.integer(sample_counts[sample_levels]), ")"),
  sample_levels
)

pdf_path <- file.path(pseudotime_dir, "Auto_parse_sample_pseudotime_combined.pdf")
pdf(pdf_path, width = 14, height = 6, onefile = TRUE)

p_samples <- plot_cells(
  cds,
  color_cells_by = "sample_label",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  cell_size = 0.6
) +
  scale_color_manual(
    values = sample_cols,
    breaks = sample_levels,
    labels = legend_labels[sample_levels],
    name = "Sample",
    na.value = "grey80",
    drop = FALSE,
    guide = guide_legend(override.aes = list(size = 4))
  ) +
  labs(
    title = paste0("Parse samples - T0-rooted trajectory (n = ", length(pt), ")"),
    color = NULL
  ) +
  theme_minimal(base_size = 11)

p_pseudotime <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  cell_size = 0.6
) +
  scale_color_viridis_c(na.value = "grey85") +
  labs(
    title = paste0("Pseudotime - root sample: ", root_sample),
    color = "Pseudotime"
  ) +
  theme_minimal(base_size = 11)

print(p_samples + p_pseudotime + plot_layout(guides = "collect"))
dev.off()

meta_out <- data.frame(
  cell = names(pt),
  sample = sample_vec[names(pt)],
  pseudotime = as.numeric(pt),
  stringsAsFactors = FALSE
)

summary_out <- meta_out %>%
  mutate(sample = factor(sample, levels = sample_levels)) %>%
  group_by(sample) %>%
  summarise(
    n_cells = n(),
    n_cells_with_pseudotime = sum(is.finite(pseudotime)),
    median_pseudotime = safe_median(pseudotime),
    mean_pseudotime = safe_mean(pseudotime),
    .groups = "drop"
  ) %>%
  arrange(sample)

saveRDS(pt, file.path(pseudotime_dir, "Auto_parse_sample_pseudotime.rds"))
saveRDS(cds, file.path(pseudotime_dir, "Auto_parse_sample_pseudotime_cds.rds"))
write_csv(meta_out, file.path(pseudotime_dir, "Auto_parse_sample_pseudotime_metadata.csv"))
write_csv(summary_out, file.path(summary_dir, "Auto_parse_sample_pseudotime_summary.csv"))

message("Saved Parse sample pseudotime outputs to: ", pseudotime_dir)
