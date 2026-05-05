suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("ggplot2")
  library("scales")
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline/analysis/summary/Auto_parse_summary.R")
}

project_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))
out_dir <- file.path(project_dir, "parse_outs")
setwd(out_dir)

summary_df <- read.csv(file.path("summary", "Auto_parse_filtering_summary.csv"), check.names = FALSE)
merged_obj <- readRDS("Auto_parse_merged.rds")

after_df <- merged_obj@meta.data %>%
  count(orig.ident, name = "Final") %>%
  rename(sample = orig.ident)

summary_df <- summary_df %>%
  left_join(after_df, by = "sample")

thresholds <- read.csv(file.path("summary", "Auto_parse_thresholds.csv"))

global_df <- summary_df %>%
  summarise(
    Raw = sum(raw),
    Doublets = sum(doublet_singlets),
    Mito = sum(`mito_DNA\npercentage < 15`),
    Genes = sum(`number of\ngenes`),
    HK = sum(`housekeeping\nexpression > 3`),
    Final = sum(Final)
  ) %>%
  pivot_longer(everything(), names_to = "step", values_to = "cell_count") %>%
  mutate(
    step = factor(
      step,
      levels = c("Raw", "Doublets", "Mito", "Genes", "HK", "Final"),
      labels = c(
        "Raw",
        "Doublet singlets",
        paste0("Mito_DNA < ", thresholds$max_mt),
        paste0("Genes ", thresholds$min_ngenes, "-", thresholds$max_ngenes),
        paste0("HK_expr > ", thresholds$min_hk_expr),
        "Final"
      )
    )
  )

global_bar <- ggplot(global_df, aes(x = step, y = cell_count, fill = step)) +
  geom_col(width = 0.7, color = "gray20", linewidth = 0.25) +
  geom_text(aes(label = comma(cell_count)), vjust = -0.4, size = 3.2) +
  scale_fill_brewer(palette = "Blues", direction = -1, guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Total Cells Remaining After Each Filter", x = NULL, y = "Number of Cells") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path("plots/Filtering_and_Summary", "global_filter_bar.png"), global_bar, width = 8, height = 4, dpi = 300)

per_sample_df <- summary_df %>%
  transmute(sample, before = raw, after = Final) %>%
  pivot_longer(c(before, after), names_to = "filter_status", values_to = "cell_count") %>%
  mutate(filter_status = factor(filter_status, levels = c("before", "after")))

per_sample_bar <- ggplot(per_sample_df, aes(x = sample, y = cell_count + 1, fill = sample, alpha = filter_status)) +
  geom_col(width = 0.75, position = position_dodge(width = 0.85), color = "gray20", linewidth = 0.15) +
  scale_y_log10(
    breaks = c(1, 11, 101, 1001, 10001, 100001),
    labels = c("0", "10", "100", "1k", "10k", "100k")
  ) +
  scale_alpha_manual(values = c(before = 1, after = 0.6), labels = c(before = "Raw", after = "Final"), name = "") +
  labs(title = "Cells Before vs After Filtering", x = NULL, y = "Number of Cells (log10 scale)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

ggsave(file.path("plots/Filtering_and_Summary", "per_sample_before_after.png"), per_sample_bar, width = 10, height = 5, dpi = 300)

step_df <- summary_df %>%
  select(sample, raw, doublet_singlets, `mito_DNA\npercentage < 15`, `number of\ngenes`, `housekeeping\nexpression > 3`, Final) %>%
  pivot_longer(-sample, names_to = "step", values_to = "cell_count") %>%
  mutate(
    step = factor(
      step,
      levels = c("raw", "doublet_singlets", "mito_DNA\npercentage < 15", "number of\ngenes", "housekeeping\nexpression > 3", "Final"),
      labels = c(
        "Raw",
        "Doublet singlets",
        paste0("Mito_DNA < ", thresholds$max_mt),
        paste0("Genes ", thresholds$min_ngenes, "-", thresholds$max_ngenes),
        paste0("HK_expr > ", thresholds$min_hk_expr),
        "Final"
      )
    )
  )

step_plot <- ggplot(step_df, aes(x = step, y = cell_count, group = sample, color = sample)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_y_continuous(labels = comma) +
  labs(title = "Per-sample filtering trajectory", x = NULL, y = "Number of Cells") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path("plots/Filtering_and_Summary", "per_sample_filtering_steps.png"), step_plot, width = 10, height = 5, dpi = 300)
