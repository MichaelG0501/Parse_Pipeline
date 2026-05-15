####################
# parse_mp_state_abundance_activity_approachB_noreg.R
#
# Description:
#   Preferred Parse state workflow. Scores scATLAS, PDO-pipeline, and
#   Parse-derived metaprograms across the 9 Parse/PDO/SUR1090 samples, then
#   assigns the PDO-pipeline state labels using Approach B with noreg
#   normalized scores. This is the active state-definition script for
#   downstream analysis.
#
# Inputs:
#   parse_outs/by_samples/<sample>/Auto_<sample>_final.rds
#   parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds
#   parse_outs/Auto_parse_metaprograms/Auto_parse_UCell_scores_filtered*.rds
#   /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   /rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds
#
# Outputs:
#   parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.{rds,csv}
#   parse_outs/cell_states/Auto_parse_PDOpipeline_mp_adj_noreg.rds
#   parse_outs/cell_states/Auto_parse_topmp_abundance_*_9samples.pdf
#   parse_outs/cell_states/Auto_parse_mp_activity_boxplots_*.pdf
#   parse_outs/logs/run_summaries/parse_mp_state_abundance_activity_approachB_noreg_*.txt
#
# Cache / replot:
#   Existing UCell score matrices and assignment RDS files are reused when all
#   sample levels are present. Plot-only changes can therefore be rerun quickly
#   from cached RDS objects without recomputing UCell scores.
#
# Methodology:
#   analysis/methodology/cell_states/mp_state_abundance_activity_approachB_noreg_methodology.md
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_helpers.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_mp_state_abundance_activity_approachB_noreg",
  parameters = list(
    state_method = parse_state_definition$preferred_method,
    normalization = parse_state_definition$normalization,
    min_group_score = parse_state_definition$min_group_score,
    hybrid_gap = parse_state_definition$hybrid_gap
  ),
  input_files = c(
    "parse_outs/by_samples/<sample>/Auto_<sample>_final.rds",
    "parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds",
    parse_reference_paths$scatlas_metaprograms,
    parse_reference_paths$pdo_metaprograms
  ),
  output_files = c(
    "parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds",
    "parse_outs/cell_states/Auto_parse_PDOpipeline_mp_adj_noreg.rds",
    "parse_outs/cell_states/Auto_parse_topmp_abundance_PDOpipeline_9samples.pdf"
  )
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

suppressPackageStartupMessages({

  library(Seurat)
  library(Matrix)
  library(UCell)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

root_dir <- parse_project_root()
qc_dir <- file.path(root_dir, "parse_outs")
out_dir <- file.path(qc_dir, "cell_states")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_order_all <- parse_all_samples
sample_order_excluding_pdo_sur1090 <- setdiff(
  sample_order_all,
  c("PDO", "SUR1090_Treated", "SUR1090_Untreated")
)
sample_cols_all <- setNames(hcl.colors(length(sample_order_all), palette = "Dark 3"), sample_order_all)

get_counts <- function(obj) {
  suppressWarnings({
    tryCatch(
      GetAssayData(obj, assay = "RNA", layer = "counts"),
      error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
    )
  })
}

label_from_map <- function(mp_vec, label_map) {
  out <- unname(label_map[mp_vec])
  names(out) <- names(mp_vec)
  out
}

label_mp <- function(mp_vec, descriptions) {
  desc <- descriptions[mp_vec]
  desc[is.na(desc)] <- mp_vec[is.na(desc)]
  out <- paste0(mp_vec, "_", desc)
  names(out) <- names(mp_vec)
  out
}

z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]
  study_sd <- clust_df %>%
    group_by(.study) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1
  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()
  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

ordered_block <- function(mps, tree_order_names) {
  out <- intersect(tree_order_names, mps)
  c(out, setdiff(mps, out))
}

make_prop_data <- function(label_vec, cell_meta, label_order, sample_order) {
  data.frame(
    sample = cell_meta[names(label_vec), "sample"],
    label = as.character(label_vec),
    stringsAsFactors = FALSE
  ) %>%
    filter(sample %in% sample_order, label %in% label_order) %>%
    count(sample, label, name = "n") %>%
    right_join(tidyr::expand_grid(sample = sample_order, label = label_order), by = c("sample", "label")) %>%
    mutate(n = replace_na(n, 0L)) %>%
    group_by(sample) %>%
    mutate(pct = 100 * n / pmax(sum(n), 1)) %>%
    ungroup()
}

plot_abundance <- function(label_vec, cell_meta, label_order, col_map, sample_order, title_text) {
  totals_df <- cell_meta %>%
    filter(sample %in% sample_order) %>%
    count(sample, name = "total_cells") %>%
    mutate(sample = factor(sample, levels = sample_order))
  scale_factor <- max(totals_df$total_cells, na.rm = TRUE) / 100
  if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1

  prop_df <- make_prop_data(label_vec, cell_meta, label_order, sample_order) %>%
    mutate(
      sample = factor(sample, levels = sample_order),
      label = factor(label, levels = rev(label_order))
    )

  ggplot(prop_df, aes(x = sample, y = pct, fill = label)) +
    geom_col(width = 0.78, colour = NA) +
    geom_point(
      data = totals_df,
      aes(x = sample, y = total_cells / scale_factor),
      inherit.aes = FALSE,
      colour = "black",
      size = 2.4
    ) +
    geom_line(
      data = totals_df,
      aes(x = sample, y = total_cells / scale_factor, group = 1),
      inherit.aes = FALSE,
      colour = "black",
      alpha = 0.45,
      linetype = "dashed",
      linewidth = 0.45
    ) +
    scale_fill_manual(values = col_map, breaks = label_order, drop = FALSE, name = NULL) +
    scale_y_continuous(
      name = "Proportion (%)",
      limits = c(0, 100),
      expand = c(0, 0),
      sec.axis = sec_axis(~ . * scale_factor, name = "Total Cells", labels = comma)
    ) +
    labs(x = NULL, title = title_text) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 17, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 15, colour = "black"),
      axis.title.y = element_text(size = 17, face = "bold"),
      axis.title.y.right = element_text(size = 17, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 13),
      legend.key.size = unit(0.55, "cm"),
      plot.margin = margin(12, 18, 12, 12)
    ) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))
}

make_activity_plot <- function(ucell_scores, cell_meta, mp_order, mp_descriptions, sample_order, title_text) {
  sample_cols <- sample_cols_all[sample_order]
  activity_long <- as.data.frame(ucell_scores[, mp_order, drop = FALSE]) %>%
    mutate(cell = rownames(.), sample = cell_meta[rownames(.), "sample"]) %>%
    filter(sample %in% sample_order) %>%
    pivot_longer(cols = all_of(mp_order), names_to = "MP", values_to = "score") %>%
    mutate(
      sample = factor(sample, levels = sample_order),
      MP = factor(MP, levels = mp_order)
    )

  activity_stats <- activity_long %>%
    group_by(MP) %>%
    summarise(
      p_value = tryCatch(kruskal.test(score ~ sample)$p.value, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significance = case_when(
        is.na(p_adj) ~ "",
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  annot_df <- activity_long %>%
    group_by(MP) %>%
    summarise(y_pos = max(score, na.rm = TRUE), .groups = "drop") %>%
    left_join(activity_stats, by = "MP") %>%
    mutate(
      y_pos = y_pos + 0.025,
      label = ifelse(!is.na(p_adj) & p_adj < 0.05, significance, "")
    )

  mp_axis_labels <- setNames(
    ifelse(!is.na(mp_descriptions[mp_order]), paste0(mp_order, "\n", mp_descriptions[mp_order]), mp_order),
    mp_order
  )

  p <- ggplot(activity_long, aes(x = MP, y = score, fill = sample, color = sample)) +
    geom_boxplot(
      position = position_dodge(width = 0.82),
      width = 0.62,
      outlier.shape = NA,
      alpha = 0.78,
      linewidth = 0.28,
      color = "black"
    ) +
    geom_text(
      data = annot_df %>% filter(label != ""),
      aes(x = MP, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    scale_fill_manual(values = sample_cols, drop = FALSE, name = "Sample") +
    scale_color_manual(values = sample_cols, guide = "none", drop = FALSE) +
    scale_x_discrete(labels = mp_axis_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
    labs(
      title = title_text,
      subtitle = "Per-cell UCell scores; stars mark BH-adjusted Kruskal-Wallis p < 0.05 across samples.",
      x = NULL,
      y = "UCell score"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 11, colour = "grey35"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, colour = "black"),
      axis.text.y = element_text(size = 11, colour = "black"),
      axis.line.x = element_blank(),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11),
      plot.margin = margin(12, 18, 12, 12)
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))

  list(plot = p, stats = activity_stats)
}

message("Loading sample count matrices")
sample_files <- file.path(
  qc_dir,
  "by_samples",
  sample_order_all,
  paste0("Auto_", sample_order_all, "_final.rds")
)
names(sample_files) <- sample_order_all
missing_files <- sample_files[!file.exists(sample_files)]
if (length(missing_files) > 0) {
  stop("Missing sample RDS files: ", paste(missing_files, collapse = ", "))
}

sample_genes <- lapply(sample_files, function(path) rownames(readRDS(path)))
common_genes <- Reduce(intersect, sample_genes)
cell_meta <- list()
need_scatlas_scoring <- !all(file.exists(
  file.path(out_dir, "Auto_parse_scATLAS_UCell_scores.rds"),
  file.path(out_dir, "Auto_parse_scATLAS_mp_adj_noreg.rds"),
  file.path(out_dir, "Auto_parse_topmp_assignments.rds")
))
need_pdo_scoring <- !all(file.exists(
  file.path(out_dir, "Auto_parse_PDOpipeline_UCell_scores.rds"),
  file.path(out_dir, "Auto_parse_PDOpipeline_mp_adj_noreg.rds"),
  file.path(out_dir, "Auto_parse_PDOpipeline_mp_adj_all_noreg.rds"),
  file.path(out_dir, "Auto_parse_PDOpipeline_topmp_assignments.rds")
))

# Parse MP setup and scoring check
parse_mp_dir <- file.path(qc_dir, "Auto_parse_metaprograms")
parse_geneNMF <- readRDS(file.path(parse_mp_dir, "Auto_parse_MP_outs_default.rds"))
parse_mp_genes <- parse_geneNMF$metaprograms.genes
bad_mps_parse <- which(parse_geneNMF$metaprograms.metrics$silhouette < 0)
parse_retained <- setdiff(names(parse_mp_genes), paste0("MP", bad_mps_parse))
coverage_tbl_parse <- parse_geneNMF$metaprograms.metrics$sampleCoverage
names(coverage_tbl_parse) <- paste0("MP", seq_along(coverage_tbl_parse))
parse_retained <- setdiff(parse_retained, names(coverage_tbl_parse)[coverage_tbl_parse < 0.25])

optimal_txt <- file.path(parse_mp_dir, "Auto_parse_optimal_nMP.txt")
optimal_nMP <- if (file.exists(optimal_txt)) readLines(optimal_txt, warn = FALSE)[1] else length(parse_mp_genes)
parse_ucell_path <- file.path(parse_mp_dir, paste0("Auto_parse_UCell_scores_filtered_nMP", optimal_nMP, ".rds"))
if (!file.exists(parse_ucell_path)) {
  parse_ucell_path <- file.path(parse_mp_dir, "Auto_parse_UCell_scores_filtered.rds")
}

need_parse_scoring <- TRUE
if (file.exists(parse_ucell_path)) {
  tmp_uc <- readRDS(parse_ucell_path)
  tmp_samples <- sub("_.*", "", rownames(tmp_uc))
  if (all(sample_order_all %in% tmp_samples)) {
    need_parse_scoring <- FALSE
  }
  rm(tmp_uc)
}

need_counts <- need_scatlas_scoring || need_pdo_scoring || need_parse_scoring
counts_list <- list()
for (sample in sample_order_all) {
  message("Loading metadata", if (need_counts) " and counts" else "", " for ", sample)
  obj <- readRDS(sample_files[[sample]])
  old_cells <- colnames(obj)
  new_cells <- paste(sample, old_cells, sep = "__")
  if (need_counts) {
    counts <- get_counts(obj)[common_genes, , drop = FALSE]
    colnames(counts) <- new_cells
    counts_list[[sample]] <- counts
  }
  cell_meta[[sample]] <- data.frame(
    cell = new_cells,
    original_cell = old_cells,
    sample = sample,
    study = "Parse_SUR1090",
    stringsAsFactors = FALSE
  )
  rm(obj)
  if (exists("counts")) rm(counts)
  gc()
}
cell_meta <- bind_rows(cell_meta)
rownames(cell_meta) <- cell_meta$cell
sample_var <- setNames(cell_meta$sample, cell_meta$cell)
study_var <- setNames(cell_meta$study, cell_meta$cell)
if (need_counts) {
  counts_all <- do.call(cbind, counts_list)
  rm(counts_list)
  gc()
}

####################
# scATLAS outputs and abundance PDF
####################
scatlas_geneNMF <- readRDS(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds"
)
scatlas_desc <- c(
  "MP1" = "G2M Cell Cycle",
  "MP9" = "G1S Cell Cycle",
  "MP2" = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5" = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8" = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7" = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
)
scatlas_state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "Immune Infiltrating" = c("MP15")
)
scatlas_state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Immune Infiltrating" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)
scatlas_mp_cols <- c(
  "MP1_G2M Cell Cycle" = "#B0B0B0",
  "MP2_MYC-related Proliferation" = "#E41A1C",
  "MP5_Epithelial IFN Resp." = "#66C2A5",
  "MP7_DNA Damage Repair" = "#999999",
  "MP8_Intestinal Diff." = "#FC8D62",
  "MP9_G1S Cell Cycle" = "#C0C0C0",
  "MP10_Columnar Diff." = "#A6D854",
  "MP12_Neuro-responsive Epi" = "#E78AC3",
  "MP13_Hypoxic Inflam. Epi." = "#984EA3",
  "MP14_Hypoxia Adapted Epi." = "#8DA0CB",
  "MP15_Immune Infiltration" = "#377EB8",
  "MP16_Secretory Diff. (Gastric)" = "#FFD92F",
  "MP17_Basal-like Transition" = "#4DAF4A",
  "MP18_Secretory Diff. (Intest.)" = "#FF7F00"
)
scatlas_retained <- rownames(scatlas_geneNMF$metaprograms.metrics)[
  scatlas_geneNMF$metaprograms.metrics$silhouette >= 0
]
scatlas_tree_order <- scatlas_geneNMF$programs.clusters[scatlas_geneNMF$programs.tree$order]
scatlas_tree_order <- paste0("MP", unique(scatlas_tree_order))
scatlas_tree_order <- scatlas_tree_order[scatlas_tree_order %in% scatlas_retained]
scatlas_cc <- intersect(c("MP1", "MP7", "MP9"), scatlas_retained)
scatlas_mp_order <- ordered_block(scatlas_cc, scatlas_tree_order)
for (state in names(scatlas_state_groups)) {
  scatlas_mp_order <- c(
    scatlas_mp_order,
    ordered_block(intersect(scatlas_state_groups[[state]], scatlas_retained), scatlas_tree_order)
  )
}
scatlas_mp_order <- unique(c(scatlas_mp_order, setdiff(scatlas_retained, scatlas_mp_order)))
scatlas_noncc_order <- setdiff(scatlas_mp_order, scatlas_cc)
scatlas_label_order <- unname(label_mp(scatlas_mp_order, scatlas_desc))
scatlas_noncc_label_order <- unname(label_mp(scatlas_noncc_order, scatlas_desc))

if (need_scatlas_scoring) {
  message("Scoring scATLAS UCell signatures")
  scatlas_features <- scatlas_geneNMF$metaprograms.genes
  scatlas_features <- lapply(scatlas_features, intersect, rownames(counts_all))
  scatlas_ucell <- UCell::ScoreSignatures_UCell(
    matrix = counts_all,
    features = scatlas_features,
    maxRank = 1500, chunk.size = 1000, ncores = 1, force.gc = TRUE
  )
  scatlas_ucell <- as.data.frame(scatlas_ucell)
  colnames(scatlas_ucell) <- sub("_UCell$", "", colnames(scatlas_ucell))
  scatlas_ucell <- as.matrix(scatlas_ucell[, intersect(scatlas_retained, colnames(scatlas_ucell)), drop = FALSE])
  saveRDS(scatlas_ucell, file.path(out_dir, "Auto_parse_scATLAS_UCell_scores.rds"))

  scatlas_mp_adj <- z_normalise(scatlas_ucell, sample_var, study_var)
  saveRDS(scatlas_mp_adj, file.path(out_dir, "Auto_parse_scATLAS_mp_adj_noreg.rds"))

  # Assignment logic Approach B
  scatlas_noncc_avail <- intersect(scatlas_noncc_order, colnames(scatlas_mp_adj))
  scatlas_group_max <- sapply(scatlas_state_groups, function(mps) {
    mps_avail <- intersect(mps, scatlas_noncc_avail)
    if (length(mps_avail) == 0) return(rep(0, nrow(scatlas_mp_adj)))
    if (length(mps_avail) == 1) return(as.numeric(scatlas_mp_adj[, mps_avail]))
    apply(scatlas_mp_adj[, mps_avail, drop = FALSE], 1, max)
  })
  scatlas_group_max <- as.matrix(scatlas_group_max)
  rownames(scatlas_group_max) <- rownames(scatlas_mp_adj)

  THRESHOLD <- parse_state_definition$min_group_score
  HYBRID_GAP_B <- parse_state_definition$hybrid_gap
  scatlas_best_idx <- max.col(scatlas_group_max, ties.method = "first")
  scatlas_best_val <- apply(scatlas_group_max, 1, max)
  scatlas_top_state <- names(scatlas_state_groups)[scatlas_best_idx]
  scatlas_top_state[scatlas_best_val < THRESHOLD] <- "Unresolved"
  
  scatlas_sorted_groups <- t(apply(scatlas_group_max, 1, sort, decreasing = TRUE))
  scatlas_gap <- scatlas_sorted_groups[, 1] - scatlas_sorted_groups[, 2]
  scatlas_top_state[(scatlas_gap < HYBRID_GAP_B) & (scatlas_top_state != "Unresolved")] <- "Hybrid"
  names(scatlas_top_state) <- rownames(scatlas_group_max)

  scatlas_top_noncc <- colnames(scatlas_mp_adj[, scatlas_noncc_avail, drop = FALSE])[max.col(scatlas_mp_adj[, scatlas_noncc_avail, drop = FALSE], ties.method = "first")]
  names(scatlas_top_noncc) <- rownames(scatlas_mp_adj)
  scatlas_top_all <- colnames(scatlas_mp_adj)[max.col(scatlas_mp_adj, ties.method = "first")]
  names(scatlas_top_all) <- rownames(scatlas_mp_adj)

  scatlas_assign <- data.frame(
    cell = rownames(scatlas_mp_adj),
    top_noncc_mp = scatlas_top_noncc,
    top_all_mp = scatlas_top_all,
    topmp_state = scatlas_top_state,
    stringsAsFactors = FALSE
  )
  saveRDS(scatlas_assign, file.path(out_dir, "Auto_parse_topmp_assignments.rds"))
  write.csv(scatlas_assign, file.path(out_dir, "Auto_parse_topmp_assignments.csv"), row.names = FALSE)
} else {
  scatlas_ucell <- readRDS(file.path(out_dir, "Auto_parse_scATLAS_UCell_scores.rds"))
  scatlas_mp_adj <- readRDS(file.path(out_dir, "Auto_parse_scATLAS_mp_adj_noreg.rds"))
  scatlas_assign <- readRDS(file.path(out_dir, "Auto_parse_topmp_assignments.rds"))
  scatlas_top_noncc <- setNames(scatlas_assign$top_noncc_mp, scatlas_assign$cell)
  scatlas_top_all <- setNames(scatlas_assign$top_all_mp, scatlas_assign$cell)
  scatlas_top_state <- setNames(scatlas_assign$topmp_state, scatlas_assign$cell)
}
scatlas_group_order_pos <- sapply(scatlas_state_groups, function(mps) {
  positions <- match(mps, scatlas_tree_order)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
scatlas_ordered_group_names <- names(sort(scatlas_group_order_pos))
scatlas_state_order <- c(scatlas_ordered_group_names, "Unresolved", "Hybrid")

scatlas_col_map <- scatlas_mp_cols[scatlas_label_order]
missing_scatlas_cols <- setdiff(scatlas_label_order, names(scatlas_col_map)[!is.na(scatlas_col_map)])
if (length(missing_scatlas_cols) > 0) {
  scatlas_col_map[missing_scatlas_cols] <- hue_pal()(length(missing_scatlas_cols))
}
scatlas_col_map <- scatlas_col_map[scatlas_label_order]
scatlas_noncc_col_map <- scatlas_col_map[scatlas_noncc_label_order]

pdf(file.path(out_dir, "Auto_parse_topmp_abundance_scATLAS_9samples.pdf"), width = 16, height = 9, useDingbats = FALSE)
print(plot_abundance(label_mp(scatlas_top_noncc, scatlas_desc), cell_meta, scatlas_noncc_label_order, scatlas_noncc_col_map, sample_order_all, "scATLAS top non-CC MP abundance"))
print(plot_abundance(label_mp(scatlas_top_all, scatlas_desc), cell_meta, scatlas_label_order, scatlas_col_map, sample_order_all, "scATLAS top MP abundance"))
print(plot_abundance(scatlas_top_state, cell_meta, scatlas_state_order, scatlas_state_cols[scatlas_state_order], sample_order_all, "scATLAS topMP state-group abundance"))
dev.off()

####################
# PDO-pipeline scoring and exact Approach B/noreg state assignment
####################
pdo_geneNMF <- readRDS(
  "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds"
)
pdo_mp_genes <- pdo_geneNMF$metaprograms.genes
bad_mps <- which(pdo_geneNMF$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
coverage_tbl <- pdo_geneNMF$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
pdo_mp_genes <- pdo_mp_genes[!names(pdo_mp_genes) %in% c(bad_mp_names, low_coverage_mps)]
pdo_retained <- names(pdo_mp_genes)

pdo_tree_order <- pdo_geneNMF$programs.clusters[pdo_geneNMF$programs.tree$order]
pdo_tree_order <- paste0("MP", rev(unique(pdo_tree_order)))
pdo_tree_order <- pdo_tree_order[pdo_tree_order %in% pdo_retained]
pdo_cc <- intersect(c("MP6", "MP7", "MP1", "MP3"), pdo_retained)
pdo_noncc <- setdiff(pdo_retained, pdo_cc)

pdo_desc_full <- c(
  "MP6" = "MP6_G2M Cell Cycle",
  "MP7" = "MP7_DNA repair",
  "MP5" = "MP5_MYC-related Proliferation",
  "MP1" = "MP1_G2M checkpoint",
  "MP3" = "MP3_G1S Cell Cycle",
  "MP8" = "MP8_Columnar Progenitor",
  "MP10" = "MP10_Inflammatory Stress Epi.",
  "MP9" = "MP9_ECM Remodeling Epi.",
  "MP4" = "MP4_Intestinal Metaplasia"
)
pdo_desc <- sub("^MP[0-9]+_", "", pdo_desc_full)
names(pdo_desc) <- names(pdo_desc_full)
pdo_state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "Stress-adaptive" = c("MP10", "MP9"),
  "SMG-like Metaplasia" = c("MP8")
)
pdo_state_groups <- lapply(pdo_state_groups, function(mps) mps[mps %in% pdo_retained])
pdo_state_groups <- pdo_state_groups[sapply(pdo_state_groups, length) > 0]
group_order_pos <- sapply(pdo_state_groups, function(mps) {
  positions <- match(mps, pdo_tree_order)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
pdo_ordered_group_names <- names(sort(group_order_pos))
pdo_state_order <- c(pdo_ordered_group_names, "Unresolved", "Hybrid")
pdo_state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)
pdo_state_cols <- pdo_state_cols[names(pdo_state_cols) %in% pdo_state_order]
pdo_mp_cols <- c(
  "MP6_G2M Cell Cycle" = "#E78AC3",
  "MP7_DNA repair" = "#999999",
  "MP5_MYC-related Proliferation" = "#E41A1C",
  "MP1_G2M checkpoint" = "#B3B3B3",
  "MP3_G1S Cell Cycle" = "#8DA0CB",
  "MP8_Columnar Progenitor" = "#FF7F00",
  "MP10_Inflammatory Stress Epi." = "#984EA3",
  "MP9_ECM Remodeling Epi." = "#A6D854",
  "MP4_Intestinal Metaplasia" = "#4DAF4A"
)

if (need_pdo_scoring) {
  message("Scoring PDO-pipeline UCell signatures")
  pdo_mp_genes <- lapply(pdo_mp_genes, intersect, rownames(counts_all))
  pdo_mp_genes <- pdo_mp_genes[lengths(pdo_mp_genes) > 0]
  pdo_ucell <- UCell::ScoreSignatures_UCell(
    matrix = counts_all,
    features = pdo_mp_genes,
    maxRank = 1500,
    chunk.size = 1000,
    ncores = 1,
    force.gc = TRUE
  )
  pdo_ucell <- as.data.frame(pdo_ucell)
  colnames(pdo_ucell) <- sub("_UCell$", "", colnames(pdo_ucell))
  pdo_ucell <- as.matrix(pdo_ucell[, intersect(pdo_retained, colnames(pdo_ucell)), drop = FALSE])

  retained_in_ucell <- intersect(pdo_retained, colnames(pdo_ucell))
  cc_in_ucell <- intersect(pdo_cc, colnames(pdo_ucell))
  non_cc_in_ucell <- intersect(pdo_noncc, colnames(pdo_ucell))
  Y_use <- pdo_ucell[, non_cc_in_ucell, drop = FALSE]
  pdo_mp_adj_noncc <- z_normalise(Y_use, sample_var, study_var)
  cc_raw <- as.matrix(pdo_ucell[, cc_in_ucell, drop = FALSE])
  pdo_mp_adj_cc <- z_normalise(cc_raw, sample_var, study_var)
  pdo_mp_adj_all <- cbind(pdo_mp_adj_noncc, pdo_mp_adj_cc)

  group_max <- sapply(pdo_state_groups, function(mps) {
    mps_avail <- intersect(mps, colnames(pdo_mp_adj_noncc))
    if (length(mps_avail) == 1) return(as.numeric(pdo_mp_adj_noncc[, mps_avail]))
    apply(pdo_mp_adj_noncc[, mps_avail, drop = FALSE], 1, max)
  })
  group_max <- as.matrix(group_max)
  rownames(group_max) <- rownames(pdo_mp_adj_noncc)

  THRESHOLD <- parse_state_definition$min_group_score
  HYBRID_GAP_B <- parse_state_definition$hybrid_gap
  best_group_idx <- max.col(group_max, ties.method = "first")
  best_group_val <- apply(group_max, 1, max)
  pdo_state <- names(pdo_state_groups)[best_group_idx]
  pdo_state[best_group_val < THRESHOLD] <- "Unresolved"
  sorted_groups <- t(apply(group_max, 1, sort, decreasing = TRUE))
  gap <- sorted_groups[, 1] - sorted_groups[, 2]
  pdo_state[(gap < HYBRID_GAP_B) & (pdo_state != "Unresolved")] <- "Hybrid"
  names(pdo_state) <- rownames(group_max)

  pdo_top_noncc <- colnames(pdo_mp_adj_noncc)[max.col(pdo_mp_adj_noncc, ties.method = "first")]
  names(pdo_top_noncc) <- rownames(pdo_mp_adj_noncc)
  pdo_top_all <- colnames(pdo_mp_adj_all)[max.col(pdo_mp_adj_all, ties.method = "first")]
  names(pdo_top_all) <- rownames(pdo_mp_adj_all)
} else {
  message("Loading cached PDO-pipeline scores and assignments")
  pdo_ucell <- readRDS(file.path(out_dir, "Auto_parse_PDOpipeline_UCell_scores.rds"))
  pdo_mp_adj_noncc <- readRDS(file.path(out_dir, "Auto_parse_PDOpipeline_mp_adj_noreg.rds"))
  pdo_mp_adj_all <- readRDS(file.path(out_dir, "Auto_parse_PDOpipeline_mp_adj_all_noreg.rds"))
  pdo_assign_cached <- readRDS(file.path(out_dir, "Auto_parse_PDOpipeline_topmp_assignments.rds"))
  pdo_state <- setNames(pdo_assign_cached$pdo_state, pdo_assign_cached$cell)
  pdo_top_noncc <- setNames(pdo_assign_cached$top_noncc_mp, pdo_assign_cached$cell)
  pdo_top_all <- setNames(pdo_assign_cached$top_all_mp, pdo_assign_cached$cell)
}

pdo_mp_order <- ordered_block(pdo_cc, pdo_tree_order)
for (state in names(pdo_state_groups)) {
  pdo_mp_order <- c(pdo_mp_order, ordered_block(intersect(pdo_state_groups[[state]], pdo_retained), pdo_tree_order))
}
pdo_mp_order <- unique(c(pdo_mp_order, setdiff(pdo_retained, pdo_mp_order)))
pdo_noncc_order <- setdiff(pdo_mp_order, pdo_cc)
pdo_label_order <- unname(pdo_desc_full[pdo_mp_order])
pdo_noncc_label_order <- unname(pdo_desc_full[pdo_noncc_order])
pdo_col_map <- pdo_mp_cols[pdo_label_order]
missing_pdo_cols <- setdiff(pdo_label_order, names(pdo_col_map)[!is.na(pdo_col_map)])
if (length(missing_pdo_cols) > 0) {
  pdo_col_map[missing_pdo_cols] <- hue_pal()(length(missing_pdo_cols))
}
pdo_col_map <- pdo_col_map[pdo_label_order]
pdo_noncc_col_map <- pdo_col_map[pdo_noncc_label_order]

pdo_assign <- data.frame(
  cell = rownames(pdo_mp_adj_all),
  original_cell = cell_meta[rownames(pdo_mp_adj_all), "original_cell"],
  sample = cell_meta[rownames(pdo_mp_adj_all), "sample"],
  top_noncc_mp = pdo_top_noncc,
  top_noncc_label = unname(pdo_desc_full[pdo_top_noncc]),
  top_all_mp = pdo_top_all,
  top_all_label = unname(pdo_desc_full[pdo_top_all]),
  pdo_state = pdo_state,
  stringsAsFactors = FALSE
)

saveRDS(pdo_ucell, file.path(out_dir, "Auto_parse_PDOpipeline_UCell_scores.rds"))
saveRDS(pdo_mp_adj_noncc, file.path(out_dir, "Auto_parse_PDOpipeline_mp_adj_noreg.rds"))
saveRDS(pdo_mp_adj_all, file.path(out_dir, "Auto_parse_PDOpipeline_mp_adj_all_noreg.rds"))
saveRDS(pdo_assign, file.path(out_dir, "Auto_parse_PDOpipeline_topmp_assignments.rds"))
write.csv(pdo_assign, file.path(out_dir, "Auto_parse_PDOpipeline_topmp_assignments.csv"), row.names = FALSE)

pdo_abundance_summary <- bind_rows(
  make_prop_data(label_from_map(pdo_top_noncc, pdo_desc_full), cell_meta, pdo_noncc_label_order, sample_order_all) %>% mutate(view = "PDO pipeline top non-CC MP"),
  make_prop_data(label_from_map(pdo_top_all, pdo_desc_full), cell_meta, pdo_label_order, sample_order_all) %>% mutate(view = "PDO pipeline top all MP"),
  make_prop_data(pdo_state, cell_meta, pdo_state_order, sample_order_all) %>% mutate(view = "PDO pipeline Approach B state")
)
write.csv(pdo_abundance_summary, file.path(out_dir, "Auto_parse_PDOpipeline_abundance_summary.csv"), row.names = FALSE)

pdf(file.path(out_dir, "Auto_parse_topmp_abundance_PDOpipeline_9samples.pdf"), width = 16, height = 9, useDingbats = FALSE)
print(plot_abundance(label_from_map(pdo_top_noncc, pdo_desc_full), cell_meta, pdo_noncc_label_order, pdo_noncc_col_map, sample_order_all, "PDO-pipeline top non-CC MP abundance"))
print(plot_abundance(label_from_map(pdo_top_all, pdo_desc_full), cell_meta, pdo_label_order, pdo_col_map, sample_order_all, "PDO-pipeline top MP abundance"))
print(plot_abundance(pdo_state, cell_meta, pdo_state_order, pdo_state_cols[pdo_state_order], sample_order_all, "PDO-pipeline Approach B noreg state abundance"))
dev.off()

####################
# Parse-derived MP scoring and abundance
####################
if (need_parse_scoring) {
  message("Scoring Parse-derived MP signatures on all samples")
  parse_features <- parse_mp_genes[parse_retained]
  parse_features <- lapply(parse_features, intersect, rownames(counts_all))
  parse_ucell_all <- UCell::ScoreSignatures_UCell(
    matrix = counts_all,
    features = parse_features,
    maxRank = 1500, chunk.size = 1000, ncores = 1, force.gc = TRUE
  )
  colnames(parse_ucell_all) <- sub("_UCell$", "", colnames(parse_ucell_all))
  parse_ucell_save <- as.data.frame(parse_ucell_all)
  rownames(parse_ucell_save) <- sub("__", "_", rownames(parse_ucell_save))
  saveRDS(parse_ucell_save, parse_ucell_path)
  parse_ucell <- as.matrix(parse_ucell_all)
} else {
  message("Loading Parse-derived UCell scores from: ", parse_ucell_path)
  parse_ucell <- readRDS(parse_ucell_path)
  rownames(parse_ucell) <- sub("^([^_]+)_", "\\1__", rownames(parse_ucell))
  parse_ucell <- as.matrix(parse_ucell[, intersect(parse_retained, colnames(parse_ucell)), drop = FALSE])
}

parse_top_all <- colnames(parse_ucell)[max.col(parse_ucell, ties.method = "first")]
names(parse_top_all) <- rownames(parse_ucell)
parse_tree_order <- paste0("MP", parse_geneNMF$programs.clusters[parse_geneNMF$programs.tree$order])
parse_tree_order <- unique(parse_tree_order[parse_tree_order %in% parse_retained])
parse_col_map <- setNames(hue_pal()(length(parse_tree_order)), parse_tree_order)

pdf(file.path(out_dir, "Auto_parse_topmp_abundance_ParseMP_9samples.pdf"), width = 16, height = 9, useDingbats = FALSE)
print(plot_abundance(parse_top_all, cell_meta, parse_tree_order, parse_col_map, sample_order_all, "Parse-derived top MP abundance"))
dev.off()

####################
# Dual MP activity PDFs
####################
write_activity_pdf <- function(sample_order, output_file, suffix_label) {
  scatlas_activity <- make_activity_plot(
    scatlas_ucell,
    cell_meta,
    intersect(scatlas_mp_order, colnames(scatlas_ucell)),
    scatlas_desc,
    sample_order,
    paste0("scATLAS MP activity by sample - ", suffix_label)
  )
  pdo_activity <- make_activity_plot(
    pdo_ucell,
    cell_meta,
    intersect(pdo_mp_order, colnames(pdo_ucell)),
    pdo_desc,
    sample_order,
    paste0("PDO-pipeline MP activity by sample - ", suffix_label)
  )
  parse_activity <- make_activity_plot(
    parse_ucell,
    cell_meta,
    intersect(parse_tree_order, colnames(parse_ucell)),
    setNames(rep(NA_character_, length(parse_tree_order)), parse_tree_order),
    sample_order,
    paste0("Parse-derived MP activity by sample - ", suffix_label)
  )
  write.csv(
    scatlas_activity$stats,
    sub("\\.pdf$", "_scATLAS_stats.csv", output_file),
    row.names = FALSE
  )
  write.csv(
    pdo_activity$stats,
    sub("\\.pdf$", "_PDOpipeline_stats.csv", output_file),
    row.names = FALSE
  )
  write.csv(
    parse_activity$stats,
    sub("\\.pdf$", "_ParseMP_stats.csv", output_file),
    row.names = FALSE
  )
  pdf(output_file, width = 20, height = 9, useDingbats = FALSE)
  print(scatlas_activity$plot)
  print(pdo_activity$plot)
  print(parse_activity$plot)
  dev.off()
}

write_activity_pdf(
  sample_order_all,
  file.path(out_dir, "Auto_parse_mp_activity_boxplots_scATLAS_and_PDOpipeline_include_PDO_SUR1090.pdf"),
  "including PDO and SUR1090 samples"
)
write_activity_pdf(
  sample_order_excluding_pdo_sur1090,
  file.path(out_dir, "Auto_parse_mp_activity_boxplots_scATLAS_and_PDOpipeline_exclude_PDO_SUR1090.pdf"),
  "excluding PDO and SUR1090 samples"
)

script_run_status <- "success"
message("Done. Outputs written to: ", out_dir)
