####################
# Description: Candidate-gene expression across six stages, eight pEMT MP activity
# distributions, baseline-relative high-cell fractions, and T4 high/low contrasts.
# Inputs (all relative to the canonical live project root):
# analysis/metaprograms/centred/T2T4hgih functional annotation.xlsx
# parse_outs/by_samples/{T0,T1,T2,T4,R4,eR4}/Auto_<sample>_final.rds
# parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/
# Auto_parse_highres_T2T4_{UCell_scores,cell_metadata,selected_mp_genes}_nMP117.rds
# Outputs: parse_outs/centred/candidate_pemt_validation/{intermediate,tables,
# figures/Auto_final,figures/Auto_supporting,reports}; shared run-summary logs.
# Cache/replot: --mode=analyse builds a compact persistent replot bundle; refuses
# to replace an existing bundle. --mode=replot reads this bundle only and rewrites
# this workflow's derived tables/figures/reports. --mode=figure1 refreshes only
# Figure 1, its source table and index. --mode=figure4 refreshes the two-page
# Figure 04 PDF, PNG companions, association tables/report and index only.
# --mode=scatter refreshes supporting per-cell scatter pages and summaries.
# Methodology: analysis/methodology/metaprograms/centred/candidate_pemt_validation_methodology.md
# Downstream: terminal exploratory candidate-prioritisation report, not validation.
# Execution registry: gnmf; PBS 4 CPUs, 48 GB, 04:00:00; from the live project root:
# qsub job_scripts/metaprograms/centred/parse_centred_candidate_pemt_validation.sh
# Direct compute-node command: Rscript analysis/metaprograms/centred/parse_centred_candidate_pemt_validation.R --mode=analyse
####################
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(UCell)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

####################
# Figure 1 uses unweighted row z-scores of six stage means. Labels retain the
# unscaled mean log1p(CP10k), not raw counts or z-scores. Stop on constant rows.
plot_candidate_expression_heatmap <- function(expression_summary, out, genes, stages) {
  plot_data <- expression_summary |>
    group_by(gene) |>
    mutate(stage_mean = mean(mean_log1p_cp10k), stage_sd = sd(mean_log1p_cp10k),
           row_z_score = (mean_log1p_cp10k - stage_mean) / stage_sd) |>
    ungroup()
  if (any(!is.finite(plot_data$row_z_score)) || any(plot_data$stage_sd <= 0)) {
    stop("Cannot z-score a candidate with zero/nonfinite variation across stages")
  }
  if (nrow(plot_data) != length(genes) * length(stages)) stop("Incomplete expression heatmap grid")
  plot_data$stage <- factor(plot_data$stage, stages)
  plot_data$gene <- factor(plot_data$gene, rev(genes))
  z_limit <- max(abs(plot_data$row_z_score))
  plot_data$label_colour <- ifelse(abs(plot_data$row_z_score) > 0.65 * z_limit, "white", "black")
  write.csv(plot_data, file.path(out, "tables/candidate_expression_heatmap_source.csv"), row.names = FALSE)
  p <- ggplot(plot_data, aes(stage, gene, fill = row_z_score)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", mean_log1p_cp10k), colour = label_colour), size = 4) +
    scale_colour_identity() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-z_limit, z_limit), name = "Row z-score") +
    scale_x_discrete(expand = c(0, 0), position = "top") +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = "Candidate expression trends across treatment and recovery",
         subtitle = "Colour: gene-wise z-score across six stages\nNumbers: mean log1p(CP10k) expression, including zero-expression cells") +
    theme_minimal(base_size = 14, base_family = "sans") +
    theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"),
          axis.text.x = element_text(face = "bold"), plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 12, margin = margin(b = 12)))
  path <- file.path(out, "figures/Auto_final/01_candidate_expression_by_stage")
  ggsave(paste0(path, ".pdf"), p, device = grDevices::cairo_pdf, width = 12, height = 10)
  ggsave(paste0(path, ".png"), p, device = grDevices::png, type = "cairo", width = 12, height = 10, dpi = 160)
}
####################

####################
# Figure 04 pairs the original binary contrast with continuous T4 associations.
# Each tested candidate is omitted from its MP score on the correlation page.
plot_t4_candidate_association <- function(bundle, out, genes, mps) {
  t4 <- which(bundle$meta$sample == "T4")
  t0 <- which(bundle$meta$sample == "T0")
  ids <- bundle$meta$cell[t4]
  if (!identical(rownames(bundle$expression), bundle$meta$cell) ||
      !identical(rownames(bundle$scores), bundle$meta$cell) ||
      !all(ids %in% rownames(bundle$rescored))) stop("T4 source data are not aligned")
  cutoff <- apply(bundle$scores[t0, mps, drop = FALSE], 2, quantile, .95, type = 7)
  flags <- sweep(bundle$scores[t4, mps, drop = FALSE], 2, cutoff, ">")
  flags <- cbind(flags, Any_pEMT_MP = rowSums(flags) > 0)
  binary <- bind_rows(lapply(colnames(flags), function(mp) {
    h <- flags[, mp]
    if (sum(h) < 3 || sum(!h) < 3) stop("Insufficient T4 binary contrast cells: ", mp)
    data.frame(MP = mp, gene = genes,
      delta_mean_log1p = colMeans(bundle$expression[t4[h], , drop = FALSE]) -
                          colMeans(bundle$expression[t4[!h], , drop = FALSE]))
  }))
  binary$gene <- factor(binary$gene, rev(genes))
  binary$MP <- factor(binary$MP, c(mps, "Any_pEMT_MP"))
  base_theme <- theme_classic(base_size = 13, base_family = "sans") +
    theme(axis.text = element_text(colour = "black"), plot.title = element_text(face = "bold"))
  page1 <- ggplot(binary, aes(MP, gene, fill = delta_mean_log1p)) +
    geom_tile(colour = "white") +
    scale_fill_gradient2(low = "#3b6b9a", mid = "white", high = "#bd5b4a", midpoint = 0,
                         name = "High minus low\nmean log1p(CP10k)") +
    labs(x = NULL, y = NULL, title = "T4 candidate expression: MP-high versus MP-low",
         subtitle = "Original signatures; descriptive cell-level association, with leave-one-gene-out checks in supporting outputs") +
    base_theme + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  associations <- bind_rows(lapply(genes, function(gene) {
    x <- bundle$expression[t4, gene]
    if (sd(x) == 0) stop("Constant T4 candidate expression: ", gene)
    bind_rows(lapply(mps, function(mp) {
      y <- bundle$rescored[ids, paste0(mp, "_loo_", gene)]
      original <- bundle$scores[t4, mp]
      if (any(!is.finite(c(x, y, original))) || sd(y) == 0 || sd(original) == 0) {
        stop("Nonfinite/constant expression or MP scores: ", gene, " / ", mp)
      }
      data.frame(gene = gene, MP = mp, n_cells = length(t4),
        pct_gene_detected = 100 * mean(x > 0),
        spearman_rho_leave_one_out = cor(x, y, method = "spearman"),
        spearman_rho_original = cor(x, original, method = "spearman"))
    }))
  }))
  summary <- associations |> group_by(gene) |>
    summarise(median_rho = median(spearman_rho_leave_one_out),
      min_rho = min(spearman_rho_leave_one_out), max_rho = max(spearman_rho_leave_one_out),
      n_positive_mps = sum(spearman_rho_leave_one_out > 0),
      pct_gene_detected = first(pct_gene_detected), .groups = "drop") |>
    arrange(desc(median_rho), gene)
  write.csv(associations, file.path(out, "tables/T4_candidate_continuous_mp_associations.csv"), row.names = FALSE)
  write.csv(summary, file.path(out, "tables/T4_candidate_continuous_association_summary.csv"), row.names = FALSE)
  plot_data <- bind_rows(
    transmute(associations, gene, MP, rho = spearman_rho_leave_one_out),
    transmute(summary, gene, MP = "Median\n(8 MPs)", rho = median_rho))
  plot_data$gene <- factor(plot_data$gene, rev(summary$gene))
  plot_data$MP <- factor(plot_data$MP, c(mps, "Median\n(8 MPs)"))
  page2 <- ggplot(plot_data, aes(MP, gene, fill = rho)) +
    geom_tile(colour = "white", linewidth = .4) +
    geom_vline(xintercept = 8.5, colour = "grey35", linewidth = .5) +
    geom_text(aes(label = sprintf("%.2f", rho), colour = abs(rho) >= .65), size = 4) +
    scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                         limits = c(-1, 1), breaks = c(-1, -.5, 0, .5, 1), name = "Spearman rho") +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = "T4 continuous gene-MP coactivation",
      subtitle = paste0("Spearman correlation across ", format(length(t4), big.mark = ","),
        " cells; tested gene omitted from each MP signature\nRows ranked by median correlation across eight MPs; zero-expression cells included"),
      caption = "Red: expression rises with MP activity; white: weak monotonic association; blue: inverse association.\nDescriptive cell-level associations; sequencing depth, shared biology and zero counts can influence correlations.") +
    theme_minimal(base_size = 14, base_family = "sans") +
    theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"),
      axis.text.x = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 12),
      plot.caption = element_text(hjust = 0, size = 10, margin = margin(t = 10)))
  path <- file.path(out, "figures/Auto_final/04_T4_candidate_high_vs_low")
  grDevices::cairo_pdf(paste0(path, ".pdf"), width = 13, height = 10, onefile = TRUE)
  print(page1)
  print(page2)
  dev.off()
  ggsave(paste0(path, ".png"), page1, device = grDevices::png, type = "cairo", width = 13, height = 10, dpi = 160)
  ggsave(paste0(path, "_page2.png"), page2, device = grDevices::png, type = "cairo", width = 13, height = 10, dpi = 160)
  # Independent rank-Pearson identity checks the reported Spearman correlations.
  check <- vapply(seq_len(nrow(associations)), function(i) {
    r <- associations[i, ]
    x <- bundle$expression[t4, r$gene]
    y <- bundle$rescored[ids, paste0(r$MP, "_loo_", r$gene)]
    abs(cor(rank(x, ties.method = "average"), rank(y, ties.method = "average")) - r$spearman_rho_leave_one_out)
  }, numeric(1))
  if (max(check) > 1e-12 || nrow(associations) != length(genes) * length(mps)) stop("Association validation failed")
  writeLines(c("status: PASS", "Figure 04: two PDF pages, with page-specific PNG companions.",
    paste0("T4 cells: ", length(t4)), paste0("Gene-MP correlations: ", nrow(associations)),
    paste0("Maximum rank-Pearson / Spearman discrepancy: ", max(check)),
    "Each correlation uses the tested gene's leave-one-out MP score; ties receive average ranks.",
    paste0("Top five by median rho: ", paste(head(summary$gene, 5), collapse = ", "))),
    file.path(out, "reports/figure4_continuous_association_validation.txt"))
  writeLines(c("# Figure 04: continuous T4 coactivation", "",
    "Page 1 retains the original MP-high versus MP-low mean expression differences. Page 2 summarises continuous Spearman correlations for every candidate gene against each of eight pEMT MPs, using all T4 cells including zero-expression cells. The candidate being tested is removed from that MP signature. No high/low cutoff is used on page 2.", "",
    "Positive rho indicates that higher expression accompanies higher MP scores; near-zero rho indicates little monotonic association; negative rho indicates an inverse trend. Near-zero correlation is not proof of no biological relationship, especially for sparse genes or nonmonotonic responses.", "",
    "Rows are ordered by median rho across the eight MPs. The final column is that median, a descriptive summary of correlated programs rather than an independent ninth MP. These correlations describe RNA coexpression, not causal activation or independent functional validation. Sequencing depth and other cell properties can confound associations; rho is not adjusted for these variables.", "",
    "Tables contain all 168 pairs with original-score and leave-one-out-score rho, cell counts and gene detection rates, plus one 21-gene summary table. No significance stars are used; effect magnitude is the visual emphasis."),
    file.path(out, "reports/figure4_continuous_association.md"))
  invisible(summary)
}
####################

####################
# Supporting per-cell views: one candidate per page, all eight MPs on that page.
# Every T4 cell is plotted; decile means describe the trend without a fitted model.
plot_t4_candidate_scatter <- function(bundle, out, genes, mps) {
  t4 <- which(bundle$meta$sample == "T4")
  ids <- bundle$meta$cell[t4]
  if (!identical(rownames(bundle$expression), bundle$meta$cell) ||
      !all(ids %in% rownames(bundle$rescored))) stop("Scatter input cells are not aligned")
  pdf_path <- file.path(out, "figures/Auto_supporting/T4_candidate_per_cell_pemt_scatter.pdf")
  grDevices::cairo_pdf(pdf_path, width = 16, height = 9, onefile = TRUE)
  device_open <- TRUE
  on.exit(if (device_open) dev.off(), add = TRUE)
  bins_all <- list()
  scatter_index <- list()
  for (gene in genes) {
    message("Per-cell scatter page: ", gene)
    dat <- bind_rows(lapply(mps, function(mp) {
      y <- bundle$expression[t4, gene]
      x <- bundle$rescored[ids, paste0(mp, "_loo_", gene)]
      if (any(!is.finite(c(x, y))) || sd(x) == 0 || sd(y) == 0) stop("Invalid scatter inputs: ", gene, "/", mp)
      data.frame(cell = ids, gene = gene, MP = mp, mp_score = x, expression = y)
    }))
    panel_stats <- dat |> group_by(MP) |>
      summarise(n_cells = n(), rho = cor(mp_score, expression, method = "spearman"), .groups = "drop")
    if (any(panel_stats$n_cells != length(t4))) stop("Scatter panel missing T4 cells")
    bins <- bind_rows(lapply(mps, function(mp) {
      d <- dat[dat$MP == mp, ]
      breaks <- unique(as.numeric(quantile(d$mp_score, seq(0, 1, .1), type = 7)))
      if (length(breaks) < 3L) stop("Too few distinct MP values for score-decile summaries")
      d$score_bin <- cut(d$mp_score, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      d |> group_by(gene, MP, score_bin) |>
        summarise(n_cells = n(), mean_mp_score = mean(mp_score),
          mean_expression = mean(expression), .groups = "drop")
    }))
    if (sum(bins$n_cells) != length(t4) * length(mps)) stop("Score-bin cell accounting mismatch")
    bins_all[[gene]] <- bins
    scatter_index[[gene]] <- data.frame(page = match(gene, genes), gene = gene,
      n_cells_per_panel = length(t4), n_MP_panels = length(mps),
      median_rho = median(panel_stats$rho))
    labels <- setNames(paste0(panel_stats$MP, "  |  rho = ", sprintf("%.2f", panel_stats$rho)), panel_stats$MP)
    dat$MP <- factor(dat$MP, mps)
    bins$MP <- factor(bins$MP, mps)
    p <- ggplot(dat, aes(mp_score, expression)) +
      geom_point(colour = "#355f7d", alpha = .16, size = .55) +
      geom_line(data = bins, aes(mean_mp_score, mean_expression), inherit.aes = FALSE,
                colour = "#c65f21", linewidth = .9) +
      geom_point(data = bins, aes(mean_mp_score, mean_expression), inherit.aes = FALSE,
                 colour = "#c65f21", size = 1.8) +
      facet_wrap(~MP, ncol = 4, scales = "free_x", labeller = labeller(MP = labels)) +
      labs(x = "MP UCell score (tested gene omitted)",
        y = paste0(gene, " expression: log1p(CP10k)"),
        title = paste0(gene, ": per-cell expression versus pEMT MP activity at T4"),
        subtitle = paste0(format(length(t4), big.mark = ","),
          " cells in every panel; orange points/line = mean expression within MP-score deciles"),
        caption = "All T4 cells included, with no jitter or subsampling; overlapping points share positions. Panel labels show Spearman rho.\nZero-expression cells are retained. Orange lines are descriptive bin means, not a fitted causal or regression model.") +
      theme_classic(base_size = 15, base_family = "sans") +
      theme(axis.text = element_text(colour = "black", size = 12),
        strip.background = element_rect(fill = "grey95", colour = "grey55"),
        strip.text = element_text(size = 13, face = "bold"),
        plot.title = element_text(face = "bold", size = 19),
        plot.subtitle = element_text(size = 14),
        plot.caption = element_text(size = 11, hjust = 0, margin = margin(t = 12)),
        panel.spacing = grid::unit(1.2, "lines"))
    print(p)
    if (identical(gene, genes[1])) {
      ggsave(file.path(out, "figures/Auto_supporting/T4_candidate_per_cell_pemt_scatter_preview.png"),
        p, device = grDevices::png, type = "cairo", width = 16, height = 9, dpi = 150)
    }
  }
  dev.off()
  device_open <- FALSE
  write.csv(bind_rows(bins_all), file.path(out, "tables/T4_candidate_scatter_decile_means.csv"), row.names = FALSE)
  write.csv(bind_rows(scatter_index), file.path(out, "tables/T4_candidate_scatter_page_index.csv"), row.names = FALSE)
  writeLines(c("status: PASS", paste0("Pages: ", length(genes)),
    paste0("Panels: ", length(genes) * length(mps)), paste0("T4 cells per panel: ", length(t4)),
    paste0("Total plotted cell-pair observations: ", length(genes) * length(mps) * length(t4)),
    "Every score-decile summary accounts for all T4 cells; no jitter or downsampling.",
    "MP scores omit the tested gene, matching Figure 04 page 2."),
    file.path(out, "reports/per_cell_scatter_validation.txt"))
  invisible(NULL)
}
####################

main <- function() {
  root <- "/rds/general/project/spatialtranscriptomics/live/Parse_Pipeline"
  Sys.setenv(AUTO_PARSE_ROOT_DIR = root)
  # This workflow runs in a four-CPU PBS allocation, not an R package check.
  # An inherited package-check flag otherwise caps BiocParallel at two workers.
  if (mode_check <- nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
    if (!nzchar(Sys.getenv("PBS_JOBID"))) stop("Inherited core-check flag: run this analysis through PBS")
    Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  }
  setwd(root)
  source("analysis/common/parse_pipeline_config.R")
  source("analysis/common/parse_pipeline_logging.R")
  mode <- sub("^--mode=", "", commandArgs(TRUE))
  if (length(mode) != 1L || !mode %in% c("analyse", "replot", "figure1", "figure4", "scatter")) stop("Specify --mode=analyse, --mode=replot, --mode=figure1, --mode=figure4 or --mode=scatter")
  stages <- c("T0", "T1", "T2", "T4", "R4", "eR4")
  mps <- c("MP48", "MP50", "MP60", "MP109", "MP16", "MP74", "MP83", "MP9")
  out <- file.path(root, "parse_outs/centred/candidate_pemt_validation")
  dirs <- file.path(out, c("intermediate", "tables", "figures/Auto_final", "figures/Auto_supporting", "reports"))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  bundle_path <- file.path(out, "intermediate/candidate_pemt_replot_bundle.rds")
  workbook <- file.path(root, "analysis/metaprograms/centred/T2T4hgih functional annotation.xlsx")
  upstream <- file.path(root, "parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter")
  paths <- setNames(file.path(upstream, paste0("Auto_parse_highres_T2T4_", c("UCell_scores", "cell_metadata", "selected_mp_genes"), "_nMP117.rds")), c("scores", "metadata", "genes"))
  sample_paths <- setNames(file.path(root, "parse_outs/by_samples", stages, paste0("Auto_", stages, "_final.rds")), stages)
  run <- parse_start_run("parse_centred_candidate_pemt_validation", parameters = list(mode = mode, MPs = mps, threshold = "T0 per-MP q95; strict >; type 7", sensitivity = "T0 q90,q97.5,q99; leave tested candidate out"), input_files = if (mode == "analyse") c(workbook, paths, sample_paths) else bundle_path, output_files = out, reused_cached = mode != "analyse")
  status <- "failed"
  on.exit(parse_finish_run(run, status = status), add = TRUE)
  write_table <- function(x, name) write.csv(x, file.path(out, "tables", paste0(name, ".csv")), row.names = FALSE, na = "NA")
  if (mode == "analyse") {
    if (file.exists(bundle_path)) stop("Bundle exists: use --mode=replot; do not overwrite source data.")
    required <- c(workbook, paths, sample_paths)
    if (any(!file.exists(required))) stop("Missing inputs: ", paste(required[!file.exists(required)], collapse = ", "))
    candidates <- as.data.frame(readxl::read_excel(workbook, sheet = "Candidate genes"))
    candidates <- candidates[!is.na(candidates$Gene) & nzchar(trimws(candidates$Gene)), , drop = FALSE]
    candidates$Gene <- trimws(candidates$Gene)
    genes <- candidates$Gene
    if (anyDuplicated(genes) || length(genes) != 21L) stop("Expected 21 unique candidate genes")
    header <- unlist(readxl::read_excel(workbook, sheet = "pEMT-Plasticity", range = "A2:H2", col_names = FALSE), use.names = FALSE)
    workbook_mps <- sub(" .*", "", trimws(header))
    if (!setequal(mps, workbook_mps)) stop("Workbook pEMT MPs do not match the agreed eight MPs")
    scores <- as.matrix(readRDS(paths[["scores"]]))
    meta <- as.data.frame(readRDS(paths[["metadata"]]))
    signatures <- readRDS(paths[["genes"]])[mps]
    if (any(lengths(signatures) == 0L) || !all(mps %in% colnames(scores))) stop("Missing MP signatures/scores")
    workbook_genes <- readxl::read_excel(workbook, sheet = "pEMT-Plasticity", skip = 3, col_names = FALSE)
    for (i in seq_along(workbook_mps)) {
      listed <- trimws(as.character(workbook_genes[[i]]))
      listed <- listed[!is.na(listed) & nzchar(listed)]
      if (!identical(listed, as.character(signatures[[workbook_mps[i]]]))) stop("Workbook signature/rank mismatch: ", workbook_mps[i])
    }
    if (!all(c("cell", "sample", "original_cell") %in% names(meta))) stop("Unexpected cached metadata schema")
    if (anyDuplicated(meta$cell) || anyDuplicated(rownames(scores))) stop("Duplicate cell IDs")
    meta <- meta[meta$sample %in% stages, , drop = FALSE]
    if (!setequal(unique(meta$sample), stages) || !all(meta$cell %in% rownames(scores))) stop("Missing stage/cell scores")
    scores <- scores[meta$cell, mps, drop = FALSE]
    if (any(!is.finite(scores))) stop("Nonfinite UCell scores")
    expression <- matrix(NA_real_, nrow(meta), length(genes), dimnames = list(meta$cell, genes))
    detected <- expression
    metadata_inventory <- list()
    feature_sets <- list()
    sensitivity_counts <- list()
    meta$nCount_RNA_measured <- NA_real_
    meta$nFeature_RNA_measured <- NA_real_
    for (stage in stages) {
      message("Extracting candidate expression: ", stage)
      obj <- readRDS(sample_paths[[stage]])
      if (!"RNA" %in% names(obj@assays)) stop("Missing RNA assay: ", stage)
      counts <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
      if (!inherits(counts, "sparseMatrix")) stop("Expected sparse RNA counts: ", stage)
      if (!all(genes %in% rownames(counts))) stop("Missing candidates: ", paste(setdiff(genes, rownames(counts)), collapse = ", "))
      idx <- which(meta$sample == stage)
      old <- meta$original_cell[idx]
      if (!setequal(old, colnames(counts))) stop("Cached and current cell populations differ: ", stage)
      counts <- counts[, old, drop = FALSE]
      md <- obj[[]][old, , drop = FALSE]
      if (!"orig.ident" %in% names(md) || !all(as.character(md$orig.ident) == stage)) stop("orig.ident mismatch: ", stage)
      lib <- Matrix::colSums(counts)
      if (any(lib <= 0)) stop("Zero-library cells")
      expression[idx, ] <- t(log1p(t(t(as.matrix(counts[genes, , drop = FALSE])) / lib) * 10000))
      detected[idx, ] <- t(as.matrix(counts[genes, , drop = FALSE] > 0))
      meta$nCount_RNA_measured[idx] <- lib
      meta$nFeature_RNA_measured[idx] <- Matrix::colSums(counts > 0)
      feature_sets[[stage]] <- rownames(counts)
      metadata_inventory[[stage]] <- data.frame(stage = stage, field = names(md), unique_values = vapply(md, function(x) length(unique(x)), integer(1)), examples = vapply(md, function(x) paste(head(unique(as.character(x)), 8), collapse = "; "), character(1)))
      if (stage %in% c("T0", "T4")) {
        colnames(counts) <- meta$cell[idx]
        sensitivity_counts[[stage]] <- counts
      }
      rm(obj, counts)
      gc()
    }
    if (any(!is.finite(expression))) stop("Nonfinite expression values")
    common_genes <- Reduce(intersect, feature_sets)
    full_signatures <- lapply(signatures, intersect, common_genes)
    excluded_signatures <- unlist(lapply(genes, function(gene) {
      setNames(lapply(full_signatures, setdiff, gene), paste0(mps, "_loo_", gene))
    }), recursive = FALSE)
    if (any(lengths(excluded_signatures) < 5L)) stop("Leave-one-out signature contains fewer than five genes")
    overlap <- bind_rows(lapply(mps, function(mp) data.frame(MP = mp, n_signature = length(signatures[[mp]]), n_scored = length(full_signatures[[mp]]), n_candidates = length(intersect(full_signatures[[mp]], genes)), candidate_genes = paste(intersect(full_signatures[[mp]], genes), collapse = "; "), n_remaining_if_all_candidates_removed = length(setdiff(full_signatures[[mp]], genes)))))
    # Recompute full and reduced signatures together to distinguish exclusion effects
    # from UCell version differences relative to the historical scores.
    features <- c(setNames(full_signatures, paste0(mps, "_full")), excluded_signatures)
    sensitivity_scores <- list()
    for (stage in c("T0", "T4")) {
      message("Scoring full/leave-one-gene-out sensitivity signatures: ", stage)
      mat <- sensitivity_counts[[stage]][common_genes, , drop = FALSE]
      sc <- as.matrix(UCell::ScoreSignatures_UCell(matrix = mat, features = features, maxRank = 1500, chunk.size = 500, ncores = 4, force.gc = TRUE))
      colnames(sc) <- sub("_UCell$", "", colnames(sc))
      if (!all(names(features) %in% colnames(sc)) || !setequal(rownames(sc), colnames(mat))) stop("Unexpected rescoring schema")
      sensitivity_scores[[stage]] <- sc[colnames(mat), names(features), drop = FALSE]
    }
    rescored <- do.call(rbind, sensitivity_scores)
    if (any(!is.finite(rescored))) stop("Nonfinite sensitivity scores")
    input_manifest <- data.frame(path = required, bytes = file.info(required)$size, modified = as.character(file.info(required)$mtime))
    bundle <- list(candidates = candidates, meta = meta, expression = expression, detected = detected, scores = scores, rescored = rescored, signatures = signatures, excluded_signatures = excluded_signatures, overlap = overlap, metadata_inventory = bind_rows(metadata_inventory), manifest = input_manifest, workbook_md5 = unname(tools::md5sum(workbook)), common_gene_count = length(common_genes), UCell_version = as.character(packageVersion("UCell")))
    saveRDS(bundle, bundle_path)
    rm(sensitivity_counts, sensitivity_scores, mat)
    gc()
  } else {
    if (!file.exists(bundle_path)) stop("Missing live replot bundle")
    bundle <- readRDS(bundle_path)
  }
  ####################
  if (mode == "scatter") {
    plot_t4_candidate_scatter(bundle, out, bundle$candidates$Gene, mps)
    index_path <- file.path(out, "reports/figure_index.md")
    if (!file.exists(index_path)) stop("Missing figure index")
    index <- readLines(index_path)
    if (!any(grepl("T4_candidate_per_cell_pemt_scatter.pdf", index, fixed = TRUE))) {
      index <- c(index, "", "## Per-cell scatter supporting evidence", "",
        "- [Per-cell gene-MP scatter plots](../figures/Auto_supporting/T4_candidate_per_cell_pemt_scatter.pdf): 21 pages in candidate-list order; eight MPs per gene, all T4 cells plotted.",
        "- [Scatter page index](../tables/T4_candidate_scatter_page_index.csv).",
        "- [Continuous association table](../tables/T4_candidate_continuous_mp_associations.csv).")
    }
    writeLines(index, index_path)
    status <- "success"
    message("COMPLETE: supporting per-cell scatter PDF")
    return(invisible(NULL))
  }
  ####################
  ####################
  # Refresh Figure 04 alone from the persistent bundle, preserving other figures.
  if (mode == "figure4") {
    plot_t4_candidate_association(bundle, out, bundle$candidates$Gene, mps)
    index_path <- file.path(out, "reports/figure_index.md")
    if (!file.exists(index_path)) stop("Missing current figure index")
    index <- readLines(index_path)
    index <- sub("candidate associations with MP-high cells.",
                 "two pages: high/low expression differences and continuous gene-MP correlations (tested gene omitted).", index, fixed = TRUE)
    writeLines(index, index_path)
    status <- "success"
    message("COMPLETE: two-page Figure 04 and continuous association tables")
    return(invisible(NULL))
  }
  ####################
  ####################
  # Scoped Figure 1 refresh reads only the live source bundle; other results stay intact.
  if (mode == "figure1") {
    summary <- bind_rows(lapply(stages, function(stage) {
      idx <- bundle$meta$sample == stage
      data.frame(stage = stage, gene = bundle$candidates$Gene, n_cells = sum(idx),
                 mean_log1p_cp10k = colMeans(bundle$expression[idx, , drop = FALSE]))
    }))
    plot_candidate_expression_heatmap(summary, out, bundle$candidates$Gene, stages)
    index_path <- file.path(out, "reports/figure_index.md")
    if (!file.exists(index_path)) stop("Missing current figure index")
    index <- readLines(index_path)
    index <- sub("mean log-normalised RNA expression and detection across six stages.",
                 "gene-wise z-score heatmap across six stages; labels show mean log1p(CP10k) expression.", index, fixed = TRUE)
    writeLines(index, index_path)
    status <- "success"
    message("COMPLETE: Figure 1 heatmap refreshed from live source data")
    return(invisible(NULL))
  }
  ####################
  genes <- bundle$candidates$Gene
  meta <- bundle$meta
  scores <- bundle$scores
  expression <- bundle$expression
  detected <- bundle$detected
  write_table(bundle$candidates, "candidate_annotations")
  write_table(bundle$manifest, "input_manifest")
  write_table(bundle$overlap, "candidate_signature_overlap")
  write_table(bundle$metadata_inventory, "metadata_inventory")
  write_table(data.frame(stage = stages, n_cells = as.integer(table(factor(meta$sample, stages)))), "cell_counts")
  expression_summary <- bind_rows(lapply(stages, function(stage) {
    idx <- which(meta$sample == stage)
    data.frame(stage = stage, gene = genes, n_cells = length(idx), mean_log1p_cp10k = colMeans(expression[idx, , drop = FALSE]), median_log1p_cp10k = apply(expression[idx, , drop = FALSE], 2, median), mean_cp10k = colMeans(expm1(expression[idx, , drop = FALSE])), pct_detected = 100 * colMeans(detected[idx, , drop = FALSE]))
  }))
  write_table(expression_summary, "candidate_expression_by_stage")
  score_long <- as.data.frame(scores) |> mutate(cell = rownames(scores), stage = meta$sample) |> pivot_longer(all_of(mps), names_to = "MP", values_to = "score")
  activity_summary <- score_long |> group_by(stage, MP) |> summarise(n_cells = n(), mean_score = mean(score), median_score = median(score), q25 = quantile(score, .25), q75 = quantile(score, .75), .groups = "drop")
  write_table(activity_summary, "mp_activity_by_stage")
  baseline <- meta$sample == "T0"
  thresholds <- bind_rows(lapply(c(.90, .95, .975, .99), function(q) data.frame(MP = mps, reference_quantile = q, threshold = apply(scores[baseline, , drop = FALSE], 2, quantile, probs = q, type = 7), n_reference_cells = sum(baseline))))
  write_table(thresholds, "baseline_thresholds")
  make_flags <- function(x, cutoff) {
    flags <- sweep(x, 2, cutoff, ">")
    cbind(flags, Any_pEMT_MP = rowSums(flags) > 0)
  }
  fractions <- bind_rows(lapply(unique(thresholds$reference_quantile), function(q) {
    cut <- thresholds$threshold[thresholds$reference_quantile == q]
    flags <- make_flags(scores, cut)
    bind_rows(lapply(stages, function(stage) {
      idx <- meta$sample == stage
      data.frame(stage = stage, MP = colnames(flags), reference_quantile = q, n_cells = sum(idx), n_high = colSums(flags[idx, , drop = FALSE]), pct_high = 100 * colMeans(flags[idx, , drop = FALSE]))
    }))
  }))
  write_table(fractions, "high_activity_fractions_threshold_sensitivity")
  write_table(filter(fractions, reference_quantile == .95), "high_activity_fractions_primary")
  write_table(filter(fractions, reference_quantile == .95, stage == "T4"), "T4_high_activity_fractions")
  primary_flags <- make_flags(scores, thresholds$threshold[thresholds$reference_quantile == .95])
  write_table(cbind(meta, setNames(as.data.frame(scores), paste0(mps, "_score")), setNames(as.data.frame(primary_flags), paste0(colnames(primary_flags), "_high"))), "per_cell_scores_and_high_activity")
  t4 <- which(meta$sample == "T4")
  compare <- function(flags, method, tested_genes = genes) {
    if (!identical(rownames(flags), meta$cell[t4])) stop("T4 flag alignment error")
    bind_rows(lapply(colnames(flags), function(mp) {
      high <- flags[, mp]
      nh <- sum(high); nl <- sum(!high)
      bind_rows(lapply(tested_genes, function(gene) {
        x <- expression[t4[high], gene]; y <- expression[t4[!high], gene]
        if (nh < 3L || nl < 3L) return(data.frame(method = method, MP = mp, gene = gene, n_high = nh, n_low = nl, status = "insufficient_cells", mean_high = NA_real_, mean_low = NA_real_, delta_mean_log1p = NA_real_, log2FC_cp10k_plus1 = NA_real_, pct_detected_high = NA_real_, pct_detected_low = NA_real_, delta_pct_detected = NA_real_, auc = NA_real_, p_cell_exploratory = NA_real_))
        auc <- (sum(rank(c(x, y), ties.method = "average")[seq_len(nh)]) - nh * (nh + 1) / 2) / (nh * nl)
        p <- if (length(unique(c(x,y))) < 2L) 1 else suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
        data.frame(method = method, MP = mp, gene = gene, n_high = nh, n_low = nl, status = "ok", mean_high = mean(x), mean_low = mean(y), delta_mean_log1p = mean(x) - mean(y), log2FC_cp10k_plus1 = log2((mean(expm1(x)) + 1) / (mean(expm1(y)) + 1)), pct_detected_high = 100 * mean(x > 0), pct_detected_low = 100 * mean(y > 0), delta_pct_detected = 100 * (mean(x > 0) - mean(y > 0)), auc = auc, p_cell_exploratory = p)
      }))
    }))
  }
  comparisons <- compare(primary_flags[t4, , drop = FALSE], "cached_original")
  rescored <- bundle$rescored
  baseline_ids <- meta$cell[baseline]
  t4_ids <- meta$cell[t4]
  baseline_sc <- rescored[baseline_ids, paste0(mps, "_full"), drop = FALSE]
  full_sc <- rescored[t4_ids, paste0(mps, "_full"), drop = FALSE]
  colnames(full_sc) <- mps
  full_cut <- apply(baseline_sc, 2, quantile, .95, type = 7)
  full_flags <- make_flags(full_sc, full_cut)
  comparisons <- bind_rows(comparisons, compare(full_flags, "rescored_full"))
  agreement <- list(data.frame(method = "full", gene = NA_character_, MP = colnames(full_flags), pct_high = 100 * colMeans(full_flags), agreement_with_cached_pct = 100 * colMeans(full_flags == primary_flags[t4, , drop = FALSE]), agreement_with_rescored_full_pct = 100))
  write_table(data.frame(MP = mps, threshold = full_cut, reference_quantile = .95), "baseline_thresholds_full")
  loo_thresholds <- list()
  loo_flags <- list()
  for (gene in genes) {
    sc <- rescored[, paste0(mps, "_loo_", gene), drop = FALSE]
    colnames(sc) <- mps
    cut <- apply(sc[baseline_ids, , drop = FALSE], 2, quantile, .95, type = 7)
    flags <- make_flags(sc[t4_ids, , drop = FALSE], cut)
    comparisons <- bind_rows(comparisons, compare(flags, "rescored_loo", tested_genes = gene))
    loo_thresholds[[gene]] <- data.frame(gene = gene, MP = mps, threshold = cut, reference_quantile = .95)
    loo_flags[[gene]] <- data.frame(gene = gene, cell = t4_ids, as.data.frame(flags))
    agreement[[gene]] <- data.frame(method = "leave_one_gene_out", gene = gene, MP = colnames(flags), pct_high = 100 * colMeans(flags), agreement_with_cached_pct = 100 * colMeans(flags == primary_flags[t4, , drop = FALSE]), agreement_with_rescored_full_pct = 100 * colMeans(flags == full_flags))
  }
  write_table(bind_rows(loo_thresholds), "baseline_thresholds_leave_one_gene_out")
  write_table(bind_rows(loo_flags), "T4_flags_leave_one_gene_out")
  comparisons <- comparisons |> group_by(method) |> mutate(p_adj_cell_exploratory = p.adjust(p_cell_exploratory, "BH")) |> ungroup()
  write_table(comparisons, "T4_candidate_high_vs_low_all_methods")
  write_table(filter(comparisons, method == "cached_original"), "T4_candidate_high_vs_low_primary")
  write_table(bind_rows(agreement), "T4_classification_sensitivity")
  version_check <- bind_rows(lapply(mps, function(mp) data.frame(MP = mp, n_cells = nrow(rescored), spearman_cached_vs_rescored = cor(scores[rownames(rescored), mp], rescored[, paste0(mp, "_full")], method = "spearman"), max_absolute_difference = max(abs(scores[rownames(rescored), mp] - rescored[, paste0(mp, "_full")])))))
  write_table(version_check, "cached_vs_rescored_score_check")
  ranking <- filter(comparisons, MP == "Any_pEMT_MP") |> select(method, gene, delta_mean_log1p, log2FC_cp10k_plus1, delta_pct_detected, auc) |> pivot_wider(names_from = method, values_from = c(delta_mean_log1p, log2FC_cp10k_plus1, delta_pct_detected, auc)) |> mutate(positive_in_both = delta_mean_log1p_cached_original > 0 & delta_mean_log1p_rescored_loo > 0) |> arrange(desc(positive_in_both), desc(auc_rescored_loo))
  ranking <- left_join(ranking, bundle$candidates, by = c("gene" = "Gene"))
  write_table(ranking, "candidate_prioritisation_any_pemt")
  qc <- bind_rows(lapply(colnames(primary_flags), function(mp) bind_rows(lapply(c(FALSE, TRUE), function(h) {
    idx <- t4[primary_flags[t4, mp] == h]
    data.frame(MP = mp, high = h, n_cells = length(idx), median_UMI = median(meta$nCount_RNA_measured[idx]), median_detected_genes = median(meta$nFeature_RNA_measured[idx]))
  }))))
  write_table(qc, "T4_high_low_sequencing_depth")

  # Exploratory quantitative figures: expression overview, activity distributions,
  # baseline-relative abundance, and within-T4 expression relationships.
  theme_set(theme_classic(base_size = 13, base_family = "sans") + theme(axis.text = element_text(colour = "black"), strip.background = element_rect(fill = "grey95"), plot.title = element_text(face = "bold"), legend.position = "right"))
  save_plot <- function(p, name, width, height, supporting = FALSE) {
    path <- file.path(out, "figures", if (supporting) "Auto_supporting" else "Auto_final", name)
    ggsave(paste0(path, ".pdf"), p, device = grDevices::cairo_pdf, width = width, height = height, limitsize = FALSE)
    # Raster companion uses the same R plot for direct review and slide placement.
    ggsave(paste0(path, ".png"), p, device = "png", type = "cairo", width = width, height = height, dpi = 160, limitsize = FALSE)
  }
  expression_summary$stage <- factor(expression_summary$stage, stages)
  expression_summary$gene <- factor(expression_summary$gene, rev(genes))
  plot_candidate_expression_heatmap(expression_summary, out, genes, stages)
  score_long$stage <- factor(score_long$stage, stages)
  score_long$MP <- factor(score_long$MP, mps)
  cut_plot <- filter(thresholds, reference_quantile == .95)
  p <- ggplot(score_long, aes(stage, score)) + geom_violin(fill = "#c8dbe7", colour = "#406780", scale = "width", linewidth = .3) + geom_boxplot(width = .13, outlier.shape = NA, fill = "white", linewidth = .3) + geom_hline(data = cut_plot, aes(yintercept = threshold), colour = "#a64836", linetype = 2) + facet_wrap(~MP, ncol = 4, scales = "free_y") + labs(x = NULL, y = "UCell score", title = "Eight pEMT/plasticity metaprograms", subtitle = "Dashed line: each MP's T0 95th percentile; box: median and IQR, whiskers: 1.5 IQR")
  save_plot(p, "02_pemt_activity_distributions", 16, 8)
  primary_frac <- filter(fractions, reference_quantile == .95)
  primary_frac$stage <- factor(primary_frac$stage, stages)
  primary_frac$MP <- factor(primary_frac$MP, c(mps, "Any_pEMT_MP"))
  p <- ggplot(primary_frac, aes(stage, pct_high, group = 1)) + geom_line(colour = "#406780") + geom_point(size = 2.5, colour = "#406780") + geom_text(aes(label = sprintf("%.1f", pct_high)), vjust = -.9, size = 3.3) + facet_wrap(~MP, ncol = 3) + scale_y_continuous(limits = c(0,108), breaks = seq(0,100,25)) + labs(x = NULL, y = "High-activity cells (%)", title = "Baseline-relative pEMT activity across stages", subtitle = "High = score > T0 95th percentile; Any_pEMT_MP = high in at least one MP (overlapping groups)")
  save_plot(p, "03_high_activity_cell_fractions", 14, 10)
  ####################
  plot_t4_candidate_association(bundle, out, genes, mps)
  plot_t4_candidate_scatter(bundle, out, genes, mps)
  ####################
  trend <- expression_summary
  trend$gene <- factor(as.character(trend$gene), genes)
  p <- ggplot(trend, aes(stage, mean_log1p_cp10k, group = 1)) + geom_line(colour = "#406780") + geom_point(colour = "#406780", size = 2) + facet_wrap(~gene, ncol = 4, scales = "free_y") + labs(x = NULL, y = "Mean log1p(CP10k)", title = "Candidate expression trajectories", subtitle = "Descriptive means; connected stages do not imply matched individual cells")
  save_plot(p, "candidate_expression_trajectories", 15, 15, TRUE)
  sens <- filter(comparisons, MP == "Any_pEMT_MP")
  sens$gene <- factor(sens$gene, rev(genes))
  p <- ggplot(sens, aes(delta_mean_log1p, gene, colour = method)) + geom_vline(xintercept = 0, colour = "grey75") + geom_point(position = position_dodge(width = .55), size = 2.5) + scale_colour_manual(values = c(cached_original = "#406780", rescored_full = "#777777", rescored_loo = "#b36b3a")) + labs(x = "High minus low mean log1p(CP10k)", y = NULL, title = "T4 association sensitivity to candidate inclusion", subtitle = "Any of eight MPs high; gene-specific sensitivity thresholds recalibrated from T0", colour = "Scores")
  save_plot(p, "leave_one_gene_out_sensitivity", 13, 10, TRUE)
  fs <- filter(fractions, stage == "T4")
  fs$MP <- factor(fs$MP, c(mps, "Any_pEMT_MP"))
  p <- ggplot(fs, aes(reference_quantile * 100, pct_high, colour = MP)) + geom_line() + geom_point() + scale_x_continuous(breaks = c(90,95,97.5,99)) + labs(x = "T0 reference percentile", y = "T4 high-activity cells (%)", title = "T4 prevalence depends on the operational cutoff")
  save_plot(p, "T4_threshold_sensitivity", 12, 7, TRUE)
  t4_fractions <- filter(fractions, reference_quantile == .95, stage == "T4")
  top <- head(ranking$gene[ranking$positive_in_both %in% TRUE], 8)
  report <- c("# Candidate expression and pEMT activity", "", "This report describes exploratory associations in the six final Parse sample objects. It is not independent functional validation.", "", "## Definitions", "", "The eight MPs are MP48, MP50, MP60, MP109, MP16, MP74, MP83 and MP9. Candidate genes and recurrence/rank annotations are read directly from the supplied workbook.", "", "Primary high activity means a cached UCell score strictly greater than that MP's T0 95th percentile (R quantile type 7). The same numerical threshold is applied at every stage. This is baseline-relative elevation, not a universal biological on/off state. T0 is a reference distribution, not a known pEMT-negative control. Low means at or below the threshold, not proven inactive. The 90th, 97.5th and 99th percentiles quantify cutoff sensitivity.", "", "Each MP is analysed separately. Any_pEMT_MP is the union of high cells, counted once per cell. Its baseline percentage can exceed 5% because eight correlated MPs are screened; this is not a family-wise 5% classification rule. Individual MP percentages must not be added.", "", "## T4 percentages", "", paste0("- ", t4_fractions$MP, ": ", sprintf("%.2f", t4_fractions$pct_high), "% (", t4_fractions$n_high, "/", t4_fractions$n_cells, ")."), "", "## Expression comparisons", "", "RNA expression is log1p(raw gene count / total RNA counts * 10,000); detection means a nonzero raw count. Tables include all-cell means, medians, detection percentages, high-minus-low expression and detection differences, log2 ratios of mean CP10k plus a pseudocount of 1, and rank AUC (0.5 means no directional separation).", "", "The original signatures are the primary analysis. Leave-one-gene-out sensitivity removes only the candidate being tested from each MP signature, not from the expression matrix or rank universe; rescoring is restricted to T0 and T4. This preserves short signatures: removing all candidates would leave MP83 with only four genes. Each candidate therefore has its own sensitivity classification, while primary high-cell prevalence always uses the original signatures. Full signatures are rescored in the same call as reduced signatures to distinguish UCell-version effects from candidate-removal effects. T0 thresholds are recalculated for each sensitivity signature. Leave-one-gene-out scores still contain other candidates and correlated genes and therefore do not provide independent validation.", "", paste0("Candidates with positive mean-expression differences in both the original and leave-one-gene-out Any_pEMT_MP comparisons, ordered by leave-one-gene-out AUC: ", paste(top, collapse = ", "), ". This is an exploratory ranking, not a significance-based selection."), "", "## Limits and interpretation", "", "Cell-level Wilcoxon p-values are exploratory only; BH adjustment spans all 21 x 9 contrasts within each scoring method. Cells are not biological replicates. No donor- or replicate-level treatment inference is claimed. The metadata inventory is provided; sequencing-depth summaries show potential high/low technical confounding. Raw per-cell source data remain in the live bundle for further analysis.", "", "MP discovery, T2/T4 enrichment selection and candidate choice used this dataset, so observed treatment-associated patterns are subject to selection bias. Gene expression in MP-high cells establishes association, not causality, druggability, protein abundance or functional necessity. Correlated MP signatures produce overlapping high groups.", "", "## Method sources", "", "- [UCell original paper](https://doi.org/10.1016/j.csbj.2021.06.043): rank-based per-cell signature scoring.", "- [UCell author guidance](https://github.com/carmonalab/UCell/issues/7): score interpretation depends on signature/data quality; no universal high-activity cutoff. The T0 q95 choice here is an operational analysis decision, not a cutoff prescribed by UCell.", "- [UCell implementation](https://github.com/carmonalab/UCell): each signature gene contributes to the rank sum; rank scoring limits but does not eliminate candidate self-contribution.", "", paste0("Sensitivity UCell version: ", bundle$UCell_version, "; common rank-universe genes: ", bundle$common_gene_count, "."))
  writeLines(report, file.path(out, "reports/current_results.md"))
  index <- c("# Current figure and table index", "", "## Selected final figures", "", "1. [Candidate expression](../figures/Auto_final/01_candidate_expression_by_stage.pdf): gene-wise z-score heatmap across six stages; labels show mean log1p(CP10k) expression.", "2. [MP activity distributions](../figures/Auto_final/02_pemt_activity_distributions.pdf): all cell scores, with baseline cutoffs.", "3. [High-cell percentages](../figures/Auto_final/03_high_activity_cell_fractions.pdf): per-MP and union percentages across stages.", "4. [T4 expression differences](../figures/Auto_final/04_T4_candidate_high_vs_low.pdf): two pages: high/low expression differences and continuous gene-MP correlations (tested gene omitted).", "", "## Supporting figures", "", "- [Gene trajectories](../figures/Auto_supporting/candidate_expression_trajectories.pdf).", "- [Leave-one-gene-out sensitivity](../figures/Auto_supporting/leave_one_gene_out_sensitivity.pdf).", "- [Threshold sensitivity](../figures/Auto_supporting/T4_threshold_sensitivity.pdf).", "", "## Key tables and source data", "", "- [T4 percentages](../tables/T4_high_activity_fractions.csv).", "- [Candidate prioritisation](../tables/candidate_prioritisation_any_pemt.csv).", "- [All high/low comparisons](../tables/T4_candidate_high_vs_low_all_methods.csv).", "- [Expression across stages](../tables/candidate_expression_by_stage.csv).", "- [Per-cell scores and classifications](../tables/per_cell_scores_and_high_activity.csv).", "- [Persistent replot bundle](../intermediate/candidate_pemt_replot_bundle.rds).", "- [Methods, interpretation and limitations](current_results.md).")
  index <- c(index, "", "## Per-cell scatter supporting evidence", "",
    "- [Per-cell gene-MP scatter plots](../figures/Auto_supporting/T4_candidate_per_cell_pemt_scatter.pdf): 21 pages in candidate-list order; eight MPs per gene, all T4 cells plotted.",
    "- [Scatter page index](../tables/T4_candidate_scatter_page_index.csv).",
    "- [Continuous association table](../tables/T4_candidate_continuous_mp_associations.csv).")
  writeLines(index, file.path(out, "reports/figure_index.md"))
  validation <- c(paste0("Cells: ", nrow(meta)), paste0("Candidates: ", length(genes)), paste0("MPs: ", length(mps)), "Unique cells and complete score/expression alignment: PASS", "Six input samples exactly match cached score cell populations: PASS", "All expression and score values finite: PASS", paste0("T4 primary contrasts with sufficient cells: ", sum(comparisons$status[comparisons$method == "cached_original"] == "ok"), "/189"), paste0("Live replot bundle bytes: ", file.info(bundle_path)$size))
  writeLines(validation, file.path(out, "reports/validation.txt"))
  status <- "success"
  message("COMPLETE: ", out)
}
main()
####################
