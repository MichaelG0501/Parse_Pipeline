####################
# parse_highres_mp_tcga_survival_volcano.R
#
# Description:
#   Terminal TCGA survival-association workflow for high-resolution Parse MPs.
#   Uses the same TCGA ESCA metadata and whole-profile GSVA reference mode as
#   scRef_Pipeline/analysis/clinical/tcga_mp_state_survival_reg_noreg.R mode
#   (b): GSVA on whole TCGA profile. Produces one volcano plot each for strict
#   increasing MPs, strict decreasing MPs, and legacy T2/T4-high MPs.
#
# Inputs:
#   parse_outs/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_selected_mp_genes_nMP117.rds
#   parse_outs/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_trend_summary_nMP117.csv
#   parse_outs/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_top_3CA_noncellcycle_nMP117.csv
#   parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_selected_mp_genes_nMP117.rds
#   parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_top_3CA_noncellcycle_nMP117.csv
#   /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/tcga_esca_meta.rds
#   /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#
# Outputs:
#   parse_outs/highres_mp_tcga_survival/tables/Auto_parse_highres_mp_tcga_survival_cox_whole_tcga.csv
#   parse_outs/highres_mp_tcga_survival/tables/Auto_parse_highres_mp_tcga_survival_gene_set_summary.csv
#   parse_outs/highres_mp_tcga_survival/figures/Auto_parse_highres_mp_tcga_survival_volcano_whole_tcga.pdf
#   parse_outs/highres_mp_tcga_survival/figures/Auto_parse_highres_mp_tcga_survival_volcano_<group>.pdf
#   parse_outs/highres_mp_tcga_survival/figures/Auto_parse_highres_mp_tcga_survival_volcano_<group>_<split>.pdf
#   parse_outs/highres_mp_tcga_survival/intermediate/Auto_parse_highres_mp_tcga_gsva_scores_whole_tcga.rds
#   parse_outs/highres_mp_tcga_survival/reports/Auto_parse_highres_mp_tcga_survival_summary.txt
#   parse_outs/logs/run_summaries/parse_highres_mp_tcga_survival_volcano_*.txt
#
# Cache / replot:
#   GSVA scores are rebuilt by default because the gene lists are small. Set
#   PARSE_TCGA_REUSE_GSVA=TRUE to reuse the cached score RDS when only plotting
#   or table formatting changes.
#
# Methodology:
#   analysis/methodology/metaprograms/highres_mp_tcga_survival_volcano_methodology.md
#
# Downstream status:
#   Terminal clinical-association figure workflow; no active downstream
#   dependency should consume its outputs without explicit review.
####################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(GSVA)
library(survival)
library(tidyr)

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

####################
# Paths and run setup
####################
project_dir <- "/rds/general/project/spatialtranscriptomics/live/Parse_Pipeline"
paths <- parse_paths(project_dir)
out_dir <- file.path(project_dir, "parse_outs", "centred", "highres_mp_tcga_survival")
tiers <- parse_output_tiers(out_dir, create = TRUE)

base_highres_dir <- file.path(project_dir, "parse_outs", "centred", "Auto_parse_highres_metaprogram_trends")

strict_gene_path <- file.path(
  base_highres_dir,
  "Auto_parse_highres_selected_mp_genes_nMP117.rds"
)
strict_trend_path <- file.path(
  base_highres_dir,
  "Auto_parse_highres_trend_summary_nMP117.csv"
)
t2t4_gene_path <- file.path(
  base_highres_dir,
  "Auto_T2T4_gt_T0eR4_filter",
  "Auto_parse_highres_T2T4_selected_mp_genes_nMP117.rds"
)
t2t4_trend_path <- file.path(
  base_highres_dir,
  "Auto_T2T4_gt_T0eR4_filter",
  "Auto_parse_highres_T2T4_filter_summary_nMP117.csv"
)
strict_3ca_label_path <- file.path(
  base_highres_dir,
  "Auto_parse_highres_top_3CA_noncellcycle_nMP117.csv"
)
t2t4_3ca_label_path <- file.path(
  base_highres_dir,
  "Auto_T2T4_gt_T0eR4_filter",
  "Auto_parse_highres_T2T4_top_3CA_noncellcycle_nMP117.csv"
)

tcga_meta_candidates <- c(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/tcga_esca_meta.rds",
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds"
)
tcga_tpm_candidates <- c(
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt",
  "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/TCGA/esca_gdc_reconstruction/tables/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt"
)

tcga_meta_path <- tcga_meta_candidates[file.exists(tcga_meta_candidates)][1]
tcga_tpm_path <- tcga_tpm_candidates[file.exists(tcga_tpm_candidates)][1]
if (is.na(tcga_meta_path)) {
  stop("Could not find TCGA metadata in expected scRef locations.")
}
if (is.na(tcga_tpm_path)) {
  stop("Could not find TCGA whole-profile TPM matrix in expected scRef locations.")
}

gsva_cache_path <- file.path(tiers$intermediate, "Auto_parse_highres_mp_tcga_gsva_scores_whole_tcga.rds")
reuse_gsva <- tolower(Sys.getenv("PARSE_TCGA_REUSE_GSVA", unset = "FALSE")) %in% c("true", "1", "yes", "y")
split_methods <- c("continuous", "median", "q1q4")

script_run <- parse_start_run(
  "parse_highres_mp_tcga_survival_volcano",
  parameters = list(
    tcga_mode = "whole_tcga_reference_gsva",
    survival_model = "Cox, EAC primary tumour only",
    split_methods = paste(split_methods, collapse = ","),
    min_gsva_genes = 5,
    reuse_gsva = reuse_gsva
  ),
  input_files = c(strict_gene_path, strict_trend_path, t2t4_trend_path, strict_3ca_label_path, t2t4_gene_path, t2t4_3ca_label_path, tcga_meta_path, tcga_tpm_path),
  output_files = c(
    file.path(tiers$tables, "Auto_parse_highres_mp_tcga_survival_cox_whole_tcga.csv"),
    file.path(tiers$tables, "Auto_parse_highres_mp_tcga_survival_gene_set_summary.csv"),
    file.path(tiers$figures, "Auto_parse_highres_mp_tcga_survival_volcano_whole_tcga.pdf"),
    gsva_cache_path
  ),
  reused_cached = reuse_gsva
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

####################
# Helpers
####################
make_tcga_page <- function(plot1, plot2, plot3, page_title) {
  gridExtra::arrangeGrob(
    plot1, plot2, plot3,
    ncol = 3,
    top = grid::textGrob(page_title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required input: ", path)
  }
  invisible(path)
}

infer_histology <- function(type_vec, detailed_vec = NA_character_) {
  joined <- tolower(paste(as.character(type_vec), as.character(detailed_vec)))
  out <- rep("Other", length(joined))
  out[grepl("adeno", joined)] <- "EAC"
  out[grepl("squamous", joined)] <- "ESCC"
  out
}

clean_gene_sets <- function(gene_sets) {
  gene_sets <- lapply(gene_sets, function(genes) {
    genes <- unique(as.character(genes))
    genes[!is.na(genes) & nzchar(genes)]
  })
  gene_sets[lengths(gene_sets) > 0]
}

make_group_gene_sets <- function(strict_genes, strict_trends, t2t4_genes, t2t4_trends) {
  t2t4_mps <- t2t4_trends$MP[t2t4_trends$retained == TRUE & t2t4_trends$MP %in% names(t2t4_genes)]
  strict_mps <- strict_trends$MP[strict_trends$retained == TRUE & strict_trends$trend_type_label == "Decrease, Consistent" & strict_trends$MP %in% names(strict_genes)]
  
  strict_mps <- setdiff(strict_mps, t2t4_mps)

  list(
    decrease_consistent = strict_genes[strict_mps],
    t2t4_high = t2t4_genes[t2t4_mps]
  ) |>
    lapply(clean_gene_sets)
}

make_unique_feature_names <- function(group_gene_sets) {
  lapply(names(group_gene_sets), function(group_name) {
    gene_sets <- group_gene_sets[[group_name]]
    names(gene_sets) <- paste(group_name, names(gene_sets), sep = "__")
    gene_sets
  }) |>
    unlist(recursive = FALSE, use.names = TRUE)
}

read_tcga_tpm <- function(path) {
  tpm_df <- data.table::fread(path)
  if (!"GeneSymbol" %in% colnames(tpm_df)) {
    stop("TCGA TPM matrix must contain a GeneSymbol column: ", path)
  }
  gene_symbols <- as.character(tpm_df$GeneSymbol)
  expr_mat <- as.matrix(tpm_df[, setdiff(colnames(tpm_df), "GeneSymbol"), with = FALSE])
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- gene_symbols
  expr_mat[is.na(expr_mat)] <- 0

  if (anyDuplicated(rownames(expr_mat)) > 0) {
    expr_mat <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE)
    gene_counts <- as.numeric(table(gene_symbols)[rownames(expr_mat)])
    expr_mat <- expr_mat / gene_counts
  }

  expr_mat
}

run_gsva_scores <- function(expr_mat, gene_sets, min_genes = 5) {
  filtered_sets <- lapply(gene_sets, function(genes) intersect(unique(genes), rownames(expr_mat)))
  filtered_sets <- filtered_sets[lengths(filtered_sets) >= min_genes]
  if (length(filtered_sets) == 0) {
    stop("No gene sets retained for GSVA after TCGA gene overlap filtering.")
  }

  if (exists("gsvaParam", where = asNamespace("GSVA"), mode = "function")) {
    param <- GSVA::gsvaParam(expr_mat, filtered_sets, kcdf = "Gaussian")
    GSVA::gsva(param)
  } else {
    GSVA::gsva(expr_mat, filtered_sets, method = "gsva", kcdf = "Gaussian")
  }
}

run_cox <- function(df, feature_cols, split_method = "continuous") {
  out <- vector("list", length(feature_cols))
  names(out) <- feature_cols
  for (feat in feature_cols) {
    d <- df |>
      dplyr::filter(
        !is.na(OS_time),
        !is.na(OS_event),
        is.finite(.data[[feat]]),
        HistologyGroup == "EAC"
      )
    if (nrow(d) < 20 || stats::var(d[[feat]], na.rm = TRUE) == 0 || sum(d$OS_event == 1, na.rm = TRUE) < 3) {
      next
    }

    if (split_method == "median") {
      med_val <- stats::median(d[[feat]], na.rm = TRUE)
      d$split_val <- factor(ifelse(d[[feat]] > med_val, "High", "Low"), levels = c("Low", "High"))
      form <- survival::Surv(OS_time, OS_event) ~ split_val
    } else if (split_method == "q1q4") {
      quants <- stats::quantile(d[[feat]], probs = c(0.25, 0.75), na.rm = TRUE)
      d <- d |>
        dplyr::filter(.data[[feat]] <= quants[1] | .data[[feat]] >= quants[2])
      if (nrow(d) < 20 || sum(d$OS_event == 1, na.rm = TRUE) < 3) {
        next
      }
      d$split_val <- factor(ifelse(d[[feat]] >= quants[2], "High", "Low"), levels = c("Low", "High"))
      form <- survival::Surv(OS_time, OS_event) ~ split_val
    } else {
      form <- stats::as.formula(paste0("survival::Surv(OS_time, OS_event) ~ `", feat, "`"))
    }

    fit <- try(survival::coxph(form, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) {
      next
    }
    ss <- summary(fit)
    out[[feat]] <- data.frame(
      feature_id = feat,
      split_method = split_method,
      HR = ss$coefficients[1, "exp(coef)"],
      P_value = ss$coefficients[1, "Pr(>|z|)"],
      n = fit$n,
      events = fit$nevent,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(out)
}

format_group_label <- function(x) {
  dplyr::case_when(
    x == "decrease_consistent" ~ "Decrease/Consistent MPs",
    x == "t2t4_high" ~ "T2/T4-high MPs",
    TRUE ~ x
  )
}

format_split_label <- function(x) {
  dplyr::case_when(
    x == "continuous" ~ "continuous Cox",
    x == "median" ~ "median high vs low Cox",
    x == "q1q4" ~ "upper vs lower quartile Cox",
    TRUE ~ x
  )
}

read_3ca_label_table <- function(path) {
  if (!file.exists(path)) {
    return(data.frame(MP = character(), top_3ca_noncc = character(), stringsAsFactors = FALSE))
  }
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) |>
    dplyr::mutate(
      top_3ca_noncc = ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc)
    )
}

make_label_lookup <- function(strict_label_path, t2t4_label_path) {
  strict_labels <- read_3ca_label_table(strict_label_path) |> dplyr::mutate(mp_group = "decrease_consistent")
  t2t4_labels <- read_3ca_label_table(t2t4_label_path) |> dplyr::mutate(mp_group = "t2t4_high")

  dplyr::bind_rows(strict_labels, t2t4_labels) |>
    dplyr::mutate(display_label = paste0(MP, "\n", top_3ca_noncc)) |>
    dplyr::select(mp_group, MP, top_3ca_noncc, display_label)
}

plot_volcano <- function(df, ttl) {
  if (nrow(df) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No Cox models available", size = 5) +
        ggplot2::labs(title = ttl)
    )
  }

  pdat <- df |>
    dplyr::mutate(
      sig = P_value < 0.05,
      log2HR = log2(HR),
      neglog10 = -log10(P_value),
      feature_label = display_label
    )

  ggplot2::ggplot(pdat, ggplot2::aes(log2HR, neglog10)) +
    ggplot2::geom_point(ggplot2::aes(color = sig), size = 2.8, alpha = 0.9) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "grey45") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey45") +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = feature_label),
      size = 2.8,
      max.overlaps = 100
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "grey70", "TRUE" = "firebrick3"),
      guide = "none"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = ttl,
      x = "log2(HR)",
      y = "-log10(p)"
    )
}

write_plot_outputs <- function(plot_obj, group_name, split_method = "continuous") {
  split_suffix <- if (split_method == "continuous") "" else paste0("_", split_method)
  stem <- paste0("Auto_parse_highres_mp_tcga_survival_volcano_", group_name, split_suffix)
  ggplot2::ggsave(
    filename = file.path(tiers$figures, paste0(stem, ".pdf")),
    plot = plot_obj,
    width = 9,
    height = 7,
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    filename = file.path(tiers$figures, paste0(stem, ".png")),
    plot = plot_obj,
    width = 9,
    height = 7,
    dpi = parse_plot_defaults$dpi
  )
}

####################
# Load inputs
####################
invisible(lapply(c(strict_gene_path, strict_trend_path, t2t4_gene_path, t2t4_trend_path), stop_if_missing))

strict_genes <- readRDS(strict_gene_path)
strict_trends <- read.csv(strict_trend_path, check.names = FALSE, stringsAsFactors = FALSE)
t2t4_genes <- readRDS(t2t4_gene_path)
t2t4_trends <- read.csv(t2t4_trend_path, check.names = FALSE, stringsAsFactors = FALSE)
group_gene_sets <- make_group_gene_sets(strict_genes, strict_trends, t2t4_genes, t2t4_trends)
all_gene_sets <- make_unique_feature_names(group_gene_sets)
label_lookup <- make_label_lookup(strict_3ca_label_path, t2t4_3ca_label_path)

if (length(all_gene_sets) == 0) {
  stop("No MP gene sets were available from high-resolution Parse outputs.")
}

tcga_meta <- readRDS(tcga_meta_path)
if (!"HistologyGroup" %in% colnames(tcga_meta)) {
  tcga_meta$HistologyGroup <- infer_histology(tcga_meta$type %||% NA_character_, tcga_meta$Cancer_Type_Detailed %||% NA_character_)
}
tcga_meta <- tcga_meta |>
  dplyr::filter(sample_type_code == "01")

tcga_tpm <- read_tcga_tpm(tcga_tpm_path)

####################
# GSVA and Cox models
####################
if (reuse_gsva && file.exists(gsva_cache_path)) {
  gsva_obj <- readRDS(gsva_cache_path)
  gsva_scores <- gsva_obj$scores
  filtered_gene_sets <- gsva_obj$filtered_gene_sets
} else {
  filtered_gene_sets <- lapply(all_gene_sets, function(genes) intersect(unique(genes), rownames(tcga_tpm)))
  filtered_gene_sets <- filtered_gene_sets[lengths(filtered_gene_sets) >= 5]
  gsva_scores <- run_gsva_scores(tcga_tpm, all_gene_sets, min_genes = 5)
  saveRDS(
    list(
      scores = gsva_scores,
      filtered_gene_sets = filtered_gene_sets,
      tcga_meta_path = tcga_meta_path,
      tcga_tpm_path = tcga_tpm_path
    ),
    gsva_cache_path
  )
}

score_df <- as.data.frame(t(gsva_scores))
score_df$sample_barcode <- rownames(score_df)
merged_df <- tcga_meta |>
  dplyr::left_join(score_df, by = "sample_barcode")

cox_res <- dplyr::bind_rows(lapply(split_methods, function(split_method) {
  run_cox(merged_df, rownames(gsva_scores), split_method = split_method)
})) |>
  tidyr::separate(feature_id, into = c("mp_group", "MP"), sep = "__", remove = FALSE) |>
  dplyr::left_join(label_lookup, by = c("mp_group", "MP")) |>
  dplyr::mutate(
    method = "whole_tcga",
    mode = "reference_gsva",
    cohort = "EAC",
    feature_type = "highres_MP",
    mp_group_label = format_group_label(mp_group),
    top_3ca_noncc = ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc),
    display_label = ifelse(is.na(display_label) | display_label == "", paste0(MP, "\n", top_3ca_noncc), display_label),
    log2HR = log2(HR),
    neglog10_p = -log10(P_value)
  ) |>
  dplyr::select(
    mp_group,
    mp_group_label,
    MP,
    method,
    mode,
    cohort,
    feature_type,
    split_method,
    HR,
    log2HR,
    P_value,
    neglog10_p,
    n,
    events,
    top_3ca_noncc,
    display_label,
    feature_id
  )

gene_set_summary <- data.frame(
  feature_id = names(all_gene_sets),
  source_gene_count = lengths(all_gene_sets),
  tcga_overlap_gene_count = lengths(lapply(all_gene_sets, function(genes) intersect(unique(genes), rownames(tcga_tpm)))),
  retained_for_gsva = names(all_gene_sets) %in% rownames(gsva_scores),
  stringsAsFactors = FALSE
) |>
  tidyr::separate(feature_id, into = c("mp_group", "MP"), sep = "__", remove = FALSE) |>
  dplyr::left_join(label_lookup, by = c("mp_group", "MP")) |>
  dplyr::mutate(mp_group_label = format_group_label(mp_group)) |>
  dplyr::select(mp_group, mp_group_label, MP, top_3ca_noncc, display_label, source_gene_count, tcga_overlap_gene_count, retained_for_gsva, feature_id)

write.csv(
  cox_res,
  file.path(tiers$tables, "Auto_parse_highres_mp_tcga_survival_cox_whole_tcga.csv"),
  row.names = FALSE
)
write.csv(
  gene_set_summary,
  file.path(tiers$tables, "Auto_parse_highres_mp_tcga_survival_gene_set_summary.csv"),
  row.names = FALSE
)

####################
# Volcano figures and report
####################
group_order <- c("t2t4_high", "decrease_consistent")
plot_list <- list()
for (split_method in split_methods) {
  for (group_name in group_order) {
    plot_key <- paste(group_name, split_method, sep = "|")
    plot_list[[plot_key]] <- plot_volcano(
      cox_res |> dplyr::filter(mp_group == group_name, split_method == .env$split_method),
      paste0("whole_tcga MP volcano (", split_method, ")")
    )
    write_plot_outputs(plot_list[[plot_key]], group_name, split_method)
  }
}

volcano_page_t2t4 <- make_tcga_page(
  plot_list[["t2t4_high|continuous"]],
  plot_list[["t2t4_high|median"]],
  plot_list[["t2t4_high|q1q4"]],
  "Centred MP volcano: T2/T4-high MPs"
)

volcano_page_decrease <- make_tcga_page(
  plot_list[["decrease_consistent|continuous"]],
  plot_list[["decrease_consistent|median"]],
  plot_list[["decrease_consistent|q1q4"]],
  "Centred MP volcano: Decrease/Consistent MPs"
)

grDevices::cairo_pdf(
  filename = file.path(tiers$figures, "Auto_parse_highres_mp_tcga_survival_volcano_whole_tcga.pdf"),
  width = 18,
  height = 8,
  onefile = TRUE
)
grid::grid.newpage()
grid::grid.draw(volcano_page_t2t4)
grid::grid.newpage()
grid::grid.draw(volcano_page_decrease)
dev.off()

summary_grid <- expand.grid(
  mp_group = group_order,
  split_method = split_methods,
  stringsAsFactors = FALSE
)
summary_lines <- c(
  "High-resolution Parse MP TCGA survival volcano workflow",
  paste0("TCGA metadata: ", tcga_meta_path),
  paste0("TCGA whole-profile TPM: ", tcga_tpm_path),
  paste0("Primary tumour EAC samples after metadata filter: ", nrow(tcga_meta |> dplyr::filter(HistologyGroup == "EAC"))),
  "",
  "Gene sets by group:",
  paste0("  - ", format_group_label(names(group_gene_sets)), ": ", lengths(group_gene_sets), " source MPs"),
  "",
  "Nominal Cox results by group and split method:",
  paste0(
    "  - ",
    format_group_label(summary_grid$mp_group),
    " / ",
    summary_grid$split_method,
    ": nominal p<0.05 = ",
    mapply(function(g, sm) {
      sum(cox_res$mp_group == g & cox_res$split_method == sm & cox_res$P_value < 0.05, na.rm = TRUE)
    }, summary_grid$mp_group, summary_grid$split_method),
    "/",
    mapply(function(g, sm) {
      sum(cox_res$mp_group == g & cox_res$split_method == sm, na.rm = TRUE)
    }, summary_grid$mp_group, summary_grid$split_method)
  )
)
writeLines(summary_lines, file.path(tiers$reports, "Auto_parse_highres_mp_tcga_survival_summary.txt"))

script_run_status <- "success"
parse_finish_run(script_run, status = script_run_status, reused_cached = reuse_gsva)
message("Saved high-resolution MP TCGA survival volcano outputs to: ", out_dir)
