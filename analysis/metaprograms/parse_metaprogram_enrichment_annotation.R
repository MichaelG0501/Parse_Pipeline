####################
# parse_metaprogram_enrichment_annotation.R
#
# Description:
#   Annotates selected Parse metaprograms against Hallmark, GO:BP, 3CA, and
#   developmental-stage reference gene sets.
#
# Inputs:
#   parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds
#   /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/*.rds
#
# Outputs:
#   parse_outs/Auto_parse_metaprograms/Auto_parse_cluster_enrich.rds
#   parse_outs/Auto_parse_metaprograms/Auto_parse_enrichment_annotation.pdf
#   parse_outs/Auto_parse_metaprograms/Auto_parse_enrich_*.png
#   parse_outs/logs/run_summaries/parse_metaprogram_enrichment_annotation_*.txt
#
# Cache / replot:
#   Uses the selected metaprogram RDS from optimal-nMP selection. Plot-only
#   changes can be regenerated from Auto_parse_cluster_enrich.rds after the
#   enrichment object exists.
#
# Methodology:
#   analysis/methodology/metaprograms/metaprogram_enrichment_annotation_methodology.md
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_metaprogram_enrichment_annotation",
  input_files = c(
    "parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds",
    parse_reference_paths$three_ca_metaprograms,
    file.path(parse_reference_paths$developmental_enrichment_dir, "*.rds")
  ),
  output_files = c(
    "parse_outs/Auto_parse_metaprograms/Auto_parse_cluster_enrich.rds",
    "parse_outs/Auto_parse_metaprograms/Auto_parse_enrichment_annotation.pdf",
    "parse_outs/Auto_parse_metaprograms/Auto_parse_enrich_*.png"
  )
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath(file.path(parse_project_root(), "analysis/metaprograms/parse_metaprogram_enrichment_annotation.R"))
}

script_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(script_dir, "..", ".."))
out_dir <- file.path(project_dir, "parse_outs")
setwd(out_dir)

parse_mp_dir <- file.path("Auto_parse_metaprograms")
if (!dir.exists(parse_mp_dir)) {
  stop("Missing Parse metaprogram output dir: ", parse_mp_dir)
}

####################
# 1. Load reference gene sets
####################

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]

MP_list <- read.csv(parse_reference_paths$three_ca_metaprograms, check.names = FALSE)
MP_list <- as.list(MP_list)
MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])
mp_term2gene <- data.frame(
  term = rep(names(MP_list), lengths(MP_list)),
  gene = unlist(MP_list),
  row.names = NULL
)
mp_term2gene$term <- sub("^MP", "3CA_mp", mp_term2gene$term)
mp_term2name <- data.frame(term = unique(mp_term2gene$term), name = unique(mp_term2gene$term))

individual_dir <- parse_reference_paths$developmental_enrichment_dir
custom_files <- list.files(individual_dir, pattern = "\\.rds$", full.names = TRUE)
custom_refs <- lapply(custom_files, readRDS)
names(custom_refs) <- sub(".*enrich_dev_", "", basename(custom_files)) %>% sub("\\.rds$", "", .)

####################
# 2. Load metaprogram results
####################

mp_default_path <- file.path(parse_mp_dir, "Auto_parse_MP_outs_default.rds")
if (!file.exists(mp_default_path)) {
  stop("Missing selected metaprogram object. Run parse_metaprogram_select_optimal_nMP.R first: ", mp_default_path)
}

geneNMF.metaprograms <- readRDS(mp_default_path)
mp_gene_lists <- geneNMF.metaprograms$metaprograms.genes
mp_assignments <- geneNMF.metaprograms$programs.clusters

bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
mp_gene_lists <- mp_gene_lists[!names(mp_gene_lists) %in% bad_mp_names]

coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
mp_gene_lists <- mp_gene_lists[!names(mp_gene_lists) %in% low_coverage_mps]

valid_cluster_ids <- as.numeric(gsub("\\D", "", names(mp_gene_lists)))
mp_assignments <- mp_assignments[mp_assignments %in% valid_cluster_ids & !is.na(mp_assignments)]

message(paste0("Filtered out by silhouette < 0: ", paste(bad_mp_names, collapse = ", ")))
message(paste0("Filtered out by sampleCoverage < 0.25: ", paste(low_coverage_mps, collapse = ", ")))
message(paste0("Using ", length(mp_gene_lists), " metaprograms after filtering"))

####################
# Placeholder: update these labels after enrichment interpretation.
####################
mp_descriptions <- setNames(names(mp_gene_lists), names(mp_gene_lists))

####################
# 3. Enrichment analysis per metaprogram
####################

cluster_enrich <- lapply(names(mp_gene_lists), function(mp_name) {
  genes <- mp_gene_lists[[mp_name]]
  mp_id <- as.numeric(gsub("\\D", "", mp_name))
  members <- names(mp_assignments)[mp_assignments == mp_id]

  message(paste0("Processing MP: ", mp_name))

  res_GO <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  res_H <- enricher(gene = genes, TERM2GENE = hallmark_term2gene, TERM2NAME = hallmark_term2name, qvalueCutoff = 0.05)
  res_M <- enricher(gene = genes, TERM2GENE = mp_term2gene, TERM2NAME = mp_term2name, qvalueCutoff = 0.05)

  res_custom_list <- lapply(names(custom_refs), function(ref_name) {
    message(paste0("  -> Running custom enrichment: ", ref_name))
    enricher(
      gene = genes,
      TERM2GENE = custom_refs[[ref_name]]$TERM2GENE,
      TERM2NAME = custom_refs[[ref_name]]$TERM2NAME,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )
  })
  names(res_custom_list) <- names(custom_refs)

  base_results <- list(
    rep_prog = mp_name,
    members = members,
    genes = genes,
    GO = res_GO,
    Hallmark = res_H,
    MPs_3CA = res_M
  )

  c(base_results, res_custom_list)
})

names(cluster_enrich) <- names(mp_gene_lists)
saveRDS(cluster_enrich, file.path(parse_mp_dir, "Auto_parse_cluster_enrich.rds"), compress = FALSE)

####################
# 4. Heatmap visualization function
####################

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order <- rev(mp_tree_order)

enrich_heatmap <- function(cluster_enrich, element,
                           top_per_program = 8, top_n = 80, cap = 7,
                           cols = viridis::magma(100, direction = -1),
                           fontsize_row = 7, fontsize_col = 9) {
  is_custom <- !element %in% c("GO", "Hallmark", "MPs_3CA")

  df_list <- lapply(names(cluster_enrich), function(prog) {
    er <- cluster_enrich[[prog]][[element]]
    if (is.null(er)) return(NULL)

    r <- tryCatch(er@result, error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0) return(NULL)

    r_sig <- r[which(r$p.adjust < 0.05 & r$p.adjust > 0), ]
    data_source <- if (is_custom) r else r_sig
    if (nrow(data_source) == 0 && !is_custom) return(NULL)

    term <- if ("Description" %in% colnames(data_source)) data_source$Description else data_source$ID
    data.frame(
      Program = prog,
      Term = term,
      padj = data_source$p.adjust,
      Overlap = data_source$GeneRatio,
      stringsAsFactors = FALSE
    )
  })

  df <- dplyr::bind_rows(df_list)
  if (is.null(df) || nrow(df) == 0) {
    df <- data.frame(Program = character(), Term = character(), padj = numeric(), Overlap = character(), stringsAsFactors = FALSE)
  }

  if (is_custom) {
    if (!element %in% names(custom_refs)) {
      message("Custom reference not found for element: ", element)
      return(invisible(NULL))
    }
    terms_use <- as.character(custom_refs[[element]]$TERM2NAME$term)
  } else {
    if (nrow(df) == 0) {
      message("No significant results found for: ", element)
      return(invisible(NULL))
    }
    terms_use <- df %>%
      dplyr::filter(padj < 0.05) %>%
      dplyr::arrange(Program, padj) %>%
      dplyr::group_by(Program) %>%
      dplyr::slice_head(n = top_per_program) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(Term) %>%
      dplyr::pull(Term)

    if (length(terms_use) > top_n) {
      terms_use <- df %>%
        dplyr::filter(Term %in% terms_use) %>%
        dplyr::group_by(Term) %>%
        dplyr::summarise(min_p = min(padj), .groups = "drop") %>%
        dplyr::arrange(min_p) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::pull(Term)
    }
  }

  ordered_mps <- paste0("MP", mp_tree_order)
  ordered_mps <- ordered_mps[ordered_mps %in% names(cluster_enrich)]
  full_grid <- expand.grid(Term = terms_use, Program = ordered_mps, stringsAsFactors = FALSE)

  final_df <- full_grid %>%
    dplyr::left_join(df, by = c("Term", "Program")) %>%
    dplyr::mutate(
      score = tidyr::replace_na(pmin(-log10(padj), cap), 0),
      display_text = tidyr::replace_na(Overlap, "")
    )

  mat <- final_df %>%
    dplyr::select(Term, Program, score) %>%
    tidyr::pivot_wider(names_from = Program, values_from = score) %>%
    as.data.frame() %>%
    { row.names(.) <- .$Term; . } %>%
    dplyr::select(-Term) %>%
    as.matrix()

  text_mat <- final_df %>%
    dplyr::select(Term, Program, display_text) %>%
    tidyr::pivot_wider(names_from = Program, values_from = display_text) %>%
    as.data.frame() %>%
    { row.names(.) <- .$Term; . } %>%
    dplyr::select(-Term) %>%
    as.matrix()

  mat <- mat[terms_use, ordered_mps[ordered_mps %in% colnames(mat)], drop = FALSE]
  text_mat <- text_mat[terms_use, colnames(mat), drop = FALSE]

  if (nrow(mat) == 0 || ncol(mat) == 0) {
    message("No matrix content for element: ", element)
    return(invisible(NULL))
  }

  mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mp_sizes <- sapply(colnames(mat), function(x) length(mp_gene_lists[[x]]))
  display_mp <- ifelse(colnames(mat) %in% names(mp_descriptions), mp_descriptions[colnames(mat)], colnames(mat))
  col_labels <- paste0(display_mp, "\nn=", mp_sizes)

  cluster_rows_param <- FALSE
  row_gaps <- NULL
  if (is_custom) {
    mat <- mat[terms_use, , drop = FALSE]
    text_mat <- text_mat[terms_use, , drop = FALSE]
  } else {
    best_mp <- colnames(mat)[max.col(mat, ties.method = "first")]
    row_order <- order(match(best_mp, colnames(mat)), -rowSums(mat))
    mat <- mat[row_order, , drop = FALSE]
    text_mat <- text_mat[row_order, , drop = FALSE]
    groups <- colnames(mat)[max.col(mat, ties.method = "first")]
    row_gaps <- which(groups[-length(groups)] != groups[-1])
  }

  breaks <- seq(0, cap, length.out = length(cols) + 1)
  pheatmap::pheatmap(
    mat,
    display_numbers = text_mat,
    number_color = "black",
    fontsize_number = fontsize_row * 1.1,
    labels_col = col_labels,
    color = cols,
    breaks = breaks,
    cluster_rows = cluster_rows_param,
    cluster_cols = FALSE,
    gaps_row = row_gaps,
    border_color = NA,
    show_colnames = TRUE,
    angle_col = 0,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    main = paste0(element, " Enrichment (-log10 padj)")
  )

  invisible(mat)
}

####################
# 5. Generate enrichment heatmaps
####################

cols_palette <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)

pdf(file.path(parse_mp_dir, "Auto_parse_enrichment_annotation.pdf"), width = 10, height = 8)
enrich_heatmap(cluster_enrich, "Hallmark", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "GO", top_per_program = 6, top_n = 60, cols = cols_palette)
enrich_heatmap(cluster_enrich, "MPs_3CA", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Early_Embryogenesis", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Normal_Development_long", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Normal_Development_short", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Organogenesis_major", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Organogenesis_sub", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Adult_Epithelium", top_per_program = 8, top_n = 80, cols = cols_palette)
enrich_heatmap(cluster_enrich, "Barretts_Oesophagus", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png(file.path(parse_mp_dir, "Auto_parse_enrich_Hallmark.png"), width = 2000, height = 1750, res = 300)
enrich_heatmap(cluster_enrich, "Hallmark", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png(file.path(parse_mp_dir, "Auto_parse_enrich_GO.png"), width = 2300, height = 2000, res = 300)
enrich_heatmap(cluster_enrich, "GO", top_per_program = 6, top_n = 60, cols = cols_palette)
dev.off()

png(file.path(parse_mp_dir, "Auto_parse_enrich_MPs_3CA.png"), width = 2000, height = 1800, res = 300)
enrich_heatmap(cluster_enrich, "MPs_3CA", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

for (element in intersect(names(custom_refs), c(
  "Early_Embryogenesis",
  "Normal_Development_long",
  "Normal_Development_short",
  "Organogenesis_major",
  "Organogenesis_sub",
  "Adult_Epithelium",
  "Barretts_Oesophagus"
))) {
  png(file.path(parse_mp_dir, paste0("Auto_parse_enrich_", element, ".png")), width = 2500, height = 1900, res = 300)
  enrich_heatmap(cluster_enrich, element, top_per_program = 8, top_n = 80, cols = cols_palette)
  dev.off()
}

script_run_status <- "success"
message("All Parse enrichment heatmaps saved to ", parse_mp_dir)
