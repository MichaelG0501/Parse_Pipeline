####################
# parse_t2t4_vs_t0_dge_pathway_response.R
#
# Description:
#   Parse treatment-response differential expression and pathway summary
#   workflow. Compares T2+T4 cells against T0 cells without batch covariates,
#   then visualises the result with volcano plots, Hallmark pathway enrichment,
#   and a PDO-FLOT-inspired pathway/metric heatmap where columns are Parse
#   samples and values are real pseudobulk scores rather than treated-untreated
#   deltas.
#
# Inputs:
#   parse_outs/by_samples/<sample>/Auto_<sample>_final.rds
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
#   /rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds
#   Hallmark gene sets from msigdbr
#
# Outputs:
#   parse_outs/t2t4_vs_t0_response/tables/parse_t2t4_vs_t0_differential_genes.csv
#   parse_outs/t2t4_vs_t0_response/tables/parse_t2t4_vs_t0_hallmark_fgsea.csv
#   parse_outs/t2t4_vs_t0_response/tables/parse_t2t4_vs_t0_sample_pathway_metric_scores.csv
#   parse_outs/t2t4_vs_t0_response/figures/parse_t2t4_vs_t0_volcano.{pdf,png}
#   parse_outs/t2t4_vs_t0_response/figures/parse_t2t4_vs_t0_hallmark_fgsea_dotplot.{pdf,png}
#   parse_outs/t2t4_vs_t0_response/figures/parse_t2t4_vs_t0_top_deg_sample_heatmap.pdf
#   parse_outs/t2t4_vs_t0_response/figures/parse_sample_real_value_pathway_metric_heatmap.pdf
#   parse_outs/t2t4_vs_t0_response/intermediate/parse_t2t4_vs_t0_response_results.rds
#   parse_outs/logs/run_summaries/parse_t2t4_vs_t0_dge_pathway_response_*.txt
#
# Cache / replot:
#   DEG, FGSEA, pseudobulk matrices, and pathway score tables are written to
#   CSV/RDS. The script recomputes them on each run so plot changes remain
#   reproducible from the current merged Seurat object.
#
# Methodology:
#   analysis/methodology/cell_states/t2t4_vs_t0_dge_pathway_response_methodology.md
#
# Downstream:
#   Terminal response figures and source tables. Not currently consumed by
#   active downstream scripts.
####################

####################
# setup
####################
source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_helpers.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_t2t4_vs_t0_dge_pathway_response",
  parameters = list(
    contrast = "T2+T4 vs T0",
    deg_test = "matrixStats approximate Wilcoxon on downsampled log-normalized cells",
    min_pct = 0.05,
    logfc_threshold = 0,
    max_cells_per_ident = 1500,
    fgsea_collection = "MSigDB Hallmark",
    heatmap_value = "sample real pseudobulk score"
  ),
  input_files = c(
    "parse_outs/by_samples/<sample>/Auto_<sample>_final.rds",
    parse_reference_paths$cell_cycle_genes,
    parse_reference_paths$pdo_metaprograms
  ),
  output_files = c(
    "parse_outs/t2t4_vs_t0_response/tables/parse_t2t4_vs_t0_differential_genes.csv",
    "parse_outs/t2t4_vs_t0_response/tables/parse_t2t4_vs_t0_hallmark_fgsea.csv",
    "parse_outs/t2t4_vs_t0_response/figures/parse_t2t4_vs_t0_volcano.pdf",
    "parse_outs/t2t4_vs_t0_response/figures/parse_sample_real_value_pathway_metric_heatmap.pdf"
  )
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(scales)
  library(grid)
})

root_dir <- parse_project_root()
paths <- parse_paths(root_dir)
out_dir <- file.path(paths$parse_outs, "t2t4_vs_t0_response")
tiers <- parse_output_tiers(out_dir, create = TRUE)
setwd(root_dir)

sample_order <- parse_samples
contrast_samples <- c("T0", "T2", "T4")
treated_group <- c("T2", "T4")
baseline_group <- "T0"
min_pct <- 0.05
logfc_threshold <- 0
max_cells_per_ident <- 1500
volcano_logfc_cut <- 0.25
volcano_fdr_cut <- 0.05

palette_contract <- c(
  neutral_dark = "#272727",
  neutral_mid = "#767676",
  neutral_light = "#D8D8D8",
  signal_blue = "#3182BD",
  signal_teal = "#33B5A5",
  accent_red = "#D24B40",
  accent_orange = "#E28E2C"
)

theme_nature_contract <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.7),
      strip.text = element_text(size = base_size - 0.3, face = "bold"),
      plot.title = element_text(size = base_size + 0.8, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "grey30"),
      panel.grid = element_blank()
    )
}
theme_set(theme_nature_contract())

save_pub_pdf_png <- function(plot, filename, width_mm = 183, height_mm = 120, dpi = 450) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
}

save_heatmap_pdf <- function(heatmap_obj, filename, width_mm = 183, height_mm = 150, ...) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  grDevices::cairo_pdf(filename, width = w, height = h, family = "Arial")
  ComplexHeatmap::draw(heatmap_obj, ...)
  dev.off()
}

get_assay_data_safe <- function(obj, assay = "RNA", layer = "counts") {
  suppressWarnings({
    tryCatch(
      SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
      error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = layer)
    )
  })
}

mean_score <- function(mat, genes) {
  genes_use <- intersect(unique(genes), rownames(mat))
  if (length(genes_use) == 0) return(rep(NA_real_, ncol(mat)))
  Matrix::colMeans(mat[genes_use, , drop = FALSE], na.rm = TRUE)
}

zscore_rows <- function(mat) {
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  mat_z[!is.finite(mat_z)] <- 0
  mat_z
}

clean_pathway_label <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  stringr::str_to_sentence(tolower(x))
}

sig_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

get_hallmark_table <- function() {
  out <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) NULL
  )
  if (is.null(out) || nrow(out) == 0) {
    out <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  }
  out
}

pathway_ids <- c(
  "HALLMARK_E2F_TARGETS" = "E2F targets",
  "HALLMARK_G2M_CHECKPOINT" = "G2M checkpoint",
  "HALLMARK_APOPTOSIS" = "Apoptosis",
  "HALLMARK_P53_PATHWAY" = "p53 pathway",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "TNFa/NFkB",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_HYPOXIA" = "Hypoxia",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" = "Unfolded protein response",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = "Oxidative phosphorylation",
  "HALLMARK_DNA_REPAIR" = "DNA repair",
  "HALLMARK_XENOBIOTIC_METABOLISM" = "Xenobiotic metabolism"
)

pathway_block_df <- data.frame(
  pathway = c(
    "E2F targets", "G2M checkpoint",
    "Apoptosis", "p53 pathway", "DNA repair",
    "TNFa/NFkB", "EMT", "Hypoxia", "Xenobiotic metabolism", "Interferon response",
    "Unfolded protein response", "Oxidative phosphorylation"
  ),
  block = c(
    rep("Cell-cycle", 2), rep("Injury", 3),
    rep("Adaptive", 5), rep("Proteostasis", 2)
  ),
  stringsAsFactors = FALSE
)
pathway_order_grouped <- pathway_block_df$pathway

composite_pathways <- list(
  "Metric: Cell-cycle / proliferation" = c("E2F targets", "G2M checkpoint"),
  "Metric: Injury / checkpoint" = c("Apoptosis", "p53 pathway", "DNA repair"),
  "Metric: Adaptive / persistence" = c("TNFa/NFkB", "EMT", "Hypoxia", "Xenobiotic metabolism", "Interferon response"),
  "Metric: Proteostasis / transition" = c("Unfolded protein response", "Oxidative phosphorylation")
)

block_cols <- c(
  "Cell-cycle" = "#F0C75E",
  "Injury" = "#E07B54",
  "Adaptive" = "#5DAA68",
  "Proteostasis" = "#5B9BD5",
  "CCSIG" = "#7A4E9D",
  "Lineage" = "#2A788E",
  "Composite" = "#6B7280"
)

sample_group_cols <- c(
  "T0 baseline" = "#0072B2",
  "T2/T4 contrast" = "#D24B40",
  "Other Parse sample" = "#9CA3AF"
)

####################
# sample file helpers
####################
sample_files <- file.path(paths$parse_outs, "by_samples", sample_order, paste0("Auto_", sample_order, "_final.rds"))
names(sample_files) <- sample_order
missing_sample_files <- sample_files[!file.exists(sample_files)]
if (length(missing_sample_files) > 0) {
  stop("Missing final per-sample RDS files: ", paste(missing_sample_files, collapse = ", "))
}

####################
# differential expression
####################
message("Running T2+T4 vs T0 differential expression ...")
set.seed(1090)
contrast_data <- lapply(contrast_samples, function(sample) {
  message("  Loading ", sample)
  obj <- readRDS(sample_files[[sample]])
  DefaultAssay(obj) <- "RNA"
  obj$orig.ident <- sample
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    obj <- tryCatch(SeuratObject::JoinLayers(obj), error = function(e) obj)
  }
  data_mat <- get_assay_data_safe(obj, assay = "RNA", layer = "data")
  if (ncol(data_mat) > max_cells_per_ident) {
    keep_cells <- sample(colnames(data_mat), max_cells_per_ident)
  } else {
    keep_cells <- colnames(data_mat)
  }
  out <- data_mat[, keep_cells, drop = FALSE]
  colnames(out) <- paste(sample, colnames(out), sep = "__")
  rm(obj, data_mat)
  gc()
  out
})
names(contrast_data) <- contrast_samples
common_de_genes <- Reduce(intersect, lapply(contrast_data, rownames))
contrast_data <- lapply(contrast_data, function(mat) mat[common_de_genes, , drop = FALSE])
expr_sparse <- do.call(cbind, contrast_data)
sample_for_cell <- sub("__.*$", "", colnames(expr_sparse))
group_for_cell <- ifelse(sample_for_cell %in% treated_group, "T2_T4", "T0")
group_for_cell <- factor(group_for_cell, levels = c("T0", "T2_T4"))
rm(contrast_data)
gc()

pct_t2t4 <- Matrix::rowMeans(expr_sparse[, group_for_cell == "T2_T4", drop = FALSE] > 0)
pct_t0 <- Matrix::rowMeans(expr_sparse[, group_for_cell == "T0", drop = FALSE] > 0)
genes_use <- names(which(pmax(pct_t2t4, pct_t0) >= min_pct))
message("  Testing ", length(genes_use), " genes after min.pct filter.")

expr_dense <- as.matrix(expr_sparse[genes_use, , drop = FALSE])
rm(expr_sparse)
gc()

group_t2t4 <- group_for_cell == "T2_T4"
group_t0 <- group_for_cell == "T0"
n_t2t4 <- sum(group_t2t4)
n_t0 <- sum(group_t0)

rank_mat <- matrixStats::rowRanks(expr_dense, ties.method = "average", preserveShape = TRUE)
rank_sum_t2t4 <- rowSums(rank_mat[, group_t2t4, drop = FALSE])
u_stat <- rank_sum_t2t4 - n_t2t4 * (n_t2t4 + 1) / 2
u_mean <- n_t2t4 * n_t0 / 2
u_sd <- sqrt(n_t2t4 * n_t0 * (n_t2t4 + n_t0 + 1) / 12)
z_score <- (u_stat - u_mean) / u_sd
p_val <- 2 * stats::pnorm(-abs(z_score))
rm(rank_mat)
gc()

mean_t2t4 <- rowMeans(expm1(expr_dense[, group_t2t4, drop = FALSE]))
mean_t0 <- rowMeans(expm1(expr_dense[, group_t0, drop = FALSE]))
log2fc <- log2((mean_t2t4 + 1e-9) / (mean_t0 + 1e-9))

deg_tbl <- data.frame(
  gene = rownames(expr_dense),
  p_val = p_val,
  avg_log2FC = log2fc,
  pct.1 = pct_t2t4[rownames(expr_dense)],
  pct.2 = pct_t0[rownames(expr_dense)],
  stringsAsFactors = FALSE
)
rm(expr_dense)
gc()

deg_tbl <- deg_tbl %>%
  mutate(
    p_val_adj = p.adjust(p_val, method = "BH"),
    log2FC = avg_log2FC,
    p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj),
    neg_log10_fdr = -log10(pmax(p_val_adj, .Machine$double.xmin)),
    pct_diff = pct.1 - pct.2,
    regulation = case_when(
      p_val_adj < volcano_fdr_cut & log2FC >= volcano_logfc_cut ~ "Higher in T2+T4",
      p_val_adj < volcano_fdr_cut & log2FC <= -volcano_logfc_cut ~ "Higher in T0",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(p_val_adj, desc(abs(log2FC)))

write.csv(deg_tbl, file.path(tiers$tables, "parse_t2t4_vs_t0_differential_genes.csv"), row.names = FALSE)
write.csv(
  deg_tbl %>% filter(regulation != "Not significant") %>% arrange(p_val_adj, desc(abs(log2FC))),
  file.path(tiers$tables, "parse_t2t4_vs_t0_significant_differential_genes.csv"),
  row.names = FALSE
)

label_genes <- deg_tbl %>%
  filter(regulation != "Not significant") %>%
  group_by(regulation) %>%
  slice_min(order_by = p_val_adj, n = 12, with_ties = FALSE) %>%
  ungroup()

volcano_cols <- c(
  "Higher in T2+T4" = unname(palette_contract["accent_red"]),
  "Higher in T0" = unname(palette_contract["signal_blue"]),
  "Not significant" = "#B8B8B8"
)

p_volcano <- ggplot(deg_tbl, aes(log2FC, neg_log10_fdr)) +
  geom_point(aes(colour = regulation), size = 0.8, alpha = 0.62, stroke = 0) +
  geom_vline(xintercept = c(-volcano_logfc_cut, volcano_logfc_cut), linetype = "dashed", linewidth = 0.25, colour = "grey55") +
  geom_hline(yintercept = -log10(volcano_fdr_cut), linetype = "dashed", linewidth = 0.25, colour = "grey55") +
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 2.0,
    min.segment.length = 0,
    max.overlaps = Inf,
    box.padding = 0.25,
    segment.size = 0.22,
    seed = 1090
  ) +
  scale_colour_manual(values = volcano_cols, drop = FALSE, name = NULL) +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  labs(
    title = "T2+T4 versus T0 differential genes",
    subtitle = "Approximate per-cell Wilcoxon test; positive log2FC indicates higher expression in T2+T4",
    x = "Average log2 fold-change",
    y = "-log10 adjusted P"
  ) +
  theme(legend.position = "top")

save_pub_pdf_png(
  p_volcano,
  file.path(tiers$figures, "parse_t2t4_vs_t0_volcano"),
  width_mm = 150,
  height_mm = 115
)

####################
# pathway enrichment
####################
message("Running Hallmark pathway enrichment ...")
hallmark_tbl <- get_hallmark_table()
hallmark_sets <- split(hallmark_tbl$gene_symbol, hallmark_tbl$gs_name)
hallmark_sets <- lapply(hallmark_sets, unique)

ranking_tbl <- deg_tbl %>%
  filter(is.finite(log2FC), !is.na(gene)) %>%
  arrange(desc(abs(log2FC))) %>%
  distinct(gene, .keep_all = TRUE)
ranking <- ranking_tbl$log2FC
names(ranking) <- ranking_tbl$gene
ranking <- sort(ranking, decreasing = TRUE)

fgsea_tbl <- fgsea::fgsea(
  pathways = hallmark_sets,
  stats = ranking,
  minSize = 15,
  maxSize = 500,
  eps = 0
) %>%
  as.data.frame() %>%
  mutate(
    pathway_label = clean_pathway_label(pathway),
    leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)),
    direction = ifelse(NES >= 0, "Enriched in T2+T4", "Enriched in T0")
  ) %>%
  arrange(padj, desc(abs(NES)))

write.csv(fgsea_tbl, file.path(tiers$tables, "parse_t2t4_vs_t0_hallmark_fgsea.csv"), row.names = FALSE)
write.csv(
  fgsea_tbl %>% select(pathway, pathway_label, direction, NES, pval, padj, size, leadingEdge),
  file.path(tiers$tables, "parse_t2t4_vs_t0_hallmark_leading_edge_summary.csv"),
  row.names = FALSE
)

fgsea_plot_tbl <- fgsea_tbl %>%
  filter(is.finite(NES), !is.na(padj)) %>%
  group_by(direction) %>%
  slice_min(order_by = padj, n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    pathway_label = factor(pathway_label, levels = pathway_label[order(NES)]),
    fdr_display = pmax(padj, 1e-10)
  )

p_fgsea <- ggplot(fgsea_plot_tbl, aes(NES, pathway_label)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(size = -log10(fdr_display), colour = NES), alpha = 0.92) +
  scale_colour_gradient2(
    low = palette_contract["signal_blue"],
    mid = "white",
    high = palette_contract["accent_red"],
    midpoint = 0,
    name = "NES"
  ) +
  scale_size_continuous(name = "-log10 FDR", range = c(1.3, 4.8)) +
  labs(
    title = "Hallmark pathway enrichment",
    subtitle = "Ranked by T2+T4 versus T0 average log2 fold-change",
    x = "Normalised enrichment score",
    y = NULL
  ) +
  theme(legend.position = "right")

save_pub_pdf_png(
  p_fgsea,
  file.path(tiers$figures, "parse_t2t4_vs_t0_hallmark_fgsea_dotplot"),
  width_mm = 170,
  height_mm = 130
)

####################
# pseudobulk sample scores
####################
message("Building sample-level pseudobulk pathway and metric scores ...")
pseudobulk_vectors <- lapply(sample_order, function(sample) {
  message("  Aggregating ", sample)
  obj <- readRDS(sample_files[[sample]])
  DefaultAssay(obj) <- "RNA"
  counts <- get_assay_data_safe(obj, assay = "RNA", layer = "counts")
  out <- Matrix::rowSums(counts)
  attr(out, "cell_n") <- ncol(counts)
  rm(obj, counts)
  gc()
  out
})
names(pseudobulk_vectors) <- sample_order
pseudobulk_cell_n <- vapply(pseudobulk_vectors, function(x) as.integer(attr(x, "cell_n")), FUN.VALUE = integer(1))
all_pb_genes <- Reduce(union, lapply(pseudobulk_vectors, names))
pseudobulk_counts <- matrix(0, nrow = length(all_pb_genes), ncol = length(sample_order), dimnames = list(all_pb_genes, sample_order))
for (sample in sample_order) {
  pseudobulk_counts[names(pseudobulk_vectors[[sample]]), sample] <- pseudobulk_vectors[[sample]]
}
storage.mode(pseudobulk_counts) <- "integer"
rm(pseudobulk_vectors)
gc()

pseudobulk_meta <- data.frame(
  sample = colnames(pseudobulk_counts),
  cell_n = as.integer(pseudobulk_cell_n[colnames(pseudobulk_counts)]),
  group = ifelse(colnames(pseudobulk_counts) == baseline_group, "T0 baseline",
    ifelse(colnames(pseudobulk_counts) %in% treated_group, "T2/T4 contrast", "Other Parse sample")
  ),
  stringsAsFactors = FALSE
)
write.csv(pseudobulk_meta, file.path(tiers$tables, "parse_sample_pseudobulk_metadata.csv"), row.names = FALSE)

pb_dge <- edgeR::DGEList(counts = pseudobulk_counts)
pb_dge <- edgeR::calcNormFactors(pb_dge)
pb_logcpm <- edgeR::cpm(pb_dge, log = TRUE, prior.count = 2)
pb_logcpm_z <- zscore_rows(pb_logcpm)

selected_gene_sets <- lapply(names(pathway_ids), function(gs_name) {
  unique(hallmark_tbl$gene_symbol[hallmark_tbl$gs_name == gs_name])
})
names(selected_gene_sets) <- unname(pathway_ids)
selected_gene_sets[["Interferon response"]] <- unique(
  hallmark_tbl$gene_symbol[
    hallmark_tbl$gs_name %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
  ]
)

pathway_score_mat <- sapply(selected_gene_sets, function(genes) mean_score(pb_logcpm_z, genes))
pathway_score_mat <- t(pathway_score_mat)
rownames(pathway_score_mat) <- names(selected_gene_sets)

cc_top50 <- character()
if (file.exists(parse_reference_paths$cell_cycle_genes)) {
  cc_genes_df <- read.csv(parse_reference_paths$cell_cycle_genes, header = TRUE, stringsAsFactors = FALSE)
  if (all(c("Gene", "Consensus") %in% colnames(cc_genes_df))) {
    cc_consensus <- cc_genes_df$Gene[cc_genes_df$Consensus == 1]
  } else {
    cc_consensus <- cc_genes_df[[1]]
  }
  cc_consensus <- intersect(unique(cc_consensus), rownames(pb_logcpm))
  if (length(cc_consensus) > 0) {
    cc_top50 <- names(sort(rowMeans(pb_logcpm[cc_consensus, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[seq_len(min(50, length(cc_consensus)))]
  }
}
if (length(cc_top50) == 0) {
  cc_top50 <- unique(c(selected_gene_sets[["E2F targets"]], selected_gene_sets[["G2M checkpoint"]]))
}

lineage_scores <- list()
pdo_geneNMF <- NULL
if (file.exists(parse_reference_paths$pdo_metaprograms)) {
  pdo_geneNMF <- readRDS(parse_reference_paths$pdo_metaprograms)
  pdo_mp_genes <- pdo_geneNMF$metaprograms.genes
  lineage_scores[["Intestinal Metaplasia"]] <- mean_score(pb_logcpm_z, pdo_mp_genes[["MP4"]])
  lineage_scores[["Columnar Progenitor"]] <- mean_score(pb_logcpm_z, pdo_mp_genes[["MP8"]])
}

extra_score_mat <- rbind(
  "CCSIG" = mean_score(pb_logcpm_z, cc_top50),
  do.call(rbind, lineage_scores)
)
colnames(extra_score_mat) <- colnames(pb_logcpm_z)

composite_score_mat <- sapply(names(composite_pathways), function(metric_name) {
  pathways <- intersect(composite_pathways[[metric_name]], rownames(pathway_score_mat))
  colMeans(pathway_score_mat[pathways, , drop = FALSE], na.rm = TRUE)
})
composite_score_mat <- t(composite_score_mat)
colnames(composite_score_mat) <- colnames(pathway_score_mat)

score_mat <- rbind(pathway_score_mat[pathway_order_grouped, sample_order, drop = FALSE], extra_score_mat[, sample_order, drop = FALSE], composite_score_mat[, sample_order, drop = FALSE])
score_blocks <- c(
  pathway_block_df$block[match(pathway_order_grouped, pathway_block_df$pathway)],
  "CCSIG",
  rep("Lineage", max(0, nrow(extra_score_mat) - 1)),
  rep("Composite", nrow(composite_score_mat))
)
names(score_blocks) <- rownames(score_mat)

score_long <- as.data.frame(score_mat) %>%
  rownames_to_column("feature") %>%
  pivot_longer(cols = all_of(sample_order), names_to = "sample", values_to = "score") %>%
  left_join(pseudobulk_meta, by = "sample") %>%
  mutate(block = score_blocks[feature])
write.csv(score_long, file.path(tiers$tables, "parse_t2t4_vs_t0_sample_pathway_metric_scores.csv"), row.names = FALSE)

contrast_score_summary <- score_long %>%
  filter(sample %in% contrast_samples) %>%
  select(feature, sample, score) %>%
  pivot_wider(names_from = sample, values_from = score) %>%
  mutate(
    T2T4_mean = rowMeans(cbind(T2, T4), na.rm = TRUE),
    T2T4_minus_T0 = T2T4_mean - T0,
    direction = case_when(
      T2T4_minus_T0 > 0 ~ "Higher in T2+T4",
      T2T4_minus_T0 < 0 ~ "Higher in T0",
      TRUE ~ "No change"
    )
  ) %>%
  arrange(desc(abs(T2T4_minus_T0)))
write.csv(contrast_score_summary, file.path(tiers$tables, "parse_t2t4_vs_t0_pathway_metric_contrast_summary.csv"), row.names = FALSE)

sample_group <- pseudobulk_meta$group[match(sample_order, pseudobulk_meta$sample)]
names(sample_group) <- sample_order
top_ha <- HeatmapAnnotation(
  Group = sample_group,
  `Cells` = anno_barplot(
    pseudobulk_meta$cell_n[match(sample_order, pseudobulk_meta$sample)],
    gp = gpar(fill = "#6B7280", col = NA),
    height = unit(14, "mm"),
    axis_param = list(gp = gpar(fontsize = 5.5))
  ),
  col = list(Group = sample_group_cols),
  annotation_name_gp = gpar(fontface = "bold", fontsize = 6.5),
  show_annotation_name = TRUE
)

row_ha <- rowAnnotation(
  Block = score_blocks,
  col = list(Block = block_cols),
  show_annotation_name = FALSE,
  annotation_legend_param = list(Block = list(title = "Feature block"))
)

score_clip <- max(0.3, quantile(abs(score_mat), 0.95, na.rm = TRUE))
score_col_fun <- colorRamp2(
  c(-score_clip, 0, score_clip),
  c("#245F7B", "white", "#B63E2F")
)

cell_fun_score <- function(j, i, x, y, w, h, fill) {
  val <- score_mat[i, j]
  if (is.finite(val)) {
    grid.text(sprintf("%.2f", val), x, y, gp = gpar(fontsize = 5.5, col = "black"))
  }
}

ht_scores <- Heatmap(
  score_mat,
  name = "Score",
  col = score_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = factor(score_blocks, levels = names(block_cols)),
  row_gap = unit(2.2, "mm"),
  left_annotation = row_ha,
  top_annotation = top_ha,
  row_names_gp = gpar(fontsize = 6.5, fontface = "bold"),
  column_names_gp = gpar(fontsize = 7.2, fontface = "bold"),
  column_title = "Parse sample pathway and metric scores",
  column_title_gp = gpar(fontface = "bold", fontsize = 8),
  heatmap_legend_param = list(title = "Real\nscore"),
  cell_fun = cell_fun_score
)

save_heatmap_pdf(
  ht_scores,
  file.path(tiers$figures, "parse_sample_real_value_pathway_metric_heatmap.pdf"),
  width_mm = 175,
  height_mm = 155,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

####################
# top DEG pseudobulk heatmap
####################
message("Writing top DEG sample heatmap ...")
top_heat_genes <- deg_tbl %>%
  filter(regulation != "Not significant") %>%
  group_by(regulation) %>%
  slice_min(order_by = p_val_adj, n = 18, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(regulation, p_val_adj) %>%
  pull(gene) %>%
  unique()

if (length(top_heat_genes) > 0) {
  heat_gene_mat <- zscore_rows(pb_logcpm[intersect(top_heat_genes, rownames(pb_logcpm)), sample_order, drop = FALSE])
  gene_direction <- deg_tbl$regulation[match(rownames(heat_gene_mat), deg_tbl$gene)]
  names(gene_direction) <- rownames(heat_gene_mat)
  direction_cols <- c(
    "Higher in T2+T4" = unname(palette_contract["accent_red"]),
    "Higher in T0" = unname(palette_contract["signal_blue"])
  )
  gene_ha <- rowAnnotation(
    Direction = gene_direction,
    col = list(Direction = direction_cols),
    annotation_name_gp = gpar(fontface = "bold", fontsize = 6.5)
  )
  heat_clip <- max(0.5, quantile(abs(heat_gene_mat), 0.95, na.rm = TRUE))
  ht_deg <- Heatmap(
    heat_gene_mat,
    name = "Gene z",
    col = colorRamp2(c(-heat_clip, 0, heat_clip), c("#245F7B", "white", "#B63E2F")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = gene_ha,
    top_annotation = top_ha,
    row_split = factor(gene_direction, levels = c("Higher in T2+T4", "Higher in T0")),
    row_names_gp = gpar(fontsize = 5.8),
    column_names_gp = gpar(fontsize = 7.2, fontface = "bold"),
    column_title = "Top differential genes across Parse samples",
    column_title_gp = gpar(fontface = "bold", fontsize = 8)
  )
  save_heatmap_pdf(
    ht_deg,
    file.path(tiers$figures, "parse_t2t4_vs_t0_top_deg_sample_heatmap.pdf"),
    width_mm = 150,
    height_mm = 150,
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
}

####################
# cache and report
####################
saveRDS(
  list(
    deg = deg_tbl,
    fgsea = fgsea_tbl,
    pseudobulk_meta = pseudobulk_meta,
    pseudobulk_logcpm = pb_logcpm,
    pathway_metric_scores = score_long,
    pathway_metric_matrix = score_mat,
    contrast_score_summary = contrast_score_summary
  ),
  file = file.path(tiers$intermediate, "parse_t2t4_vs_t0_response_results.rds")
)

writeLines(
  c(
    "# Parse T2+T4 vs T0 response summary",
    "",
    "## Contrast",
    "- T2 and T4 cells were pooled and compared against T0 cells.",
    "- No batch covariate was included, matching the sequencing-together assumption.",
    "",
    "## Main outputs",
    "- `figures/parse_t2t4_vs_t0_volcano.pdf`",
    "- `figures/parse_t2t4_vs_t0_hallmark_fgsea_dotplot.pdf`",
    "- `figures/parse_sample_real_value_pathway_metric_heatmap.pdf`",
    "- `tables/parse_t2t4_vs_t0_differential_genes.csv`",
    "- `tables/parse_t2t4_vs_t0_hallmark_fgsea.csv`",
    "- `tables/parse_t2t4_vs_t0_sample_pathway_metric_scores.csv`"
  ),
  con = file.path(tiers$reports, "parse_t2t4_vs_t0_response_summary.md")
)

script_run_status <- "success"
parse_finish_run(script_run, status = script_run_status)
message("T2+T4 vs T0 response workflow complete.")
