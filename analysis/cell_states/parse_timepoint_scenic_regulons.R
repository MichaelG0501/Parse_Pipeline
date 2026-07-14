####################
# parse_timepoint_scenic_regulons.R
#
# Description:
#   Run SCENIC regulon discovery and regulon activity analysis for Parse
#   treatment-response timepoints. This workflow ignores cell-state and
#   metaprogram concepts completely: the biological entities are the six
#   Parse timepoints T0, T1, T2, T4, R4, and eR4.
#
# Inputs:
#   parse_outs/Auto_parse_merged.rds
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/cistarget_databases_rcistarget_mc9nr/
#     hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
#     hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
#
# Outputs:
#   parse_outs/cell_states/timepoint_scenic/intermediate/scenic_combined/
#   parse_outs/cell_states/timepoint_scenic/intermediate/*regulon_auc.rds
#   parse_outs/cell_states/timepoint_scenic/intermediate/*rss.rds
#   parse_outs/cell_states/timepoint_scenic/tables/*selected_cells.csv
#   parse_outs/cell_states/timepoint_scenic/tables/*differential_regulon_activity.csv
#   parse_outs/cell_states/timepoint_scenic/tables/*network_edges.csv
#   parse_outs/cell_states/timepoint_scenic/tables/*regulon_targets.csv
#   parse_outs/cell_states/timepoint_scenic/figures/*regulon_heatmap.pdf
#   parse_outs/cell_states/timepoint_scenic/figures/*network.pdf
#   parse_outs/cell_states/timepoint_scenic/figures/*overlap*.pdf
#   parse_outs/logs/run_summaries/parse_timepoint_scenic_regulons_<timestamp>.txt
#
# Cache/replot behavior:
#   SCENIC intermediate files under intermediate/scenic_* are reused when
#   present. Tables and figures are regenerated from cached SCENIC objects.
#
# Methodology:
#   analysis/methodology/cell_states/timepoint_scenic_regulons_methodology.md
#
# Downstream status:
#   Active terminal workflow for timepoint-level regulon discovery and
#   differential regulon activity comparison. No active downstream consumers.
#
# Usage:
#   Rscript analysis/cell_states/parse_timepoint_scenic_regulons.R run_mode=preflight
#   Rscript analysis/cell_states/parse_timepoint_scenic_regulons.R run_mode=combined
#   Rscript analysis/cell_states/parse_timepoint_scenic_regulons.R run_mode=combined genie3_nparts=100
####################

####################
# Package setup and shared configuration
####################
root_dir <- Sys.getenv(
  "AUTO_PARSE_ROOT_DIR",
  unset = "/rds/general/project/spatialtranscriptomics/ephemeral/Parse_Pipeline"
)
if (!dir.exists(root_dir)) {
  root_dir <- "/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline"
}
if (!dir.exists(root_dir)) stop("Parse_Pipeline root is missing. Ask the user what to do.")
setwd(root_dir)

source(file.path(root_dir, "analysis/common/parse_pipeline_config.R"))
source(file.path(root_dir, "analysis/common/parse_pipeline_helpers.R"))
source(file.path(root_dir, "analysis/common/parse_pipeline_logging.R"))

parse_load_or_stop(c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ComplexHeatmap",
  "circlize",
  "Matrix",
  "data.table",
  "scales",
  "igraph",
  "ggraph",
  "tidygraph",
  "grid",
  "SCENIC",
  "AUCell",
  "RcisTarget",
  "GENIE3",
  "doRNG",
  "doMC"
))

####################
# Small utilities
####################
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || !nzchar(as.character(x[1]))) return(y)
  x[1]
}

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) next
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    out[[parts[1]]] <- paste(parts[-1], collapse = "=")
  }
  out
}

format_regulon_name <- function(x) {
  x <- gsub("_extended$", "", x)
  x <- gsub(" \\([0-9]+g\\)$", "", x)
  gsub(" \\([0-9]+ genes\\)$", "", x)
}

check_file_or_stop <- function(path, label = path) {
  if (!file.exists(path)) stop("Missing required input file: ", label, ". Ask the user what to do.")
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}

check_dir_or_stop <- function(path, label = path) {
  if (!dir.exists(path)) stop("Missing required input directory: ", label, ". Ask the user what to do.")
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}

write_matrix_csv <- function(mat, path, row_label = "regulon") {
  df <- as.data.frame(mat, check.names = FALSE)
  df[[row_label]] <- rownames(mat)
  df <- df[, c(row_label, setdiff(colnames(df), row_label)), drop = FALSE]
  utils::write.csv(df, path, row.names = FALSE)
}

extract_regulon_targets <- function(x) {
  if (requireNamespace("GSEABase", quietly = TRUE) && methods::is(x, "GeneSet")) return(unique(GSEABase::geneIds(x)))
  if (is.character(x)) return(unique(x))
  if (is.list(x) && !is.null(x$gene)) return(unique(as.character(x$gene)))
  if (!is.null(names(x))) return(unique(names(x)))
  unique(as.character(x))
}

extract_counts_for_cells <- function(seurat_obj, cells, assay = "RNA") {
  counts_mat <- tryCatch(
    parse_get_counts(seurat_obj, assay = assay),
    error = function(e) NULL
  )
  if (!is.null(counts_mat)) {
    missing_cells <- setdiff(cells, colnames(counts_mat))
    if (length(missing_cells) > 0) stop("Counts matrix is missing selected cells. Ask the user what to do.")
    counts_mat <- counts_mat[, cells, drop = FALSE]
    if (!inherits(counts_mat, "dgCMatrix")) counts_mat <- methods::as(counts_mat, "dgCMatrix")
    return(counts_mat)
  }

  assay_obj <- seurat_obj[[assay]]
  count_layers <- SeuratObject::Layers(assay_obj, search = "^counts")
  if (length(count_layers) == 0) stop("No counts layers found in assay ", assay, ". Ask the user what to do.")
  layer_mats <- lapply(count_layers, function(layer_name) {
    mat <- SeuratObject::LayerData(seurat_obj, assay = assay, layer = layer_name)
    keep <- intersect(cells, colnames(mat))
    if (length(keep) == 0) return(NULL)
    mat[, keep, drop = FALSE]
  })
  layer_mats <- layer_mats[!vapply(layer_mats, is.null, logical(1))]
  if (length(layer_mats) == 0) stop("Selected cells were not found in any counts layer. Ask the user what to do.")

  found_cells <- unlist(lapply(layer_mats, colnames), use.names = FALSE)
  missing_cells <- setdiff(cells, found_cells)
  if (length(missing_cells) > 0) {
    stop("Counts layers are missing selected cells: ", paste(head(missing_cells, 20), collapse = ", "), ". Ask the user what to do.")
  }

  row_sets_identical <- all(vapply(layer_mats, function(mat) identical(rownames(mat), rownames(layer_mats[[1]])), logical(1)))
  if (row_sets_identical) {
    counts_mat <- do.call(cbind, layer_mats)
  } else {
    all_genes <- Reduce(union, lapply(layer_mats, rownames))
    layer_mats <- lapply(layer_mats, function(mat) {
      if (identical(rownames(mat), all_genes)) return(mat)
      expanded <- Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        dims = c(length(all_genes), ncol(mat)),
        dimnames = list(all_genes, colnames(mat))
      )
      expanded[rownames(mat), colnames(mat)] <- mat
      expanded
    })
    counts_mat <- do.call(cbind, layer_mats)
  }
  counts_mat <- counts_mat[, cells, drop = FALSE]
  if (!inherits(counts_mat, "dgCMatrix")) counts_mat <- methods::as(counts_mat, "dgCMatrix")
  counts_mat
}

detect_db_files <- function(db_dir) {
  check_dir_or_stop(db_dir, "SCENIC cisTarget database directory")
  db_files <- list.files(db_dir, pattern = "\\.feather$", full.names = FALSE)
  db_files <- db_files[grepl("hg38|refseq-r80|hgnc", db_files, ignore.case = TRUE)]
  preferred <- db_files[grepl("mc9nr|refseq-r80", db_files, ignore.case = TRUE)]
  if (length(preferred) > 0) db_files <- preferred
  db_500 <- db_files[grepl("500bp", db_files, ignore.case = TRUE)][1]
  db_10k <- db_files[grepl("10kb", db_files, ignore.case = TRUE)][1]
  db_files <- unique(c(db_500, db_10k, db_files))
  db_files <- db_files[!is.na(db_files)]
  if (length(db_files) == 0) stop("No human mc9nr cisTarget feather databases found in ", db_dir, ". Ask the user what to do.")
  invisible(vapply(file.path(db_dir, db_files), check_file_or_stop, character(1)))
  db_files
}

patch_scenic_annotation_lookup <- function() {
  scenic_ns <- asNamespace("SCENIC")
  original_fun <- get("getDbAnnotations", envir = scenic_ns)
  patched_fun <- original_fun
  body(patched_fun) <- quote({
    dbAnnotFiles <- scenicOptions@settings$db_annotFiles
    if (!is.null(dbAnnotFiles)) {
      motifAnnotations <- NULL
      for (annotPath in dbAnnotFiles) {
        motifAnnot <- data.table::fread(annotPath)
        motifAnnot$annotationSource <- factor(motifAnnot$annotationSource)
        colnames(motifAnnot)[1] <- "motif"
        levels(motifAnnot$annotationSource) <- c(
          levels(motifAnnot$annotationSource),
          c("directAnnotation", "inferredBy_Orthology", "inferredBy_MotifSimilarity", "inferredBy_MotifSimilarity_n_Orthology")
        )
        motifAnnotations <- rbind(motifAnnotations, motifAnnot)
      }
    } else {
      if (is.na(getDatasetInfo(scenicOptions, "org"))) stop("Please provide an organism (scenicOptions@inputDatasetInfo$org).")
      org <- getDatasetInfo(scenicOptions, "org")
      if (org == "hgnc") motifAnnotName <- "motifAnnotations_hgnc"
      if (org == "mgi") motifAnnotName <- "motifAnnotations_mgi"
      if (org == "dmel") motifAnnotName <- "motifAnnotations_dmel"
      if (!is.null(scenicOptions@settings$db_mcVersion) && scenicOptions@settings$db_mcVersion == "v8") motifAnnotName <- paste0(motifAnnotName, "_v8")
      annot_env <- new.env(parent = baseenv())
      data(list = motifAnnotName, package = "RcisTarget", envir = annot_env, verbose = FALSE)
      if (!exists(motifAnnotName, envir = annot_env, inherits = FALSE)) {
        v9_name <- paste0(motifAnnotName, "_v9")
        data(list = v9_name, package = "RcisTarget", envir = annot_env, verbose = FALSE)
        if (exists(v9_name, envir = annot_env, inherits = FALSE)) assign(motifAnnotName, get(v9_name, envir = annot_env), envir = annot_env)
      }
      motifAnnotations <- get(motifAnnotName, envir = annot_env, inherits = FALSE)
    }
    return(motifAnnotations)
  })
  unlockBinding("getDbAnnotations", scenic_ns)
  assign("getDbAnnotations", patched_fun, envir = scenic_ns)
  lockBinding("getDbAnnotations", scenic_ns)
  invisible(TRUE)
}

scenic_gene_filtering_sparse <- function(exprMat, scenicOptions, minCountsPerGene, minSamples) {
  dbFilePath <- getDatabases(scenicOptions)[[1]]
  outFile_genesKept <- getIntName(scenicOptions, "genesKept")
  if (is.null(dbFilePath)) stop("dbFilePath")
  if (is.data.frame(exprMat)) stop("data.frame expression matrices are not supported")
  if (any(table(rownames(exprMat)) > 1)) stop("Expression matrix rownames should be unique")
  if (inherits(exprMat, "Matrix") || inherits(exprMat, "sparseMatrix")) {
    nCountsPerGene <- Matrix::rowSums(exprMat, na.rm = TRUE)
    nCellsPerGene <- Matrix::rowSums(exprMat > 0, na.rm = TRUE)
  } else {
    nCountsPerGene <- rowSums(exprMat, na.rm = TRUE)
    nCellsPerGene <- rowSums(exprMat > 0, na.rm = TRUE)
  }
  genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minCountsPerGene)]
  nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
  genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
  motifRankings <- RcisTarget::importRankings(dbFilePath)
  genesInDatabase <- colnames(RcisTarget::getRanking(motifRankings))
  genesKept <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
  if (!is.null(outFile_genesKept)) saveRDS(genesKept, file = outFile_genesKept)
  genesKept
}

make_long_regulon_table <- function(mean_auc_mat, rss_mat, label_order) {
  rows <- lapply(label_order, function(label) {
    other_labels <- setdiff(label_order, label)
    data.frame(
      timepoint = label,
      regulon = rownames(mean_auc_mat),
      regulon_label = format_regulon_name(rownames(mean_auc_mat)),
      mean_auc = as.numeric(mean_auc_mat[, label]),
      mean_auc_other = if (length(other_labels) > 0) rowMeans(mean_auc_mat[, other_labels, drop = FALSE], na.rm = TRUE) else NA_real_,
      rss = as.numeric(rss_mat[, label]),
      rss_other = if (length(other_labels) > 0) rowMeans(rss_mat[, other_labels, drop = FALSE], na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  df <- dplyr::bind_rows(rows)
  df |>
    dplyr::mutate(
      delta_mean_auc = mean_auc - mean_auc_other,
      delta_rss = rss - rss_other
    ) |>
    dplyr::group_by(regulon) |>
    dplyr::mutate(
      mean_auc_z_across_timepoints = as.numeric(scale(mean_auc)),
      rss_z_across_timepoints = as.numeric(scale(rss))
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(timepoint) |>
    dplyr::arrange(dplyr::desc(rss), dplyr::desc(delta_mean_auc), .by_group = TRUE) |>
    dplyr::mutate(timepoint_rank_by_rss = dplyr::row_number()) |>
    dplyr::ungroup()
}

####################
# SCENIC execution
####################
run_scenic_for_cells <- function(counts_mat,
                                 label_map,
                                 label_order,
                                 label_colours,
                                 run_dir,
                                 out_prefix,
                                 dataset_title,
                                 db_dir,
                                 db_files,
                                 n_cores,
                                 genie3_nparts,
                                 genie3_resume,
                                 top_regulons_per_label,
                                 tiers) {
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  setwd(run_dir)
  on.exit(setwd(old_wd), add = TRUE)

  label_map <- label_map[colnames(counts_mat)]
  if (any(is.na(label_map))) stop("Some SCENIC cells lack timepoint labels in ", out_prefix, ". Ask the user what to do.")
  label_order <- label_order[label_order %in% unique(as.character(label_map))]
  if (length(label_order) == 0) stop("No labels left for SCENIC run ", out_prefix, ". Ask the user what to do.")

  scenicOptions <- SCENIC::initializeScenic(
    org = "hgnc",
    dbDir = db_dir,
    dbs = db_files,
    datasetTitle = dataset_title,
    nCores = n_cores
  )
  min_counts_per_gene <- max(3 * 0.01 * ncol(counts_mat), 20)
  min_samples <- max(0.01 * ncol(counts_mat), 20)
  genes_kept_path <- file.path("int", "1.1_genesKept.Rds")
  if (file.exists(genes_kept_path)) {
    genes_kept <- readRDS(genes_kept_path)
  } else {
    genes_kept <- scenic_gene_filtering_sparse(
      counts_mat,
      scenicOptions = scenicOptions,
      minCountsPerGene = min_counts_per_gene,
      minSamples = min_samples
    )
  }
  if (length(genes_kept) < 500) stop("SCENIC gene filtering retained fewer than 500 genes for ", out_prefix, ". Ask the user what to do.")

  expr_mat_use <- counts_mat[genes_kept, , drop = FALSE]
  if (!is.matrix(expr_mat_use)) expr_mat_use <- as.matrix(expr_mat_use)

  if (!file.exists(file.path("int", "1.2_corrMat.Rds"))) SCENIC::runCorrelation(expr_mat_use, scenicOptions)
  if (!file.exists(file.path("int", "1.4_GENIE3_linkList.Rds"))) {
    SCENIC::runGenie3(
      expr_mat_use,
      scenicOptions,
      nParts = genie3_nparts,
      resumePreviousRun = genie3_resume
    )
  }
  if (!file.exists(file.path("int", "1.6_tfModules_asDF.Rds"))) scenicOptions <- SCENIC::runSCENIC_1_coexNetwork2modules(scenicOptions)
  if (!(file.exists(file.path("int", "2.6_regulons_asGeneSet.Rds")) && file.exists(file.path("int", "2.6_regulons_asIncidMat.Rds")))) {
    scenicOptions <- SCENIC::runSCENIC_2_createRegulons(scenicOptions)
  }
  if (!file.exists(file.path("int", "3.4_regulonAUC.Rds"))) SCENIC::runSCENIC_3_scoreCells(scenicOptions, exprMat = counts_mat)

  regulon_auc <- SCENIC::loadInt(scenicOptions, "aucell_regulonAUC")
  regulons <- SCENIC::loadInt(scenicOptions, "regulons")
  auc_mat <- AUCell::getAUC(regulon_auc)
  label_map <- label_map[colnames(auc_mat)]

  mean_auc_mat <- as.matrix(sapply(label_order, function(label) {
    cells <- names(label_map)[label_map == label]
    rowMeans(auc_mat[, cells, drop = FALSE], na.rm = TRUE)
  }))
  colnames(mean_auc_mat) <- label_order

  if (length(label_order) >= 2) {
    rss_mat <- tryCatch(SCENIC::calcRSS(AUC = auc_mat, cellAnnotation = label_map), error = function(e) NULL)
    if (is.null(rss_mat)) rss_mat <- mean_auc_mat
    rss_mat <- as.matrix(rss_mat)
    rss_mat <- rss_mat[, label_order, drop = FALSE]
  } else {
    rss_mat <- mean_auc_mat
  }
  mean_auc_mat <- mean_auc_mat[, label_order, drop = FALSE]

  saveRDS(regulon_auc, file.path(tiers$intermediate, paste0(out_prefix, "_regulon_auc.rds")))
  saveRDS(regulons, file.path(tiers$intermediate, paste0(out_prefix, "_regulons.rds")))
  saveRDS(mean_auc_mat, file.path(tiers$intermediate, paste0(out_prefix, "_mean_auc.rds")))
  saveRDS(rss_mat, file.path(tiers$intermediate, paste0(out_prefix, "_rss.rds")))
  write_matrix_csv(mean_auc_mat, file.path(tiers$tables, paste0(out_prefix, "_mean_auc.csv")))
  write_matrix_csv(rss_mat, file.path(tiers$tables, paste0(out_prefix, "_rss.csv")))

  differential_df <- make_long_regulon_table(mean_auc_mat, rss_mat, label_order)
  utils::write.csv(
    differential_df,
    file.path(tiers$tables, paste0(out_prefix, "_differential_regulon_activity.csv")),
    row.names = FALSE
  )

  top_regulons <- unique(unlist(lapply(label_order, function(label) {
    vals <- sort(rss_mat[, label], decreasing = TRUE)
    names(vals)[seq_len(min(top_regulons_per_label, length(vals)))]
  })))
  top_regulons <- top_regulons[!is.na(top_regulons)]
  if (length(top_regulons) == 0) stop("No regulons available for plotting in ", out_prefix, ". Ask the user what to do.")

  plot_rss_mat <- rss_mat[top_regulons, label_order, drop = FALSE]
  plot_rss_scaled <- t(scale(t(plot_rss_mat)))
  plot_rss_scaled[!is.finite(plot_rss_scaled)] <- 0
  rownames(plot_rss_scaled) <- format_regulon_name(rownames(plot_rss_scaled))

  annotation_df <- data.frame(Timepoint = factor(label_order, levels = label_order))
  rownames(annotation_df) <- label_order
  ha_cols <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_df,
    col = list(Timepoint = label_colours[label_order]),
    show_annotation_name = FALSE
  )
  rss_col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
  grDevices::pdf(file.path(tiers$figures, paste0(out_prefix, "_regulon_heatmap.pdf")), width = 12, height = 10, useDingbats = FALSE)
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(
      plot_rss_scaled,
      name = "Scaled RSS",
      col = rss_col_fun,
      top_annotation = ha_cols,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_column_dend = FALSE,
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 10),
      column_names_rot = 45,
      heatmap_legend_param = list(title = "Scaled RSS")
    ),
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  grid::grid.text(
    paste0("SCENIC regulon specificity: ", gsub("_", " ", out_prefix)),
    x = grid::unit(4, "mm"),
    y = grid::unit(1, "npc") - grid::unit(4, "mm"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  grDevices::dev.off()

  ####################
  # Top-20 per-timepoint RSS-specific heatmap with fixed timepoint order
  ####################
  top20_specific_df <- dplyr::bind_rows(lapply(label_order, function(label) {
    vals <- sort(rss_mat[, label], decreasing = TRUE)
    keep <- names(vals)[seq_len(min(20, length(vals)))]
    data.frame(
      regulon = keep,
      regulon_label = format_regulon_name(keep),
      selected_timepoint = label,
      selected_rss = as.numeric(vals[keep]),
      stringsAsFactors = FALSE
    )
  })) |>
    dplyr::mutate(selected_timepoint = factor(selected_timepoint, levels = label_order)) |>
    dplyr::arrange(selected_timepoint, dplyr::desc(selected_rss), regulon_label)
  utils::write.csv(
    top20_specific_df,
    file.path(tiers$tables, paste0(out_prefix, "_top20_specific_regulons_by_timepoint.csv")),
    row.names = FALSE
  )
  top20_rss_mat <- rss_mat[top20_specific_df$regulon, label_order, drop = FALSE]
  write_matrix_csv(top20_rss_mat, file.path(tiers$tables, paste0(out_prefix, "_top20_specific_regulon_rss_matrix.csv")))
  top20_rss_scaled <- t(scale(t(top20_rss_mat)))
  top20_rss_scaled[!is.finite(top20_rss_scaled)] <- 0
  rownames(top20_rss_scaled) <- paste0(top20_specific_df$regulon_label, " [", top20_specific_df$selected_timepoint, "]")
  row_split <- factor(as.character(top20_specific_df$selected_timepoint), levels = label_order)
  grDevices::pdf(file.path(tiers$figures, paste0(out_prefix, "_top20_specific_regulon_heatmap.pdf")), width = 12, height = 16, useDingbats = FALSE)
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(
      top20_rss_scaled,
      name = "Scaled RSS",
      col = rss_col_fun,
      top_annotation = ha_cols,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_column_dend = FALSE,
      row_split = row_split,
      row_title_rot = 0,
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = 7),
      column_names_gp = grid::gpar(fontsize = 10),
      column_names_rot = 45,
      heatmap_legend_param = list(title = "Scaled RSS")
    ),
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  grid::grid.text(
    "Top 20 SCENIC regulons per timepoint ordered by RSS specificity",
    x = grid::unit(4, "mm"),
    y = grid::unit(1, "npc") - grid::unit(4, "mm"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  grDevices::dev.off()


  ####################
  # Balanced-count top-RSS and RSS-gap heatmaps
  ####################
  balanced_counts <- table(label_map)[label_order]
  balanced_n <- min(balanced_counts)
  set.seed(1090)
  balanced_cells <- unlist(lapply(label_order, function(label) {
    cells <- names(label_map)[label_map == label]
    if (length(cells) > balanced_n) cells <- sample(cells, balanced_n)
    cells
  }), use.names = FALSE)
  balanced_label_map <- label_map[balanced_cells]
  balanced_auc_mat <- auc_mat[, names(balanced_label_map), drop = FALSE]
  balanced_mean_auc_mat <- as.matrix(sapply(label_order, function(label) {
    cells <- names(balanced_label_map)[balanced_label_map == label]
    rowMeans(balanced_auc_mat[, cells, drop = FALSE], na.rm = TRUE)
  }))
  colnames(balanced_mean_auc_mat) <- label_order
  balanced_rss_mat <- tryCatch(
    SCENIC::calcRSS(AUC = balanced_auc_mat, cellAnnotation = balanced_label_map, cellTypes = label_order),
    error = function(e) NULL
  )
  if (is.null(balanced_rss_mat)) balanced_rss_mat <- balanced_mean_auc_mat
  balanced_rss_mat <- as.matrix(balanced_rss_mat)[, label_order, drop = FALSE]
  saveRDS(balanced_mean_auc_mat, file.path(tiers$intermediate, paste0(out_prefix, "_balanced", balanced_n, "_mean_auc.rds")))
  saveRDS(balanced_rss_mat, file.path(tiers$intermediate, paste0(out_prefix, "_balanced", balanced_n, "_rss.rds")))
  write_matrix_csv(balanced_mean_auc_mat, file.path(tiers$tables, paste0(out_prefix, "_balanced", balanced_n, "_mean_auc.csv")))
  write_matrix_csv(balanced_rss_mat, file.path(tiers$tables, paste0(out_prefix, "_balanced", balanced_n, "_rss.csv")))
  utils::write.csv(
    data.frame(timepoint = label_order, n_cells = as.integer(table(balanced_label_map)[label_order]), stringsAsFactors = FALSE),
    file.path(tiers$tables, paste0(out_prefix, "_balanced", balanced_n, "_selected_cells_summary.csv")),
    row.names = FALSE
  )

  balanced_top20_df <- dplyr::bind_rows(lapply(label_order, function(label) {
    vals <- sort(balanced_rss_mat[, label], decreasing = TRUE)
    keep <- names(vals)[seq_len(min(20, length(vals)))]
    data.frame(
      regulon = keep,
      regulon_label = format_regulon_name(keep),
      selected_timepoint = label,
      selected_rss = as.numeric(vals[keep]),
      stringsAsFactors = FALSE
    )
  })) |>
    dplyr::mutate(selected_timepoint = factor(selected_timepoint, levels = label_order)) |>
    dplyr::arrange(selected_timepoint, dplyr::desc(selected_rss), regulon_label)
  utils::write.csv(
    balanced_top20_df,
    file.path(tiers$tables, paste0(out_prefix, "_balanced", balanced_n, "_top20_specific_regulons_by_timepoint.csv")),
    row.names = FALSE
  )
  balanced_top20_mat <- balanced_rss_mat[balanced_top20_df$regulon, label_order, drop = FALSE]
  write_matrix_csv(balanced_top20_mat, file.path(tiers$tables, paste0(out_prefix, "_balanced", balanced_n, "_top20_specific_regulon_rss_matrix.csv")))
  balanced_top20_scaled <- t(scale(t(balanced_top20_mat)))
  balanced_top20_scaled[!is.finite(balanced_top20_scaled)] <- 0
  rownames(balanced_top20_scaled) <- paste0(balanced_top20_df$regulon_label, " [", balanced_top20_df$selected_timepoint, "]")
  balanced_row_split <- factor(as.character(balanced_top20_df$selected_timepoint), levels = label_order)
  grDevices::pdf(file.path(tiers$figures, paste0(out_prefix, "_balanced", balanced_n, "_top20_specific_regulon_heatmap.pdf")), width = 12, height = 16, useDingbats = FALSE)
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(
      balanced_top20_scaled,
      name = "Scaled RSS",
      col = rss_col_fun,
      top_annotation = ha_cols,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_column_dend = FALSE,
      row_split = balanced_row_split,
      row_title_rot = 0,
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = 7),
      column_names_gp = grid::gpar(fontsize = 10),
      column_names_rot = 45,
      heatmap_legend_param = list(title = "Scaled RSS")
    ),
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  grid::grid.text(
    paste0("Balanced top 20 SCENIC regulons per timepoint (n=", balanced_n, ")"),
    x = grid::unit(4, "mm"),
    y = grid::unit(1, "npc") - grid::unit(4, "mm"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  grDevices::dev.off()

  balanced_gap_df <- dplyr::bind_rows(lapply(label_order, function(label) {
    rows <- lapply(rownames(balanced_rss_mat), function(regulon) {
      vals <- balanced_rss_mat[regulon, label_order]
      ordered_vals <- sort(vals, decreasing = TRUE)
      best_label <- names(ordered_vals)[1]
      if (!identical(best_label, label)) return(NULL)
      data.frame(
        regulon = regulon,
        regulon_label = format_regulon_name(regulon),
        selected_timepoint = label,
        selected_rss = as.numeric(ordered_vals[1]),
        next_best_timepoint = names(ordered_vals)[2],
        next_best_rss = as.numeric(ordered_vals[2]),
        rss_gap = as.numeric(ordered_vals[1] - ordered_vals[2]),
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(rows) |>
      dplyr::arrange(dplyr::desc(rss_gap), dplyr::desc(selected_rss), regulon_label) |>
      dplyr::slice_head(n = 20)
  })) |>
    dplyr::mutate(selected_timepoint = factor(selected_timepoint, levels = label_order)) |>
    dplyr::arrange(selected_timepoint, dplyr::desc(rss_gap), dplyr::desc(selected_rss), regulon_label)
  utils::write.csv(
    balanced_gap_df,
    file.path(tiers$tables, paste0(out_prefix, "_balanced", balanced_n, "_top20_gap_regulons_by_timepoint.csv")),
    row.names = FALSE
  )
  balanced_gap_mat <- balanced_rss_mat[balanced_gap_df$regulon, label_order, drop = FALSE]
  write_matrix_csv(balanced_gap_mat, file.path(tiers$tables, paste0(out_prefix, "_balanced", balanced_n, "_top20_gap_regulon_rss_matrix.csv")))
  balanced_gap_scaled <- t(scale(t(balanced_gap_mat)))
  balanced_gap_scaled[!is.finite(balanced_gap_scaled)] <- 0
  rownames(balanced_gap_scaled) <- paste0(balanced_gap_df$regulon_label, " [", balanced_gap_df$selected_timepoint, "]")
  balanced_gap_row_split <- factor(as.character(balanced_gap_df$selected_timepoint), levels = label_order)
  grDevices::pdf(file.path(tiers$figures, paste0(out_prefix, "_balanced", balanced_n, "_top20_gap_regulon_heatmap.pdf")), width = 12, height = 16, useDingbats = FALSE)
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(
      balanced_gap_scaled,
      name = "Scaled RSS",
      col = rss_col_fun,
      top_annotation = ha_cols,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_column_dend = FALSE,
      row_split = balanced_gap_row_split,
      row_title_rot = 0,
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = 7),
      column_names_gp = grid::gpar(fontsize = 10),
      column_names_rot = 45,
      heatmap_legend_param = list(title = "Scaled RSS")
    ),
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  grid::grid.text(
    paste0("Balanced top 20 regulons by RSS gap vs next timepoint (n=", balanced_n, ")"),
    x = grid::unit(4, "mm"),
    y = grid::unit(1, "npc") - grid::unit(4, "mm"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  grDevices::dev.off()
  ####################


  edge_df <- dplyr::bind_rows(lapply(label_order, function(label) {
    vals <- sort(rss_mat[, label], decreasing = TRUE)
    keep <- names(vals)[seq_len(min(top_regulons_per_label, length(vals)))]
    data.frame(
      regulon = keep,
      regulon_label = format_regulon_name(keep),
      timepoint = label,
      weight = as.numeric(vals[keep]),
      stringsAsFactors = FALSE
    )
  })) |>
    dplyr::distinct(regulon_label, timepoint, .keep_all = TRUE) |>
    dplyr::filter(is.finite(weight), weight > 0)
  utils::write.csv(edge_df, file.path(tiers$tables, paste0(out_prefix, "_network_edges.csv")), row.names = FALSE)

  if (nrow(edge_df) > 0) {
    node_df <- data.frame(name = unique(c(edge_df$regulon_label, edge_df$timepoint)), stringsAsFactors = FALSE) |>
      dplyr::mutate(
        node_type = ifelse(name %in% label_order, "Timepoint", "Regulon"),
        timepoint = ifelse(node_type == "Timepoint", name, "Regulon")
      )
    network_graph <- tidygraph::tbl_graph(
      nodes = node_df,
      edges = edge_df |> dplyr::transmute(from = regulon_label, to = timepoint, weight = weight),
      directed = FALSE
    )
    network_fill <- c(label_colours[label_order], Regulon = "grey35")
    grDevices::pdf(file.path(tiers$figures, paste0(out_prefix, "_network.pdf")), width = 15, height = 10, useDingbats = FALSE)
    print(
      ggraph::ggraph(network_graph, layout = "stress") +
        ggraph::geom_edge_link(ggplot2::aes(width = weight, alpha = weight), colour = "grey70") +
        ggraph::scale_edge_width(range = c(0.4, 2.2)) +
        ggraph::scale_edge_alpha(range = c(0.3, 0.9)) +
        ggraph::geom_node_point(ggplot2::aes(fill = timepoint, shape = node_type), size = 5, colour = "black", stroke = 0.3) +
        ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 3, max.overlaps = 50) +
        ggplot2::scale_shape_manual(values = c(Timepoint = 21, Regulon = 22)) +
        ggplot2::scale_fill_manual(values = network_fill, drop = FALSE) +
        ggplot2::theme_void(base_size = 12) +
        ggplot2::labs(title = paste0("SCENIC timepoint regulatory network: ", gsub("_", " ", out_prefix))) +
        ggplot2::guides(edge_width = "none", edge_alpha = "none", shape = "none", fill = "none")
    )
    grDevices::dev.off()
  }

  regulon_target_df <- dplyr::bind_rows(lapply(seq_len(nrow(edge_df)), function(i) {
    reg_name <- edge_df$regulon[i]
    reg_targets <- extract_regulon_targets(regulons[[reg_name]])
    data.frame(
      timepoint = edge_df$timepoint[i],
      regulon = reg_name,
      regulon_label = edge_df$regulon_label[i],
      rss_weight = edge_df$weight[i],
      n_targets = length(reg_targets),
      targets_preview = paste(head(reg_targets, 40), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(regulon_target_df, file.path(tiers$tables, paste0(out_prefix, "_regulon_targets.csv")), row.names = FALSE)

  list(
    out_prefix = out_prefix,
    run_dir = run_dir,
    label_order = label_order,
    n_cells = ncol(counts_mat),
    n_genes_kept = length(genes_kept),
    n_regulons = nrow(auc_mat),
    top_regulons = top_regulons,
    mean_auc_mat = mean_auc_mat,
    rss_mat = rss_mat
  )
}

####################
# Main workflow
####################
arg_list <- parse_args(commandArgs(trailingOnly = TRUE))
run_mode <- arg_list[["run_mode"]] %||% Sys.getenv("RUN_MODE", unset = "combined")
sample_id <- arg_list[["sample_id"]] %||% Sys.getenv("SAMPLE_ID", unset = "")
n_cores <- as.integer(arg_list[["n_cores"]] %||% Sys.getenv("N_CORES", unset = "12"))
cells_per_timepoint <- as.integer(arg_list[["cells_per_timepoint"]] %||% Sys.getenv("CELLS_PER_TIMEPOINT", unset = "0"))
genie3_nparts <- as.integer(arg_list[["genie3_nparts"]] %||% Sys.getenv("GENIE3_NPARTS", unset = "100"))
genie3_resume <- tolower(arg_list[["genie3_resume"]] %||% Sys.getenv("GENIE3_RESUME", unset = "true")) %in% c("true", "1", "yes", "y")
top_regulons_per_label <- as.integer(arg_list[["top_regulons_per_label"]] %||% "12")
random_seed <- as.integer(arg_list[["random_seed"]] %||% "1090")
db_dir <- arg_list[["db_dir"]] %||% Sys.getenv(
  "SCENIC_DB_DIR",
  unset = "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/cistarget_databases_rcistarget_mc9nr"
)
db_dir <- normalizePath(db_dir, winslash = "/", mustWork = FALSE)

if (!run_mode %in% c("preflight", "combined")) {
  stop("run_mode must be one of preflight, combined. Per-timepoint SCENIC runs are intentionally disabled because they create non-comparable regulon dictionaries.")
}

paths <- parse_paths(root_dir)
out_base <- file.path(paths$parse_outs, "cell_states", "timepoint_scenic")
tiers <- parse_output_tiers(out_base, create = TRUE)

input_files <- c(
  file.path(paths$parse_outs, "Auto_parse_merged.rds"),
  db_dir
)
expected_outputs <- c(
  file.path(tiers$tables, "timepoint_scenic_database_confirmation.csv"),
  file.path(tiers$tables, "combined_timepoint_differential_regulon_activity.csv")
)
run <- parse_start_run(
  script_name = paste0(
    "parse_timepoint_scenic_regulons_",
    run_mode,
    if (nzchar(sample_id)) paste0("_", sample_id) else ""
  ),
  parameters = list(
    run_mode = run_mode,
    sample_id = sample_id,
    n_cores = n_cores,
    genie3_nparts = genie3_nparts,
    genie3_resume = genie3_resume,
    cells_per_timepoint = cells_per_timepoint,
    top_regulons_per_label = top_regulons_per_label,
    random_seed = random_seed,
    db_dir = db_dir
  ),
  input_files = input_files,
  output_files = expected_outputs,
  reused_cached = TRUE
)
run_status <- "failed"
run_notes <- character()
on.exit(parse_finish_run(run, status = run_status, output_files = expected_outputs, reused_cached = TRUE, notes = run_notes), add = TRUE)

merged_path <- check_file_or_stop(file.path(paths$parse_outs, "Auto_parse_merged.rds"), "merged Parse Seurat object")
db_files <- detect_db_files(db_dir)
utils::write.csv(
  data.frame(
    selected_database_dir = db_dir,
    selected_database_files = paste(db_files, collapse = " | "),
    reason = "PDO SCENIC workflow defaults to cistarget_databases_rcistarget_mc9nr and these are the hg38 refseq-r80 mc9nr RcisTarget databases.",
    stringsAsFactors = FALSE
  ),
  file.path(tiers$tables, "timepoint_scenic_database_confirmation.csv"),
  row.names = FALSE
)

patch_scenic_annotation_lookup()

set.seed(random_seed)
message("Loading merged Parse Seurat object: ", merged_path)
merged_obj <- readRDS(merged_path)
if (!"orig.ident" %in% colnames(merged_obj@meta.data)) stop("orig.ident is missing from merged object metadata. Ask the user what to do.")
missing_samples <- setdiff(parse_samples, unique(as.character(merged_obj$orig.ident)))
if (length(missing_samples) > 0) stop("Missing required Parse timepoint(s): ", paste(missing_samples, collapse = ", "), ". Ask the user what to do.")

all_timepoint_cells <- unlist(lapply(parse_samples, function(sample_name) {
  cells <- colnames(merged_obj)[as.character(merged_obj$orig.ident) == sample_name]
  if (length(cells) == 0) stop("No cells found for ", sample_name, ". Ask the user what to do.")
  if (cells_per_timepoint > 0 && length(cells) > cells_per_timepoint) {
    cells <- sample(cells, cells_per_timepoint)
  }
  stats::setNames(cells, rep(sample_name, length(cells)))
}), use.names = FALSE)
selected_meta <- data.frame(
  cell = all_timepoint_cells,
  timepoint = as.character(merged_obj$orig.ident[all_timepoint_cells]),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(timepoint = factor(timepoint, levels = parse_samples)) |>
  dplyr::arrange(timepoint, cell)
utils::write.csv(selected_meta, file.path(tiers$tables, "timepoint_scenic_selected_cells.csv"), row.names = FALSE)

if (run_mode == "preflight") {
  preflight_counts <- extract_counts_for_cells(merged_obj, selected_meta$cell, assay = "RNA")
  preflight_summary <- selected_meta |>
    dplyr::count(timepoint, name = "n_cells") |>
    dplyr::mutate(
      n_genes_in_counts = nrow(preflight_counts),
      db_dir = db_dir,
      db_files = paste(db_files, collapse = " | ")
    )
  utils::write.csv(preflight_summary, file.path(tiers$tables, "timepoint_scenic_preflight_summary.csv"), row.names = FALSE)
  expected_outputs <- c(
    file.path(tiers$tables, "timepoint_scenic_database_confirmation.csv"),
    file.path(tiers$tables, "timepoint_scenic_selected_cells.csv"),
    file.path(tiers$tables, "timepoint_scenic_preflight_summary.csv")
  )
  run_notes <- c(
    "Preflight completed without starting SCENIC.",
    "The selected cisTarget directory is cistarget_databases_rcistarget_mc9nr.",
    "State and metaprogram assignments are not loaded or used."
  )
  run_status <- "success"
  message("Preflight checks passed. No SCENIC run was started.")
  parse_finish_run(run, status = run_status, output_files = expected_outputs, reused_cached = TRUE, notes = run_notes)
  quit(save = "no")
}

counts_all <- extract_counts_for_cells(merged_obj, selected_meta$cell, assay = "RNA")
label_map_all <- stats::setNames(as.character(selected_meta$timepoint), selected_meta$cell)

run_results <- list()
if (run_mode == "combined") {
  message("Running combined six-timepoint SCENIC analysis.")
  run_results[["combined"]] <- run_scenic_for_cells(
    counts_mat = counts_all,
    label_map = label_map_all,
    label_order = parse_samples,
    label_colours = parse_sample_colours,
    run_dir = file.path(tiers$intermediate, "scenic_combined"),
    out_prefix = "combined_timepoint",
    dataset_title = "Parse_timepoint_scenic_combined",
    db_dir = db_dir,
    db_files = db_files,
    n_cores = n_cores,
    genie3_nparts = genie3_nparts,
    genie3_resume = genie3_resume,
    top_regulons_per_label = top_regulons_per_label,
    tiers = tiers
  )
}

if (length(run_results) > 0) {
  summary_df <- dplyr::bind_rows(lapply(names(run_results), function(name) {
    res <- run_results[[name]]
    data.frame(
      run = name,
      run_dir = res$run_dir,
      labels = paste(res$label_order, collapse = " | "),
      n_cells = res$n_cells,
      n_genes_kept = res$n_genes_kept,
      n_regulons = res$n_regulons,
      stringsAsFactors = FALSE
    )
  }))
  run_summary_table <- file.path(tiers$tables, "timepoint_scenic_run_summary.csv")
  utils::write.csv(summary_df, run_summary_table, row.names = FALSE)
}

expected_outputs <- c(
  file.path(tiers$tables, "timepoint_scenic_database_confirmation.csv"),
  file.path(tiers$tables, "timepoint_scenic_selected_cells.csv"),
  file.path(tiers$tables, "timepoint_scenic_run_summary.csv")
)
if (run_mode == "combined") {
  expected_outputs <- c(
    expected_outputs,
    file.path(tiers$tables, "combined_timepoint_differential_regulon_activity.csv"),
    file.path(tiers$tables, "combined_timepoint_top20_specific_regulons_by_timepoint.csv"),
    file.path(tiers$tables, "combined_timepoint_top20_specific_regulon_rss_matrix.csv"),
    file.path(tiers$tables, "combined_timepoint_balanced2600_rss.csv"),
    file.path(tiers$tables, "combined_timepoint_balanced2600_top20_specific_regulons_by_timepoint.csv"),
    file.path(tiers$tables, "combined_timepoint_balanced2600_top20_gap_regulons_by_timepoint.csv"),
    file.path(tiers$figures, "combined_timepoint_regulon_heatmap.pdf"),
    file.path(tiers$figures, "combined_timepoint_top20_specific_regulon_heatmap.pdf"),
    file.path(tiers$figures, "combined_timepoint_balanced2600_top20_specific_regulon_heatmap.pdf"),
    file.path(tiers$figures, "combined_timepoint_balanced2600_top20_gap_regulon_heatmap.pdf"),
    file.path(tiers$figures, "combined_timepoint_network.pdf")
  )
}
run_notes <- c(
  "The selected cisTarget directory is cistarget_databases_rcistarget_mc9nr.",
  "State and metaprogram assignments are not loaded or used.",
    "SCENIC is run only in combined mode so all timepoints share a single regulon dictionary for clean AUC/RSS comparison.",
    paste0("GENIE3 resumePreviousRun is ", genie3_resume, " with nParts=", genie3_nparts, "."),
    if (cells_per_timepoint > 0) paste0("Cells were capped at ", cells_per_timepoint, " per timepoint.") else "All available cells from the six Parse timepoints were used."
  )
run_status <- "success"
message("Saved Parse timepoint SCENIC outputs under: ", out_base)
