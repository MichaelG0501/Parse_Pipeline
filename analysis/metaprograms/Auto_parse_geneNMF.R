####################
# Auto_parse_geneNMF.R
# GeneNMF metaprogram analysis for Parse samples only.
# Excludes PDO and SUR1090 samples, and keeps T0/T1/T2/T4/R4/eR4.
####################

suppressPackageStartupMessages({
  library(GeneNMF)
  library(RColorBrewer)
  library(msigdbr)
  library(fgsea)
  library(UCell)
  library(Seurat)
  library(SeuratObject)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline/analysis/metaprograms/Auto_parse_geneNMF.R")
}

script_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(script_dir, "..", ".."))
out_dir <- file.path(project_dir, "parse_outs")
setwd(out_dir)

parse_mp_dir <- file.path("Auto_parse_metaprograms")
mp_results_dir <- file.path(parse_mp_dir, "Metaprogrammes_Results")
dir.create(parse_mp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mp_results_dir, recursive = TRUE, showWarnings = FALSE)

requested_samples <- c("T0", "T1", "T2", "T3", "T4", "R4", "eR4")
parse_samples <- c("T0", "T1", "T2", "T4", "R4", "eR4")
excluded_samples <- c("PDO", "SUR1090_Treated", "SUR1090_Untreated")

message("Requested samples: ", paste(requested_samples, collapse = ", "))
message("Using Parse samples: ", paste(parse_samples, collapse = ", "))
message("Excluding samples: ", paste(excluded_samples, collapse = ", "))
missing_requested <- setdiff(requested_samples, parse_samples)
if (length(missing_requested) > 0) {
  message("Requested samples absent from current Parse QC outputs: ", paste(missing_requested, collapse = ", "))
}

####################
# Load input data
####################

parse_list_path <- file.path(parse_mp_dir, "Auto_parse_list_geneNMF.rds")
parse_merged_path <- file.path(parse_mp_dir, "Auto_parse_merged_geneNMF.rds")

sample_files <- setNames(
  file.path("by_samples", parse_samples, paste0("Auto_", parse_samples, "_final.rds")),
  parse_samples
)
missing_files <- names(sample_files)[!file.exists(sample_files)]
if (length(missing_files) > 0) {
  stop("Missing final sample RDS files: ", paste(missing_files, collapse = ", "))
}

if (file.exists(parse_list_path)) {
  message("Loading existing ", parse_list_path)
  parse.list <- readRDS(parse_list_path)
} else {
  message("Building Parse GeneNMF sample list...")
  parse.list <- lapply(sample_files, readRDS)
  names(parse.list) <- parse_samples
  parse.list <- parse.list[parse_samples]
  saveRDS(parse.list, parse_list_path, compress = FALSE)
}

if (file.exists(parse_merged_path)) {
  message("Loading existing ", parse_merged_path)
  parse_obj <- readRDS(parse_merged_path)
} else {
  message("Building merged Parse object for UCell scoring...")
  parse_obj <- merge(parse.list[[1]], y = parse.list[-1], add.cell.ids = parse_samples, project = "parse_geneNMF")
  parse_obj <- NormalizeData(parse_obj, verbose = FALSE)
  saveRDS(parse_obj, parse_merged_path, compress = FALSE)
}

####################
# Step 1: multiNMF
####################

geneNMF_out_path <- file.path(parse_mp_dir, "Auto_parse_geneNMF_outs.rds")
if (file.exists(geneNMF_out_path)) {
  message("Loading existing ", geneNMF_out_path)
  geneNMF.programs <- readRDS(geneNMF_out_path)
} else {
  message("Running multiNMF...")
  geneNMF.programs <- multiNMF(parse.list, assay = "RNA", k = 4:9, min.exp = 0.05)
  saveRDS(geneNMF.programs, file = geneNMF_out_path, compress = FALSE)
}

####################
# Step 2: Loop over nMP values.
# Range follows find_optimal_nmf.R so the optimum is data-driven.
####################

k_vals <- 4:35
for (k in k_vals) {
  rds_path <- file.path(mp_results_dir, paste0("Auto_parse_geneNMF_metaprograms_nMP_", k, ".rds"))
  png_path <- file.path(mp_results_dir, paste0("Auto_parse_metaprograms_heatmap_nMP_", k, ".png"))

  if (file.exists(rds_path) && file.exists(png_path)) {
    message(paste0("nMP = ", k, " already exists, skipping."))
    next
  }

  message(paste0("Running getMetaPrograms with nMP = ", k))
  geneNMF.metaprograms <- getMetaPrograms(
    geneNMF.programs,
    metric = "cosine",
    specificity.weight = 5,
    weight.explained = 0.5,
    nMP = k,
    min.confidence = 0.5
  )
  saveRDS(geneNMF.metaprograms, file = rds_path, compress = FALSE)

  n_colors <- min(k, 12)
  anno_colors <- brewer.pal(n = max(n_colors, 3), name = "Paired")
  anno_colors <- anno_colors[seq_len(length(geneNMF.metaprograms$metaprograms.genes))]
  names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

  png(png_path, width = 3000, height = 2500, res = 300)
  plotMetaPrograms(
    geneNMF.metaprograms,
    annotation_colors = anno_colors,
    similarity.cutoff = c(0, 1)
  )
  dev.off()
  message(paste0("Saved nMP = ", k))
}

####################
# Step 3: Select optimal nMP if available.
# Auto_parse_find_optimal_nmf.R writes Auto_parse_optimal_nMP.txt.
####################

optimal_txt <- file.path(parse_mp_dir, "Auto_parse_optimal_nMP.txt")
default_nMP <- if (file.exists(optimal_txt)) {
  as.integer(readLines(optimal_txt, warn = FALSE)[1])
} else {
  NA_integer_
}

if (is.na(default_nMP)) {
  message("Optimal nMP has not been selected yet. Run Auto_parse_find_optimal_nmf.R next.")
} else {
  rds_default <- file.path(mp_results_dir, paste0("Auto_parse_geneNMF_metaprograms_nMP_", default_nMP, ".rds"))
  if (!file.exists(rds_default)) {
    stop(paste0("Selected nMP file not found: ", rds_default))
  }

  message(paste0("Setting nMP = ", default_nMP, " as default."))
  geneNMF.metaprograms <- readRDS(rds_default)
  saveRDS(geneNMF.metaprograms, file = file.path(parse_mp_dir, "Auto_parse_MP_outs_default.rds"), compress = FALSE)

  anno_colors <- brewer.pal(n = min(default_nMP, 12), name = "Paired")
  anno_colors <- anno_colors[seq_len(length(geneNMF.metaprograms$metaprograms.genes))]
  names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

  png(file.path(parse_mp_dir, "Auto_parse_metaprograms_heatmap.png"), width = 3000, height = 2500, res = 300)
  plotMetaPrograms(
    geneNMF.metaprograms,
    annotation_colors = anno_colors,
    similarity.cutoff = c(0, 1)
  )
  dev.off()
}

####################
# Step 4: UCell scoring if optimal nMP exists.
####################

if (exists("geneNMF.metaprograms")) {
  bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
  bad_mp_names <- paste0("MP", bad_mps)

  coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
  names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
  low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]

  exclude_mps <- unique(c(bad_mp_names, low_coverage_mps))
  if (length(exclude_mps) > 0) {
    message(paste0("Excluding MPs: ", paste(exclude_mps, collapse = ", ")))
  }

  mp.genes <- geneNMF.metaprograms$metaprograms.genes
  mp.genes <- mp.genes[!names(mp.genes) %in% exclude_mps]

  final_path <- file.path(parse_mp_dir, "Auto_parse_final_geneNMF.rds")
  ucell_path <- file.path(parse_mp_dir, paste0("Auto_parse_UCell_scores_filtered_nMP", default_nMP, ".rds"))
  legacy_ucell_path <- file.path(parse_mp_dir, "Auto_parse_UCell_scores_filtered.rds")

  if (file.exists(final_path) && file.exists(ucell_path)) {
    message("UCell outputs already exist, skipping scoring.")
    parse_obj <- readRDS(final_path)
  } else if (file.exists(legacy_ucell_path)) {
    legacy_ucell_scores <- readRDS(legacy_ucell_path)
    if (all(names(mp.genes) %in% colnames(legacy_ucell_scores))) {
      message("Using matching legacy UCell score file for selected nMP.")
      ucell_scores <- legacy_ucell_scores[, names(mp.genes), drop = FALSE]
      saveRDS(ucell_scores, ucell_path, compress = FALSE)
    } else {
      message("Legacy UCell score file does not match selected nMP; rescoring.")
      if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
        parse_obj <- SeuratObject::JoinLayers(parse_obj, assay = "RNA")
      }
      parse_obj <- AddModuleScore_UCell(parse_obj, features = mp.genes, ncores = 4, name = "")
      ucell_scores <- parse_obj@meta.data[, names(mp.genes), drop = FALSE]
      saveRDS(ucell_scores, ucell_path, compress = FALSE)
      saveRDS(parse_obj, final_path, compress = FALSE)
    }
  } else {
    message("Running UCell scoring...")
    if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
      parse_obj <- SeuratObject::JoinLayers(parse_obj, assay = "RNA")
    }
    parse_obj <- AddModuleScore_UCell(parse_obj, features = mp.genes, ncores = 4, name = "")
    ucell_scores <- parse_obj@meta.data[, names(mp.genes), drop = FALSE]
    saveRDS(ucell_scores, ucell_path, compress = FALSE)
    saveRDS(parse_obj, final_path, compress = FALSE)
    message("Saved ", ucell_path)
  }

  png(file.path(parse_mp_dir, "Auto_parse_vln_origident.png"), width = 5000, height = 2500, res = 300)
  print(VlnPlot(parse_obj, features = names(mp.genes), group.by = "orig.ident", pt.size = 0, ncol = 5))
  dev.off()
}

message("Auto_parse_geneNMF.R completed successfully.")
