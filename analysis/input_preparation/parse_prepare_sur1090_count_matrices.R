#!/usr/bin/env Rscript
####################
# parse_prepare_sur1090_count_matrices.R
#
# Description:
#   Converts SUR1090 treated/untreated PDO CSV count matrices into the Parse
#   pipeline DGE_filtered matrix/metadata layout expected by the QC workflow.
#
# Inputs:
#   /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/00_counts_matrix_all/SUR1090_Treated_PDO.csv
#   /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/00_counts_matrix_all/SUR1090_Untreated_PDO.csv
#
# Outputs:
#   parse_outs/input/output_combined/<sample>/DGE_filtered/count_matrix.mtx.gz
#   parse_outs/input/output_combined/<sample>/DGE_filtered/all_genes.csv.gz
#   parse_outs/input/output_combined/<sample>/DGE_filtered/cell_metadata.csv.gz
#   parse_outs/logs/run_summaries/parse_prepare_sur1090_count_matrices_*.txt
#
# Methodology:
#   analysis/methodology/input_preparation/prepare_sur1090_count_matrices_methodology.md
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_prepare_sur1090_count_matrices",
  input_files = "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/00_counts_matrix_all/SUR1090_*_PDO.csv",
  output_files = "parse_outs/input/output_combined/<sample>/DGE_filtered/{count_matrix.mtx.gz,all_genes.csv.gz,cell_metadata.csv.gz}"
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

suppressPackageStartupMessages({
  library("Matrix")
  library("data.table")
})

# Paths
input_dir <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/00_counts_matrix_all"
output_base <- file.path(parse_project_root(), "parse_outs/input/output_combined")

# Samples to process
samples <- c("SUR1090_Treated_PDO", "SUR1090_Untreated_PDO")

for (sample_name in samples) {
  output_name <- sub("_PDO$", "", sample_name)
  output_dir <- file.path(output_base, output_name, "DGE_filtered")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  input_file <- file.path(input_dir, paste0(sample_name, ".csv"))
  message("Processing: ", sample_name)
  
  # Read CSV - first row is header with cell barcodes, first column is gene names
  count_mat <- fread(input_file, header = TRUE)
  
  # First column contains gene names
  gene_names <- count_mat[[1]]
  count_mat[[1]] <- NULL
  
  # Convert to matrix (genes x cells)
  counts <- as.matrix(count_mat)
  rownames(counts) <- gene_names
  
  # Cell IDs are in the column names
  cell_ids <- colnames(count_mat)
  n_cells <- ncol(counts)
  n_genes <- nrow(counts)
  
  message("  Matrix: ", n_genes, " genes x ", n_cells, " cells")
  
  # TRANSPOSE: need cells as rows (first index), genes as columns
  counts_t <- t(counts)  # now cells x genes
  counts_sparse <- as(counts_t, "dgCMatrix")
  
  message("  Sparse matrix: ", nrow(counts_sparse), " cells x ", ncol(counts_sparse), " genes")
  
  # Save count matrix in MTX format (cells x genes)
  writeMM(counts_sparse, file.path(output_dir, "count_matrix.mtx.gz"))
  
  # Gene metadata
  gene_df <- data.frame(
    gene_name = gene_names,
    gene_id = NA,
    stringsAsFactors = FALSE
  )
  fwrite(gene_df, file.path(output_dir, "all_genes.csv.gz"))
  
  # Cell metadata - calculate using base R (for transposed matrix, use rowSums)
  ncount_rna <- as.numeric(colSums(counts))  # for genes x cells matrix, use colSums
  nfeature_rna <- as.numeric(colSums(counts > 0))
  
  meta_df <- data.frame(
    bc_wells = cell_ids,
    nCount_RNA = ncount_rna,
    nFeature_RNA = nfeature_rna,
    percent.mt = NA,
    stringsAsFactors = FALSE
  )
  fwrite(meta_df, file.path(output_dir, "cell_metadata.csv.gz"))
  
  message("  Saved ", n_cells, " cells x ", n_genes, " genes")
}

script_run_status <- "success"
message("Done!")
