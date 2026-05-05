#!/usr/bin/env Rscript
# Script to convert SUR1090 CSV files to pipeline format

suppressPackageStartupMessages({
  library("Matrix")
  library("data.table")
})

# Paths
input_dir <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/00_counts_matrix_all"
output_base <- "/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline/parse_outs/input/output_combined"

# Samples to process
samples <- c("SUR1090_Treated_PDO", "SUR1090_Untreated_PDO")

for (sample_name in samples) {
  output_dir <- file.path(output_base, sample_name, "DGE_filtered")
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

message("Done!")