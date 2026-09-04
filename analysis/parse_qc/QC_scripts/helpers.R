suppressPackageStartupMessages({
  library("Matrix")
  library("data.table")
})

sample_order <- c("PDO", "T0", "T1", "T2", "T4", "R4", "eR4")

max_mt <- 15
min_ngenes <- 2500
max_ngenes <- 13000
min_hk_expr <- 3
n_pcs <- 30

housekeeping_genes <- c(
  "ACTB", "GAPDH", "RPS11", "RPS13", "RPS14", "RPS15", "RPS16", "RPS18",
  "RPS19", "RPS20", "RPL10", "RPL13", "RPL15", "RPL18"
)

marker_panel <- list(
  b.cell = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
  dendritic = c("CLEC10A", "CCR7", "CD86"),
  endothelial = c("ENG", "CLEC14A", "CLDN5", "VWF", "CDH5"),
  epithelial = c("KRT7", "MUC1", "KRT19", "EPCAM"),
  fibroblast = c("COL3A1", "COL1A2", "LUM", "COL1A1", "COL6A3", "DCN"),
  macrophage = c("CSF1R", "TYROBP", "CD14", "CD163", "AIF1", "CD68"),
  mast = c("MS4A2", "CPA3", "TPSB2", "TPSAB1"),
  nk.cell = c("GNLY", "NKG7", "PRF1", "GZMB", "KLRB1"),
  plasma = c("MZB1", "JCHAIN", "DERL3"),
  t.cell = c("CD3E", "CD3D", "CD2", "CD3G"),
  housekeeping = housekeeping_genes
)

marker_genes <- unique(unlist(marker_panel, use.names = FALSE))

read_parse_sample <- function(sample_dir, sample_name) {
  meta <- fread(file.path(sample_dir, "cell_metadata.csv.gz"))
  genes <- fread(file.path(sample_dir, "all_genes.csv.gz"))
  counts <- t(readMM(gzfile(file.path(sample_dir, "count_matrix.mtx.gz"))))
  rownames(counts) <- make.unique(genes$gene_name)
  colnames(counts) <- make.unique(meta$bc_wells)

  meta <- as.data.frame(meta)
  rownames(meta) <- colnames(counts)
  meta$orig.ident <- sample_name
  meta$sample <- sample_name

  list(
    counts = counts,
    meta = meta,
    genes = genes
  )
}

compute_hk_mean <- function(counts) {
  totals <- Matrix::colSums(counts)
  hk_present <- intersect(housekeeping_genes, rownames(counts))
  if (length(hk_present) == 0) {
    return(rep(NA_real_, ncol(counts)))
  }
  hk_counts <- counts[hk_present, , drop = FALSE]
  hk_cpm <- t(t(hk_counts) / totals) * 1e6
  hk_cpm[is.na(hk_cpm)] <- 0
  Matrix::colMeans(log2(as.matrix(hk_cpm) / 10 + 1))
}

get_marker_expr <- function(counts) {
  keep_genes <- intersect(marker_genes, rownames(counts))
  expr <- matrix(
    0,
    nrow = length(marker_genes),
    ncol = ncol(counts),
    dimnames = list(marker_genes, colnames(counts))
  )
  if (length(keep_genes) == 0) {
    return(expr)
  }
  totals <- Matrix::colSums(counts)
  marker_counts <- counts[keep_genes, , drop = FALSE]
  marker_cpm <- t(t(marker_counts) / totals) * 1e6
  marker_cpm[is.na(marker_cpm)] <- 0
  expr[keep_genes, ] <- log2(as.matrix(marker_cpm) / 10 + 1)
  expr
}

get_doublet_rate <- function(ncells) {
  if (ncells <= 1000) {
    return(0.008)
  }
  if (ncells <= 5000) {
    return(0.04)
  }
  if (ncells <= 10000) {
    return(0.08)
  }
  if (ncells <= 20000) {
    return(0.16)
  }
  if (ncells <= 30000) {
    return(0.24)
  }
  0.24
}

