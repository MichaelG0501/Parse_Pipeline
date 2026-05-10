####################
# Shared helpers for Parse sample-level Monocle3 pseudotime workflows.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(igraph)
  library(patchwork)
  library(scales)
})

set.seed(12345)

get_parse_root_dir <- function() {
  env_root <- Sys.getenv("AUTO_PARSE_ROOT_DIR", unset = "")
  if (nzchar(env_root)) return(normalizePath(env_root, mustWork = FALSE))

  candidates <- c(
    "/rds/general/project/spatialtranscriptomics/ephemeral/Parse_Pipeline",
    "/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline"
  )
  candidates[file.exists(candidates)][1]
}

root_dir <- get_parse_root_dir()
setwd(root_dir)
qc_dir <- file.path(root_dir, "parse_outs")
pseudotime_dir <- file.path(qc_dir, "pseudotime_samples")
summary_dir <- file.path(qc_dir, "summary")
dir.create(pseudotime_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

####################
# Six-sample trajectory configuration. The previous all-sample run including
# PDO and SUR1090 is intentionally ignored; this default keeps only the Parse
# treatment-response samples requested for the rerun.
####################
analysis_samples <- strsplit(
  Sys.getenv("AUTO_PARSE_PSEUDOTIME_SAMPLES", unset = "T0,T1,T2,T4,R4,eR4"),
  ",",
  fixed = TRUE
)[[1]]
analysis_samples <- trimws(analysis_samples)
analysis_samples <- analysis_samples[nzchar(analysis_samples)]

sample_order_all <- analysis_samples

root_sample <- Sys.getenv("AUTO_PARSE_ROOT_SAMPLE", unset = "T0")
min_cells_per_sample <- as.integer(Sys.getenv("AUTO_PARSE_MIN_CELLS_PER_SAMPLE", unset = "30"))
sample_cols_all <- c(
  "T0" = "#0072B2",
  "T1" = "#E69F00",
  "T2" = "#009E73",
  "T4" = "#D55E00",
  "R4" = "#CC79A7",
  "eR4" = "#56B4E9"
)
missing_config_cols <- setdiff(sample_order_all, names(sample_cols_all))
if (length(missing_config_cols) > 0) {
  sample_cols_all <- c(
    sample_cols_all,
    setNames(hue_pal()(length(missing_config_cols)), missing_config_cols)
  )
}
sample_cols_all <- sample_cols_all[sample_order_all]

load_parse_object <- function() {
  obj_path <- file.path(qc_dir, "Auto_parse_merged.rds")
  if (!file.exists(obj_path)) stop("Missing Parse Seurat object: ", obj_path)
  readRDS(obj_path)
}

detect_sample_col <- function(seurat_obj) {
  candidates <- c("sample", "orig.ident", "Sample", "sample_id")
  sample_col <- candidates[candidates %in% colnames(seurat_obj@meta.data)][1]
  if (is.na(sample_col)) stop("No sample column found in Seurat metadata.")
  sample_col
}

ordered_samples <- function(sample_vec) {
  present <- unique(as.character(sample_vec))
  c(intersect(sample_order_all, present), sort(setdiff(present, sample_order_all)))
}

get_sample_colours <- function(samples) {
  samples <- as.character(samples)
  missing <- setdiff(samples, names(sample_cols_all))
  if (length(missing) > 0) {
    extra_cols <- setNames(hue_pal()(length(missing)), missing)
    return(c(sample_cols_all[intersect(names(sample_cols_all), samples)], extra_cols)[samples])
  }
  sample_cols_all[samples]
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

prepare_parse_sample_trajectory <- function(seurat_obj) {
  if (inherits(seurat_obj[["RNA"]], "Assay5")) {
    seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
  }

  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

  n_pcs <- min(30, ncol(seurat_obj) - 1)
  if (n_pcs < 2) stop("Too few cells for PCA.")

  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = n_pcs, verbose = FALSE)
  dims_use <- min(15, n_pcs)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims_use, verbose = FALSE)

  cds <- as.cell_data_set(seurat_obj)
  cds <- cluster_cells(cds, verbose = FALSE)
  cds <- learn_graph(cds, verbose = FALSE, use_partition = FALSE)
  cds
}

build_parse_trajectory <- function() {
  message("Loading Parse Seurat object ...")
  parse_obj <- load_parse_object()
  sample_col <- detect_sample_col(parse_obj)
  parse_obj$sample_label <- as.character(parse_obj@meta.data[[sample_col]])
  missing_samples <- setdiff(analysis_samples, unique(parse_obj$sample_label))
  if (length(missing_samples) > 0) {
    stop("Requested samples missing from Parse object: ", paste(missing_samples, collapse = ", "))
  }
  parse_obj <- subset(parse_obj, subset = sample_label %in% analysis_samples)

  sample_levels <- ordered_samples(parse_obj$sample_label)
  sample_counts <- table(factor(parse_obj$sample_label, levels = sample_levels))
  keep_samples <- names(sample_counts)[sample_counts >= min_cells_per_sample]
  keep_samples <- intersect(analysis_samples, keep_samples)
  dropped_samples <- setdiff(analysis_samples, keep_samples)
  if (length(dropped_samples) > 0) {
    message("Dropping samples below minimum cell threshold: ", paste(dropped_samples, collapse = ", "))
  }
  if (!root_sample %in% keep_samples) {
    stop("Root sample '", root_sample, "' is not available after filtering.")
  }
  parse_obj <- subset(parse_obj, subset = sample_label %in% keep_samples)
  parse_obj$sample_label <- factor(as.character(parse_obj$sample_label), levels = keep_samples)

  root_cells <- colnames(parse_obj)[as.character(parse_obj$sample_label) == root_sample]
  if (length(root_cells) == 0) {
    stop("No root cells found for root sample '", root_sample, "'.")
  }

  message("Preparing Monocle3 trajectory using sample labels; root sample: ", root_sample)
  cds <- prepare_parse_sample_trajectory(parse_obj)
  cds <- order_cells(cds, root_cells = intersect(root_cells, colnames(cds)))

  pt <- pseudotime(cds)
  pt[is.infinite(pt)] <- NA_real_

  list(
    seurat_obj = parse_obj,
    cds = cds,
    pseudotime = pt,
    sample_col = sample_col,
    sample_levels = keep_samples,
    root_cells = root_cells
  )
}

get_or_build_parse_trajectory <- function(rebuild = FALSE) {
  cds_path <- file.path(pseudotime_dir, "Auto_parse_sample_pseudotime_cds.rds")
  pt_path <- file.path(pseudotime_dir, "Auto_parse_sample_pseudotime.rds")
  meta_path <- file.path(pseudotime_dir, "Auto_parse_sample_pseudotime_metadata.csv")

  if (!rebuild && all(file.exists(cds_path, pt_path, meta_path))) {
    message("Loading existing Parse pseudotime assets from: ", pseudotime_dir)
    cds <- readRDS(cds_path)
    pt <- readRDS(pt_path)
    meta <- read_csv(meta_path, show_col_types = FALSE)
    sample_levels <- ordered_samples(meta$sample)
    root_cells <- meta$cell[meta$sample == root_sample]
    return(list(
      seurat_obj = NULL,
      cds = cds,
      pseudotime = pt,
      metadata = meta,
      sample_col = "sample",
      sample_levels = sample_levels,
      root_cells = root_cells
    ))
  }

  build_parse_trajectory()
}

extract_graph_structure <- function(cds) {
  graph_obj <- principal_graph(cds)[["UMAP"]]
  graph_coords <- cds@principal_graph_aux[["UMAP"]]$dp_mst

  edge_df <- igraph::as_data_frame(graph_obj, what = "edges") %>%
    mutate(
      x = graph_coords[1, from],
      y = graph_coords[2, from],
      xend = graph_coords[1, to],
      yend = graph_coords[2, to]
    )

  node_df <- data.frame(
    node = colnames(graph_coords),
    x = graph_coords[1, ],
    y = graph_coords[2, ],
    stringsAsFactors = FALSE
  )

  list(graph = graph_obj, graph_coords = graph_coords, edges = edge_df, nodes = node_df)
}

get_root_label_position <- function(cds, graph_bits, root_cells) {
  aux <- cds@principal_graph_aux[["UMAP"]]
  root_node <- NA_character_

  if (!is.null(aux$pr_graph_cell_proj_closest_vertex)) {
    closest_vertex <- flatten_closest_vertex(aux$pr_graph_cell_proj_closest_vertex)
    common_root <- intersect(root_cells, names(closest_vertex))
    if (length(common_root) > 0) {
      vals <- as.character(closest_vertex[common_root])
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) {
        root_node <- names(sort(table(vals), decreasing = TRUE))[1]
        root_node <- coerce_graph_vertex_name(root_node, graph_bits$graph)
      }
    }
  }

  if (!is.na(root_node) && root_node %in% graph_bits$nodes$node) {
    return(graph_bits$nodes %>% filter(node == root_node) %>% slice(1) %>% mutate(label = "ROOT"))
  }

  umap_mat <- reducedDims(cds)$UMAP
  common_root <- intersect(root_cells, rownames(umap_mat))
  if (length(common_root) > 0) {
    xy <- umap_mat[common_root, , drop = FALSE]
    return(data.frame(
      node = "root_centroid",
      x = median(xy[, 1], na.rm = TRUE),
      y = median(xy[, 2], na.rm = TRUE),
      label = "ROOT",
      stringsAsFactors = FALSE
    ))
  }

  NULL
}

project_point_to_segment <- function(px, py, ax, ay, bx, by) {
  abx <- bx - ax
  aby <- by - ay
  ab2 <- abx * abx + aby * aby

  if (!is.finite(ab2) || ab2 == 0) {
    return(list(x = ax, y = ay, t = 0, dist2 = (px - ax)^2 + (py - ay)^2))
  }

  t <- ((px - ax) * abx + (py - ay) * aby) / ab2
  t <- max(0, min(1, t))
  proj_x <- ax + t * abx
  proj_y <- ay + t * aby
  list(x = proj_x, y = proj_y, t = t, dist2 = (px - proj_x)^2 + (py - proj_y)^2)
}

project_cells_to_graph <- function(cds, edges_df, sample_map) {
  umap_mat <- reducedDims(cds)$UMAP
  cell_names <- rownames(umap_mat)
  pt_vec <- pseudotime(cds)
  sample_levels <- ordered_samples(sample_map[cell_names])

  bind_rows(lapply(seq_along(cell_names), function(i) {
    px <- umap_mat[i, 1]
    py <- umap_mat[i, 2]
    best_dist2 <- Inf
    best_p <- NULL
    best_e <- 1L

    for (e in seq_len(nrow(edges_df))) {
      pr <- project_point_to_segment(px, py, edges_df$x[e], edges_df$y[e], edges_df$xend[e], edges_df$yend[e])
      if (pr$dist2 < best_dist2) {
        best_dist2 <- pr$dist2
        best_p <- pr
        best_e <- e
      }
    }

    data.frame(
      cell = cell_names[i],
      umap_x = px,
      umap_y = py,
      graph_x = best_p$x,
      graph_y = best_p$y,
      projection_distance = sqrt(best_dist2),
      edge_index = best_e,
      sample = factor(sample_map[cell_names[i]], levels = sample_levels),
      pseudotime = as.numeric(pt_vec[cell_names[i]]),
      stringsAsFactors = FALSE
    )
  }))
}

compute_graph_weights <- function(graph_obj, graph_coords) {
  edge_df <- igraph::as_data_frame(graph_obj, what = "edges")
  if (nrow(edge_df) == 0) return(graph_obj)

  edge_weights <- purrr::map2_dbl(edge_df$from, edge_df$to, function(from_node, to_node) {
    from_xy <- graph_coords[, from_node]
    to_xy <- graph_coords[, to_node]
    sqrt(sum((from_xy - to_xy) ^ 2))
  })

  igraph::E(graph_obj)$weight <- edge_weights
  graph_obj
}

flatten_closest_vertex <- function(closest_vertex) {
  if (is.null(closest_vertex)) return(character())

  if (is.matrix(closest_vertex) || is.data.frame(closest_vertex)) {
    out <- as.character(closest_vertex[, 1, drop = TRUE])
    names(out) <- rownames(closest_vertex)
    return(out)
  }

  out <- as.character(closest_vertex)
  out
}

coerce_graph_vertex_name <- function(raw_vertex, graph_obj) {
  graph_nodes <- igraph::V(graph_obj)$name
  raw_vertex <- as.character(raw_vertex)[1]

  if (length(raw_vertex) == 0 || is.na(raw_vertex) || raw_vertex == "") {
    return(NA_character_)
  }

  if (raw_vertex %in% graph_nodes) {
    return(raw_vertex)
  }

  raw_num <- suppressWarnings(as.numeric(raw_vertex))
  if (!is.na(raw_num)) {
    raw_idx <- as.integer(raw_num)
    if (raw_idx >= 1 && raw_idx <= length(graph_nodes)) {
      return(graph_nodes[raw_idx])
    }
  }

  NA_character_
}

nearest_graph_vertex <- function(point_xy, graph_coords) {
  vertex_names <- colnames(graph_coords)
  dists <- apply(graph_coords, 2, function(node_xy) sqrt(sum((point_xy - node_xy) ^ 2)))
  vertex_names[which.min(dists)]
}

get_sample_medoid_cell <- function(umap_df) {
  if (nrow(umap_df) == 1) return(umap_df$cell[1])
  dist_mat <- as.matrix(dist(umap_df[, c("UMAP_1", "UMAP_2"), drop = FALSE]))
  umap_df$cell[which.min(rowSums(dist_mat))]
}

make_empty_matrix <- function(labels) {
  mat <- matrix(NA_real_, nrow = length(labels), ncol = length(labels))
  rownames(mat) <- labels
  colnames(mat) <- labels
  diag(mat) <- 0
  mat
}

matrix_to_long <- function(mat, method_name) {
  as.data.frame(as.table(mat), stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rename(sample_a = Var1, sample_b = Var2, distance = Freq) %>%
    mutate(method = method_name)
}
