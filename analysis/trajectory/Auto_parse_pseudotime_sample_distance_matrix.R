####################
# Auto_parse_pseudotime_sample_distance_matrix.R
#
# Monocle3-based sample-to-sample distance comparison for Parse data.
# Samples replace states throughout the reference workflow. The trajectory is
# rooted on T0, and the primary output includes T0-to-sample distances such as
# T0 to T1 and T0 to R4.
####################

helper_candidates <- c(
  "analysis/trajectory/Auto_parse_pseudotime_helpers.R",
  "analysis/cell_states/Auto_parse_pseudotime_helpers.R",
  "Auto_parse_pseudotime_helpers.R"
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) stop("Could not find Auto_parse_pseudotime_helpers.R")
source(helper_path)

distance_dir <- file.path(pseudotime_dir, "sample_distance_pseudotime")
dir.create(distance_dir, recursive = TRUE, showWarnings = FALSE)

build_heatmap_pdf <- function(long_df, sample_levels, file_path) {
  if (nrow(long_df) == 0) return(invisible(NULL))

  plot_df <- long_df %>%
    group_by(method) %>%
    mutate(
      method_range = max(distance, na.rm = TRUE) - min(distance, na.rm = TRUE),
      distance_scaled = ifelse(
        is.finite(method_range) & method_range > 0,
        (distance - min(distance, na.rm = TRUE)) / method_range,
        0
      )
    ) %>%
    ungroup() %>%
    mutate(
      sample_a = factor(sample_a, levels = sample_levels),
      sample_b = factor(sample_b, levels = sample_levels)
    )

  p <- ggplot(plot_df, aes(x = sample_a, y = sample_b, fill = distance_scaled)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(
      aes(label = ifelse(is.na(distance), "NA", sprintf("%.2f", distance))),
      size = 4.1,
      fontface = "bold",
      colour = "black"
    ) +
    scale_fill_gradient(
      low = "#fff5f0",
      high = "#99000d",
      na.value = "grey90",
      name = "Within-method\nscaled distance"
    ) +
    facet_wrap(~method) +
    theme_minimal(base_size = 15) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 13, face = "bold", colour = "black"),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 13),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 12, colour = "grey30")
    ) +
    labs(
      title = "Parse Sample Distance Method Comparison",
      subtitle = "Tile colour is scaled within each method; labels show raw distances."
    )

  ggsave(file_path, p, width = 15, height = 10)
}

####################
# Interconnected Node Plot (Network style)
# Each node is a sample. Layout is determined by MDS (Multi-Dimensional Scaling)
# such that Euclidean distance between nodes reflects the pseudotime distance.
# All samples are interconnected with edges.
####################
build_interconnected_network_pdf <- function(matrix_list, sample_levels, sample_cols, file_path) {
  plot_list <- list()
  
  for (method_name in names(matrix_list)) {
    mat <- matrix_list[[method_name]]
    present_samples <- intersect(sample_levels, rownames(mat))
    mat_sub <- mat[present_samples, present_samples, drop = FALSE]
    
    # Remove rows/cols with all NA or zero (except for the root)
    keep <- apply(mat_sub, 1, function(x) any(is.finite(x) & x > 0)) | 
           apply(mat_sub, 2, function(x) any(is.finite(x) & x > 0)) |
           (rownames(mat_sub) == root_sample)
    
    mat_sub <- mat_sub[keep, keep, drop = FALSE]
    if (nrow(mat_sub) < 2) next
    
    # Symmetrize and ensure no NAs for MDS
    mat_sub[is.na(mat_sub)] <- 0
    dist_obj <- as.dist(mat_sub)
    
    # MDS layout
    set.seed(12345)
    mds_res <- tryCatch({
      cmdscale(dist_obj, k = 2, eig = TRUE)
    }, error = function(e) NULL)
    
    if (is.null(mds_res) || is.null(mds_res$points)) {
      message("MDS failed for method: ", method_name)
      next
    }
    
    layout_df <- as.data.frame(mds_res$points)
    colnames(layout_df) <- c("x", "y")
    layout_df$sample <- rownames(layout_df)
    layout_df$sample <- factor(layout_df$sample, levels = sample_levels)
    
    # Create edges for all pairs
    edge_df <- as.data.frame(as.table(mat_sub)) %>%
      rename(sample_a = Var1, sample_b = Var2, distance = Freq) %>%
      mutate(sample_a = as.character(sample_a), sample_b = as.character(sample_b)) %>%
      filter(sample_a < sample_b) %>%
      filter(distance > 0) %>%
      left_join(layout_df %>% select(sample, x, y), by = c("sample_a" = "sample")) %>%
      left_join(layout_df %>% select(sample, xend = x, yend = y), by = c("sample_b" = "sample"))
    
    # Add jitter to ensure valid plotting even for 1D MDS layouts
    set.seed(123)
    x_range <- diff(range(layout_df$x, na.rm=TRUE))
    if (x_range == 0) x_range <- 1
    layout_df$x <- layout_df$x + runif(nrow(layout_df), -1e-2 * x_range, 1e-2 * x_range)
    layout_df$y <- layout_df$y + runif(nrow(layout_df), -1e-2 * x_range, 1e-2 * x_range)

    # Re-project edges to jittered coords
    edge_df <- as.data.frame(as.table(mat_sub)) %>%
      rename(sample_a = Var1, sample_b = Var2, distance = Freq) %>%
      mutate(sample_a = as.character(sample_a), sample_b = as.character(sample_b)) %>%
      filter(sample_a < sample_b) %>%
      filter(distance > 0) %>%
      left_join(layout_df %>% select(sample, x, y), by = c("sample_a" = "sample")) %>%
      left_join(layout_df %>% select(sample, xend = x, yend = y), by = c("sample_b" = "sample"))

    # Normalize distance for aesthetics (edge strength)
    if (nrow(edge_df) > 0) {
      edge_df <- edge_df %>%
        mutate(
          inv_dist = 1 / (distance + 1e-6),
          rel_inv_dist = (inv_dist - min(inv_dist)) / (max(inv_dist) - min(inv_dist) + 1e-6)
        )
    }

    p <- ggplot() +
      geom_segment(
        data = edge_df,
        aes(x = x, y = y, xend = xend, yend = yend, alpha = rel_inv_dist, linewidth = rel_inv_dist),
        colour = "grey70",
        lineend = "round"
      ) +
      geom_point(
        data = layout_df,
        aes(x = x, y = y, fill = sample),
        shape = 21,
        size = 11,
        colour = "black",
        stroke = 1.2
      ) +
      ggrepel::geom_text_repel(
        data = layout_df,
        aes(x = x, y = y, label = sample),
        fontface = "bold",
        size = 8,
        box.padding = 0.6,
        point.padding = 0.4
      ) +
      scale_fill_manual(values = sample_cols) +
      scale_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
      scale_linewidth_continuous(range = c(1.5, 5.0), guide = "none") +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        plot.margin = margin(20, 20, 20, 20),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()
      ) +
      labs(title = gsub("_", " ", method_name))
    
    plot_list[[method_name]] <- p
  }
  
  if (length(plot_list) > 0) {
    combined <- wrap_plots(plot_list, ncol = 3) +
      plot_annotation(
        title = "Sample Distance Network",
        theme = theme(
          plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
          plot.subtitle = element_blank()
        )
      )
    
    ggsave(file_path, combined, width = 24, height = 14)
  }
}

message("=== Parse sample distance matrix ===")
result <- get_or_build_parse_trajectory(rebuild = FALSE)
cds <- result$cds
sample_levels <- result$sample_levels

sample_map <- as.character(colData(cds)$sample_label)
names(sample_map) <- colnames(cds)

pt_vec <- pseudotime(cds)
pt_vec[is.infinite(pt_vec)] <- NA_real_

umap_mat <- reducedDims(cds)$UMAP
umap_df <- data.frame(
  cell = rownames(umap_mat),
  UMAP_1 = umap_mat[, 1],
  UMAP_2 = umap_mat[, 2],
  sample = sample_map[rownames(umap_mat)],
  pseudotime = as.numeric(pt_vec[rownames(umap_mat)]),
  stringsAsFactors = FALSE
) %>%
  filter(sample %in% sample_levels)

graph_bits <- extract_graph_structure(cds)
graph_obj <- compute_graph_weights(graph_bits$graph, graph_bits$graph_coords)
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- flatten_closest_vertex(closest_vertex)

sample_summary <- umap_df %>%
  mutate(sample = factor(sample, levels = sample_levels)) %>%
  group_by(sample) %>%
  summarise(
    n_cells = n(),
    n_cells_with_pseudotime = sum(is.finite(pseudotime)),
    median_pseudotime = safe_median(pseudotime),
    mean_pseudotime = safe_mean(pseudotime),
    medoid_cell = get_sample_medoid_cell(cur_data()),
    centroid_x = mean(UMAP_1, na.rm = TRUE),
    centroid_y = mean(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    medoid_vertex = coerce_graph_vertex_name(closest_vertex[medoid_cell], graph_obj),
    centroid_vertex = nearest_graph_vertex(c(centroid_x, centroid_y), graph_bits$graph_coords)
  ) %>%
  ungroup()

present_samples <- as.character(sample_summary$sample)
directed_median_mat <- make_empty_matrix(sample_levels)
directed_mean_mat <- make_empty_matrix(sample_levels)
geodesic_medoid_mat <- make_empty_matrix(sample_levels)
geodesic_centroid_mat <- make_empty_matrix(sample_levels)
umap_centroid_mat <- make_empty_matrix(sample_levels)

for (sample_a in present_samples) {
  for (sample_b in present_samples) {
    row_a <- sample_summary %>% filter(sample == sample_a)
    row_b <- sample_summary %>% filter(sample == sample_b)
    if (nrow(row_a) == 0 || nrow(row_b) == 0) next

    directed_median_mat[sample_a, sample_b] <- abs(row_b$median_pseudotime[1] - row_a$median_pseudotime[1])
    directed_mean_mat[sample_a, sample_b] <- abs(row_b$mean_pseudotime[1] - row_a$mean_pseudotime[1])

    medoid_distance <- suppressWarnings(
      igraph::distances(
        graph_obj,
        v = row_a$medoid_vertex[1],
        to = row_b$medoid_vertex[1],
        weights = igraph::E(graph_obj)$weight
      )[1, 1]
    )
    geodesic_medoid_mat[sample_a, sample_b] <- if (is.finite(medoid_distance)) as.numeric(medoid_distance) else NA_real_

    centroid_distance <- suppressWarnings(
      igraph::distances(
        graph_obj,
        v = row_a$centroid_vertex[1],
        to = row_b$centroid_vertex[1],
        weights = igraph::E(graph_obj)$weight
      )[1, 1]
    )
    geodesic_centroid_mat[sample_a, sample_b] <- if (is.finite(centroid_distance)) as.numeric(centroid_distance) else NA_real_

    umap_centroid_mat[sample_a, sample_b] <- sqrt(
      (row_a$centroid_x[1] - row_b$centroid_x[1]) ^ 2 +
        (row_a$centroid_y[1] - row_b$centroid_y[1]) ^ 2
    )
  }
}

matrix_list <- list(
  directed_pseudotime_median = directed_median_mat,
  directed_pseudotime_mean = directed_mean_mat,
  principal_graph_geodesic_medoid = geodesic_medoid_mat,
  principal_graph_geodesic_centroid = geodesic_centroid_mat,
  umap_centroid_euclidean = umap_centroid_mat
)

saveRDS(matrix_list, file.path(distance_dir, "Auto_parse_sample_distance_matrices.rds"))
write_csv(sample_summary, file.path(distance_dir, "Auto_parse_sample_distance_sample_summary.csv"))

for (method_name in names(matrix_list)) {
  write.csv(
    matrix_list[[method_name]],
    file.path(distance_dir, paste0("Auto_parse_", method_name, "_sample_matrix.csv")),
    quote = FALSE
  )
}

long_df <- bind_rows(lapply(names(matrix_list), function(method_name) {
  matrix_to_long(matrix_list[[method_name]], method_name)
}))
write_csv(long_df, file.path(distance_dir, "Auto_parse_sample_distance_long.csv"))

root_distance_df <- bind_rows(lapply(names(matrix_list), function(method_name) {
  mat <- matrix_list[[method_name]]
  data.frame(
    root_sample = root_sample,
    target_sample = colnames(mat),
    method = method_name,
    distance = as.numeric(mat[root_sample, ]),
    stringsAsFactors = FALSE
  )
})) %>%
  filter(target_sample != root_sample) %>%
  arrange(method, target_sample)

write_csv(root_distance_df, file.path(distance_dir, "Auto_parse_sample_distance_from_T0.csv"))
write_csv(root_distance_df, file.path(summary_dir, "Auto_parse_sample_distance_from_T0.csv"))

build_heatmap_pdf(
  long_df,
  sample_levels,
  file.path(distance_dir, "Auto_parse_sample_distance_method_comparison_heatmap.pdf")
)

build_interconnected_network_pdf(
  matrix_list,
  sample_levels,
  get_sample_colours(sample_levels),
  file.path(distance_dir, "Auto_parse_sample_distance_from_T0_nodeplot.pdf")
)

message("Saved Parse sample distance outputs to: ", distance_dir)
