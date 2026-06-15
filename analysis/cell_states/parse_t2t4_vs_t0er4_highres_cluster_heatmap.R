####################
# parse_t2t4_vs_t0er4_highres_cluster_heatmap.R
#
# Description:
#   Visualise high-resolution NMF MP cluster score changes per cell state
#   comparing T2+T4 vs T0+eR4 in Parse samples.
#   Uses manually annotated functional clusters of high-resolution metaprograms
#   derived from the PDO dataset (nMP156).
#
# Inputs:
#   parse_outs/by_samples/<sample>/Auto_<sample>_final.rds (T0, T2, T4, eR4)
#   parse_outs/Auto_parse_all_meta.rds
#   /rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_pdo_flot_highres_metaprogram_trends/Auto_pdo_flot_highres_geneNMF_metaprograms_nMP156.rds
#
# Outputs:
#   parse_outs/cell_states/t2t4_vs_t0er4_highres_clusters/...
#
# Cache / replot:
#   Saves calculated UCell scores to avoid re-running if intermediate file exists.
#
# Methodology:
#   analysis/methodology/cell_states/t2t4_vs_t0er4_highres_cluster_heatmap_methodology.md
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(stringr)
  library(UCell)
})

setwd("/rds/general/project/spatialtranscriptomics/ephemeral/Parse_Pipeline/parse_outs")
out_dir <- "cell_states/t2t4_vs_t0er4_highres_clusters"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "intermediate"), showWarnings = FALSE)

state_levels <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation"
)

state_labels_split <- c(
  "Classic\nProliferative",
  "Basal to\nIntest. Meta",
  "SMG-like\nMetaplasia",
  "Stress-\nadaptive",
  "3CA EMT &\nProt. mat."
)
names(state_labels_split) <- state_levels

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Stress-adaptive"       = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8"
)

increase_clusters <- list(
  "Replication stress &\ngenome maintenance" = c("MP31", "MP48", "MP52", "MP64", "MP39", "MP155"),
  "Chemotherapy-induced stress &\ninflammatory injury" = c("MP145", "MP56", "MP26", "MP47", "MP84"),
  "Wound-response &\nEMT-like plasticity" = c("MP128", "MP113"),
  "Mitotic /\nproliferative recovery" = c("MP33"),
  "Inflammatory-metabolic\nepithelial reprogramming" = c("MP49")
)

decrease_clusters <- list(
  "Differentiated epithelial /\nBarrett's lineage" = c("MP28", "MP32", "MP38", "MP46", "MP61"),
  "Lipid, xenobiotic &\ndetox metabolism" = c("MP55", "MP17", "MP75", "MP85"),
  "Immune modulation /\ninterferon response" = c("MP78", "MP62", "MP69", "MP70", "MP65"),
  "Stem / progenitor\nidentity & quiescence" = c("MP37", "MP21", "MP28"), 
  "ECM, adhesion &\nstromal interaction" = c("MP76", "MP44", "MP73")
)

all_clusters <- c(increase_clusters, decrease_clusters)
cluster_direction <- c(
  rep("T2T4 > T0eR4", length(increase_clusters)),
  rep("T2T4 < T0eR4", length(decrease_clusters))
)
names(cluster_direction) <- names(all_clusters)

ucell_cache <- file.path(out_dir, "intermediate", "parse_t2t4_vs_t0er4_ucell_scores.rds")

if (file.exists(ucell_cache)) {
  message("Loading cached UCell scores...")
  ucell <- readRDS(ucell_cache)
} else {
  message("Loading PDO nMP156 metaprograms...")
  pdo_nmp_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_pdo_flot_highres_metaprogram_trends/Auto_pdo_flot_highres_geneNMF_metaprograms_nMP156.rds"
  pdo_mps <- readRDS(pdo_nmp_path)
  mp.genes <- pdo_mps$metaprograms.genes
  names(mp.genes) <- paste0("MP", seq_along(mp.genes))
  
  # Only need MPs present in our clusters
  needed_mps <- unique(unlist(all_clusters))
  mp.genes <- mp.genes[names(mp.genes) %in% needed_mps]
  
  samples_to_process <- c("T0", "T2", "T4", "eR4")
  ucell_list <- list()
  
  for (smpl in samples_to_process) {
    message("Scoring sample: ", smpl)
    obj_path <- file.path("by_samples", smpl, paste0("Auto_", smpl, "_final.rds"))
    if (!file.exists(obj_path)) {
      stop("Sample file missing: ", obj_path)
    }
    obj <- readRDS(obj_path)
    # Ensure default assay is RNA
    DefaultAssay(obj) <- "RNA"
    
    # Calculate UCell
    obj <- AddModuleScore_UCell(obj, features = mp.genes, ncores = 4, name = "")
    
    # Extract matrix
    scores <- obj@meta.data[, names(mp.genes)]
    
    # Fix rownames to match global metadata (<sample>__<cell>)
    if (!any(grepl(paste0("^", smpl, "__"), rownames(scores)))) {
      rownames(scores) <- paste0(smpl, "__", rownames(scores))
    }
    
    ucell_list[[smpl]] <- scores
  }
  
  message("Combining UCell scores...")
  ucell <- do.call(rbind, unname(ucell_list))
  saveRDS(ucell, ucell_cache)
}

# Fix any prepended names from old cache or previous runs
rownames(ucell) <- sub("^[A-Za-z0-9]+\\.(?=[A-Za-z0-9]+__)", "", rownames(ucell), perl = TRUE)

message("Loading cell states...")
states <- readRDS("cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds")

states <- states %>%
  mutate(state = pdo_state)

emt_prot_states <- c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III")
states$state[states$state %in% emt_prot_states] <- "3CA_EMT_and_Protein_maturation"

common_cells <- intersect(rownames(ucell), states$cell)
cat("Common cells:", length(common_cells), "\n")

cell_df <- states %>%
  filter(cell %in% common_cells) %>%
  mutate(
    Treatment = ifelse(sample %in% c("T2", "T4"), "T2T4", 
                       ifelse(sample %in% c("T0", "eR4"), "T0eR4", NA))
  ) %>%
  filter(!is.na(Treatment)) %>%
  filter(state %in% state_levels)

cat("Cells with valid states and treatments:", nrow(cell_df), "\n")
print(table(cell_df$state, cell_df$Treatment))

# Filter state_levels to only those with valid cells to avoid downstream out-of-bounds errors
state_levels <- intersect(state_levels, unique(cell_df$state))

####################
# Compute cluster scores per cell
####################
message("Computing cluster scores per cell ...")

all_mps_needed <- unique(unlist(all_clusters))
available_mps <- intersect(all_mps_needed, colnames(ucell))
missing_mps <- setdiff(all_mps_needed, colnames(ucell))
if (length(missing_mps) > 0) {
  warning("MPs not found in UCell scores (will be skipped): ", paste(missing_mps, collapse = ", "))
}

cluster_scores <- sapply(names(all_clusters), function(cluster_name) {
  mps <- intersect(all_clusters[[cluster_name]], colnames(ucell))
  if (length(mps) == 0) return(rep(NA_real_, nrow(cell_df)))
  if (length(mps) == 1) {
    return(ucell[cell_df$cell, mps])
  }
  rowMeans(ucell[cell_df$cell, mps, drop = FALSE], na.rm = TRUE)
})
rownames(cluster_scores) <- cell_df$cell

####################
# Aggregate: mean cluster score per state × Treatment
####################
message("Aggregating per state × Treatment ...")

cluster_long <- as.data.frame(cluster_scores, check.names = FALSE) %>%
  mutate(cell = cell_df$cell, state = cell_df$state, Treatment = cell_df$Treatment) %>%
  pivot_longer(cols = all_of(names(all_clusters)), names_to = "cluster", values_to = "score")

# Note: Since Parse doesn't have multiple patients natively separated in the metadata in the same paired way as PDO,
# we aggregate directly by Treatment and compute the delta. 
# We don't do paired t-tests per patient here, we just use a two-sample t-test comparing all T2T4 cells vs T0eR4 cells for each state.

agg <- cluster_long %>%
  group_by(state, Treatment, cluster) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), n_cells = n(), .groups = "drop")

delta_df <- agg %>%
  pivot_wider(names_from = Treatment, values_from = c(mean_score, n_cells)) %>%
  filter(!is.na(mean_score_T0eR4), !is.na(mean_score_T2T4)) %>%
  mutate(mean_delta = mean_score_T2T4 - mean_score_T0eR4)

# Calculate p-values comparing cells directly
p_vals <- cluster_long %>%
  group_by(state, cluster) %>%
  summarise(
    p_value = tryCatch(
      t.test(score ~ Treatment)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

mean_delta <- delta_df %>%
  left_join(p_vals, by = c("state", "cluster")) %>%
  mutate(
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

mean_scores_abs <- agg

# Save CSVs
write.csv(mean_delta, file.path(out_dir, "tables", "parse_t2t4_vs_t0er4_cluster_delta_scores.csv"), row.names = FALSE)
write.csv(mean_scores_abs, file.path(out_dir, "tables", "parse_t2t4_vs_t0er4_cluster_absolute_scores.csv"), row.names = FALSE)

####################
# PLOT 1: Delta heatmap (T2T4 - T0eR4)
####################
message("Generating delta heatmap ...")

cluster_order <- names(all_clusters)
direction_labels <- cluster_direction

delta_mat <- mean_delta %>%
  select(cluster, state, mean_delta) %>%
  pivot_wider(names_from = state, values_from = mean_delta) %>%
  column_to_rownames("cluster") %>%
  as.matrix()
delta_mat <- delta_mat[cluster_order, state_levels, drop = FALSE]

sig_mat <- mean_delta %>%
  select(cluster, state, sig_label) %>%
  pivot_wider(names_from = state, values_from = sig_label) %>%
  column_to_rownames("cluster") %>%
  as.matrix()
sig_mat <- sig_mat[cluster_order, state_levels, drop = FALSE]
sig_mat[is.na(sig_mat)] <- ""

clip_val <- max(0.002, quantile(abs(delta_mat), 0.95, na.rm = TRUE))
col_fun <- colorRamp2(c(-clip_val, 0, clip_val), c("#245F7B", "white", "#B63E2F"))

row_direction <- factor(direction_labels[cluster_order],
                         levels = c("T2T4 > T0eR4", "T2T4 < T0eR4"))
ha_row <- rowAnnotation(
  Direction = row_direction,
  col = list(Direction = c("T2T4 > T0eR4" = "#B63E2F", "T2T4 < T0eR4" = "#245F7B")),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontface = "bold", fontsize = 9)
)

ha_top <- HeatmapAnnotation(
  State = state_levels,
  col = list(State = state_cols[state_levels]),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontface = "bold", fontsize = 9)
)

col_labels_delta <- state_labels_split[colnames(delta_mat)]

cell_fun_delta <- function(j, i, x, y, w, h, fill) {
  val <- delta_mat[i, j]
  lbl <- sig_mat[i, j]
  grid.text(sprintf("%.4f", val), x, y - unit(1, "mm"), gp = gpar(fontsize = 7))
  if (!is.na(lbl) && lbl != "") {
    grid.text(lbl, x, y + unit(2.5, "mm"), gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
  }
}

ht_delta <- Heatmap(
  delta_mat,
  name = "Mean \u0394 score\n(T2T4 \u2212 T0eR4)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_labels = col_labels_delta,
  column_names_rot = 0,
  column_names_centered = TRUE,
  show_column_names = TRUE,
  row_split = row_direction,
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_gap = unit(5, "mm"),
  left_annotation = ha_row,
  top_annotation = ha_top,
  cell_fun = cell_fun_delta,
  width = unit(ncol(delta_mat) * 2.5, "cm"),
  height = unit(nrow(delta_mat) * 1.2, "cm")
)

delta_pdf <- file.path(out_dir, "figures", "parse_t2t4_vs_t0er4_cluster_delta_heatmap.pdf")
message("Writing: ", delta_pdf)
pdf(delta_pdf, width = 10, height = 8)
draw(ht_delta,
     column_title = "High-Res MP Cluster Score Change (T2T4 vs T0eR4 per state)",
     column_title_gp = gpar(fontface = "bold", fontsize = 14),
     merge_legend = TRUE)
dev.off()

png(sub(".pdf$", ".png", delta_pdf), width = 10, height = 8, units = "in", res = 300)
draw(ht_delta,
     column_title = "High-Res MP Cluster Score Change (T2T4 vs T0eR4 per state)",
     column_title_gp = gpar(fontface = "bold", fontsize = 14),
     merge_legend = TRUE)
dev.off()

####################
# PLOT 2: Absolute scores heatmap (T0eR4 | T2T4 per state)
####################
message("Generating absolute score heatmap ...")

abs_wide <- mean_scores_abs %>%
  mutate(col_label = paste0(state, "\n", Treatment)) %>%
  select(cluster, col_label, mean_score) %>%
  pivot_wider(names_from = col_label, values_from = mean_score) %>%
  column_to_rownames("cluster")

col_order <- unlist(lapply(state_levels, function(st) {
  c(paste0(st, "\nT0eR4"), paste0(st, "\nT2T4"))
}))
col_order <- col_order[col_order %in% colnames(abs_wide)]
abs_mat <- as.matrix(abs_wide[cluster_order, col_order, drop = FALSE])

norm_mat <- t(apply(abs_mat, 1, scale))
rownames(norm_mat) <- rownames(abs_mat)
colnames(norm_mat) <- colnames(abs_mat)

col_state_split <- factor(
  rep(state_labels_split[state_levels], each = 2)[seq_along(col_order)],
  levels = state_labels_split[state_levels]
)

col_short <- gsub(".*\n", "", col_order)

norm_clip <- max(1.5, quantile(abs(norm_mat), 0.98, na.rm = TRUE))
col_fun_norm <- colorRamp2(c(-norm_clip, 0, norm_clip), c("#245F7B", "white", "#B63E2F"))

ha_row2 <- rowAnnotation(
  Direction = row_direction,
  col = list(Direction = c("T2T4 > T0eR4" = "#B63E2F", "T2T4 < T0eR4" = "#245F7B")),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontface = "bold", fontsize = 9)
)

ha_top2 <- HeatmapAnnotation(
  State = rep(state_levels, each = 2)[seq_along(col_order)],
  col = list(State = state_cols),
  show_annotation_name = FALSE
)

cell_fun_abs <- function(j, i, x, y, w, h, fill) {
  grid.text(sprintf("%.4f", abs_mat[i, j]), x, y, gp = gpar(fontsize = 6.5))
}

ht_abs <- Heatmap(
  norm_mat,
  name = "Row Z-score\n(UCell score)",
  col = col_fun_norm,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 30,
  show_column_names = TRUE,
  column_labels = col_short,
  column_split = col_state_split,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
  column_gap = unit(4, "mm"),
  row_split = row_direction,
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_gap = unit(5, "mm"),
  left_annotation = ha_row2,
  top_annotation = ha_top2,
  cell_fun = cell_fun_abs,
  width = unit(length(col_order) * 1.3, "cm"),
  height = unit(nrow(norm_mat) * 1.2, "cm")
)

abs_pdf <- file.path(out_dir, "figures", "parse_t2t4_vs_t0er4_cluster_absolute_heatmap.pdf")
message("Writing: ", abs_pdf)
pdf(abs_pdf, width = 12, height = 8)
draw(ht_abs,
     column_title = "High-Res MP Cluster Scores (T0eR4 & T2T4 per state)",
     column_title_gp = gpar(fontface = "bold", fontsize = 14),
     merge_legend = TRUE)
dev.off()

png(sub(".pdf$", ".png", abs_pdf), width = 12, height = 8, units = "in", res = 300)
draw(ht_abs,
     column_title = "High-Res MP Cluster Scores (T0eR4 & T2T4 per state)",
     column_title_gp = gpar(fontface = "bold", fontsize = 14),
     merge_legend = TRUE)
dev.off()

####################
# PLOT 3: MP-level Delta Heatmap (Side by Side)
####################
message("Generating MP-level delta heatmap ...")

mp_to_cluster <- data.frame(
  MP = unname(unlist(all_clusters)),
  Cluster = rep(names(all_clusters), lengths(all_clusters)),
  stringsAsFactors = FALSE
)
mp_to_cluster$Direction <- cluster_direction[mp_to_cluster$Cluster]

available_mps_plot3 <- intersect(mp_to_cluster$MP, colnames(ucell))
mp_to_cluster <- mp_to_cluster[mp_to_cluster$MP %in% available_mps_plot3, ]

unique_mps <- unique(mp_to_cluster$MP)

mp_long <- as.data.frame(ucell[cell_df$cell, unique_mps, drop=FALSE]) %>%
  mutate(cell = cell_df$cell, state = cell_df$state, Treatment = cell_df$Treatment) %>%
  pivot_longer(cols = all_of(unique_mps), names_to = "MP", values_to = "score")

mp_agg <- mp_long %>%
  group_by(state, Treatment, MP) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

mp_mean_delta <- mp_agg %>%
  pivot_wider(names_from = Treatment, values_from = mean_score) %>%
  filter(!is.na(T0eR4), !is.na(T2T4)) %>%
  mutate(mean_delta = T2T4 - T0eR4) %>%
  left_join(
    mp_long %>%
      group_by(state, MP) %>%
      summarise(
        p_value = tryCatch(t.test(score ~ Treatment)$p.value, error = function(e) NA_real_),
        .groups = "drop"
      ),
    by = c("state", "MP")
  ) %>%
  mutate(
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  left_join(mp_to_cluster, by = "MP", relationship = "many-to-many")

mp_inc <- mp_mean_delta %>% filter(Direction == "T2T4 > T0eR4")
mp_dec <- mp_mean_delta %>% filter(Direction == "T2T4 < T0eR4")

create_mp_ht <- function(df, title) {
  if (nrow(df) == 0) return(NULL)
  
  mat <- df %>% select(MP, state, mean_delta) %>% distinct() %>% pivot_wider(names_from = state, values_from = mean_delta) %>% column_to_rownames("MP") %>% as.matrix()
  mat <- mat[, state_levels[state_levels %in% colnames(mat)], drop = FALSE]
  
  sig <- df %>% select(MP, state, sig_label) %>% distinct() %>% pivot_wider(names_from = state, values_from = sig_label) %>% column_to_rownames("MP") %>% as.matrix()
  sig <- sig[rownames(mat), colnames(mat), drop = FALSE]
  sig[is.na(sig)] <- ""
  
  df_unique <- df %>% select(MP, Cluster) %>% distinct()
  cluster_lvls <- names(all_clusters)
  df_unique <- df_unique %>% mutate(Cluster = factor(Cluster, levels = cluster_lvls)) %>% arrange(Cluster, MP)
  
  mat <- mat[df_unique$MP, , drop = FALSE]
  sig <- sig[df_unique$MP, , drop = FALSE]
  rownames(mat) <- df_unique$MP
  rownames(sig) <- df_unique$MP
  
  row_split_fac <- factor(df_unique$Cluster, levels = cluster_lvls[cluster_lvls %in% df_unique$Cluster])
  
  clip <- max(0.002, quantile(abs(mat), 0.95, na.rm = TRUE))
  col_f <- colorRamp2(c(-clip, 0, clip), c("#245F7B", "white", "#B63E2F"))
  
  cluster_colors <- rainbow(length(cluster_lvls))
  names(cluster_colors) <- cluster_lvls
  
  ha_row <- rowAnnotation(
    Cluster = df_unique$Cluster,
    col = list(Cluster = cluster_colors),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  
  cell_f <- function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.4f", mat[i, j]), x, y, gp = gpar(fontsize = 7))
  }
  
  Heatmap(
    mat,
    name = "\u0394 Score",
    col = col_f,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 8, fontface = "bold"),
    column_names_gp = gpar(fontsize = 8, fontface = "bold"),
    column_labels = paste0("\n", state_labels_split[colnames(mat)]),
    column_names_rot = 0,
    column_names_centered = TRUE,
    row_split = row_split_fac,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 8, fontface = "bold"),
    row_gap = unit(3, "mm"),
    cell_fun = cell_f,
    left_annotation = ha_row,
    top_annotation = HeatmapAnnotation(
      State = colnames(mat),
      col = list(State = state_cols[colnames(mat)]),
      show_annotation_name = FALSE,
      show_legend = FALSE
    ),
    show_heatmap_legend = (title == "T2T4 > T0eR4"),
    heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"),
    width = unit(ncol(mat) * 2.5, "cm"),
    height = unit(nrow(mat) * 0.8, "cm"),
    column_title = paste(title, "MPs")
  )
}

ht_inc <- create_mp_ht(mp_inc, "T2T4 > T0eR4")
ht_dec <- create_mp_ht(mp_dec, "T2T4 < T0eR4")

mp_pdf <- file.path(out_dir, "figures", "parse_t2t4_vs_t0er4_MP_delta_heatmap.pdf")
message("Writing: ", mp_pdf)
pdf(mp_pdf, width = 20, height = 12)
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
if (!is.null(ht_inc)) draw(ht_inc, newpage = FALSE, heatmap_legend_side = "bottom")
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
if (!is.null(ht_dec)) draw(ht_dec, newpage = FALSE, heatmap_legend_side = "bottom")
popViewport()
dev.off()

png(sub(".pdf$", ".png", mp_pdf), width = 20, height = 12, units = "in", res = 300)
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
if (!is.null(ht_inc)) draw(ht_inc, newpage = FALSE, heatmap_legend_side = "bottom")
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
if (!is.null(ht_dec)) draw(ht_dec, newpage = FALSE, heatmap_legend_side = "bottom")
popViewport()
dev.off()

message("=== parse_t2t4_vs_t0er4_highres_cluster_heatmap.R completed successfully ===")
