####################
# parse_pipeline_config.R
#
# Shared Parse_Pipeline configuration.
#
# Inputs:
#   Environment variables only:
#     AUTO_PARSE_ROOT_DIR: optional project-root override.
#
# Outputs:
#   No files are written unless parse_output_tiers(..., create = TRUE) is used.
#
# Methodology:
#   analysis/methodology/common/shared_configuration_and_logging_methodology.md
#
# Notes:
#   This file centralizes constants that were previously repeated across
#   downstream scripts. Workflow-specific scripts should source this file before
#   defining local paths, sample orders, plotting dimensions, or state thresholds.
####################

parse_project_root <- function() {
  env_root <- Sys.getenv("AUTO_PARSE_ROOT_DIR", unset = "")
  if (nzchar(env_root)) return(normalizePath(env_root, mustWork = FALSE))

  candidates <- c(
    "/rds/general/project/spatialtranscriptomics/ephemeral/Parse_Pipeline",
    "/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline"
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("Could not locate Parse_Pipeline root. Set AUTO_PARSE_ROOT_DIR.")
  }
  normalizePath(existing[1], mustWork = FALSE)
}

parse_paths <- function(root_dir = parse_project_root()) {
  parse_outs <- file.path(root_dir, "parse_outs")
  list(
    root_dir = root_dir,
    analysis_dir = file.path(root_dir, "analysis"),
    parse_outs = parse_outs,
    intermediate = file.path(parse_outs, "intermediate"),
    tables = file.path(parse_outs, "tables"),
    figures = file.path(parse_outs, "figures"),
    logs = file.path(parse_outs, "logs"),
    reports = file.path(parse_outs, "reports"),
    methodology = file.path(root_dir, "analysis", "methodology")
  )
}

parse_output_tiers <- function(base_dir, create = TRUE) {
  tiers <- list(
    intermediate = file.path(base_dir, "intermediate"),
    tables = file.path(base_dir, "tables"),
    figures = file.path(base_dir, "figures"),
    logs = file.path(base_dir, "logs"),
    reports = file.path(base_dir, "reports")
  )
  if (isTRUE(create)) {
    invisible(lapply(tiers, dir.create, recursive = TRUE, showWarnings = FALSE))
  }
  tiers
}

parse_samples <- c("T0", "T1", "T2", "T4", "R4", "eR4")
parse_all_samples <- c(
  parse_samples,
  "PDO",
  "SUR1090_Untreated",
  "SUR1090_Treated"
)
parse_external_samples <- c("PDO", "SUR1090_Untreated", "SUR1090_Treated")
parse_root_sample <- "T0"

parse_metadata_columns <- list(
  sample = "orig.ident",
  sample_fallbacks = c("sample", "orig.ident", "Sample", "sample_id"),
  cell = "cell",
  study = "study"
)

parse_state_definition <- list(
  preferred_method = "Approach B",
  normalization = "noreg",
  unresolved_label = "Unresolved",
  hybrid_label = "Hybrid",
  min_group_score = 0.5,
  hybrid_gap = 0.3,
  primary_assignment = "parse_outs/cell_states/Auto_parse_PDOpipeline_topmp_assignments.rds",
  primary_adjusted_scores = "parse_outs/cell_states/Auto_parse_PDOpipeline_mp_adj_noreg.rds"
)

parse_sample_colours <- c(
  "T0" = "#0072B2",
  "T1" = "#E69F00",
  "T2" = "#009E73",
  "T4" = "#D55E00",
  "R4" = "#CC79A7",
  "eR4" = "#56B4E9",
  "PDO" = "#4D4D4D",
  "SUR1090_Untreated" = "#0099A8",
  "SUR1090_Treated" = "#C44E52"
)

parse_plot_defaults <- list(
  dpi = 300,
  slide_base_size = 16,
  slide_axis_text_size = 13,
  slide_legend_text_size = 13,
  slide_legend_title_size = 14,
  slide_point_size = 2.7,
  pdf_width_wide = 16,
  pdf_height_wide = 9,
  pdf_width_tall = 12,
  pdf_height_tall = 10
)

parse_thresholds <- list(
  min_metaprogram_sample_coverage = 0.25,
  ucell_max_rank = 1500,
  highres_score_ncores = 2,
  pseudotime_min_cells_per_sample = 30
)

parse_reference_paths <- list(
  scatlas_epithelial = "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/EAC_Ref_epi.rds",
  scatlas_metaprograms = "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds",
  scatlas_ucell_scores = "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/UCell_nMP19_filtered.rds",
  pdo_metaprograms = "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds",
  three_ca_metaprograms = "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv",
  cell_cycle_genes = "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  developmental_enrichment_dir = "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage",
  carroll_reference = "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Carroll_2023_reference.rds",
  gene_order = "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt",
  velocity_gtf = "/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
)
