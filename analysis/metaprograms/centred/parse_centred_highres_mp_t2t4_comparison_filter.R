####################
# parse_centred_highres_mp_t2t4_comparison_filter.R
#
# Description:
#   Runs the existing T2/T4-high high-resolution MP comparison filter against
#   centred GeneNMF programmes, with all outputs isolated under parse_outs/centred/.
#   The only method change relative to the source workflow is that the upstream
#   NMF programmes were generated with multiNMF(center = TRUE).
#
# Inputs:
#   parse_outs/centred/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_geneNMF_metaprograms_nMP<half_total>.rds
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_UCell_scores_nMP<half_total>.rds
#   parse_outs/by_samples/<sample>/Auto_<sample>_final.rds
#   3CA, cell-cycle, and developmental enrichment reference files listed in AGENTS.md
#
# Outputs:
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/*
#   parse_outs/logs/run_summaries/parse_centred_highres_mp_t2t4_comparison_filter_*.txt
#
# Cache / replot:
#   Uses the same --mode argument as the uncentred workflow: score, enrich,
#   excel, or all. Cached centred nMP117 metaprograms, UCell scores, trend
#   tables, selected genes, and annotation outputs are reused when present.
#
# Methodology:
#   analysis/methodology/metaprograms/legacy_highres_mp_t2t4_comparison_filter_methodology.md
#
# Downstream status:
#   Feeds centred publication-style T2/T4-high MP heatmaps and centred-vs-
#   uncentred T2/T4 best-match comparison only.
####################

source("analysis/common/parse_pipeline_config.R")

source_path <- file.path(
  parse_project_root(),
  "analysis",
  "metaprograms",
  "parse_highres_mp_t2t4_comparison_filter.R"
)
if (!file.exists(source_path)) {
  stop("Missing source T2/T4 comparison workflow: ", source_path)
}

script_text <- readLines(source_path, warn = FALSE)

replace_fixed <- function(x, pattern, replacement) {
  gsub(pattern, replacement, x, fixed = TRUE)
}

script_text <- replace_fixed(
  script_text,
  "legacy_parse_highres_mp_t2t4_comparison_filter",
  "parse_centred_highres_mp_t2t4_comparison_filter"
)
script_text <- replace_fixed(
  script_text,
  "parse_highres_mp_t2t4_comparison_filter",
  "parse_centred_highres_mp_t2t4_comparison_filter"
)
script_text <- replace_fixed(
  script_text,
  "parse_outs/Auto_parse_metaprograms",
  "parse_outs/centred/Auto_parse_metaprograms"
)
script_text <- replace_fixed(
  script_text,
  "parse_outs/Auto_parse_highres_metaprogram_trends",
  "parse_outs/centred/Auto_parse_highres_metaprogram_trends"
)
script_text <- replace_fixed(
  script_text,
  'parse_mp_dir <- file.path(qc_dir, "Auto_parse_metaprograms")',
  'parse_mp_dir <- file.path(qc_dir, "centred", "Auto_parse_metaprograms")'
)
script_text <- replace_fixed(
  script_text,
  'base_highres_dir <- file.path(qc_dir, "Auto_parse_highres_metaprogram_trends")',
  'base_highres_dir <- file.path(qc_dir, "centred", "Auto_parse_highres_metaprogram_trends")'
)

eval(parse(text = script_text), envir = globalenv())
