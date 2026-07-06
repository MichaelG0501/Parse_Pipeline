####################
# parse_centred_highres_mp_strict_mean_median_trend_filter.R
#
# Description:
#   Runs the existing high-resolution strict mean/median trend and annotation
#   workflow against centred GeneNMF programmes, with all outputs isolated under
#   parse_outs/centred/.
#
# Inputs:
#   parse_outs/centred/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds
#   parse_outs/by_samples/<sample>/Auto_<sample>_final.rds
#   3CA, cell-cycle, and developmental enrichment reference files listed in AGENTS.md
#
# Outputs:
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_selected_mp_genes_nMP<half_total>.{rds,csv}
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_trend_summary_nMP<half_total>.csv
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_top_3CA_noncellcycle_nMP<half_total>.csv
#   parse_outs/centred/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_enrichment_annotation_nMP<half_total>.pdf
#   parse_outs/logs/run_summaries/parse_centred_highres_mp_strict_mean_median_trend_filter_*.txt
#
# Cache / replot:
#   Uses the same --mode argument as the uncentred workflow: score, enrich,
#   excel, or all. Cached centred metaprograms, UCell scores, annotations, and
#   reports under parse_outs/centred/ are reused when present.
#
# Methodology:
#   analysis/methodology/metaprograms/highres_mp_strict_mean_median_trend_filter_methodology.md
#
# Downstream status:
#   Feeds centred MP heatmaps and centred-vs-uncentred comparison only.
####################

source("analysis/common/parse_pipeline_config.R")

source_path <- file.path(
  parse_project_root(),
  "analysis",
  "metaprograms",
  "parse_highres_mp_strict_mean_median_trend_filter.R"
)
if (!file.exists(source_path)) {
  stop("Missing source high-resolution workflow: ", source_path)
}

script_text <- readLines(source_path, warn = FALSE)

replace_fixed <- function(x, pattern, replacement) {
  gsub(pattern, replacement, x, fixed = TRUE)
}

script_text <- replace_fixed(
  script_text,
  "parse_highres_mp_strict_mean_median_trend_filter",
  "parse_centred_highres_mp_strict_mean_median_trend_filter"
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
  'out_dir <- file.path(qc_dir, "Auto_parse_highres_metaprogram_trends")',
  'out_dir <- file.path(qc_dir, "centred", "Auto_parse_highres_metaprogram_trends")'
)

eval(parse(text = script_text), envir = globalenv())
