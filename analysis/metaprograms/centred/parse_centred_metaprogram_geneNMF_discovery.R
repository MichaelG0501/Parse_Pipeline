####################
# parse_centred_metaprogram_geneNMF_discovery.R
#
# Description:
#   Runs the Parse six-sample GeneNMF programme sweep with multiNMF(center = TRUE).
#   This is an isolated comparison workflow for the centred-data NMF method.
#
# Inputs:
#   parse_outs/by_samples/<sample>/Auto_<sample>_final.rds
#   Optional cache: parse_outs/centred/Auto_parse_metaprograms/Auto_parse_list_geneNMF.rds
#   Optional cache: parse_outs/centred/Auto_parse_metaprograms/Auto_parse_merged_geneNMF.rds
#
# Outputs:
#   parse_outs/centred/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds
#   parse_outs/centred/Auto_parse_metaprograms/Auto_parse_geneNMF_program_config.csv
#   parse_outs/logs/run_summaries/parse_centred_metaprogram_geneNMF_discovery_*.txt
#
# Cache / replot:
#   Cached sample-list, merged object, and centred GeneNMF outputs are reused
#   when present. Delete or move centred cache files manually before rerunning
#   if a fresh centred NMF sweep is required.
#
# Methodology:
#   analysis/methodology/metaprograms/metaprogram_geneNMF_discovery_methodology.md
#
# Downstream status:
#   Feeds the centred high-resolution trend, annotation, heatmap, and comparison
#   workflows only. Does not replace the canonical uncentred GeneNMF outputs.
####################

source("analysis/common/parse_pipeline_config.R")
source("analysis/common/parse_pipeline_logging.R")

script_run <- parse_start_run(
  "parse_centred_metaprogram_geneNMF_discovery",
  parameters = list(
    samples = paste(parse_samples, collapse = ","),
    k_range = "4:9",
    center = TRUE,
    min.exp = 0.05
  ),
  input_files = "parse_outs/by_samples/<sample>/Auto_<sample>_final.rds",
  output_files = c(
    "parse_outs/centred/Auto_parse_metaprograms/Auto_parse_geneNMF_outs.rds",
    "parse_outs/centred/Auto_parse_metaprograms/Auto_parse_geneNMF_program_config.csv"
  )
)
script_run_status <- "failed"
on.exit(parse_finish_run(script_run, status = script_run_status), add = TRUE)

suppressPackageStartupMessages({
  library(GeneNMF)
  library(Seurat)
  library(SeuratObject)
})

project_dir <- parse_project_root()
paths <- parse_paths(project_dir)
centred_dir <- file.path(paths$parse_outs, "centred")
parse_mp_dir <- file.path(centred_dir, "Auto_parse_metaprograms")
dir.create(parse_mp_dir, recursive = TRUE, showWarnings = FALSE)

sample_files <- setNames(
  file.path(paths$parse_outs, "by_samples", parse_samples, paste0("Auto_", parse_samples, "_final.rds")),
  parse_samples
)
missing_files <- names(sample_files)[!file.exists(sample_files)]
if (length(missing_files) > 0) {
  stop("Missing final sample RDS files: ", paste(missing_files, collapse = ", "))
}

parse_list_path <- file.path(parse_mp_dir, "Auto_parse_list_geneNMF.rds")
parse_merged_path <- file.path(parse_mp_dir, "Auto_parse_merged_geneNMF.rds")
geneNMF_out_path <- file.path(parse_mp_dir, "Auto_parse_geneNMF_outs.rds")
config_path <- file.path(parse_mp_dir, "Auto_parse_geneNMF_program_config.csv")

if (file.exists(parse_list_path)) {
  message("Loading existing centred sample list: ", parse_list_path)
  parse.list <- readRDS(parse_list_path)
} else {
  message("Building centred GeneNMF sample list.")
  parse.list <- lapply(sample_files, readRDS)
  names(parse.list) <- parse_samples
  parse.list <- parse.list[parse_samples]
  saveRDS(parse.list, parse_list_path, compress = FALSE)
}

if (file.exists(parse_merged_path)) {
  message("Loading existing centred merged object: ", parse_merged_path)
  parse_obj <- readRDS(parse_merged_path)
} else {
  message("Building merged object cache for centred GeneNMF comparison.")
  parse_obj <- merge(parse.list[[1]], y = parse.list[-1], add.cell.ids = parse_samples, project = "parse_centred_geneNMF")
  parse_obj <- NormalizeData(parse_obj, verbose = FALSE)
  saveRDS(parse_obj, parse_merged_path, compress = FALSE)
}

if (file.exists(geneNMF_out_path)) {
  message("Loading existing centred GeneNMF programmes: ", geneNMF_out_path)
  geneNMF.programs <- readRDS(geneNMF_out_path)
} else {
  message("Running centred multiNMF with center = TRUE.")
  geneNMF.programs <- multiNMF(
    parse.list,
    assay = "RNA",
    k = 4:9,
    min.exp = 0.05,
    center = TRUE
  )
  saveRDS(geneNMF.programs, file = geneNMF_out_path, compress = FALSE)
}

program_counts <- vapply(geneNMF.programs, function(x) ncol(x$w), numeric(1))
config <- data.frame(
  sample = names(program_counts),
  nmf_programmes = as.integer(program_counts),
  center = TRUE,
  k_range = "4:9",
  min_exp = 0.05,
  stringsAsFactors = FALSE
)
write.csv(config, config_path, row.names = FALSE)

message("Total centred NMF programmes: ", sum(program_counts))
script_run_status <- "success"
message("parse_centred_metaprogram_geneNMF_discovery.R completed successfully.")
