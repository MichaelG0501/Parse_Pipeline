####################
# Compatibility entry point for the Parse sample-level distance matrix workflow.
# Samples are used in place of states.
####################

script_candidates <- c(
  "analysis/trajectory/Auto_parse_pseudotime_sample_distance_matrix.R",
  "analysis/cell_states/Auto_parse_pseudotime_sample_distance_matrix.R",
  "Auto_parse_pseudotime_sample_distance_matrix.R"
)
script_path <- script_candidates[file.exists(script_candidates)][1]
if (is.na(script_path)) stop("Could not find Auto_parse_pseudotime_sample_distance_matrix.R")
source(script_path)
