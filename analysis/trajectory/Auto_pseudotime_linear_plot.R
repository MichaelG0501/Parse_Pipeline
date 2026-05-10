####################
# Compatibility entry point for the Parse sample-level linear trajectory report.
####################

script_candidates <- c(
  "analysis/trajectory/Auto_parse_pseudotime_linear_plot.R",
  "analysis/cell_states/Auto_parse_pseudotime_linear_plot.R",
  "Auto_parse_pseudotime_linear_plot.R"
)
script_path <- script_candidates[file.exists(script_candidates)][1]
if (is.na(script_path)) stop("Could not find Auto_parse_pseudotime_linear_plot.R")
source(script_path)
