####################
# DELETE CANDIDATE.
#
# Compatibility entry point for the old Parse sample-level linear trajectory
# report name. The active script is:
#   analysis/trajectory/parse_pseudotime_linear_sample_report.R
#
# This wrapper is retained only so the user can delete it manually after
# confirming no old PBS wrapper still calls it.
####################

script_candidates <- c(
  "analysis/trajectory/parse_pseudotime_linear_sample_report.R"
)
script_path <- script_candidates[file.exists(script_candidates)][1]
if (is.na(script_path)) stop("Could not find parse_pseudotime_linear_sample_report.R")
source(script_path)
