####################
# DELETE CANDIDATE.
#
# Compatibility entry point for the old state-distance script name. The active
# sample-distance workflow is:
#   analysis/trajectory/parse_pseudotime_sample_distance_maps.R
#
# This wrapper is retained only so the user can delete it manually after
# confirming no old PBS wrapper still calls it.
####################

script_candidates <- c(
  "analysis/trajectory/parse_pseudotime_sample_distance_maps.R"
)
script_path <- script_candidates[file.exists(script_candidates)][1]
if (is.na(script_path)) stop("Could not find parse_pseudotime_sample_distance_maps.R")
source(script_path)
