####################
# parse_pipeline_logging.R
#
# Lightweight run-summary logging for Parse_Pipeline scripts.
#
# Inputs:
#   Optional parse_start_run() fields supplied by each workflow script.
#
# Outputs:
#   parse_outs/logs/run_summaries/<script>_<timestamp>.txt
#
# Methodology:
#   analysis/methodology/common/shared_configuration_and_logging_methodology.md
####################

if (!exists("parse_project_root")) {
  source("analysis/common/parse_pipeline_config.R")
}

parse_log_timestamp <- function(time = Sys.time()) {
  format(time, "%Y%m%d_%H%M%S")
}

parse_format_vector <- function(x) {
  if (is.null(x) || length(x) == 0) return("  - none")
  paste0("  - ", as.character(x), collapse = "\n")
}

parse_format_named_list <- function(x) {
  if (is.null(x) || length(x) == 0) return("  - none")
  paste0("  - ", names(x), ": ", vapply(x, paste, collapse = ", ", FUN.VALUE = character(1)), collapse = "\n")
}

parse_start_run <- function(script_name,
                            parameters = list(),
                            input_files = character(),
                            output_files = character(),
                            reused_cached = FALSE,
                            log_dir = file.path(parse_paths()$logs, "run_summaries")) {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  start_time <- Sys.time()
  run <- list(
    script_name = script_name,
    start_time = start_time,
    parameters = parameters,
    input_files = input_files,
    output_files = output_files,
    reused_cached = reused_cached,
    log_path = file.path(log_dir, paste0(gsub("[^A-Za-z0-9_]+", "_", script_name), "_", parse_log_timestamp(start_time), ".txt"))
  )
  message("Run summary will be written to: ", run$log_path)
  run
}

parse_finish_run <- function(run,
                             status = "success",
                             output_files = run$output_files,
                             reused_cached = run$reused_cached,
                             notes = character(),
                             include_session = TRUE) {
  end_time <- Sys.time()
  lines <- c(
    paste0("script: ", run$script_name),
    paste0("status: ", status),
    paste0("start_time: ", format(run$start_time, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("end_time: ", format(end_time, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("elapsed_minutes: ", round(as.numeric(difftime(end_time, run$start_time, units = "mins")), 3)),
    "",
    "input_files:",
    parse_format_vector(run$input_files),
    "",
    "output_files:",
    parse_format_vector(output_files),
    "",
    "parameters:",
    parse_format_named_list(run$parameters),
    "",
    paste0("reused_cached: ", paste(reused_cached, collapse = ", ")),
    "",
    "notes:",
    parse_format_vector(notes)
  )
  if (isTRUE(include_session)) {
    lines <- c(lines, "", "session_info:", paste0("  ", utils::capture.output(utils::sessionInfo())))
  }
  writeLines(lines, run$log_path)
  invisible(run$log_path)
}
