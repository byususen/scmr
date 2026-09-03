#!/usr/bin/env Rscript
# Usage: Rscript inst/scripts/02_run_article_simulation.R smoke /path/to/results
# Install this checkout first: R CMD INSTALL .
args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args)) args[1] else "smoke"
output <- if (length(args) > 1L) args[2] else "scmr-results"
if (!requireNamespace("scmr", quietly = TRUE) || utils::packageVersion("scmr") < "0.3.0") {
  stop("Install scmr >= 0.3.0 from this checkout before running the study.")
}
study_dir <- system.file("simulation", package = "scmr", mustWork = TRUE)
for (name in c("config", "io", "evaluate", "summarize", "run")) source(file.path(study_dir, paste0(name, ".R")))
run_article_simulation(article_config(profile, output_dir = output))
