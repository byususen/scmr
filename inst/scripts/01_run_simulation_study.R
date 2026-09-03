#!/usr/bin/env Rscript
# Compatibility entry point: uses the maintained article runner (smoke by default).
# Usage: Rscript inst/scripts/01_run_simulation_study.R smoke /path/to/results
if (!requireNamespace("scmr", quietly = TRUE) || utils::packageVersion("scmr") < "0.3.0") {
  stop("Install scmr >= 0.3.0 from this checkout first.")
}
source(system.file("scripts", "02_run_article_simulation.R", package = "scmr", mustWork = TRUE))
