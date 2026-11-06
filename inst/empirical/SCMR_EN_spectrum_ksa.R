# =============================================================================
# SCORE-SPECTRUM (PHASE-TRANSITION) CHECK ON THE FULL-DATA KSA FITS
#
# For every regime of the saved full-data fits (interpret/fits), computes the
# eigenvalues of the whitened segment score vectors and counts the spikes above
# the noise bulk, calibrated by a parametric bootstrap from the fitted
# probabilities (statistic R_k = lambda_k / median(lambda), invariant to a common
# inflation of the score variance). A regime with no spike is "stable": it shows
# no sign of a hidden sub-regime. No refitting; only aggregates are written.
#
# Run in RStudio (about 10-20 minutes):
#   source("C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax/SCMR_EN_spectrum_ksa.R")
# Output: interpret/spectrum_regimes.csv and interpret/spectrum_eigen.csv
# =============================================================================
syntax_dir0 <- Sys.getenv("SCMR_SYNTAX_DIR", "C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax")
v4_script <- Sys.getenv("SCMR_V4_SCRIPT", file.path(syntax_dir0, "SCMR_EN_empirical_article_v4.R"))
spec_funs <- Sys.getenv("SCMR_SPECTRUM_FUNS", file.path(syntax_dir0, "spectrum_funs.R"))
G_list <- as.integer(strsplit(Sys.getenv("SCMR_SPECTRUM_G", "1,3,6,8,10,12,15"), ",")[[1]])
n_boot <- as.integer(Sys.getenv("SCMR_SPECTRUM_BOOT", "200"))
designs_env <- Sys.getenv("SCMR_DESIGNS")
Sys.setenv(SCMR_DESIGNS = "none")
source(v4_script)                       # data and design matrices only
Sys.setenv(SCMR_DESIGNS = designs_env)
source(spec_funs)
suppressWarnings(loadNamespace("glmnet"))

int_dir <- file.path(output_dir, "interpret")
banner("SCORE SPECTRUM OF THE FULL-DATA FITS")
all_rows <- seq_len(nrow(dat))
global <- global_for("FULL", 0L, "S", all_rows, outer_seed + 999L)
x <- X_base[, global$x_colnames, drop = FALSE]
y <- dat$phase
units <- sort(unique(seg_id))
unit_num <- match(seg_id, units)

reg_rows <- list(); eig_rows <- list()
for (G in G_list) {
  path <- file.path(int_dir, "fits", sprintf("SCMR_EN_S_full_G%02d.rds", G))
  if (!file.exists(path)) { cat("  missing", basename(path), "\n"); next }
  f <- readRDS(path)
  ug <- as.integer(f$group_unit[as.character(units)])
  for (g in sort(unique(ug))) {
    rows_g <- which(ug[unit_num] == g)
    m_g <- length(unique(unit_num[rows_g]))
    if (m_g < 5L) next
    t0 <- proc.time()[["elapsed"]]
    sb <- spectrum_boot(x, y, f$prob_train, unit_num, rows_g, B = n_boot, seed = 100L * G + g)
    reg_rows[[length(reg_rows) + 1L]] <- data.frame(G = G, Regime = g, Segments = m_g, Rows = length(rows_g),
      d = sb$d, R1 = sb$R[1], R2 = sb$R[2], R3 = sb$R[3], Rcrit95 = sb$crit,
      Rcrit99 = stats::quantile(sb$boot, 0.99, names = FALSE), Spikes = sb$spikes,
      Seconds = round(proc.time()[["elapsed"]] - t0, 1))
    eig_rows[[length(eig_rows) + 1L]] <- data.frame(G = G, Regime = g, k = seq_len(min(20L, length(sb$R))),
                                                    R = sb$R[seq_len(min(20L, length(sb$R)))])
    cat(sprintf("  G = %2d regime %2d: %3d segments, R1 = %.2f (crit %.2f), spikes %d\n",
                G, g, m_g, sb$R[1], sb$crit, sb$spikes))
  }
  rm(f); invisible(gc())
}
reg <- do.call(rbind, reg_rows)
utils::write.csv(reg, file.path(int_dir, "spectrum_regimes.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, eig_rows), file.path(int_dir, "spectrum_eigen.csv"), row.names = FALSE)
cat("\nRegimes with at least one spike, by G:\n")
print(aggregate(cbind(Unstable = Spikes > 0, Regimes = 1) ~ G, reg, sum))
cat("SPECTRUM DONE\n")
