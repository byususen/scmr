# =============================================================================
# PHASE-COEXISTENCE CHECK ON THE FULL-DATA KSA FITS (after SCMR_EN_spectrum_ksa.R)
#
# Uses the spike counts already in interpret/spectrum_regimes.csv (no bootstrap
# is repeated). For each spike eigenvector of each regime:
#   (a) spatial order   : Moran's I of the segment loadings, 9,999 permutations;
#   (b) two phases      : bimodality of the loadings against continuous
#                         spatially autocorrelated fields with the same Moran's I.
# A regime holds a discrete sub-regime only if (a) and (b) both reject; p-values
# are Bonferroni-adjusted over the spikes of a regime and Holm-adjusted over the
# regimes of each G. Aggregates only.
#
# Run in RStudio (a few minutes):
#   source("C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax/SCMR_EN_phase_ksa.R")
# Output: interpret/phase_regimes.csv, interpret/phase_spikes.csv
# =============================================================================
syntax_dir0 <- Sys.getenv("SCMR_SYNTAX_DIR", "C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax")
v4_script <- Sys.getenv("SCMR_V4_SCRIPT", file.path(syntax_dir0, "SCMR_EN_empirical_article_v4.R"))
spec_funs <- Sys.getenv("SCMR_SPECTRUM_FUNS", file.path(syntax_dir0, "spectrum_funs.R"))
n_perm <- as.integer(Sys.getenv("SCMR_PHASE_PERM", "9999"))
n_null <- as.integer(Sys.getenv("SCMR_PHASE_NULL", "999"))
designs_env <- Sys.getenv("SCMR_DESIGNS")
Sys.setenv(SCMR_DESIGNS = "none")
source(v4_script)
Sys.setenv(SCMR_DESIGNS = designs_env)
source(spec_funs)
suppressWarnings(loadNamespace("glmnet"))

int_dir <- file.path(output_dir, "interpret")
spec <- utils::read.csv(file.path(int_dir, "spectrum_regimes.csv"))
banner("PHASE-COEXISTENCE CHECK")
all_rows <- seq_len(nrow(dat))
global <- global_for("FULL", 0L, "S", all_rows, outer_seed + 999L)
x <- X_base[, global$x_colnames, drop = FALSE]
y <- dat$phase
units <- sort(unique(seg_id))
unit_num <- match(seg_id, units)
uxy <- as.matrix(seg_xy)[match(units, seg_id), , drop = FALSE]
knn_adj <- function(xy, k = 8L) {
  D <- as.matrix(stats::dist(xy)); diag(D) <- Inf
  A <- matrix(0, nrow(xy), nrow(xy))
  for (i in seq_len(nrow(xy))) A[i, order(D[i, ])[seq_len(min(k, nrow(xy) - 1L))]] <- 1
  pmax(A, t(A))
}
reg_rows <- list(); spk_rows <- list()
for (G in sort(unique(spec$G))) {
  f <- readRDS(file.path(int_dir, "fits", sprintf("SCMR_EN_S_full_G%02d.rds", G)))
  ug <- as.integer(f$group_unit[as.character(units)])
  sg <- spec[spec$G == G, ]
  for (r in seq_len(nrow(sg))) {
    g <- sg$Regime[r]
    rows_g <- which(ug[unit_num] == g)
    g_units <- sort(unique(unit_num[rows_g]))
    si <- unit_scores_info(x, y, f$prob_train, unit_num, rows_g)
    ph <- phase_order_test(si$U, si$n_b, si$J, knn_adj(uxy[g_units, , drop = FALSE]), sg$Spikes[r],
                           k = 5L, n_perm = n_perm, n_null = n_null, seed = 100L * G + g)
    if (nrow(ph)) spk_rows[[length(spk_rows) + 1L]] <- data.frame(G = G, Regime = g, ph)
    reg_rows[[length(reg_rows) + 1L]] <- data.frame(G = G, Regime = g, Segments = length(g_units),
      SpikesTested = nrow(ph), OrderedSpikes = sum(ph$p_order < 0.05), CoexistSpikes = sum(ph$p_coexist < 0.05),
      p_regime = regime_phase_p(ph))
    cat(sprintf("  G = %2d regime %2d: %3d segments, spikes %d, ordered %d, two-phase %d, p = %.4f\n",
                G, g, length(g_units), nrow(ph), sum(ph$p_order < 0.05), sum(ph$p < 0.05), regime_phase_p(ph)))
  }
  rm(f); invisible(gc())
}
reg <- do.call(rbind, reg_rows)
reg$p_holm <- stats::ave(reg$p_regime, reg$G, FUN = function(p) stats::p.adjust(p, "holm"))
reg$Split <- reg$p_holm < 0.05
utils::write.csv(reg, file.path(int_dir, "phase_regimes.csv"), row.names = FALSE)
if (length(spk_rows)) utils::write.csv(do.call(rbind, spk_rows), file.path(int_dir, "phase_spikes.csv"), row.names = FALSE)
cat("\nRegimes holding a discrete sub-regime (Holm 5% within G):\n")
print(aggregate(cbind(Split, Regimes = 1) ~ G, reg, sum))
cat("PHASE DONE\n")
