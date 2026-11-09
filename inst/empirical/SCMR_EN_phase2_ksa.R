# =============================================================================
# PHASE-COUNT AND PROJECTED-DIP CHECK ON THE FULL-DATA KSA FITS (after SCMR_EN_spectrum_ksa.R)
#
# Uses the spike counts already in interpret/spectrum_regimes.csv (no bootstrap
# is repeated). The spikes whose segment loadings are spatially ordered (Moran's
# I, 9,999 permutations) are analysed jointly:
#   phase count : k-means clusters of the ordered loadings against continuous
#                 ordered (SAR) fields with matched Moran's I;
#   projected dip: Hartigan's dip test of multimodality of each ordered loading
#                 and each principal axis of the ordered loadings.
# p-values are Holm-adjusted over the regimes of each G. Aggregates only.
#
# Run in RStudio (a few minutes):
#   source("C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax/SCMR_EN_phase2_ksa.R")
# Output: interpret/phase2_regimes.csv, interpret/phase2_spikes.csv
# =============================================================================
syntax_dir0 <- Sys.getenv("SCMR_SYNTAX_DIR", "C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax")
v4_script <- Sys.getenv("SCMR_V4_SCRIPT", file.path(syntax_dir0, "SCMR_EN_empirical_article_v4.R"))
spec_funs <- Sys.getenv("SCMR_SPECTRUM_FUNS", file.path(syntax_dir0, "spectrum_funs.R"))
n_perm <- as.integer(Sys.getenv("SCMR_PHASE_PERM", "9999"))
n_null <- as.integer(Sys.getenv("SCMR_PHASE_NULL", "199"))
if (!requireNamespace("diptest", quietly = TRUE)) install.packages("diptest")
designs_env <- Sys.getenv("SCMR_DESIGNS")
Sys.setenv(SCMR_DESIGNS = "none")
source(v4_script)
Sys.setenv(SCMR_DESIGNS = designs_env)
source(spec_funs)
suppressWarnings(loadNamespace("glmnet"))

int_dir <- file.path(output_dir, "interpret")
spec <- utils::read.csv(file.path(int_dir, "spectrum_regimes.csv"))
banner("PHASE-COUNT AND PROJECTED-DIP CHECK")
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
    A_g <- knn_adj(uxy[g_units, , drop = FALSE])
    pc <- phase_count_test(si$U, si$n_b, si$J, A_g, sg$Spikes[r], n_perm = n_perm, n_null = n_null, seed = 100L * G + g)
    pd <- phase_dip_test(si$U, si$n_b, si$J, A_g, sg$Spikes[r], n_perm = n_perm, seed = 100L * G + g)
    if (nrow(pd$spikes)) spk_rows[[length(spk_rows) + 1L]] <- data.frame(G = G, Regime = g, pd$spikes)
    reg_rows[[length(reg_rows) + 1L]] <- data.frame(G = G, Regime = g, Segments = length(g_units),
      SpikesTested = nrow(pd$spikes), OrderedSpikes = pd$r, p_count = pc$p, Phases = pc$phases,
      p_dip = pd$p, DipProjections = pd$n_proj)
    cat(sprintf("  G = %2d regime %2d: %3d segments, spikes %d, ordered %d, phases %d (p %.4f), dip p %.4f\n",
                G, g, length(g_units), nrow(pd$spikes), pd$r, pc$phases, pc$p, pd$p))
  }
  rm(f); invisible(gc())
}
reg <- do.call(rbind, reg_rows)
holm <- function(p) stats::ave(p, reg$G, FUN = function(q) stats::p.adjust(q, "holm"))
reg$p_count_holm <- holm(reg$p_count); reg$p_dip_holm <- holm(reg$p_dip)
reg$SplitCount <- reg$p_count_holm < 0.05; reg$SplitDip <- reg$p_dip_holm < 0.05
utils::write.csv(reg, file.path(int_dir, "phase2_regimes.csv"), row.names = FALSE)
if (length(spk_rows)) utils::write.csv(do.call(rbind, spk_rows), file.path(int_dir, "phase2_spikes.csv"), row.names = FALSE)
cat("\nRegimes holding discrete sub-regimes (Holm 5% within G):\n")
print(aggregate(cbind(SplitCount, SplitDip, Regimes = 1) ~ G, reg, sum))
cat("PHASE2 DONE\n")
