# =============================================================================
# RE-SCORING THE FULL-DATA KSA FITS WITH SIZE-AWARE PLIC (no refitting)
#
# Reads the saved fits interpret/fits/SCMR_EN_S_full_GXX.rds written by
# SCMR_EN_interpret_v4.R and recomputes, from the per-regime components already
# stored in fit$criteria$cluster_components:
#
#   PLIC-BIC (current)  = -2 l - 2 logPL + log(N)   * (sum_g df_g + d_phi)
#   PLIC-BIC-rows       = -2 l - 2 logPL + sum_g df_g * log(n_g) + d_phi * log(B)
#   PLIC-BIC-blocks     = -2 l - 2 logPL + sum_g df_g * log(m_g) + d_phi * log(B)
#   PLIC-AIC            = -2 l - 2 logPL + 2 * (sum_g df_g + d_phi)
#
# n_g = rows (segment-months) in regime g, m_g = segments in regime g,
# N = rows, B = segments, d_phi = 1 when phi is estimated.
# The size-aware forms follow from a regime-wise Laplace approximation: each
# regime's coefficients are informed only by that regime's own observations.
#
# Output (aggregates only; no segment identifiers or coordinates):
#   interpret/size_aware_plic.csv        one row per G
#   interpret/size_aware_regime_parts.csv one row per G x regime
#   interpret/size_aware_starts.csv      multi-start objectives per G
#
# Run in RStudio:  source(".../Syntax/SCMR_EN_rescore_size_aware.R")
# Takes a few seconds per G (only reads the RDS files).
# =============================================================================
results_dir <- Sys.getenv("SCMR_OUTPUT_DIR",
  "C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Results/Real data v4")
int_dir <- file.path(results_dir, "interpret")
fit_files <- sort(list.files(file.path(int_dir, "fits"), pattern = "^SCMR_EN_S_full_G[0-9]+\\.rds$", full.names = TRUE))
if (!length(fit_files)) stop("No fits found in ", file.path(int_dir, "fits"))
suppressWarnings(loadNamespace("glmnet"))

rows <- list(); parts_all <- list(); starts_all <- list()
for (path in fit_files) {
  f <- readRDS(path)
  cr <- f$criteria
  pc <- cr$cluster_components
  G <- if (!is.null(f$G)) f$G else nrow(pc)
  N <- cr$n_rows; B <- cr$n_units
  ll <- cr$logLik; logpl <- cr$PottsLogPseudoLik; d_phi <- cr$df_label
  df_g <- pc$DFEffective_g; n_g <- pc$NRowsCluster; m_g <- pc$NUnitsCluster
  base <- -2 * ll - 2 * logpl
  rows[[length(rows) + 1L]] <- data.frame(
    G = G, Phi = if (is.null(f$phi)) NA_real_ else f$phi,
    N = N, B = B, Neg2LogLik = -2 * ll, Neg2LogPL = -2 * logpl,
    DF_total = sum(df_g), d_phi = d_phi,
    MinRegimeSegments = min(m_g), MaxRegimeSegments = max(m_g),
    Pen_BIC_N = log(N) * (sum(df_g) + d_phi),
    Pen_BIC_rows = sum(df_g * log(n_g)) + d_phi * log(B),
    Pen_BIC_blocks = sum(df_g * log(m_g)) + d_phi * log(B),
    PLIC_AIC = base + 2 * (sum(df_g) + d_phi),
    PLIC_BIC = base + log(N) * (sum(df_g) + d_phi),
    PLIC_BIC_rows = base + sum(df_g * log(n_g)) + d_phi * log(B),
    PLIC_BIC_blocks = base + sum(df_g * log(m_g)) + d_phi * log(B),
    Minutes = if (is.null(f$fit_minutes)) NA_real_ else f$fit_minutes)
  parts_all[[length(parts_all) + 1L]] <- data.frame(
    G = G, Regime = pc$Cluster, Segments = m_g, Rows = n_g,
    ShareSegments = round(m_g / B, 4), DFEffective = df_g,
    Neg2LogLik = pc$ClusterNeg2LogLik, Neg2LogLikPerRow = pc$ClusterNeg2LogLik / n_g,
    ActivePredictors = pc$ActivePredictors_g, Lambda = pc$LambdaUsed_g)
  st <- f$diagnostics$starts
  if (!is.null(st)) starts_all[[length(starts_all) + 1L]] <- data.frame(G = G, st, row.names = NULL)
  rm(f); invisible(gc())
}
tab <- do.call(rbind, rows); tab <- tab[order(tab$G), ]
for (k in c("PLIC_AIC", "PLIC_BIC", "PLIC_BIC_rows", "PLIC_BIC_blocks"))
  tab[[paste0("Delta_", k)]] <- tab[[k]] - min(tab[[k]])
chosen <- vapply(c("PLIC_AIC", "PLIC_BIC", "PLIC_BIC_rows", "PLIC_BIC_blocks"),
                 function(k) tab$G[which.min(tab[[k]])], numeric(1))
parts <- do.call(rbind, parts_all)
utils::write.csv(tab, file.path(int_dir, "size_aware_plic.csv"), row.names = FALSE)
utils::write.csv(parts, file.path(int_dir, "size_aware_regime_parts.csv"), row.names = FALSE)
if (length(starts_all)) utils::write.csv(do.call(rbind, starts_all), file.path(int_dir, "size_aware_starts.csv"), row.names = FALSE)

cat("\n==== Size-aware PLIC on the full-data fits ====\n")
print(tab[, c("G", "Phi", "MinRegimeSegments", "DF_total", "PLIC_AIC", "PLIC_BIC", "PLIC_BIC_rows", "PLIC_BIC_blocks")],
      row.names = FALSE, digits = 7)
cat("\nChosen G:\n"); print(chosen)
cat("\nFiles written to", int_dir, ":\n  size_aware_plic.csv, size_aware_regime_parts.csv, size_aware_starts.csv\n")
cat("RESCORE DONE\n")
