# =============================================================================
# INTERPRETATION OF THE KSA APPLICATION (after SCMR_EN_empirical_article_v4.R)
#
# 1. Final SCMR-EN-S fits on ALL segments (no hold-out) for the G grid below,
#    with the same tuning as the cross-validation (global CV of the static
#    specification, inner validation of the local penalty, segment labels,
#    Potts phi estimated). PLIC-AIC and PLIC-BIC are reported over the grid.
# 2. Regime outputs per fitted G: segment memberships with soft (Potts)
#    weights, cluster x regency table, monthly phase calendar per regime,
#    standardized coefficients per regime and class, and maps (one panel per
#    regency, colour-blind-safe palette, PNG 300 dpi).
# 3. Segment bootstrap of the paired E2 kappa differences (SCMR-EN-S with G
#    chosen per fold by PLIC-AIC and by PLIC-BIC versus each comparator), using
#    the pooled out-of-fold predictions of the v4 run.
#
# Resume-safe: one RDS per fitted G under <output_dir>/interpret/fits.
# Environment: SCMR_INTERP_G (default "1,3,6,8,10,12,15"), SCMR_INTERP_STAGE
# ("all", "fit" or "report"), SCMR_BOOT (default 1000), SCMR_INTERP_WORKERS
# (default 1; > 1 fits the G values in parallel R worker processes, each
# loading the data once), plus the v4 variables.
# =============================================================================
interp_G <- suppressWarnings(as.integer(strsplit(Sys.getenv("SCMR_INTERP_G", "1,3,6,8,10,12,15"), ",")[[1]]))
interp_G <- sort(unique(interp_G[!is.na(interp_G) & interp_G >= 1]), decreasing = TRUE)
interp_stage <- Sys.getenv("SCMR_INTERP_STAGE", "all")
n_boot <- as.integer(Sys.getenv("SCMR_BOOT", "1000"))
n_workers <- max(1L, as.integer(Sys.getenv("SCMR_INTERP_WORKERS", "1")))

syntax_dir0 <- Sys.getenv("SCMR_SYNTAX_DIR", "C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax")
v4_script <- Sys.getenv("SCMR_V4_SCRIPT", file.path(syntax_dir0, "SCMR_EN_empirical_article_v4.R"))
designs_env <- Sys.getenv("SCMR_DESIGNS")
Sys.setenv(SCMR_DESIGNS = "none")
source(v4_script)                       # data, design matrices and helpers only
Sys.setenv(SCMR_DESIGNS = designs_env)

int_dir <- file.path(output_dir, "interpret")
dir.create(file.path(int_dir, "fits"), recursive = TRUE, showWarnings = FALSE)
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000",
           "#999999", "#882255", "#44AA99", "#117733", "#332288", "#AA4499", "#DDCC77")
regency_names <- c("3205" = "Garut", "3212" = "Indramayu", "6171" = "Pontianak")

# -----------------------------------------------------------------------------
# 1. Final fits on all segments
# -----------------------------------------------------------------------------
all_rows <- seq_len(nrow(dat))
seed_full <- outer_seed + 999L
fit_file <- function(G) file.path(int_dir, "fits", sprintf("SCMR_EN_S_full_G%02d.rds", G))
if (interp_stage %in% c("all", "fit")) {
  banner("FINAL FITS ON ALL SEGMENTS")
  global <- global_for("FULL", 0L, "S", all_rows, seed_full)
  cvfit <- global$fits[[1]]$fit
  lam_base <- if (inherits(cvfit, "cv.glmnet")) cvfit[[local_lambda_rule]] else global$lambda_used[1]
  tuned <- tune_local("FULL", 0L, "S", all_rows, lam_base, global$alpha, seed_full, "direct")
  lam <- lam_base * tuned$multiplier
  features <- local_features("FULL", 0L, "S", all_rows, global, seed_full)
  x_cols <- global$x_colnames
  alpha_full <- global$alpha
  cat("  alpha", alpha_full, "| local lambda", signif(lam, 3), "(x", tuned$multiplier, ")\n")
  # One fit; uses only objects created by the v4 script plus the arguments, so
  # it runs unchanged in a worker that has sourced the v4 script.
  fit_one <- function(G, x_cols, alpha_full, lam, features, seed_full, path) {
    tryCatch({
      cat(sprintf("  [%s] SCMR-EN-S G = %d on all segments\n", format(Sys.time(), "%H:%M"), G))
      t0 <- proc.time()[["elapsed"]]
      fit <- fit_scmr(X_base[, x_cols, drop = FALSE], dat$phase, model = "scmr", G = G, unit_id = seg_id,
                      coords = seg_xy, alpha = alpha_full, lambda = lam, control = control_for("SCMR-EN"),
                      seed = seed_full, init_features = features)
      fit$fit_minutes <- (proc.time()[["elapsed"]] - t0) / 60
      saveRDS(fit, path)
      sprintf("G = %2d done in %.1f min | phi %.3f | PLIC-AIC %.0f | PLIC-BIC %.0f", G, fit$fit_minutes,
              if (is.null(fit$phi)) NA_real_ else fit$phi, fit$criteria$CriterionPLIC_AIC, fit$criteria$CriterionPLIC_BIC)
    }, error = function(e) sprintf("G = %2d FAILED: %s", G, conditionMessage(e)))
  }
  todo <- interp_G[!file.exists(vapply(interp_G, fit_file, ""))]
  for (G in setdiff(interp_G, todo)) cat("  G =", G, "exists\n")
  if (length(todo) && n_workers > 1L && length(todo) > 1L) {
    nw <- min(n_workers, length(todo))
    cat("  fitting G =", paste(todo, collapse = ", "), "in", nw, "worker processes (progress in",
        file.path(int_dir, "workers_log.txt"), ")\n")
    cl <- parallel::makePSOCKcluster(nw, outfile = file.path(int_dir, "workers_log.txt"))
    msgs <- tryCatch({
      parallel::clusterCall(cl, function(script) {
        Sys.setenv(SCMR_DESIGNS = "none")
        sink(tempfile()); on.exit(sink())   # keep the workers' data-loading output out of the log
        source(script)
        TRUE
      }, v4_script)
      paths <- stats::setNames(vapply(todo, fit_file, ""), todo)
      parallel::clusterApplyLB(cl, todo, fun = function(G, one, paths, ...) one(G, ..., path = paths[[as.character(G)]]),
                               one = fit_one, paths = paths, x_cols = x_cols, alpha_full = alpha_full, lam = lam,
                               features = features, seed_full = seed_full)
    }, finally = parallel::stopCluster(cl))
  } else msgs <- NULL
  if (length(todo) && is.null(msgs)) msgs <- lapply(todo, function(G) fit_one(G, x_cols, alpha_full, lam, features, seed_full, fit_file(G)))
  for (m in msgs) cat("   ", m, "\n")
}

# -----------------------------------------------------------------------------
# 2. Regime outputs
# -----------------------------------------------------------------------------
if (interp_stage %in% c("all", "report")) {
  banner("REGIME OUTPUTS")
  loadNamespace("glmnet")                # coef() methods for the saved fits (report-only runs fit nothing)
  fits <- list()
  for (G in interp_G) if (file.exists(fit_file(G))) fits[[as.character(G)]] <- readRDS(fit_file(G))
  if (length(fits)) {
    crit <- do.call(rbind, lapply(fits, function(f) data.frame(G = f$G, Phi = f$phi %||% NA_real_,
      PLIC_AIC = f$criteria$CriterionPLIC_AIC, PLIC_BIC = f$criteria$CriterionPLIC_BIC,
      df_effective = f$criteria$df_effective_total, logPL = f$criteria$PottsLogPseudoLik,
      Minutes = f$fit_minutes %||% NA_real_)))
    crit$Chosen <- trimws(paste(ifelse(crit$PLIC_AIC == min(crit$PLIC_AIC), "PLIC-AIC", ""),
                                ifelse(crit$PLIC_BIC == min(crit$PLIC_BIC), "PLIC-BIC", "")))
    utils::write.csv(crit, file.path(int_dir, "full_fit_criteria.csv"), row.names = FALSE)
    print(crit, row.names = FALSE, digits = 6)

    seg_tab2 <- unique(data.frame(Segment = seg_id, Regency = dat$regency, Lat = seg_xy[, 1], Lon = seg_xy[, 2]))
    x_sd <- apply(X_base, 2, stats::sd)
    # Potts conditional weight of each segment's own regime given the regimes of
    # its k nearest segments (the weight a segment would receive with its own
    # data removed; the posterior weights are ~1 because every segment carries
    # many observations).
    potts_weight <- function(mem, phi, G, k = 8L) {
      D <- as.matrix(stats::dist(as.matrix(mem[, c("Lat", "Lon")]))); diag(D) <- Inf
      nb <- t(apply(D, 1, function(d) order(d)[seq_len(k)]))
      cnt <- t(apply(nb, 1, function(j) tabulate(mem$Cluster[j], G)))
      w <- exp(phi * cnt); w <- w / rowSums(w)
      own <- cbind(seq_len(nrow(mem)), mem$Cluster)
      data.frame(PottsWeight = w[own], NeighbourShare = cnt[own] / k,
                 MajorityAgree = cnt[own] == apply(cnt, 1, max))
    }
    coherence <- list()
    for (f in fits) {
      G <- f$G
      if (G < 2) next                    # G = 1 enters the criterion table only
      tag <- sprintf("G%02d", G)
      # Memberships with soft weights (Potts conditional of the 8 nearest segments).
      pp <- f$PP_unit
      mem <- data.frame(Segment = names(f$group_unit), Cluster = as.integer(f$group_unit),
                        PosteriorWeight = apply(pp[names(f$group_unit), , drop = FALSE], 1, max))
      mem <- merge(seg_tab2, mem, by = "Segment")
      mem <- cbind(mem, potts_weight(mem, f$phi, G))
      coherence[[tag]] <- data.frame(G = G, Phi = f$phi, MajorityAgree = mean(mem$MajorityAgree),
        MixedNeighbours = mean(mem$NeighbourShare < 1), PottsWeight_median = stats::median(mem$PottsWeight),
        PottsWeight_below_half = mean(mem$PottsWeight < 0.5), PosteriorWeight_min = min(mem$PosteriorWeight))
      utils::write.csv(mem, file.path(int_dir, sprintf("memberships_%s.csv", tag)), row.names = FALSE)
      tab <- table(Cluster = mem$Cluster, Regency = regency_names[mem$Regency])
      utils::write.csv(as.data.frame.matrix(tab), file.path(int_dir, sprintf("cluster_by_regency_%s.csv", tag)))
      cat("\n== G =", G, "| clusters x regency (segments)\n"); print(tab)
      # Monthly phase calendar per regime: share of each phase by calendar month.
      cl_row <- mem$Cluster[match(seg_id, mem$Segment)]
      cal <- as.data.frame(prop.table(table(Cluster = cl_row, Month = dat$month, Phase = dat$phase), c(1, 2)))
      names(cal)[4] <- "Share"
      utils::write.csv(cal, file.path(int_dir, sprintf("phase_calendar_%s.csv", tag)), row.names = FALSE)
      # Standardized coefficients (per SD of the predictor), centred over classes.
      cf <- coef(f, "sum_to_zero")
      feats <- dimnames(cf)[[3]][-1]
      long <- do.call(rbind, lapply(seq_len(G), function(g) {
        m <- cf[g, , feats, drop = FALSE][1, , ]
        data.frame(Cluster = g, Class = rep(rownames(m), times = ncol(m)), Feature = rep(colnames(m), each = nrow(m)),
                   Coef = as.vector(m), CoefStd = as.vector(sweep(m, 2, x_sd[colnames(m)], "*")))
      }))
      long <- long[long$Coef != 0, ]
      utils::write.csv(long, file.path(int_dir, sprintf("coefficients_%s.csv", tag)), row.names = FALSE)
      # Top 10 features per regime by the range of standardized effects over classes.
      top <- do.call(rbind, lapply(split(long, long$Cluster), function(z) {
        r <- tapply(z$CoefStd, z$Feature, function(v) diff(range(c(0, v))))
        r <- sort(r, decreasing = TRUE)[seq_len(min(10, length(r)))]
        data.frame(Cluster = z$Cluster[1], Feature = names(r), Range = as.numeric(r))
      }))
      utils::write.csv(top, file.path(int_dir, sprintf("top_features_%s.csv", tag)), row.names = FALSE)
      # Maps: one panel per regency (study areas only), colour-blind-safe palette.
      cols <- okabe[seq_len(G)]
      grDevices::png(file.path(int_dir, sprintf("map_regimes_%s.png", tag)), width = 3300, height = 1500, res = 300)
      graphics::layout(matrix(c(1, 2, 3, 4, 4, 4), 2, byrow = TRUE), heights = c(1, 0.22))
      op <- graphics::par(mar = c(4, 4, 3, 1))
      for (rg in names(regency_names)) {
        z <- mem[mem$Regency == rg, ]
        if (!nrow(z)) { graphics::plot.new(); next }
        graphics::plot(z$Lon, z$Lat, pch = 21, bg = cols[z$Cluster], col = "grey20",
                       cex = 0.5 + 1.5 * z$PottsWeight, asp = 1, xlab = "Longitude", ylab = "Latitude",
                       main = regency_names[[rg]])
      }
      graphics::par(mar = c(0, 0, 0, 0))
      graphics::plot.new()
      graphics::legend("top", legend = paste("Regime", seq_len(G)), pt.bg = cols, pch = 21, horiz = TRUE,
                       bty = "n", cex = 0.9)
      graphics::mtext(paste("Point size: Potts conditional weight of the segment's own regime given its 8 nearest segments.",
                            "Map lines delineate study areas and do not necessarily depict accepted national boundaries."),
                      side = 1, line = -1, cex = 0.55)
      graphics::par(op)
      grDevices::dev.off()
    }
    if (length(coherence)) {
      coh <- do.call(rbind, coherence); rownames(coh) <- NULL
      utils::write.csv(coh, file.path(int_dir, "spatial_coherence.csv"), row.names = FALSE)
      cat("\nSpatial coherence of the fitted regimes (8 nearest segments):\n"); print(coh, row.names = FALSE, digits = 3)
    }
  }

  # ---------------------------------------------------------------------------
  # 3. Segment bootstrap of paired E2 differences
  # ---------------------------------------------------------------------------
  banner("SEGMENT BOOTSTRAP (E2)")
  rp_file <- file.path(output_dir, "row_predictions.csv")
  res_file <- file.path(output_dir, "results.csv")
  if (file.exists(rp_file) && file.exists(res_file)) {
    res <- utils::read.csv(res_file, stringsAsFactors = FALSE)
    rp <- utils::read.csv(rp_file, stringsAsFactors = FALSE, colClasses = c(Segment = "character", Phase = "character"))
    rp <- rp[rp$Design == "E2_mapping", ]
    e2 <- res[res$Design == "E2_mapping", ]
    pick <- function(model, crit) {
      z <- e2[e2$Model == model & e2$G > 1 & is.finite(e2[[crit]]), ]
      do.call(rbind, lapply(split(z, z$Fold), function(d) d[which.min(d[[crit]]), c("Fold", "G")]))
    }
    sel <- list(`SCMR-EN-S (PLIC-AIC)` = cbind(pick("SCMR-EN-S", "CriterionPLIC_AIC"), Model = "SCMR-EN-S"),
                `SCMR-EN-S (PLIC-BIC)` = cbind(pick("SCMR-EN-S", "CriterionPLIC_BIC"), Model = "SCMR-EN-S"),
                `SCR-original-S (BIC asli)` = cbind(pick("SCR-original-S", "CriterionSCR_BIC_original"), Model = "SCR-original-S"))
    fixed <- c("Global-EN-S" = 1L, "Regency-EN-S" = 3L, "GW-EN-S" = 0L, "LightGBM-S" = 0L)
    for (m in names(fixed)) sel[[m]] <- data.frame(Fold = sort(unique(rp$Fold)), G = fixed[[m]], Model = m)
    preds <- lapply(sel, function(s) {
      z <- merge(rp, s, by = c("Fold", "Model", "G"))
      z[order(z$Row), c("Row", "Segment", "Phase", "Pred")]
    })
    ok <- vapply(preds, nrow, integer(1)) > 0
    preds <- preds[ok]
    rows <- Reduce(intersect, lapply(preds, `[[`, "Row"))
    preds <- lapply(preds, function(z) z[match(rows, z$Row), ])
    truth <- factor(preds[[1]]$Phase, levels = class_levels)
    seg <- preds[[1]]$Segment
    kappa_of <- function(pred, idx) {
      tb <- table(factor(truth[idx], levels = class_levels), factor(pred[idx], levels = class_levels))
      n <- sum(tb); po <- sum(diag(tb)) / n; pe <- sum(rowSums(tb) * colSums(tb)) / n^2
      (po - pe) / (1 - pe)
    }
    pred_chr <- lapply(preds, function(z) as.character(z$Pred))
    seg_rows <- split(seq_along(seg), seg)
    point <- vapply(pred_chr, kappa_of, numeric(1), idx = seq_along(seg))
    set.seed(outer_seed + 77L)
    boot <- t(vapply(seq_len(n_boot), function(b) {
      idx <- unlist(seg_rows[sample.int(length(seg_rows), replace = TRUE)], use.names = FALSE)
      vapply(pred_chr, kappa_of, numeric(1), idx = idx)
    }, numeric(length(pred_chr))))
    out <- list()
    for (ref in intersect(c("SCMR-EN-S (PLIC-AIC)", "SCMR-EN-S (PLIC-BIC)"), names(pred_chr))) {
      for (m in setdiff(names(pred_chr), ref)) {
        d <- boot[, ref] - boot[, m]
        out[[length(out) + 1L]] <- data.frame(Reference = ref, Comparator = m,
          Kappa_reference = point[[ref]], Kappa_comparator = point[[m]], Difference = point[[ref]] - point[[m]],
          CI_low = stats::quantile(d, 0.025), CI_high = stats::quantile(d, 0.975), Share_positive = mean(d > 0))
      }
    }
    bt <- do.call(rbind, out); rownames(bt) <- NULL
    utils::write.csv(bt, file.path(int_dir, "bootstrap_E2_pooled_kappa.csv"), row.names = FALSE)
    cat("Pooled out-of-fold kappa and 95% segment-bootstrap intervals of differences (", n_boot, "resamples ):\n")
    print(bt, row.names = FALSE, digits = 3)
  } else cat("row_predictions.csv or results.csv not found in", output_dir, "\n")
}
banner("INTERPRETATION DONE")
cat("Outputs in:", int_dir, "\n")
