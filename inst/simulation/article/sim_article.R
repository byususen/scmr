# Simulation study of the SCMR-EN article.
#
# One method (SCMR-EN: multinomial EN per cluster, Potts prior on block labels
# with phi estimated by pseudo-likelihood, PLIC for G) and two predictor
# specifications; every method in a cell uses the same specification as the
# data-generating process:
#   spec = static   y_it | x_it, g ~ multinomial logit (no previous class)
#   spec = dynamic  y_it | y_i,t-1, x_it, g ~ Markov multinomial logit
#                   (previous class as one-hot predictors; at new locations the
#                   previous class is observed, as in a nowcast)
#
# Experiments
#   exp=S1  core benchmark on the domain of Sugasawa and Murakami (2021):
#           patterns global, grid (2 x 3), irregular (disconnected), smooth GP
#   exp=S4  blocks of 3 x 3 units sharing the regime (like KSA segments):
#           labels per block (proposed) versus labels per unit
#   exp=S3  oracle check of the recovery theorem: misassignment versus the
#           evidence per label (T x block size) and phi; phi-hat on unit and
#           block graphs (degeneracy lemma)
#
# Usage: Rscript sim_article.R exp=S1 spec=static pattern=grid reps=50 rep_start=1 cores=4 out=dir
suppressPackageStartupMessages(library(scmr))
`%||%` <- function(x, y) if (is.null(x)) y else x

opt <- list(exp = "S1", spec = "static", pattern = "grid", reps = 2L, rep_start = 1L, cores = 1L,
            out = "sim_out", n_units = 150L, n_time = 12L, block_size = 1L, G = 3L, delta = 1,
            G_grid = "1,2,3,4,5,6,7,8", p_active = 5L, p_inactive = 10L, coef_k = 10L, n_starts = 6L,
            k_neighbors = 8L, new_share = 0.2, scr = TRUE, gw = TRUE, unit_labels = FALSE,
            multipliers = "0.0625,0.25,1,4")
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2L && kv[1] %in% names(opt)) opt[[kv[1]]] <- utils::type.convert(kv[2], as.is = TRUE)
}
G_grid <- as.integer(strsplit(as.character(opt$G_grid), ",", fixed = TRUE)[[1]])
dynamic <- identical(opt$spec, "dynamic")
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
cat("Options:\n"); str(opt)

ari <- function(a, b) {
  tab <- table(a, b); p2 <- function(x) sum(x * (x - 1) / 2)
  n2 <- p2(sum(tab)); r <- p2(rowSums(tab)); k <- p2(colSums(tab)); e <- r * k / n2
  if (((r + k) / 2 - e) == 0) return(NA_real_)
  (p2(tab) - e) / ((r + k) / 2 - e)
}
metrics <- function(y, prob, classes) {
  pred <- factor(classes[max.col(prob, ties.method = "first")], levels = classes)
  calc_overall_metrics(y, pred, prob, classes)[, c("Accuracy", "Kappa", "LogLoss", "BrierScore", "MacroF1")]
}
pattern_id <- function(p) match(p, c("global", "grid", "irregular", "smooth", "blocks"))

simulate <- function(rep, n_time = opt$n_time, block_size = opt$block_size, pattern = opt$pattern,
                     spec_dynamic = dynamic, seed = NULL) {
  seed <- seed %||% (100000L * (1L + spec_dynamic) + 1000L * rep + 10L * pattern_id(pattern) + (block_size > 1L))
  sim <- simulate_scmr_panel(n_units = opt$n_units, n_time = n_time, pattern = pattern, G = opt$G,
                             p_active = opt$p_active, p_inactive = opt$p_inactive, delta = opt$delta,
                             seed = seed, domain = "scr", grid = c(2, 3), smooth_type = "gp",
                             dynamic = spec_dynamic, block_size = block_size)
  sim$seed <- seed
  sim
}

# Mean squared error of the unit-level active slopes (centred over classes).
unit_beta_mse <- function(f, sim, train_units, label_of_unit) {
  feats <- paste0("x", seq_len(opt$p_active))
  ui <- match(train_units, unique(sim$unit_id))
  truth <- sim$unit_beta[ui, , , drop = FALSE]
  est <- array(0, dim(truth))
  if (inherits(f, "scmr_gw")) {
    uc <- sim$unit_coords[ui, , drop = FALSE]
    near <- apply(uc, 1L, function(s) which.min(colSums((t(f$anchor_coords) - s)^2)))
    for (a in unique(near)) {
      cf <- scmr:::engine_coef(f$engines[[a]], "sum_to_zero")
      est[near == a, , ] <- rep(t(cf[, feats, drop = FALSE]), each = sum(near == a))
    }
  } else {
    cf <- coef(f, "sum_to_zero")
    grp <- if (f$G == 1L) rep(1L, length(train_units)) else if (f$model == "fixed_clusters") {
      sim$unit_regime[ui]
    } else as.integer(f$group_unit[label_of_unit[train_units]])
    for (g in unique(grp)) est[grp == g, , ] <- rep(t(cf[g, , feats, drop = TRUE]), each = sum(grp == g))
  }
  mean((est - truth)^2)
}

run_dataset <- function(rep) {
  sim <- simulate(rep)
  seed <- sim$seed
  pattern <- opt$pattern
  classes <- levels(sim$y)
  X <- if (dynamic) cbind(sim$x, sim$lag_x) else sim$x
  blocks <- unique(sim$block_id)
  units <- unique(sim$unit_id)
  set.seed(seed)
  new_blocks <- sample(blocks, round(opt$new_share * length(blocks)))
  tr <- !sim$block_id %in% new_blocks
  nw <- !tr
  train_units <- unique(sim$unit_id[tr])
  bs <- opt$block_size
  # Proposed labels: one per block (a block of one unit when block_size = 1).
  lab <- sim$block_id; lab_xy <- sim$block_coords
  label_of_unit <- stats::setNames(sim$unit_block, units)
  Gt <- if (pattern == "global") 1L else if (pattern == "smooth") NA_integer_ else length(sim$beta)
  base <- list(lambda_scale = "sum", type_multinomial = "ungrouped", min_units = if (bs > 1L) 4L else 8L,
               min_per_class = 3L, max_iter = 30L, n_starts = as.integer(opt$n_starts),
               coef_init_k = if (bs > 1L) max(2L, round(opt$coef_k / bs * 3)) else as.integer(opt$coef_k),
               coef_init_anchors = 200L, coef_init_alpha = 0.1, coef_init_lambda = 0.0625,
               k_neighbors = as.integer(opt$k_neighbors), tiny_movement_max_units = 0, tiny_movement_rate_tol = 0,
               tiny_movement_revert = FALSE, alpha_grid = 0.5, lambda_rule = "lambda.min")
  ctrl <- function(...) do.call(scmr_control, utils::modifyList(base, list(...)))
  safe <- function(expr) tryCatch(expr, error = function(e) e)
  timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- safe(expr); attr(v, "sec") <- proc.time()[["elapsed"]] - t0; v }
  glob <- timed(fit_scmr(X[tr, ], sim$y[tr], model = "global", alpha = 0.5, lambda_rule = "lambda.min",
                         unit_id = lab[tr], control = ctrl(), seed = seed))
  if (inherits(glob, "error")) stop("global fit failed: ", conditionMessage(glob))
  # Local penalty, as in the application: lambda.1se of the global CV times a
  # multiplier chosen by inner validation (20% of the training blocks held out,
  # two-stage fit at G = true G or 3, log-loss at the held-out blocks).
  cvfit <- glob$fits[[1]]$fit
  lam_base <- if (inherits(cvfit, "cv.glmnet")) cvfit$lambda.1se else glob$lambda_used[1]
  mults <- as.numeric(strsplit(as.character(opt$multipliers), ",", fixed = TRUE)[[1]])
  tr_blocks <- unique(lab[tr])
  set.seed(seed + 7L)
  val_blocks <- sample(tr_blocks, round(0.2 * length(tr_blocks)))
  fi <- tr & !lab %in% val_blocks; vi <- tr & lab %in% val_blocks
  G_tune <- if (!is.na(Gt) && Gt > 1L) Gt else 3L
  tune_tab <- do.call(rbind, lapply(mults, function(mm) {
    f <- safe(fit_scmr(X[fi, ], sim$y[fi], G = G_tune, unit_id = lab[fi], coords = lab_xy[fi, ], alpha = 0.5,
                       lambda = lam_base * mm, control = ctrl(update_memberships = FALSE, n_starts = 1L,
                                                              init_method = "kmeans"), seed = seed))
    ll <- if (inherits(f, "error")) NA_real_ else {
      pv <- predict(f, X[vi, , drop = FALSE], new_unit_id = lab[vi], new_coords = lab_xy[vi, ], membership = "potts")
      metrics(sim$y[vi], pv, classes)$LogLoss
    }
    data.frame(Multiplier = mm, LogLoss = ll)
  }))
  lam <- if (all(is.na(tune_tab$LogLoss))) lam_base else lam_base * tune_tab$Multiplier[which.min(tune_tab$LogLoss)]
  fitG <- function(G, ids = lab, xy = lab_xy, ...) {
    if (G == 1L) return(glob)
    timed(fit_scmr(X[tr, ], sim$y[tr], G = G, unit_id = ids[tr], coords = xy[tr, ], alpha = 0.5,
                   lambda = lam, control = ctrl(...), seed = seed))
  }
  crit_names <- c(PLIC_AIC = "CriterionPLIC_AIC", PLIC_BIC = "CriterionPLIC_BIC",
                  SCR_AIC_eff = "CriterionSCR_AIC_effective", SCR_BIC_eff = "CriterionSCR_BIC_effective",
                  CB_BIC = "CriterionCB_BIC")
  grid_fits <- lapply(G_grid, function(G) fitG(G, phi_update = "pl"))
  sel <- bind(lapply(seq_along(G_grid), function(j) {
    f <- grid_fits[[j]]
    if (inherits(f, "error")) return(data.frame(G = G_grid[j], Error = conditionMessage(f)))
    v <- vapply(crit_names, function(cn) f$criteria[[cn]] %||% NA_real_, numeric(1))
    cbind(data.frame(G = G_grid[j], Error = ""), as.data.frame(as.list(v)),
          data.frame(Phi = f$phi %||% NA_real_, LogPL = f$criteria$PottsLogPseudoLik %||% NA_real_,
                     Gfit = f$G))
  }))
  pick <- function(crit) {
    ok <- is.finite(sel[[crit]] %||% NA)
    if (!any(ok)) NA_integer_ else sel$G[ok][which.min(sel[[crit]][ok])]
  }
  G_sel <- pick("PLIC_AIC")
  G_star <- if (is.na(Gt)) G_sel else Gt
  models <- list(`Global-EN` = glob)
  models$`SCMR-EN` <- if (!is.na(G_sel)) grid_fits[[match(G_sel, G_grid)]] else simpleError("no G selected")
  if (!is.na(Gt) && Gt > 1L) {
    models$`SCMR-EN@trueG` <- if (Gt %in% G_grid) grid_fits[[match(Gt, G_grid)]] else fitG(Gt, phi_update = "pl")
  }
  if (isTRUE(G_star > 1L)) {
    models$`SCMR-EN-phi1` <- fitG(G_star, phi_update = "fixed", phi = 1)
    models$`TwoStage-EN` <- fitG(G_star, update_memberships = FALSE, n_starts = 1L, init_method = "kmeans")
    if (isTRUE(opt$unit_labels) && bs > 1L) {
      models$`SCMR-EN-unitlabels` <- fitG(G_star, ids = sim$unit_id, xy = sim$coords, phi_update = "pl",
                                          min_units = 8L, coef_init_k = as.integer(opt$coef_k))
    }
  }
  if (!is.na(Gt) && Gt > 1L) {
    models$`TCMR-EN (oracle)` <- timed(fit_scmr(X[tr, ], sim$y[tr], model = "fixed_clusters", cluster = sim$regime[tr],
                                                unit_id = lab[tr], alpha = 0.5, lambda = lam, control = ctrl(), seed = seed))
  }
  if (isTRUE(opt$scr)) {
    scr_fit <- function(G) timed(fit_scr_original(X[tr, ], sim$y[tr], lab[tr], lab_xy[tr, ], G = G, phi = 1,
                                                  k_neighbors = 5, seed = seed))
    scr_grid <- lapply(G_grid, scr_fit)
    sel$SCR_BIC_orig <- vapply(scr_grid, function(f) if (inherits(f, "error")) NA_real_ else f$criteria$CriterionSCR_BIC_original, numeric(1))
    sel$SCR_Geff <- vapply(scr_grid, function(f) if (inherits(f, "error")) NA_real_ else f$G, numeric(1))
    G_scr <- pick("SCR_BIC_orig")
    models$`SCR-original` <- if (!is.na(G_scr)) scr_grid[[match(G_scr, G_grid)]] else simpleError("no G selected")
    if (!is.na(Gt) && Gt > 1L) models$`SCR-original@trueG` <- if (Gt %in% G_grid) scr_grid[[match(Gt, G_grid)]] else scr_fit(Gt)
  }
  if (isTRUE(opt$gw)) {
    models$`GW-EN` <- timed(fit_gw_multinom_en(X[tr, ], sim$y[tr], sim$unit_id[tr], sim$coords[tr, ],
                                               k_grid = if (bs > 1L) c(50, 100, 200) else c(25, 50, 100),
                                               alpha = 0.5, lambda = lam, anchors = 80, seed = seed))
  }
  sel$Rep <- rep; sel$TrueG <- Gt
  # Evaluation at the held-out blocks (new locations).
  xn <- X[nw, , drop = FALSE]; yn <- sim$y[nw]
  rows <- list()
  for (nm in names(models)) {
    f <- models[[nm]]
    if (inherits(f, "error")) { rows[[length(rows) + 1L]] <- data.frame(Model = nm, Error = conditionMessage(f)); next }
    is_gw <- inherits(f, "scmr_gw")
    unit_lab <- identical(nm, "SCMR-EN-unitlabels")
    prob <- if (is_gw) predict(f, xn, sim$coords[nw, ]) else if (f$model == "fixed_clusters") {
      predict(f, xn, cluster = sim$regime[nw])
    } else if (f$model == "global" || f$G == 1L) predict(f, xn) else if (unit_lab) {
      predict(f, xn, new_unit_id = sim$unit_id[nw], new_coords = sim$coords[nw, ], membership = "potts")
    } else predict(f, xn, new_unit_id = lab[nw], new_coords = lab_xy[nw, ], membership = "potts")
    m <- metrics(yn, prob, classes)
    beta_mse <- tryCatch(unit_beta_mse(f, sim, train_units, if (unit_lab) stats::setNames(train_units, train_units) else label_of_unit),
                         error = function(e) NA_real_)
    unit_ari <- NA_real_
    if (!is_gw && f$model == "scmr" && f$G > 1L && !is.na(Gt) && Gt > 1L) {
      est <- if (unit_lab) f$group_unit[train_units] else f$group_unit[label_of_unit[train_units]]
      unit_ari <- ari(sim$unit_regime[match(train_units, units)], as.integer(est))
    }
    monotone <- NA
    it <- if (!is_gw) f$diagnostics$iterations else NULL
    if (!is.null(it$PenalizedObjective) && nrow(it)) {
      monotone <- all(tapply(it$PenalizedObjective, it$Start, function(v) all(diff(v) >= -1e-6 * max(1, abs(v)))))
    }
    rows[[length(rows) + 1L]] <- cbind(data.frame(Model = nm, Error = ""), m,
      data.frame(UnitARI = unit_ari, BetaMSE = beta_mse, Phi = if (!is_gw) f$phi %||% NA_real_ else NA_real_,
                 G = if (is_gw) NA_integer_ else f$G, Seconds = attr(f, "sec") %||% NA_real_, Monotone = monotone))
  }
  cols <- unique(unlist(lapply(rows, names)))
  res <- do.call(rbind, lapply(rows, function(r) { for (c in setdiff(cols, names(r))) r[[c]] <- NA; r[cols] }))
  res$Rep <- rep; res$Seed <- seed; res$G_PLIC_AIC <- G_sel; res$TrueG <- Gt
  res$LambdaBase <- lam_base; res$Multiplier <- lam / lam_base
  list(results = res, selection = sel)
}

# S3: oracle assignment error versus evidence per label (T x block size) and
# phi, with the true parameters; and phi-hat at the true labels on the unit
# graph and on the block graph (degeneracy lemma).
theory_check <- function(rep) {
  out <- list()
  for (bs in c(1L, 9L)) {
    sim <- simulate(rep, n_time = 24L, block_size = bs, pattern = "irregular",
                    seed = 700000L + 1000L * rep + 100L * dynamic + bs)
    C <- nlevels(sim$y); G <- length(sim$beta)
    yi <- as.integer(sim$y); li <- as.integer(sim$lag)
    ll <- sapply(seq_len(G), function(g) {
      eta <- (if (dynamic) sim$gamma[[g]][li, ] else matrix(sim$alpha[[g]], length(yi), C, byrow = TRUE)) +
        sim$x[, seq_len(opt$p_active)] %*% sim$beta[[g]]
      eta[cbind(seq_along(yi), yi)] - apply(eta, 1, function(e) max(e) + log(sum(exp(e - max(e)))))
    })
    bl <- unique(sim$block_id)
    bxy <- sim$block_coords[match(bl, sim$block_id), , drop = FALSE]
    btruth <- sim$unit_regime[match(bl, sim$unit_block)]
    wb <- scmr:::build_unit_weights(bxy, k = opt$k_neighbors)$W
    Sb <- scmr:::potts_support(wb, btruth, G)
    wu <- scmr:::build_unit_weights(sim$unit_coords, k = opt$k_neighbors)$W
    phi_unit <- as.numeric(scmr_potts_phi(wu, sim$unit_regime, G, phi_max = 20))
    phi_block <- as.numeric(scmr_potts_phi(wb, btruth, G, phi_max = 20))
    within <- mean(vapply(seq_len(nrow(wu)), function(i) {
      nb <- which(wu[i, ] > 0); mean(sim$unit_block[nb] == sim$unit_block[i])
    }, numeric(1)))
    for (Tn in c(1L, 2L, 4L, 8L, 12L, 24L)) {
      keep <- sim$time <= Tn
      U <- rowsum(ll[keep, , drop = FALSE], sim$block_id[keep])[bl, , drop = FALSE]
      for (phi in c(0, 0.5, 1, 2)) {
        score <- U + phi * Sb
        out[[length(out) + 1L]] <- data.frame(Rep = rep, BlockSize = bs, T = Tn, Evidence = Tn * bs, Phi = phi,
          Misassigned = mean(max.col(score, ties.method = "first") != btruth),
          PhiHatUnitGraph = phi_unit, PhiHatBlockGraph = phi_block, WithinBlockNeighbourShare = within)
      }
    }
  }
  do.call(rbind, out)
}

reps <- seq.int(opt$rep_start, length.out = opt$reps)
t_start <- Sys.time()
par_apply <- function(X, FUN) {
  if (opt$cores > 1L && .Platform$OS.type == "unix") parallel::mclapply(X, FUN, mc.cores = opt$cores, mc.preschedule = FALSE) else lapply(X, FUN)
}
bind <- function(lst) {
  lst <- Filter(function(d) is.data.frame(d) && nrow(d), lst)
  if (!length(lst)) return(data.frame())
  nm <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) { for (c in setdiff(nm, names(d))) d[[c]] <- NA; d[nm] }))
}
tag <- data.frame(Exp = opt$exp, Spec = opt$spec, Pattern = opt$pattern, BlockSize = opt$block_size,
                  NUnits = opt$n_units, NTime = opt$n_time)
if (opt$exp == "S3") {
  th <- bind(par_apply(reps, function(r) tryCatch(theory_check(r), error = function(e) data.frame(Rep = r, Error = conditionMessage(e)))))
  th <- cbind(tag[rep(1L, nrow(th)), c("Exp", "Spec")], th)
  utils::write.csv(th, file.path(opt$out, "theory.csv"), row.names = FALSE)
  print(stats::aggregate(cbind(Misassigned, PhiHatUnitGraph, PhiHatBlockGraph, WithinBlockNeighbourShare) ~ BlockSize + T + Phi,
                         th, mean), digits = 3)
} else {
  out <- par_apply(reps, function(r) tryCatch(run_dataset(r), error = function(e) {
    list(results = data.frame(Model = "ALL", Error = conditionMessage(e), Rep = r), selection = NULL)
  }))
  results <- bind(lapply(out, `[[`, "results"))
  selection <- bind(lapply(out, `[[`, "selection"))
  results <- cbind(tag[rep(1L, nrow(results)), ], results)
  if (nrow(selection)) selection <- cbind(tag[rep(1L, nrow(selection)), ], selection)
  utils::write.csv(results, file.path(opt$out, "results.csv"), row.names = FALSE)
  utils::write.csv(selection, file.path(opt$out, "selection.csv"), row.names = FALSE)
  ok <- results[results$Error == "" & !is.na(results$Kappa), ]
  if (nrow(ok)) {
    cat("\n== Mean over replications (new locations) ==\n")
    print(stats::aggregate(cbind(Kappa, LogLoss, UnitARI, BetaMSE, Phi, G, Seconds) ~ Model, ok, mean, na.action = na.pass), digits = 3)
    cat("\nMonotone objective in every SCMR fit:", all(ok$Monotone[grepl("^SCMR", ok$Model)], na.rm = TRUE), "\n")
  }
  err <- results[results$Error != "", ]
  if (nrow(err)) { cat("\n== Errors ==\n"); print(err[, c("Rep", "Model", "Error")]) }
  if (nrow(selection)) {
    cat("\n== Selected G (true G =", selection$TrueG[1], ") ==\n")
    for (crit in intersect(c("PLIC_AIC", "PLIC_BIC", "SCR_AIC_eff", "SCR_BIC_eff", "CB_BIC", "SCR_BIC_orig"), names(selection))) {
      best <- vapply(split(selection, selection$Rep), function(d) {
        ok <- is.finite(d[[crit]]); if (!any(ok)) NA_integer_ else d$G[ok][which.min(d[[crit]][ok])]
      }, integer(1))
      cat(sprintf("  %-13s", crit)); print(table(factor(best, levels = G_grid)))
    }
  }
}
cat("\nElapsed:", format(Sys.time() - t_start), "\n")
