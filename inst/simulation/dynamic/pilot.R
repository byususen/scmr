# Dynamic SCMR-EN (DSCMR-EN) pilot study on the panel Markov generator.
#
# Usage: Rscript pilot.R reps=10 n_units=150 n_time=24 patterns=irregular,blocks,global \
#                         G=3 delta=1 out=dynamic_pilot cores=4 theory=TRUE
#
# Per dataset: 80% of units train, 20% are new locations. Models:
#   DSCMR-EN       lag state, estimated phi (pseudo-likelihood), G = true (and G grid for selection)
#   DSCMR-EN-phi1  same with phi fixed at 1 (ablation)
#   TwoStage-EN    k-means partition of the coordinates, no membership updates
#   SCMR-EN-static no lag state, estimated phi
#   Global-EN      lag state, one model
#   GW-EN          geographically weighted multinomial EN with lag state
#   TCMR-EN        true partition (oracle)
# New-location prediction modes:
#   lag_observed   previous class known (nowcast at survey points)
#   filter         previous class unknown, exact forward filter
#   plugin         previous class unknown, argmax of the previous prediction
suppressPackageStartupMessages(library(scmr))

opt <- list(reps = 3L, n_units = 150L, n_time = 24L, patterns = "irregular,blocks,global",
            G = 3L, delta = 1, out = "dynamic_pilot", cores = 1L, theory = TRUE,
            G_grid = "1,2,3,4,5", p_active = 5L, p_inactive = 10L, coef_k = 30L, n_starts = 3L,
            domain = "square", smooth_type = "linear", scr = TRUE, membership = "potts")
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2L && kv[1] %in% names(opt)) opt[[kv[1]]] <- utils::type.convert(kv[2], as.is = TRUE)
}
patterns <- strsplit(opt$patterns, ",", fixed = TRUE)[[1]]
G_grid <- as.integer(strsplit(as.character(opt$G_grid), ",", fixed = TRUE)[[1]])
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
  m <- calc_overall_metrics(y, pred, prob, classes)
  m[, c("Accuracy", "Kappa", "LogLoss", "BrierScore", "MacroF1")]
}

# Generic forward filter / plug-in for any row-wise probability function.
set_lag <- function(x, lagc, k) { x[, lagc] <- 0; if (length(k) == 1L) x[, lagc[k]] <- 1 else x[cbind(seq_len(nrow(x)), match(lagc[k], colnames(x)))] <- 1; x }
generic_filter <- function(prob_fun, x, unit, time, lagc, a0) {
  C <- length(lagc)
  P <- lapply(seq_len(C), function(k) prob_fun(set_lag(x, lagc, k)))
  out <- matrix(0, nrow(x), C)
  for (u in unique(unit)) {
    idx <- which(unit == u); idx <- idx[order(time[idx])]
    a <- a0
    for (r in idx) {
      a <- vapply(seq_len(C), function(c) sum(vapply(seq_len(C), function(k) P[[k]][r, c], 0) * a), 0)
      a <- a / sum(a); out[r, ] <- a
    }
  }
  out
}
generic_plugin <- function(prob_fun, x, unit, time, lagc, start_class) {
  out <- matrix(0, nrow(x), length(lagc))
  steps <- sort(unique(time))
  cur <- setNames(rep(start_class, length(unique(unit))), unique(unit))
  for (t in steps) {
    idx <- which(time == t)
    xt <- x
    xt[idx, ] <- set_lag(x[idx, , drop = FALSE], lagc, cur[unit[idx]])
    p <- prob_fun(xt)[idx, , drop = FALSE]
    out[idx, ] <- p
    cur[unit[idx]] <- max.col(p, ties.method = "first")
  }
  out
}

# Mean squared error of the unit-level active slopes (centred over classes)
# at the training units, the main accuracy measure of the SCR paper.
unit_beta_mse <- function(f, sim, train_units) {
  pa <- opt$p_active
  feats <- paste0("x", seq_len(pa))
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
    grp <- if (f$G == 1L) rep(1L, length(train_units)) else as.integer(f$group_unit[train_units])
    for (g in unique(grp)) {
      est[grp == g, , ] <- rep(t(cf[g, , feats, drop = TRUE]), each = sum(grp == g))
    }
  }
  mean((est - truth)^2)
}

run_dataset <- function(pattern, rep) {
  seed <- 1000L * rep + match(pattern, c("irregular", "blocks", "global", "smooth", "grid"))
  sim <- simulate_scmr_panel(n_units = opt$n_units, n_time = opt$n_time, pattern = pattern,
                             G = opt$G, p_active = opt$p_active, p_inactive = opt$p_inactive,
                             delta = opt$delta, seed = seed, domain = opt$domain,
                             smooth_type = opt$smooth_type)
  classes <- levels(sim$y)
  lagc <- colnames(sim$lag_x)
  X <- cbind(sim$x, sim$lag_x)
  units <- unique(sim$unit_id)
  set.seed(seed)
  new_units <- sample(units, round(0.2 * length(units)))
  tr <- !sim$unit_id %in% new_units
  nw <- !tr
  Gt <- if (pattern == "global") 1L else if (pattern == "smooth") opt$G else length(sim$beta)
  base <- list(lambda_scale = "sum", type_multinomial = "ungrouped", min_units = 8, min_per_class = 3,
               max_iter = 30L, n_starts = as.integer(opt$n_starts), coef_init_k = as.integer(opt$coef_k), coef_init_anchors = 200L,
               coef_init_alpha = 0.1, coef_init_lambda = 0.0625, k_neighbors = 8,
               tiny_movement_max_units = 0, tiny_movement_rate_tol = 0, tiny_movement_revert = FALSE,
               alpha_grid = 0.5, lambda_rule = "lambda.min")
  ctrl <- function(...) do.call(scmr_control, utils::modifyList(base, list(...)))
  t0 <- proc.time()[["elapsed"]]
  glob <- fit_scmr(X[tr, ], sim$y[tr], model = "global", alpha = 0.5, lambda_rule = "lambda.min",
                   control = ctrl(), seed = seed)
  lam <- glob$lambda_used[1]
  fitG <- function(G, x = X, ...) {
    if (G == 1L) return(glob)
    fit_scmr(x[tr, ], sim$y[tr], G = G, unit_id = sim$unit_id[tr], coords = sim$coords[tr, ],
             alpha = 0.5, lambda = lam, control = ctrl(...), seed = seed)
  }
  # G grid for selection with the proposed model.
  grid_fits <- lapply(G_grid, function(G) tryCatch(fitG(G, phi_update = "pl"), error = function(e) e))
  sel <- do.call(rbind, lapply(seq_along(G_grid), function(j) {
    f <- grid_fits[[j]]
    if (inherits(f, "error")) return(data.frame(G = G_grid[j], PLIC_AIC = NA, PLIC_BIC = NA, SCR_AIC_eff = NA, Phi = NA, Error = conditionMessage(f)))
    data.frame(G = G_grid[j], PLIC_AIC = f$criteria$CriterionPLIC_AIC, PLIC_BIC = f$criteria$CriterionPLIC_BIC,
               SCR_AIC_eff = f$criteria$CriterionSCR_AIC_effective, Phi = f$phi %||% NA_real_, Error = "")
  }))
  sel$Pattern <- pattern; sel$Rep <- rep; sel$TrueG <- if (pattern == "smooth") NA else Gt
  main <- if (Gt %in% G_grid) grid_fits[[match(Gt, G_grid)]] else tryCatch(fitG(Gt, phi_update = "pl"), error = function(e) e)
  # Original SCR (unpenalized, simultaneous updates, phi = 1, k-means start):
  # G chosen by its own BIC over the same grid, and also fitted at the true G.
  scr_fit <- function(G) tryCatch(fit_scr_original(X[tr, ], sim$y[tr], sim$unit_id[tr], sim$coords[tr, ], G = G,
                                                   phi = 1, k_neighbors = 5, seed = seed), error = function(e) e)
  if (isTRUE(opt$scr)) {
    scr_grid <- lapply(G_grid, scr_fit)
    sel$SCR_orig_BIC <- vapply(scr_grid, function(f) if (inherits(f, "error")) NA_real_ else f$criteria$CriterionSCR_BIC_original, numeric(1))
    sel$SCR_orig_Geff <- vapply(scr_grid, function(f) if (inherits(f, "error")) NA_real_ else f$G, numeric(1))
  }
  models <- list(`DSCMR-EN` = main, `Global-EN` = glob)
  if (isTRUE(opt$scr)) {
    models$`SCR-original` <- if (Gt %in% G_grid) scr_grid[[match(Gt, G_grid)]] else scr_fit(Gt)
  }
  if (Gt > 1L) {
    models$`DSCMR-EN-phi1` <- tryCatch(fitG(Gt, phi_update = "fixed", phi = 1), error = function(e) e)
    models$`TwoStage-EN` <- tryCatch(fitG(Gt, update_memberships = FALSE, n_starts = 1L, init_method = "kmeans"), error = function(e) e)
    models$`SCMR-EN-static` <- tryCatch(fitG(Gt, x = sim$x, phi_update = "pl"), error = function(e) e)
    if (pattern != "smooth") {
      models$`TCMR-EN` <- tryCatch(fit_scmr(X[tr, ], sim$y[tr], model = "fixed_clusters", cluster = sim$regime[tr],
                                            unit_id = sim$unit_id[tr], alpha = 0.5, lambda = lam,
                                            control = ctrl(), seed = seed), error = function(e) e)
    }
  }
  models$`GW-EN` <- tryCatch(fit_gw_multinom_en(X[tr, ], sim$y[tr], sim$unit_id[tr], sim$coords[tr, ],
                                                k_grid = c(25, 50, 100), alpha = 0.5, lambda = lam,
                                                anchors = 80, seed = seed), error = function(e) e)
  fit_seconds <- proc.time()[["elapsed"]] - t0
  a0 <- as.numeric(table(sim$y[tr])) / sum(tr)
  start_class <- which.max(a0)
  xn <- X[nw, ]; un <- sim$unit_id[nw]; tn <- sim$time[nw]; cn <- sim$coords[nw, ]; yn <- sim$y[nw]
  rows <- list()
  monotone <- NA
  for (nm in names(models)) {
    f <- models[[nm]]
    if (inherits(f, "error")) { rows[[length(rows) + 1L]] <- data.frame(Model = nm, Mode = "error", Error = conditionMessage(f)); next }
    static <- nm == "SCMR-EN-static"
    is_gw <- inherits(f, "scmr_gw")
    prob_fun <- if (is_gw) function(x) predict(f, x, cn) else if (f$model == "fixed_clusters") {
      function(x) predict(f, x, cluster = sim$regime[nw])
    } else if (static) function(x) predict(f, x[, colnames(sim$x)], new_unit_id = un, new_coords = cn, membership = opt$membership)
    else function(x) predict(f, x, new_unit_id = un, new_coords = cn, membership = opt$membership)
    preds <- list(lag_observed = prob_fun(xn))
    if (!static) {
      preds$filter <- if (!is_gw && f$model %in% c("scmr", "global")) {
        scmr_filter_predict(f, xn, un, tn, cn, lagc, membership = opt$membership)
      } else generic_filter(prob_fun, xn, un, tn, lagc, a0)
      preds$plugin <- generic_plugin(prob_fun, xn, un, tn, lagc, start_class)
    }
    beta_mse <- tryCatch(unit_beta_mse(f, sim, units[!units %in% new_units]), error = function(e) NA_real_)
    train_ari <- if (!is_gw && f$model == "scmr" && f$G > 1L && pattern %in% c("irregular", "blocks", "grid")) {
      ari(sim$unit_regime[match(f$unit_levels, units)], f$group_unit)
    } else NA_real_
    if (nm == "DSCMR-EN" && !is.null(f$diagnostics$iterations$PenalizedObjective)) {
      tr_obj <- f$diagnostics$iterations
      monotone <- all(tapply(tr_obj$PenalizedObjective, tr_obj$Start, function(v) all(diff(v) >= -1e-6 * max(1, abs(v)))))
    }
    for (md in names(preds)) {
      m <- metrics(yn, preds[[md]], classes)
      rows[[length(rows) + 1L]] <- cbind(data.frame(Model = nm, Mode = md, Error = ""), m,
        data.frame(TrainARI = train_ari, BetaMSE = beta_mse, Phi = if (!is_gw) f$phi %||% NA_real_ else NA_real_,
                   G = if (is_gw) NA_integer_ else f$G))
    }
  }
  res <- do.call(rbind, lapply(rows, function(r) { for (c in setdiff(c("Accuracy", "Kappa", "LogLoss", "BrierScore", "MacroF1", "TrainARI", "BetaMSE", "Phi", "G"), names(r))) r[[c]] <- NA; r }))
  res$Pattern <- pattern; res$Rep <- rep; res$Seed <- seed; res$FitSeconds <- fit_seconds
  res$MonotoneDSCMR <- monotone
  list(results = res, selection = sel)
}

# Theory check (Theorem 4): oracle assignment error versus T and phi.
theory_check <- function(pattern, rep, Tmax = 24L) {
  seed <- 5000L * rep + match(pattern, c("irregular", "blocks"))
  sim <- simulate_scmr_panel(n_units = opt$n_units, n_time = Tmax, pattern = pattern, G = opt$G,
                             p_active = opt$p_active, p_inactive = opt$p_inactive, delta = opt$delta, seed = seed)
  C <- nlevels(sim$y); G <- length(sim$beta)
  yi <- as.integer(sim$y); li <- as.integer(sim$lag)
  # Row log-likelihood under each regime with the true parameters.
  ll <- sapply(seq_len(G), function(g) {
    eta <- sim$gamma[[g]][li, ] + sim$x[, seq_len(opt$p_active)] %*% sim$beta[[g]]
    eta[cbind(seq_along(yi), yi)] - apply(eta, 1, function(e) max(e) + log(sum(exp(e - max(e)))))
  })
  w <- scmr:::build_unit_weights(sim$unit_coords, k = 8)$W
  truth <- sim$unit_regime
  S <- scmr:::potts_support(w, truth, G)
  out <- list()
  for (Tn in c(1L, 2L, 4L, 8L, 12L, 24L)) {
    keep <- sim$time <= Tn
    U <- rowsum(ll[keep, , drop = FALSE], sim$unit_id[keep])[unique(sim$unit_id), , drop = FALSE]
    for (phi in c(0, 0.5, 1, 2)) {
      score <- U + phi * S
      out[[length(out) + 1L]] <- data.frame(Pattern = pattern, Rep = rep, T = Tn, Phi = phi,
                                            Misassigned = mean(max.col(score, ties.method = "first") != truth))
    }
  }
  do.call(rbind, out)
}

tasks <- expand.grid(rep = seq_len(opt$reps), pattern = patterns, stringsAsFactors = FALSE)
runner <- function(i) tryCatch(run_dataset(tasks$pattern[i], tasks$rep[i]),
                               error = function(e) list(results = data.frame(Model = "ALL", Mode = "error",
                                 Error = conditionMessage(e), Pattern = tasks$pattern[i], Rep = tasks$rep[i]),
                                 selection = NULL))
t_start <- Sys.time()
out <- if (opt$cores > 1L && .Platform$OS.type == "unix") parallel::mclapply(seq_len(nrow(tasks)), runner, mc.cores = opt$cores) else lapply(seq_len(nrow(tasks)), runner)
bind <- function(lst) { nm <- unique(unlist(lapply(lst, names))); do.call(rbind, lapply(lst, function(d) { for (c in setdiff(nm, names(d))) d[[c]] <- NA; d[nm] })) }
results <- bind(lapply(out, `[[`, "results"))
selection <- bind(Filter(Negate(is.null), lapply(out, `[[`, "selection")))
utils::write.csv(results, file.path(opt$out, "results.csv"), row.names = FALSE)
utils::write.csv(selection, file.path(opt$out, "selection.csv"), row.names = FALSE)
if (isTRUE(opt$theory)) {
  th <- do.call(rbind, lapply(intersect(patterns, c("irregular", "blocks")), function(p) do.call(rbind, lapply(seq_len(opt$reps), function(r) theory_check(p, r)))))
  utils::write.csv(th, file.path(opt$out, "theory.csv"), row.names = FALSE)
  cat("\n== Oracle misassignment rate by T and phi ==\n")
  print(stats::aggregate(Misassigned ~ Pattern + T + Phi, th, mean))
}
cat("\nElapsed:", format(Sys.time() - t_start), "\n")
ok <- results[results$Mode != "error" & !is.na(results$Kappa), ]
cat("\n== Mean metrics at new locations ==\n")
print(stats::aggregate(cbind(Kappa, LogLoss, TrainARI, BetaMSE) ~ Pattern + Model + Mode, ok, mean, na.action = na.pass), digits = 3)
err <- results[results$Mode == "error", ]
if (nrow(err)) { cat("\n== Errors ==\n"); print(err[, c("Pattern", "Rep", "Model", "Error")]) }
if (nrow(selection)) {
  cat("\n== Selected G ==\n")
  for (crit in intersect(c("PLIC_AIC", "PLIC_BIC", "SCR_AIC_eff", "SCR_orig_BIC"), names(selection))) {
    best <- do.call(rbind, lapply(split(selection, list(selection$Pattern, selection$Rep), drop = TRUE), function(d) {
      d <- d[is.finite(d[[crit]]), ]; if (!nrow(d)) return(NULL)
      data.frame(Pattern = d$Pattern[1], Rep = d$Rep[1], Selected = d$G[which.min(d[[crit]])], TrueG = d$TrueG[1])
    }))
    cat(crit, "\n"); print(table(best$Pattern, best$Selected))
  }
}
cat("\nDSCMR-EN monotone in every dataset:", all(ok$MonotoneDSCMR[ok$Model == "DSCMR-EN"], na.rm = TRUE), "\n")
