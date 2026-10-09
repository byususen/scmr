dyn_fixture <- function(C = 3, n_units = 40, n_time = 4, pattern = "blocks", G = 2, seed = 7) {
  sim <- simulate_scmr_panel(n_units = n_units, n_time = n_time, pattern = pattern, G = G,
                             n_classes = C, p_active = 2, p_inactive = 1, delta = 1.5, seed = seed)
  x <- cbind(sim$x, sim$lag_x)
  c(sim, list(xx = x, lag_cols = colnames(sim$lag_x)))
}

test_that("lag design marks the previous class of the same unit", {
  y <- factor(c("a", "b", "b", "a", "c"), levels = c("a", "b", "c"))
  d <- scmr_lag_design(y, unit_id = c(1, 1, 1, 2, 2), time = c(1, 2, 4, 1, 2))
  expect_identical(d$available, c(FALSE, TRUE, FALSE, FALSE, TRUE))
  expect_equal(unname(d$x[2, ]), c(1, 0, 0))
  expect_equal(unname(d$x[5, ]), c(1, 0, 0))
  expect_true(all(is.na(d$x[1, ])))
  expect_identical(d$columns, c("lag_a", "lag_b", "lag_c"))
})

test_that("local change of the Potts pseudo-likelihood is exact", {
  set.seed(3)
  coords <- cbind(runif(30), runif(30))
  w <- scmr:::build_unit_weights(coords, k = 4)$W
  g <- sample(1:3, 30, replace = TRUE)
  S <- scmr:::potts_support(w, g, 3)
  nb <- scmr:::potts_neighbors(w)
  for (u in c(1, 7, 19)) for (h in 1:3) {
    g2 <- g; g2[u] <- h
    brute <- scmr:::potts_log_pseudolikelihood(w, g2, 0.8, 3) - scmr:::potts_log_pseudolikelihood(w, g, 0.8, 3)
    expect_equal(scmr:::potts_delta_logpl(u, h, g, S, nb, 0.8), brute, tolerance = 1e-10)
  }
  S2 <- scmr:::potts_update_support(S, 7, g[7], 1, nb)
  g3 <- g; g3[7] <- 1
  expect_equal(S2, scmr:::potts_support(w, g3, 3))
})

test_that("phi estimate maximizes the concave pseudo-likelihood", {
  set.seed(4)
  coords <- cbind(runif(80), runif(80))
  w <- scmr:::build_unit_weights(coords, k = 5)$W
  g <- ifelse(coords[, 1] + 0.3 * rnorm(80) < 0.5, 1L, 2L)
  est <- scmr_potts_phi(w, g, 2, phi_start = 0.1)
  f <- function(phi) scmr:::potts_log_pseudolikelihood(w, g, phi, 2)
  ref <- optimize(f, c(0, 20), maximum = TRUE, tol = 1e-10)$maximum
  expect_equal(as.numeric(est), ref, tolerance = 1e-5)
  grid <- seq(0, 5, by = 0.25)
  vals <- vapply(grid, f, numeric(1))
  expect_true(all(diff(diff(vals)) <= 1e-9))
})

test_that("estimated-phi and size-adaptive alternation never decreases the objective", {
  d <- dyn_fixture(n_units = 50, n_time = 5)
  for (opt in list(list(phi_update = "pl"), list(penalty_size = "adaptive"),
                   list(phi_update = "pl", penalty_size = "adaptive"))) {
    ctrl <- do.call(scmr_control, c(list(lambda_scale = "sum", type_multinomial = "ungrouped",
                                         min_units = 6, min_per_class = 2, max_iter = 15,
                                         tiny_movement_max_units = 0, tiny_movement_rate_tol = 0,
                                         tiny_movement_revert = FALSE), opt))
    f <- fit_scmr(d$xx, d$y, G = 2, unit_id = d$unit_id, coords = d$coords, alpha = 0.5,
                  lambda = 0.02, control = ctrl, seed = 9)
    tr <- f$diagnostics$iterations$PenalizedObjective
    expect_true(all(diff(tr) >= -1e-6 * max(1, abs(tr))))
    expect_gte(f$criteria$PenalizedObjective, max(tr) - 1e-6 * max(1, abs(tr)))
    expect_true(is.finite(f$criteria$CriterionPLIC_AIC))
    if (identical(opt$phi_update, "pl")) {
      expect_true(f$criteria$phi_estimated)
      expect_equal(f$criteria$df_label, 1)
      expect_equal(f$phi, tail(f$diagnostics$iterations$Phi, 1))
    }
  }
})

brute_marginals <- function(trans, a0, y_obs) {
  # trans: T x C x C array (c, k); a0 initial distribution of y_0.
  Tn <- dim(trans)[1]
  C <- dim(trans)[2]
  paths <- as.matrix(expand.grid(rep(list(seq_len(C)), Tn + 1L)))
  prob <- a0[paths[, 1]]
  for (t in seq_len(Tn)) prob <- prob * trans[cbind(t, paths[, t + 1L], paths[, t])]
  keep <- rep(TRUE, nrow(paths))
  for (t in seq_len(Tn)) if (!is.na(y_obs[t])) keep <- keep & paths[, t + 1L] == y_obs[t]
  prob <- prob * keep
  t(vapply(seq_len(Tn), function(t) vapply(seq_len(C), function(c) sum(prob[paths[, t + 1L] == c]), numeric(1)),
           numeric(C))) / sum(prob)
}

test_that("forward filter and wave imputation equal brute-force enumeration", {
  d <- dyn_fixture(C = 3, n_units = 40, n_time = 4, pattern = "global")
  fit <- fit_scmr(d$xx, d$y, model = "global", alpha = 0.5, lambda = 0.01,
                  control = scmr_control(type_multinomial = "ungrouped"))
  rows <- which(d$unit_id == d$unit_id[1])
  p_f <- scmr_filter_predict(fit, d$xx[rows, ], d$unit_id[rows], d$time[rows], lag_columns = d$lag_cols)
  trans <- scmr:::transition_array(fit$fits[[1]], d$xx[rows, ], d$lag_cols, fit$class_levels)
  a0 <- as.numeric(table(fit$y_train)) / length(fit$y_train)
  expect_equal(unname(p_f), brute_marginals(trans, a0, rep(NA, 4)), tolerance = 1e-10)
  yo <- as.character(d$y[rows]); yo[c(2, 3)] <- NA
  p_i <- scmr_impute_waves(fit, d$xx[rows, ], yo, d$unit_id[rows], d$time[rows], lag_columns = d$lag_cols)
  expect_equal(unname(p_i[1:4, ]), brute_marginals(trans, a0, match(yo, fit$class_levels)),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("dynamic spatial fit predicts valid mixtures at new units", {
  d <- dyn_fixture(C = 3, n_units = 60, n_time = 5, pattern = "blocks", G = 2)
  train <- d$unit_id %in% sprintf("U%04d", 1:48)
  ctrl <- scmr_control(lambda_scale = "sum", type_multinomial = "ungrouped", min_units = 6,
                       min_per_class = 2, phi_update = "pl")
  fit <- fit_scmr(d$xx[train, ], d$y[train], G = 2, unit_id = d$unit_id[train], coords = d$coords[train, ],
                  alpha = 0.5, lambda = 0.02, control = ctrl, seed = 2)
  new <- !train
  for (rule in c("potts", "proportion", "majority")) {
    p <- predict(fit, d$xx[new, ], new_unit_id = d$unit_id[new], new_coords = d$coords[new, ], membership = rule)
    expect_equal(rowSums(p), rep(1, sum(new)), tolerance = 1e-10)
  }
  pf <- scmr_filter_predict(fit, d$xx[new, ], d$unit_id[new], d$time[new], d$coords[new, ], d$lag_cols)
  expect_equal(rowSums(pf), rep(1, sum(new)), tolerance = 1e-10)
  yo <- as.character(d$y[new]); yo[seq(2, length(yo), by = 2)] <- NA
  pi <- scmr_impute_waves(fit, d$xx[new, ], yo, d$unit_id[new], d$time[new], d$coords[new, ], d$lag_cols)
  expect_equal(rowSums(attr(pi, "cluster_posterior")), rep(1, 12), tolerance = 1e-10, ignore_attr = TRUE)
  obs <- which(!is.na(yo))
  expect_equal(pi[cbind(obs, match(yo[obs], fit$class_levels))], rep(1, length(obs)), tolerance = 1e-10)
})

test_that("panel generator, block folds and GW comparator behave", {
  sim <- simulate_scmr_panel(n_units = 90, n_time = 3, pattern = "irregular", G = 3, seed = 5)
  expect_setequal(unique(sim$unit_regime), 1:3)
  expect_equal(length(sim$y), 270L)
  expect_equal(rowSums(sim$lag_x), rep(1, 270))
  blocks <- substr(sim$unit_id, 1, 4)
  folds <- scmr_block_folds(sim$unit_id, nfolds = 3, seed = 1)
  expect_true(all(tapply(folds, sim$unit_id, function(z) length(unique(z))) == 1L))
  expect_true(max(table(tapply(folds, sim$unit_id, `[`, 1))) - min(table(tapply(folds, sim$unit_id, `[`, 1))) <= 1)
  gw <- fit_gw_multinom_en(cbind(sim$x, sim$lag_x), sim$y, sim$unit_id, sim$coords,
                           k_grid = c(20, 40), anchors = 15, seed = 2)
  p <- predict(gw, cbind(sim$x, sim$lag_x)[1:20, ], sim$coords[1:20, ])
  expect_equal(rowSums(p), rep(1, 20), tolerance = 1e-10, ignore_attr = TRUE)
  expect_true(gw$k %in% c(20, 40))
})

test_that("original SCR port runs, predicts and reports its BIC", {
  d <- dyn_fixture(C = 3, n_units = 60, n_time = 4, pattern = "blocks", G = 2)
  f <- fit_scr_original(d$xx, d$y, d$unit_id, d$coords, G = 2, kmeans_starts = 5, maxitr = 20)
  expect_true(isTRUE(f$scr_original))
  expect_lte(f$G, 2L)
  expect_true(is.finite(f$criteria$CriterionSCR_BIC_original))
  expect_identical(f$penalty, "none")
  p <- scmr_filter_predict(f, d$xx, d$unit_id, d$time, d$coords, d$lag_cols)
  expect_equal(rowSums(p), rep(1, nrow(d$xx)), tolerance = 1e-10)
})

test_that("panel generator supports the SCR domain, grid regions and GP smooth fields", {
  g <- simulate_scmr_panel(n_units = 120, n_time = 2, pattern = "grid", domain = "scr", seed = 4)
  expect_equal(length(g$beta), 6L)
  expect_true(all(g$unit_coords[, 1]^2 + 0.5 * g$unit_coords[, 2]^2 > 0.25))
  s <- simulate_scmr_panel(n_units = 80, n_time = 2, pattern = "smooth", smooth_type = "gp", seed = 4)
  expect_equal(dim(s$unit_beta), c(80L, 5L, 6L))
  expect_equal(apply(s$unit_beta, c(1, 2), sum), matrix(0, 80, 5), tolerance = 1e-10)
  expect_gt(stats::sd(s$unit_beta[, 1, 1]), 0)
})
