ext_fixture <- function() {
  set.seed(23)
  n_unit <- 60
  coords_u <- cbind(runif(n_unit), runif(n_unit))
  regime <- ifelse(coords_u[, 1] < 0.5, 1L, 2L)
  id <- rep(sprintf("U%02d", seq_len(n_unit)), each = 3)
  coords <- coords_u[rep(seq_len(n_unit), each = 3), ]
  x <- matrix(rnorm(3 * n_unit * 3), ncol = 3, dimnames = list(NULL, c("a", "b", "c")))
  b <- rbind(c(1.5, -1, 0), c(-1.5, 1, 0))[regime[rep(seq_len(n_unit), each = 3)], ]
  eta <- cbind(x[, 1] * b[, 1], x[, 2] * b[, 2], 0)
  pr <- exp(eta) / rowSums(exp(eta))
  y <- factor(apply(pr, 1, function(p) sample(c("k1", "k2", "k3"), 1, prob = p)))
  list(x = x, y = y, id = id, coords = coords)
}

test_that("Potts pseudo-likelihood matches a hand calculation", {
  w <- Matrix::Matrix(matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), 3), sparse = TRUE)
  expected <- -log(2) + 1 - 2 * log(exp(1) + 1)
  expect_equal(scmr:::potts_log_pseudolikelihood(w, c(1, 1, 2), phi = 1), expected)
  expect_equal(scmr:::potts_log_pseudolikelihood(w, c(1, 1, 1), phi = 1, G = 1), 0)
})

test_that("sum-scale alternation never decreases the penalized objective", {
  d <- ext_fixture()
  ctrl <- scmr_control(lambda_scale = "sum", type_multinomial = "ungrouped", min_units = 8,
                       min_per_class = 3, max_iter = 15, tiny_movement_max_units = 0,
                       tiny_movement_rate_tol = 0, tiny_movement_revert = FALSE)
  f <- fit_scmr(d$x, d$y, G = 2, unit_id = d$id, coords = d$coords, penalty = "elastic_net",
                alpha = 0.5, lambda = 0.02, control = ctrl, seed = 5)
  tr <- f$diagnostics$iterations$PenalizedObjective
  expect_true(all(diff(tr) >= -1e-6 * max(1, abs(tr))))
  expect_gte(f$criteria$PenalizedObjective, max(tr) - 1e-6 * max(1, abs(tr)))
  expect_identical(f$lambda_scale, "sum")
  expect_true(all(is.finite(c(f$criteria$CriterionPLIC_BIC, f$criteria$CriterionPLIC_AIC))))
  expect_lte(f$criteria$PottsLogPseudoLik, 0)
})

test_that("globally standardized engines back-transform coefficients exactly", {
  d <- ext_fixture()
  labels <- ifelse(d$coords[, 1] < 0.5, "west", "east")
  ctrl <- scmr_control(lambda_scale = "sum", type_multinomial = "ungrouped", min_units = 8, min_per_class = 3)
  f <- fit_scmr(d$x * 7 + 3, d$y, model = "fixed_clusters", cluster = labels, unit_id = d$id,
                penalty = "elastic_net", alpha = 0.5, lambda = 0.02, control = ctrl)
  b <- coef(f)
  for (g in seq_len(f$G)) {
    idx <- which(f$row_group == g)
    eta <- cbind(1, (d$x * 7 + 3)[idx, ]) %*% t(b[g, , ])
    pr <- exp(eta - apply(eta, 1, max)); pr <- pr / rowSums(pr)
    expect_equal(unname(f$prob_train[idx, ]), unname(pr), tolerance = 1e-8)
  }
})

test_that("two-stage and multi-start variants are recorded", {
  d <- ext_fixture()
  base <- list(type_multinomial = "ungrouped", min_units = 8, min_per_class = 3,
               coef_init_k = 10, coef_init_anchors = 30)
  two <- fit_scmr(d$x, d$y, G = 2, unit_id = d$id, coords = d$coords, alpha = 0.5, lambda = 0.02,
                  control = do.call(scmr_control, c(base, update_memberships = FALSE)), seed = 3)
  expect_true(two$converged)
  expect_identical(two$convergence_reason, "two_stage_initial_partition")
  ms <- fit_scmr(d$x, d$y, G = 2, unit_id = d$id, coords = d$coords, alpha = 0.5, lambda = 0.02,
                 control = do.call(scmr_control, c(base, lambda_scale = "sum", n_starts = 2L)), seed = 3)
  st <- ms$diagnostics$starts
  expect_equal(nrow(st), 2L)
  expect_equal(sum(st$Selected), 1L)
  expect_equal(st$Objective[st$Selected], max(st$Objective, na.rm = TRUE))
  feats <- scmr_local_coefficients(d$x, d$y, d$id, d$coords, k = 10, anchors = 20)
  expect_equal(dim(feats), c(60L, 3L * 4L))
})

test_that("irregular scenario produces disconnected Latin-square regimes", {
  dat <- simulate_scmr_data(n_obs = 300, n_new = 30, p = 4, active = 1:2,
                            scenario = "clustered_irregular", seed = 11)
  expect_setequal(unique(dat$observed$true_cluster), c("Q1", "Q2", "Q3"))
  ctrl <- scmr_simulation_control()
  centers <- as.matrix(expand.grid(s1 = c(-2, 0, 2) / 3, s2 = c(1, 3, 5) / 3))
  reg <- scmr:::irregular_regime_index(centers, ctrl)
  expect_true(all(table(reg) == 3L))
  expect_equal(apply(dat$observed$beta, c(1, 3), sum), matrix(0, 300, 5), tolerance = 1e-12,
               ignore_attr = TRUE)
})
