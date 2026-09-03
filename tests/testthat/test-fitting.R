fit_fixture <- function() {
  set.seed(17)
  list(x = matrix(rnorm(360), 120, 3, dimnames = list(NULL, c("a", "b", "c"))),
    y = factor(rep(c("first", "second", "third"), 40), levels = c("third", "first", "second")),
    id = rep(sprintf("U%02d", 1:40), each = 3),
    coords = cbind(rep(seq_len(40), each = 3), rep(rep(1:2, 20), each = 3)))
}

test_that("all penalties use exactly the same G1 path as global fitting", {
  d <- fit_fixture()
  ctrl <- scmr_control(type_multinomial = "ungrouped", min_units = 3, min_per_class = 2,
                       alpha_grid = c(.3, .7))
  for (penalty in c("none", "ridge", "lasso", "elastic_net")) {
    args <- list(x = d$x, y = d$y, unit_id = d$id, coords = d$coords,
                 penalty = penalty, control = ctrl, seed = 99, nfolds = 3)
    global <- do.call(fit_scmr, c(args, list(model = "global")))
    spatial <- do.call(fit_scmr, c(args, list(model = "scmr", G = 1)))
    expect_equal(spatial$prob_train, global$prob_train, tolerance = 1e-12)
    expect_equal(coef(spatial), coef(global), tolerance = 1e-12)
    expect_equal(spatial$lambda_used, global$lambda_used)
    expect_equal(spatial$criteria, global$criteria)
    expect_true(global$converged)
    expect_equal(dim(predict(global, d$x[1, ])), c(1L, 3L))
    expect_equal(dim(predict(spatial, d$x[1:3, ])), c(3L, 3L))
    if (penalty != "none") {
      expect_equal(rowSums(global$prob_oof_train), rep(1, 120), tolerance = 1e-10)
      # glmnet's prevalidated predictions are logits, not probabilities.
      j <- which.min(abs(global$fit$lambda - global$lambda))
      link <- global$fit$fit.preval[, , j]
      expected <- exp(link - apply(link, 1, max))
      expected <- expected / rowSums(expected)
      expect_equal(unname(global$prob_oof_train), unname(expected), tolerance = 1e-10)
    }
  }
})

test_that("unpenalized coefficients reproduce predictions on the input scale", {
  d <- fit_fixture()
  for (k in c(2, 3)) {
    keep <- d$y %in% levels(d$y)[seq_len(k)]
    x <- d$x[keep, 1, drop = FALSE] * 20 + 100
    f <- fit_scmr(x, droplevels(d$y[keep]), model = "global", penalty = "none")
    b <- matrix(coef(f)[1, , ], nrow = k)
    eta <- cbind(1, x) %*% t(b)
    pr <- exp(eta - apply(eta, 1, max)); pr <- pr / rowSums(pr)
    expect_equal(unname(predict(f, x)), unname(pr), tolerance = 1e-8)
    expect_equal(colSums(b), rep(0, 2), tolerance = 1e-12)
    expect_equal(as.numeric(attr(logLik(f), "df")), (k - 1) * 2)
  }
})

test_that("fixed clusters preserve labels and spatial moves preserve unit feasibility", {
  d <- fit_fixture()
  labels <- ifelse(d$coords[, 1] <= 20, "west", "east")
  ctrl <- scmr_control(min_units = 15, min_per_class = 2, max_iter = 2,
    type_multinomial = "ungrouped", tiny_movement_max_units = 0, tiny_movement_rate_tol = 0)
  for (penalty in c("none", "ridge", "lasso", "elastic_net")) {
    extra <- if (penalty == "none") list() else list(lambda = c(.03, .05))
    f <- do.call(fit_scmr, c(list(x = d$x, y = d$y, unit_id = d$id,
      model = "fixed_clusters", cluster = labels, penalty = penalty, control = ctrl), extra))
    expect_identical(f$fixed_cluster_levels[f$row_group], labels)
    expect_error(predict(f, d$x[1, , drop = FALSE], cluster = "unknown"), "unknown")
    expect_equal(rowSums(predict(f, d$x, cluster = labels)), rep(1, 120), tolerance = 1e-10)
    expect_equal(f$criteria$CBBIC_WeightedN_Check, 120)
  }
  f <- fit_scmr(d$x, d$y, G = 2, unit_id = d$id, coords = d$coords,
    penalty = "elastic_net", lambda = .03, control = ctrl, seed = 43)
  expect_true(all(table(f$group_unit) >= 15))
  expect_true(all(vapply(split(f$row_group, d$id), function(z) length(unique(z)) == 1L, logical(1))))
  expect_true(all(table(f$row_group, d$y) >= 2))
  expect_equal(predict(f, d$x, d$id, d$coords), f$prob_train, tolerance = 1e-10)
  expect_equal(rowSums(predict(f, d$x[1:3, ], c("new1", "new2", "new3"), d$coords[1:3, ])), rep(1, 3))
  expect_false(any(f$diagnostics$tuning$Stage == "final"))
})

test_that("EDF and likelihood criteria obey independently computed identities", {
  d <- fit_fixture()
  f <- fit_scmr(d$x, d$y, model = "global", penalty = "none", unit_id = d$id)
  ll <- sum(log(f$prob_train[cbind(seq_along(d$y), as.integer(d$y))]))
  expect_equal(as.numeric(logLik(f)), ll)
  expect_equal(f$criteria$CriterionSCR_BIC_original, -2 * ll + log(120) * 8)
  expect_equal(f$criteria$CriterionCB_BIC, -2 * ll + log(40) * 8)
  expect_equal(f$criteria$CriterionCB_AIC, -2 * ll + 16)
  for (type in c("grouped", "ungrouped")) {
    df <- vapply(c(.001, 10), function(lam) {
      z <- fit_scmr(d$x, d$y, model = "global", penalty = "ridge", lambda = lam,
                    control = scmr_control(type_multinomial = type))
      z$criteria$df_effective_total
    }, numeric(1))
    expect_true(all(df >= 2 & df <= 8 + 1e-7))
    expect_lt(df[2], df[1])
  }
  w <- Matrix::Matrix(matrix(c(0, 2, 1, 2, 0, 3, 1, 3, 0), 3), sparse = TRUE)
  expect_equal(scmr:::spatial_bonus(w, c(1, 1, 2)), 2)
})

test_that("seeding and validation protect caller state", {
  d <- fit_fixture(); before <- .Random.seed
  f <- fit_scmr(d$x, d$y, model = "global", penalty = "lasso", nfolds = 3)
  expect_identical(.Random.seed, before)
  expect_error(fit_scmr(d$x, d$y, model = "global", lambda = -1), "lambda")
  expect_error(fit_scmr(d$x, factor(d$y, levels = c(levels(d$y), "absent")), model = "global"), "every declared")
  expect_error(fit_scmr(d$x, d$y, G = 2, unit_id = d$id, coords = matrix(0, 2, 2)), "coords")
  expect_error(fit_scmr(d$x, d$y, model = "fixed_clusters", unit_id = d$id,
                       cluster = rep(1:2, 60)), "multiple fixed")
})
