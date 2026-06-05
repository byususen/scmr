test_that("simulation returns expected objects", {
  dat <- simulate_scmr_data(n_obs = 50, n_new = 10, p = 6, active = 1:3, seed = 1)
  expect_equal(nrow(dat$observed$X), 50)
  expect_equal(nrow(dat$new$X), 10)
  expect_true(all(levels(dat$observed$y) == c("C1", "C2", "C3")))
})

test_that("metrics are finite for simple probabilities", {
  y <- factor(c("a", "b", "a"), levels = c("a", "b"))
  prob <- matrix(c(.8, .2, .3, .7, .6, .4), ncol = 2, byrow = TRUE)
  colnames(prob) <- c("a", "b")
  pred <- factor(colnames(prob)[max.col(prob)], levels = levels(y))
  m <- calc_overall_metrics(y, pred, prob, levels(y))
  expect_true(is.finite(m$Accuracy))
  expect_true(is.finite(m$LogLoss))
})

test_that("penalty alpha resolver handles ridge/lasso/enet", {
  expect_equal(scmr:::resolve_penalty_alpha("ridge")$alpha_default, 0)
  expect_equal(scmr:::resolve_penalty_alpha("lasso")$alpha_default, 1)
  expect_equal(scmr:::resolve_penalty_alpha("elastic_net", alpha_grid = c(0.4, 0.8))$alpha_grid, c(0.4, 0.8))
})
