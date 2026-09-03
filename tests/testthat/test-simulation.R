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

test_that("new data generators preserve region proportions and identification", {
  set.seed(21); rng <- .Random.seed
  for (scenario in c("global", "clustered_balanced", "clustered_imbalanced", "smooth")) {
    d <- simulate_scmr_data(90, 18, p = 4, active = 1:2, scenario = scenario, seed = 73)
    expect_equal(rowSums(d$observed$prob), rep(1, 90), tolerance = 1e-12)
    expect_equal(as.vector(apply(d$observed$beta, c(1, 3), sum)), rep(0, 90 * 5), tolerance = 1e-12)
    expect_identical(.Random.seed, rng)
    if (scenario == "clustered_balanced") expect_equal(as.numeric(d$observed_region_counts), rep(15, 6))
    if (scenario == "clustered_imbalanced") expect_equal(as.numeric(d$observed_region_counts), c(25, 20, 20, 15, 5, 5))
    if (scenario == "smooth") expect_true(all(is.na(d$observed$true_cluster)))
  }
  zero <- simulate_scmr_data(30, 0, p = 2, active = integer(), scenario = "clustered_balanced")
  expect_equal(dim(zero$new$beta), c(0L, 3L, 3L))
  expect_true(all(zero$observed$beta[, , -1] == 0))
  d <- simulate_scmr_data(36, 0, p = 2, active = 1, scenario = "clustered_balanced", heterogeneity_strength = 0)
  expect_true(all(d$observed$beta[, 1, 2] == .5))
  expect_error(scmr_simulation_control(domain_exclusion_radius = 100), "entire")
  expect_error(simulate_scmr_data(10, p = 2, active = 3), "active")
  expect_error(simulate_scmr_data(10, class_levels = c("a", "b")), "three")
})

test_that("stratified splitting retains training rows without changing RNG", {
  y <- factor(c("a", "a", "b", "b", "c"))
  set.seed(17); before <- .Random.seed
  test <- make_stratified_split(y, test_prop = .99)
  expect_equal(table(y[-test]), table(factor(c("a", "b", "c"))))
  expect_identical(.Random.seed, before)
  expect_error(make_stratified_split(y, region = "x"), "region")
})
