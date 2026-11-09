phase_sim <- function(side = 12, n_b = 40, seed = 1, smooth = FALSE, regimes = TRUE) {
  set.seed(seed)
  gx <- rep(seq_len(side), side); gy <- rep(seq_len(side), each = side)
  B <- side^2
  reg <- if (regimes) ifelse(gx <= side / 2, 1L, 2L) else rep(1L, B)
  x <- matrix(stats::rnorm(B * n_b * 3), ncol = 3, dimnames = list(NULL, paste0("x", 1:3)))
  r <- rep(reg, each = n_b)
  f <- rep((gx + gy - side - 1) / (side - 1), each = n_b)
  eta2 <- ifelse(r == 1L, 1.5, -1.5) * x[, 1] + ifelse(r == 1L, 0.8, -0.8)
  if (smooth) eta2 <- 2 * f * x[, 1] + 1.2 * f
  eta <- cbind(0, eta2, 0.5 * x[, 2])
  pr <- exp(eta - apply(eta, 1, max)); pr <- pr / rowSums(pr)
  y <- factor(apply(pr, 1, function(q) sample(c("a", "b", "c"), 1, prob = q)), levels = c("a", "b", "c"))
  list(x = x, y = y, id = rep(sprintf("U%03d", seq_len(B)), each = n_b),
       coords = cbind(rep(gx, each = n_b), rep(gy, each = n_b)), reg = reg)
}
phase_ctrl <- function() scmr_control(lambda_scale = "sum", phi_update = "pl", min_units = 5, min_per_class = 2,
                                      k_neighbors = 8, max_iter = 20, tiny_movement_max_units = 0,
                                      tiny_movement_rate_tol = 0, tiny_movement_revert = FALSE, n_starts = 1)

test_that("phase test flags discrete regimes but not their correct partition", {
  skip_if_not_installed("diptest")
  d <- phase_sim()
  f1 <- fit_scmr(d$x, d$y, G = 1, unit_id = d$id, coords = d$coords, alpha = 0.5, lambda = 0.005,
                 control = phase_ctrl(), seed = 1)
  t1 <- scmr_phase_test(f1, d$x, d$y, d$id, d$coords, B = 60, n_perm = 499, seed = 2)
  expect_equal(nrow(t1), 1L)
  expect_identical(t1$Phase, "discrete")
  expect_s3_class(attr(t1, "spikes"), "data.frame")
  units <- sort(unique(d$id))
  truth <- d$reg[match(units, sprintf("U%03d", seq_len(144)))]
  f2 <- fit_scmr(d$x, d$y, G = 2, unit_id = d$id, coords = d$coords, alpha = 0.5, lambda = 0.005,
                 control = phase_ctrl(), seed = 1, initial_cluster = truth)
  t2 <- scmr_phase_test(f2, d$x, d$y, d$id, d$coords, B = 60, n_perm = 499, seed = 2)
  expect_equal(nrow(t2), 2L)
  expect_false(any(t2$Phase == "discrete"))
})

test_that("phase test calls a smooth gradient continuous, not discrete", {
  skip_if_not_installed("diptest")
  d <- phase_sim(smooth = TRUE, regimes = FALSE, seed = 3)
  f1 <- fit_scmr(d$x, d$y, G = 1, unit_id = d$id, coords = d$coords, alpha = 0.5, lambda = 0.005,
                 control = phase_ctrl(), seed = 1)
  t1 <- scmr_phase_test(f1, d$x, d$y, d$id, d$coords, B = 60, n_perm = 499, seed = 2)
  expect_false(t1$Phase == "discrete")
  expect_true(t1$OrderedSpikes >= 1)
})
