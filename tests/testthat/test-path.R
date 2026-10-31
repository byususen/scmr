sim_block_regimes <- function(side = 10, n_b = 10, seed = 1) {
  set.seed(seed)
  gx <- rep(seq_len(side), side)
  gy <- rep(seq_len(side), each = side)
  reg <- ifelse(gx <= 6, 1L, 2L)
  reg[gx >= 2 & gx <= 4 & gy >= 2 & gy <= 4] <- 3L          # 9-block regime
  B <- side^2
  beta <- list(cbind(c(0, 0, 0), rbind(c(1.5, -1, 0, 0), c(-1, 1.5, 0, 0), 0)),
               cbind(c(0, .5, -.5), rbind(c(-1.5, 1, 0, 0), c(1, -1.5, 0, 0), 0)),
               cbind(c(1, 0, 0), rbind(c(0, 0, 2, -2), c(0, 0, -2, 2), 0)))
  x <- matrix(stats::rnorm(B * n_b * 4), ncol = 4, dimnames = list(NULL, paste0("x", 1:4)))
  r <- rep(reg, each = n_b)
  eta <- t(vapply(seq_len(nrow(x)), function(i) drop(beta[[r[i]]] %*% c(1, x[i, ])), numeric(3)))
  pr <- exp(eta - apply(eta, 1, max))
  pr <- pr / rowSums(pr)
  y <- factor(apply(pr, 1, function(q) sample(c("a", "b", "c"), 1, prob = q)), levels = c("a", "b", "c"))
  list(x = x, y = y, id = rep(sprintf("B%03d", seq_len(B)), each = n_b),
       coords = cbind(rep(gx, each = n_b), rep(gy, each = n_b)), reg = reg)
}

path_control <- function(...) {
  scmr_control(lambda_scale = "sum", phi_update = "pl", min_units = 5, min_per_class = 2,
               k_neighbors = 8, max_iter = 20, tiny_movement_max_units = 0,
               tiny_movement_rate_tol = 0, tiny_movement_revert = FALSE, ...)
}

test_that("size-aware PLIC-BIC matches its regime-wise definition", {
  d <- sim_block_regimes()
  f <- fit_scmr(d$x, d$y, G = 3, unit_id = d$id, coords = d$coords, alpha = 0.5, lambda = 0.01,
                control = path_control(), seed = 3)
  cr <- f$criteria
  pc <- cr$cluster_components
  base <- -2 * cr$logLik - 2 * cr$PottsLogPseudoLik
  expect_equal(cr$CriterionPLIC_BIC_rows,
               base + sum(pc$DFEffective_g * log(pc$NRowsCluster)) + cr$df_label * log(cr$n_units))
  expect_equal(cr$CriterionPLIC_BIC_units,
               base + sum(pc$DFEffective_g * log(pc$NUnitsCluster)) + cr$df_label * log(cr$n_units))
  # Every regime is smaller than the sample, so the size-aware penalty is lighter.
  expect_lt(cr$CriterionPLIC_BIC_rows, cr$CriterionPLIC_BIC)
})

test_that("phi warm-up holds phi, then estimates it, with a monotone objective", {
  expect_error(scmr_control(phi_warmup_sweeps = -1), "nonnegative")
  expect_error(scmr_control(phi_warmup_sweeps = 1.5), "nonnegative")
  d <- sim_block_regimes()
  f <- fit_scmr(d$x, d$y, G = 3, unit_id = d$id, coords = d$coords, alpha = 0.5, lambda = 0.01,
                control = path_control(phi_warmup_sweeps = 2L, phi_warmup_value = 0), seed = 3)
  tr <- f$obj_trace
  expect_gt(nrow(tr), 2L)
  expect_equal(tr$Phi[1], 0)
  expect_gt(f$phi, 0)
  expect_true(all(diff(tr$PenalizedObjective) >= -1e-6))
  # The alternation does not stop at the sweep that releases phi.
  expect_gt(nrow(tr), min(which(tr$Phi > 0)))
})

test_that("split-merge path recovers a small regime and records every candidate", {
  d <- sim_block_regimes()
  p <- scmr_split_merge_path(d$x, d$y, d$id, d$coords, G_max = 4, alpha = 0.5, lambda = 0.01,
                             control = path_control(), seed = 7, n_split = 2L, n_merge = 1L)
  expect_s3_class(p, "scmr_path")
  expect_equal(p$table$G, 1:4)
  expect_true(all(c("split", "merge") %in% p$log$Step))
  expect_true(all(is.finite(p$table$Objective)))
  tab <- table(d$reg, p$fits[[3]]$group_unit[p$units])
  # Each true regime maps to one estimated regime (up to a few boundary units).
  expect_true(all(apply(tab, 1, max) / rowSums(tab) >= 0.9))
  expect_equal(unname(p$selected["PLIC_BIC_rows"]), 3)
  # keep_fits = FALSE keeps memberships and criteria only.
  q <- scmr_split_merge_path(d$x, d$y, d$id, d$coords, G_max = 2, alpha = 0.5, lambda = 0.01,
                             control = path_control(), seed = 7, merge_back = FALSE, keep_fits = FALSE)
  expect_null(q$fits[[2]]$fits)
  expect_length(q$fits[[2]]$group_unit, 100L)
})
