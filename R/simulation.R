inside_irregular_domain <- function(coords) {
  coords[, 1]^2 + 0.5 * coords[, 2]^2 > 0.5^2
}

sample_domain <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- matrix(NA_real_, nrow = n, ncol = 2)
  colnames(out) <- c("s1", "s2")
  got <- 0L
  while (got < n) {
    m <- max(1000L, 2L * (n - got))
    cand <- cbind(stats::runif(m, -1, 1), stats::runif(m, 0, 2))
    cand <- cand[inside_irregular_domain(cand), , drop = FALSE]
    take <- min(n - got, nrow(cand))
    if (take > 0) {
      out[(got + 1):(got + take), ] <- cand[seq_len(take), , drop = FALSE]
      got <- got + take
    }
  }
  out
}

region_levels6 <- c("R11", "R12", "R13", "R21", "R22", "R23")

region_label_6 <- function(coords) {
  s1 <- coords[, 1]
  s2 <- coords[, 2]
  j <- ifelse(s1 <= 0, 1L, 2L)
  k <- cut(s2, breaks = c(0, 2/3, 4/3, 2), include.lowest = TRUE, labels = FALSE)
  paste0("R", j, k)
}

region_index_6 <- function(coords) {
  s1 <- coords[, 1]
  s2 <- coords[, 2]
  j <- ifelse(s1 <= 0, 0L, 1L)
  k <- cut(s2, breaks = c(0, 2/3, 4/3, 2), include.lowest = TRUE, labels = FALSE) - 1L
  list(j = j, k = k)
}

exp_cov_matrix <- function(coords, range_param, nugget = 1e-8, kernel_power = 1) {
  D <- as.matrix(stats::dist(coords))
  Sigma <- exp(- (D / range_param)^kernel_power)
  diag(Sigma) <- diag(Sigma) + nugget
  Sigma
}

sim_gp <- function(coords, range_param, tau2 = 1, seed = NULL, kernel_power = 1) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- tau2 * exp_cov_matrix(coords, range_param, kernel_power = kernel_power)
  as.numeric(MASS::mvrnorm(n = 1, mu = rep(0, nrow(coords)), Sigma = Sigma))
}

simulate_predictors <- function(coords, p = 25, active = 1:5,
                                rho_active = 0.50, rho_inactive = 0.10,
                                eta = 0.2, predictor_eta_scale = 0.30,
                                kernel_power = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(coords)
  X <- matrix(NA_real_, nrow = n, ncol = p)
  eta_eff <- eta * predictor_eta_scale
  common_active <- sim_gp(coords, eta_eff, tau2 = 1, kernel_power = kernel_power)
  common_inactive <- sim_gp(coords, eta_eff, tau2 = 1, kernel_power = kernel_power)
  active_set <- sort(unique(active))
  inactive_set <- setdiff(seq_len(p), active_set)
  for (m in active_set) {
    zm <- sim_gp(coords, eta_eff, tau2 = 1, kernel_power = kernel_power)
    X[, m] <- rho_active * common_active + sqrt(1 - rho_active^2) * zm
  }
  for (m in inactive_set) {
    zm <- sim_gp(coords, eta_eff, tau2 = 1, kernel_power = kernel_power)
    X[, m] <- rho_inactive * common_inactive + sqrt(1 - rho_inactive^2) * zm
  }
  X <- as.matrix(scale(X))
  X[!is.finite(X)] <- 0
  colnames(X) <- paste0("x", seq_len(p))
  X
}

simulate_beta_global <- function(coords, p = 25, active = 1:5, class_levels = c("C1", "C2", "C3")) {
  n <- nrow(coords)
  C <- length(class_levels)
  beta <- array(0, dim = c(n, C, p + 1), dimnames = list(NULL, class_levels, c("(Intercept)", paste0("x", seq_len(p)))))
  a_base <- c(0.50, -0.45, 0.40, -0.35, 0.30)
  b_base <- c(-0.45, 0.40, -0.35, 0.30, -0.25)
  beta[, 1, 1] <- 0.30
  beta[, 2, 1] <- -0.30
  beta[, 3, 1] <- -(beta[, 1, 1] + beta[, 2, 1])
  for (jj in seq_along(active)) {
    m <- active[jj]
    beta[, 1, m + 1] <- a_base[jj]
    beta[, 2, m + 1] <- b_base[jj]
    beta[, 3, m + 1] <- -(beta[, 1, m + 1] + beta[, 2, m + 1])
  }
  beta
}

simulate_beta_clustered <- function(coords, p = 25, active = 1:5, class_levels = c("C1", "C2", "C3")) {
  n <- nrow(coords)
  C <- length(class_levels)
  beta <- array(0, dim = c(n, C, p + 1), dimnames = list(NULL, class_levels, c("(Intercept)", paste0("x", seq_len(p)))))
  idx <- region_index_6(coords)
  g1_vals <- c(-1, 0)
  g2_vals <- c(0, 2/3, 4/3)
  a_base <- c(0.9, -0.8, 0.7, -0.6, 0.5)
  b_base <- c(-0.8, 0.7, -0.6, 0.5, -0.4)
  for (i in seq_len(n)) {
    g1 <- g1_vals[idx$j[i] + 1]
    g2 <- g2_vals[idx$k[i] + 1]
    q1 <- g1 + g2 - 0.25
    q2 <- g1 - g2 + 0.25
    beta[i, 1, 1] <- 0.8 * q1
    beta[i, 2, 1] <- -0.7 * q2
    beta[i, 3, 1] <- -(beta[i, 1, 1] + beta[i, 2, 1])
    for (jj in seq_along(active)) {
      m <- active[jj]
      beta[i, 1, m + 1] <- a_base[jj] * q1
      beta[i, 2, m + 1] <- b_base[jj] * q2
      beta[i, 3, m + 1] <- -(beta[i, 1, m + 1] + beta[i, 2, m + 1])
    }
  }
  beta
}

softmax_from_beta <- function(X, beta) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  C <- dim(beta)[2]
  eta <- matrix(0, nrow = n, ncol = C)
  for (cc in seq_len(C)) {
    beta_c <- matrix(beta[, cc, -1, drop = TRUE], nrow = n, ncol = p)
    eta[, cc] <- beta[, cc, 1] + rowSums(X * beta_c)
  }
  P <- softmax_rows(eta)
  colnames(P) <- dimnames(beta)[[2]]
  P
}

sample_multiclass <- function(P, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  classes <- colnames(P)
  apply(P, 1, function(pr) sample(classes, size = 1, prob = pr))
}

#' Simulate data for SCMR experiments
#'
#' @param n_obs Number of observed samples.
#' @param n_new Number of independent new samples.
#' @param p Number of predictors.
#' @param active Active predictor indices.
#' @param eta Spatial range scale for predictors.
#' @param scenario One of `"global"`, `"clustered_balanced"`, or `"clustered_imbalanced"`.
#' @param class_levels Class labels.
#' @param seed Random seed.
#' @return A list with `observed` and `new` data lists.
#' @export
simulate_scmr_data <- function(n_obs = 3000, n_new = 300, p = 25, active = 1:5,
                               eta = 0.2,
                               scenario = c("global", "clustered_balanced", "clustered_imbalanced"),
                               class_levels = c("C1", "C2", "C3"), seed = 1000) {
  scenario <- match.arg(scenario)
  coords_all <- sample_domain(n_obs + n_new, seed = seed + 1)
  region_all <- region_label_6(coords_all)
  X_all <- simulate_predictors(coords_all, p = p, active = active, eta = eta, seed = seed + 2)
  beta_all <- switch(
    scenario,
    global = simulate_beta_global(coords_all, p = p, active = active, class_levels = class_levels),
    clustered_balanced = simulate_beta_clustered(coords_all, p = p, active = active, class_levels = class_levels),
    clustered_imbalanced = simulate_beta_clustered(coords_all, p = p, active = active, class_levels = class_levels)
  )
  P_all <- softmax_from_beta(X_all, beta_all)
  Y_all <- factor(sample_multiclass(P_all, seed = seed + 4), levels = class_levels)
  obs_idx <- seq_len(n_obs)
  new_idx <- (n_obs + 1):(n_obs + n_new)
  list(
    observed = list(
      X = X_all[obs_idx, , drop = FALSE],
      y = Y_all[obs_idx],
      coords = coords_all[obs_idx, , drop = FALSE],
      beta = beta_all[obs_idx, , , drop = FALSE],
      prob = P_all[obs_idx, , drop = FALSE],
      region = region_all[obs_idx],
      unit = paste0("U", obs_idx)
    ),
    new = list(
      X = X_all[new_idx, , drop = FALSE],
      y = Y_all[new_idx],
      coords = coords_all[new_idx, , drop = FALSE],
      beta = beta_all[new_idx, , , drop = FALSE],
      prob = P_all[new_idx, , drop = FALSE],
      region = region_all[new_idx],
      unit = paste0("N", seq_along(new_idx))
    ),
    scenario = scenario,
    active = active,
    class_levels = class_levels
  )
}
