# Score-spectrum (phase-transition) diagnostics for SCMR regimes.
#
# For a fitted regime with units b = 1..m, let u_b be the unit score vector
# (sum over the unit's rows of (e_y - p) (x) (1, x), classes 1..C-1) and J-bar
# the average per-row Fisher information times the mean unit size. Under a
# homogeneous regime the whitened scores w_b = (nbar / n_b)^{1/2} Jbar^{-1/2} u_b
# are approximately mean-zero with covariance sigma^2 I_d (sigma^2 = 1 for a
# correct model; > 1 under within-unit dependence or overdispersion), so the
# eigenvalues of S = (1/m) sum_b w_b w_b' follow the Marchenko-Pastur law with
# ratio gamma = d / m and upper edge sigma^2 (1 + sqrt(gamma))^2. A hidden
# sub-regime adds a spike; by the BBP phase transition it separates from the
# bulk only when its strength exceeds sigma^2 sqrt(gamma).

mp_quantile <- function(prob, gamma, sigma2 = 1, grid = 4000L) {
  a <- sigma2 * (1 - sqrt(gamma))^2; b <- sigma2 * (1 + sqrt(gamma))^2
  x <- seq(a, b, length.out = grid + 2L)[-c(1L, grid + 2L)]
  dens <- sqrt(pmax((b - x) * (x - a), 0)) / (2 * pi * gamma * sigma2 * x)
  cdf <- cumsum(dens); cdf <- cdf / cdf[length(cdf)]
  stats::approx(cdf, x, xout = prob, ties = "ordered", rule = 2)$y
}

# Tracy-Widom (beta = 1) upper quantile of the largest eigenvalue (Johnstone 2001).
tw_crit <- function(m, d, sigma2, level = 0.95) {
  q <- c(`0.9` = 0.4501, `0.95` = 0.9793, `0.99` = 2.0234)[as.character(level)]
  mu <- (sqrt(m - 1) + sqrt(d))^2
  s <- (sqrt(m - 1) + sqrt(d)) * (1 / sqrt(m - 1) + 1 / sqrt(d))^(1 / 3)
  sigma2 * (mu + s * q) / m
}

# Unit scores and pooled Fisher information for rows `rows_g` of a fit.
# Optional reduction of the score dimension: with q_proj > 0 the standardized
# predictors are replaced by their first q_proj principal components (computed on
# all rows), so d = (C - 1) (q_proj + 1). A smaller d / m keeps the spike
# eigenvectors close to the signal direction (the eigenvector overlap in the
# spiked model falls as d / m grows), which the loading tests need.
unit_scores_info <- function(x, y, prob, unit_num, rows_g = seq_len(nrow(x)), q_proj = 0L) {
  Xs <- scale(x); Xs[!is.finite(Xs)] <- 0
  if (q_proj > 0L && q_proj < ncol(Xs)) {
    pcs <- .scmr_pc_cache(Xs, q_proj)
    Xs <- Xs %*% pcs
  }
  Z <- cbind(1, Xs)
  Z <- Z[rows_g, , drop = FALSE]
  P <- prob[rows_g, , drop = FALSE]
  C <- ncol(P); K <- C - 1L; p1 <- ncol(Z)
  Y <- matrix(0, length(rows_g), C); Y[cbind(seq_along(rows_g), as.integer(y[rows_g]))] <- 1
  u_rows <- do.call(cbind, lapply(seq_len(K), function(c) Z * (Y[, c] - P[, c])))
  un <- factor(unit_num[rows_g])
  U <- rowsum(u_rows, un, reorder = TRUE)
  n_b <- as.numeric(table(un))
  J <- matrix(0, K * p1, K * p1)
  for (c1 in seq_len(K)) for (c2 in c1:K) {
    w <- if (c1 == c2) P[, c1] * (1 - P[, c1]) else -P[, c1] * P[, c2]
    blk <- crossprod(Z * w, Z)
    J[(c1 - 1) * p1 + seq_len(p1), (c2 - 1) * p1 + seq_len(p1)] <- blk
    J[(c2 - 1) * p1 + seq_len(p1), (c1 - 1) * p1 + seq_len(p1)] <- t(blk)
  }
  list(U = U, n_b = n_b, J = J / length(rows_g))     # J per row
}

.scmr_pc_env <- new.env()
.scmr_pc_cache <- function(Xs, q) {
  key <- paste(nrow(Xs), ncol(Xs), q, sum(Xs[seq_len(min(50L, nrow(Xs))), 1]))
  if (!identical(.scmr_pc_env$key, key)) {
    sv <- if (nrow(Xs) > 20000L) svd(Xs[sort(sample.int(nrow(Xs), 20000L)), ], nu = 0, nv = q) else svd(Xs, nu = 0, nv = q)
    assign("rot", sv$v[, seq_len(q), drop = FALSE], envir = .scmr_pc_env)
    assign("key", key, envir = .scmr_pc_env)
  }
  .scmr_pc_env$rot
}

# Spectrum of whitened unit scores and the number of spikes above the edge.
#   ridge: relative ridge on the information before whitening
#   proj: if not NULL, number of leading information-weighted directions kept
score_spectrum <- function(si, ridge = 1e-3, level = 0.95, max_spikes = 10L) {
  m <- nrow(si$U); d <- ncol(si$U)
  nbar <- mean(si$n_b)
  Jb <- si$J * nbar
  ev <- eigen(Jb, symmetric = TRUE)
  keep <- ev$values > ridge * max(ev$values)
  Wh <- ev$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(ev$values[keep]), sum(keep))
  W <- (si$U * sqrt(nbar / si$n_b)) %*% Wh
  d <- ncol(W)
  lam <- if (d <= m) eigen(crossprod(W) / m, symmetric = TRUE, only.values = TRUE)$values else
    eigen(tcrossprod(W) / m, symmetric = TRUE, only.values = TRUE)$values
  gamma <- d / m
  # Noise level from the bulk, re-estimated after removing detected spikes.
  k <- 0L
  repeat {
    bulk <- lam[(k + 1L):length(lam)]
    if (gamma > 1) bulk <- bulk[seq_len(max(1L, min(m, d) - k))]
    sigma2 <- stats::median(bulk) / mp_quantile(0.5, gamma)
    crit <- tw_crit(m, d, sigma2, level)
    k_new <- sum(lam[seq_len(min(max_spikes, length(lam)))] > crit)
    if (k_new == k) break
    k <- k_new
    if (k >= max_spikes) break
  }
  list(lambda = lam, m = m, d = d, gamma = gamma, sigma2 = sigma2, crit = crit,
       edge = sigma2 * (1 + sqrt(gamma))^2, spikes = k, ratio1 = lam[1] / crit)
}

# Scale-free version calibrated by a parametric bootstrap from the fitted
# probabilities: R_k = lambda_k / median(lambda). A common inflation of the
# score variance (within-unit dependence, overdispersion) cancels in R; the
# bootstrap supplies the finite-sample shape (unequal unit sizes, rank of the
# unit information, d / m) that the asymptotic Marchenko-Pastur law misses.
whitened_eigen <- function(U, n_b, J, ridge = 1e-3) {
  m <- nrow(U); nbar <- mean(n_b)
  ev <- eigen(J * nbar, symmetric = TRUE)
  keep <- ev$values > ridge * max(ev$values)
  Wh <- ev$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(ev$values[keep]), sum(keep))
  W <- (U * sqrt(nbar / n_b)) %*% Wh
  lam <- if (ncol(W) <= m) eigen(crossprod(W) / m, symmetric = TRUE, only.values = TRUE)$values else
    eigen(tcrossprod(W) / m, symmetric = TRUE, only.values = TRUE)$values
  lam[seq_len(min(m - 1L, ncol(W)))]          # drop the centring zero
}
spectrum_boot <- function(x, y, prob, unit_num, rows_g = seq_len(nrow(x)), B = 200L, level = 0.95,
                          max_spikes = 10L, seed = 1L, q_proj = 0L) {
  si <- unit_scores_info(x, y, prob, unit_num, rows_g, q_proj)
  lam <- whitened_eigen(si$U, si$n_b, si$J)
  R <- lam / stats::median(lam)
  P <- prob[rows_g, , drop = FALSE]
  cum <- t(apply(P, 1, cumsum))
  set.seed(seed)
  lev <- levels(y)
  Rb <- replicate(B, {
    ys <- factor(lev[1 + rowSums(stats::runif(nrow(P)) > cum)], levels = lev)
    yy <- y; yy[rows_g] <- ys
    sb <- unit_scores_info(x, yy, prob, unit_num, rows_g, q_proj)
    lb <- whitened_eigen(sb$U, sb$n_b, si$J)
    lb[1] / stats::median(lb)
  })
  crit <- stats::quantile(Rb, level, names = FALSE)
  list(R = R, crit = crit, spikes = sum(R[seq_len(min(max_spikes, length(R)))] > crit),
       m = nrow(si$U), d = ncol(si$U), boot = Rb)
}

# Spatial order of the leading directions. A sub-regime is a spatially coherent
# set of units, so the unit loadings of its spike are spatially autocorrelated
# ("ordered phase"); unstructured unit effects (e.g. segment random effects)
# also inflate eigenvalues but give loadings without spatial order. Moran's I of
# the loadings of the top k eigenvectors, with a permutation p-value.
#   A: m x m symmetric adjacency of the regime's units (same order as U rows)
spike_moran <- function(U, n_b, J, A, k = 3L, n_perm = 499L, ridge = 1e-3, seed = 1L) {
  m <- nrow(U); nbar <- mean(n_b)
  ev <- eigen(J * nbar, symmetric = TRUE)
  keep <- ev$values > ridge * max(ev$values)
  Wh <- ev$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(ev$values[keep]), sum(keep))
  W <- (U * sqrt(nbar / n_b)) %*% Wh
  e <- eigen(tcrossprod(W) / m, symmetric = TRUE)
  L <- e$vectors[, seq_len(min(k, m - 1L)), drop = FALSE]
  A <- as.matrix(A); diag(A) <- 0; S0 <- sum(A)
  moran <- function(z) { z <- z - mean(z); (m / S0) * sum(z * (A %*% z)) / sum(z^2) }
  set.seed(seed)
  out <- t(apply(L, 2, function(z) {
    I <- moran(z)
    perm <- replicate(n_perm, moran(sample(z)))
    c(I = I, p = (1 + sum(perm >= I)) / (n_perm + 1))
  }))
  data.frame(k = seq_len(nrow(out)), MoranI = out[, "I"], p = out[, "p"])
}

# ---------------------------------------------------------------------------
# Phase-coexistence test (first-order versus continuous transition).
# A discrete sub-regime makes the unit loadings of its spike two-valued plus
# noise ("two coexisting phases"); a smooth spatial gradient gives spatially
# ordered but unimodal loadings. For each spike eigenvector we test
#   (a) spatial order: Moran's I of the loadings, permutation p-value;
#   (b) coexistence: Sarle's bimodality coefficient of the loadings, against a
#       null of a *continuous* spatially autocorrelated field (Gaussian SAR
#       field on the same graph with the same Moran's I), so spatial order
#       alone does not count as a discrete phase.
# The spike is a discrete sub-regime only if both reject (intersection-union:
# p = max(p_a, p_b)). Within a regime the spike p-values are Bonferroni
# adjusted; across the regimes of one G they are Holm adjusted.
bimodality_coef <- function(z) {
  n <- length(z); z <- z - mean(z); s2 <- mean(z^2)
  g <- mean(z^3) / s2^1.5; k <- mean(z^4) / s2^2 - 3
  (g^2 + 1) / (k + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}
moran_I <- function(Z, A) {          # Z: m x B matrix of (centred) columns
  Z <- sweep(Z, 2, colMeans(Z))
  (nrow(Z) / sum(A)) * colSums(Z * (A %*% Z)) / colSums(Z^2)
}
sar_fields <- function(A, rho, B, seed) {
  Wr <- A / pmax(rowSums(A), 1)
  M <- diag(nrow(A)) - rho * Wr
  set.seed(seed)
  solve(M, matrix(stats::rnorm(nrow(A) * B), nrow(A)))
}
phase_order_test <- function(U, n_b, J, A, n_spikes, k = 5L, n_perm = 9999L, n_null = 999L,
                             ridge = 1e-3, seed = 1L) {
  m <- nrow(U); nbar <- mean(n_b)
  A <- as.matrix(A); diag(A) <- 0
  ev <- eigen(J * nbar, symmetric = TRUE)
  keep <- ev$values > ridge * max(ev$values)
  Wh <- ev$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(ev$values[keep]), sum(keep))
  W <- (U * sqrt(nbar / n_b)) %*% Wh
  e <- eigen(tcrossprod(W) / m, symmetric = TRUE)
  kk <- min(k, max(n_spikes, 0L), m - 1L)
  if (kk < 1L) return(data.frame(k = integer(), MoranI = numeric(), p_order = numeric(),
                                 BC = numeric(), p_coexist = numeric(), p = numeric()))
  # Null fields: a grid of SAR strengths, matched later on Moran's I.
  rhos <- c(0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95)
  fields <- lapply(seq_along(rhos), function(i) sar_fields(A, rhos[i], n_null, seed + i))
  field_I <- lapply(fields, function(F) moran_I(F, A))
  field_BC <- lapply(fields, function(F) apply(F, 2, bimodality_coef))
  set.seed(seed)
  out <- lapply(seq_len(kk), function(j) {
    z <- e$vectors[, j]
    I0 <- moran_I(matrix(z), A)
    perm <- moran_I(sapply(seq_len(n_perm), function(b) sample(z)), A)
    p_order <- (1 + sum(perm >= I0)) / (n_perm + 1)
    bc0 <- bimodality_coef(z)
    # Null BCs from continuous fields whose Moran's I is near I0.
    allI <- unlist(field_I); allBC <- unlist(field_BC)
    near <- order(abs(allI - I0))[seq_len(n_null)]
    p_coexist <- (1 + sum(allBC[near] >= bc0)) / (n_null + 1)
    data.frame(k = j, MoranI = I0, p_order = p_order, BC = bc0, p_coexist = p_coexist,
               p = max(p_order, p_coexist))
  })
  do.call(rbind, out)
}
# Regime-level p-value: smallest spike p, Bonferroni over the spikes tested.
regime_phase_p <- function(tab) if (!nrow(tab)) 1 else min(1, min(tab$p) * nrow(tab))

# ---------------------------------------------------------------------------
# Number of coexisting phases. The spatially ordered spikes of a regime are
# analysed jointly: with r ordered spikes, the m x r matrix of unit loadings
# (columns scaled to unit variance) is split by k-means for k = 2..k_max and the
# between-cluster share of variance T_k is compared with its distribution under
# *continuous* ordered fields (independent Gaussian SAR fields on the same graph,
# each matched to the Moran's I of one loading column). Discrete sub-regimes
# (q coexisting phases) give well-separated clusters; a gradient or wave does
# not. p = min_k p_k * (k_max - 1); the estimated number of phases is the k with
# the largest standardized excess among the significant k (1 if none).
phase_count_test <- function(U, n_b, J, A, n_spikes, k_spikes = 5L, k_max = 4L, n_perm = 9999L,
                             n_null = 199L, alpha_order = 0.05, ridge = 1e-3, seed = 1L) {
  m <- nrow(U); nbar <- mean(n_b)
  A <- as.matrix(A); diag(A) <- 0
  empty <- list(r = 0L, p = 1, phases = 1L, spikes = data.frame(), Tk = data.frame())
  kk <- min(k_spikes, max(n_spikes, 0L), m - 1L)
  if (kk < 1L) return(empty)
  ev <- eigen(J * nbar, symmetric = TRUE)
  keep <- ev$values > ridge * max(ev$values)
  Wh <- ev$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(ev$values[keep]), sum(keep))
  W <- (U * sqrt(nbar / n_b)) %*% Wh
  e <- eigen(tcrossprod(W) / m, symmetric = TRUE)
  set.seed(seed)
  sp <- do.call(rbind, lapply(seq_len(kk), function(j) {
    z <- e$vectors[, j]
    I0 <- moran_I(matrix(z), A)
    perm <- moran_I(sapply(seq_len(n_perm), function(b) sample(z)), A)
    data.frame(k = j, MoranI = I0, p_order = (1 + sum(perm >= I0)) / (n_perm + 1))
  }))
  ord <- sp$k[sp$p_order < alpha_order]
  out <- list(r = length(ord), p = 1, phases = 1L, spikes = sp, Tk = data.frame())
  if (!length(ord)) return(out)
  L <- scale(e$vectors[, ord, drop = FALSE])
  ks <- 2:min(k_max, m - 1L)
  between <- function(X) vapply(ks, function(k) {
    km <- stats::kmeans(X, centers = k, nstart = 10L, iter.max = 50L)
    km$betweenss / km$totss
  }, numeric(1))
  T_obs <- between(L)
  # Continuous null: SAR fields matched column by column on Moran's I.
  rhos <- c(0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98)
  pool <- do.call(cbind, lapply(seq_along(rhos), function(i) sar_fields(A, rhos[i], n_null, seed + 10L * i)))
  pool_I <- moran_I(pool, A)
  near <- lapply(sp$MoranI[match(ord, sp$k)], function(I0) order(abs(pool_I - I0))[seq_len(2L * n_null)])
  T_null <- t(vapply(seq_len(n_null), function(b) {
    X <- scale(sapply(near, function(ix) pool[, ix[sample.int(length(ix), 1L)]]))
    between(matrix(X, nrow = m))
  }, numeric(length(ks))))
  if (length(ks) == 1L) T_null <- matrix(T_null, ncol = 1L)
  p_k <- vapply(seq_along(ks), function(i) (1 + sum(T_null[, i] >= T_obs[i])) / (n_null + 1), numeric(1))
  z_k <- (T_obs - colMeans(T_null)) / apply(T_null, 2, stats::sd)
  sig <- which(p_k < 0.05)
  out$p <- min(1, min(p_k) * length(ks))
  out$phases <- if (out$p < 0.05 && length(sig)) ks[sig[which.max(z_k[sig])]] else 1L
  out$Tk <- data.frame(k = ks, T_obs = T_obs, T_null_mean = colMeans(T_null), z = z_k, p = p_k)
  out
}

# ---------------------------------------------------------------------------
# Projected dip rule. Among the spatially ordered spikes of a regime, test the
# unit loadings for multimodality with Hartigan's dip test, which is calibrated
# against the uniform, the least favourable unimodal law (Hartigan & Hartigan
# 1985). A smooth gradient gives ordered but unimodal (at most uniform-like)
# loadings; coexisting discrete phases give multimodal loadings. Projections
# tested: each ordered loading column and, when r >= 2, each principal axis of
# the ordered loadings (cf. projected dip-means, Kalogeratos & Likas 2012;
# Chamalis & Likas 2018). p = min p * number of projections.
phase_dip_test <- function(U, n_b, J, A, n_spikes, k_spikes = 5L, n_perm = 9999L,
                           alpha_order = 0.05, ridge = 1e-3, seed = 1L) {
  m <- nrow(U); nbar <- mean(n_b)
  A <- as.matrix(A); diag(A) <- 0
  kk <- min(k_spikes, max(n_spikes, 0L), m - 1L)
  res <- list(r = 0L, p = 1, n_proj = 0L, spikes = data.frame())
  if (kk < 1L || m < 6L) return(res)
  ev <- eigen(J * nbar, symmetric = TRUE)
  keep <- ev$values > ridge * max(ev$values)
  Wh <- ev$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(ev$values[keep]), sum(keep))
  W <- (U * sqrt(nbar / n_b)) %*% Wh
  e <- eigen(tcrossprod(W) / m, symmetric = TRUE)
  set.seed(seed)
  sp <- do.call(rbind, lapply(seq_len(kk), function(j) {
    z <- e$vectors[, j]
    I0 <- moran_I(matrix(z), A)
    perm <- moran_I(sapply(seq_len(n_perm), function(b) sample(z)), A)
    data.frame(k = j, MoranI = I0, p_order = (1 + sum(perm >= I0)) / (n_perm + 1),
               p_dip = diptest::dip.test(z)$p.value)
  }))
  ord <- sp$k[sp$p_order < alpha_order]
  res$spikes <- sp; res$r <- length(ord)
  if (!length(ord)) return(res)
  L <- e$vectors[, ord, drop = FALSE]
  projs <- list(); for (j in seq_along(ord)) projs[[length(projs) + 1L]] <- L[, j]
  if (length(ord) >= 2L) {
    pc <- stats::prcomp(L, center = TRUE, scale. = TRUE)$x
    for (j in seq_len(ncol(pc))) projs[[length(projs) + 1L]] <- pc[, j]
  }
  pd <- vapply(projs, function(z) diptest::dip.test(z)$p.value, numeric(1))
  res$n_proj <- length(pd)
  res$p <- min(1, min(pd) * length(pd))
  res
}
