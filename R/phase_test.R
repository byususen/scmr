# Two-stage phase test for the regimes of a fitted SCMR model.
#
# For regime g with units b = 1..m, the unit score vector is
#   u_b = sum_{i in b} (e_{y_i} - p_i)_{1..C-1} (x) (1, x_i)
# at the fitted probabilities p_i, and J is the average per-row Fisher
# information of the regime. The whitened scores
#   w_b = (nbar / n_b)^{1/2} (nbar J)^{-1/2} u_b
# have covariance close to sigma^2 I under a homogeneous regime, so the
# eigenvalues of S = (1/m) sum_b w_b w_b' form a Marchenko-Pastur bulk; a hidden
# sub-regime adds an eigenvalue that leaves the bulk (the BBP phase transition).
#
# Stage 0 (spikes): R_k = lambda_k / median(lambda) is compared with the
#   bootstrap distribution of R_1 under the fitted probabilities; a common
#   inflation of the score variance (within-unit dependence) cancels in R.
# Stage 1 (order): each spike eigenvector gives unit loadings; Moran's I on the
#   unit graph with a permutation p-value. Ordered spikes indicate spatial
#   non-stationarity inside the regime; unstructured unit effects do not order.
# Stage 2 (discreteness): Hartigan's dip test of multimodality on each ordered
#   loading vector and on each principal axis of the ordered loadings. The dip
#   test is calibrated against the uniform, the least favourable unimodal law,
#   so smooth gradients (ordered but unimodal loadings) are not declared
#   discrete; coexisting discrete phases give multimodal loadings.
# Within a regime the stage-2 p-value is Bonferroni-adjusted over projections;
# across the regimes of the fit it is Holm-adjusted.

phase_unit_scores <- function(x, y, prob, unit_num, rows) {
  Z <- cbind(1, x[rows, , drop = FALSE])
  P <- prob[rows, , drop = FALSE]
  C <- ncol(P); K <- C - 1L; p1 <- ncol(Z)
  Y <- matrix(0, length(rows), C)
  Y[cbind(seq_along(rows), as.integer(y[rows]))] <- 1
  un <- factor(unit_num[rows])
  U <- rowsum(do.call(cbind, lapply(seq_len(K), function(c) Z * (Y[, c] - P[, c]))), un, reorder = TRUE)
  J <- matrix(0, K * p1, K * p1)
  for (c1 in seq_len(K)) for (c2 in c1:K) {
    w <- if (c1 == c2) P[, c1] * (1 - P[, c1]) else -P[, c1] * P[, c2]
    blk <- crossprod(Z * w, Z)
    J[(c1 - 1) * p1 + seq_len(p1), (c2 - 1) * p1 + seq_len(p1)] <- blk
    J[(c2 - 1) * p1 + seq_len(p1), (c1 - 1) * p1 + seq_len(p1)] <- t(blk)
  }
  list(U = U, n_b = as.numeric(table(un)), J = J / length(rows), units = as.integer(levels(un)))
}

phase_whitener <- function(J, nbar, ridge = 1e-3) {
  ev <- eigen(J * nbar, symmetric = TRUE)
  keep <- ev$values > ridge * max(ev$values)
  ev$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(ev$values[keep]), sum(keep))
}

phase_eigen <- function(U, n_b, Wh, vectors = FALSE) {
  m <- nrow(U); nbar <- mean(n_b)
  W <- (U * sqrt(nbar / n_b)) %*% Wh
  e <- eigen(tcrossprod(W) / m, symmetric = TRUE, only.values = !vectors)
  r <- seq_len(min(m - 1L, ncol(W)))            # drop the centring zero
  list(values = e$values[r], vectors = if (vectors) e$vectors[, r, drop = FALSE] else NULL)
}

phase_moran <- function(Z, A) {
  Z <- sweep(as.matrix(Z), 2, colMeans(as.matrix(Z)))
  (nrow(Z) / sum(A)) * colSums(Z * (A %*% Z)) / colSums(Z^2)
}

phase_knn <- function(xy, k) {
  D <- as.matrix(stats::dist(xy)); diag(D) <- Inf
  A <- matrix(0, nrow(xy), nrow(xy))
  for (i in seq_len(nrow(xy))) A[i, order(D[i, ])[seq_len(min(k, nrow(xy) - 1L))]] <- 1
  pmax(A, t(A))
}

phase_regime <- function(x, y, prob, unit_num, rows, xy, B, n_perm, k_neighbors, max_spikes,
                         alpha, ridge, seed) {
  sc <- phase_unit_scores(x, y, prob, unit_num, rows)
  m <- nrow(sc$U)
  out <- list(Units = m, Spikes = 0L, OrderedSpikes = 0L, R1 = NA_real_, Rcrit = NA_real_,
              MaxMoran = NA_real_, p_dip = 1, spikes = data.frame())
  if (m < 6L) return(out)
  Wh <- phase_whitener(sc$J, mean(sc$n_b), ridge)
  e <- phase_eigen(sc$U, sc$n_b, Wh, vectors = TRUE)
  R <- e$values / stats::median(e$values)
  # Stage 0: bootstrap calibration of the largest ratio under the fitted model.
  P <- prob[rows, , drop = FALSE]
  cum <- t(apply(P, 1, cumsum))
  lev <- levels(y)
  Rb <- with_scmr_seed(seed, replicate(B, {
    yy <- y
    yy[rows] <- factor(lev[1 + rowSums(stats::runif(nrow(P)) > cum)], levels = lev)
    sb <- phase_unit_scores(x, yy, prob, unit_num, rows)
    lb <- phase_eigen(sb$U, sb$n_b, Wh)$values
    lb[1] / stats::median(lb)
  }))
  crit <- stats::quantile(Rb, 1 - alpha, names = FALSE)
  n_spk <- sum(R[seq_len(min(max_spikes, length(R)))] > crit)
  out$Spikes <- n_spk; out$R1 <- R[1]; out$Rcrit <- crit
  if (!n_spk) return(out)
  # Stage 1: spatial order of the spike loadings.
  A <- phase_knn(xy[sc$units, , drop = FALSE], k_neighbors)
  L <- e$vectors[, seq_len(n_spk), drop = FALSE]
  I0 <- phase_moran(L, A)
  p_order <- with_scmr_seed(seed + 1L, vapply(seq_len(n_spk), function(j) {
    perm <- phase_moran(sapply(seq_len(n_perm), function(b) sample(L[, j])), A)
    (1 + sum(perm >= I0[j])) / (n_perm + 1)
  }, numeric(1)))
  out$spikes <- data.frame(k = seq_len(n_spk), R = R[seq_len(n_spk)], MoranI = I0, p_order = p_order)
  ord <- which(p_order < alpha)
  out$OrderedSpikes <- length(ord); out$MaxMoran <- max(I0)
  if (!length(ord)) return(out)
  # Stage 2: multimodality of the ordered loadings (projected dip).
  projs <- lapply(ord, function(j) L[, j])
  if (length(ord) >= 2L) {
    pc <- stats::prcomp(L[, ord, drop = FALSE], center = TRUE, scale. = TRUE)$x
    projs <- c(projs, lapply(seq_len(ncol(pc)), function(j) pc[, j]))
  }
  pd <- vapply(projs, function(z) diptest::dip.test(z)$p.value, numeric(1))
  out$p_dip <- min(1, min(pd) * length(pd))
  out
}

#' Two-stage phase test: are the regimes of a fit internally homogeneous,
#' continuously varying, or hiding discrete sub-regimes?
#'
#' For every regime of a fitted SCMR model, the unit score vectors at the
#' fitted probabilities are whitened by the regime's Fisher information and
#' their eigenvalues compared with the noise bulk (bootstrap-calibrated spikes).
#' Spikes whose unit loadings are spatially ordered (Moran's I, permutation
#' test) indicate spatial non-stationarity inside the regime; if those loadings
#' are also multimodal (Hartigan's dip test, calibrated against the uniform),
#' the regime contains discrete sub-regimes. See Details for the statistics.
#'
#' Each regime is classified as `"homogeneous"` (no spatially ordered spike),
#' `"continuous"` (ordered spikes with unimodal loadings: smooth variation that
#' further regimes would only approximate) or `"discrete"` (Holm-adjusted dip
#' p-value below `alpha`: a split is warranted). The test needs enough rows per
#' unit for the loadings to resolve; with few rows per unit its power is low.
#'
#' @param fit A fitted `scmr_fit` (model `"scmr"`, or G = 1) on `x`, `y`.
#' @param x,y,unit_id,coords The data the model was fitted to (row-level).
#' @param B Bootstrap replicates for the spike threshold.
#' @param n_perm Permutations for Moran's I.
#' @param k_neighbors Nearest neighbours for the unit graph.
#' @param max_spikes Largest number of spikes examined per regime.
#' @param alpha Level of every stage.
#' @param seed Random seed.
#' @return A data frame with one row per regime: units, spikes, ordered
#'   spikes, largest spike ratio and its threshold, largest Moran's I, dip
#'   p-value (raw and Holm-adjusted) and the classification. The spike-level
#'   results are attached as attribute `"spikes"`.
#' @export
scmr_phase_test <- function(fit, x, y, unit_id, coords, B = 200L, n_perm = 9999L, k_neighbors = 8L,
                            max_spikes = 5L, alpha = 0.05, seed = 1L) {
  if (!requireNamespace("diptest", quietly = TRUE)) {
    stop("Package 'diptest' is required: install.packages(\"diptest\").", call. = FALSE)
  }
  x <- as.matrix(x)
  xs <- scale(x); xs[!is.finite(xs)] <- 0
  y <- factor(y, levels = colnames(fit$prob_train) %||% levels(factor(y)))
  units <- sort(unique(as.character(unit_id)))
  unit_num <- match(as.character(unit_id), units)
  coords <- as.matrix(coords)
  xy <- if (nrow(coords) == length(unit_id)) coords[match(units, as.character(unit_id)), , drop = FALSE] else coords
  ug <- as.integer(fit$group_unit[units])
  if (anyNA(ug)) ug <- rep(1L, length(units))
  regimes <- sort(unique(ug))
  res <- lapply(regimes, function(g) {
    rows <- which(ug[unit_num] == g)
    phase_regime(xs, y, fit$prob_train, unit_num, rows, xy, B, n_perm, k_neighbors, max_spikes,
                 alpha, 1e-3, seed + 100L * g)
  })
  tab <- data.frame(Regime = regimes,
                    Units = vapply(res, function(r) r$Units, numeric(1)),
                    Spikes = vapply(res, function(r) r$Spikes, numeric(1)),
                    OrderedSpikes = vapply(res, function(r) r$OrderedSpikes, numeric(1)),
                    R1 = vapply(res, function(r) r$R1, numeric(1)),
                    Rcrit = vapply(res, function(r) r$Rcrit, numeric(1)),
                    MaxMoran = vapply(res, function(r) r$MaxMoran, numeric(1)),
                    p_dip = vapply(res, function(r) r$p_dip, numeric(1)))
  tab$p_dip_holm <- stats::p.adjust(tab$p_dip, "holm")
  tab$Phase <- ifelse(tab$p_dip_holm < alpha, "discrete",
                      ifelse(tab$OrderedSpikes > 0, "continuous", "homogeneous"))
  attr(tab, "spikes") <- do.call(rbind, lapply(seq_along(res), function(i)
    if (nrow(res[[i]]$spikes)) data.frame(Regime = regimes[i], res[[i]]$spikes) else NULL))
  tab
}
