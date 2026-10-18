# Potts label model: neighbour support, exact local changes of the Besag
# pseudo-log-likelihood, and the maximum pseudo-likelihood estimate of phi.
#
#   log PL(g; phi) = sum_r [ phi * S_r(g_r) - log sum_h exp(phi * S_r(h)) ],
#   S_r(h)        = sum_q w_rq I(g_q = h).
#
# The label term of the penalized objective is either phi * (sum of agreeing
# edge weights) (phi fixed, the original SCR penalty) or log PL(g; phi) (phi
# estimated). Both changes caused by moving one unit only involve that unit and
# its neighbours, so membership updates stay exact and cheap.

# Neighbour support matrix S (units x G).
potts_support <- function(w, groups, G) {
  onehot <- matrix(0, length(groups), G)
  onehot[cbind(seq_along(groups), groups)] <- 1
  as.matrix(w %*% onehot)
}

log_sum_exp <- function(a) {
  m <- max(a)
  m + log(sum(exp(a - m)))
}

# Row-wise sparse neighbour lists of a symmetric weight matrix.
potts_neighbors <- function(w) {
  w <- methods::as(w, "generalMatrix")
  s <- Matrix::summary(w)
  s <- s[s$x != 0 & s$i != s$j, , drop = FALSE]
  idx <- split(s$j, factor(s$i, levels = seq_len(nrow(w))))
  val <- split(s$x, factor(s$i, levels = seq_len(nrow(w))))
  list(index = unname(idx), weight = unname(val))
}

# Change of log PL when unit u moves from groups[u] to h, given support S.
potts_delta_logpl <- function(u, h, groups, S, nb, phi) {
  old <- groups[u]
  if (h == old) return(0)
  # Unit u's own term: its support row does not depend on its own label.
  delta <- phi * (S[u, h] - S[u, old])
  # (the normaliser of unit u's own term is unchanged)
  jj <- nb$index[[u]]
  if (!length(jj)) return(delta)
  ww <- nb$weight[[u]]
  for (a in seq_along(jj)) {
    r <- jj[a]
    s_old <- S[r, ]
    s_new <- s_old
    s_new[old] <- s_new[old] - ww[a]
    s_new[h] <- s_new[h] + ww[a]
    gr <- groups[r]
    delta <- delta + phi * (s_new[gr] - s_old[gr]) -
      (log_sum_exp(phi * s_new) - log_sum_exp(phi * s_old))
  }
  delta
}

# Update support S after unit u moved from old to h.
potts_update_support <- function(S, u, old, h, nb) {
  jj <- nb$index[[u]]
  if (length(jj)) {
    S[jj, old] <- S[jj, old] - nb$weight[[u]]
    S[jj, h] <- S[jj, h] + nb$weight[[u]]
  }
  S
}

potts_logpl_from_support <- function(S, groups, phi) {
  a <- phi * S
  amax <- apply(a, 1L, max)
  sum(a[cbind(seq_along(groups), groups)] - amax - log(rowSums(exp(a - amax))))
}

# Score and observed information of log PL in phi.
potts_logpl_derivatives <- function(S, groups, phi) {
  a <- phi * S
  a <- a - apply(a, 1L, max)
  p <- exp(a)
  p <- p / rowSums(p)
  m1 <- rowSums(p * S)
  m2 <- rowSums(p * S^2)
  list(score = sum(S[cbind(seq_along(groups), groups)] - m1),
       info = sum(m2 - m1^2))
}

#' Maximum pseudo-likelihood estimate of the Potts interaction
#'
#' log PL(g; phi) is concave in phi (its second derivative is minus the sum of
#' the conditional variances of the neighbour support), so the maximiser on
#' `[0, phi_max]` is unique and is found by safeguarded Newton steps.
#' @param w Symmetric spatial weight matrix between units.
#' @param groups Integer unit labels in `1:G`.
#' @param G Number of clusters.
#' @param phi_start Starting value.
#' @param phi_max Upper bound (the estimate diverges when every unit agrees
#'   with all of its neighbours).
#' @param tol Convergence tolerance on the score.
#' @return The estimate, with attributes `logpl` and `iterations`.
#' @export
scmr_potts_phi <- function(w, groups, G = max(groups), phi_start = 1, phi_max = 20, tol = 1e-8) {
  groups <- as.integer(groups)
  if (G <= 1L) return(structure(0, logpl = 0, iterations = 0L))
  S <- potts_support(w, groups, G)
  f <- function(phi) potts_logpl_from_support(S, groups, phi)
  phi <- min(max(phi_start, 0), phi_max)
  it <- 0L
  for (it in seq_len(100L)) {
    d <- potts_logpl_derivatives(S, groups, phi)
    if (abs(d$score) < tol * max(1, length(groups))) break
    step <- if (d$info > 1e-12) d$score / d$info else sign(d$score)
    cand <- min(max(phi + step, 0), phi_max)
    # Concavity: back-track until the objective does not decrease.
    f0 <- f(phi)
    while (f(cand) < f0 - 1e-12 && abs(cand - phi) > 1e-10) cand <- (cand + phi) / 2
    if (abs(cand - phi) < 1e-10) break
    phi <- cand
  }
  structure(phi, logpl = f(phi), iterations = it)
}
