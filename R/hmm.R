# Hybrid HMM smoothing of static class probabilities.
#
# A static model gives p(y_t | x_t). Under the hidden Markov assumption
# (x_t independent of the rest given y_t), p(x_t | y_t) is proportional to
# p(y_t | x_t) / pi(y_t), so with a phase transition matrix A the posterior of
# the whole sequence is
#   P(y_1:T | x_1:T) proportional to pi(y_1) e_1(y_1) prod_t A(y_t-1, y_t) e_t(y_t),
#   e_t(c) = {p(y_t = c | x_t) / pi(c)}^tau.
# tau < 1 tempers the emissions when the features of neighbouring months share
# information (for example lagged satellite indices). Marginals come from the
# forward-backward recursion ("smooth") or the forward recursion only
# ("filter", uses x up to t). The result can be tempered (temperature) and
# pooled linearly with the static probabilities (weight).

#' Phase transition matrix from observed class sequences
#'
#' Counts transitions between consecutive time points of the same unit.
#' @param y Classes (factor or character).
#' @param unit_id,time Unit and integer time index per row.
#' @param classes Class levels; defaults to `levels(y)`.
#' @param pseudo Pseudo-count added to every cell.
#' @return A list with `transition` (rows: previous class, columns: current
#'   class, rows sum to one) and `prior` (class frequencies).
#' @export
scmr_transition_matrix <- function(y, unit_id, time, classes = NULL, pseudo = 0.5) {
  y <- as.character(y)
  if (is.null(classes)) classes <- if (is.factor(y)) levels(y) else sort(unique(y))
  classes <- as.character(classes)
  n <- length(y)
  if (length(unit_id) != n || length(time) != n) stop("y, unit_id and time must have equal length.", call. = FALSE)
  key <- paste(unit_id, time, sep = "\r")
  prev <- match(paste(unit_id, time - 1, sep = "\r"), key)
  ok <- !is.na(prev) & !is.na(y) & !is.na(y[prev])
  A <- table(factor(y[prev[ok]], levels = classes), factor(y[ok], levels = classes))
  A <- unclass(A) + pseudo
  A <- A / rowSums(A)
  dimnames(A) <- list(previous = classes, current = classes)
  prior <- as.numeric(table(factor(y, levels = classes)) + pseudo)
  list(transition = A, prior = stats::setNames(prior / sum(prior), classes))
}

#' Smooth static class probabilities over time with a hybrid HMM
#'
#' @param prob Probability matrix (rows, classes) from a static model.
#' @param unit_id,time Unit and integer time index per row; a gap in time
#'   restarts the chain.
#' @param transition Transition matrix (previous class in rows), for example
#'   from [scmr_transition_matrix()].
#' @param prior Class prior of the static model's training data.
#' @param tau Emission exponent (1: plain hybrid HMM).
#' @param temperature Temperature applied to the HMM marginals.
#' @param weight Weight of the HMM marginals in the linear pool with `prob`
#'   (1: HMM only, 0: `prob` unchanged).
#' @param direction `"smooth"` (forward-backward) or `"filter"` (forward only).
#' @return A probability matrix like `prob`.
#' @export
scmr_hmm_smooth <- function(prob, unit_id, time, transition, prior, tau = 1, temperature = 1,
                            weight = 1, direction = c("smooth", "filter")) {
  direction <- match.arg(direction)
  prob <- as.matrix(prob)
  n <- nrow(prob)
  C <- ncol(prob)
  if (!all(dim(transition) == C) || length(prior) != C) {
    stop("transition must be C x C and prior of length C.", call. = FALSE)
  }
  if (length(unit_id) != n || length(time) != n) stop("unit_id and time must match the rows of prob.", call. = FALSE)
  A <- unname(as.matrix(transition))
  prior <- as.numeric(prior)
  E <- sweep(pmax(prob, 1e-300), 2, prior, "/")^tau
  out <- matrix(0, n, C)
  for (idx in panel_order(unit_id, time)) {
    Tn <- length(idx)
    restart <- c(TRUE, diff(time[idx]) != 1)
    f <- matrix(0, Tn, C)
    for (t in seq_len(Tn)) {
      pr <- if (restart[t]) prior else as.numeric(f[t - 1, ] %*% A)
      a <- pr * E[idx[t], ]
      f[t, ] <- a / sum(a)
    }
    m <- f
    if (direction == "smooth" && Tn > 1) {
      b <- matrix(1, Tn, C)
      for (t in (Tn - 1):1) {
        if (restart[t + 1]) next
        v <- as.numeric(A %*% (E[idx[t + 1], ] * b[t + 1, ]))
        b[t, ] <- v / sum(v)
      }
      m <- f * b
      m <- m / rowSums(m)
    }
    out[idx, ] <- m
  }
  if (temperature != 1) {
    l <- log(pmax(out, 1e-300)) / temperature
    l <- l - apply(l, 1, max)
    out <- exp(l) / rowSums(exp(l))
  }
  out <- (1 - weight) * prob + weight * out
  out <- out / rowSums(out)
  dimnames(out) <- dimnames(prob)
  out
}
