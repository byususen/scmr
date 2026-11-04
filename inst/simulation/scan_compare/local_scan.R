# Copied unchanged from the Stage 26 working package (Bahan/Syntax/scmr_pkg/R) for the scan comparison.
# Local offset likelihood used only to rank candidate starting regions.

#' Add local multinomial offset likelihood gains to a candidate pool
#'
#' For a fixed global fit and region W, maximizes the multinomial log
#' likelihood after adding a K=C-1 dimensional intercept shift inside W.
#' A small ridge term stabilizes sparse regions. This score is not calibrated
#' for the search over regions and must not be used as a G-selection test.
#'
#' @param scan Pool returned by `scmr_flexible_scan_starts()`.
#' @param y Ordered categorical response used by the global fit.
#' @param unit_id Ordered spatial unit ID for each response.
#' @param global_fit One-cluster fit on the same ordered data.
#' @param ridge Nonnegative ridge term on the local intercept shift.
#' @param max_iter Maximum Newton steps per region.
#' @param tol Absolute convergence tolerance for the Newton step.
#' @return Candidate pool with an added `intercept_local` score column.
#' @export
scmr_add_local_likelihood_scores <- function(scan, y, unit_id, global_fit,
                                              ridge = 0.1, max_iter = 15L,
                                              tol = 1e-7) {
  if (!is.list(scan) || is.null(scan$regions) ||
      is.null(scan$all_scores) || is.null(scan$unit_ids))
    stop("scan must be a flexible candidate pool.", call. = FALSE)
  if (!inherits(global_fit, "scmr_fit") || global_fit$G != 1L ||
      !identical(as.character(global_fit$y_train), as.character(y)) ||
      !identical(as.character(global_fit$unit_id), as.character(unit_id)))
    stop("global_fit, y and unit_id must refer to the same ordered data.",
         call. = FALSE)
  if (length(ridge) != 1L || !is.finite(ridge) || ridge < 0 ||
      length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter < 1L || max_iter != floor(max_iter) ||
      length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop("Invalid optimizer controls.", call. = FALSE)
  p <- global_fit$prob_train
  if (nrow(p) != length(y) || ncol(p) < 3L ||
      !identical(colnames(p), levels(factor(y))))
    stop("Global probabilities must match the response classes.", call. = FALSE)
  p <- pmax(p, 1e-12)
  p <- p / rowSums(p)
  logp <- log(p)
  C <- ncol(p)
  K <- C - 1L
  yi <- match(as.character(y), colnames(p))
  if (anyNA(yi)) stop("Response has unknown class.", call. = FALSE)
  rows_by_unit <- split(seq_along(unit_id), as.character(unit_id))
  if (!all(scan$unit_ids %in% names(rows_by_unit)))
    stop("Candidate IDs do not match unit_id.", call. = FALSE)
  local_gain <- function(region) {
    idx <- unlist(rows_by_unit[as.character(region)], use.names = FALSE)
    if (!length(idx)) return(NA_real_)
    a <- logp[idx, , drop = FALSE]
    yloc <- yi[idx]
    observed <- tabulate(yloc, nbins = C)[seq_len(K)]
    calc <- function(delta) {
      z <- cbind(sweep(a[, seq_len(K), drop = FALSE], 2L, delta, "+"),
                 a[, C])
      shift <- apply(z, 1L, max)
      e <- exp(z - shift)
      den <- rowSums(e)
      q <- e[, seq_len(K), drop = FALSE] / den
      value <- sum(observed * delta) - sum(shift + log(den)) -
        0.5 * ridge * sum(delta^2)
      gradient <- observed - colSums(q) - ridge * delta
      information <- diag(colSums(q), nrow = K) - crossprod(q) +
        diag(ridge, nrow = K)
      list(value = value, gradient = gradient, information = information)
    }
    delta <- numeric(K)
    current <- calc(delta)
    baseline <- current$value
    for (iter in seq_len(max_iter)) {
      step <- tryCatch(solve(current$information, current$gradient),
                       error = function(e) rep(NA_real_, K))
      if (any(!is.finite(step)) || max(abs(step)) < tol) break
      factor <- 1
      repeat {
        proposal <- calc(delta + factor * step)
        if (is.finite(proposal$value) &&
            proposal$value >= current$value - 1e-10) break
        factor <- factor / 2
        if (factor < 1 / 128) break
      }
      if (factor < 1 / 128 || proposal$value < current$value) break
      delta <- delta + factor * step
      current <- proposal
      if (max(abs(factor * step)) < tol) break
    }
    2 * max(0, current$value - baseline)
  }
  gain <- vapply(scan$regions, local_gain, numeric(1L))
  if (any(!is.finite(gain)))
    stop("Local likelihood optimization returned nonfinite scores.",
         call. = FALSE)
  scan$all_scores <- cbind(scan$all_scores, intercept_local = gain)
  scan$local_likelihood <- list(ridge = ridge, max_iter = max_iter, tol = tol)
  scan
}
