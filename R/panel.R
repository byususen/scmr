# Panel (unit x time) support: lagged-state design, forward filtering of class
# probabilities when the previous class is not observed, and forward-backward
# imputation of missing survey waves.
#
# The dynamic model is p_g(y_t = c | y_{t-1} = k, x_t): the previous class
# enters the local multinomial models as one-hot "lag" columns, so the
# transition intercepts gamma_{g,k,c} are ordinary (penalized) coefficients of
# those columns.

#' Lagged-class design for panel data
#'
#' Builds one indicator column per class for the class observed at the previous
#' time point of the same unit.
#' @param y Class labels (factor or vector), one per row.
#' @param unit_id Unit identifier per row.
#' @param time Integer-valued time index per row (for example months since start).
#' @param classes Class levels; defaults to `levels(factor(y))`.
#' @param prefix Prefix of the column names.
#' @return A list with `x` (rows x classes indicator matrix; rows without an
#'   observed previous time point are `NA`), `available` (logical: previous
#'   class observed at time - 1), `lag` (factor of previous classes) and
#'   `columns` (column names, in class order).
#' @export
scmr_lag_design <- function(y, unit_id, time, classes = NULL, prefix = "lag_") {
  y <- if (is.factor(y)) y else factor(y)
  classes <- classes %||% levels(y)
  y <- factor(as.character(y), levels = classes)
  n <- length(y)
  if (length(unit_id) != n || length(time) != n) stop("y, unit_id and time must have equal length.", call. = FALSE)
  if (any(!is.finite(time)) || any(time != floor(time))) stop("time must be integer valued.", call. = FALSE)
  key <- paste(unit_id, time, sep = "\r")
  if (anyDuplicated(key)) stop("Each unit may appear once per time point.", call. = FALSE)
  prev <- match(paste(unit_id, time - 1, sep = "\r"), key)
  lag <- factor(rep(NA_character_, n), levels = classes)
  ok <- !is.na(prev)
  lag[ok] <- y[prev[ok]]
  available <- ok & !is.na(lag)
  cols <- paste0(prefix, make.names(classes))
  x <- matrix(NA_real_, n, length(classes), dimnames = list(NULL, cols))
  x[available, ] <- 0
  x[cbind(which(available), as.integer(lag[available]))] <- 1
  list(x = x, available = available, lag = lag, columns = cols)
}

# Per-regime transition arrays P[row, c, k] = p_g(c | k, x_row).
transition_array <- function(engine, newx, lag_columns, classes) {
  C <- length(classes)
  out <- array(0, c(nrow(newx), C, C))
  for (k in seq_len(C)) {
    xk <- newx
    xk[, lag_columns] <- 0
    xk[, lag_columns[k]] <- 1
    out[, , k] <- engine_predict(engine, xk)
  }
  out
}

panel_order <- function(unit_id, time) {
  ord <- order(as.character(unit_id), time)
  split(ord, factor(as.character(unit_id)[ord], levels = unique(as.character(unit_id)[ord])))
}

check_lag_columns <- function(fit, lag_columns) {
  if (length(lag_columns) != length(fit$class_levels) || !all(lag_columns %in% fit$x_colnames)) {
    stop("lag_columns must name one fitted column per class, in class order.", call. = FALSE)
  }
}

# Starting class distribution: supplied per unit, else training class shares.
initial_distribution <- function(fit, init, units) {
  C <- length(fit$class_levels)
  base <- as.numeric(table(factor(fit$y_train, levels = fit$class_levels))) / length(fit$y_train)
  out <- matrix(base, length(units), C, byrow = TRUE, dimnames = list(units, fit$class_levels))
  if (is.null(init)) return(out)
  if (is.matrix(init)) {
    if (!all(units %in% rownames(init))) stop("init must have one row per prediction unit.", call. = FALSE)
    out[] <- init[units, fit$class_levels, drop = FALSE]
  } else {
    init <- init[!is.na(init)]
    hit <- intersect(names(init), units)
    for (u in hit) { out[u, ] <- 0; out[u, match(as.character(init[[u]]), fit$class_levels)] <- 1 }
  }
  out / rowSums(out)
}

#' Forward-filtered class probabilities when the previous class is unobserved
#'
#' For a dynamic fit (lag columns among the predictors) and a unit observed
#' over consecutive time points with features x_1, ..., x_T but no classes,
#' the exact predictive probabilities under the fitted model are obtained by
#' the forward recursion
#' alpha_t(c) = sum_k p_g(c | k, x_t) alpha_(t-1)(k), computed per cluster and
#' mixed with the unit's cluster weights. A gap in time restarts the recursion
#' from the initial distribution.
#' @param fit A fitted `scmr_fit` whose predictors include `lag_columns`.
#' @param newx Predictor matrix with the fitted columns; values in the lag
#'   columns are ignored.
#' @param unit_id,time Unit and integer time index per row.
#' @param new_coords Coordinates per row, needed for units not in training.
#' @param lag_columns Names of the lag indicator columns, in class order.
#' @param init Optional initial class per unit (named vector) or a matrix of
#'   initial probabilities (rows named by unit); defaults to training shares.
#' @param membership Cluster-weight rule for unseen units (see `predict()`).
#' @return A probability matrix (rows of `newx`, classes in columns).
#' @export
scmr_filter_predict <- function(fit, newx, unit_id, time, new_coords = NULL, lag_columns,
                                init = NULL, membership = c("potts", "proportion", "majority")) {
  membership <- match.arg(membership)
  check_lag_columns(fit, lag_columns)
  newx <- check_newx(newx, fit$x_colnames)
  n <- nrow(newx)
  C <- length(fit$class_levels)
  pp <- if (fit$model == "scmr") assign_test_PP(fit, unit_id, new_coords, n, membership) else
    matrix(1, n, 1L)
  trans <- lapply(fit$fits, transition_array, newx = newx, lag_columns = lag_columns,
                  classes = fit$class_levels)
  seqs <- panel_order(unit_id, time)
  a0 <- initial_distribution(fit, init, names(seqs))
  out <- matrix(0, n, C, dimnames = list(NULL, fit$class_levels))
  for (u in names(seqs)) {
    idx <- seqs[[u]]
    for (g in seq_along(trans)) {
      wg <- pp[idx[1], g]
      if (wg <= 0) next
      a <- a0[u, ]
      prev_t <- NA
      for (r in idx) {
        if (!is.na(prev_t) && time[r] != prev_t + 1) a <- a0[u, ]
        a <- as.numeric(trans[[g]][r, , ] %*% a)
        a <- a / sum(a)
        out[r, ] <- out[r, ] + wg * a
        prev_t <- time[r]
      }
    }
  }
  normalize_prob_matrix(out, fit$class_levels, n)
}

#' Impute missing survey waves by forward-backward smoothing
#'
#' Rows with `NA` in `y` are imputed from the observed classes of the same unit
#' (before and after the gap) and the features, under the fitted dynamic
#' model. Cluster weights are updated by the likelihood of the observed classes,
#' so the result is the exact posterior P(y_t | observed classes, x) of the
#' fitted model.
#' @inheritParams scmr_filter_predict
#' @param y Observed classes with `NA` for missing waves.
#' @return A probability matrix with one row per input row (observed rows are
#'   indicator rows) and attribute `cluster_posterior` (units x clusters).
#' @export
scmr_impute_waves <- function(fit, newx, y, unit_id, time, new_coords = NULL, lag_columns,
                              init = NULL, membership = c("potts", "proportion", "majority")) {
  membership <- match.arg(membership)
  check_lag_columns(fit, lag_columns)
  newx <- check_newx(newx, fit$x_colnames)
  n <- nrow(newx)
  C <- length(fit$class_levels)
  yi <- match(as.character(y), fit$class_levels)
  if (any(is.na(yi) & !is.na(y))) stop("y contains unknown classes.", call. = FALSE)
  pp <- if (fit$model == "scmr") assign_test_PP(fit, unit_id, new_coords, n, membership) else
    matrix(1, n, 1L)
  trans <- lapply(fit$fits, transition_array, newx = newx, lag_columns = lag_columns,
                  classes = fit$class_levels)
  seqs <- panel_order(unit_id, time)
  a0 <- initial_distribution(fit, init, names(seqs))
  G <- length(trans)
  out <- matrix(0, n, C, dimnames = list(NULL, fit$class_levels))
  post_g <- matrix(0, length(seqs), G, dimnames = list(names(seqs), paste0("G", seq_len(G))))
  evidence <- function(r) { e <- rep(1, C); if (!is.na(yi[r])) { e[] <- 0; e[yi[r]] <- 1 }; e }
  for (u in names(seqs)) {
    idx <- seqs[[u]]
    Tn <- length(idx)
    restart <- c(TRUE, diff(time[idx]) != 1)
    logw <- rep(-Inf, G)
    marg <- vector("list", G)
    for (g in seq_len(G)) {
      wg <- pp[idx[1], g]
      if (wg <= 0) next
      f <- matrix(0, Tn, C)
      loglik <- 0
      for (t in seq_len(Tn)) {
        prior <- if (restart[t]) a0[u, ] else f[t - 1, ]
        a <- as.numeric(trans[[g]][idx[t], , ] %*% prior) * evidence(idx[t])
        z <- sum(a)
        loglik <- loglik + log(max(z, .Machine$double.xmin))
        f[t, ] <- if (z > 0) a / z else rep(1 / C, C)
      }
      b <- matrix(1, Tn, C)
      if (Tn > 1) for (t in (Tn - 1):1) {
        if (restart[t + 1]) { b[t, ] <- 1; next }
        v <- as.numeric(crossprod(trans[[g]][idx[t + 1], , ], evidence(idx[t + 1]) * b[t + 1, ]))
        b[t, ] <- v / max(sum(v), .Machine$double.xmin)
      }
      m <- f * b
      marg[[g]] <- m / rowSums(m)
      logw[g] <- log(wg) + loglik
    }
    r <- exp(logw - max(logw))
    r <- r / sum(r)
    post_g[u, ] <- r
    for (g in which(r > 0)) out[idx, ] <- out[idx, ] + r[g] * marg[[g]]
  }
  out <- normalize_prob_matrix(out, fit$class_levels, n)
  attr(out, "cluster_posterior") <- post_g
  out
}
