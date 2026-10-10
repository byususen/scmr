#' Original spatially clustered regression (benchmark)
#'
#' A multinomial port of the SCR algorithm of Sugasawa and Murakami (2021),
#' following their reference code (github.com/sshonosuke/SCR, `SCR()`):
#' unpenalized local models (`nnet::multinom` here, `lm`/`glm` there), a fixed
#' spatial weight `phi`, initial groups from the best of `kmeans_starts`
#' k-means runs on the coordinates (smallest total within-cluster sum of
#' squares), and a *simultaneous* update of all labels given the previous
#' labels, \eqn{g_i \leftarrow \arg\max_g \{\ell_i(g) + \phi \sum_j w_{ij} I(g_j = g)\}}.
#' A group whose size is at most p + 1 (or that lacks a class) keeps its
#' previous coefficients. Iterations stop when the objective
#' \eqn{\sum_i \max_g Q_i(g)} changes by less than `tol` (relative) or returns
#' to its running maximum, as in the reference code; the simultaneous update is
#' not guaranteed to be monotone.
#'
#' Differences from the reference code, needed for panel data and multinomial
#' responses: labels are assigned per spatial unit (all rows of a unit share a
#' label, so the unit score sums its rows' log-likelihoods); groups that become
#' empty are dropped. With one row per unit the update is the original one.
#' The original BIC is `criteria$CriterionSCR_BIC_original`.
#' @inheritParams fit_scmr
#' @param phi Fixed spatial weight (the reference code uses 1).
#' @param k_neighbors Number of nearest units in the binary symmetric weights.
#' @param kmeans_starts Number of k-means runs for the initial groups.
#' @param maxitr Maximum number of iterations.
#' @param tol Relative convergence tolerance.
#' @return A `scmr_fit` object (`model = "scmr"`, `penalty = "none"`) with
#'   `scr_original = TRUE`; prediction, filtering and imputation functions apply.
#' @export
fit_scr_original <- function(x, y, unit_id, coords, G, phi = 1, k_neighbors = 5,
                             kmeans_starts = 20, maxitr = 100, tol = 1e-5,
                             control = scmr_control(), seed = 123) {
  data <- prepare_scmr_data(x, y, unit_id, coords)
  control <- utils::modifyList(control, list(phi = phi, k_neighbors = k_neighbors,
                                             weight_type = "binary", phi_update = "fixed",
                                             row_standardize_weights = FALSE))
  classes <- levels(data$y)
  units <- sort(unique(data$unit_id))
  unit_num <- match(data$unit_id, units)
  ucoords <- get_unit_coords(data$unit_id, data$coords)
  w <- build_unit_weights(ucoords, k_neighbors, "binary", NULL, TRUE, FALSE)$W
  p <- ncol(data$x) + 1L
  fit_engine <- function(idx) tryCatch(
    fit_local_engine(data$x[idx, , drop = FALSE], factor(data$y[idx], levels = classes),
                     "none", NA_real_, NA_real_, control),
    error = function(e) NULL)
  t0 <- proc.time()[["elapsed"]]
  out <- with_scmr_seed(seed, {
    runs <- lapply(seq_len(kmeans_starts), function(k) stats::kmeans(ucoords, G))
    groups <- runs[[which.min(vapply(runs, function(r) r$tot.withinss, numeric(1)))]]$cluster
    global <- fit_engine(seq_len(nrow(data$x)))
    if (is.null(global)) stop("The unpenalized global multinomial model failed.", call. = FALSE)
    engines <- rep(list(global), G)
    val <- 0
    mval <- -Inf
    iterations <- list()
    converged <- FALSE
    for (it in seq_len(maxitr)) {
      cval <- val
      # Parameter step: refit every group that is large enough.
      for (g in seq_len(G)) {
        idx <- which(groups[unit_num] == g)
        if (length(idx) > p + 1L && all(table(factor(data$y[idx], levels = classes)) > 0L)) {
          e <- fit_engine(idx)
          if (!is.null(e)) engines[[g]] <- e
        }
      }
      # Simultaneous assignment step given the previous labels.
      row_ll <- vapply(engines, function(e) true_class_logp(engine_predict(e, data$x), data$y, classes),
                       numeric(nrow(data$x)))
      row_ll <- matrix(row_ll, nrow(data$x), G)
      unit_ll <- rowsum(row_ll, unit_num, reorder = TRUE)
      Q <- unit_ll + phi * potts_support(w, groups, G)
      new_groups <- max.col(Q, ties.method = "first")
      val <- sum(apply(Q, 1L, max))
      iterations[[it]] <- data.frame(Start = 1L, Iter = it, LogLik = sum(row_ll[cbind(seq_len(nrow(data$x)), new_groups[unit_num])]),
                                     SpatialBonus = spatial_bonus(w, new_groups), Phi = phi,
                                     LabelTerm = phi * spatial_bonus(w, new_groups), Objective = val,
                                     PenalizedObjective = val, Changes = sum(new_groups != groups),
                                     AcceptedChanges = sum(new_groups != groups), Reverted = FALSE)
      groups <- new_groups
      dd <- abs(cval - val) / abs(val)
      mval <- max(mval, cval)
      if (dd < tol || abs(mval - val) < tol) { converged <- TRUE; break }
    }
    list(groups = groups, engines = engines, iterations = iterations, converged = converged, it = it)
  })
  # Drop groups left empty by the simultaneous update.
  used <- sort(unique(out$groups))
  groups_unit <- match(out$groups, used)
  fits <- out$engines[used]
  for (g in seq_along(fits)) fits[[g]]$n_fit <- sum(groups_unit[unit_num] == g)
  ans <- new_scmr_fit("scmr", "none", fits, data, groups_unit[unit_num], control, seed,
                      iterations = bind_diagnostics(out$iterations), membership_converged = out$converged,
                      reason = if (out$converged) "scr_objective_converged" else "max_iter_reached",
                      w = w, bandwidth = NA_real_)
  ans$scr_original <- TRUE
  ans$G_requested <- G
  ans$phi <- phi
  ans$lambda_rule <- "none"
  ans$runtime <- list(total_sec = proc.time()[["elapsed"]] - t0, tuning_sec = 0)
  ans
}
