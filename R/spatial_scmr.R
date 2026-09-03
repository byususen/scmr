# One membership algorithm for penalized and unpenalized local engines.
membership_sweep <- function(fits, data, groups, unit_levels, w, control, seed) {
  G <- length(fits)
  unit_num <- match(data$unit_id, unit_levels)
  rows <- split(seq_len(nrow(data$x)), factor(unit_num, levels = seq_along(unit_levels)))
  contrib <- precompute_unit_contributions(rows, data$y, levels(data$y))
  counts <- compute_cluster_counts(groups, G, contrib$unit_total, contrib$unit_class)
  min_class <- resolve_min_class(control$min_per_class, levels(data$y))
  row_ll <- vapply(fits, function(f) true_class_logp(engine_predict(f, data$x), data$y, levels(data$y)),
                   numeric(nrow(data$x)))
  row_ll <- matrix(row_ll, nrow(data$x), G)
  unit_ll <- matrix(0, length(unit_levels), G)
  # Aggregate explicitly in unit-level order; input rows need not be sorted.
  for (u in seq_along(rows)) unit_ll[u, ] <- colSums(row_ll[rows[[u]], , drop = FALSE])
  order_u <- if (control$update_order == "random") with_scmr_seed(seed, sample.int(length(groups))) else seq_along(groups)
  changes <- 0L
  for (u in order_u) {
    old <- groups[u]
    weights <- as.numeric(w[u, ])
    nn <- which(weights != 0)
    neighbor <- vapply(seq_len(G), function(g) sum(weights[nn][groups[nn] == g]), numeric(1))
    score <- unit_ll[u, ] + control$phi * neighbor
    for (candidate in order(score, decreasing = TRUE)) {
      if (score[candidate] <= score[old] + 1e-12) break
      if (!move_is_feasible(u, candidate, groups, counts$unit_count, counts$class_count,
                             contrib$unit_class, control$min_units, min_class)) next
      groups[u] <- candidate
      counts$unit_count[old] <- counts$unit_count[old] - 1L
      counts$unit_count[candidate] <- counts$unit_count[candidate] + 1L
      counts$total_count[old] <- counts$total_count[old] - contrib$unit_total[u]
      counts$total_count[candidate] <- counts$total_count[candidate] + contrib$unit_total[u]
      counts$class_count[old, ] <- counts$class_count[old, ] - contrib$unit_class[u, ]
      counts$class_count[candidate, ] <- counts$class_count[candidate, ] + contrib$unit_class[u, ]
      changes <- changes + 1L
      break
    }
  }
  list(groups = groups, changes = changes, row_ll = row_ll)
}

fit_partition_core <- function(data, model, penalty, G, groups, fixed_levels,
                                alpha_info, alpha_fixed, lambda, control,
                                initial_cluster, seed) {
  units <- sort(unique(data$unit_id))
  unit_num <- match(data$unit_id, units)
  rows <- split(seq_len(nrow(data$x)), factor(unit_num, levels = seq_along(units)))
  contrib <- precompute_unit_contributions(rows, data$y, levels(data$y))
  min_class <- resolve_min_class(control$min_per_class, levels(data$y))
  w <- NULL
  bandwidth <- NULL
  if (model == "scmr") {
    coords <- get_unit_coords(data$unit_id, data$coords)
    wobj <- build_unit_weights(coords, control$k_neighbors, control$weight_type,
                               control$weight_bandwidth, TRUE, control$row_standardize_weights)
    w <- wobj$W
    bandwidth <- wobj$bandwidth_used
    unit_groups <- initialize_feasible_groups(coords, G, contrib$unit_total, contrib$unit_class,
      control$min_units, min_class, control$init_tries, scmr_seed(seed, 1000 + G),
      control$init_method, initial_cluster)
    groups <- unit_groups[unit_num]
  } else {
    unit_groups <- vapply(rows, function(idx) as.integer(groups[idx[1]]), integer(1))
    if (!is_feasible_partition(unit_groups, G, contrib$unit_total, contrib$unit_class,
                                control$min_units, min_class)) {
      stop("Fixed clusters do not satisfy min_units/min_per_class.", call. = FALSE)
    }
  }
  tuning <- list()
  tune_seconds <- 0
  alpha <- alpha_info$alpha_default
  tune_lambda <- penalty != "none" && is.null(lambda)
  allow_alpha <- tune_lambda && !alpha_fixed && penalty == "elastic_net" &&
    control$tune_alpha_after_convergence
  if (penalty == "none") {
    lambda <- rep(NA_real_, G)
  } else if (!is.null(lambda)) {
    lambda <- rep(lambda, length.out = G)
  } else {
    t0 <- proc.time()[["elapsed"]]
    init <- tune_initial_partition(data$x, data$y, groups, alpha_info$alpha_grid, control, seed)
    alpha <- init$alpha
    lambda <- init$lambda
    tuning[[length(tuning) + 1L]] <- init$trace
    tune_seconds <- tune_seconds + proc.time()[["elapsed"]] - t0
  }
  fit_groups <- function(groups, alpha, lambda) {
    lapply(seq_len(G), function(g) {
      idx <- which(groups == g)
      fit_local_engine(data$x[idx, , drop = FALSE], data$y[idx], penalty, alpha, lambda[g], control)
    })
  }
  iterations <- list()
  membership_converged <- model == "fixed_clusters"
  reason <- if (membership_converged) "fixed_memberships" else "max_iter_reached"
  tiny <- FALSE
  if (model == "scmr") for (iter in seq_len(control$max_iter)) {
    fits <- fit_groups(groups, alpha, lambda)
    old <- unit_groups
    sweep <- membership_sweep(fits, data, unit_groups, units, w, control, scmr_seed(seed, 100000 + iter))
    unit_groups <- sweep$groups
    tiny <- sweep$changes > 0 && (sweep$changes <= control$tiny_movement_max_units ||
                                   sweep$changes / length(units) <= control$tiny_movement_rate_tol)
    reverted <- tiny && control$tiny_movement_revert
    if (reverted) unit_groups <- old
    groups <- unit_groups[unit_num]
    ll <- sum(sweep$row_ll[cbind(seq_len(nrow(data$x)), groups)])
    bonus <- spatial_bonus(w, unit_groups)
    iterations[[iter]] <- data.frame(Iter = iter, LogLik = ll, SpatialBonus = bonus,
      Objective = ll + control$phi * bonus, Changes = sweep$changes,
      AcceptedChanges = if (reverted) 0L else sweep$changes, Reverted = reverted)
    if (control$verbose) message("SCMR G=", G, " iteration=", iter, " changes=", sweep$changes)
    if (sweep$changes == 0L || tiny) {
      membership_converged <- TRUE
      reason <- if (tiny) "tiny_membership_change" else "zero_membership_change"
      break
    }
    if (tune_lambda && iter < control$max_iter) {
      t0 <- proc.time()[["elapsed"]]
      update <- tune_local_partition(data$x, data$y, groups, alpha, lambda, control, seed, iter)
      lambda <- update$lambda
      tuning[[length(tuning) + 1L]] <- update$trace
      tune_seconds <- tune_seconds + proc.time()[["elapsed"]] - t0
    }
  }
  alpha_changed <- FALSE
  if (allow_alpha && membership_converged && !tiny) {
    t0 <- proc.time()[["elapsed"]]
    update <- tune_final_partition(data$x, data$y, groups, alpha, lambda, control, seed)
    alpha <- update$alpha
    lambda <- update$lambda
    alpha_changed <- update$changed
    tuning[[length(tuning) + 1L]] <- update$trace
    tune_seconds <- tune_seconds + proc.time()[["elapsed"]] - t0
  }
  fits <- fit_groups(groups, alpha, lambda)
  if (model == "scmr" && alpha_changed) {
    check <- membership_sweep(fits, data, unit_groups, units, w, control, scmr_seed(seed, 200000))
    if (check$changes > 0L) {
      membership_converged <- FALSE
      reason <- "membership_not_stationary_after_final_tuning"
    }
  }
  ans <- new_scmr_fit(model, penalty, fits, data, groups, control, seed,
                      bind_diagnostics(tuning), bind_diagnostics(iterations),
                      membership_converged, reason, w, bandwidth, fixed_levels)
  ans$lambda_rule <- if (penalty == "none") "none" else if (tune_lambda) "holdout_tolerance" else "fixed"
  ans$runtime_tuning <- tune_seconds
  ans
}

#' Fit spatially clustered multinomial regression
#' @inheritParams fit_scmr
#' @return A scmr_spatial object inheriting from scmr_fit.
#' @export
fit_spatial_scmr <- function(x, y, unit_id, coords, G,
                             penalty = c("elastic_net", "lasso", "ridge", "none"),
                             alpha = NULL, alpha_grid = NULL, lambda = NULL,
                             control = scmr_control(), initial_cluster = NULL, seed = 123) {
  fit_scmr(x, y, model = "scmr", penalty = match.arg(penalty), G = G, unit_id = unit_id,
           coords = coords, alpha = alpha, alpha_grid = alpha_grid, lambda = lambda,
           control = control, initial_cluster = initial_cluster, seed = seed)
}
