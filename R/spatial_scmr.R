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

# Penalty term sum_g s_g * P(B_g). In sum-scale mode s_g is fixed (glmnet lambda
# times the cluster size at fitting time); in mean-scale mode s_g = n_g * lambda_g
# changes with the current cluster size n_g.
partition_penalty <- function(fits, groups, sum_scale = FALSE) {
  sum(vapply(seq_along(fits), function(g) {
    if (fits[[g]]$kind != "glmnet") return(0)
    size <- if (sum_scale) fits[[g]]$n_fit else sum(groups == g)
    size * fits[[g]]$lambda * engine_penalty(fits[[g]])
  }, numeric(1)))
}

fit_partition_core <- function(data, model, penalty, G, groups, fixed_levels,
                                alpha_info, alpha_fixed, lambda, control,
                                initial_cluster, seed, init_features = NULL) {
  units <- sort(unique(data$unit_id))
  unit_num <- match(data$unit_id, units)
  rows <- split(seq_len(nrow(data$x)), factor(unit_num, levels = seq_along(units)))
  contrib <- precompute_unit_contributions(rows, data$y, levels(data$y))
  min_class <- resolve_min_class(control$min_per_class, levels(data$y))
  # Sum-scale mode: fixed lambda_g * n_g, global standardization, no inner tuning.
  sum_scale <- identical(control$lambda_scale, "sum") && penalty != "none"
  pre <- if (sum_scale && control$standardize) make_pre_transform(data$x) else NULL
  xt <- apply_pre_transform(data$x, pre)
  tctrl <- control
  if (!is.null(pre)) tctrl$standardize <- FALSE
  w <- NULL
  bandwidth <- NULL
  coords <- NULL
  if (model == "scmr") {
    coords <- get_unit_coords(data$unit_id, data$coords)
    wobj <- build_unit_weights(coords, control$k_neighbors, control$weight_type,
                               control$weight_bandwidth, TRUE, control$row_standardize_weights)
    w <- wobj$W
    bandwidth <- wobj$bandwidth_used
  } else {
    fixed_unit_groups <- vapply(rows, function(idx) as.integer(groups[idx[1]]), integer(1))
    if (!is_feasible_partition(fixed_unit_groups, G, contrib$unit_total, contrib$unit_class,
                                control$min_units, min_class)) {
      stop("Fixed clusters do not satisfy min_units/min_per_class.", call. = FALSE)
    }
  }
  tune_lambda <- penalty != "none" && is.null(lambda)
  allow_alpha <- tune_lambda && !alpha_fixed && penalty == "elastic_net" &&
    control$tune_alpha_after_convergence && !sum_scale
  local_retune <- tune_lambda && !sum_scale
  update <- model == "scmr" && control$update_memberships

  fit_groups <- function(groups, alpha, lambda) {
    lapply(seq_len(G), function(g) {
      idx <- which(groups == g)
      lam <- if (penalty == "none") NA_real_ else if (sum_scale) lambda[g] / length(idx) else lambda[g]
      engine <- fit_local_engine(data$x[idx, , drop = FALSE], data$y[idx], penalty, alpha, lam, control, pre)
      engine$n_fit <- length(idx)
      engine
    })
  }
  objective_of <- function(fits, groups, unit_groups) {
    ll <- sum(vapply(seq_len(G), function(g) {
      idx <- which(groups == g)
      sum(true_class_logp(engine_predict(fits[[g]], data$x[idx, , drop = FALSE]),
                          data$y[idx], levels(data$y)))
    }, numeric(1)))
    ll + control$phi * spatial_bonus(w, unit_groups) - partition_penalty(fits, groups, sum_scale)
  }

  run_start <- function(unit_groups, start_id) {
    groups <- unit_groups[unit_num]
    tuning <- list()
    tune_seconds <- 0
    alpha <- alpha_info$alpha_default
    lam <- if (penalty == "none") rep(NA_real_, G) else if (!is.null(lambda)) rep(lambda, length.out = G) else NULL
    if (tune_lambda) {
      t0 <- proc.time()[["elapsed"]]
      init <- tune_initial_partition(xt, data$y, groups, alpha_info$alpha_grid, tctrl, seed)
      alpha <- init$alpha
      lam <- init$lambda
      tuning[[length(tuning) + 1L]] <- init$trace
      tune_seconds <- tune_seconds + proc.time()[["elapsed"]] - t0
    }
    # Convert the mean-scale lambda chosen on the initial partition to sum scale.
    if (sum_scale) lam <- lam * tabulate(groups, G)
    iterations <- list()
    membership_converged <- !update
    reason <- if (model == "fixed_clusters") "fixed_memberships" else
      if (!update) "two_stage_initial_partition" else "max_iter_reached"
    tiny <- FALSE
    if (update) for (iter in seq_len(control$max_iter)) {
      fits <- fit_groups(groups, alpha, lam)
      old <- unit_groups
      sweep <- membership_sweep(fits, data, unit_groups, units, w, control,
                                scmr_seed(seed, 100000 + iter + 7919L * (start_id - 1L)))
      unit_groups <- sweep$groups
      tiny <- sweep$changes > 0 && (sweep$changes <= control$tiny_movement_max_units ||
                                     sweep$changes / length(units) <= control$tiny_movement_rate_tol)
      reverted <- tiny && control$tiny_movement_revert
      if (reverted) unit_groups <- old
      groups <- unit_groups[unit_num]
      ll <- sum(sweep$row_ll[cbind(seq_len(nrow(data$x)), groups)])
      bonus <- spatial_bonus(w, unit_groups)
      # Q evaluated at (current coefficients, updated memberships).
      iterations[[iter]] <- data.frame(Start = start_id, Iter = iter, LogLik = ll, SpatialBonus = bonus,
        Objective = ll + control$phi * bonus,
        PenalizedObjective = ll + control$phi * bonus - partition_penalty(fits, groups, sum_scale),
        Changes = sweep$changes, AcceptedChanges = if (reverted) 0L else sweep$changes,
        Reverted = reverted)
      if (control$verbose) message("SCMR G=", G, " start=", start_id, " iteration=", iter,
                                   " changes=", sweep$changes)
      if (sweep$changes == 0L || tiny) {
        membership_converged <- TRUE
        reason <- if (tiny) "tiny_membership_change" else "zero_membership_change"
        break
      }
      if (local_retune && iter < control$max_iter) {
        t0 <- proc.time()[["elapsed"]]
        upd <- tune_local_partition(xt, data$y, groups, alpha, lam, tctrl, seed, iter)
        lam <- upd$lambda
        tuning[[length(tuning) + 1L]] <- upd$trace
        tune_seconds <- tune_seconds + proc.time()[["elapsed"]] - t0
      }
    }
    alpha_changed <- FALSE
    if (allow_alpha && membership_converged && !tiny) {
      t0 <- proc.time()[["elapsed"]]
      upd <- tune_final_partition(xt, data$y, groups, alpha, lam, tctrl, seed)
      alpha <- upd$alpha
      lam <- upd$lambda
      alpha_changed <- upd$changed
      tuning[[length(tuning) + 1L]] <- upd$trace
      tune_seconds <- tune_seconds + proc.time()[["elapsed"]] - t0
    }
    fits <- fit_groups(groups, alpha, lam)
    if (update && alpha_changed) {
      check <- membership_sweep(fits, data, unit_groups, units, w, control, scmr_seed(seed, 200000))
      if (check$changes > 0L) {
        membership_converged <- FALSE
        reason <- "membership_not_stationary_after_final_tuning"
      }
    }
    list(unit_groups = unit_groups, groups = groups, fits = fits, tuning = tuning,
         iterations = iterations, membership_converged = membership_converged,
         reason = reason, tune_seconds = tune_seconds,
         objective = objective_of(fits, groups, unit_groups))
  }

  starts <- data.frame()
  if (model == "scmr") {
    methods <- control$start_methods %||%
      c(control$init_method, rep(c("coefficient", "kmeans"), length.out = max(0L, control$n_starts - 1L)))
    methods <- rep(methods, length.out = control$n_starts)
    if (!is.null(initial_cluster)) methods[1] <- "supplied"
    if (any(methods == "coefficient")) {
      if (is.null(init_features)) {
        init_features <- scmr_local_coefficients(data$x, data$y, data$unit_id, data$coords,
          k = control$coef_init_k, alpha = control$coef_init_alpha, lambda = control$coef_init_lambda,
          anchors = control$coef_init_anchors, type_multinomial = control$type_multinomial,
          seed = scmr_seed(seed, 4242))
      }
      init_features <- as.matrix(init_features)
      if (!is.null(rownames(init_features))) init_features <- init_features[units, , drop = FALSE]
      if (nrow(init_features) != length(units)) stop("init_features needs one row per spatial unit.", call. = FALSE)
    }
    results <- list()
    for (s in seq_along(methods)) {
      ug <- tryCatch(initialize_feasible_groups(coords, G, contrib$unit_total, contrib$unit_class,
          control$min_units, min_class, control$init_tries, scmr_seed(seed, 1000 + G + 31L * (s - 1L)),
          if (methods[s] == "supplied") control$init_method else methods[s],
          if (methods[s] == "supplied") initial_cluster else NULL, init_features),
        error = function(e) e)
      if (inherits(ug, "error")) {
        if (s == 1L && length(methods) == 1L) stop(conditionMessage(ug), call. = FALSE)
        starts <- rbind(starts, data.frame(Start = s, Method = methods[s], Objective = NA_real_,
          MembershipConverged = NA, Reason = paste("initialization_failed:", conditionMessage(ug))))
        next
      }
      res <- run_start(ug, s)
      res$start <- s
      results[[length(results) + 1L]] <- res
      starts <- rbind(starts, data.frame(Start = s, Method = methods[s], Objective = res$objective,
        MembershipConverged = res$membership_converged, Reason = res$reason))
    }
    if (!length(results)) stop("Every initialization failed.", call. = FALSE)
    best <- results[[which.max(vapply(results, function(z) z$objective, numeric(1)))]]
    starts$Selected <- starts$Start == best$start
  } else {
    best <- run_start(fixed_unit_groups, 1L)
  }
  ans <- new_scmr_fit(model, penalty, best$fits, data, best$groups, control, seed,
                      bind_diagnostics(best$tuning), bind_diagnostics(best$iterations),
                      best$membership_converged, best$reason, w, bandwidth, fixed_levels)
  ans$lambda_rule <- if (penalty == "none") "none" else if (tune_lambda) "holdout_tolerance" else "fixed"
  ans$lambda_scale <- if (sum_scale) "sum" else "mean"
  ans$lambda_sum <- if (penalty == "none") rep(NA_real_, G) else ans$lambda_used * tabulate(best$groups, G)
  ans$runtime_tuning <- sum(vapply(if (model == "scmr") results else list(best),
                                   function(z) z$tune_seconds, numeric(1)))
  ans$diagnostics$starts <- starts
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
