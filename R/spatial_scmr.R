#' Fit spatially clustered multinomial regression
#'
#' This function fits the spatial SCMR model. Rows sharing the same `unit_id`
#' are assigned to the same spatial cluster. Each cluster has a local penalized
#' multinomial regression model, and membership updates combine local likelihood
#' with a Potts-style spatial smoothness term.
#'
#' @param x Numeric matrix or sparse matrix of predictors.
#' @param y Factor or character vector of class labels.
#' @param unit_id Spatial unit identifier. Rows sharing a unit move together.
#' @param coords Two-column coordinate matrix for rows in `x`.
#' @param G Number of spatial clusters.
#' @param penalty Penalty type: `"ridge"`, `"lasso"`, or `"elastic_net"`.
#' @param alpha Optional fixed alpha value.
#' @param alpha_grid Optional candidate alpha values for elastic-net tuning.
#' @param lambda Optional fixed lambda value. Use a scalar or a vector of length `G`.
#' @param control List from [scmr_control()].
#' @param initial_cluster Optional initial cluster labels for unique spatial units.
#' @param seed Random seed.
#' @return An object of class `scmr_spatial`.
#' @export
fit_spatial_scmr <- function(x, y, unit_id, coords, G,
                             penalty = c("elastic_net", "lasso", "ridge"),
                             alpha = NULL,
                             alpha_grid = NULL,
                             lambda = NULL,
                             control = scmr_control(),
                             initial_cluster = NULL,
                             seed = 123) {
  penalty <- match.arg(penalty)
  if (!is.null(seed)) set.seed(seed)

  alpha_info <- resolve_penalty_alpha(penalty, alpha = alpha, alpha_grid = alpha_grid %||% control$alpha_grid)
  control$alpha_grid <- alpha_info$alpha_grid
  if (penalty %in% c("ridge", "lasso") || !is.null(alpha)) {
    control$tune_alpha_after_convergence <- FALSE
  }
  if (!is.null(lambda)) {
    control$tune_lambda <- FALSE
    control$tune_alpha_after_convergence <- FALSE
  }

  y <- factor(y)
  class_levels <- levels(y)
  K <- length(class_levels)
  unit_factor <- factor(unit_id)
  unit_levels <- levels(unit_factor)
  unit_row_index <- split(seq_len(nrow(x)), unit_factor)
  unit_num <- as.integer(unit_factor)
  m <- length(unit_levels)
  unit_coords <- get_unit_coords(unit_id, coords)
  contrib <- precompute_unit_contributions(unit_row_index, y, class_levels)
  unit_total <- contrib$unit_total
  unit_class <- contrib$unit_class
  min_class_vec <- resolve_min_class(control$min_per_class, class_levels)

  wobj <- build_unit_weights(
    unit_coords,
    k = control$k_neighbors,
    weight_type = control$weight_type,
    bandwidth = control$weight_bandwidth,
    symmetrize = TRUE,
    row_standardize = control$row_standardize_weights
  )
  W <- wobj$W

  group_unit <- initialize_feasible_groups(
    unit_coords, G, unit_total, unit_class,
    min_units = control$min_units,
    min_class_vec = min_class_vec,
    tries = control$init_tries,
    seed = seed,
    method = control$init_method,
    initial_cluster = initial_cluster
  )

  counts <- compute_cluster_counts(group_unit, G, unit_total, unit_class)
  unit_count <- counts$unit_count
  total_count <- counts$total_count
  class_count <- counts$class_count

  alpha_current <- if (!is.null(alpha)) as.numeric(alpha)[1] else alpha_info$alpha_default
  lambda_vec <- rep(NA_real_, G)
  tuning_initial <- list()

  if (!is.null(lambda)) {
    lambda_vec <- rep(as.numeric(lambda), length.out = G)
  } else {
    row_group_init <- group_unit[unit_num]
    for (g in seq_len(G)) {
      idx_g <- which(row_group_init == g)
      tune_g <- fast_holdout_tune_alpha_lambda(
        x[idx_g, , drop = FALSE], y[idx_g], control$alpha_grid,
        class_levels = class_levels, control = control,
        seed = seed + 8300 + G + 17 * g
      )
      lambda_vec[g] <- tune_g$lambda
      if (isTRUE(control$tune_alpha_after_convergence) && !is.null(tune_g$trace) && nrow(tune_g$trace) > 0) {
        tune_g$trace$Cluster <- g
        tuning_initial[[length(tuning_initial) + 1L]] <- tune_g$trace
      }
    }

    all_init <- if (length(tuning_initial) > 0) do.call(rbind, tuning_initial) else data.frame()
    if (nrow(all_init) > 0 && is.null(alpha)) {
      agg <- stats::aggregate(HoldoutSelectedLoss ~ Alpha, all_init, mean, na.rm = TRUE)
      agg <- agg[order(agg$HoldoutSelectedLoss), , drop = FALSE]
      alpha_current <- as.numeric(agg$Alpha[1])
      for (g in seq_len(G)) {
        rr_g <- all_init[
          all_init$Cluster == g &
            abs(all_init$Alpha - alpha_current) < 1e-12 &
            is.finite(all_init$LambdaSelected), , drop = FALSE]
        if (nrow(rr_g) > 0) {
          lambda_vec[g] <- rr_g$LambdaSelected[which.min(rr_g$HoldoutSelectedLoss)]
        }
      }
    }
  }

  if (!is.finite(alpha_current)) alpha_current <- alpha_info$alpha_default
  if (any(!is.finite(lambda_vec) | lambda_vec <= 0)) {
    stop("All SCMR cluster lambda values must be positive and finite.", call. = FALSE)
  }

  make_row_loglik_from_fits <- function(fits_current) {
    row_loglik_by_g <- matrix(NA_real_, nrow = nrow(x), ncol = G)
    for (gg in seq_len(G)) {
      lambda_gg <- get_s_value(fits_current[[gg]], NULL)
      p_all_gg <- predict_multinom_prob(fits_current[[gg]], newx = x, s = lambda_gg)
      p_all_gg <- normalize_prob_matrix(p_all_gg, class_levels, nrow(x))
      row_loglik_by_g[, gg] <- true_class_logp(p_all_gg, y, class_levels)
    }
    row_loglik_by_g
  }

  update_membership_once <- function(row_loglik_by_g) {
    unit_loglik_by_g <- matrix(0, nrow = m, ncol = G)
    for (u in seq_len(m)) {
      idx_u <- unit_row_index[[u]]
      unit_loglik_by_g[u, ] <- colSums(row_loglik_by_g[idx_u, , drop = FALSE])
    }
    order_u <- seq_len(m)
    if (identical(control$update_order, "random")) order_u <- sample(order_u)
    changes <- 0L
    if (G > 1) {
      for (u in order_u) {
        g_old <- group_unit[u]
        membership_mat <- stats::model.matrix(~ factor(group_unit, levels = seq_len(G)) - 1)
        neigh_score <- as.numeric(W[u, ] %*% membership_mat)
        score <- unit_loglik_by_g[u, ] + control$phi * neigh_score
        candidate_order <- order(score, decreasing = TRUE)
        chosen <- g_old
        for (g_new in candidate_order) {
          if (move_is_feasible(u, g_new, group_unit, unit_count, class_count,
                               unit_class, control$min_units, min_class_vec)) {
            chosen <- g_new
            break
          }
        }
        if (chosen != g_old) {
          group_unit[u] <<- chosen
          unit_count[g_old] <<- unit_count[g_old] - 1L
          unit_count[chosen] <<- unit_count[chosen] + 1L
          total_count[g_old] <<- total_count[g_old] - unit_total[u]
          total_count[chosen] <<- total_count[chosen] + unit_total[u]
          class_count[g_old, ] <<- class_count[g_old, ] - unit_class[u, ]
          class_count[chosen, ] <<- class_count[chosen, ] + unit_class[u, ]
          changes <- changes + 1L
        }
      }
    }
    changes
  }

  obj_trace <- data.frame()
  converged <- FALSE
  convergence_reason <- "max_iter_reached"
  skip_final_alpha <- FALSE

  for (iter in seq_len(control$max_iter)) {
    row_group <- group_unit[unit_num]
    fits <- vector("list", G)
    for (g in seq_len(G)) {
      idx_g <- which(row_group == g)
      fits[[g]] <- fit_fixed_lambda_glmnet(
        x[idx_g, , drop = FALSE], y[idx_g], alpha_current, lambda_vec[g],
        class_levels = class_levels, standardize = control$standardize,
        type_multinomial = control$type_multinomial
      )
    }

    row_loglik_by_g <- make_row_loglik_from_fits(fits)
    old_group_unit <- group_unit
    changes <- update_membership_once(row_loglik_by_g)
    row_group_new <- group_unit[unit_num]
    current_loglik <- sum(row_loglik_by_g[cbind(seq_len(nrow(x)), row_group_new)])
    Wsum <- summary(W)
    spatial_bonus <- if (G == 1) {
      0.5 * sum(Wsum$x)
    } else {
      0.5 * sum(Wsum$x[group_unit[Wsum$i] == group_unit[Wsum$j]])
    }
    obj_trace <- rbind(obj_trace, data.frame(
      Iter = iter,
      LogLik = current_loglik,
      SpatialBonus = spatial_bonus,
      Objective = current_loglik + control$phi * spatial_bonus,
      Changes = changes
    ))

    if (changes == 0) {
      converged <- TRUE
      convergence_reason <- "zero_membership_change"
      break
    }

    is_tiny <- changes <= control$tiny_movement_max_units ||
      (changes / m) <= control$tiny_movement_rate_tol
    if (is_tiny) {
      if (isTRUE(control$tiny_movement_revert)) {
        group_unit <- old_group_unit
        counts_revert <- compute_cluster_counts(group_unit, G, unit_total, unit_class)
        unit_count <- counts_revert$unit_count
        total_count <- counts_revert$total_count
        class_count <- counts_revert$class_count
      }
      converged <- TRUE
      convergence_reason <- "tiny_membership_change_skip_retune"
      skip_final_alpha <- TRUE
      break
    }

    if (isTRUE(control$tune_lambda)) {
      row_group_after <- group_unit[unit_num]
      for (g in seq_len(G)) {
        idx_g <- which(row_group_after == g)
        upd <- local_lambda_update(
          x[idx_g, , drop = FALSE], y[idx_g], alpha_current, lambda_vec[g],
          class_levels = class_levels, control = control,
          seed = seed + 9100 + 1000 * G + 17 * g + iter
        )
        lambda_vec[g] <- upd$lambda
      }
    }
  }

  if (isTRUE(control$tune_lambda) && isTRUE(control$tune_alpha_after_convergence) && !isTRUE(skip_final_alpha)) {
    alpha_update <- final_tune_alpha_after_convergence(
      x, y, group_unit[unit_num], alpha_current, lambda_vec,
      class_levels, control, seed = seed + 12000 + 1000 * G
    )
    alpha_current <- alpha_update$alpha
    lambda_vec <- alpha_update$lambda_vec
  }

  final_row_group <- group_unit[unit_num]
  final_fits <- vector("list", G)
  for (g in seq_len(G)) {
    idx_g <- which(final_row_group == g)
    final_fits[[g]] <- fit_fixed_lambda_glmnet(
      x[idx_g, , drop = FALSE], y[idx_g], alpha_current, lambda_vec[g],
      class_levels = class_levels, standardize = control$standardize,
      type_multinomial = control$type_multinomial
    )
  }

  prob_train <- matrix(NA_real_, nrow = nrow(x), ncol = K)
  colnames(prob_train) <- class_levels
  row_loglik_by_g_final <- matrix(NA_real_, nrow = nrow(x), ncol = G)
  active_per_group <- integer(G)
  df_effective_group <- numeric(G)

  for (g in seq_len(G)) {
    idx_g <- which(final_row_group == g)
    lambda_g <- get_s_value(final_fits[[g]], NULL)
    active_per_group[g] <- active_predictor_count(final_fits[[g]], s = lambda_g)
    df_effective_group[g] <- elastic_net_effective_df(
      final_fits[[g]], x[idx_g, , drop = FALSE], s = lambda_g,
      alpha = alpha_current, class_levels = class_levels,
      standardize = control$standardize
    )
    p_all_g <- predict_multinom_prob(final_fits[[g]], newx = x, s = lambda_g)
    p_all_g <- normalize_prob_matrix(p_all_g, class_levels, nrow(x))
    row_loglik_by_g_final[, g] <- true_class_logp(p_all_g, y, class_levels)
    prob_train[idx_g, ] <- p_all_g[idx_g, , drop = FALSE]
  }

  criteria <- compute_criteria_from_row_loglik(
    row_loglik_by_g = row_loglik_by_g_final,
    group_unit = group_unit,
    unit_num = unit_num,
    unit_id = unit_id,
    y = y,
    x_ncol = ncol(x),
    W = W,
    active_per_group = active_per_group,
    G = G,
    lambda_used = lambda_vec,
    class_levels = class_levels,
    unit_total = unit_total,
    unit_class = unit_class,
    df_effective_group = df_effective_group,
    alpha = alpha_current,
    phi = control$phi,
    criterion_balance_gamma = control$criterion_balance_gamma
  )

  PP_unit <- matrix(0, nrow = m, ncol = G)
  PP_unit[cbind(seq_len(m), group_unit)] <- 1
  colnames(PP_unit) <- paste0("G", seq_len(G))
  rownames(PP_unit) <- unit_levels

  out <- list(
    model = "scmr",
    penalty = penalty,
    G = G,
    group_unit = stats::setNames(group_unit, unit_levels),
    PP_unit = PP_unit,
    unit_levels = unit_levels,
    train_unit_coords = unit_coords,
    class_levels = class_levels,
    fits = final_fits,
    prob_train = normalize_prob_matrix(prob_train, class_levels, nrow(x)),
    row_loglik_by_g = row_loglik_by_g_final,
    y_train = y,
    row_group = final_row_group,
    active_per_group = active_per_group,
    lambda_used = lambda_vec,
    df_effective_group = df_effective_group,
    W = W,
    control = control,
    obj_trace = obj_trace,
    converged = converged,
    convergence_reason = convergence_reason,
    criteria = criteria,
    x_colnames = colnames(x),
    alpha = alpha_current,
    alpha_final = alpha_current,
    initial_tuning = if (exists("all_init")) all_init else data.frame()
  )
  class(out) <- "scmr_spatial"
  out
}

assign_test_PP <- function(fit, new_unit_id, new_coords) {
  new_unit_id <- as.character(new_unit_id)
  train_units <- rownames(fit$PP_unit)
  train_coords <- fit$train_unit_coords
  train_PP <- fit$PP_unit
  new_coord_df <- unique(data.frame(
    unit = new_unit_id,
    xcoord = new_coords[, 1],
    ycoord = new_coords[, 2],
    stringsAsFactors = FALSE
  ))
  row_unit_factor <- factor(new_unit_id, levels = new_coord_df$unit)
  new_units <- new_coord_df$unit
  G <- ncol(train_PP)
  PP_new_unit <- matrix(0, nrow = length(new_units), ncol = G)
  colnames(PP_new_unit) <- colnames(train_PP)
  rownames(PP_new_unit) <- new_units

  for (i in seq_along(new_units)) {
    u <- new_units[i]
    if (u %in% train_units) {
      PP_new_unit[i, ] <- train_PP[u, ]
    } else {
      coord_i <- as.numeric(new_coord_df[i, c("xcoord", "ycoord")])
      d <- sqrt(rowSums((sweep(train_coords, 2, coord_i, "-"))^2))
      kk <- min(fit$control$k_neighbors, length(d))
      nn <- order(d)[seq_len(kk)]
      if (fit$control$weight_type == "binary") {
        ww <- rep(1, kk)
      } else {
        bw <- fit$control$weight_bandwidth
        if (is.null(bw) || !is.finite(bw) || bw <= 0) {
          bw <- stats::median(d[nn], na.rm = TRUE)
        }
        ww <- exp(-(d[nn]^2) / (bw^2))
      }
      ww <- ww / sum(ww)
      PP_new_unit[i, ] <- colSums(train_PP[nn, , drop = FALSE] * ww)
    }
  }
  PP_new_row <- PP_new_unit[as.character(row_unit_factor), , drop = FALSE]
  rownames(PP_new_row) <- NULL
  normalize_prob_matrix(PP_new_row, colnames(train_PP), length(new_unit_id))
}

#' @export
predict.scmr_spatial <- function(object, newx, new_unit_id, new_coords,
                                 type = c("prob", "class", "membership"), ...) {
  type <- match.arg(type)
  PP_new <- assign_test_PP(object, new_unit_id, new_coords)
  if (type == "membership") return(PP_new)

  G <- length(object$fits)
  classes <- object$class_levels
  prob <- matrix(0, nrow = nrow(newx), ncol = length(classes))
  colnames(prob) <- classes
  for (g in seq_len(G)) {
    p_g <- predict_multinom_prob(object$fits[[g]], newx = newx, s = object$lambda_used[g])
    p_g <- normalize_prob_matrix(p_g, classes, nrow(newx))
    prob <- prob + p_g * PP_new[, g]
  }
  prob <- normalize_prob_matrix(prob, classes, nrow(newx))
  if (type == "prob") return(prob)
  factor(classes[max.col(prob)], levels = classes)
}

#' @export
print.scmr_spatial <- function(x, ...) {
  cat("Spatially clustered multinomial regression (SCMR)\n")
  cat("  Penalty:", x$penalty, "\n")
  cat("  Clusters:", x$G, "\n")
  cat("  Classes:", paste(x$class_levels, collapse = ", "), "\n")
  cat("  Alpha:", x$alpha_final, "\n")
  cat("  Converged:", x$converged, " (", x$convergence_reason, ")\n", sep = "")
  invisible(x)
}
