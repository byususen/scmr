make_holdout_validation_split <- function(y, val_prop = 0.30, seed = 123) {
  set.seed(seed)
  y <- factor(y)
  idx_all <- seq_along(y)
  split_idx <- split(idx_all, y)
  val_idx <- unlist(lapply(split_idx, function(idx) {
    n <- length(idx)
    if (n <= 3) return(idx[1])
    n_val <- max(1, floor(val_prop * n))
    n_val <- min(n_val, n - 2)
    sample(idx, size = n_val)
  }))
  val_idx <- sort(unique(val_idx))
  list(train_idx = setdiff(idx_all, val_idx), val_idx = val_idx)
}

path_logloss_by_lambda <- function(fit, x_val, y_val, lambda_values, class_levels) {
  arr <- stats::predict(fit, newx = x_val, s = lambda_values, type = "response")
  out <- rep(Inf, length(lambda_values))
  if (length(dim(arr)) == 3) {
    for (li in seq_along(lambda_values)) {
      prob <- arr[, , li, drop = FALSE]
      prob <- prob[, , 1]
      if (is.null(colnames(prob)) && !is.null(dimnames(arr)[[2]])) colnames(prob) <- dimnames(arr)[[2]]
      prob <- normalize_prob_matrix(prob, class_levels, length(y_val))
      out[li] <- multiclass_logloss(prob, factor(y_val, levels = class_levels))
    }
  } else {
    prob <- normalize_prob_matrix(as.matrix(arr), class_levels, length(y_val))
    out[1] <- multiclass_logloss(prob, factor(y_val, levels = class_levels))
  }
  out
}

select_lambda_1se_like <- function(lambda_values, losses, tol = 0.002) {
  lambda_values <- as.numeric(lambda_values)
  losses <- as.numeric(losses)
  ok <- is.finite(lambda_values) & lambda_values > 0 & is.finite(losses)
  if (!any(ok)) return(list(lambda = NA_real_, loss = Inf, min_loss = Inf))
  idx_ok <- which(ok)
  min_pos <- idx_ok[which.min(losses[idx_ok])]
  min_loss <- losses[min_pos]
  eligible <- idx_ok[losses[idx_ok] <= min_loss + tol]
  sel_pos <- eligible[which.max(lambda_values[eligible])]
  list(lambda = lambda_values[sel_pos], loss = losses[sel_pos], min_loss = min_loss)
}

fast_holdout_tune_alpha_lambda <- function(x, y, alpha_grid, class_levels, control,
                                           seed = 123) {
  y <- factor(y, levels = class_levels)
  split <- make_holdout_validation_split(y, val_prop = control$holdout_val_prop, seed = seed)
  tr_idx <- split$train_idx
  val_idx <- split$val_idx
  rows <- list()
  best <- list(alpha = NA_real_, lambda = NA_real_, loss = Inf, min_loss = Inf)
  for (aa in alpha_grid) {
    fit <- tryCatch({
      glmnet::glmnet(
        x = x[tr_idx, , drop = FALSE],
        y = factor(y[tr_idx], levels = class_levels),
        family = "multinomial", alpha = aa, nlambda = control$holdout_nlambda,
        standardize = control$standardize, type.multinomial = control$type_multinomial,
        maxit = 100000
      )
    }, error = function(e) e)
    if (inherits(fit, "error")) next
    lambda_values <- as.numeric(fit$lambda)
    lambda_values <- lambda_values[is.finite(lambda_values) & lambda_values > 0]
    losses <- tryCatch(
      path_logloss_by_lambda(fit, x[val_idx, , drop = FALSE], y[val_idx], lambda_values, class_levels),
      error = function(e) rep(Inf, length(lambda_values))
    )
    sel <- select_lambda_1se_like(lambda_values, losses, tol = control$lambda_1se_tol)
    rows[[length(rows) + 1L]] <- data.frame(
      Alpha = aa, LambdaSelected = sel$lambda,
      HoldoutSelectedLoss = sel$loss, HoldoutMinLoss = sel$min_loss,
      NumLambda = length(lambda_values), stringsAsFactors = FALSE
    )
    if (is.finite(sel$loss) && sel$loss < best$loss) {
      best <- list(alpha = aa, lambda = sel$lambda, loss = sel$loss, min_loss = sel$min_loss)
    }
  }
  if (!is.finite(best$lambda) || !is.finite(best$alpha)) {
    stop("Holdout alpha/lambda tuning failed.", call. = FALSE)
  }
  best$trace <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()
  best
}

local_lambda_update <- function(x, y, alpha, lambda_old, class_levels, control,
                                seed = 123) {
  y <- factor(y, levels = class_levels)
  lambda_old <- as.numeric(lambda_old)[1]
  if (!is.finite(lambda_old) || lambda_old <= 0) {
    return(list(lambda = lambda_old, loss = Inf, action = "reuse_bad_old_lambda"))
  }
  lambda_grid <- sort(unique(as.numeric(lambda_old * control$local_lambda_factors)), decreasing = TRUE)
  lambda_grid <- lambda_grid[is.finite(lambda_grid) & lambda_grid > 0]
  split <- make_holdout_validation_split(y, val_prop = control$holdout_val_prop, seed = seed)
  tr_idx <- split$train_idx
  val_idx <- split$val_idx
  fit <- tryCatch({
    lambda_path <- sort(unique(as.numeric(c(lambda_grid * 10, lambda_grid, lambda_grid / 10))), decreasing = TRUE)
    lambda_path <- lambda_path[is.finite(lambda_path) & lambda_path > 0]
    glmnet::glmnet(
      x = x[tr_idx, , drop = FALSE], y = factor(y[tr_idx], levels = class_levels),
      family = "multinomial", alpha = alpha, lambda = lambda_path,
      standardize = control$standardize, type.multinomial = control$type_multinomial,
      maxit = 100000
    )
  }, error = function(e) e)
  if (inherits(fit, "error")) return(list(lambda = lambda_old, loss = Inf, action = "reuse_after_fit_error"))
  losses <- tryCatch(
    path_logloss_by_lambda(fit, x[val_idx, , drop = FALSE], y[val_idx], lambda_grid, class_levels),
    error = function(e) rep(Inf, length(lambda_grid))
  )
  sel <- select_lambda_1se_like(lambda_grid, losses, tol = control$lambda_1se_tol)
  old_idx <- which.min(abs(log(lambda_grid) - log(lambda_old)))
  old_loss <- losses[old_idx]
  if (is.finite(sel$loss) && (!is.finite(old_loss) || sel$loss < old_loss - control$lambda_update_tol)) {
    return(list(lambda = sel$lambda, loss = sel$loss, action = "update_lambda"))
  }
  list(lambda = lambda_old, loss = old_loss, action = "reuse_lambda")
}

final_tune_alpha_after_convergence <- function(x, y, row_group, alpha_current,
                                               lambda_vec_current, class_levels,
                                               control, seed = 123) {
  G <- max(row_group)
  alpha_candidates <- unique(pmin(pmax(alpha_current + control$final_alpha_delta, 0.05), 1.0))
  alpha_candidates <- sort(alpha_candidates)
  result <- list()
  for (aa in alpha_candidates) {
    lambda_vec <- numeric(G)
    loss_vec <- rep(Inf, G)
    for (g in seq_len(G)) {
      idx_g <- which(row_group == g)
      upd <- local_lambda_update(
        x[idx_g, , drop = FALSE], y[idx_g], aa, lambda_vec_current[g],
        class_levels, control, seed = seed + 29 * g
      )
      lambda_vec[g] <- upd$lambda
      loss_vec[g] <- upd$loss
    }
    result[[as.character(aa)]] <- list(
      alpha = aa, lambda_vec = lambda_vec, loss_vec = loss_vec,
      agg_loss = mean(loss_vec[is.finite(loss_vec)], na.rm = TRUE)
    )
  }
  agg <- vapply(result, function(z) z$agg_loss, numeric(1))
  if (!any(is.finite(agg))) {
    return(list(alpha = alpha_current, lambda_vec = lambda_vec_current, action = "reuse_after_alpha_tuning_failure"))
  }
  best <- result[[names(agg)[which.min(agg)]]]
  current_loss <- if (as.character(alpha_current) %in% names(result)) result[[as.character(alpha_current)]]$agg_loss else Inf
  if (is.finite(best$agg_loss) && (!is.finite(current_loss) || best$agg_loss < current_loss - control$alpha_update_tol)) {
    return(list(alpha = best$alpha, lambda_vec = best$lambda_vec, action = "update_alpha_after_convergence"))
  }
  list(alpha = alpha_current, lambda_vec = lambda_vec_current, action = "reuse_alpha")
}
