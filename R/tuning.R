# Holdout tuning returns tables; the study adds labels and writes files.
make_holdout_validation_split <- function(y, val_prop = 0.30, seed = 123) {
  with_scmr_seed(seed, {
    strata <- split(seq_along(y), factor(y))
    val <- unlist(lapply(strata, function(idx) {
      if (length(idx) < 3L) return(integer())
      nval <- min(max(1L, floor(val_prop * length(idx))), length(idx) - 2L)
      idx[sample.int(length(idx), nval)]
    }), use.names = FALSE)
    if (!length(val)) stop("Insufficient class support for a holdout split.", call. = FALSE)
    list(train_idx = setdiff(seq_along(y), val), val_idx = sort(val))
  })
}

path_logloss_by_lambda <- function(fit, x_val, y_val, lambda_values, class_levels) {
  arr <- stats::predict(fit, newx = glmnet_design(x_val), s = lambda_values, type = "response")
  if (length(dim(arr)) != 3L) stop("Unexpected multinomial prediction dimensions.", call. = FALSE)
  vapply(seq_along(lambda_values), function(i) {
    p <- matrix(arr[, , i, drop = FALSE], nrow = dim(arr)[1], ncol = dim(arr)[2],
                dimnames = dimnames(arr)[1:2])
    multiclass_logloss(p, y_val, class_levels)
  }, numeric(1))
}

select_lambda_1se_like <- function(lambda_values, losses, tol = 0.002) {
  ok <- which(is.finite(lambda_values) & lambda_values > 0 & is.finite(losses))
  if (!length(ok)) return(list(lambda = NA_real_, loss = Inf, min_loss = Inf, index = NA_integer_))
  best <- ok[which.min(losses[ok])]
  eligible <- ok[losses[ok] <= losses[best] + tol]
  chosen <- eligible[which.max(lambda_values[eligible])]
  list(lambda = lambda_values[chosen], loss = losses[chosen], min_loss = losses[best], index = chosen)
}

holdout_path <- function(x, y, alpha, control, seed, lambda_values = NULL) {
  tryCatch({
    split <- make_holdout_validation_split(y, control$holdout_val_prop, seed)
    tr <- split$train_idx
    va <- split$val_idx
    warnings <- character()
    path <- if (is.null(lambda_values)) NULL else
      sort(unique(c(lambda_values * 10, lambda_values, lambda_values / 10)), decreasing = TRUE)
    raw <- withCallingHandlers(glmnet::glmnet(
      x = glmnet_design(x[tr, , drop = FALSE]), y = y[tr],
      family = "multinomial", alpha = alpha, lambda = path,
      nlambda = control$holdout_nlambda, standardize = control$standardize,
      type.multinomial = control$type_multinomial, maxit = control$glmnet_maxit),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    lambdas <- if (is.null(lambda_values)) raw$lambda else lambda_values
    loss <- path_logloss_by_lambda(raw, x[va, , drop = FALSE], y[va], lambdas, levels(y))
    if (!is.null(raw$jerr) && raw$jerr != 0L) loss[] <- Inf
    sel <- select_lambda_1se_like(lambdas, loss, control$lambda_1se_tol)
    trace <- data.frame(Alpha = alpha, Lambda = lambdas, Loss = loss,
                        Selected = seq_along(lambdas) == sel$index,
                        NTrain = length(tr), NValidation = length(va),
                        Status = if (is.finite(sel$loss)) "ok" else "solver_failure",
                        Message = paste(warnings, collapse = ";"), stringsAsFactors = FALSE)
    list(lambda = sel$lambda, loss = sel$loss, min_loss = sel$min_loss,
         lambdas = lambdas, losses = loss, trace = trace)
  }, error = function(e) list(lambda = NA_real_, loss = Inf, min_loss = Inf,
      lambdas = numeric(), losses = numeric(),
      trace = data.frame(Alpha = alpha, Lambda = NA_real_, Loss = Inf, Selected = FALSE,
                         NTrain = NA_integer_, NValidation = NA_integer_, Status = "error",
                         Message = conditionMessage(e), stringsAsFactors = FALSE)))
}

tag_tuning_trace <- function(trace, stage, cluster, iter = 0L) {
  trace$Stage <- stage
  trace$Cluster <- cluster
  trace$Iter <- iter
  trace
}

tune_initial_partition <- function(x, y, groups, alpha_grid, control, seed) {
  G <- max(groups)
  losses <- matrix(Inf, G, length(alpha_grid))
  lambdas <- matrix(NA_real_, G, length(alpha_grid))
  trace <- list()
  for (g in seq_len(G)) for (a in seq_along(alpha_grid)) {
    idx <- which(groups == g)
    ans <- holdout_path(x[idx, , drop = FALSE], y[idx], alpha_grid[a], control,
                        scmr_seed(seed, 8300 + G + 17 * g))
    losses[g, a] <- ans$loss
    lambdas[g, a] <- ans$lambda
    trace[[length(trace) + 1L]] <- tag_tuning_trace(ans$trace, "initial", g)
  }
  # Every cluster must have a valid result for a candidate shared alpha.
  score <- colMeans(losses)
  if (!any(is.finite(score))) stop("Initial tuning failed for at least one cluster at every alpha.", call. = FALSE)
  best <- which.min(score)
  list(alpha = alpha_grid[best], lambda = lambdas[, best],
       trace = do.call(rbind, trace))
}

tune_local_partition <- function(x, y, groups, alpha, lambda, control, seed, iter) {
  trace <- list()
  for (g in seq_along(lambda)) {
    idx <- which(groups == g)
    grid <- sort(unique(lambda[g] * c(1, control$local_lambda_factors)), decreasing = TRUE)
    ans <- holdout_path(x[idx, , drop = FALSE], y[idx], alpha, control,
                        scmr_seed(seed, 9100 + 1000 * length(lambda) + 17 * g + iter), grid)
    old <- if (length(ans$losses)) ans$losses[which.min(abs(grid - lambda[g]))] else Inf
    updated <- is.finite(ans$loss) && ans$loss < old - control$lambda_update_tol
    tr <- tag_tuning_trace(ans$trace, "local_lambda", g, iter)
    tr$LambdaBefore <- lambda[g]
    if (updated) lambda[g] <- ans$lambda
    tr$LambdaAfter <- lambda[g]
    trace[[g]] <- tr
  }
  list(lambda = lambda, trace = do.call(rbind, trace))
}

tune_final_partition <- function(x, y, groups, alpha, lambda, control, seed) {
  candidates <- sort(unique(c(alpha, pmin(pmax(alpha + control$final_alpha_delta, .05), 1))))
  losses <- matrix(Inf, length(lambda), length(candidates))
  selected <- matrix(NA_real_, length(lambda), length(candidates))
  trace <- list()
  for (g in seq_along(lambda)) for (a in seq_along(candidates)) {
    idx <- which(groups == g)
    grid <- sort(unique(lambda[g] * c(1, control$local_lambda_factors)), decreasing = TRUE)
    ans <- holdout_path(x[idx, , drop = FALSE], y[idx], candidates[a], control,
                        scmr_seed(seed, 12000 + 1000 * length(lambda) + 29 * g), grid)
    losses[g, a] <- ans$loss
    selected[g, a] <- ans$lambda
    trace[[length(trace) + 1L]] <- tag_tuning_trace(ans$trace, "final_alpha", g)
  }
  score <- colMeans(losses)
  old <- match(alpha, candidates)
  changed <- FALSE
  if (any(is.finite(score))) {
    best <- which.min(score)
    changed <- score[best] < score[old] - control$alpha_update_tol
    if (changed) {
      alpha <- candidates[best]
      lambda <- selected[, best]
    }
  }
  list(alpha = alpha, lambda = lambda, changed = changed,
       trace = do.call(rbind, trace))
}

bind_diagnostics <- function(rows) {
  rows <- Filter(function(z) !is.null(z) && nrow(z) > 0, rows)
  if (!length(rows)) return(data.frame())
  cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(z) {
    for (nm in setdiff(cols, names(z))) z[[nm]] <- NA
    z[, cols, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
