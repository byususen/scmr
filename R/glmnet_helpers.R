make_stratified_foldid <- function(y, nfolds = 5, seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  y <- factor(y)
  tab <- table(y)
  tab <- tab[tab > 0]
  if (length(tab) < 2) stop("Need at least two classes for cv.glmnet().", call. = FALSE)
  nfolds_use <- min(as.integer(nfolds), min(tab))
  if (!is.finite(nfolds_use) || nfolds_use < 3L) {
    stop("Cross-validation needs at least three folds and three observations per class.", call. = FALSE)
  }
  foldid <- integer(length(y))
  for (cl in levels(y)) {
    idx <- which(y == cl)
    if (length(idx) == 0) next
    foldid[idx] <- sample(rep(seq_len(nfolds_use), length.out = length(idx)))
  }
  foldid
}

get_s_value <- function(fit, s_rule = NULL) {
  if (!is.null(s_rule) && is.numeric(s_rule)) return(as.numeric(s_rule)[1])
  if (inherits(fit, "cv.glmnet")) {
    if (is.null(s_rule)) s_rule <- "lambda.1se"
    if (identical(s_rule, "lambda.1se")) return(as.numeric(fit$lambda.1se))
    if (identical(s_rule, "lambda.min")) return(as.numeric(fit$lambda.min))
    return(as.numeric(s_rule)[1])
  }
  if (is.null(s_rule) && !is.null(fit$selected_lambda)) return(as.numeric(fit$selected_lambda)[1])
  if (is.null(s_rule) && !is.null(fit$fixed_lambda)) return(as.numeric(fit$fixed_lambda)[1])
  if (!is.null(fit$lambda) && length(fit$lambda) > 0) return(as.numeric(fit$lambda[1]))
  stop("Cannot determine lambda value.", call. = FALSE)
}

get_s_index <- function(cvfit, s_rule) {
  s_val <- get_s_value(cvfit, s_rule)
  which.min(abs(cvfit$lambda - s_val))
}

predict_multinom_prob <- function(fit, newx, s = NULL) {
  s_use <- if (inherits(fit, "cv.glmnet")) {
    if (is.null(s)) "lambda.1se" else s
  } else {
    get_s_value(fit, s)
  }
  # Linear predictors and a max-shifted softmax: glmnet's "response" can
  # overflow to NaN for large linear predictors (extrapolated local fits).
  arr <- stats::predict(fit, newx = glmnet_design(newx), s = s_use, type = "link")
  if (length(dim(arr)) == 3) {
    link <- matrix(arr[, , 1, drop = FALSE], nrow = dim(arr)[1],
                   ncol = dim(arr)[2], dimnames = dimnames(arr)[1:2])
  } else {
    link <- as.matrix(arr)
  }
  prob <- softmax_rows(link)
  dimnames(prob) <- dimnames(link)
  if (is.null(colnames(prob))) {
    cls <- NULL
    if (!is.null(fit$classnames)) cls <- fit$classnames
    if (is.null(cls) && !is.null(fit$glmnet.fit$classnames)) cls <- fit$glmnet.fit$classnames
    if (!is.null(cls) && length(cls) == ncol(prob)) colnames(prob) <- cls
  }
  prob
}

extract_cv_preval_prob_matrix <- function(cvfit, y_levels, s_rule, n_test) {
  if (is.null(cvfit$fit.preval)) stop("cv.glmnet object does not contain fit.preval. Use keep = TRUE.", call. = FALSE)
  idx <- get_s_index(cvfit, s_rule)
  preval <- cvfit$fit.preval
  if (is.list(preval) && !is.data.frame(preval)) {
    out <- do.call(cbind, lapply(preval, function(mat) as.matrix(mat)[, idx]))
    if (!is.null(names(preval))) colnames(out) <- names(preval)
  } else if (length(dim(preval)) == 3) {
    out <- matrix(preval[, , idx, drop = FALSE], nrow = dim(preval)[1],
                  ncol = dim(preval)[2], dimnames = dimnames(preval)[1:2])
  } else {
    stop("Unsupported fit.preval structure.", call. = FALSE)
  }
  # cv.glmnet keeps prevalidated linear predictors for multinomial models.
  normalize_prob_matrix(softmax_rows(out), y_levels = y_levels, n_test = n_test)
}

glmnet_design <- function(x) {
  if (ncol(x) == 1L) cbind(x, `.scmr_padding` = 0) else x
}

fit_fixed_lambda_glmnet <- function(x, y, alpha, lambda, class_levels,
                                    standardize = TRUE,
                                    type_multinomial = "grouped", maxit = 100000,
                                    weights = NULL) {
  y <- factor(y, levels = class_levels)
  lambda <- as.numeric(lambda)[1]
  if (!is.finite(lambda) || lambda <= 0) stop("lambda must be positive.", call. = FALSE)
  # A decreasing warm-start path that ends exactly at the requested lambda.
  # Nearly separable local data (strong lag-state predictors, tiny lambda) can
  # stop glmnet before the end of a short path, or make it fail while naming a
  # truncated path; a longer path with more warm starts is tried next.
  attempt <- function(n_steps, top) {
    path <- sort(unique(c(exp(seq(log(lambda * top), log(lambda), length.out = n_steps)), lambda)),
                 decreasing = TRUE)
    fit <- tryCatch(glmnet::glmnet(
      x = glmnet_design(x), y = y, family = "multinomial", alpha = alpha, lambda = path,
      standardize = standardize, type.multinomial = type_multinomial, maxit = maxit,
      weights = weights %||% rep(1, nrow(x))), error = function(e) e)
    if (inherits(fit, "error")) return(fit)
    if (!any(abs(fit$lambda - lambda) <= abs(lambda) * 1e-8)) {
      err <- simpleError("The glmnet path did not reach the requested lambda; inspect solver limits.")
      err$fit <- fit
      return(err)
    }
    fit
  }
  fit <- attempt(4L, 100)
  if (inherits(fit, "error")) fit <- attempt(25L, 1000)
  if (inherits(fit, "error")) fit <- attempt(60L, 1e4)
  used <- lambda
  if (inherits(fit, "error")) {
    # Last resort: the smallest lambda glmnet reached (a slightly stronger
    # penalty), reported through the fit and a warning.
    if (is.null(fit$fit) || !length(fit$fit$lambda)) stop(conditionMessage(fit), call. = FALSE)
    partial <- fit$fit
    used <- min(partial$lambda)
    warning(sprintf("glmnet stopped before lambda = %.3g; using lambda = %.3g.", lambda, used), call. = FALSE)
    fit <- partial
  }
  fit$fixed_lambda <- used
  fit$selected_lambda <- used
  fit$requested_lambda <- lambda
  fit
}

active_predictor_count <- function(fit, s = NULL, tol = 1e-8) {
  s_use <- if (is.null(s)) get_s_value(fit, NULL) else s
  cf <- stats::coef(fit, s = s_use)
  active <- character(0)
  if (is.list(cf)) {
    for (mat in cf) {
      mat <- as.matrix(mat)
      vals <- as.numeric(mat[, 1])
      names(vals) <- rownames(mat)
      active <- union(active, names(vals)[abs(vals) > tol])
    }
  }
  length(setdiff(active, "(Intercept)"))
}
