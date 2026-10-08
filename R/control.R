#' Control parameters for SCMR fitting
#'
#' @param phi Spatial smoothness coefficient for the Potts-style clustering term.
#' @param k_neighbors Number of nearest neighbors used to form spatial weights.
#' @param weight_type Spatial weight type: either `"binary"` or `"exp"`.
#' @param weight_bandwidth Optional bandwidth for exponential spatial weights.
#' @param row_standardize_weights If `TRUE`, row-standardize the spatial weight matrix.
#' @param max_iter Maximum number of SCMR membership-update iterations.
#' @param init_method Initial spatial partition method: `"kmeans"` (coordinates),
#'   `"random"`, or `"coefficient"` (k-means on local coefficient vectors; see
#'   [scmr_local_coefficients()]).
#' @param init_tries Number of random/k-means initialization attempts for feasible clusters.
#' @param update_order Unit update order in each SCMR iteration, either `"random"` or `"fixed"`.
#' @param alpha_grid Candidate alpha values used when tuning elastic-net models.
#' @param lambda_rule Selection rule for global `cv.glmnet()`, usually `"lambda.1se"` or `"lambda.min"`.
#' @param holdout_val_prop Validation proportion for cluster-local holdout tuning.
#' @param holdout_nlambda Number of lambda values in the initial holdout path.
#' @param lambda_1se_tol Fixed loss tolerance for choosing the largest acceptable holdout lambda; not a standard error.
#' @param lambda_update_tol Minimum validation-loss improvement needed to update local lambda.
#' @param alpha_update_tol Minimum validation-loss improvement needed to update alpha after convergence.
#' @param local_lambda_factors Multipliers around current lambda for local lambda updates.
#' @param final_alpha_delta Candidate alpha changes around the current alpha after convergence.
#' @param tune_lambda If `TRUE`, tune cluster-local lambdas during SCMR fitting.
#' @param tune_alpha_after_convergence If `TRUE`, optionally retune alpha after membership convergence.
#' @param min_units Minimum number of unique spatial units per cluster.
#' @param min_per_class Minimum number of observations per class per cluster.
#' @param standardize Passed to `glmnet::glmnet()` and `glmnet::cv.glmnet()`.
#' @param type_multinomial Passed to glmnet for multinomial models.
#' @param tiny_movement_max_units Number of moved units treated as a tiny update.
#' @param tiny_movement_rate_tol Relative moved-unit rate treated as a tiny update.
#' @param tiny_movement_revert If `TRUE`, revert the final tiny movement before stopping.
#' @param criterion_balance_gamma Cluster-balance exponent used in CB-BIC.
#' @param glmnet_maxit Maximum iterations for a glmnet fit.
#' @param multinom_maxit Maximum iterations for an unpenalized multinomial fit.
#' @param multinom_maxnwts Maximum weights for the unpenalized multinomial engine.
#' @param selection_tol Absolute native-coefficient support threshold.
#' @param lambda_scale `"mean"` (glmnet scale, tuned during iterations; legacy)
#'   or `"sum"`. With `"sum"`, one sum-scale penalty strength
#'   kappa = mean(lambda) * n / G, tuned on the first initial partition, is used
#'   for every cluster and start; predictors are standardized once globally and
#'   no tuning happens inside the alternation, so the penalized Potts objective
#'   is non-decreasing (monotone algorithm) and comparable across starts.
#' @param update_memberships If `FALSE`, keep the initial partition (two-stage
#'   benchmark: spatial clustering followed by local elastic-net fits).
#' @param n_starts Number of initial partitions; the fit with the largest
#'   penalized objective is returned.
#' @param start_methods Optional character vector of initialization methods for
#'   the starts, recycled to `n_starts`; defaults to `init_method` then
#'   alternating `"coefficient"` and `"kmeans"`.
#' @param coef_init_k Number of nearest spatial units in each local fit used by
#'   coefficient-space initialization.
#' @param coef_init_alpha,coef_init_lambda Elastic-net settings of the local fits.
#' @param coef_init_anchors Maximum number of anchor units at which local fits
#'   are computed; other units take the features of their nearest anchor.
#' @param verbose Print iteration progress when true.
#' @return A named list of control settings.
#' @export
scmr_control <- function(
    phi = 1,
    k_neighbors = 5,
    weight_type = c("binary", "exp"),
    weight_bandwidth = NULL,
    row_standardize_weights = FALSE,
    max_iter = 10,
    init_method = c("kmeans", "random", "coefficient"),
    init_tries = 200,
    update_order = c("random", "fixed"),
    alpha_grid = c(0.6, 0.7, 0.8),
    lambda_rule = "lambda.1se",
    holdout_val_prop = 0.30,
    holdout_nlambda = 60,
    lambda_1se_tol = 0.002,
    lambda_update_tol = 0.001,
    alpha_update_tol = 0.003,
    local_lambda_factors = c(3, 2, 1, 1/2, 1/3),
    final_alpha_delta = c(-0.1, 0, 0.1),
    tune_lambda = TRUE,
    tune_alpha_after_convergence = TRUE,
    min_units = 10,
    min_per_class = 5,
    standardize = TRUE,
    type_multinomial = "grouped",
    tiny_movement_max_units = 5,
    tiny_movement_rate_tol = 0.001,
    tiny_movement_revert = TRUE,
    criterion_balance_gamma = 1,
    glmnet_maxit = 100000,
    multinom_maxit = 1000,
    multinom_maxnwts = 200000,
    selection_tol = 1e-8,
    lambda_scale = c("mean", "sum"),
    update_memberships = TRUE,
    n_starts = 1L,
    start_methods = NULL,
    coef_init_k = 50L,
    coef_init_alpha = 0.2,
    coef_init_lambda = 0.05,
    coef_init_anchors = 300L,
    verbose = FALSE) {
  weight_type <- match.arg(weight_type)
  init_method <- match.arg(init_method)
  lambda_scale <- match.arg(lambda_scale)
  update_order <- match.arg(update_order)

  out <- list(
    phi = phi,
    k_neighbors = k_neighbors,
    weight_type = weight_type,
    weight_bandwidth = weight_bandwidth,
    row_standardize_weights = row_standardize_weights,
    max_iter = max_iter,
    init_method = init_method,
    init_tries = init_tries,
    update_order = update_order,
    alpha_grid = alpha_grid,
    lambda_rule = lambda_rule,
    holdout_val_prop = holdout_val_prop,
    holdout_nlambda = holdout_nlambda,
    lambda_1se_tol = lambda_1se_tol,
    lambda_update_tol = lambda_update_tol,
    alpha_update_tol = alpha_update_tol,
    local_lambda_factors = local_lambda_factors,
    final_alpha_delta = final_alpha_delta,
    tune_lambda = tune_lambda,
    tune_alpha_after_convergence = tune_alpha_after_convergence,
    min_units = min_units,
    min_per_class = min_per_class,
    standardize = standardize,
    type_multinomial = type_multinomial,
    tiny_movement_max_units = tiny_movement_max_units,
    tiny_movement_rate_tol = tiny_movement_rate_tol,
    tiny_movement_revert = tiny_movement_revert,
    criterion_balance_gamma = criterion_balance_gamma,
    glmnet_maxit = glmnet_maxit,
    multinom_maxit = multinom_maxit,
    multinom_maxnwts = multinom_maxnwts,
    selection_tol = selection_tol,
    lambda_scale = lambda_scale,
    update_memberships = update_memberships,
    n_starts = n_starts,
    start_methods = start_methods,
    coef_init_k = coef_init_k,
    coef_init_alpha = coef_init_alpha,
    coef_init_lambda = coef_init_lambda,
    coef_init_anchors = coef_init_anchors,
    verbose = verbose
  )
  validate_scmr_control(out)
  out
}

validate_scmr_control <- function(control) {
  for (nm in c("standardize", "row_standardize_weights", "tune_lambda", "tune_alpha_after_convergence",
               "tiny_movement_revert", "update_memberships", "verbose")) {
    if (!is.logical(control[[nm]]) || length(control[[nm]]) != 1L || is.na(control[[nm]])) {
      stop(nm, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  choices <- list(weight_type = c("binary", "exp"), init_method = c("kmeans", "random", "coefficient"),
    update_order = c("random", "fixed"), lambda_rule = c("lambda.min", "lambda.1se"),
    type_multinomial = c("grouped", "ungrouped"))
  for (nm in names(choices)) {
    if (length(control[[nm]]) != 1L || !control[[nm]] %in% choices[[nm]]) stop("Invalid ", nm, call. = FALSE)
  }
  integer_fields <- c("k_neighbors", "max_iter", "init_tries", "holdout_nlambda",
                      "min_units", "glmnet_maxit", "multinom_maxit", "multinom_maxnwts",
                      "n_starts", "coef_init_k", "coef_init_anchors")
  for (nm in integer_fields) {
    z <- control[[nm]]
    if (length(z) != 1L || !is.finite(z) || z < 1 || z != floor(z)) {
      stop(nm, " must be a positive integer.", call. = FALSE)
    }
  }
  for (nm in c("phi", "lambda_1se_tol", "lambda_update_tol", "alpha_update_tol",
               "tiny_movement_max_units", "tiny_movement_rate_tol",
               "criterion_balance_gamma", "selection_tol")) {
    z <- control[[nm]]
    if (length(z) != 1L || !is.finite(z) || z < 0) {
      stop(nm, " must be finite and nonnegative.", call. = FALSE)
    }
  }
  if (!is.null(control$start_methods) &&
      (!is.character(control$start_methods) || !length(control$start_methods) ||
       any(!control$start_methods %in% c("kmeans", "random", "coefficient")))) {
    stop("start_methods must contain kmeans, random, or coefficient.", call. = FALSE)
  }
  if (length(control$coef_init_alpha) != 1L || !is.finite(control$coef_init_alpha) ||
      control$coef_init_alpha < 0 || control$coef_init_alpha > 1) stop("coef_init_alpha must lie in [0, 1].", call. = FALSE)
  if (length(control$coef_init_lambda) != 1L || !is.finite(control$coef_init_lambda) ||
      control$coef_init_lambda <= 0) stop("coef_init_lambda must be positive.", call. = FALSE)
  if (length(control$holdout_val_prop) != 1L || !is.finite(control$holdout_val_prop) || control$holdout_val_prop <= 0 ||
      control$holdout_val_prop >= 1) stop("holdout_val_prop must lie in (0, 1).", call. = FALSE)
  if (!length(control$min_per_class) || any(!is.finite(control$min_per_class)) ||
      any(control$min_per_class < 1 | control$min_per_class != floor(control$min_per_class))) {
    stop("min_per_class must contain positive integers.", call. = FALSE)
  }
  if (!control$type_multinomial %in% c("grouped", "ungrouped")) {
    stop("type_multinomial must be grouped or ungrouped.", call. = FALSE)
  }
  if (!length(control$alpha_grid) || any(!is.finite(control$alpha_grid)) ||
      any(control$alpha_grid < 0 | control$alpha_grid > 1)) stop("Invalid alpha_grid.", call. = FALSE)
  if (!length(control$local_lambda_factors) || any(!is.finite(control$local_lambda_factors)) ||
      any(control$local_lambda_factors <= 0)) stop("Invalid local_lambda_factors.", call. = FALSE)
  if (!length(control$final_alpha_delta) || any(!is.finite(control$final_alpha_delta))) {
    stop("Invalid final_alpha_delta.", call. = FALSE)
  }
  if (!is.null(control$weight_bandwidth) &&
      (length(control$weight_bandwidth) != 1L || !is.finite(control$weight_bandwidth) ||
       control$weight_bandwidth <= 0)) stop("weight_bandwidth must be positive.", call. = FALSE)
  invisible(control)
}

resolve_penalty_alpha <- function(penalty = c("elastic_net", "lasso", "ridge", "none"),
                                  alpha = NULL,
                                  alpha_grid = NULL) {
  penalty <- match.arg(penalty)
  if (!is.null(alpha) && (length(alpha) != 1L || !is.finite(alpha))) stop("alpha must be a finite scalar.", call. = FALSE)
  if (penalty == "none") {
    if (!is.null(alpha)) stop("alpha is not used with penalty = 'none'.", call. = FALSE)
    return(list(penalty = penalty, alpha_grid = numeric(), alpha_default = NA_real_))
  }

  if (penalty == "ridge") {
    if (!is.null(alpha) && any(abs(as.numeric(alpha)) > .Machine$double.eps^0.5)) {
      stop("For penalty = 'ridge', alpha must be 0 or NULL.", call. = FALSE)
    }
    return(list(penalty = penalty, alpha_grid = 0, alpha_default = 0))
  }

  if (penalty == "lasso") {
    if (!is.null(alpha) && any(abs(as.numeric(alpha) - 1) > .Machine$double.eps^0.5)) {
      stop("For penalty = 'lasso', alpha must be 1 or NULL.", call. = FALSE)
    }
    return(list(penalty = penalty, alpha_grid = 1, alpha_default = 1))
  }

  if (!is.null(alpha)) {
    ag <- as.numeric(alpha)
  } else if (!is.null(alpha_grid)) {
    ag <- as.numeric(alpha_grid)
  } else {
    ag <- c(0.6, 0.7, 0.8)
  }

  if (!length(ag) || any(!is.finite(ag) | ag <= 0 | ag > 1)) {
    stop("Elastic-net alpha values must be in (0, 1].", call. = FALSE)
  }
  ag <- sort(unique(ag))
  if (!is.null(alpha) && length(ag) != 1L) stop("alpha must be a scalar.", call. = FALSE)
  list(penalty = penalty, alpha_grid = ag, alpha_default = stats::median(ag))
}
