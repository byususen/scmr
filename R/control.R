#' Control parameters for SCMR fitting
#'
#' @param phi Spatial smoothness coefficient for the Potts-style clustering term.
#' @param k_neighbors Number of nearest neighbors used to form spatial weights.
#' @param weight_type Spatial weight type: either `"binary"` or `"exp"`.
#' @param weight_bandwidth Optional bandwidth for exponential spatial weights.
#' @param row_standardize_weights If `TRUE`, row-standardize the spatial weight matrix.
#' @param max_iter Maximum number of SCMR membership-update iterations.
#' @param init_method Initial spatial partition method, either `"kmeans"` or `"random"`.
#' @param init_tries Number of random/k-means initialization attempts for feasible clusters.
#' @param update_order Unit update order in each SCMR iteration, either `"random"` or `"fixed"`.
#' @param alpha_grid Candidate alpha values used when tuning elastic-net models.
#' @param lambda_rule Selection rule for global `cv.glmnet()`, usually `"lambda.1se"` or `"lambda.min"`.
#' @param holdout_val_prop Validation proportion for cluster-local holdout tuning.
#' @param holdout_nlambda Number of lambda values in the initial holdout path.
#' @param lambda_1se_tol Tolerance for the largest-lambda 1-SE-like holdout rule.
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
#' @return A named list of control settings.
#' @export
scmr_control <- function(
    phi = 1,
    k_neighbors = 5,
    weight_type = c("binary", "exp"),
    weight_bandwidth = NULL,
    row_standardize_weights = FALSE,
    max_iter = 10,
    init_method = c("kmeans", "random"),
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
    criterion_balance_gamma = 1) {
  weight_type <- match.arg(weight_type)
  init_method <- match.arg(init_method)
  update_order <- match.arg(update_order)

  list(
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
    criterion_balance_gamma = criterion_balance_gamma
  )
}

resolve_penalty_alpha <- function(penalty = c("elastic_net", "lasso", "ridge"),
                                  alpha = NULL,
                                  alpha_grid = NULL) {
  penalty <- match.arg(penalty)

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

  ag <- sort(unique(ag[is.finite(ag) & ag > 0 & ag <= 1]))
  if (length(ag) == 0) stop("Elastic-net alpha values must be in (0, 1].", call. = FALSE)
  list(penalty = penalty, alpha_grid = ag, alpha_default = stats::median(ag))
}
