#' Fit spatially clustered or global multinomial regression
#'
#' The package fits global, estimated spatial-cluster, and fixed-membership
#' multinomial models. Study labels such as TCMR belong to the simulation.
#' @param x Finite numeric predictor matrix, optionally sparse.
#' @param y Categorical response with every declared class present in training.
#' @param model One of "scmr", "global", or "fixed_clusters".
#' @param penalty One of "elastic_net", "lasso", "ridge", or "none".
#' @param G Number of estimated clusters; optional consistency check for fixed clusters.
#' @param unit_id One spatial unit identifier per row. Repeated units move together.
#' @param coords Numeric two-column coordinates, required for estimated spatial clusters.
#' @param alpha Optional fixed scalar alpha. Ridge uses zero; lasso uses one.
#' @param alpha_grid Optional elastic-net alpha candidates.
#' @param lambda Optional fixed positive lambda, scalar or one per cluster.
#' @param control Fitting settings from scmr_control().
#' @param initial_cluster Optional initial labels for sorted unique spatial units.
#' @param nfolds Number of global CV folds, at least three.
#' @param lambda_rule Global CV rule: "lambda.1se" or "lambda.min".
#' @param seed Root random seed, or NULL to use the current RNG state.
#' @param cluster Fixed cluster label for each training row.
#' @param ... Reserved for extensions.
#' @return A scmr_fit object with local fits, training probabilities, memberships,
#'   criteria, cluster components, and diagnostic tables. No files are written.
#' @details A spatial model with G=1 uses the global tuning and fitting path.
#'   Information criteria exclude the spatial bonus. Penalized EDF is a local
#'   curvature approximation conditional on the active set and memberships.
#' @export
fit_scmr <- function(x, y, model = c("scmr", "global", "fixed_clusters"),
                     penalty = c("elastic_net", "lasso", "ridge", "none"),
                     G = NULL, unit_id = NULL, coords = NULL, alpha = NULL,
                     alpha_grid = NULL, lambda = NULL, control = scmr_control(),
                     initial_cluster = NULL, nfolds = 5,
                     lambda_rule = c("lambda.1se", "lambda.min"), seed = 123,
                     cluster = NULL, ...) {
  call <- match.call()
  model <- match.arg(model)
  penalty <- match.arg(penalty)
  if (!is.list(control) || any(!names(control) %in% names(scmr_control()))) stop("Unknown control setting.", call. = FALSE)
  control <- utils::modifyList(scmr_control(), control)
  validate_scmr_control(control)
  if (missing(lambda_rule)) lambda_rule <- control$lambda_rule
  lambda_rule <- match.arg(lambda_rule, c("lambda.1se", "lambda.min"))
  if (!is.null(alpha_grid)) control$alpha_grid <- alpha_grid
  alpha_info <- resolve_penalty_alpha(penalty, alpha, alpha_grid %||% control$alpha_grid)
  if (penalty == "none" && !is.null(lambda)) stop("lambda is not used with penalty = 'none'.", call. = FALSE)
  if (penalty != "none" && !control$tune_lambda && is.null(lambda)) {
    stop("Supply lambda when tune_lambda is FALSE.", call. = FALSE)
  }
  if (!is.null(lambda) && (any(!is.finite(lambda)) || any(lambda <= 0))) {
    stop("lambda must contain finite positive values.", call. = FALSE)
  }
  if (model == "scmr" && (is.null(unit_id) || is.null(coords) || is.null(G))) {
    stop("Spatial fitting requires G, unit_id, and coords.", call. = FALSE)
  }
  data <- prepare_scmr_data(x, y, unit_id, coords)
  if (!is.null(G) && (length(G) != 1L || !is.finite(G) || G < 1 || G != floor(G))) {
    stop("G must be a positive integer.", call. = FALSE)
  }
  if (model == "fixed_clusters") {
    if (length(cluster) != nrow(data$x) || anyNA(cluster)) {
      stop("cluster must supply one nonmissing fixed label per row.", call. = FALSE)
    }
    fixed_levels <- if (is.factor(cluster)) levels(droplevels(cluster)) else sort(unique(as.character(cluster)))
    groups <- match(as.character(cluster), fixed_levels)
    if (!is.null(G) && G != length(fixed_levels)) stop("G does not match fixed labels.", call. = FALSE)
    G <- length(fixed_levels)
    per_unit <- split(groups, data$unit_id)
    if (any(vapply(per_unit, function(z) length(unique(z)), integer(1)) != 1L)) {
      stop("A spatial unit cannot belong to multiple fixed clusters.", call. = FALSE)
    }
  } else {
    fixed_levels <- NULL
    groups <- NULL
    if (model == "global") G <- 1L
  }
  if (!is.null(lambda) && !length(lambda) %in% c(1L, G)) {
    stop("lambda must be scalar or contain one value per cluster.", call. = FALSE)
  }
  t0 <- proc.time()[["elapsed"]]
  out <- with_scmr_seed(seed, {
    if (model == "global" || (model == "scmr" && G == 1L)) {
      if (model == "scmr") {
        if (length(unique(data$unit_id)) < control$min_units ||
            any(table(data$y) < resolve_min_class(control$min_per_class, levels(data$y)))) {
          stop("G=1 does not satisfy the spatial feasibility constraints.", call. = FALSE)
        }
      }
      ans <- fit_global_core(data, penalty, alpha_info, lambda, control, seed, nfolds, lambda_rule)
      if (model == "scmr") {
        ans$model <- "scmr"
        class(ans) <- c("scmr_spatial", "scmr_fit")
        ans$convergence_reason <- if (ans$converged) "global_G1_path" else ans$convergence_reason
      }
      ans
    } else {
      fit_partition_core(data, model, penalty, G, groups, fixed_levels, alpha_info,
                         !is.null(alpha), lambda, control, initial_cluster, seed)
    }
  })
  out$call <- call
  out$runtime <- list(total_sec = proc.time()[["elapsed"]] - t0,
                      tuning_sec = out$runtime_tuning)
  out
}

#' Predict from an SCMR package model
#' @param object A fitted SCMR object.
#' @param ... Arguments passed to predict().
#' @return Predicted probabilities, classes, or memberships.
#' @export
predict_scmr <- function(object, ...) stats::predict(object, ...)
