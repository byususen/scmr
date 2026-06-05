#' Fit SCMR or a global penalized multinomial benchmark
#'
#' `fit_scmr()` is the main user-facing function. Use `model = "global"`
#' for a global penalized multinomial model, or `model = "scmr"` for
#' spatially clustered multinomial regression. Use `penalty` to choose ridge,
#' lasso, or elastic net.
#'
#' @param x Numeric matrix or sparse matrix of predictors.
#' @param y Factor or character vector of class labels.
#' @param model Either `"global"` or `"scmr"`.
#' @param penalty One of `"ridge"`, `"lasso"`, or `"elastic_net"`.
#' @param G Number of spatial clusters. Required when `model = "scmr"`.
#' @param unit_id Spatial unit identifier. Required when `model = "scmr"`.
#' @param coords Two-column coordinate matrix. Required when `model = "scmr"`.
#' @param alpha Optional fixed alpha value. Ridge uses 0 and lasso uses 1.
#' @param alpha_grid Optional elastic-net alpha grid.
#' @param lambda Optional fixed lambda. For SCMR, use a scalar or one value per cluster.
#' @param control Control list from [scmr_control()].
#' @param initial_cluster Optional initial cluster labels for unique spatial units.
#' @param nfolds Number of CV folds for the global model.
#' @param lambda_rule `"lambda.1se"` or `"lambda.min"` when global lambda is tuned.
#' @param seed Random seed.
#' @param ... Reserved for future extensions.
#' @return A fitted `scmr_global` or `scmr_spatial` object.
#' @export
fit_scmr <- function(x, y,
                     model = c("scmr", "global"),
                     penalty = c("elastic_net", "lasso", "ridge"),
                     G = NULL,
                     unit_id = NULL,
                     coords = NULL,
                     alpha = NULL,
                     alpha_grid = NULL,
                     lambda = NULL,
                     control = scmr_control(),
                     initial_cluster = NULL,
                     nfolds = 5,
                     lambda_rule = c("lambda.1se", "lambda.min"),
                     seed = 123,
                     ...) {
  model <- match.arg(model)
  penalty <- match.arg(penalty)
  lambda_rule <- match.arg(lambda_rule)

  if (model == "global") {
    return(fit_global_scmr(
      x = x,
      y = y,
      penalty = penalty,
      alpha = alpha,
      alpha_grid = alpha_grid,
      lambda = lambda,
      nfolds = nfolds,
      lambda_rule = lambda_rule,
      seed = seed,
      standardize = control$standardize,
      type_multinomial = control$type_multinomial
    ))
  }

  if (is.null(G)) stop("G must be supplied when model = 'scmr'.", call. = FALSE)
  if (is.null(unit_id)) stop("unit_id must be supplied when model = 'scmr'.", call. = FALSE)
  if (is.null(coords)) stop("coords must be supplied when model = 'scmr'.", call. = FALSE)

  fit_spatial_scmr(
    x = x,
    y = y,
    unit_id = unit_id,
    coords = coords,
    G = G,
    penalty = penalty,
    alpha = alpha,
    alpha_grid = alpha_grid,
    lambda = lambda,
    control = control,
    initial_cluster = initial_cluster,
    seed = seed
  )
}

#' Predict from an SCMR package model
#'
#' Thin wrapper around the package's S3 prediction methods.
#' @param object A fitted `scmr_global` or `scmr_spatial` object.
#' @param ... Arguments passed to `predict()`.
#' @export
predict_scmr <- function(object, ...) {
  stats::predict(object, ...)
}
