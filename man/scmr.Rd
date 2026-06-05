\name{scmr}
\alias{scmr}
\alias{fit_scmr}
\alias{fit_global_scmr}
\alias{fit_spatial_scmr}
\alias{scmr_control}
\alias{simulate_scmr_data}
\title{Spatially clustered multinomial regression}
\description{
The \pkg{scmr} package fits global penalized multinomial regression and
spatially clustered multinomial regression for multiclass categorical outcomes.
The main public interface is \code{fit_scmr()}, with \code{model = "global"}
or \code{model = "scmr"}. Penalties can be selected with \code{penalty =
"ridge"}, \code{"lasso"}, or \code{"elastic_net"}.
}
\usage{
fit_scmr(x, y, model = c("scmr", "global"),
         penalty = c("elastic_net", "lasso", "ridge"),
         G = NULL, unit_id = NULL, coords = NULL,
         alpha = NULL, alpha_grid = NULL, lambda = NULL,
         control = scmr_control(), initial_cluster = NULL,
         nfolds = 5, lambda_rule = c("lambda.1se", "lambda.min"),
         seed = 123, ...)

scmr_control(phi = 1, k_neighbors = 5, weight_type = c("binary", "exp"),
             weight_bandwidth = NULL, row_standardize_weights = FALSE,
             max_iter = 10, init_method = c("kmeans", "random"),
             init_tries = 200, update_order = c("random", "fixed"),
             alpha_grid = c(0.6, 0.7, 0.8), lambda_rule = "lambda.1se",
             holdout_val_prop = 0.30, holdout_nlambda = 60,
             lambda_1se_tol = 0.002, lambda_update_tol = 0.001,
             alpha_update_tol = 0.003,
             local_lambda_factors = c(3, 2, 1, 1/2, 1/3),
             final_alpha_delta = c(-0.1, 0, 0.1), tune_lambda = TRUE,
             tune_alpha_after_convergence = TRUE, min_units = 10,
             min_per_class = 5, standardize = TRUE,
             type_multinomial = "grouped", tiny_movement_max_units = 5,
             tiny_movement_rate_tol = 0.001, tiny_movement_revert = TRUE,
             criterion_balance_gamma = 1)
}
\arguments{
\item{x}{Numeric matrix or sparse matrix of predictors.}
\item{y}{Factor or character response vector.}
\item{model}{Either \code{"global"} or \code{"scmr"}.}
\item{penalty}{One of \code{"ridge"}, \code{"lasso"}, or \code{"elastic_net"}.}
\item{G}{Number of spatial clusters for \code{model = "scmr"}.}
\item{unit_id}{Spatial unit identifier for \code{model = "scmr"}.}
\item{coords}{Two-column coordinate matrix for \code{model = "scmr"}.}
\item{alpha}{Optional fixed alpha value.}
\item{alpha_grid}{Optional alpha grid for elastic-net tuning.}
\item{lambda}{Optional fixed lambda; for SCMR, scalar or length \code{G}.}
\item{control}{A list from \code{scmr_control()}.}
\item{initial_cluster}{Optional initial labels for unique spatial units.}
\item{nfolds}{Number of folds for the global model.}
\item{lambda_rule}{Global lambda selection rule when lambda is not fixed.}
\item{seed}{Random seed.}
\item{...}{Reserved for future use.}
}
\value{
A fitted \code{scmr_global} or \code{scmr_spatial} object.
}
