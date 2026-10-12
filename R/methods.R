#' Methods for fitted SCMR models
#' @param object A fitted scmr_fit object.
#' @param newx Predictor matrix for prediction.
#' @param new_unit_id Unit identifiers for spatial prediction.
#' @param new_coords Coordinates for spatial prediction.
#' @param type Return probabilities, classes, or cluster membership weights.
#' @param cluster Fixed cluster labels for fixed-membership prediction.
#' @param membership Rule for the cluster weights of unseen units:
#'   `"proportion"` (share of neighbours), `"potts"` (Potts full conditional
#'   with the fitted phi) or `"majority"` (hard majority, as in SCR).
#' @param parameterization Return native or sum-to-zero coefficients.
#' @param ... Reserved arguments.
#' @return predict returns a matrix or factor; coef returns a cluster-by-class-
#'   by-term array; logLik returns a logLik object with approximate effective df.
#' @name scmr-methods
NULL

# Membership weights of prediction rows. Training units keep their estimated
# label. An unseen unit gets weights from its k nearest training units:
#   "proportion": share of (kernel-weighted) neighbours in each cluster (default);
#   "potts":      Potts full conditional softmax(phi * S_0(h)), S_0 = neighbour
#                 support; the Bayes predictive mixture under the label model;
#   "majority":   all weight on the cluster with the largest support (SCR rule).
assign_test_PP <- function(fit, new_unit_id, new_coords, n = length(new_unit_id),
                           rule = c("proportion", "potts", "majority")) {
  rule <- match.arg(rule)
  G <- fit$G
  if (G == 1L) return(matrix(1, n, 1L, dimnames = list(NULL, "G1")))
  if (length(new_unit_id) != n || anyNA(new_unit_id)) stop("Supply new_unit_id for every prediction row.", call. = FALSE)
  new_unit_id <- as.character(new_unit_id)
  new_units <- unique(new_unit_id)
  unseen <- !new_units %in% fit$unit_levels
  coords <- NULL
  if (any(unseen)) {
    new_coords <- as.matrix(new_coords)
    if (!is.numeric(new_coords) || nrow(new_coords) != n || ncol(new_coords) != 2L ||
        any(!is.finite(new_coords))) stop("Supply finite new_coords for unseen units.", call. = FALSE)
    coords <- get_unit_coords(new_unit_id, new_coords)
    coords <- coords[match(new_units, sort(unique(new_unit_id))), , drop = FALSE]
  }
  phi <- fit$phi %||% fit$control$phi
  pp <- matrix(0, length(new_units), G, dimnames = list(new_units, paste0("G", seq_len(G))))
  for (i in seq_along(new_units)) {
    if (!unseen[i]) {
      pp[i, ] <- fit$PP_unit[new_units[i], ]
    } else {
      distance <- sqrt(rowSums(sweep(fit$train_unit_coords, 2, coords[i, ], "-")^2))
      nn <- order(distance)[seq_len(min(fit$control$k_neighbors, length(distance)))]
      weights <- if (fit$control$weight_type == "binary") rep(1, length(nn)) else {
        # Preserve the training bandwidth and normalize in log space.
        logw <- -(distance[nn] / fit$bandwidth_used)^2
        exp(logw - max(logw))
      }
      support <- colSums(fit$PP_unit[nn, , drop = FALSE] * weights)
      pp[i, ] <- switch(rule,
        proportion = support / sum(support),
        potts = { a <- phi * support; e <- exp(a - max(a)); e / sum(e) },
        majority = { z <- numeric(G); z[which.max(support)] <- 1; z })
    }
  }
  normalize_prob_matrix(pp[new_unit_id, , drop = FALSE], colnames(pp), n)
}

#' @rdname scmr-methods
#' @export
predict.scmr_fit <- function(object, newx, new_unit_id = NULL, new_coords = NULL,
                              type = c("prob", "class", "membership"), cluster = NULL,
                              membership = c("proportion", "potts", "majority"), ...) {
  type <- match.arg(type)
  membership <- match.arg(membership)
  newx <- check_newx(newx, object$x_colnames)
  n <- nrow(newx)
  if (object$model == "fixed_clusters") {
    if (length(cluster) != n || anyNA(cluster)) stop("Supply one fixed cluster label per prediction row.", call. = FALSE)
    group <- match(as.character(cluster), object$fixed_cluster_levels)
    if (anyNA(group)) stop("Prediction contains unknown fixed cluster labels.", call. = FALSE)
    pp <- matrix(0, n, object$G, dimnames = list(NULL, paste0("G", seq_len(object$G))))
    pp[cbind(seq_len(n), group)] <- 1
  } else if (object$model == "global") {
    pp <- matrix(1, n, 1L, dimnames = list(NULL, "G1"))
  } else {
    pp <- assign_test_PP(object, new_unit_id, new_coords, n, membership)
  }
  if (type == "membership") return(pp)
  prob <- matrix(0, n, length(object$class_levels), dimnames = list(NULL, object$class_levels))
  for (g in seq_along(object$fits)) prob <- prob + engine_predict(object$fits[[g]], newx) * pp[, g]
  prob <- normalize_prob_matrix(prob, object$class_levels, n)
  if (type == "prob") return(prob)
  factor(object$class_levels[max.col(prob, ties.method = "first")], levels = object$class_levels)
}

#' @rdname scmr-methods
#' @export
coef.scmr_fit <- function(object, parameterization = c("sum_to_zero", "native"), ...) {
  parameterization <- match.arg(parameterization)
  out <- array(0, c(object$G, length(object$class_levels), length(object$x_colnames) + 1L),
               dimnames = list(paste0("G", seq_len(object$G)), object$class_levels,
                               c("(Intercept)", object$x_colnames)))
  for (g in seq_len(object$G)) out[g, , ] <- engine_coef(object$fits[[g]], parameterization)
  out
}

#' @rdname scmr-methods
#' @export
logLik.scmr_fit <- function(object, ...) {
  structure(object$criteria$logLik, class = "logLik",
            df = object$criteria$df_effective_total, nobs = object$criteria$n_rows)
}

#' @rdname scmr-methods
#' @param x A fitted scmr_fit object.
#' @export
print.scmr_fit <- function(x, ...) {
  cat("SCMR multinomial fit\n")
  cat("  Model:", x$model, "| penalty:", x$penalty, "| clusters:", x$G, "\n")
  cat("  Classes:", paste(x$class_levels, collapse = ", "), "\n")
  cat("  Converged:", x$converged, "|", x$convergence_reason, "\n")
  invisible(x)
}

# Backward-compatible subclass dispatch.
#' @export
predict.scmr_global <- function(object, newx, type = c("prob", "class", "membership"), ...) {
  predict.scmr_fit(object, newx, type = match.arg(type), ...)
}
#' @export
predict.scmr_spatial <- function(object, newx, new_unit_id = NULL, new_coords = NULL,
                                 type = c("prob", "class", "membership"), ...) {
  predict.scmr_fit(object, newx, new_unit_id, new_coords, match.arg(type), ...)
}
#' @export
print.scmr_global <- function(x, ...) print.scmr_fit(x, ...)
#' @export
print.scmr_spatial <- function(x, ...) print.scmr_fit(x, ...)
