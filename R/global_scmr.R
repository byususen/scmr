#' Fit a global penalized multinomial regression model
#'
#' This is the non-spatial benchmark used by [fit_scmr()] when
#' `model = "global"`.
#'
#' @param x Numeric matrix or sparse matrix of predictors.
#' @param y Factor or character vector of class labels.
#' @param penalty Penalty type: `"ridge"`, `"lasso"`, or `"elastic_net"`.
#' @param alpha Fixed alpha value. For ridge this is 0; for lasso this is 1.
#' @param alpha_grid Candidate alpha values for elastic-net tuning.
#' @param lambda Optional fixed lambda value. If `NULL`, lambda is chosen by `cv.glmnet()`.
#' @param nfolds Number of stratified folds for cross-validation.
#' @param lambda_rule `"lambda.1se"` or `"lambda.min"` when `lambda = NULL`.
#' @param seed Random seed for fold assignment.
#' @param standardize Passed to `glmnet::glmnet()` or `glmnet::cv.glmnet()`.
#' @param type_multinomial Passed to glmnet for multinomial models.
#' @return An object of class `scmr_global`.
#' @export
fit_global_scmr <- function(x, y,
                            penalty = c("elastic_net", "lasso", "ridge"),
                            alpha = NULL,
                            alpha_grid = NULL,
                            lambda = NULL,
                            nfolds = 5,
                            lambda_rule = c("lambda.1se", "lambda.min"),
                            seed = 123,
                            standardize = TRUE,
                            type_multinomial = "grouped") {
  penalty <- match.arg(penalty)
  lambda_rule <- match.arg(lambda_rule)
  alpha_info <- resolve_penalty_alpha(penalty, alpha = alpha, alpha_grid = alpha_grid)
  y <- factor(y)
  class_levels <- levels(y)

  if (!is.null(lambda)) {
    alpha_use <- if (!is.null(alpha)) as.numeric(alpha)[1] else alpha_info$alpha_default
    lambda_use <- as.numeric(lambda)[1]
    fit <- fit_fixed_lambda_glmnet(
      x = x, y = y, alpha = alpha_use, lambda = lambda_use,
      class_levels = class_levels, standardize = standardize,
      type_multinomial = type_multinomial
    )
    prob_train <- predict_multinom_prob(fit, x, lambda_use)
    prob_train <- normalize_prob_matrix(prob_train, class_levels, length(y))
    criteria <- compute_global_criteria(
      fit, x, y, prob_train, alpha_use, lambda_use, class_levels,
      standardize = standardize
    )
    out <- list(
      model = "global",
      penalty = penalty,
      fit = fit,
      alpha = alpha_use,
      lambda = lambda_use,
      lambda_rule = "fixed",
      class_levels = class_levels,
      tuning = data.frame(Alpha = alpha_use, LambdaSelected = lambda_use,
                          LambdaRule = "fixed", stringsAsFactors = FALSE),
      prob_train = prob_train,
      criteria = criteria,
      x_colnames = colnames(x)
    )
    class(out) <- "scmr_global"
    return(out)
  }

  foldid <- make_stratified_foldid(y, nfolds = nfolds, seed = seed)
  best <- NULL
  tuning <- list()

  for (aa in alpha_info$alpha_grid) {
    t0 <- Sys.time()
    cvfit <- glmnet::cv.glmnet(
      x = x,
      y = y,
      family = "multinomial",
      alpha = aa,
      foldid = foldid,
      type.measure = "deviance",
      standardize = standardize,
      type.multinomial = type_multinomial,
      maxit = 100000,
      keep = TRUE
    )
    lambda_selected <- get_s_value(cvfit, lambda_rule)
    selected_idx <- get_s_index(cvfit, lambda_rule)
    selected_cvm <- as.numeric(cvfit$cvm[selected_idx])
    tuning[[length(tuning) + 1L]] <- data.frame(
      Alpha = aa,
      LambdaRule = lambda_rule,
      LambdaSelected = lambda_selected,
      LambdaMin = cvfit$lambda.min,
      Lambda1SE = cvfit$lambda.1se,
      CVMSelected = selected_cvm,
      CVMMin = min(cvfit$cvm, na.rm = TRUE),
      RuntimeSec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
      stringsAsFactors = FALSE
    )
    if (is.null(best) || selected_cvm < best$cvm_selected) {
      best <- list(fit = cvfit, alpha = aa, lambda = lambda_selected,
                   cvm_selected = selected_cvm)
    }
  }

  prob_train <- predict_multinom_prob(best$fit, x, lambda_rule)
  prob_train <- normalize_prob_matrix(prob_train, class_levels, length(y))
  criteria <- compute_global_criteria(
    best$fit, x, y, prob_train, best$alpha, lambda_rule, class_levels,
    standardize = standardize
  )

  out <- list(
    model = "global",
    penalty = penalty,
    fit = best$fit,
    alpha = best$alpha,
    lambda = best$lambda,
    lambda_rule = lambda_rule,
    class_levels = class_levels,
    tuning = do.call(rbind, tuning),
    prob_train = prob_train,
    criteria = criteria,
    x_colnames = colnames(x)
  )
  class(out) <- "scmr_global"
  out
}

#' @export
predict.scmr_global <- function(object, newx, type = c("prob", "class"), ...) {
  type <- match.arg(type)
  s_use <- if (identical(object$lambda_rule, "fixed")) object$lambda else object$lambda_rule
  prob <- predict_multinom_prob(object$fit, newx, s_use)
  prob <- normalize_prob_matrix(prob, object$class_levels, nrow(newx))
  if (type == "prob") return(prob)
  factor(object$class_levels[max.col(prob)], levels = object$class_levels)
}

#' @export
print.scmr_global <- function(x, ...) {
  cat("Global penalized multinomial regression\n")
  cat("  Penalty:", x$penalty, "\n")
  cat("  Classes:", paste(x$class_levels, collapse = ", "), "\n")
  cat("  Alpha:", x$alpha, "\n")
  cat("  Lambda:", signif(x$lambda, 5), "\n")
  invisible(x)
}
