article_subset <- function(data, index) {
  lapply(data, function(z) {
    if (length(dim(z)) == 3L) z[index, , , drop = FALSE]
    else if (length(dim(z)) == 2L) z[index, , drop = FALSE] else z[index]
  })
}

article_ari <- function(truth, predicted) {
  if (anyNA(truth) || length(truth) < 2L) return(NA_real_)
  tab <- table(truth, predicted)
  pairs <- function(x) sum(x * (x - 1) / 2)
  n2 <- pairs(sum(tab)); a <- pairs(rowSums(tab)); b <- pairs(colSums(tab))
  expected <- a * b / n2; denom <- (a + b) / 2 - expected
  if (abs(denom) < 1e-12) return(if (a == b && pairs(tab) == a) 1 else 0)
  (pairs(tab) - expected) / denom
}

article_local_coef <- function(fit, pp) {
  b <- stats::coef(fit, parameterization = "sum_to_zero")
  out <- array(0, c(nrow(pp), dim(b)[2:3]),
    dimnames = list(NULL, dimnames(b)[[2]], dimnames(b)[[3]]))
  for (cc in seq_len(dim(b)[2])) for (j in seq_len(dim(b)[3])) out[, cc, j] <- pp %*% b[, cc, j]
  out
}

article_support <- function(selected, truth) {
  tp <- sum(selected & truth); fp <- sum(selected & !truth)
  fn <- sum(!selected & truth); tn <- sum(!selected & !truth)
  div <- function(a, b) if (b > 0) a / b else NA_real_
  data.frame(TP = tp, FP = fp, FN = fn, TN = tn, TPR = div(tp, tp + fn),
    FPR = div(fp, fp + tn), Precision = div(tp, tp + fp), F1 = div(2 * tp, 2 * tp + fp + fn))
}

article_parameters <- function(truth, estimate, region, active, per_observation = FALSE) {
  groups <- if (per_observation) as.character(seq_len(dim(truth)[1])) else as.character(region)
  result <- list()
  for (rg in unique(groups)) {
    ii <- which(groups == rg)
    grid <- expand.grid(Class = dimnames(truth)[[2]], Term = dimnames(truth)[[3]], stringsAsFactors = FALSE)
    grid$Group <- rg; grid$N <- length(ii)
    grid$Role <- c("intercept", ifelse(seq_len(dim(truth)[3] - 1L) %in% active, "active", "inactive"))[rep(seq_len(dim(truth)[3]), each = dim(truth)[2])]
    grid$TrueMean <- as.vector(apply(truth[ii, , , drop = FALSE], c(2, 3), mean))
    grid$EstimatedMean <- as.vector(apply(estimate[ii, , , drop = FALSE], c(2, 3), mean))
    diff <- estimate[ii, , , drop = FALSE] - truth[ii, , , drop = FALSE]
    grid$Bias <- as.vector(apply(diff, c(2, 3), mean))
    grid$MSE <- as.vector(apply(diff^2, c(2, 3), mean))
    result[[length(result) + 1L]] <- grid
  }
  article_bind(result)
}

article_evaluate <- function(fit, datasets, meta, cfg, reused = FALSE) {
  result <- meta
  result$Converged <- fit$converged
  result$ConvergenceReason <- fit$convergence_reason
  result$Eligible <- !cfg$require_convergence || fit$converged
  result$Alpha <- fit$alpha
  result$LambdaMean <- if (all(is.na(fit$lambda_used))) NA_real_ else mean(fit$lambda_used)
  result$LambdaMin <- if (all(is.na(fit$lambda_used))) NA_real_ else min(fit$lambda_used)
  result$LambdaMax <- if (all(is.na(fit$lambda_used))) NA_real_ else max(fit$lambda_used)
  result$FitTimeSec <- if (reused) 0 else fit$runtime$total_sec
  result$SourceFitTimeSec <- fit$runtime$total_sec
  result$ReusedGlobalG1 <- reused
  result$Warnings <- paste(fit$diagnostics$warnings, collapse = " | ")
  scalar_criteria <- Filter(function(z) is.atomic(z) && length(z) == 1L, fit$criteria)
  for (nm in names(scalar_criteria)) result[[nm]] <- scalar_criteria[[nm]]
  result$EDFMethod <- paste(fit$criteria$edf_method, collapse = ";")
  tables <- list(overall = list(), class = list(), region = list(), region_class = list(),
    parameters_region = list(), parameters_observation = list(), memberships = list(),
    selection = list(), coefficients_cluster = list(), parameters_cluster = list())
  for (dataset_name in names(datasets)) {
    dat <- datasets[[dataset_name]]
    if (!nrow(dat$X)) next
    predict_args <- list(object = fit, newx = dat$X, new_unit_id = dat$unit, new_coords = dat$coords)
    if (fit$model == "fixed_clusters") predict_args$cluster <- dat$true_cluster
    prob <- if (dataset_name == "Train") fit$prob_train else do.call(stats::predict, predict_args)
    pp <- do.call(stats::predict, c(predict_args, list(type = "membership")))
    predicted <- fit$class_levels[max.col(prob, ties.method = "first")]
    beta_hat <- article_local_coef(fit, pp)
    residual <- beta_hat - dat$beta
    active <- seq_len(meta$NActive); inactive <- setdiff(seq_len(meta$P), active)
    overall <- scmr::calc_overall_metrics(dat$y, predicted, prob, fit$class_levels)
    overall$ProbabilityMSE <- mean((prob - dat$prob)^2)
    overall$BetaMSE <- mean(residual^2)
    overall$BetaMSEIntercept <- mean(residual[, , 1, drop = FALSE]^2)
    overall$BetaMSEActive <- mean(residual[, , active + 1L, drop = FALSE]^2)
    overall$BetaMSEInactive <- if (length(inactive)) mean(residual[, , inactive + 1L, drop = FALSE]^2) else NA_real_
    overall$ARI <- article_ari(dat$true_cluster, max.col(pp, ties.method = "first"))
    overall$Dataset <- dataset_name
    tables$overall[[dataset_name]] <- article_meta(overall, meta)
    for (nm in setdiff(names(overall), "Dataset")) result[[paste(dataset_name, nm, sep = "_")]] <- overall[[nm]]
    cls <- scmr::calc_class_metrics(dat$y, predicted, fit$class_levels)
    cls$Dataset <- dataset_name
    tables$class[[dataset_name]] <- article_meta(cls, meta)
    for (rg in sort(unique(dat$region))) {
      ii <- which(dat$region == rg)
      reg <- scmr::calc_overall_metrics(dat$y[ii], predicted[ii], prob[ii, , drop = FALSE], fit$class_levels)
      reg$ProbabilityMSE <- mean((prob[ii, , drop = FALSE] - dat$prob[ii, , drop = FALSE])^2)
      reg$BetaMSE <- mean(residual[ii, , , drop = FALSE]^2)
      reg$Region <- rg; reg$Dataset <- dataset_name
      tables$region[[length(tables$region) + 1L]] <- article_meta(reg, meta)
      rc <- scmr::calc_class_metrics(dat$y[ii], predicted[ii], fit$class_levels)
      rc$Region <- rg; rc$Dataset <- dataset_name
      tables$region_class[[length(tables$region_class) + 1L]] <- article_meta(rc, meta)
    }
    pars <- article_parameters(dat$beta, beta_hat, dat$region, active)
    pars$Dataset <- dataset_name
    tables$parameters_region[[dataset_name]] <- article_meta(pars, meta)
    if (cfg$export_per_observation) {
      pars <- article_parameters(dat$beta, beta_hat, dat$region, active, TRUE)
      pars$Unit <- dat$unit[as.integer(pars$Group)]; pars$Dataset <- dataset_name
      tables$parameters_observation[[dataset_name]] <- article_meta(pars, meta)
    }
    if (cfg$export_membership) {
      mem <- expand.grid(Row = seq_len(nrow(pp)), Cluster = seq_len(ncol(pp)))
      mem$Unit <- dat$unit[mem$Row]; mem$Probability <- as.vector(pp)
      mem$Dataset <- dataset_name
      tables$memberships[[dataset_name]] <- article_meta(mem, meta)
    }
    if (dataset_name == "Train") {
      pars <- article_parameters(dat$beta, beta_hat, paste0("G", fit$row_group), active)
      pars$Dataset <- dataset_name
      tables$parameters_cluster[[dataset_name]] <- article_meta(pars, meta)
      centered <- stats::coef(fit)
      native <- stats::coef(fit, parameterization = "native")
      for (g in seq_len(fit$G)) {
        ii <- which(fit$row_group == g)
        truth_mask <- apply(abs(dat$beta[ii, , -1, drop = FALSE]) > cfg$fit_control$selection_tol, c(2, 3), any)
        est_mask <- matrix(abs(centered[g, , -1, drop = FALSE]) > cfg$fit_control$selection_tol, nrow = length(fit$class_levels))
        pred_native <- matrix(abs(native[g, , -1, drop = FALSE]) > cfg$fit_control$selection_tol, nrow = length(fit$class_levels))
        sel <- rbind(article_support(colSums(pred_native) > 0, colSums(truth_mask) > 0),
                     article_support(est_mask, truth_mask))
        sel$Level <- c("predictor_native", "coefficient_sum_to_zero")
        sel$Cluster <- g
        tables$selection[[g]] <- article_meta(sel, meta)
        cf <- expand.grid(Class = fit$class_levels, Term = dimnames(native)[[3]], stringsAsFactors = FALSE)
        cf$Cluster <- g; cf$Native <- as.vector(native[g, , ])
        cf$SumToZero <- as.vector(centered[g, , ])
        cf$NativeSelected <- abs(cf$Native) > cfg$fit_control$selection_tol
        tables$coefficients_cluster[[g]] <- article_meta(cf, meta)
      }
    }
  }
  tables <- lapply(tables, article_bind)
  for (nm in c("tuning", "iterations", "solver")) tables[[paste0("diagnostics_", nm)]] <- article_meta(fit$diagnostics[[nm]], meta)
  tables$cluster_components <- article_meta(fit$cluster_components, meta)
  list(result = result, tables = tables, fit = if (cfg$save_fits) fit else NULL)
}
