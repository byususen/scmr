new_scmr_fit <- function(model, penalty, fits, data, groups, control, seed,
                         tuning = data.frame(), iterations = data.frame(),
                         membership_converged = TRUE, reason = "fit_complete",
                         w = NULL, bandwidth = NULL, fixed_levels = NULL) {
  unit_levels <- sort(unique(data$unit_id))
  unit_num <- match(data$unit_id, unit_levels)
  group_unit <- vapply(seq_along(unit_levels), function(u) {
    z <- unique(groups[unit_num == u])
    if (length(z) != 1L) stop("A spatial unit has conflicting memberships.", call. = FALSE)
    as.integer(z)
  }, integer(1))
  G <- length(fits)
  pp <- matrix(0, length(unit_levels), G,
               dimnames = list(unit_levels, paste0("G", seq_len(G))))
  pp[cbind(seq_along(group_unit), group_unit)] <- 1
  prob <- matrix(0, nrow(data$x), nlevels(data$y), dimnames = list(NULL, levels(data$y)))
  for (g in seq_len(G)) {
    idx <- which(groups == g)
    prob[idx, ] <- engine_predict(fits[[g]], data$x[idx, , drop = FALSE])
  }
  criteria <- scmr_fit_criteria(fits, data$x, data$y, groups, data$unit_id, control, w, group_unit)
  solver <- data.frame(Cluster = seq_len(G),
    Converged = vapply(fits, function(z) isTRUE(z$converged), logical(1)),
    Message = vapply(fits, function(z) paste(z$warnings, collapse = ";"), character(1)))
  converged <- membership_converged && all(solver$Converged)
  if (!all(solver$Converged)) reason <- paste(reason, "local_solver_not_converged", sep = ";")
  parts <- criteria$cluster_components
  out <- list(
    schema_version = 1L, model = model, penalty = penalty, G = G,
    class_levels = levels(data$y), x_colnames = colnames(data$x),
    fits = fits, fit = if (G == 1L) fits[[1]]$fit else NULL,
    alpha = fits[[1]]$alpha, alpha_final = fits[[1]]$alpha,
    lambda_used = vapply(fits, function(z) z$lambda, numeric(1)),
    lambda = if (G == 1L) fits[[1]]$lambda else NULL,
    prob_train = prob, prob_oof_train = NULL, y_train = data$y,
    row_group = groups, group_unit = stats::setNames(group_unit, unit_levels),
    PP_unit = pp, unit_levels = unit_levels, unit_id = data$unit_id,
    train_unit_coords = if (is.null(data$coords)) NULL else get_unit_coords(data$unit_id, data$coords),
    fixed_cluster_levels = fixed_levels, control = control, seed = seed,
    W = w, bandwidth_used = bandwidth, criteria = criteria, cluster_components = parts,
    active_per_group = parts$ActivePredictors_g, df_effective_group = parts$DFEffective_g,
    converged = converged, convergence_reason = reason, obj_trace = iterations,
    initial_tuning = tuning[tuning$Stage %in% "initial", , drop = FALSE],
    diagnostics = list(tuning = tuning, iterations = iterations, solver = solver,
                       membership_converged = membership_converged,
                       warnings = unique(unlist(lapply(fits, function(z) z$warnings)))),
    no_enet = penalty == "none")
  class(out) <- c(switch(model, global = "scmr_global", scmr = "scmr_spatial",
                         fixed_clusters = "scmr_fixed"), "scmr_fit")
  out
}

fit_global_core <- function(data, penalty, alpha_info, lambda, control, seed, nfolds, lambda_rule) {
  tuning <- data.frame()
  oof <- NULL
  tune_seconds <- 0
  if (penalty == "none" || !is.null(lambda)) {
    engine <- fit_local_engine(data$x, data$y, penalty, alpha_info$alpha_default,
                               if (penalty == "none") NA_real_ else lambda, control)
    rule <- if (penalty == "none") "none" else "fixed"
  } else {
    t0 <- proc.time()[["elapsed"]]
    foldid <- make_stratified_foldid(data$y, nfolds, scmr_seed(seed, 7100))
    candidates <- list()
    trace <- list()
    for (a in alpha_info$alpha_grid) {
      warnings <- character()
      raw <- tryCatch(withCallingHandlers(glmnet::cv.glmnet(
        x = glmnet_design(data$x), y = data$y, family = "multinomial", alpha = a,
        foldid = foldid, type.measure = "deviance", standardize = control$standardize,
        type.multinomial = control$type_multinomial, maxit = control$glmnet_maxit, keep = TRUE),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }), error = function(e) e)
      if (inherits(raw, "error")) {
        trace[[length(trace) + 1L]] <- data.frame(Stage = "global_cv", Cluster = 1L, Iter = 0L,
          Alpha = a, Lambda = NA_real_, Loss = Inf, Selected = FALSE,
          Status = "error", Message = conditionMessage(raw))
        next
      }
      idx <- get_s_index(raw, lambda_rule)
      score <- raw$cvm[idx]
      if (!is.null(raw$glmnet.fit$jerr) && raw$glmnet.fit$jerr != 0L) score <- Inf
      candidates[[length(candidates) + 1L]] <- list(fit = raw, alpha = a,
        lambda = get_s_value(raw, lambda_rule), score = score, warnings = warnings)
      trace[[length(trace) + 1L]] <- data.frame(Stage = "global_cv", Cluster = 1L, Iter = 0L,
        Alpha = a, Lambda = raw$lambda, Loss = raw$cvm,
        Selected = seq_along(raw$lambda) == idx,
        Status = if (is.finite(score)) "ok" else "solver_failure", Message = paste(warnings, collapse = ";"))
    }
    scores <- vapply(candidates, function(z) z$score, numeric(1))
    if (!any(is.finite(scores))) stop("Global cross-validation failed for every alpha.", call. = FALSE)
    best <- candidates[[which.min(scores)]]
    engine <- engine_from_glmnet(best$fit, data$x, levels(data$y), best$alpha, best$lambda, control, best$warnings)
    oof <- tryCatch(extract_cv_preval_prob_matrix(best$fit, levels(data$y), lambda_rule, nrow(data$x)),
                    error = function(e) NULL)
    tuning <- bind_diagnostics(trace)
    tune_seconds <- proc.time()[["elapsed"]] - t0
    rule <- lambda_rule
  }
  out <- new_scmr_fit("global", penalty, list(engine), data, rep(1L, nrow(data$x)),
                      control, seed, tuning, reason = "global_fit")
  out$prob_oof_train <- oof
  out$lambda_rule <- rule
  out$runtime_tuning <- tune_seconds
  out
}

#' Fit a global multinomial benchmark
#' @inheritParams fit_scmr
#' @param standardize Standardize predictors within each numerical fit.
#' @param type_multinomial Either grouped or ungrouped multinomial penalization.
#' @return A fitted object inheriting from scmr_fit and scmr_global.
#' @export
fit_global_scmr <- function(x, y, penalty = c("elastic_net", "lasso", "ridge", "none"),
                            alpha = NULL, alpha_grid = NULL, lambda = NULL,
                            nfolds = 5, lambda_rule = c("lambda.1se", "lambda.min"),
                            seed = 123, standardize = TRUE, type_multinomial = "grouped",
                            unit_id = NULL, control = NULL) {
  if (is.null(control)) control <- scmr_control(standardize = standardize, type_multinomial = type_multinomial)
  fit_scmr(x, y, model = "global", penalty = match.arg(penalty), alpha = alpha,
           alpha_grid = alpha_grid, lambda = lambda, control = control, nfolds = nfolds,
           lambda_rule = match.arg(lambda_rule), seed = seed, unit_id = unit_id)
}
