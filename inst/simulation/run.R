article_dgp_diagnostics <- function(dat, meta) {
  corr <- stats::cor(dat$observed$X)
  active <- seq_len(meta$NActive); inactive <- setdiff(seq_len(meta$P), active)
  pair_mean <- function(idx) {
    if (length(idx) < 2) return(NA_real_)
    z <- corr[idx, idx, drop = FALSE]; mean(z[upper.tri(z)])
  }
  counts <- data.frame(Region = names(dat$observed_region_counts),
    ObservedCount = as.integer(dat$observed_region_counts), NewCount = as.integer(dat$new_region_counts),
    TargetProportion = as.numeric(dat$target_region_probs),
    ActiveCorrelation = pair_mean(active), InactiveCorrelation = pair_mean(inactive),
    CrossCorrelation = if (length(inactive)) mean(corr[active, inactive, drop = FALSE]) else NA_real_)
  counts$ObservedProportion <- counts$ObservedCount / sum(counts$ObservedCount)
  counts$NewProportion <- if (sum(counts$NewCount)) counts$NewCount / sum(counts$NewCount) else NA_real_
  article_meta(counts, meta)
}

run_article_simulation <- function(cfg = article_config()) {
  cfg <- do.call(article_config, cfg)
  out <- article_open(cfg)
  grid <- article_grid(cfg)
  grid$DatasetID <- vapply(seq_len(nrow(grid)), function(i) article_hash(grid[i, , drop = FALSE]), character(1))
  expected <- article_bind(lapply(seq_len(nrow(grid)), function(i) {
    models <- article_models(grid$Scenario[i], cfg)
    z <- article_meta(models, grid[i, , drop = FALSE])
    z$RunID <- vapply(seq_len(nrow(z)), function(j) article_hash(z[j, , drop = FALSE]), character(1))
    z
  }))
  article_csv(expected, file.path(out, "expected_runs.csv"))
  message("Results: ", out)
  for (i in seq_len(nrow(grid))) {
    meta <- grid[i, , drop = FALSE]
    models <- expected[expected$DatasetID == meta$DatasetID, , drop = FALSE]
    completed <- file.exists(file.path(out, "runs", paste0(models$RunID, ".rds")))
    if (all(completed)) { message("Already complete: ", meta$Scenario, " repeat ", meta$Repeat); next }
    message("Generating ", meta$Scenario, " repeat ", meta$Repeat)
    data_path <- file.path(out, "data", paste0(meta$DatasetID, ".rds"))
    data_bundle <- tryCatch({
      if (file.exists(data_path)) readRDS(data_path) else {
        dat <- scmr::simulate_scmr_data(n_obs = cfg$n_obs, n_new = cfg$n_new, p = meta$P,
          active = seq_len(meta$NActive), eta = meta$Eta, scenario = meta$Scenario,
          class_levels = cfg$class_levels, seed = meta$Seed, class_balance = meta$ClassBalance,
          heterogeneity_strength = meta$Heterogeneity, control = cfg$dgp_control)
        test_idx <- scmr::make_stratified_split(dat$observed$y, dat$observed$region, cfg$test_prop, meta$Seed)
        train_idx <- setdiff(seq_len(cfg$n_obs), test_idx)
        z <- list(datasets = list(Train = article_subset(dat$observed, train_idx),
          Test = article_subset(dat$observed, test_idx), New = dat$new),
          diagnostics = article_dgp_diagnostics(dat, meta))
        article_save(z, data_path)
        z
      }
    }, error = function(e) e)
    if (inherits(data_bundle, "error")) {
      for (j in which(!completed)) article_save(article_meta(data.frame(Stage = "data_generation",
        Message = conditionMessage(data_bundle)), models[j, , drop = FALSE]),
        file.path(out, "errors", paste0(models$RunID[j], ".rds")))
      next
    }
    datasets <- data_bundle$datasets; train <- datasets$Train
    globals <- list()
    # Local coefficient features are computed once per dataset and reused across G.
    feature_cache <- NULL
    local_features <- function() {
      if (is.null(feature_cache)) {
        mc <- article_run_control(cfg, "SCMR-EN-MS")
        feature_cache <<- scmr::scmr_local_coefficients(train$X, train$y, train$unit, train$coords,
          k = mc$coef_init_k, alpha = mc$coef_init_alpha, lambda = mc$coef_init_lambda,
          anchors = mc$coef_init_anchors, type_multinomial = mc$type_multinomial, seed = meta$Seed)
      }
      feature_cache
    }
    for (j in seq_len(nrow(models))) {
      run <- models[j, , drop = FALSE]
      path <- file.path(out, "runs", paste0(run$RunID, ".rds"))
      if (file.exists(path)) next
      message("  ", run$Model, " G=", run$G)
      error_path <- file.path(out, "errors", paste0(run$RunID, ".rds"))
      tryCatch({
        can_reuse <- length(cfg$fit_control$min_per_class) == 1L &&
          length(unique(train$unit)) >= cfg$fit_control$min_units &&
          all(table(train$y) >= cfg$fit_control$min_per_class)
        reused <- can_reuse && run$Mode == "scmr" && run$G == 1L && !is.null(globals[[run$Penalty]])
        fit <- if (reused) globals[[run$Penalty]] else scmr::fit_scmr(
          x = train$X, y = train$y, model = run$Mode, penalty = run$Penalty, G = run$G,
          unit_id = train$unit, coords = train$coords,
          cluster = if (run$Mode == "fixed_clusters") train$true_cluster else NULL,
          control = article_run_control(cfg, run$Model), nfolds = cfg$nfolds, seed = meta$Seed,
          init_features = if (run$Model == "SCMR-EN-MS") local_features() else NULL)
        if (run$Mode == "global") globals[[run$Penalty]] <- fit
        bundle <- article_evaluate(fit, datasets, run, cfg, reused)
        article_save(bundle, path)
        if (file.exists(error_path)) unlink(error_path)
      }, error = function(e) {
        article_save(article_meta(data.frame(Stage = "fit_or_evaluate", Message = conditionMessage(e)), run), error_path)
        message("    Failed: ", conditionMessage(e))
      })
    }
    # Derived CSV reports can always be rebuilt from committed model bundles.
    article_summarize(out, expected, cfg)
  }
  dgp <- lapply(list.files(file.path(out, "data"), pattern = "\\.rds$", full.names = TRUE), function(p) readRDS(p)$diagnostics)
  article_csv(article_bind(dgp), file.path(out, "dgp_diagnostics.csv"))
  summary <- article_summarize(out, expected, cfg)
  message("Completed fits: ", nrow(summary$results), "/", nrow(expected),
    "; errors: ", nrow(summary$errors), "; not converged: ", sum(summary$coverage$Status == "not_converged"))
  invisible(c(list(output_dir = out), summary))
}
