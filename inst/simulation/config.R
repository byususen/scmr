# Article settings live here; numerical algorithms live in the scmr package.
article_config <- function(profile = c("smoke", "pilot", "main"), output_dir = "scmr-results", ...) {
  profile <- match.arg(profile)
  fit <- scmr::scmr_control(type_multinomial = "ungrouped", alpha_grid = c(.3, .5, .7, .9))
  cfg <- list(profile = profile, output_dir = output_dir, n_obs = 3000L, n_new = 300L,
    n_repeats = 10L, test_prop = .30, n_active = 5L, n_inactive = 5L,
    scenarios = c("global", "clustered_balanced", "clustered_imbalanced", "smooth"),
    eta = .2, class_balance = "balanced", heterogeneity_strength = 1,
    class_levels = c("C1", "C2", "C3"), G = 1:10, nfolds = 5L, seed = 1000L,
    fit_control = fit, dgp_control = scmr::scmr_simulation_control(),
    resume = TRUE, save_fits = FALSE, export_per_observation = FALSE,
    export_membership = FALSE, require_convergence = TRUE,
    # Optional extra models: "TwoStage-EN" (k-means partition, no membership
    # updates) and "SCMR-EN-MS" (monotone sum-scale algorithm with multi-start
    # coefficient-space initialization). Both are evaluated over the full G grid.
    extra_models = character(), ms_control = list(lambda_scale = "sum", n_starts = 3L,
      tiny_movement_max_units = 0, tiny_movement_rate_tol = 0, tiny_movement_revert = FALSE,
      max_iter = 30L))
  if (profile == "pilot") {
    cfg$n_obs <- 600L; cfg$n_new <- 60L; cfg$n_repeats <- 2L; cfg$G <- c(1L, 3L, 6L)
  }
  if (profile == "smoke") {
    cfg$n_obs <- 360L; cfg$n_new <- 36L; cfg$n_repeats <- 1L
    cfg$n_active <- 2L; cfg$n_inactive <- 2L; cfg$G <- c(1L, 2L)
    cfg$fit_control$alpha_grid <- c(.3, .7)
    cfg$fit_control$holdout_nlambda <- 15L
    cfg$fit_control$min_units <- 3L
    cfg$fit_control$min_per_class <- 2L
    cfg$fit_control$max_iter <- 4L
  }
  extra <- list(...)
  if (length(extra) && (is.null(names(extra)) || any(!names(extra) %in% names(cfg)))) stop("Unknown or unnamed study setting.")
  cfg[names(extra)] <- extra
  if (length(cfg$n_obs) != 1L || cfg$n_obs < 3 || cfg$n_new < 0 || cfg$n_repeats < 1 ||
      cfg$test_prop <= 0 || cfg$test_prop >= 1 || any(cfg$n_active < 1) || any(cfg$n_inactive < 0) ||
      any(cfg$G < 1) || any(cfg$eta <= 0) || any(cfg$heterogeneity_strength < 0)) stop("Invalid study dimensions.")
  for (nm in c("n_obs", "n_new", "n_repeats", "n_active", "n_inactive", "G", "seed")) {
    if (any(!is.finite(cfg[[nm]])) || any(cfg[[nm]] != floor(cfg[[nm]]))) stop(nm, " must contain integers.")
  }
  if (any(!cfg$scenarios %in% c("global", "clustered_balanced", "clustered_imbalanced", "smooth",
                                 "clustered_irregular")) ||
      any(!cfg$extra_models %in% c("TwoStage-EN", "SCMR-EN-MS")) ||
      any(!cfg$class_balance %in% c("balanced", "imbalanced"))) stop("Unknown scenario or class balance.")
  if (anyDuplicated(cfg$G)) stop("G must contain distinct candidates.")
  cfg$fit_control <- do.call(scmr::scmr_control, cfg$fit_control)
  cfg$dgp_control <- do.call(scmr::scmr_simulation_control, cfg$dgp_control)
  cfg
}

article_grid <- function(cfg) {
  grid <- expand.grid(Repeat = seq_len(cfg$n_repeats), Scenario = cfg$scenarios,
    Eta = cfg$eta, ClassBalance = cfg$class_balance, NActive = cfg$n_active,
    NInactive = cfg$n_inactive, Heterogeneity = cfg$heterogeneity_strength,
    stringsAsFactors = FALSE)
  grid$Heterogeneity[!grepl("^clustered", grid$Scenario)] <- 1
  grid <- unique(grid)
  grid$P <- grid$NActive + grid$NInactive
  grid$Seed <- with(grid, cfg$seed + Repeat * 1000000 + round(Eta * 1000) + P * 100 +
    NActive * 10000 + round(Heterogeneity * 100000) +
    ifelse(Scenario == "smooth", 5000000, ifelse(Scenario == "clustered_balanced", 7000000,
    ifelse(Scenario == "clustered_imbalanced", 9000000, ifelse(Scenario == "clustered_irregular", 11000000, 0)))) +
    ifelse(ClassBalance == "imbalanced", 13000000, 0))
  if (any(grid$Seed > .Machine$integer.max - 10000)) stop("Study seed grid exceeds R's seed range.")
  # Detect, rather than silently permit, collisions in the source script's seed formula.
  if (anyDuplicated(grid$Seed)) stop("Study settings produce duplicate seeds; use separate configurations.")
  rownames(grid) <- NULL
  grid
}

article_true_G <- function(scenario) {
  switch(scenario, global = 1L, clustered_balanced = 6L, clustered_imbalanced = 6L,
         clustered_irregular = 3L, smooth = NA_integer_)
}

article_models <- function(scenario, cfg) {
  models <- data.frame(Model = c("Global", "Global-EN", rep("SCMR-EN", length(cfg$G))),
    Mode = c("global", "global", rep("scmr", length(cfg$G))),
    Penalty = c("none", "elastic_net", rep("elastic_net", length(cfg$G))),
    G = c(1L, 1L, cfg$G), stringsAsFactors = FALSE)
  true_G <- article_true_G(scenario)
  if (scenario != "smooth") models <- rbind(models,
    data.frame(Model = "SCMR", Mode = "scmr", Penalty = "none", G = true_G))
  if (grepl("^clustered", scenario)) models <- rbind(models,
    data.frame(Model = c("TCMR", "TCMR-EN"), Mode = "fixed_clusters",
      Penalty = c("none", "elastic_net"), G = true_G))
  for (extra in cfg$extra_models) {
    g_extra <- setdiff(cfg$G, 1L)
    if (length(g_extra)) models <- rbind(models,
      data.frame(Model = extra, Mode = "scmr", Penalty = "elastic_net", G = g_extra))
  }
  models
}

# Fitting controls for each study model; extra models modify the article control.
article_run_control <- function(cfg, model) {
  ctrl <- cfg$fit_control
  if (identical(model, "TwoStage-EN")) ctrl$update_memberships <- FALSE
  if (identical(model, "SCMR-EN-MS")) ctrl[names(cfg$ms_control)] <- cfg$ms_control
  do.call(scmr::scmr_control, ctrl)
}
