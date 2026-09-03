#' Control the synthetic data generator
#'
#' @param ... Named overrides to the returned defaults. Unknown names are errors.
#' @return A named list of domain, covariance, coefficient and region settings.
#' @details Region probabilities specify location counts, independently of the
#' response class balance. The six region labels and rectangular partition must
#' remain consistent. Gaussian fields are generated jointly at observed and new
#' locations; new locations assess spatial interpolation.
#' @export
scmr_simulation_control <- function(...) {
  defaults <- simulation_defaults()
  supplied <- list(...)
  if (length(supplied) && (is.null(names(supplied)) || any(names(supplied) == "") ||
                          anyDuplicated(names(supplied)))) {
    stop("Simulation settings must have unique names.", call. = FALSE)
  }
  unknown <- setdiff(names(supplied), names(defaults))
  if (length(unknown)) stop("Unknown simulation setting: ", paste(unknown, collapse = ", "), call. = FALSE)
  defaults[names(supplied)] <- supplied
  for (nm in setdiff(names(defaults), c("region_levels6", "smooth_center_gp"))) {
    z <- defaults[[nm]]
    if (!is.numeric(z) || !length(z) || any(!is.finite(z))) {
      stop(nm, " must contain finite numbers.", call. = FALSE)
    }
  }
  for (nm in c("domain_s1_range", "domain_s2_range", "class_imbalance_shift_range_general",
               "beta_cluster_g1_values", "smooth_psi_intercept")) {
    if (length(defaults[[nm]]) != 2L) stop(nm, " must have length two.", call. = FALSE)
  }
  for (nm in c("domain_s1_range", "domain_s2_range")) {
    if (diff(defaults[[nm]]) <= 0) stop(nm, " must be increasing.", call. = FALSE)
  }
  if (!identical(defaults$region_levels6, paste0("R", rep(1:2, each = 3), rep(1:3, 2)))) {
    stop("The generator uses region labels R11 through R23.", call. = FALSE)
  }
  br <- defaults$region_s2_breaks
  if (length(br) != 4L || any(diff(br) <= 0) ||
      !isTRUE(all.equal(br[c(1, 4)], defaults$domain_s2_range)) ||
      length(defaults$region_s1_split) != 1L ||
      defaults$region_s1_split <= defaults$domain_s1_range[1] ||
      defaults$region_s1_split >= defaults$domain_s1_range[2]) {
    stop("Region boundaries must partition the domain.", call. = FALSE)
  }
  for (nm in c("cluster_member_probs_balanced", "cluster_member_probs_imbalanced")) {
    z <- defaults[[nm]]
    if (length(z) != 6L || any(z <= 0)) stop(nm, " must contain six positive weights.", call. = FALSE)
    if (is.null(names(z))) names(z) <- defaults$region_levels6
    if (!setequal(names(z), defaults$region_levels6) || anyDuplicated(names(z))) {
      stop(nm, " has incompatible region names.", call. = FALSE)
    }
    defaults[[nm]] <- z[defaults$region_levels6] / sum(z)
  }
  scalar_positive <- c("gp_nugget", "predictor_gp_tau2", "predictor_eta_scale",
    "predictor_kernel_power", "smooth_psi_active_min", "smooth_psi_active_max")
  for (nm in scalar_positive) {
    if (length(defaults[[nm]]) != 1L || defaults[[nm]] <= 0) stop(nm, " must be positive.", call. = FALSE)
  }
  if (defaults$predictor_kernel_power > 2 || any(defaults$smooth_psi_intercept <= 0) ||
      defaults$smooth_psi_active_max < defaults$smooth_psi_active_min) {
    stop("Invalid Gaussian field ranges or kernel power.", call. = FALSE)
  }
  for (nm in c("predictor_rho_active", "predictor_rho_inactive")) {
    if (length(defaults[[nm]]) != 1L || abs(defaults[[nm]]) > 1) stop(nm, " must be in [-1, 1].", call. = FALSE)
  }
  if (length(defaults$beta_cluster_g2_values) != 3L || length(defaults$class_imbalance_shift_three) != 3L ||
      abs(sum(defaults$class_imbalance_shift_three)) > 1e-10) {
    stop("Three-class offsets must have length three and class shifts must sum to zero.", call. = FALSE)
  }
  if (!is.logical(defaults$smooth_center_gp) || length(defaults$smooth_center_gp) != 1L ||
      is.na(defaults$smooth_center_gp)) stop("smooth_center_gp must be TRUE or FALSE.", call. = FALSE)
  for (nm in c("domain_exclusion_s1_weight", "domain_exclusion_s2_weight", "domain_exclusion_radius",
               "smooth_beta_tau2_intercept", "smooth_beta_tau2_active")) {
    if (length(defaults[[nm]]) != 1L || defaults[[nm]] < 0) stop(nm, " must be nonnegative.", call. = FALSE)
  }
  bounds <- region_bounds_6(defaults)
  for (i in seq_len(nrow(bounds))) {
    corners <- as.matrix(expand.grid(c(bounds$s1_min[i], bounds$s1_max[i]), c(bounds$s2_min[i], bounds$s2_max[i])))
    if (!any(inside_irregular_domain(corners, defaults))) stop("The exclusion removes an entire sampling region.", call. = FALSE)
  }
  defaults
}

#' Simulate data for SCMR experiments
#'
#' @param n_obs Number of observed locations.
#' @param n_new Number of new locations, possibly zero.
#' @param p Number of predictors.
#' @param active Distinct active predictor indices, possibly empty.
#' @param eta Spatial range scale for predictors.
#' @param scenario Global, balanced or imbalanced six-region, or smooth coefficients.
#' @param class_levels Three distinct response labels. Model fitting itself allows
#'   any number of classes greater than one.
#' @param seed Random seed. The caller's random state is restored.
#' @param class_balance Whether to add an intercept shift favouring some classes.
#'   `"balanced"` means no extra shift, not exactly equal response counts.
#' @param heterogeneity_strength Nonnegative multiplier of clustered deviations
#'   from the global coefficients; zero produces global coefficients.
#' @param control Settings from [scmr_simulation_control()].
#' @return A list with `observed` and `new` data lists containing `X`, `y`, `coords`,
#'   `beta`, `prob`, `region`, `true_cluster`, and `unit`, plus scenario metadata.
#' @details Coefficients use a sum-to-zero class representation. Predictors are
#'   scaled over all generated locations as part of the data-generating process.
#'   Observed and new locations share the same Gaussian fields. Region counts in
#'   clustered scenarios follow the supplied proportions with integer rounding.
#' @export
simulate_scmr_data <- function(n_obs = 3000, n_new = 300, p = 25, active = 1:5,
                               eta = 0.2,
                               scenario = c("global", "clustered_balanced", "clustered_imbalanced", "smooth"),
                               class_levels = c("C1", "C2", "C3"), seed = 1000,
                               class_balance = c("balanced", "imbalanced"),
                               heterogeneity_strength = 1, control = scmr_simulation_control()) {
  scenario <- match.arg(scenario)
  class_balance <- match.arg(class_balance)
  control <- do.call(scmr_simulation_control, control)
  for (nm in c("n_obs", "n_new", "p")) {
    v <- get(nm)
    if (length(v) != 1L || !is.finite(v) || v != as.integer(v) || v < if (nm == "n_new") 0 else 1) {
      stop(nm, " must be a valid integer count.", call. = FALSE)
    }
  }
  if (n_obs + n_new < 2L) stop("Generate at least two locations.", call. = FALSE)
  if (anyNA(active) || anyDuplicated(active) || any(active != as.integer(active)) || any(active < 1 | active > p)) {
    stop("active must contain distinct predictor indices.", call. = FALSE)
  }
  if (length(class_levels) != 3L || anyNA(class_levels) || anyDuplicated(class_levels)) {
    stop("This data generator requires three distinct class labels.", call. = FALSE)
  }
  if (length(eta) != 1L || !is.finite(eta) || eta <= 0 || length(heterogeneity_strength) != 1L ||
      !is.finite(heterogeneity_strength) || heterogeneity_strength < 0) stop("Invalid eta or heterogeneity strength.", call. = FALSE)
  if (length(seed) != 1L || !is.finite(seed) || seed != floor(seed) || seed < 0 || seed > .Machine$integer.max - 20) {
    stop("seed must be a nonnegative integer below the integer limit minus 20.", call. = FALSE)
  }
  out <- with_scmr_seed(seed, make_simulation_dataset(
    n_obs, n_new, p, active, eta, scenario, class_balance, seed, class_levels,
    smooth_setup = list(psi_intercept = control$smooth_psi_intercept,
      beta_tau2_intercept = control$smooth_beta_tau2_intercept,
      beta_tau2_active = control$smooth_beta_tau2_active, center_gp = control$smooth_center_gp),
    predictor_eta_scale = control$predictor_eta_scale,
    predictor_kernel_power = control$predictor_kernel_power,
    predictor_rho_active = control$predictor_rho_active,
    predictor_rho_inactive = control$predictor_rho_inactive,
    cluster_member_probs = choose_cluster_member_probs(scenario, control),
    heterogeneity_strength = heterogeneity_strength, control = control))
  out$scenario <- scenario
  out$active <- active
  out$class_levels <- class_levels
  out$control <- control
  out
}
