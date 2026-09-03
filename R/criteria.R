# The Potts objective is reported separately from the IC likelihood.
spatial_bonus <- function(w, groups) {
  if (is.null(w)) return(0)
  # General sparse storage exposes both triangles consistently.
  entries <- Matrix::summary(methods::as(w, "generalMatrix"))
  sum(entries$x[groups[entries$i] == groups[entries$j]]) / 2
}

scmr_fit_criteria <- function(fits, x, y, row_group, unit_id, control,
                              w = NULL, group_unit = NULL) {
  g_total <- length(fits)
  n <- nrow(x)
  assigned <- numeric(n)
  parts <- vector("list", g_total)
  for (g in seq_len(g_total)) {
    idx <- which(row_group == g)
    if (!length(idx)) stop("Cannot evaluate an empty cluster.", call. = FALSE)
    lp <- true_class_logp(engine_predict(fits[[g]], x[idx, , drop = FALSE]),
                         y[idx], fits[[g]]$classes)
    assigned[idx] <- lp
    edf <- engine_effective_df(fits[[g]], x[idx, , drop = FALSE], control$selection_tol)
    parts[[g]] <- data.frame(
      Cluster = g, NRowsCluster = length(idx), NUnitsCluster = length(unique(unit_id[idx])),
      ClusterLogLik = sum(lp), ClusterNeg2LogLik = -2 * sum(lp),
      DFOriginal_g = edf$nominal, DFEffective_g = edf$df,
      DFNominalActive_g = edf$df_nominal_active,
      ActivePredictors_g = edf$active_predictors,
      ActiveCoefficients_g = edf$active_coefficients,
      LambdaUsed_g = fits[[g]]$lambda, EDFMethod = edf$method,
      stringsAsFactors = FALSE)
  }
  parts <- do.call(rbind, parts)
  m <- parts$NUnitsCluster
  a <- (mean(m) / m)^control$criterion_balance_gamma
  weight <- a * n / sum(a * parts$NRowsCluster)
  parts$CBBIC_n_c <- parts$NRowsCluster
  parts$CBBIC_m_c <- m
  parts$CBBIC_m_bar <- mean(m)
  parts$CBBIC_a_c <- a
  parts$CBBIC_w_c <- weight
  parts$CBBIC_weighted_N_c <- weight * parts$NRowsCluster
  parts$CBBIC_WeightedNeg2LogLik_c <- weight * parts$ClusterNeg2LogLik
  parts$CBBIC_UnitPenalty_c <- log(pmax(m, 2)) * parts$DFEffective_g
  parts$CBBIC_Component_c <- parts$CBBIC_WeightedNeg2LogLik_c + parts$CBBIC_UnitPenalty_c
  parts$CBAIC_UnitPenalty_c <- 2 * parts$DFEffective_g
  parts$CBAIC_Component_c <- parts$CBBIC_WeightedNeg2LogLik_c + parts$CBAIC_UnitPenalty_c
  ll <- sum(assigned)
  df_orig <- sum(parts$DFOriginal_g)
  df_eff <- sum(parts$DFEffective_g)
  cb_fit <- sum(parts$CBBIC_WeightedNeg2LogLik_c)
  cb_pen <- sum(parts$CBBIC_UnitPenalty_c)
  bonus <- spatial_bonus(w, group_unit)
  criteria <- list(
    logLik = ll, neg2_logLik = -2 * ll, spatial_bonus = bonus,
    objective_value = ll + control$phi * bonus,
    neg2_objective = -2 * (ll + control$phi * bonus),
    n_rows = n, n_units = length(unique(unit_id)), n_predictors = ncol(x),
    df_original_total = df_orig, df_effective_total = df_eff,
    df_nominal_active_total = sum(parts$DFNominalActive_g),
    edf_method = unique(parts$EDFMethod),
    edf_is_approximation = any(vapply(fits, function(z) z$kind == "glmnet", logical(1))),
    active_predictors_total = sum(parts$ActivePredictors_g),
    penalty_scr_bic_original = log(n) * df_orig,
    penalty_scr_bic_effective = log(n) * df_eff,
    penalty_scr_aic_effective = 2 * df_eff,
    CriterionSCR_BIC_original = -2 * ll + log(n) * df_orig,
    CriterionSCR_BIC_effective = -2 * ll + log(n) * df_eff,
    CriterionSCR_AIC_effective = -2 * ll + 2 * df_eff,
    CBBIC_Fit = cb_fit, CBBIC_Penalty = cb_pen,
    CBBIC_WeightedN_Check = sum(parts$CBBIC_weighted_N_c),
    CBBIC_m_bar = mean(m), CBBIC_EmptyClusterCount = 0L,
    CBBIC_InfeasibleClusterCount = 0L, CBBIC_FeasibilityPenalty = 0,
    CriterionCB_BIC = cb_fit + cb_pen,
    CriterionCB_AIC = cb_fit + 2 * df_eff,
    cluster_components = parts)
  criteria$CriterionV3_CB_uBIC <- criteria$CriterionCB_BIC
  criteria
}
