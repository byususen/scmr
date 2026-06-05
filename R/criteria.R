compute_criteria_from_row_loglik <- function(row_loglik_by_g, group_unit, unit_num,
                                             unit_id, y, x_ncol, W,
                                             active_per_group, G, lambda_used,
                                             class_levels, unit_total, unit_class,
                                             df_effective_group, alpha,
                                             phi = 1,
                                             criterion_balance_gamma = 1) {
  k <- length(class_levels)
  n_rows <- length(y)
  n_units <- length(unique(unit_id))
  group_unit <- as.integer(group_unit)
  row_group <- group_unit[unit_num]

  assigned_loglik <- row_loglik_by_g[cbind(seq_len(n_rows), row_group)]
  logLik <- sum(assigned_loglik)
  neg2_logLik <- -2 * logLik

  Wsum <- summary(W)
  if (G == 1) {
    spatial_bonus <- 0.5 * sum(Wsum$x)
  } else {
    same_cluster <- group_unit[Wsum$i] == group_unit[Wsum$j]
    spatial_bonus <- 0.5 * sum(Wsum$x[same_cluster])
  }
  objective_value <- logLik + phi * spatial_bonus

  df_original_group <- rep(k * (x_ncol + 1), G)
  df_original_total <- sum(df_original_group)
  df_effective_group <- as.numeric(df_effective_group)
  bad <- !is.finite(df_effective_group)
  if (any(bad)) df_effective_group[bad] <- max(k - 1, 1) * (active_per_group[bad] + 1)
  df_effective_total <- sum(df_effective_group)

  cluster_logLik <- numeric(G)
  cluster_neg2 <- numeric(G)
  cluster_n_rows <- integer(G)
  cluster_n_units <- integer(G)
  cluster_components <- vector("list", G)

  for (g in seq_len(G)) {
    idx_g_rows <- which(row_group == g)
    idx_g_units <- which(group_unit == g)
    cluster_n_rows[g] <- length(idx_g_rows)
    cluster_n_units[g] <- length(idx_g_units)
    cluster_logLik[g] <- if (length(idx_g_rows) > 0) sum(row_loglik_by_g[idx_g_rows, g]) else 0
    cluster_neg2[g] <- -2 * cluster_logLik[g]
    cluster_components[[g]] <- data.frame(
      Cluster = g,
      NRowsCluster = cluster_n_rows[g],
      NUnitsCluster = cluster_n_units[g],
      ClusterLogLik = cluster_logLik[g],
      ClusterNeg2LogLik = cluster_neg2[g],
      DFOriginal_g = df_original_group[g],
      DFEffective_g = df_effective_group[g],
      ActivePredictors_g = active_per_group[g],
      LambdaUsed_g = lambda_used[g],
      stringsAsFactors = FALSE
    )
  }
  cluster_components <- do.call(rbind, cluster_components)

  n_c <- cluster_n_rows
  m_c <- cluster_n_units
  m_bar <- mean(m_c[m_c > 0], na.rm = TRUE)
  if (!is.finite(m_bar) || is.na(m_bar)) m_bar <- 1

  a_c <- (m_bar / pmax(m_c, 1)) ^ criterion_balance_gamma
  a_c[!is.finite(a_c)] <- 1
  row_scale <- n_rows / sum(a_c * n_c, na.rm = TRUE)
  if (!is.finite(row_scale) || is.na(row_scale)) row_scale <- 1
  w_c <- row_scale * a_c

  cb_fit_c <- w_c * cluster_neg2
  cb_fit <- sum(cb_fit_c, na.rm = TRUE)
  cb_penalty_c <- log(pmax(m_c, 2)) * df_effective_group
  cb_penalty <- sum(cb_penalty_c, na.rm = TRUE)

  criterion_cb_bic <- cb_fit + cb_penalty
  criterion_scr_bic <- neg2_logLik + log(max(n_rows, 2)) * df_original_total

  cluster_components$CBuBIC_n_c <- n_c
  cluster_components$CBuBIC_m_c <- m_c
  cluster_components$CBuBIC_m_bar <- m_bar
  cluster_components$CBuBIC_a_c <- a_c
  cluster_components$CBuBIC_w_c <- w_c
  cluster_components$CBuBIC_weighted_N_c <- w_c * n_c
  cluster_components$CBuBIC_WeightedNeg2LogLik_c <- cb_fit_c
  cluster_components$CBuBIC_UnitPenalty_c <- cb_penalty_c
  cluster_components$CBuBIC_Component_c <- cb_fit_c + cb_penalty_c

  list(
    logLik = logLik,
    neg2_logLik = neg2_logLik,
    spatial_bonus = spatial_bonus,
    objective_value = objective_value,
    n_rows = n_rows,
    n_units = n_units,
    n_predictors = x_ncol,
    df_original_total = df_original_total,
    df_effective_total = df_effective_total,
    active_predictors_total = sum(active_per_group),
    CriterionSCR_BIC_original = criterion_scr_bic,
    CBuBIC_Fit = cb_fit,
    CBuBIC_Penalty = cb_penalty,
    CBuBIC_m_bar = m_bar,
    CBuBIC_WeightedN_Check = sum(w_c * n_c, na.rm = TRUE),
    CriterionV3_CB_uBIC = criterion_cb_bic,
    cluster_components = cluster_components
  )
}

compute_global_criteria <- function(fit_global, X_train, y_train, prob_train,
                                    alpha, s_rule, class_levels,
                                    standardize = TRUE) {
  logLik_global <- sum(true_class_logp(prob_train, y_train, class_levels))
  neg2_logLik_global <- -2 * logLik_global
  df_original_global <- length(class_levels) * (ncol(X_train) + 1)
  df_effective_global <- elastic_net_effective_df(
    fit_global, X_train, s = s_rule, alpha = alpha,
    class_levels = class_levels, standardize = standardize
  )
  list(
    logLik = logLik_global,
    neg2_logLik = neg2_logLik_global,
    n_rows = nrow(X_train),
    n_units = NA_real_,
    n_predictors = ncol(X_train),
    active_predictors_total = NA_real_,
    df_original_total = df_original_global,
    df_effective_total = df_effective_global,
    CriterionSCR_BIC_original = neg2_logLik_global + log(max(nrow(X_train), 2)) * df_original_global,
    CriterionV3_CB_uBIC = neg2_logLik_global + log(max(nrow(X_train), 2)) * df_effective_global,
    lambda_used = get_s_value(fit_global, s_rule)
  )
}
