# Public reproducibility script for the SCMR simulation study.
# This script uses only synthetic data and can be committed to GitHub.
# Increase n_repeats, n_obs, n_new, p_grid, and k_grid for the full study.

library(scmr)
library(Matrix)

output_dir <- file.path(getwd(), "simulation_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Small defaults for quick public verification.
# Suggested full study settings:
# n_obs = 3000, n_new = 300, n_repeats = 50,
# scenarios = c("global", "clustered_balanced", "clustered_imbalanced"),
# p_grid = c(25), eta_grid = c(0.2), k_grid = 1:10.
n_obs <- 300
n_new <- 60
n_repeats <- 2
scenarios <- c("global", "clustered_balanced")
eta_grid <- c(0.2)
p_grid <- c(10)
k_grid <- 1:3
active_count <- 5
penalties <- c("elastic_net")

ctrl <- scmr_control(
  max_iter = 3,
  init_tries = 20,
  init_method = "kmeans",
  min_units = 5,
  min_per_class = 2,
  alpha_grid = c(0.6, 0.8),
  holdout_nlambda = 30,
  tune_lambda = TRUE,
  tune_alpha_after_convergence = TRUE
)

results <- list()
row_id <- 1L

for (rr in seq_len(n_repeats)) {
  for (scenario in scenarios) {
    for (eta in eta_grid) {
      for (p in p_grid) {
        for (penalty in penalties) {
          seed <- 1000 + rr * 10000 + round(eta * 1000) + p
          message("Repeat ", rr, " | scenario = ", scenario, " | penalty = ", penalty, " | p = ", p)

          dat <- simulate_scmr_data(
            n_obs = n_obs,
            n_new = n_new,
            p = p,
            active = seq_len(active_count),
            eta = eta,
            scenario = scenario,
            seed = seed
          )

          obs <- dat$observed
          test_idx <- make_stratified_split(obs$y, obs$region, test_prop = 0.30, seed = seed)
          train_idx <- setdiff(seq_along(obs$y), test_idx)

          X_train <- Matrix(obs$X[train_idx, , drop = FALSE], sparse = TRUE)
          X_test <- Matrix(obs$X[test_idx, , drop = FALSE], sparse = TRUE)
          y_train <- factor(obs$y[train_idx], levels = levels(obs$y))
          y_test <- factor(obs$y[test_idx], levels = levels(obs$y))

          # Global benchmark. Choose model = "global".
          fit_g <- fit_scmr(
            X_train,
            y_train,
            model = "global",
            penalty = penalty,
            control = ctrl,
            seed = seed
          )
          prob_g <- predict(fit_g, X_test, type = "prob")
          pred_g <- factor(colnames(prob_g)[max.col(prob_g)], levels = levels(y_train))
          met_g <- calc_overall_metrics(y_test, pred_g, prob_g, levels(y_train))
          results[[row_id]] <- cbind(
            data.frame(
              Repeat = rr,
              Scenario = scenario,
              Eta = eta,
              P = p,
              Penalty = penalty,
              Model = "global",
              FixedK = 1,
              Alpha = fit_g$alpha,
              Lambda = fit_g$lambda,
              CriterionSCR_BIC_original = fit_g$criteria$CriterionSCR_BIC_original,
              CriterionV3_CB_uBIC = fit_g$criteria$CriterionV3_CB_uBIC
            ),
            met_g
          )
          row_id <- row_id + 1L

          # Spatially clustered SCMR. Choose model = "scmr" and set G.
          for (G in k_grid) {
            fit_s <- tryCatch({
              fit_scmr(
                X_train,
                y_train,
                model = "scmr",
                penalty = penalty,
                G = G,
                unit_id = obs$unit[train_idx],
                coords = obs$coords[train_idx, , drop = FALSE],
                control = ctrl,
                seed = seed + G
              )
            }, error = function(e) {
              message("SCMR failed for G = ", G, ": ", conditionMessage(e))
              NULL
            })

            if (is.null(fit_s)) next

            prob_s <- predict(
              fit_s,
              X_test,
              new_unit_id = obs$unit[test_idx],
              new_coords = obs$coords[test_idx, , drop = FALSE],
              type = "prob"
            )
            pred_s <- factor(colnames(prob_s)[max.col(prob_s)], levels = levels(y_train))
            met_s <- calc_overall_metrics(y_test, pred_s, prob_s, levels(y_train))
            results[[row_id]] <- cbind(
              data.frame(
                Repeat = rr,
                Scenario = scenario,
                Eta = eta,
                P = p,
                Penalty = penalty,
                Model = "scmr",
                FixedK = G,
                Alpha = fit_s$alpha_final,
                Lambda = mean(fit_s$lambda_used),
                CriterionSCR_BIC_original = fit_s$criteria$CriterionSCR_BIC_original,
                CriterionV3_CB_uBIC = fit_s$criteria$CriterionV3_CB_uBIC,
                Converged = fit_s$converged,
                ConvergenceReason = fit_s$convergence_reason
              ),
              met_s
            )
            row_id <- row_id + 1L
          }
        }
      }
    }
  }
}

results_df <- do.call(rbind, results)
write.csv(results_df, file.path(output_dir, "scmr_simulation_results_by_k.csv"), row.names = FALSE)

best_cb <- do.call(rbind, lapply(
  split(results_df, list(results_df$Repeat, results_df$Scenario, results_df$Penalty, results_df$Model), drop = TRUE),
  function(z) z[which.min(z$CriterionV3_CB_uBIC), , drop = FALSE]
))
write.csv(best_cb, file.path(output_dir, "scmr_best_by_cb_bic.csv"), row.names = FALSE)

message("Done. Results written to: ", output_dir)
