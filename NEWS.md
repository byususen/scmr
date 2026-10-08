# scmr 0.4.0 (development; not yet checked with R CMD check)

- Monotone algorithm: `scmr_control(lambda_scale = "sum")` fixes one common
  sum-scale penalty strength kappa = mean(lambda) * n / G, standardizes predictors
  once globally and skips tuning inside the alternation. Every membership move
  and every refit then increases the penalized Potts objective, which is
  recorded per iteration as `PenalizedObjective`.
- Coefficient-space initialization (`init_method = "coefficient"`,
  `scmr_local_coefficients()`) and multi-start fitting (`n_starts`,
  `start_methods`); the start with the largest penalized objective is returned
  and all starts are listed in `diagnostics$starts`.
- Two-stage benchmark: `update_memberships = FALSE` keeps the initial partition.
- New criteria `CriterionPLIC_BIC` and `CriterionPLIC_AIC`: ordinary likelihood,
  Potts pseudo-likelihood label cost, and EDF penalty. `PottsLogPseudoLik` and
  `PenalizedObjective` are reported in `criteria`.
- New data-generating scenario `"clustered_irregular"` (three disconnected
  regimes on a Latin-square tiling); study support for the extra models
  `"TwoStage-EN"` and `"SCMR-EN-MS"` and for PLIC selection reports.
- Default behaviour is unchanged.

# scmr 0.3.0

- Separate model fitting in `R/` from the article study in `inst/simulation/`.
- Add `penalty = "none"` via `nnet::multinom` and `model = "fixed_clusters"`.
  Retain ridge, lasso, and elastic net for all three model modes.
- Add a common `scmr_fit` object, prediction, coefficient, likelihood and print
  methods, solver diagnostics, and shared fitting paths for spatial G=1/global.
- Consolidate five criteria and their per-cluster components. Penalized EDF is
  an explicitly labelled active-Hessian approximation. Grouped penalties include
  nonzero group-norm curvature; native coefficients are separate from centered
  coefficients used for truth comparison.
- Correct probability orientation for square and single-observation matrices,
  absent-class log-loss, F1 with zero true positives, and glmnet CV logits.
- Update source/destination unit counts after every membership move. Use an
  undirected spatial bonus consistently, including after row standardization.
  Exponential prediction weights preserve the training bandwidth.
- Retain fixed alpha/lambda choices, capture tuning failures, prevent final alpha
  updates after membership iteration exhaustion, and flag changed memberships
  following final tuning. Reject paths that never reach a requested fixed lambda.
- Add explicit DGP controls, smooth coefficient fields, separately controlled
  class imbalance and region imbalance, and clustered heterogeneity strength.
  New locations share the generated Gaussian fields. Correct zero-new-row handling
  and validate generator settings. Seeded fits/generation/splits restore caller RNG.
- Export article evaluation, parameter recovery, diagnostics, and all five
  selection rules from resumable per-run checkpoints. Preserve failure coverage
  and require a complete eligible G grid for selection by default.
- Correct unit-major compact membership export. Replace the old quick-study
  script with an entry point to the maintained article runner.
- Regenerate help from source and remove obsolete combined help aliases.

## Migration

Refit models saved under 0.2.0. Existing entry-point names and global/spatial
subclasses remain, but `fits` now contains engine wrappers under the common
object contract. Use supported methods instead of indexing solver internals.
`fit`, `lambda` (G=1), `alpha_final`, and `obj_trace` remain compatibility fields.
The default package penalty remains grouped; article profiles request ungrouped.

Results are not expected to be bit-for-bit identical to the original script:
several calculations have been corrected, random streams are explicitly scoped,
failed fits are visible, and nonconverged/incomplete G searches are excluded from
selection. The holdout lambda tolerance is a fixed loss margin, not a true 1-SE rule.
