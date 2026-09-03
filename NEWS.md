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
