# Article simulation guide

The supplied standalone script has been split into reusable package code and
study orchestration. This directory contains no fitting, probability, or EDF
implementation. It calls `scmr::simulate_scmr_data()`, `scmr::fit_scmr()`,
`predict()`, `coef()`, and the package metric functions.

## First run

From the repository root:

```sh
R CMD INSTALL .
Rscript inst/scripts/02_run_article_simulation.R smoke ./results
```

The default entry point is `smoke`. Choose `pilot` or `main` explicitly.
The old `01_run_simulation_study.R` entry point delegates to this runner.

| Profile | Observed / new locations | Repeats | Active / inactive predictors | SCMR-EN G |
|---|---:|---:|---:|---|
| smoke | 360 / 36 | 1 | 2 / 2 | 1, 2 |
| pilot | 600 / 60 | 2 | 5 / 5 | 1, 3, 6 |
| main | 3000 / 300 | 10 | 5 / 5 | 1 through 10 |

All profiles include global, balanced-cluster, imbalanced-cluster, and smooth
scenarios. Global and Global-EN run everywhere. Unpenalized SCMR uses G=1 in
the global scenario, G=6 in clustered scenarios, and is omitted for smooth
truth. TCMR and TCMR-EN run only for clustered truth. Thus smoke requests 23
model runs, pilot 54, and main 550. G=1 study rows can reuse global fits; the
package independently implements and tests the same G=1/global fitting path.

## Customize settings

Load the helpers from the installed package:

```r
for (name in c("config", "io", "evaluate", "summarize", "run")) {
  source(system.file("simulation", paste0(name, ".R"), package = "scmr"))
}
cfg <- article_config("pilot", output_dir = "results",
                      n_repeats = 3L, class_balance = c("balanced", "imbalanced"))
cfg$fit_control$max_iter <- 20L
cfg$export_membership <- TRUE
result <- run_article_simulation(cfg)
result$output_dir
```

Prefer overrides through `article_config()`; if editing nested controls,
validate them again with `do.call(scmr_control, cfg$fit_control)` and
`do.call(scmr_simulation_control, cfg$dgp_control)` before running.

`n_active`, `n_inactive`, `eta`, `class_balance`, and `heterogeneity_strength`
can be vectors. Their combinations define the study grid. Heterogeneity
strength applies only to clustered coefficients; nonclustered scenarios run
once per other configuration with strength recorded as one.

The main profile preserves the attachment's sample sizes, 30% region-by-class
test split, seed formula, predictor covariance settings, six-region proportions,
coefficient patterns, class shifts, and tuning defaults. It explicitly uses
ungrouped multinomial EN and alpha candidates .3, .5, .7, .9. The package default
remains grouped. The smaller profiles change fitting constraints/path length
only where specified in `config.R`.

## Reproducibility and resume

A fingerprint of resolved configuration, package/study function definitions,
engine versions, and R version creates a subdirectory such as
`smoke-012345abcdef`. Output location and the resume flag do not affect the
fingerprint. Results generated with changed settings or implementations use a
new directory. Timing values may differ across machines; numerical reproducibility
also depends on numerical libraries and the R random-number implementation.

- `manifest.rds`, `config.R`, `session-info.txt`: settings, source-script hash,
  package versions, and runtime information.
- `data/<dataset-id>.rds`: generated/split datasets and DGP diagnostics. This
  preserves the expensive Gaussian fields across restarts.
- `runs/<run-id>.rds`: a completed model's result and output tables, committed
  after fitting and evaluation succeed. `save_fits = TRUE` also stores fits.
- `errors/<run-id>.rds`: failures with study metadata and error messages.

Repeating the command skips complete bundles, retries errors, and reconstructs
CSV reports. Nonconverged fits are completed computations and are not silently
retried. Change solver settings to give them a new configuration. `resume = FALSE`
refuses to overwrite an existing matching study. No working directory changes,
Windows-specific paths, or cleanup of old study outputs occur.

Checkpoints are intended for one writer per output directory. Run separate
configurations in separate directories when launching concurrent jobs.

## Reports

| Output | Contents |
|---|---|
| `expected_runs.csv`, `run_coverage.csv`, `errors.csv` | Requested runs, actual status, failures |
| `result_by_k.csv` | One row per fitted model, criteria, convergence, runtimes, Train/Test/New metrics |
| `overall.csv`, `class.csv`, `region.csv`, `region_class.csv` | Classification and truth-based evaluation |
| `parameters_region.csv` | True/estimated coefficient means, bias, MSE by region/class/term |
| `parameters_cluster.csv`, `coefficients_cluster.csv` | Coefficient recovery and native/centered coefficients by fitted cluster |
| `selection.csv` | Predictor support recovery and centered coefficient support recovery |
| `cluster_components.csv`, `criterion_audit.csv` | Per-cluster likelihood, EDF, penalties, weighted sample-size checks |
| `diagnostics_tuning.csv`, `diagnostics_iterations.csv`, `diagnostics_solver.csv` | Tuning choices/failures, membership changes, solver status |
| `dgp_diagnostics.csv` | Region counts/proportions and realized predictor correlations |
| `best_G_<rule>.csv`, `best_G_frequency_<rule>.csv` | Best eligible G and replicate selection frequencies for each of five rules |
| `selection_coverage_<rule>.csv` | Whether every requested G was available and eligible |
| `result_summary.csv`, `*_summary.csv`, `summary_selected_<rule>.csv` | Means, finite replicate counts, Monte Carlo SEs and t-based 95% intervals |

Optional `parameters_observation.csv` and `memberships.csv` are controlled by
`export_per_observation` and `export_membership`; otherwise their reports are
empty. Large per-location exports are disabled by default. The new schemas
replace the standalone script's append-only CSVs; downstream analyses should
use column names documented by the generated reports.

`FitTimeSec` is zero for reused G=1 rows; `SourceFitTimeSec` records the original
fit time and `ReusedGlobalG1` identifies reuse. With `save_fits = TRUE`, a reused
row stores its global source fit. Data-generation and report-writing time are
not included in fit time.

## Evaluation conventions

Predictions mix local class probabilities using membership weights. Coefficient
recovery at a location uses the same weights applied to centered coefficients;
this weighted coefficient field is a recovery summary, not a claim that the
probability mixture equals a single multinomial logit with those coefficients.

The true coefficients and `coef(fit)` both sum to zero across classes. This
makes coefficient MSE comparable between `nnet`'s baseline representation and
`glmnet`'s symmetric representation. Predictor support uses native coefficients;
coefficient-level recovery uses centered coefficients on both sides. Native
zeros remain separately available and are the basis of penalty/EDF calculations.

ARI compares hard predicted memberships against discrete true clusters. It is
undefined for the smooth scenario. Class metrics preserve all model classes;
F1 is zero when TP=0 and FP+FN>0, and undefined when TP=FP=FN=0. Macro metrics
omit undefined entries. Undefined precision contributes zero to the existing
support-weighted precision convention. Brier score is the average sum across
classes; probability MSE averages across rows and classes.

`class_balance = "balanced"` means no extra intercept shift. It does not force
equal class counts. Balanced/imbalanced clustered scenarios refer to location
counts in the six regions and can be crossed independently with class imbalance.
Observed and new locations share Gaussian fields; new-location evaluation is
spatial interpolation. Predictor standardization over all generated locations
is part of the DGP. Model/validation standardization is fitted within each
training fit.

## Criteria and EDF

Let `n` be the number of training rows, `m_g` the unique units in cluster `g`,
`n_g` its rows, `L_g` its ordinary multinomial log-likelihood, and `d_g` its EDF.
Write `L = sum(L_g)`, `d = sum(d_g)`, and `d0 = G (K - 1) (p + 1)`.

- SCR-BIC original: `-2 L + log(n) d0`.
- SCR-BIC effective: `-2 L + log(n) d`.
- SCR-AIC effective: `-2 L + 2 d`.
- With `a_g = (mean(m_g)/m_g)^gamma` and
  `w_g = n a_g / sum(a_h n_h)`, CB-BIC is
  `sum(-2 w_g L_g + log(max(m_g, 2)) d_g)`.
- CB-AIC replaces the last penalty by `2 d_g`.

The weights satisfy `sum(w_g n_g) = n`. The spatial fitting objective is
`L + phi * sum_{u<v} W_uv I(z_u=z_v)` and is reported separately. G=1 uses the
global path, which reports zero spatial bonus and identical global criteria.

For unpenalized fits, EDF is `(K-1) rank([1, X])` within each cluster. For EN,
EDF is `tr((H + P)^+ H)` on native active coordinates, where H is the full
multinomial likelihood Hessian and P is penalty curvature on the solver's
predictor scale. Intercepts are unpenalized. Ungrouped L1 has zero curvature
within a fixed nonzero support. Grouped penalties include the vector-norm
curvature of each nonzero group, in addition to ridge curvature. A pseudoinverse
handles nonidentifiability; reported EDF is bounded to the numerical rank of H.

These are conditional local approximations. They do not estimate the additional
degrees of freedom from choosing memberships, alpha, or lambda. The study
implements these criteria for comparison; it does not establish their statistical
validity for every sampling design.

## Failures and intentional changes

All fits remain in `result_by_k.csv`, including those that did not converge.
By default only converged fits enter summaries. A replicate can select G only
if every requested SCMR-EN candidate has a finite criterion and meets the
eligibility rule. This avoids silently selecting among a reduced grid after
failures. The count of usable replicates is reported alongside every summary.
With one replicate, standard errors and confidence intervals are unavailable.

The source script is preserved by its SHA-256 in the manifest:
`921bd7679b8e659824df38dcae3a57e10abe44f2776222375b5ca0b59acd211d`.
This refactor intentionally corrects its probability/metric, membership-count,
spatial-bonus, and stopping/reporting defects. It is not expected to reproduce
those defects or produce identical output files. See the repository's `NEWS.md`
for migration details.
