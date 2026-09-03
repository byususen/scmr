# SCMR package and article simulation: refactor specification

Status: implemented on `refactor/scmr-simulation-v6`; this document records the
accepted design. See `dev/validation.md` for verification and remaining limits.

## Starting point

- Package baseline: version 0.2.0, commit `036f4af291961e57b4b23924bdf51212bcbbfef7`.
- Source: `Pasted text(20260903-053441).txt`, 9,140 lines.
- Source SHA-256: `921bd7679b8e659824df38dcae3a57e10abe44f2776222375b5ca0b59acd211d`.
- The source is identical to `Pasted text(20260903-051643).txt`.
- This work concerns multinomial SCMR and its simulation study. The broader GLM package remains a later project.

The goal is one implementation of each statistical operation in the package.
The article study will call the package to generate data, fit models, and
predict, then evaluate results against simulated truth and export them.

## 1. Public fitting interface

Retain `fit_scmr()` and its existing arguments. Add `model = "fixed_clusters"`,
`penalty = "none"`, and a `cluster` argument for fixed memberships. Append new
arguments without changing existing positional argument meanings.

| Study label | `model` | `penalty` | Membership treatment |
|---|---|---|---|
| Global | `"global"` | `"none"` | One group |
| Global-EN | `"global"` | `"elastic_net"` | One group |
| SCMR | `"scmr"` | `"none"` | Estimated |
| SCMR-EN | `"scmr"` | `"elastic_net"` | Estimated |
| TCMR | `"fixed_clusters"` | `"none"` | Supplied and held fixed |
| TCMR-EN | `"fixed_clusters"` | `"elastic_net"` | Supplied and held fixed |

Keep the package's ridge and lasso options. The six study labels belong to
simulation metadata; the numerical code uses model and penalty arguments.

Input rules:

- `x`: numeric matrix or supported sparse matrix, with stable predictor names.
- `y`: categorical response with an explicit, preserved class order.
- Estimated spatial membership: require `unit_id`, coordinates, and `G`.
- Fixed membership: `cluster` supplies one label per training row. Rows sharing
  a spatial unit must have the same fixed label. Infer `G` from these labels,
  or check consistency when the caller supplies it. Coordinates are optional.
- Global models may receive `unit_id` so criteria that use spatial-unit counts
  can distinguish units from observation rows. Without IDs, treat rows as units.
- Fixed-membership prediction requires a known cluster label for each new row;
  reject unknown labels explicitly.
- Validate training class support, dimensions, finite inputs, and tuning
  settings before fitting. Preserve the model class order during evaluation,
  including when a test subset omits a class.

`scmr_control()` contains algorithm settings, including spatial weights,
feasibility constraints, iteration limits, standardization, solver limits, and
tuning controls. It does not contain scenario names, replicate numbers, paths,
or output filenames.

For ridge, alpha remains zero; for lasso, alpha remains one. Explicit alpha and
lambda choices are respected. Unpenalized fitting uses the multinomial engine
with no regularization and bypasses alpha/lambda tuning.

The article configuration explicitly sets `type_multinomial = "ungrouped"` and
its alpha grid. Changes to package defaults are a separate documented decision.

## 2. Common fitted-object contract

Keep existing `scmr_global` and `scmr_spatial` subclasses and add a common
`scmr_fit` parent. Introduce `scmr_fixed` for fixed-membership models. These are
the implemented classes; existing prediction and print dispatch must remain compatible.

All fitted objects expose the following core fields:

| Field | Meaning |
|---|---|
| `schema_version` | Version of the object structure |
| `call`, `model`, `penalty`, `G` | Model identity and call |
| `class_levels`, `x_colnames` | Ordered response classes and predictors |
| `fits` | List of local engine fits, length `G`, including global `G = 1` |
| `alpha` | Shared fitted alpha; `NA_real_` when unpenalized |
| `lambda_used` | Numeric vector of length `G`; `NA` when unpenalized |
| `prob_train` | Training probability matrix, observations by response classes |
| `row_group` | Integer membership aligned to the original training row order |
| `group_unit`, `PP_unit`, `unit_levels` | Unit memberships with explicit ordering |
| `control`, `seed` | Resolved fitting settings and root seed |
| `criteria` | Named criteria, log-likelihood, dimensions, and EDF method metadata |
| `cluster_components` | One row per cluster with counts, likelihood, EDF, and penalties |
| `converged`, `convergence_reason` | Actual optimization status and stopping reason |
| `diagnostics` | Iteration history, tuning tables, solver status, and warnings |

Retain existing fields such as global `fit`, global `lambda`, spatial
`alpha_final`, and `obj_trace` as documented compatibility aliases where needed.
An object's common fields must not depend on which study wrapper called it.

Implemented methods:

- `predict()`: probabilities or response classes; spatial models also return
  membership weights. Predictions preserve row and class order for every batch
  size, including one row and square probability matrices.
- `coef()`: cluster-by-class-by-term array on the input predictor scale, with an
  explicit choice of native or sum-to-zero parameterization.
- `logLik()`: ordinary training multinomial log-likelihood, with the meaning of
  any supplied degrees of freedom documented.
- `print()`: model, penalty, cluster count, and convergence summary.

Keep native coefficient support for penalty/EDF calculations separate from
coefficients transformed for truth comparison. Centering a coefficient vector
can change which individual entries are zero.

## 3. Package and study boundary

| Operation in the attachment | Destination |
|---|---|
| Penalized and unpenalized local fitting | Package engine helpers |
| Shared iteration and constrained membership updates | Package spatial fitter |
| Alpha/lambda tuning | Package tuning functions |
| Hessian, active-set EDF, likelihood, and all five criteria | Package |
| Probability prediction and coefficient extraction | Package methods |
| Probability validation and general classification metrics | Package utilities |
| Reusable spatial data generation | `simulate_scmr_data()` with explicit arguments |
| Article coefficient values, sample sizes, scenario grids, and repeats | Study configuration |
| Coefficient MSE, support recovery against truth, and ARI reporting | Study evaluation |
| Replicate loops and model/cluster grids | Study runner |
| CSV exports, run keys, checkpoints, and resume logic | Study I/O |
| Best-G selection, replicate summaries, and report tables | Study post-processing |

Fitting functions return diagnostic tables. The runner attaches `Repeat`,
`Scenario`, and other study metadata and saves them. Fitting does not call
`append_row()`, change directories, create output folders, or delete files.

Keep only the active definition of each duplicated function from the attachment.
Replace dependencies on global variables with explicit arguments or resolved
control settings. Internal random splits and iteration order must derive from
documented seed inputs, independent of plotting and exporting operations.

Implemented study layout:

- `inst/simulation/config.R`: pilot/full configurations and study parameters.
- `inst/simulation/run.R`: data generation, splitting, and package calls.
- `inst/simulation/evaluate.R`: truth-based and regional evaluation.
- `inst/simulation/io.R`: fixed table schemas, writes, manifests, and resume.
- `inst/simulation/summarize.R`: best-G selection and replicate summaries.
- `inst/scripts/02_run_article_simulation.R`: short entry point.

The study loads `scmr` and calls supported interfaces. Its helpers must not
redefine package fitting, prediction, or criterion functions. The current
three-class restriction belongs to the supplied DGP; fitting constraints must
be validated separately and described accurately.

## 4. Numerical behavior to preserve or correct explicitly

1. Correct probability orientation, absent-class log-loss, zero-TP F1, and the
   sparse symmetric spatial-bonus calculation. These four defects were
   reproduced in isolated checks of the supplied implementation.
2. Check unit-count updates during membership moves. The attachment changes
   labels and class counts but does not update `unit_count` after each move.
3. Consolidate all five criteria in the package. The current script computes
   some only during CSV post-processing. EDF remains an approximation with an
   explicit method label. The ungrouped active-set Hessian calculation must
   not silently serve as a grouped-penalty EDF implementation.
4. Define the G=1 spatial fit to use the global fitting/tuning path. Compare
   probabilities, coefficients, likelihood, and applicable criteria using
   identical data, unit IDs, controls, and seeds. Copying CSV rows is not an
   independent check of fitting equivalence. Record fit reuse in runtime data.
5. Distinguish membership stopping from local-solver convergence. A fixed
   partition does not imply successful coefficient optimization. Check whether
   final alpha retuning changes the conditions behind a convergence claim.
6. Align coefficient parameterization before MSE comparison. Keep prediction
   mixtures and membership-weighted coefficient summaries explicitly distinct;
   they do not generally define the same probability model.
7. Version study output schemas and include a configuration fingerprint in
   resume checks. Mark a run complete only after its required outputs succeed.

Record deliberate numerical changes in `NEWS.md`. Migration checks should
compare against explicitly corrected reference behavior, without requiring
known defects to be reproduced.

## 5. First implementation increment and acceptance checks

The first code increment corrected package probability and metric helpers
and added focused regression checks. Those checks establish the baseline for
the remaining statistical code.

Its acceptance cases are:

- Preserve an already valid non-symmetric K-by-K probability matrix.
- Handle a single predicted observation without losing dimensions or labels.
- Obtain log-loss `-log(0.6)` when the observed class has predicted probability
  0.6, even if other model classes are absent from the evaluation subset.
- Obtain macro-F1 `1/3` when one of three equally supported classes is perfectly
  predicted and the other two are completely confused.
- Preserve existing ordinary-case metric behavior except where a documented
  defect requires correction.

The object/validation contract, local engines, tuning, EDF, criteria, clustering,
fixed memberships and study wrappers have now been implemented. A complete
small smoke study has run; the main 550-fit article study remains to be run.
The verification record is in `dev/validation.md`.
