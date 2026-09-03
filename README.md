# scmr: Spatially Clustered Multinomial Regression

`scmr` fits multinomial regression models with global, estimated spatial-cluster,
or supplied fixed-cluster memberships. Each supports no penalty, ridge, lasso,
and elastic net. Model fitting and the article simulation use the same package
functions.

## Install and try a model

From a local checkout:

```sh
R CMD INSTALL .
```

Or install the published GitHub version:

```r
# install.packages("remotes")
remotes::install_github("byususen/scmr")
```

The following example needs the 0.3.0 refactor:

```r
library(scmr)
dat <- simulate_scmr_data(n_obs = 120, n_new = 12, p = 4, active = 1:2)
train <- dat$observed

fit <- fit_scmr(train$X, train$y, model = "global",
                penalty = "elastic_net", alpha = 0.5, seed = 123)
predict(fit, dat$new$X, type = "prob")
coef(fit)               # cluster x class x term; sum-to-zero coefficients
logLik(fit)             # ordinary training log-likelihood
fit$criteria
fit$diagnostics
```

## Fit the six study models

| Study model | `model` | `penalty` | Additional arguments |
|---|---|---|---|
| Global | `"global"` | `"none"` | — |
| Global-EN | `"global"` | `"elastic_net"` | — |
| SCMR | `"scmr"` | `"none"` | `G`, `unit_id`, `coords` |
| SCMR-EN | `"scmr"` | `"elastic_net"` | `G`, `unit_id`, `coords` |
| TCMR | `"fixed_clusters"` | `"none"` | `cluster` |
| TCMR-EN | `"fixed_clusters"` | `"elastic_net"` | `cluster` |

Use `penalty = "ridge"` for alpha zero or `"lasso"` for alpha one. Elastic-net
alpha can be fixed or tuned over `alpha_grid`. A supplied positive `lambda`
fixes the penalty strength; spatial/fixed models also accept one lambda per
cluster. With fixed lambda and no explicit alpha, elastic net uses the median
of the resolved alpha grid.

```r
ctrl <- scmr_control(type_multinomial = "ungrouped",
                     alpha_grid = c(.3, .5, .7, .9),
                     min_units = 5, min_per_class = 2)
spatial <- fit_scmr(train$X, train$y, model = "scmr", G = 2,
                    unit_id = train$unit, coords = train$coords,
                    control = ctrl, seed = 123)
predict(spatial, dat$new$X, new_unit_id = dat$new$unit,
        new_coords = dat$new$coords)
spatial$converged
spatial$convergence_reason
```

For fixed memberships, `cluster` supplies a label per training row. Prediction
requires the corresponding known labels through `predict(fit, newx,
cluster = ...)`. TCMR uses simulated true cluster labels; it is an oracle benchmark.

## Inputs and predictions

- `x` is a finite numeric predictor matrix. Sparse input is accepted and
  converted to a dense matrix internally. Predictor names must be unique.
- `y` contains every declared training class. Fitting supports two or more
  classes; the supplied data generator has three classes.
- Coordinates have two columns. Repeated rows for a spatial unit share one
  coordinate pair and move together. Initial labels follow sorted unique unit IDs.
- Every fitted cluster must meet `min_units` and `min_per_class`.
- Probabilities always have one row per observation and one column per model
  class, even for one-row batches or test sets missing a response class.
- Spatial predictions reuse memberships for known units. Unseen units receive
  neighbor-based membership weights, and predictions mix local probabilities.
  Exponential weights reuse the training bandwidth.

All fits inherit from `scmr_fit`. Inspect `prob_train`, `row_group`, `PP_unit`,
`alpha`, `lambda_used`, `criteria`, `cluster_components`, and `diagnostics`.
A spatial model with `G = 1` uses the global fitting and tuning path.

`coef(fit)` returns coefficients on the original predictor scale, centered
across classes. `coef(fit, parameterization = "native")` returns the solver's
representation. Native zeros determine penalized support/EDF; centering can
change which individual coefficients are zero.

## Run the article simulation

Install this checkout, then run from the repository root:

```sh
Rscript inst/scripts/02_run_article_simulation.R smoke ./results
Rscript inst/scripts/02_run_article_simulation.R pilot ./results
Rscript inst/scripts/02_run_article_simulation.R main ./results
```

Start with `smoke`. The main profile uses 3,000 observed and 300 new locations,
10 repetitions, four scenarios, and SCMR-EN `G = 1:10`: 550 requested model runs.
Gaussian-field generation uses dense covariance matrices, so the main study
requires substantially more time and memory than the smoke profile.

The same command resumes matching completed runs. Settings and code identify
a separate output directory; changing a configuration never appends to old
results. Per-run RDS checkpoints are the source of truth, and CSV reports can
be regenerated from them.

The study files have separate responsibilities:

| File | Responsibility |
|---|---|
| `inst/simulation/config.R` | Profiles, scenarios, model grids, DGP and fitting settings |
| `inst/simulation/run.R` | Generate data, split samples, call package functions |
| `inst/simulation/evaluate.R` | Compare predictions, coefficients, and memberships with truth |
| `inst/simulation/io.R` | Checkpoints, manifests, and CSV output |
| `inst/simulation/summarize.R` | Five best-G rules, summaries, selection frequencies, coverage |

See [the simulation guide](inst/simulation/README.md) for custom settings,
outputs, and reproducibility details.

## Statistical conventions

The package default remains grouped multinomial penalization; the article
profile explicitly requests **ungrouped** penalization. Global models use
stratified cross-validation. Multi-cluster models use class-stratified holdout
loss to choose a shared alpha and cluster-specific lambdas. The local
`lambda_1se_tol` is a fixed loss tolerance, not an estimated standard error.

The five criteria are SCR-BIC with nominal or effective degrees of freedom,
SCR-AIC with effective degrees of freedom, CB-BIC, and CB-AIC. They use ordinary
multinomial likelihood; the spatial bonus is reported separately. Penalized EDF
is a local active-Hessian approximation conditional on the selected support,
penalties, and memberships. It does not account for the search over partitions
or tuning parameters. Formula details appear in the simulation guide.

Check convergence before interpreting a fit. Reaching the iteration limit,
local solver failure, or membership changes after final tuning is recorded.
Tiny-movement stopping remains a configurable approximation. Study reports
retain unsuccessful convergence statuses and exclude those fits from summaries
by default. A replicate selects G only when all requested candidates are eligible.

## Changes from 0.2.0

This refactor adds unpenalized and fixed-membership fits, common fitted-object
methods, the smooth data scenario, and actual balanced/imbalanced region sizes.
It also corrects probability orientation, absent-class log-loss, zero-TP F1,
CV probability extraction, membership export ordering, spatial bonus evaluation,
and unit counts during membership updates.

Results should be regenerated. Serialized 0.2.0 fitted objects should be refitted
with 0.3.0; low-level engine storage and study output schemas have changed. See
[NEWS](NEWS.md) and [the refactor specification](dev/scmr-refactor-spec.md).
