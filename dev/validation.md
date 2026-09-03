# Validation of the 0.3.0 refactor

Date: 2026-09-03. Branch: `refactor/scmr-simulation-v6`.
Baseline: `036f4af291961e57b4b23924bdf51212bcbbfef7` (0.2.0).

## Environment and package check

- Native R 4.3.3 on Ubuntu 24.04.3, x86_64 Linux.
- glmnet 4.1-8, nnet, Matrix, MASS, testthat and roxygen2 available locally.
- `R CMD INSTALL .`: successful.
- `R CMD build .`: successful.
- `R CMD check --no-manual scmr_0.3.0.tar.gz`: **Status: OK**.
- Tests during package checking: **136 passed, 0 failed, 0 warnings, 0 skipped**.
- Generated Rd files have no duplicate aliases or signature mismatches.

The manual/PDF build and other R/platform versions were not checked. The
repository's existing GitHub Actions workflow remains available for remote checks.

## Numerical and interface checks

- Square probability matrices retain orientation; single rows preserve class
  dimensions; invalid probabilities and labels are rejected.
- Absent-class log-loss preserves the model's full class probabilities.
  Zero-TP F1 and associated macro/weighted metrics use the corrected formula.
- Actual spatial G=1 and global fits agree in training probabilities,
  coefficients, selected lambdas and criteria for none, ridge, lasso and EN,
  including tuned penalized fits.
- Multinomial CV prevalidated logits transform to the expected probabilities.
- Centered unpenalized coefficients reproduce predictions on the original
  predictor scale for binary and three-class, one-predictor problems.
- Fixed memberships and prediction labels remain aligned. Estimated membership
  fits retain unit and class-count feasibility and keep repeated units together.
- Likelihood, nominal BIC, unit-count CB-BIC and CB-AIC match independent
  calculations. Ridge EDF falls as penalization increases in the checked example.
  Sparse spatial bonus matches a three-node hand calculation.
- All four DGP scenarios return normalized probabilities and sum-to-zero
  coefficients. Clustered region counts match the intended integer proportions.
  Zero new rows, no active predictors and zero heterogeneity are handled.
- Fits, generation and stratified splitting preserve the caller's RNG state.

## Study integration and smoke run

An installed-package integration test generates data, fits four global/G=1
study rows, writes reports, resumes in the same R session without rewriting
completed run bundles, and rejects incompatible overwrite requests. It also
checks that a missing G grid cannot silently select a best G.

The final smoke profile used 360 observed and 36 new locations, two active and
two inactive predictors, all four scenarios, one replicate, and SCMR-EN G=1,2.
It completed **23 of 23 requested fits with zero execution errors**. Twenty-two
fits converged; one fit is retained with a nonconvergence flag. That fit is
excluded from summaries by the default policy. All five criteria were finite,
run IDs were unique, and every cluster-balanced weighted sample-size audit
agreed with the training row count to within 1e-8.

All five best-G reports contained one selected row for each of the four
scenarios. The default policy requires every candidate G to be eligible.
The small integration study and full smoke outputs are verification artifacts,
not estimates for the article. The main 550-fit study has **not** been run.

## Review limits

The tests establish implementation behavior on the checked cases. They do not
establish theoretical validity of the conditional EDF/IC approximations, global
optimality of spatial membership search, or absence of separation in unpenalized
fits. Inspect solver diagnostics and run the pilot/main studies before drawing
scientific conclusions. Regenerate results instead of combining 0.2.0 or old
standalone-script CSVs with the new schemas.
