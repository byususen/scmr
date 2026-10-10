# =============================================================================
# EMPIRICAL ARTICLE CODE v4 - SCMR-EN on KSA + Sentinel-1/2
#
# Requires scmr >= 0.5.0 (install from the scmr_pkg folder next to this file).
#
# One method (SCMR-EN: multinomial elastic net per cluster, Potts prior on
# segment labels with phi estimated by pseudo-likelihood, PLIC for G), two
# predictor specifications. Every method within a design uses the SAME
# specification, training rows, penalty tuning and prediction rule, so that
# differences come only from how spatial heterogeneity is modelled.
#   S (static)  month dummies + Sentinel features
#   D (dynamic) S + one-hot previous phase (Markov transition intercepts)
#
# Designs
#   E2_mapping    MAIN, spec S. 5-fold CV over whole KSA segments stratified
#                 by regency; test segments are new locations (no survey).
#   E1_nowcast    Spec D. Train on months 1..18, predict months 19..24 at the
#                 survey points with the observed previous phase.
#   E3_imputation Spec D. 15% of monthly records hidden; imputed by exact
#                 forward-backward under the fitted model (scmr_impute_waves),
#                 versus carry-forward.
# Methods (suffix -S / -D = specification)
#   Global-EN     one multinomial EN (CV lambda.min)
#   SCMR-EN       proposed; G in the grid, chosen by PLIC-AIC
#   TwoStage-EN   k-means of segment centroids, no label updates
#   Regency-EN    one model per regency (administrative strata)
#   SCR-original  port of the SCR reference code (unpenalized, phi = 1, BIC)
#   GW-EN         geographically weighted multinomial EN
#   LightGBM      gradient boosting with the same features plus coordinates
#   ablations at the G chosen by PLIC-AIC: SCMR-EN-<spec>-phi1 (phi fixed at 1)
#   and SCMR-EN-S-unit (labels per sub-segment instead of per segment)
# Resume-safe: one RDS per (design, fold, method, G) under output_dir/runs.
# Fitted dynamic runs of v3c (same settings) are reused when their penalty
# matches. Environment overrides: SCMR_SYNTAX_DIR, SCMR_DATA_FILE,
# SCMR_OUTPUT_DIR, SCMR_PROFILE ("main", "pilot", "ci"), SCMR_DESIGNS,
# SCMR_FOLDS, SCMR_G_GRID, SCMR_USE_GW, SCMR_USE_SCR.
# =============================================================================

# =============================================================================
# 0. SETTINGS
# =============================================================================
syntax_dir <- Sys.getenv("SCMR_SYNTAX_DIR", "C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax")
data_file <- Sys.getenv("SCMR_DATA_FILE", file.path(syntax_dir, "df_month_rec_final.csv"))
output_dir <- Sys.getenv("SCMR_OUTPUT_DIR", file.path(dirname(syntax_dir), "Results", "Real data v4"))
results_root <- file.path(dirname(syntax_dir), "Results")
reuse_cache_dirs <- file.path(results_root, c("Real data v3c", "Real data v3b", "Real data v3"), "cache")
reuse_v3c <- file.path(results_root, "Real data v3c")
run_profile <- Sys.getenv("SCMR_PROFILE", "main")

designs_all <- c("E2_mapping", "E1_nowcast", "E3_imputation")
designs_to_run <- designs_all
if (nzchar(Sys.getenv("SCMR_DESIGNS"))) designs_to_run <- strsplit(Sys.getenv("SCMR_DESIGNS"), ",")[[1]]
use_gw <- TRUE
if (nzchar(Sys.getenv("SCMR_USE_GW"))) use_gw <- as.logical(Sys.getenv("SCMR_USE_GW"))
use_scr_original <- TRUE
if (nzchar(Sys.getenv("SCMR_USE_SCR"))) use_scr_original <- as.logical(Sys.getenv("SCMR_USE_SCR"))
use_lightgbm <- requireNamespace("lightgbm", quietly = TRUE)
scr_maxitr <- 30L

profiles <- list(
  ci    = list(folds = 1L, G_grid = c(1L, 2L), n_starts = 2L, max_iter = 10L, gw_anchors = 30L,
               gw_k = c(10, 20), coef_k = 3L, coef_k_unit = 10L, coef_anchors = 60L),
  pilot = list(folds = 1L, G_grid = c(1L, 3L, 6L), n_starts = 3L, max_iter = 20L, gw_anchors = 300L,
               gw_k = c(50, 100, 200), coef_k = 3L, coef_k_unit = 12L, coef_anchors = 600L),
  main  = list(folds = 1:5, G_grid = c(1L, 3L, 6L, 8L), n_starts = 3L, max_iter = 20L, gw_anchors = 300L,
               gw_k = c(50, 100, 200), coef_k = 3L, coef_k_unit = 12L, coef_anchors = 600L))
if (!run_profile %in% names(profiles)) stop("Unknown profile: ", run_profile)
prof <- profiles[[run_profile]]
if (nzchar(Sys.getenv("SCMR_FOLDS"))) prof$folds <- as.integer(strsplit(Sys.getenv("SCMR_FOLDS"), ",")[[1]])
if (nzchar(Sys.getenv("SCMR_G_GRID"))) prof$G_grid <- as.integer(strsplit(Sys.getenv("SCMR_G_GRID"), ",")[[1]])
n_folds <- 5L
outer_seed <- 2026L
train_months_E1 <- 18L
mask_prop_E3 <- 0.15

keep_vi_optical <- c("CI_red_edge", "EVI", "EVI2", "GNDVI", "MNDWI", "MSAVI",
                     "MSI", "NDRE", "NDVI", "NDWI", "NIRv", "PSRI")
keep_radar <- c("rvi", "vh", "vv")
phase_labels <- c("1" = "Early vegetative", "2" = "Late vegetative", "3" = "Generative",
                  "4" = "Harvest", "5" = "Land preparation", "6" = "Other")
alpha_grid <- c(0.5, 0.9)
min_units <- 5L                  # segments per cluster (10 sub-segments for unit labels)
min_per_class <- 5L
k_neighbors <- 8L                # Potts graph: 8 nearest segments
membership_rule <- "potts"       # new locations: Potts conditional of the 8 nearest segments
# Local penalty = lambda.1se of the global model of the same specification
# times a multiplier chosen by inner validation (held-out training segments,
# log-loss of the design's own prediction mode, two-stage fit with G = 3).
local_lambda_rule <- "lambda.1se"
tune_multipliers <- c(1, 4, 16, 64)
tune_G <- 3L
show_progress <- TRUE

# =============================================================================
# 1. PACKAGES
# =============================================================================
if (!requireNamespace("scmr", quietly = TRUE) || utils::packageVersion("scmr") < "0.5.0") {
  stop("Install scmr >= 0.5.0 first (install.packages(<scmr_pkg folder>, repos = NULL, type = 'source')).")
}
suppressPackageStartupMessages(library(scmr))
cat("scmr", as.character(utils::packageVersion("scmr")), "| profile", run_profile,
    "| designs", paste(designs_to_run, collapse = ","), "| LightGBM", use_lightgbm, "\n")
dir.create(file.path(output_dir, "runs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "errors"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "cache"), recursive = TRUE, showWarnings = FALSE)
banner <- function(...) cat("\n==== ", ..., " ====\n", sep = "")

# =============================================================================
# 2. DATA
# =============================================================================
is_model_feature <- function(nm) {
  if (nm %in% c("id_segmen", "id_subsegmen", "id_month", "phase", "STRATA", "year", "month",
                "month_end", "lati_mean", "long_mean", "s2r_med_lati_mean", "s2r_med_long_mean")) return(FALSE)
  grepl("^(s2r_med_|s1_med_|s2r_month_|s1_month_|s2r_w2_|s1_w2_)", nm)
}
is_selected <- function(nm, keep) any(vapply(keep, function(v) grepl(paste0("(^|_)", v, "(_|$)"), nm), logical(1)))
keep_feature_reduced <- function(nm) {
  if (grepl("^(s2r_med_|s1_med_)", nm)) {
    if (grepl("_lag[4-9]+$", nm)) return(FALSE)
    if (grepl("^s2r_", nm)) return(is_selected(nm, keep_vi_optical))
    if (grepl("^s1_", nm)) return(is_selected(nm, keep_radar))
  }
  if (grepl("^(s2r_month_|s1_month_|s2r_w2_|s1_w2_)", nm)) {
    if (grepl("^(s2r_w2_|s1_w2_)", nm) && grepl("_max_|_min_", nm)) return(FALSE)
    if (grepl("^s2r_", nm)) return(is_selected(nm, keep_vi_optical))
    if (grepl("^s1_", nm)) return(is_selected(nm, keep_radar))
  }
  FALSE
}

banner("LOAD DATA")
raw <- if (grepl("\\.rds$", data_file, ignore.case = TRUE)) readRDS(data_file) else
  utils::read.csv(data_file, stringsAsFactors = FALSE,
                  colClasses = c(id_segmen = "character", id_subsegmen = "character", id_month = "character"))
raw$id_segmen <- as.character(raw$id_segmen); raw$id_subsegmen <- as.character(raw$id_subsegmen)
raw <- raw[, !grepl("last", names(raw)), drop = FALSE]
raw$phase <- as.character(as.integer(round(as.numeric(raw$phase))))
raw <- raw[!is.na(raw$phase) & raw$phase != "12", , drop = FALSE]
raw$phase[raw$phase %in% c("6", "7", "8")] <- "6"
raw$regency <- substr(raw$id_segmen, 1, 4)
raw$time <- (raw$year - min(raw$year)) * 12L + raw$month
x_vars <- names(raw)[vapply(names(raw), is_model_feature, logical(1))]
x_vars <- x_vars[vapply(raw[x_vars], is.numeric, logical(1))]
x_vars <- x_vars[vapply(x_vars, keep_feature_reduced, logical(1))]
if ("s1_w2_min_rvi" %in% names(raw)) raw$s1_w2_min_rvi[raw$s1_w2_min_rvi < -10] <- NA
dat <- raw[, unique(c("id_segmen", "id_subsegmen", "phase", "regency", "month", "time",
                      "lati_mean", "long_mean", x_vars)), drop = FALSE]
dat <- dat[stats::complete.cases(dat), , drop = FALSE]
dat <- dat[order(dat$id_subsegmen, dat$time), , drop = FALSE]
rm(raw); invisible(gc())

class_levels <- as.character(1:6)
class_levels <- class_levels[class_levels %in% unique(dat$phase)]
dat$phase <- factor(dat$phase, levels = class_levels)
lagd <- scmr_lag_design(dat$phase, dat$id_subsegmen, dat$time, class_levels)
lag_cols <- lagd$columns
unit_xy <- stats::aggregate(cbind(lati_mean, long_mean) ~ id_subsegmen, data = dat, FUN = mean)
dat$unit_lat <- unit_xy$lati_mean[match(dat$id_subsegmen, unit_xy$id_subsegmen)]
dat$unit_lon <- unit_xy$long_mean[match(dat$id_subsegmen, unit_xy$id_subsegmen)]
coords_all <- cbind(dat$unit_lat, dat$unit_lon)
# Labels are carried by survey segments (blocks); the ablation "unit" labels
# sub-segments. Coordinates: segment centroid / sub-segment centroid per row.
seg_id <- dat$id_segmen
seg_tab <- stats::aggregate(cbind(unit_lat, unit_lon) ~ id_segmen,
                            data = unique(dat[, c("id_segmen", "id_subsegmen", "unit_lat", "unit_lon")]), FUN = mean)
seg_xy <- cbind(seg_tab$unit_lat[match(seg_id, seg_tab$id_segmen)], seg_tab$unit_lon[match(seg_id, seg_tab$id_segmen)])
sub_id <- dat$id_subsegmen

month_mm <- stats::model.matrix(~ factor(month, levels = 1:12) - 1, dat)
colnames(month_mm) <- paste0("month", 1:12)
X_base <- cbind(month_mm, as.matrix(dat[, x_vars]))
X_full <- cbind(X_base, lagd$x)                     # lag columns NA where unavailable
X_lagfill <- X_full; X_lagfill[is.na(X_lagfill)] <- 0
XY <- cbind(lat = dat$unit_lat, lon = dat$unit_lon)   # coordinates as LightGBM features
cat("Rows:", nrow(dat), "| units:", length(unique(dat$id_subsegmen)), "| segments:",
    length(unique(dat$id_segmen)), "| satellite predictors:", length(x_vars),
    "| rows with previous phase:", sum(lagd$available), "\n")
print(table(Phase = dat$phase, Regency = dat$regency))

# =============================================================================
# 3. HELPERS
# =============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x
metric_block <- function(y, prob, prefix) {
  pred <- factor(class_levels[max.col(prob, ties.method = "first")], levels = class_levels)
  m <- calc_overall_metrics(factor(y, levels = class_levels), pred, prob, class_levels)
  names(m) <- paste0(prefix, names(m))
  m
}
class_block <- function(y, prob, meta) {
  pred <- factor(class_levels[max.col(prob, ties.method = "first")], levels = class_levels)
  cm <- calc_class_metrics(factor(y, levels = class_levels), pred, class_levels)
  cm$PhaseLabel <- unname(phase_labels[cm$Class])
  cbind(meta[rep(1L, nrow(cm)), , drop = FALSE], cm)
}
scalar_criteria <- function(fit) {
  if (is.null(fit$criteria)) return(data.frame(row.names = 1L))
  z <- Filter(function(v) is.atomic(v) && length(v) == 1L, fit$criteria)
  as.data.frame(z, stringsAsFactors = FALSE)
}
file_key <- function(design, fold, model, G) sprintf("%s_fold%02d_%s_G%02d.rds", design, fold, gsub("[^A-Za-z0-9]", "", model), G)
run_file <- function(design, fold, model, G) file.path(output_dir, "runs", file_key(design, fold, model, G))
err_file <- function(design, fold, model, G) file.path(output_dir, "errors", file_key(design, fold, model, G))
keep_cols <- function(X) apply(X, 2, function(v) stats::sd(v) > 0)
# Model name = <family>-<S|D>[-<variant>]
parse_model <- function(model) {
  m <- regmatches(model, regexec("^(.*)-([SD])(-(.*))?$", model))[[1]]
  list(family = m[2], spec = m[3], variant = if (nzchar(m[5])) m[5] else "")
}
design_seed <- function(design, fold) outer_seed + 100L * fold + match(design, designs_all)

base_control <- list(type_multinomial = "ungrouped", alpha_grid = alpha_grid, min_units = min_units,
                     min_per_class = min_per_class, k_neighbors = k_neighbors, weight_type = "binary",
                     lambda_rule = "lambda.min", lambda_scale = "sum", max_iter = prof$max_iter,
                     tiny_movement_max_units = 0, tiny_movement_rate_tol = 0, tiny_movement_revert = FALSE,
                     coef_init_k = prof$coef_k, coef_init_anchors = prof$coef_anchors,
                     coef_init_alpha = 0.1, coef_init_lambda = 0.0625, verbose = show_progress)
control_for <- function(family, variant = "") {
  extra <- switch(family,
    "SCMR-EN" = if (variant == "phi1") list(phi_update = "fixed", phi = 1, n_starts = prof$n_starts) else
      list(phi_update = "pl", n_starts = prof$n_starts),
    "TwoStage-EN" = list(update_memberships = FALSE, n_starts = 1L, init_method = "kmeans"),
    list())
  if (variant == "unit") extra <- c(extra, list(min_units = 2L * min_units, coef_init_k = prof$coef_k_unit))
  do.call(scmr_control, utils::modifyList(base_control, extra))
}
fit_lightgbm <- function(X, y, seed) {
  set.seed(seed)
  ds <- lightgbm::lgb.Dataset(X, label = as.integer(y) - 1L)
  lightgbm::lgb.train(params = list(objective = "multiclass", num_class = nlevels(y), learning_rate = 0.05,
                                    num_leaves = 63, min_data_in_leaf = 50, feature_fraction = 0.8,
                                    bagging_fraction = 0.8, bagging_freq = 1, verbose = -1, seed = seed),
                      data = ds, nrounds = 400)
}
predict_lightgbm <- function(model, X) {
  p <- predict(model, X)
  if (is.null(dim(p))) p <- matrix(p, ncol = length(class_levels), byrow = TRUE)
  colnames(p) <- class_levels
  p
}
X_of <- function(spec, rows, lag_filled = FALSE) {
  if (spec == "S") X_base[rows, , drop = FALSE] else if (lag_filled) X_lagfill[rows, , drop = FALSE] else X_full[rows, , drop = FALSE]
}

# Fitted v3c runs (dynamic models, segment labels, identical settings) reused
# for E1/E3 when the tuned local penalty equals the one v3c used.
v3c_names <- c("Global-EN-D" = "Global-EN", "SCMR-EN-D" = "DSCMR-EN", "TwoStage-EN-D" = "TwoStage-EN",
               "SCMR-EN-D-phi1" = "DSCMR-EN-phi1", "SCR-original-D" = "SCR-original")
reuse_from_v3c <- function(design, fold, model, G, ctx) {
  if (!isTRUE(ctx$v3c_ok) || !model %in% names(v3c_names)) return(FALSE)
  old <- file.path(reuse_v3c, "runs", file_key(design, fold, v3c_names[[model]], G))
  if (!file.exists(old)) return(FALSE)
  b <- readRDS(old)
  pm <- parse_model(model)
  ren <- function(d) {
    if (is.null(d) || !nrow(d)) return(d)
    d$Model <- model; d$Spec <- pm$spec; d$Family <- pm$family; d$Variant <- pm$variant; d$Level <- NULL
    d
  }
  b$result <- ren(b$result); b$class <- ren(b$class); b$membership <- ren(b$membership)
  saveRDS(b, run_file(design, fold, model, G))
  cat(sprintf("  %-18s G=%d | reused from v3c\n", model, G))
  TRUE
}

# Fit one method (or reuse), evaluate in the design's prediction mode, save.
run_model <- function(design, fold, model, G, ctx) {
  path <- run_file(design, fold, model, G)
  if (!file.exists(path)) reuse_from_v3c(design, fold, model, G, ctx)
  if (file.exists(path)) return(invisible(readRDS(path)))
  pm <- parse_model(model)
  spec <- pm$spec; fam <- pm$family; var <- pm$variant
  sc <- ctx$spec[[spec]]
  train <- ctx$train; test <- ctx$test; seed <- ctx$seed
  meta <- data.frame(Design = design, Fold = fold, Seed = seed, Model = model, Family = fam, Spec = spec,
                     Variant = var, G = G)
  Xtr <- X_of(spec, train)
  kc <- keep_cols(Xtr); if (spec == "D") kc[lag_cols[lag_cols %in% names(kc)]] <- TRUE
  cols <- colnames(Xtr)[kc]
  Xtr <- Xtr[, cols, drop = FALSE]
  ytr <- droplevels(dat$phase[train])
  if (!identical(levels(ytr), class_levels)) stop("A class is absent from training rows.")
  unit_lab <- var == "unit"
  lab <- if (unit_lab) sub_id else seg_id
  lab_xy <- if (unit_lab) coords_all else seg_xy
  t0 <- proc.time()[["elapsed"]]
  fit <- switch(fam,
    "Global-EN" = sc$global,
    "LightGBM" = fit_lightgbm(cbind(Xtr, XY[train, ]), ytr, seed),
    "GW-EN" = fit_gw_multinom_en(Xtr, ytr, sub_id[train], coords_all[train, ], k_grid = prof$gw_k,
                                 alpha = sc$alpha, lambda = sc$lambda, anchors = prof$gw_anchors, seed = seed),
    "SCR-original" = fit_scr_original(Xtr, ytr, seg_id[train], seg_xy[train, ], G = G, phi = 1,
                                      k_neighbors = 5, maxitr = scr_maxitr, seed = seed),
    "Regency-EN" = fit_scmr(Xtr, ytr, model = "fixed_clusters", cluster = dat$regency[train],
                            unit_id = seg_id[train], alpha = sc$alpha, lambda = sc$lambda,
                            control = control_for(fam), seed = seed),
    fit_scmr(Xtr, ytr, model = "scmr", G = G, unit_id = lab[train], coords = lab_xy[train, ],
             alpha = sc$alpha, lambda = sc$lambda, control = control_for(fam, var), seed = seed,
             init_features = if (unit_lab) sc$features_unit else sc$features))
  fit_sec <- proc.time()[["elapsed"]] - t0
  te_lab <- lab[test]; te_xy <- lab_xy[test, , drop = FALSE]
  prob_fun <- function(x, idx) {
    if (fam == "LightGBM") return(predict_lightgbm(fit, cbind(x, XY[test[idx], , drop = FALSE])))
    if (fam == "GW-EN") return(predict(fit, x, coords_all[test[idx], , drop = FALSE]))
    if (fam == "Regency-EN") return(predict(fit, x, cluster = dat$regency[test[idx]]))
    if (fit$model == "global" || fit$G == 1L) return(predict(fit, x))
    predict(fit, x, new_unit_id = te_lab[idx], new_coords = te_xy[idx, , drop = FALSE], membership = membership_rule)
  }
  preds <- list()
  if (ctx$mode == "direct") {
    preds$direct <- list(rows = seq_along(test), prob = prob_fun(X_of(spec, test)[, cols, drop = FALSE], seq_along(test)))
  } else if (ctx$mode == "lag_observed") {
    ok <- which(lagd$available[test])
    preds$lag_observed <- list(rows = ok, prob = prob_fun(X_full[test[ok], cols, drop = FALSE], ok))
  } else if (ctx$mode == "impute") {
    p <- scmr_impute_waves(fit, X_lagfill[test, cols, drop = FALSE], ctx$y_obs, sub_id[test], dat$time[test],
                           te_xy, lag_cols, membership = membership_rule, cluster_id = te_lab)
    rows <- which(is.na(ctx$y_obs))
    preds$impute <- list(rows = rows, prob = p[rows, , drop = FALSE])
  }
  results <- list(); classes <- list(); rowpred <- list()
  for (md in names(preds)) {
    rows <- preds[[md]]$rows; pmat <- preds[[md]]$prob
    y <- dat$phase[test[rows]]
    m <- cbind(meta, Mode = md, FitSec = fit_sec, metric_block(y, pmat, "Test"))
    if (inherits(fit, "scmr_fit")) {
      m <- cbind(m, data.frame(Phi = fit$phi %||% NA_real_, Converged = fit$converged %||% NA,
                               Iterations = if (NROW(fit$diagnostics$iterations)) max(fit$diagnostics$iterations$Iter) else 0L),
                 scalar_criteria(fit))
    }
    results[[md]] <- m
    classes[[md]] <- cbind(class_block(y, pmat, meta), Mode = md)
    rowpred[[md]] <- data.frame(Row = test[rows], Mode = md, Pred = class_levels[max.col(pmat, ties.method = "first")],
                                pmat, check.names = FALSE)
  }
  is_clustered <- inherits(fit, "scmr_fit") && identical(fit$model, "scmr") && fit$G > 1L
  bundle <- list(result = do.call(rbind, results), class = do.call(rbind, classes),
                 predictions = do.call(rbind, rowpred),
                 membership = if (is_clustered) data.frame(meta[rep(1L, length(fit$group_unit)), ],
                   Unit = names(fit$group_unit), Cluster = as.integer(fit$group_unit), row.names = NULL) else NULL,
                 iterations = if (inherits(fit, "scmr_fit")) fit$diagnostics$iterations else NULL)
  saveRDS(bundle, path)
  if (file.exists(err_file(design, fold, model, G))) unlink(err_file(design, fold, model, G))
  r1 <- bundle$result[1, ]
  cat(sprintf("  %-18s G=%d | %s kappa %.4f | %.1f min\n", model, G, r1$Mode, r1$TestKappa, fit_sec / 60))
  invisible(bundle)
}
safe_run <- function(design, fold, model, G, ctx) {
  if (!file.exists(run_file(design, fold, model, G))) cat(sprintf("  [%s] start %s G=%d\n", format(Sys.time(), "%H:%M"), model, G))
  calls <- NULL
  tryCatch(withCallingHandlers(run_model(design, fold, model, G, ctx), error = function(e) {
    calls <<- vapply(sys.calls(), function(cl) paste(deparse(cl, nlines = 1L), collapse = ""), character(1))
  }), error = function(e) {
    saveRDS(list(message = conditionMessage(e), calls = calls), err_file(design, fold, model, G))
    cat("  ", model, "G =", G, "FAILED:", conditionMessage(e), "\n")
    keep <- calls[!grepl("^(tryCatch|withCallingHandlers|doTryCatch|tryCatchList|tryCatchOne|simpleError|stop|h\\()", calls)]
    cat("    call stack (innermost last):\n", paste0("      ", utils::tail(substr(keep, 1, 110), 12), collapse = "\n"), "\n")
    NULL
  })
}

# Global model of a specification (CV over alpha_grid, lambda.min), cached.
global_for <- function(design, fold, spec, train, seed) {
  name <- if (spec == "D") "global" else "globalS"
  gfile <- file.path(output_dir, "cache", sprintf("%s_fold%02d_%s.rds", design, fold, name))
  if (!file.exists(gfile) && spec == "D") {
    old <- file.path(reuse_cache_dirs, basename(gfile)); old <- old[file.exists(old)]
    if (length(old)) file.copy(old[1], gfile)
  }
  if (file.exists(gfile)) return(readRDS(gfile))
  cat(sprintf("  [%s] global cross-validation (spec %s)\n", format(Sys.time(), "%H:%M"), spec))
  Xtr <- X_of(spec, train)
  kc <- keep_cols(Xtr); if (spec == "D") kc[lag_cols] <- TRUE
  g <- fit_scmr(Xtr[, kc, drop = FALSE], droplevels(dat$phase[train]), model = "global",
                unit_id = sub_id[train], control = control_for("Global-EN"), nfolds = 5L,
                lambda_rule = "lambda.min", seed = seed)
  saveRDS(g, gfile)
  g
}

# Inner validation of the local penalty multiplier: one fold of whole training
# segments held out; a two-stage fit (G = tune_G) on the rest is scored by the
# log-loss of the design's prediction mode on the held-out segments.
tune_local <- function(design, fold, spec, train, lam_base, alp, seed, mode) {
  tfile <- file.path(output_dir, "cache", sprintf("%s_fold%02d_tuning%s.rds", design, fold, spec))
  if (file.exists(tfile)) { tuned <- readRDS(tfile); print(tuned$table, row.names = FALSE, digits = 4); return(tuned) }
  cat(sprintf("  [%s] inner validation of the local penalty (spec %s)\n", format(Sys.time(), "%H:%M"), spec))
  inner <- scmr_block_folds(dat$id_segmen[train], dat$regency[train], nfolds = 5, seed = seed + 17L)
  fit_rows <- train[inner != 1L]
  val_rows <- train[inner == 1L]
  if (mode != "direct") val_rows <- val_rows[lagd$available[val_rows]]
  Xfit <- X_of(spec, fit_rows)
  kc <- keep_cols(Xfit); if (spec == "D") kc[lag_cols] <- TRUE
  cols <- colnames(Xfit)[kc]
  yfit <- droplevels(dat$phase[fit_rows])
  tab <- data.frame()
  for (m in tune_multipliers) {
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(fit_scmr(Xfit[, cols, drop = FALSE], yfit, G = tune_G, unit_id = seg_id[fit_rows],
                             coords = seg_xy[fit_rows, ], alpha = alp, lambda = lam_base * m,
                             control = utils::modifyList(control_for("TwoStage-EN"), list(verbose = FALSE)),
                             seed = seed), error = function(e) e)
    if (inherits(fit, "error")) { cat("    multiplier", m, "failed:", conditionMessage(fit), "\n"); next }
    p <- predict(fit, X_of(spec, val_rows)[, cols, drop = FALSE], new_unit_id = seg_id[val_rows],
                 new_coords = seg_xy[val_rows, , drop = FALSE], membership = membership_rule)
    mb <- metric_block(dat$phase[val_rows], p, "")
    tab <- rbind(tab, data.frame(Multiplier = m, Lambda = lam_base * m, LogLoss = mb$LogLoss, Kappa = mb$Kappa,
                                 Minutes = (proc.time()[["elapsed"]] - t0) / 60))
  }
  if (!nrow(tab)) stop("Inner validation failed for every multiplier.")
  tuned <- list(multiplier = tab$Multiplier[which.min(tab$LogLoss)], table = tab)
  print(tab, row.names = FALSE, digits = 4)
  saveRDS(tuned, tfile)
  tuned
}

local_features <- function(design, fold, spec, train, global, seed, unit = FALSE) {
  name <- paste0(if (spec == "D") "features" else "featuresS", if (unit) "_unit" else "")
  ffile <- file.path(output_dir, "cache", sprintf("%s_fold%02d_%s.rds", design, fold, name))
  if (!file.exists(ffile) && spec == "D" && !unit) {
    old <- file.path(reuse_v3c, "cache", basename(ffile))
    if (file.exists(old)) file.copy(old, ffile)
  }
  if (file.exists(ffile)) return(readRDS(ffile))
  cat(sprintf("  [%s] local coefficient features (spec %s%s)\n", format(Sys.time(), "%H:%M"), spec, if (unit) ", sub-segments" else ""))
  Xtr <- X_of(spec, train)[, global$x_colnames, drop = FALSE]
  f <- scmr_local_coefficients(Xtr, droplevels(dat$phase[train]), if (unit) sub_id[train] else seg_id[train],
                               if (unit) coords_all[train, ] else seg_xy[train, ],
                               k = if (unit) prof$coef_k_unit else prof$coef_k, alpha = 0.1, lambda = 0.0625,
                               anchors = prof$coef_anchors, type_multinomial = "ungrouped", seed = seed)
  saveRDS(f, ffile)
  f
}

# G with the smallest PLIC-AIC among the SCMR-EN fits (G > 1) of a split.
selected_G <- function(design, fold, model, G_grid) {
  Gs <- setdiff(G_grid, 1L)
  crit <- vapply(Gs, function(G) {
    f <- run_file(design, fold, model, G)
    if (!file.exists(f)) return(NA_real_)
    v <- readRDS(f)$result$CriterionPLIC_AIC
    if (is.null(v)) NA_real_ else v[1]
  }, numeric(1))
  if (all(is.na(crit))) NA_integer_ else Gs[which.min(crit)]
}

run_split <- function(design, fold, train, test, spec, mode, G_grid, y_obs = NULL) {
  seed <- design_seed(design, fold)
  banner(design, " | fold ", fold, " | spec ", spec, " | train rows ", length(train), " | test rows ", length(test))
  global <- global_for(design, fold, spec, train, seed)
  alp <- global$alpha
  cvfit <- global$fits[[1]]$fit
  lam_base <- if (inherits(cvfit, "cv.glmnet")) cvfit[[local_lambda_rule]] else global$lambda_used[1]
  tuned <- tune_local(design, fold, spec, train, lam_base, alp, seed, mode)
  lam <- lam_base * tuned$multiplier
  # v3c runs are reusable when v3c used the same local penalty.
  v3c_ok <- FALSE
  old_t <- file.path(reuse_v3c, "cache", sprintf("%s_fold%02d_tuning.rds", design, fold))
  if (spec == "D" && file.exists(old_t)) {
    ot <- readRDS(old_t)
    v3c_ok <- isTRUE(all.equal(ot$multiplier, tuned$multiplier)) && identical(ot$rule, membership_rule)
  }
  cat("  Global-EN alpha", alp, "| lambda.min", signif(global$lambda_used[1], 3), "| local lambda", signif(lam, 3),
      "(", local_lambda_rule, "x", tuned$multiplier, ") | reuse v3c runs:", v3c_ok, "\n")
  features <- if (any(G_grid > 1L)) local_features(design, fold, spec, train, global, seed) else NULL
  ctx <- list(train = train, test = test, seed = seed, mode = mode, y_obs = y_obs, v3c_ok = v3c_ok, spec = list())
  ctx$spec[[spec]] <- list(global = global, alpha = alp, lambda = lam, features = features, features_unit = NULL)
  M <- function(f, v = "") paste0(f, "-", spec, if (nzchar(v)) paste0("-", v) else "")
  run <- function(model, G) safe_run(design, fold, model, G, ctx)
  run(M("Global-EN"), 1L)
  for (G in setdiff(G_grid, 1L)) { run(M("SCMR-EN"), G); run(M("TwoStage-EN"), G) }
  if (mode != "impute") run(M("Regency-EN"), 3L)
  G_sel <- selected_G(design, fold, M("SCMR-EN"), G_grid)
  if (!is.na(G_sel)) {
    cat("  PLIC-AIC selects G =", G_sel, "for the ablations\n")
    if (mode != "impute") run(M("SCMR-EN", "phi1"), G_sel)
    if (design == "E2_mapping") {
      ctx$spec[[spec]]$features_unit <- local_features(design, fold, spec, train, global, seed, unit = TRUE)
      run(M("SCMR-EN", "unit"), G_sel)
    }
  }
  if (use_scr_original) for (G in setdiff(G_grid, 1L)) run(M("SCR-original"), G)
  if (mode != "impute") {
    if (use_gw) run(M("GW-EN"), 0L)
    if (use_lightgbm) run(M("LightGBM"), 0L)
  }
}

# =============================================================================
# 4. DESIGNS
# =============================================================================
if ("E2_mapping" %in% designs_to_run) {
  folds <- scmr_block_folds(dat$id_segmen, dat$regency, nfolds = n_folds, seed = outer_seed)
  for (f in prof$folds) run_split("E2_mapping", f, which(folds != f), which(folds == f), "S", "direct", prof$G_grid)
}
if ("E1_nowcast" %in% designs_to_run) {
  T0 <- min(train_months_E1, max(dat$time) - 1L)
  run_split("E1_nowcast", 1L, which(dat$time <= T0 & lagd$available), which(dat$time > T0 & lagd$available),
            "D", "lag_observed", prof$G_grid)
}
if ("E3_imputation" %in% designs_to_run) {
  set.seed(outer_seed + 3L)
  cand <- which(dat$time > min(dat$time))
  masked <- sort(sample(cand, round(mask_prop_E3 * length(cand))))
  prev_row <- match(paste(dat$id_subsegmen, dat$time - 1L), paste(dat$id_subsegmen, dat$time))
  train <- which(lagd$available & !(seq_len(nrow(dat)) %in% masked) & !(prev_row %in% masked))
  y_obs <- as.character(dat$phase); y_obs[masked] <- NA
  run_split("E3_imputation", 1L, train, seq_len(nrow(dat)), "D", "impute", prof$G_grid, y_obs = y_obs)
  # Carry-forward baseline: last observed phase of the same sub-segment.
  cf <- character(length(masked))
  for (j in seq_along(masked)) {
    r <- masked[j]; u <- dat$id_subsegmen[r]
    prev <- which(dat$id_subsegmen == u & dat$time < dat$time[r] & !is.na(y_obs))
    cf[j] <- if (length(prev)) y_obs[prev[which.max(dat$time[prev])]] else names(which.max(table(y_obs)))
  }
  pm <- matrix(0, length(masked), length(class_levels), dimnames = list(NULL, class_levels))
  pm[cbind(seq_along(masked), match(cf, class_levels))] <- 1
  pm <- 0.98 * pm + 0.02 / length(class_levels)
  cfres <- cbind(data.frame(Design = "E3_imputation", Fold = 1L, Seed = NA, Model = "CarryForward", Family = "CarryForward",
                            Spec = "D", Variant = "", G = 0L, Mode = "impute", FitSec = 0),
                 metric_block(dat$phase[masked], pm, "Test"))
  saveRDS(list(result = cfres), run_file("E3_imputation", 1L, "CarryForward", 0L))
}

# =============================================================================
# 5. REPORTS
# =============================================================================
banner("REPORTS")
bind_fill <- function(rows) {
  rows <- Filter(function(z) is.data.frame(z) && nrow(z) > 0L, rows)
  if (!length(rows)) return(data.frame())
  cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(z) { for (nm in setdiff(cols, names(z))) z[[nm]] <- NA; z[, cols] }))
}
paths <- list.files(file.path(output_dir, "runs"), pattern = "\\.rds$", full.names = TRUE)
bundles <- lapply(paths, readRDS)
results <- bind_fill(lapply(bundles, `[[`, "result"))
utils::write.csv(results, file.path(output_dir, "results.csv"), row.names = FALSE)
utils::write.csv(bind_fill(lapply(bundles, `[[`, "class")), file.path(output_dir, "eval_per_class.csv"), row.names = FALSE)
utils::write.csv(bind_fill(lapply(bundles, `[[`, "membership")), file.path(output_dir, "memberships.csv"), row.names = FALSE)
preds <- bind_fill(lapply(seq_along(bundles), function(i) {
  p <- bundles[[i]]$predictions
  if (is.null(p)) return(NULL)
  cbind(bundles[[i]]$result[rep(1L, nrow(p)), c("Design", "Fold", "Model", "G")], p,
        Unit = dat$id_subsegmen[p$Row], Segment = dat$id_segmen[p$Row], Regency = dat$regency[p$Row],
        Time = dat$time[p$Row], Phase = as.character(dat$phase[p$Row]))
}))
utils::write.csv(preds, file.path(output_dir, "row_predictions.csv"), row.names = FALSE)
utils::write.csv(bind_fill(lapply(bundles, function(b) {
  it <- b$iterations
  if (is.null(it) || !nrow(it)) return(NULL)
  cbind(b$result[rep(1L, nrow(it)), c("Design", "Fold", "Model", "G")], it)
})), file.path(output_dir, "iterations.csv"), row.names = FALSE)
if (nrow(results)) {
  perf <- stats::aggregate(results[, c("TestKappa", "TestAccuracy", "TestMacroF1", "TestLogLoss")],
                           results[, c("Design", "Model", "G", "Mode")], function(v) mean(v, na.rm = TRUE))
  utils::write.csv(perf, file.path(output_dir, "performance.csv"), row.names = FALSE)
  print(perf[order(perf$Design, perf$Mode, -perf$TestKappa), ], digits = 3, row.names = FALSE)
  # Methods with their own choice of G, per fold: SCMR-EN by PLIC-AIC, SCR-original
  # by its original BIC, TwoStage-EN at the G of SCMR-EN; paired with Global-EN.
  chosen <- do.call(rbind, lapply(split(results, list(results$Design, results$Fold), drop = TRUE), function(d) {
    pick <- function(model, crit) {
      z <- d[d$Model == model & d$G > 1, ]
      if (!nrow(z) || !crit %in% names(z) || all(!is.finite(z[[crit]]))) return(NULL)
      z[which.min(z[[crit]]), ]
    }
    out <- list()
    for (sp in unique(d$Spec)) {
      s <- pick(paste0("SCMR-EN-", sp), "CriterionPLIC_AIC")
      if (!is.null(s)) {
        out$s <- transform(s, Choice = "PLIC-AIC")
        tw <- d[d$Model == paste0("TwoStage-EN-", sp) & d$G == s$G, ]
        if (nrow(tw)) out$t <- transform(tw[1, ], Choice = "G of SCMR-EN")
      }
      r <- pick(paste0("SCR-original-", sp), "CriterionSCR_BIC_original")
      if (!is.null(r)) out$r <- transform(r, Choice = "SCR BIC")
    }
    others <- d[!grepl("^(SCMR-EN-[SD]|TwoStage-EN|SCR-original)", d$Model) | grepl("-(phi1|unit)$", d$Model), ]
    if (nrow(others)) out$o <- transform(others, Choice = "fixed")
    bind_fill(out)
  }))
  if (!is.null(chosen) && nrow(chosen)) {
    utils::write.csv(chosen, file.path(output_dir, "chosen_models.csv"), row.names = FALSE)
    summ <- stats::aggregate(chosen[, c("TestKappa", "TestMacroF1", "TestLogLoss")], chosen[, c("Design", "Mode", "Model", "Choice")],
                             function(v) mean(v, na.rm = TRUE))
    nf <- stats::aggregate(list(Folds = chosen$Fold), chosen[, c("Design", "Mode", "Model", "Choice")], length)
    summ <- merge(summ, nf)
    cat("\nMethods with their own G (mean over folds):\n")
    print(summ[order(summ$Design, -summ$TestKappa), ], digits = 3, row.names = FALSE)
    # Paired differences against Global-EN of the same specification, per fold.
    glob <- chosen[grepl("^Global-EN-", chosen$Model), c("Design", "Fold", "TestKappa", "TestLogLoss")]
    names(glob)[3:4] <- c("GlobalKappa", "GlobalLogLoss")
    pd <- merge(chosen, glob)
    pd$dKappa <- pd$TestKappa - pd$GlobalKappa; pd$dLogLoss <- pd$TestLogLoss - pd$GlobalLogLoss
    pds <- do.call(rbind, lapply(split(pd, list(pd$Design, pd$Model), drop = TRUE), function(z) data.frame(
      Design = z$Design[1], Model = z$Model[1], Folds = nrow(z), dKappa = mean(z$dKappa),
      SE_dKappa = if (nrow(z) > 1) stats::sd(z$dKappa) / sqrt(nrow(z)) else NA_real_,
      dLogLoss = mean(z$dLogLoss), SE_dLogLoss = if (nrow(z) > 1) stats::sd(z$dLogLoss) / sqrt(nrow(z)) else NA_real_)))
    utils::write.csv(pds, file.path(output_dir, "paired_vs_global.csv"), row.names = FALSE)
    cat("\nPaired difference to Global-EN (mean over folds, SE across folds):\n")
    print(pds[order(pds$Design, -pds$dKappa), ], digits = 3, row.names = FALSE)
  }
}
banner("DONE")
cat("Outputs in:", output_dir, "| failed runs:", length(list.files(file.path(output_dir, "errors"))), "\n")
