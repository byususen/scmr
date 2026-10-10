# =============================================================================
# EMPIRICAL ARTICLE CODE v3 - dynamic SCMR-EN (DSCMR-EN) on KSA + Sentinel-1/2
#
# Requires scmr >= 0.5.0 (install from the scmr_pkg folder next to this file).
#
# Designs (see the research design document, section 7):
#   E2_mapping    MAIN. 5-fold CV over whole KSA segments, stratified by regency.
#                 Test units are new locations; their previous phase is NOT used:
#                 probabilities come from the exact forward filter
#                 (scmr_filter_predict) or, for comparison, plug-in of the
#                 previous predicted phase.
#   E1_nowcast    Temporal hold-out: train on months 1..18, predict months 19..24
#                 at the survey points with the observed previous phase.
#   E3_imputation 15% of monthly records hidden at random; imputed by
#                 forward-backward smoothing (scmr_impute_waves) and compared
#                 with carry-forward of the last observed phase.
# Models: Global-EN, DSCMR-EN (estimated phi), DSCMR-EN-phi1 (ablation),
#         TwoStage-EN, SCMR-EN-static (no previous phase), GW-EN (optional),
#         LightGBM (optional, if the lightgbm package is installed),
#         SCR-original (benchmark port of Sugasawa & Murakami's reference code).
# Cluster labels are fitted per survey segment by default (SCMR_MEMBERSHIP_LEVEL).
# Resume-safe: one RDS per (design, fold, model, G) under output_dir/runs.
# Settings can be overridden with environment variables SCMR_SYNTAX_DIR,
# SCMR_DATA_FILE, SCMR_OUTPUT_DIR and SCMR_PROFILE ("pilot", "main", "ci").
# =============================================================================

# =============================================================================
# 0. SETTINGS
# =============================================================================
syntax_dir <- Sys.getenv("SCMR_SYNTAX_DIR", "C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax")
data_file <- Sys.getenv("SCMR_DATA_FILE", file.path(syntax_dir, "df_month_rec_final.csv"))
output_dir <- Sys.getenv("SCMR_OUTPUT_DIR", file.path(dirname(syntax_dir), "Results", "Real data v3c"))
# Cached global fits of an earlier run are reused when present (the global
# model does not depend on the membership level).
reuse_cache_dirs <- file.path(dirname(syntax_dir), "Results", c("Real data v3b", "Real data v3"), "cache")
# Spatial level of the cluster labels. "segment": one label per survey segment
# (its 3 x 3 sub-segments share it) with the Potts graph among segments;
# "subsegment": one label per sub-segment (runs v3 and v3b). In v3b the
# sub-segment graph linked almost only sub-segments of the same segment, so the
# estimated phi saturated and the Potts term gave no coherence between
# segments. Class sequences (lags, filter, imputation) stay per sub-segment.
membership_level <- Sys.getenv("SCMR_MEMBERSHIP_LEVEL", "segment")
run_profile <- Sys.getenv("SCMR_PROFILE", "pilot")   # "pilot" first, then "main"

designs_to_run <- c("E2_mapping", "E1_nowcast", "E3_imputation")
use_gw <- TRUE                   # geographically weighted comparator (slow on full data)
# Optional overrides (used to split the work into parallel cloud jobs):
#   SCMR_DESIGNS="E2_mapping"  SCMR_FOLDS="1,2"  SCMR_USE_GW="FALSE"
if (nzchar(Sys.getenv("SCMR_DESIGNS"))) designs_to_run <- strsplit(Sys.getenv("SCMR_DESIGNS"), ",")[[1]]
if (nzchar(Sys.getenv("SCMR_USE_GW"))) use_gw <- as.logical(Sys.getenv("SCMR_USE_GW"))
use_lightgbm <- requireNamespace("lightgbm", quietly = TRUE)
use_scr_original <- TRUE         # benchmark: original SCR algorithm (unpenalized; may be slow)
scr_maxitr <- 30L
if (nzchar(Sys.getenv("SCMR_USE_SCR"))) use_scr_original <- as.logical(Sys.getenv("SCMR_USE_SCR"))

profiles <- list(
  ci    = list(folds = 1L, G_grid = c(1L, 2L), n_starts = 2L, max_iter = 10L, gw_anchors = 30L,
               gw_k = c(10, 20), coef_k = 10L, coef_anchors = 60L, E3_G = 2L),
  # Local coefficient initialisation over about one segment (k = 12 units) with
  # six starts: in the panel pilot this raised the ARI of irregular regimes at
  # T = 12 from 0.23 to 0.78.
  pilot = list(folds = 1L, G_grid = c(1L, 3L, 6L), n_starts = 3L, max_iter = 20L, gw_anchors = 300L,
               gw_k = c(50, 100, 200), coef_k = 12L, coef_anchors = 600L, E3_G = 3L),
  main  = list(folds = 1:5, G_grid = 1:8, n_starts = 6L, max_iter = 30L, gw_anchors = 300L,
               gw_k = c(50, 100, 200, 400), coef_k = 12L, coef_anchors = 600L, E3_G = NA))
if (!run_profile %in% names(profiles)) stop("Unknown profile: ", run_profile)
if (!membership_level %in% c("segment", "subsegment")) stop("Unknown membership level: ", membership_level)
prof <- profiles[[run_profile]]
if (membership_level == "segment") prof$coef_k <- 3L   # about 27 sub-segments, as 12 sub-segments before
if (nzchar(Sys.getenv("SCMR_FOLDS"))) prof$folds <- as.integer(strsplit(Sys.getenv("SCMR_FOLDS"), ",")[[1]])
if (nzchar(Sys.getenv("SCMR_G_GRID"))) {
  prof$G_grid <- as.integer(strsplit(Sys.getenv("SCMR_G_GRID"), ",")[[1]])
  prof$E3_G <- NA
}
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
min_units <- if (membership_level == "segment") 5L else 10L   # labelled units per cluster
min_per_class <- 5L
k_neighbors <- 8L                # Potts / new-location neighbours among labelled units
# Penalty of the local (clustered and GW) models: the global CV "lambda.1se".
# On the KSA data lambda.min is tiny (about 4e-5): the previous phase nearly
# separates the classes inside clusters, glmnet then converges very slowly and
# the clustered fits take hours. Global-EN itself still uses lambda.min.
local_lambda_rule <- "lambda.1se"
# The local penalty is lambda.1se times a multiplier, and the cluster weights of
# new locations follow a membership rule; both are chosen by an inner
# validation on whole training segments (filter log-loss of a two-stage fit),
# before any test data are touched. The first pilot showed near-separation
# inside clusters (the previous phase is very predictive) at lambda.1se.
tune_multipliers <- c(1, 4, 16, 64)
tune_rules <- c("proportion", "potts")
tune_G <- 3L
show_progress <- TRUE            # print every membership iteration

# =============================================================================
# 1. PACKAGES
# =============================================================================
if (!requireNamespace("scmr", quietly = TRUE) || utils::packageVersion("scmr") < "0.5.0") {
  stop("Install scmr >= 0.5.0 first (install.packages(<scmr_pkg folder>, repos = NULL, type = 'source')).")
}
suppressPackageStartupMessages(library(scmr))
cat("scmr", as.character(utils::packageVersion("scmr")), "| profile", run_profile,
    "| membership level", membership_level, "| LightGBM", use_lightgbm, "\n")
dir.create(file.path(output_dir, "runs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "errors"), recursive = TRUE, showWarnings = FALSE)
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
# Labelled units (mem_id) and their coordinates: segment centroids or the
# sub-segment centroids.
if (membership_level == "segment") {
  mem_id <- dat$id_segmen
  seg_xy <- stats::aggregate(cbind(unit_lat, unit_lon) ~ id_segmen,
                             data = unique(dat[, c("id_segmen", "id_subsegmen", "unit_lat", "unit_lon")]), FUN = mean)
  mem_xy <- cbind(seg_xy$unit_lat[match(mem_id, seg_xy$id_segmen)], seg_xy$unit_lon[match(mem_id, seg_xy$id_segmen)])
} else {
  mem_id <- dat$id_subsegmen
  mem_xy <- coords_all
}

month_mm <- stats::model.matrix(~ factor(month, levels = 1:12) - 1, dat)
colnames(month_mm) <- paste0("month", 1:12)
X_base <- cbind(month_mm, as.matrix(dat[, x_vars]))
X_full <- cbind(X_base, lagd$x)                     # lag columns NA where unavailable
X_lagfill <- X_full; X_lagfill[is.na(X_lagfill)] <- 0
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
run_file <- function(design, fold, model, G) file.path(output_dir, "runs",
  sprintf("%s_fold%02d_%s_G%02d.rds", design, fold, gsub("[^A-Za-z0-9]", "", model), G))
err_file <- function(design, fold, model, G) file.path(output_dir, "errors",
  sprintf("%s_fold%02d_%s_G%02d.rds", design, fold, gsub("[^A-Za-z0-9]", "", model), G))
keep_cols <- function(X) apply(X, 2, function(v) stats::sd(v) > 0)

base_control <- list(type_multinomial = "ungrouped", alpha_grid = alpha_grid, min_units = min_units,
                     min_per_class = min_per_class, k_neighbors = k_neighbors, weight_type = "binary",
                     lambda_rule = "lambda.min", lambda_scale = "sum", max_iter = prof$max_iter,
                     tiny_movement_max_units = 0, tiny_movement_rate_tol = 0, tiny_movement_revert = FALSE,
                     coef_init_k = prof$coef_k, coef_init_anchors = prof$coef_anchors,
                     coef_init_alpha = 0.1, coef_init_lambda = 0.0625, verbose = show_progress)
control_for <- function(model) {
  extra <- switch(model,
    "Global-EN" = list(),
    "DSCMR-EN" = list(phi_update = "pl", n_starts = prof$n_starts),
    "SCMR-EN-static" = list(phi_update = "pl", n_starts = prof$n_starts),
    "DSCMR-EN-phi1" = list(phi_update = "fixed", phi = 1, n_starts = prof$n_starts),
    "TwoStage-EN" = list(update_memberships = FALSE, n_starts = 1L, init_method = "kmeans"))
  do.call(scmr_control, utils::modifyList(base_control, extra))
}

# Exact filter for scmr/global fits; plug-in of the previous predicted phase
# for any row-wise probability function.
set_lag <- function(x, k) {
  x[, lag_cols] <- 0
  x[cbind(seq_len(nrow(x)), match(lag_cols[k], colnames(x)))] <- 1
  x
}
plugin_predict <- function(prob_fun, x, unit, time, start_class) {
  out <- matrix(0, nrow(x), length(class_levels), dimnames = list(NULL, class_levels))
  cur <- stats::setNames(rep(start_class, length(unique(unit))), unique(unit))
  last_t <- stats::setNames(rep(-Inf, length(unique(unit))), unique(unit))
  for (t in sort(unique(time))) {
    idx <- which(time == t)
    gap <- last_t[unit[idx]] != t - 1
    cur[unit[idx][gap]] <- start_class
    p <- prob_fun(set_lag(x[idx, , drop = FALSE], cur[unit[idx]]), idx)
    out[idx, ] <- p
    cur[unit[idx]] <- max.col(p, ties.method = "first")
    last_t[unit[idx]] <- t
  }
  out
}
generic_filter <- function(prob_fun, x, unit, time, a0) {
  C <- length(class_levels)
  P <- lapply(seq_len(C), function(k) prob_fun(set_lag(x, rep(k, nrow(x))), seq_len(nrow(x))))
  out <- matrix(0, nrow(x), C, dimnames = list(NULL, class_levels))
  for (u in unique(unit)) {
    idx <- which(unit == u); idx <- idx[order(time[idx])]
    a <- a0; prev <- NA
    for (r in idx) {
      if (!is.na(prev) && time[r] != prev + 1) a <- a0
      a <- as.numeric(vapply(seq_len(C), function(c) sum(vapply(P, function(Pk) Pk[r, c], 0) * a), 0))
      a <- a / sum(a); out[r, ] <- a; prev <- time[r]
    }
  }
  out
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

# Fit one model (or reuse), evaluate under the design's prediction modes, save.
run_model <- function(design, fold, model, G, train, test, seed, global, lam, alp, features = NULL,
                      modes = c("filter", "plugin", "lag_observed"), y_obs_impute = NULL) {
  path <- run_file(design, fold, model, G)
  if (file.exists(path)) return(invisible(readRDS(path)))
  meta <- data.frame(Design = design, Fold = fold, Seed = seed, Model = model, G = G,
                     Level = if (model %in% c("GW-EN", "LightGBM", "Global-EN")) "none" else membership_level)
  static <- model == "SCMR-EN-static"
  cols_all <- if (static) colnames(X_base) else colnames(X_full)
  tr_rows <- train
  Xtr <- X_full[tr_rows, cols_all, drop = FALSE]
  kc <- keep_cols(Xtr); kc[lag_cols[lag_cols %in% names(kc)]] <- TRUE
  cols <- cols_all[kc]
  Xtr <- Xtr[, cols, drop = FALSE]
  ytr <- droplevels(dat$phase[tr_rows])
  if (!identical(levels(ytr), class_levels)) stop("A class is absent from training rows.")
  t0 <- proc.time()[["elapsed"]]
  fit <- if (model == "Global-EN") global else if (model == "LightGBM") {
    fit_lightgbm(Xtr, ytr, seed)
  } else if (model == "GW-EN") {
    fit_gw_multinom_en(Xtr, ytr, dat$id_subsegmen[tr_rows], coords_all[tr_rows, ], k_grid = prof$gw_k,
                       alpha = alp, lambda = lam, anchors = prof$gw_anchors, seed = seed)
  } else if (model == "SCR-original") {
    # Benchmark: multinomial port of the reference SCR code (unpenalized,
    # simultaneous updates, phi = 1, k-means start, 5 nearest neighbours).
    fit_scr_original(Xtr, ytr, mem_id[tr_rows], mem_xy[tr_rows, ], G = G, phi = 1,
                     k_neighbors = 5, maxitr = scr_maxitr, seed = seed)
  } else if (G == 1L && !static) {
    global
  } else {
    fit_scmr(Xtr, ytr, model = if (G == 1L) "global" else "scmr", G = G, unit_id = mem_id[tr_rows],
             coords = mem_xy[tr_rows, ], alpha = alp, lambda = lam, control = control_for(model),
             seed = seed, init_features = if (!static && G > 1L) features else NULL)
  }
  fit_sec <- proc.time()[["elapsed"]] - t0
  te_unit <- dat$id_subsegmen[test]; te_time <- dat$time[test]; te_xy <- coords_all[test, , drop = FALSE]
  te_mem <- mem_id[test]; te_mem_xy <- mem_xy[test, , drop = FALSE]
  Xte <- X_lagfill[test, cols, drop = FALSE]
  prob_fun <- function(x, idx) {
    if (inherits(fit, "lgb.Booster")) return(predict_lightgbm(fit, x))
    if (inherits(fit, "scmr_gw")) return(predict(fit, x, te_xy[idx, , drop = FALSE]))
    if (fit$model == "global") return(predict(fit, x))
    predict(fit, x, new_unit_id = te_mem[idx], new_coords = te_mem_xy[idx, , drop = FALSE], membership = new_membership)
  }
  a0 <- as.numeric(table(ytr)) / length(ytr)
  preds <- list()
  if (static) {
    preds$direct <- prob_fun(Xte, seq_along(test))
  } else {
    if ("lag_observed" %in% modes) {
      ok <- lagd$available[test]
      if (any(ok)) preds$lag_observed <- list(rows = which(ok), prob = prob_fun(X_full[test[ok], cols, drop = FALSE], which(ok)))
    }
    if ("filter" %in% modes) {
      preds$filter <- if (inherits(fit, c("scmr_fit"))) {
        scmr_filter_predict(fit, Xte, te_unit, te_time, te_mem_xy, lag_cols, membership = new_membership,
                            cluster_id = te_mem)
      } else generic_filter(prob_fun, Xte, te_unit, te_time, a0)
    }
    if ("plugin" %in% modes) preds$plugin <- plugin_predict(prob_fun, Xte, te_unit, te_time, which.max(a0))
    if ("impute" %in% modes) {
      preds$impute <- scmr_impute_waves(fit, Xte, y_obs_impute, te_unit, te_time, te_mem_xy, lag_cols,
                                        membership = new_membership, cluster_id = te_mem)
    }
  }
  results <- list(); classes <- list(); rowpred <- list()
  for (md in names(preds)) {
    p <- preds[[md]]
    rows <- if (is.list(p)) p$rows else seq_along(test)
    pm <- if (is.list(p)) p$prob else p
    if (md == "impute") { rows <- which(is.na(y_obs_impute)); pm <- pm[rows, , drop = FALSE] }
    y <- dat$phase[test[rows]]
    m <- cbind(meta, Mode = md, FitSec = fit_sec, metric_block(y, pm, "Test"))
    if (inherits(fit, "scmr_fit")) {
      m <- cbind(m, data.frame(Phi = fit$phi %||% NA_real_, Converged = fit$converged,
                               Iterations = if (nrow(fit$diagnostics$iterations)) max(fit$diagnostics$iterations$Iter) else 0L),
                 scalar_criteria(fit))
    }
    results[[md]] <- m
    classes[[md]] <- cbind(class_block(y, pm, meta), Mode = md)
    rowpred[[md]] <- data.frame(Row = test[rows], Mode = md, Pred = class_levels[max.col(pm, ties.method = "first")],
                                pm, check.names = FALSE)
  }
  bundle <- list(result = do.call(rbind, lapply(results, function(z) z)), class = do.call(rbind, classes),
                 predictions = do.call(rbind, rowpred),
                 membership = if (inherits(fit, "scmr_fit") && fit$G > 1L) data.frame(meta[rep(1L, length(fit$group_unit)), ],
                   Unit = names(fit$group_unit), Cluster = as.integer(fit$group_unit), row.names = NULL) else NULL,
                 coefficients = if (inherits(fit, "scmr_fit")) coef(fit) else NULL,
                 iterations = if (inherits(fit, "scmr_fit")) fit$diagnostics$iterations else NULL)
  saveRDS(bundle, path)
  if (file.exists(err_file(design, fold, model, G))) unlink(err_file(design, fold, model, G))
  r1 <- bundle$result[1, ]
  cat(sprintf("  %-15s G=%d | %s kappa %.4f | %.1f min\n", model, G, r1$Mode, r1$TestKappa, fit_sec / 60))
  invisible(bundle)
}
safe_run <- function(design, fold, model, G, ...) {
  if (!file.exists(run_file(design, fold, model, G))) {
    cat(sprintf("  [%s] start %s G=%d\n", format(Sys.time(), "%H:%M"), model, G))
  }
  calls <- NULL
  tryCatch(withCallingHandlers(run_model(design, fold, model, G, ...), error = function(e) {
    calls <<- vapply(sys.calls(), function(cl) paste(deparse(cl, nlines = 1L), collapse = ""), character(1))
  }), error = function(e) {
    saveRDS(list(message = conditionMessage(e), calls = calls), err_file(design, fold, model, G))
    cat("  ", model, "G =", G, "FAILED:", conditionMessage(e), "\n")
    keep <- calls[!grepl("^(tryCatch|withCallingHandlers|doTryCatch|tryCatchList|tryCatchOne|simpleError|stop|h\\()", calls)]
    cat("    call stack (innermost last):\n", paste0("      ", utils::tail(substr(keep, 1, 110), 12), collapse = "\n"), "\n")
    NULL
  })
}

global_fit <- function(train, seed) {
  Xtr <- X_full[train, , drop = FALSE]
  kc <- keep_cols(Xtr); kc[lag_cols] <- TRUE
  fit_scmr(Xtr[, kc, drop = FALSE], droplevels(dat$phase[train]), model = "global",
           unit_id = dat$id_subsegmen[train], control = control_for("Global-EN"),
           nfolds = 5L, lambda_rule = "lambda.min", seed = seed)
}

new_membership <- "proportion"
# Inner validation: one fold of whole training segments is held out; a
# two-stage fit (G = tune_G) on the rest is scored by the forward-filter
# log-loss on the held-out segments for every multiplier and membership rule.
tune_local <- function(design, fold, train, lam_base, alp, seed, cache_dir) {
  tfile <- file.path(cache_dir, sprintf("%s_fold%02d_tuning.rds", design, fold))
  if (file.exists(tfile)) {
    tuned <- readRDS(tfile)
    print(tuned$table, row.names = FALSE, digits = 4)
    return(tuned)
  }
  cat(sprintf("  [%s] inner validation of the local penalty and membership rule\n", format(Sys.time(), "%H:%M")))
  inner <- scmr_block_folds(dat$id_segmen[train], dat$regency[train], nfolds = 5, seed = seed + 17L)
  fit_rows <- train[inner != 1L]
  val_units <- unique(dat$id_subsegmen[train[inner == 1L]])
  val_rows <- which(dat$id_subsegmen %in% val_units)
  Xfit <- X_full[fit_rows, , drop = FALSE]
  kc <- keep_cols(Xfit); kc[lag_cols] <- TRUE
  cols <- colnames(X_full)[kc]
  yfit <- droplevels(dat$phase[fit_rows])
  tab <- data.frame()
  for (m in tune_multipliers) {
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(fit_scmr(Xfit[, cols, drop = FALSE], yfit, G = tune_G, unit_id = mem_id[fit_rows],
                             coords = mem_xy[fit_rows, ], alpha = alp, lambda = lam_base * m,
                             control = utils::modifyList(control_for("TwoStage-EN"), list(verbose = FALSE)),
                             seed = seed), error = function(e) e)
    if (inherits(fit, "error")) { cat("    multiplier", m, "failed:", conditionMessage(fit), "\n"); next }
    for (rule in tune_rules) {
      p <- scmr_filter_predict(fit, X_lagfill[val_rows, cols, drop = FALSE], dat$id_subsegmen[val_rows],
                               dat$time[val_rows], mem_xy[val_rows, , drop = FALSE], lag_cols, membership = rule,
                               cluster_id = mem_id[val_rows])
      mb <- metric_block(dat$phase[val_rows], p, "")
      tab <- rbind(tab, data.frame(Multiplier = m, Lambda = lam_base * m, Rule = rule,
                                   FilterLogLoss = mb$LogLoss, FilterKappa = mb$Kappa,
                                   TrainPenalizedObjective = fit$criteria$PenalizedObjective,
                                   Minutes = (proc.time()[["elapsed"]] - t0) / 60))
    }
  }
  if (!nrow(tab)) stop("Inner validation failed for every multiplier.")
  best <- tab[which.min(tab$FilterLogLoss), ]
  tuned <- list(multiplier = best$Multiplier, rule = as.character(best$Rule), table = tab)
  print(tab, row.names = FALSE, digits = 4)
  saveRDS(tuned, tfile)
  tuned
}

run_split <- function(design, fold, train, test, modes, G_grid, y_obs_impute = NULL) {
  seed <- outer_seed + 100L * fold + match(design, designs_to_run)
  banner(design, " | fold ", fold, " | train rows ", length(train), " | test rows ", length(test))
  # The global fit and the local coefficient features are cached per split, so
  # a rerun after an interruption does not repeat them.
  cache_dir <- file.path(output_dir, "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  gfile <- file.path(cache_dir, sprintf("%s_fold%02d_global.rds", design, fold))
  old_g <- file.path(reuse_cache_dirs, basename(gfile))
  old_g <- old_g[file.exists(old_g)]
  if (!file.exists(gfile) && length(old_g)) file.copy(old_g[1], gfile)
  if (file.exists(gfile)) global <- readRDS(gfile) else {
    cat(sprintf("  [%s] global cross-validation\n", format(Sys.time(), "%H:%M")))
    global <- global_fit(train, seed)
    saveRDS(global, gfile)
  }
  alp <- global$alpha
  cvfit <- global$fits[[1]]$fit
  lam_base <- if (inherits(cvfit, "cv.glmnet")) cvfit[[local_lambda_rule]] else global$lambda_used[1]
  tuned <- if (any(G_grid > 1L) || use_gw) tune_local(design, fold, train, lam_base, alp, seed, cache_dir) else
    list(multiplier = 1, rule = "proportion")
  lam <- lam_base * tuned$multiplier
  new_membership <<- tuned$rule
  cat("  Global-EN alpha", alp, "| lambda.min", signif(global$lambda_used[1], 3),
      "| local lambda", signif(lam, 3), "(", local_lambda_rule, "x", tuned$multiplier, ")",
      "| new-location membership", new_membership, "\n")
  features <- NULL
  if (any(G_grid > 1L)) {
    ffile <- file.path(cache_dir, sprintf("%s_fold%02d_features.rds", design, fold))
    if (file.exists(ffile)) features <- readRDS(ffile) else {
      cat(sprintf("  [%s] local coefficient features\n", format(Sys.time(), "%H:%M")))
      Xtr <- X_full[train, global$x_colnames, drop = FALSE]
      features <- scmr_local_coefficients(Xtr, droplevels(dat$phase[train]), mem_id[train],
                                          mem_xy[train, ], k = prof$coef_k, alpha = 0.1, lambda = 0.0625,
                                          anchors = prof$coef_anchors, type_multinomial = "ungrouped", seed = seed)
      saveRDS(features, ffile)
    }
  }
  args <- list(train = train, test = test, seed = seed, global = global, lam = lam, alp = alp,
               features = features, modes = modes, y_obs_impute = y_obs_impute)
  do.call(safe_run, c(list(design, fold, "Global-EN", 1L), args))
  for (G in setdiff(G_grid, 1L)) {
    do.call(safe_run, c(list(design, fold, "DSCMR-EN", G), args))
    if (design != "E3_imputation") do.call(safe_run, c(list(design, fold, "TwoStage-EN", G), args))
  }
  if (design == "E2_mapping") {
    for (G in setdiff(G_grid, 1L)) {
      do.call(safe_run, c(list(design, fold, "DSCMR-EN-phi1", G), args))
      do.call(safe_run, c(list(design, fold, "SCMR-EN-static", G), args))
    }
    if (use_gw) do.call(safe_run, c(list(design, fold, "GW-EN", 0L), args))
    if (use_lightgbm) do.call(safe_run, c(list(design, fold, "LightGBM", 0L), args))
    if (use_scr_original) for (G in setdiff(G_grid, 1L)) do.call(safe_run, c(list(design, fold, "SCR-original", G), args))
  }
  if (design == "E1_nowcast" && use_lightgbm) do.call(safe_run, c(list(design, fold, "LightGBM", 0L), args))
}

# =============================================================================
# 4. DESIGNS
# =============================================================================
if ("E2_mapping" %in% designs_to_run) {
  folds <- scmr_block_folds(dat$id_segmen, dat$regency, nfolds = n_folds, seed = outer_seed)
  for (f in prof$folds) {
    test <- which(folds == f)
    train <- which(folds != f & lagd$available)
    run_split("E2_mapping", f, train, test, c("filter", "plugin", "lag_observed"), prof$G_grid)
  }
}
if ("E1_nowcast" %in% designs_to_run) {
  T0 <- min(train_months_E1, max(dat$time) - 1L)
  train <- which(dat$time <= T0 & lagd$available)
  test <- which(dat$time > T0 & lagd$available)
  run_split("E1_nowcast", 1L, train, test, "lag_observed", prof$G_grid)
}
if ("E3_imputation" %in% designs_to_run) {
  set.seed(outer_seed + 3L)
  cand <- which(dat$time > min(dat$time))
  masked <- sort(sample(cand, round(mask_prop_E3 * length(cand))))
  prev_row <- match(paste(dat$id_subsegmen, dat$time - 1L), paste(dat$id_subsegmen, dat$time))
  train <- which(lagd$available & !(seq_len(nrow(dat)) %in% masked) & !(prev_row %in% masked))
  test <- seq_len(nrow(dat))
  y_obs <- as.character(dat$phase); y_obs[masked] <- NA
  G_E3 <- if (is.na(prof$E3_G)) prof$G_grid else c(1L, prof$E3_G)
  run_split("E3_imputation", 1L, train, test, "impute", G_E3, y_obs_impute = y_obs)
  # Carry-forward baseline: last observed phase of the same unit.
  cf <- character(length(masked))
  for (j in seq_along(masked)) {
    r <- masked[j]; u <- dat$id_subsegmen[r]
    prev <- which(dat$id_subsegmen == u & dat$time < dat$time[r] & !is.na(y_obs))
    cf[j] <- if (length(prev)) y_obs[prev[which.max(dat$time[prev])]] else names(which.max(table(y_obs)))
  }
  pm <- matrix(0, length(masked), length(class_levels), dimnames = list(NULL, class_levels))
  pm[cbind(seq_along(masked), match(cf, class_levels))] <- 1
  pm <- 0.98 * pm + 0.02 / length(class_levels)
  cfres <- cbind(data.frame(Design = "E3_imputation", Fold = 1L, Seed = NA, Model = "CarryForward", G = 0L,
                            Mode = "impute", FitSec = 0), metric_block(dat$phase[masked], pm, "Test"))
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
  crit <- intersect(c("CriterionPLIC_AIC", "CriterionPLIC_BIC", "CriterionSCR_AIC_effective"), names(results))
  sc <- results[results$Model %in% c("DSCMR-EN", "Global-EN") & results$Mode %in% c("filter", "lag_observed", "impute"), ]
  if (length(crit) && nrow(sc)) {
    cat("\nSelected G (DSCMR-EN grid, Global-EN as G = 1):\n")
    for (cr in crit) for (key in unique(paste(sc$Design, sc$Fold, sc$Mode))) {
      z <- sc[paste(sc$Design, sc$Fold, sc$Mode) == key & is.finite(sc[[cr]]), ]
      if (nrow(z)) cat(sprintf("  %-28s %-30s G = %d\n", cr, key, z$G[which.min(z[[cr]])]))
    }
  }
}
banner("DONE")
cat("Outputs in:", output_dir, "| failed runs:", length(list.files(file.path(output_dir, "errors"))), "\n")
