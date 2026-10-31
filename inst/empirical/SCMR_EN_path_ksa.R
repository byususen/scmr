# =============================================================================
# SPLIT-MERGE PATH ON THE KSA DATA (after SCMR_EN_empirical_article_v4.R)
#
# Fits SCMR-EN-S on all segments for G = 1..G_max with scmr_split_merge_path():
# each G starts from the best fit at G - 1 with one regime split along its block
# score vectors; a backward pass tries merges. Same tuning as the full-data
# interpretation fits (global CV of the static specification, inner validation
# of the local penalty, segment labels, Potts phi estimated), so the objectives
# can be compared with the multi-start fits in interpret/fits.
#
# Requires the scmr version on branch feature/size-aware-init:
#   remotes::install_github("byususen/scmr@feature/size-aware-init")
#
# Environment:
#   SCMR_PATH_GMAX   largest G (default 10)
#   SCMR_PATH_BASE   "none" (default) or "global" (shared-base regimes)
#   SCMR_PATH_WARM   phi warm-up sweeps (default 0)
#   SCMR_PATH_SPLIT  split candidates per G (default 2)
#   SCMR_PATH_MERGE  merge candidates per G in the backward pass (default 1)
#   plus the v4 variables (SCMR_SYNTAX_DIR, SCMR_OUTPUT_DIR, ...).
#
# Output in <output_dir>/path_<base>_w<warm>/ (aggregates only):
#   path_table.csv  one row per G: objective, phi, smallest regime, EDF, criteria
#   path_log.csv    every candidate fit: G, step, objective, seconds
#   path.rds        the path object (memberships and criteria per G; stays local)
# =============================================================================
env_chr <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
G_max <- as.integer(env_chr("SCMR_PATH_GMAX", "10"))
base_model <- env_chr("SCMR_PATH_BASE", "none")
warm <- as.integer(env_chr("SCMR_PATH_WARM", "0"))
n_split <- as.integer(env_chr("SCMR_PATH_SPLIT", "2"))
n_merge <- as.integer(env_chr("SCMR_PATH_MERGE", "1"))

syntax_dir0 <- Sys.getenv("SCMR_SYNTAX_DIR", "C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax")
v4_script <- Sys.getenv("SCMR_V4_SCRIPT", file.path(syntax_dir0, "SCMR_EN_empirical_article_v4.R"))
designs_env <- Sys.getenv("SCMR_DESIGNS")
Sys.setenv(SCMR_DESIGNS = "none")
source(v4_script)                       # data, design matrices and helpers only
Sys.setenv(SCMR_DESIGNS = designs_env)
if (!exists("scmr_split_merge_path")) stop("Install scmr from branch feature/size-aware-init first.")

path_dir <- file.path(output_dir, sprintf("path_%s_w%d", base_model, warm))
dir.create(path_dir, recursive = TRUE, showWarnings = FALSE)
banner(sprintf("SPLIT-MERGE PATH | G = 1..%d | base %s | warm-up %d", G_max, base_model, warm))
all_rows <- seq_len(nrow(dat))
seed_full <- outer_seed + 999L          # as in SCMR_EN_interpret_v4.R
global <- global_for("FULL", 0L, "S", all_rows, seed_full)
cvfit <- global$fits[[1]]$fit
lam_base <- if (inherits(cvfit, "cv.glmnet")) cvfit[[local_lambda_rule]] else global$lambda_used[1]
tuned <- tune_local("FULL", 0L, "S", all_rows, lam_base, global$alpha, seed_full, "direct")
lam <- lam_base * tuned$multiplier
cat("  alpha", global$alpha, "| local lambda", signif(lam, 3), "(x", tuned$multiplier, ")\n")

ctrl <- utils::modifyList(control_for("SCMR-EN"),
                          list(phi_warmup_sweeps = warm, base_model = base_model, verbose = FALSE))
ctrl <- do.call(scmr_control, ctrl[intersect(names(ctrl), names(formals(scmr_control)))])
t0 <- proc.time()[["elapsed"]]
pth <- scmr_split_merge_path(X_base[, global$x_colnames, drop = FALSE], dat$phase, seg_id, seg_xy,
                             G_max = G_max, alpha = global$alpha, lambda = lam, control = ctrl,
                             seed = seed_full, n_split = n_split, n_merge = n_merge,
                             keep_fits = FALSE, verbose = TRUE)
cat(sprintf("  path done in %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))
saveRDS(pth, file.path(path_dir, "path.rds"))
utils::write.csv(pth$table, file.path(path_dir, "path_table.csv"), row.names = FALSE)
utils::write.csv(pth$log, file.path(path_dir, "path_log.csv"), row.names = FALSE)
print(pth)
cat("PATH DONE\n")
