# Scan versus split-merge path on the Stage 26/28 8x8 design (G = 2 fixed).
# Methods (same SCMR engine, penalty and Potts settings for every final fit):
#   kmeans5      five k-means coordinate starts, best penalized objective
#   scan_current Stage 26 flexible scan (BFS + rectangles + beam + motifs, sizes
#                5..14, local-likelihood intercept ranking, within-size ranks), 5 fits
#   path_kmeans  scmr_split_merge_path, 2-means score split (1 fit at G = 2)
#   path_scan    scmr_split_merge_path, scan split (sizes 5..32, 5 fits at G = 2)
#   truth        one fit started at the simulated partition (diagnostic only)
# Env: CMP_REP_FROM, CMP_REP_TO, CMP_CASES (max cells, 18 incl. 2 null),
#      CMP_METHODS (comma list), CMP_OUT, CMP_CORES.
suppressMessages(library(scmr))
here <- Sys.getenv("CMP_DIR", "inst/simulation/scan_compare")
ns <- asNamespace("scmr")
for (nm in ls(ns, all.names = TRUE)) if (!exists(nm, envir = globalenv(), inherits = FALSE))
  assign(nm, get(nm, envir = ns), envir = globalenv())
source(file.path(here, "flexible_scan.R"))
source(file.path(here, "local_scan.R"))
prepare_scmr_data <- function(x, y, unit_id, coords, component = NULL) {
  if (!is.null(component)) stop("Only one spatial component is supported")
  get("prepare_scmr_data", envir = ns)(x, y, unit_id, coords)
}
build_unit_weights <- function(coords, k = 5L, component = NULL, ...) {
  get("build_unit_weights", envir = ns)(coords, k = k, ...)
}
source(file.path(here, "design_stage26.R"))
env_int <- function(k, d) as.integer(Sys.getenv(k, d))
rep_from <- env_int("CMP_REP_FROM", "1"); rep_to <- env_int("CMP_REP_TO", "1")
max_cases <- env_int("CMP_CASES", "18"); cores <- env_int("CMP_CORES", "1")
methods <- strsplit(Sys.getenv("CMP_METHODS", "kmeans5,scan_current,path_kmeans,path_scan,truth"), ",")[[1]]
out <- Sys.getenv("CMP_OUT", "compare_grid8.csv")

cases <- expand.grid(P = c(24L, 105L), Shape = c("rectangle", "L"), Size = c(6L, 10L),
                     Signal = c("intercept", "slope"), stringsAsFactors = FALSE)
cases <- rbind(cases, data.frame(P = c(24L, 105L), Shape = "L", Size = 6L, Signal = "null"))
cases <- cases[seq_len(max_cases), , drop = FALSE]

fit_g2 <- function(d, init = NULL, seed, method = "kmeans") {
  ctl <- d$control; ctl$n_starts <- 1L; ctl$init_method <- method; ctl$start_methods <- method
  tryCatch(suppressWarnings(fit_scmr(d$x, d$y, model = "scmr", G = 2L, unit_id = d$block,
    coords = d$coords, penalty = "elastic_net", alpha = .5, lambda = .02, control = ctl,
    seed = seed, initial_cluster = init)), error = function(e) NULL)
}
objective <- function(f) if (is.null(f)) -Inf else f$criteria$PenalizedObjective
minor <- function(f) {
  lab <- as.integer(f$group_unit[as.character(seq_len(64L))])
  if (anyNA(lab) || length(unique(lab)) != 2L) return(integer())
  which(lab == which.min(tabulate(lab, 2L)))
}
n_components <- function(region) {
  if (!length(region)) return(0L)
  rem <- region; n <- 0L
  while (length(rem)) {
    n <- n + 1L; seen <- rem[1L]; fr <- seen
    while (length(fr)) {
      u <- fr[1L]; fr <- fr[-1L]
      near <- rem[abs((rem - 1L) %% 8L - (u - 1L) %% 8L) + abs((rem - 1L) %/% 8L - (u - 1L) %/% 8L) == 1L]
      add <- setdiff(near, seen); seen <- union(seen, add); fr <- union(fr, add)
    }
    rem <- setdiff(rem, seen)
  }
  n
}
truth_labels <- function(d) { lab <- rep(1L, 64L); lab[d$truth] <- 2L; lab[order(as.character(seq_len(64L)))] }

run_method <- function(d, method) {
  t0 <- proc.time()[["elapsed"]]
  fits <- list(); n_fits <- 0L
  if (method == "kmeans5") {
    fits <- lapply(1:5, function(j) fit_g2(d, seed = d$seed + 1000L * j)); n_fits <- 5L
  } else if (method == "scan_current") {
    one <- suppressWarnings(fit_scmr(d$x, d$y, model = "scmr", G = 1L, unit_id = d$block, coords = d$coords,
      penalty = "elastic_net", alpha = .5, lambda = .02, control = d$control, seed = 11L))
    pool <- scmr_flexible_scan_starts(d$x, d$y, d$block, d$coords, one, k_neighbors = 4L, sizes = 5:14,
      top_k = 5L, min_per_class = 2L, beam_width = 4L, max_jaccard = .6, include_motifs = TRUE,
      motif_neighbors = 9L, motif_sizes = 5:8)
    pool <- scmr_add_local_likelihood_scores(pool, d$y, d$block, one, ridge = .1, max_iter = 15L)
    st <- suppressWarnings(scmr_select_flexible_starts(pool, top_k = 5L, max_jaccard = .6, slope_score = "raw",
      intercept_score = "local_likelihood", size_calibration = "within_size")$starts)
    fits <- lapply(seq_along(st), function(j) fit_g2(d, st[[j]], d$seed + 1000L * j)); n_fits <- length(st)
  } else if (method %in% c("path_kmeans", "path_scan")) {
    p <- tryCatch(suppressWarnings(scmr_split_merge_path(d$x, d$y, d$block, d$coords, G_max = 2L,
      penalty = "elastic_net", alpha = .5, lambda = .02, control = d$control, seed = d$seed + 1000L,
      n_split = 1L, merge_back = FALSE, split_rule = if (method == "path_scan") "scan" else "kmeans",
      scan_n = 5L, scan_max_size = env_int("CMP_SCAN_MAX", "32"), scan_anchors = 64L, scan_beam = 4L)), error = function(e) NULL)
    if (!is.null(p) && length(p$fits) >= 2L) fits <- list(p$fits[[2L]])
    n_fits <- if (is.null(p)) 0L else sum(p$log$G == 2L)
  } else if (method == "truth") {
    fits <- list(fit_g2(d, truth_labels(d), d$seed + 1000L)); n_fits <- 1L
  }
  q <- vapply(fits, objective, numeric(1))
  win <- if (length(q) && any(is.finite(q))) fits[[which.max(q)]] else NULL
  pred <- if (is.null(win)) integer() else minor(win)
  tp <- length(intersect(pred, d$truth))
  data.frame(P = d$p, Shape = d$shape, Size = d$size, Signal = d$signal, Rep = d$rep, Seed = d$seed,
    Method = method, FitsG2 = n_fits, Objective = if (is.null(win)) NA_real_ else max(q),
    IoU = if (d$signal == "null" || !length(pred)) NA_real_ else tp / length(union(pred, d$truth)),
    Exact = if (d$signal == "null" || !length(pred)) NA_integer_ else as.integer(setequal(pred, d$truth)),
    Components = if (length(pred)) n_components(pred) else NA_integer_,
    Seconds = proc.time()[["elapsed"]] - t0)
}
jobs <- expand.grid(i = seq_len(nrow(cases)), rep = rep_from:rep_to)
res <- parallel::mclapply(seq_len(nrow(jobs)), function(k) {
  cc <- cases[jobs$i[k], ]
  d <- design_case(cc$P, cc$Shape, cc$Size, cc$Signal, jobs$rep[k])
  do.call(rbind, lapply(methods, function(m) tryCatch(run_method(d, m), error = function(e) {
    message("FAILED ", m, " case ", jobs$i[k], " rep ", jobs$rep[k], ": ", conditionMessage(e)); NULL })))
}, mc.cores = cores)
res <- do.call(rbind, res)
if (nzchar(Sys.getenv("CMP_LABEL"))) res$Method <- paste0(res$Method, Sys.getenv("CMP_LABEL"))
utils::write.csv(res, out, row.names = FALSE)
alt <- res[res$Signal != "null", ]
print(aggregate(cbind(IoU, Exact, Components, Seconds) ~ Method, alt, mean), digits = 3)
print(aggregate(IoU ~ Method + Signal + Shape, alt, mean), digits = 3)
