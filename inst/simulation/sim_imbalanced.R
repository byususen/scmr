# =============================================================================
# KSA-like imbalanced regimes: does the criterion find G, and does the
# split-merge path beat the multi-start fit?
#
# Blocks on a side x side grid, n_b rows per block, p predictors (first 6
# active), K = 4 classes. Regimes are contiguous; shares are imbalanced.
# Scenarios (env SIM_SCEN, comma list):
#   null : G = 1
#   s3   : shares 0.60 / 0.30 / 0.10
#   s4   : shares 0.55 / 0.30 / 0.10 / 0.05
#   s4t  : shares 0.60 / 0.28 / 0.10 / 0.02 (tiny regime)
# Regime deviations from the base coefficients (env SIM_DEV):
#   dense  : 3 active slopes x 2 classes (+/- delta) and intercept shifts
#   sparse : intercept shifts (+/- delta) for 2 classes and one slope for 1 class
# Other env: SIM_SIDE, SIM_NB, SIM_P, SIM_REPS, SIM_REP0 (first replicate),
# SIM_GMAX, SIM_DELTA, SIM_WARM (phi warm-up sweeps), SIM_BASE (none/global),
# SIM_CORES, SIM_OUT, SIM_STANDARD (1 = also run the multi-start fits),
# SIM_SPLIT (path split rule: kmeans, scan or both).
# Methods, each over G = 1..G_max:
#   standard : fit_scmr with n_starts = 3 (kmeans, coefficient, kmeans)
#   path     : scmr_split_merge_path (score splits + backward merges)
# Output: one row per replicate x scenario x method x G (criteria, objective,
# ARI vs truth, seconds) in SIM_OUT (default sim_imbalanced.csv).
# =============================================================================
suppressMessages(library(scmr))
env_int <- function(k, d) as.integer(Sys.getenv(k, d))
side <- env_int("SIM_SIDE", "20"); n_b <- env_int("SIM_NB", "12"); p <- env_int("SIM_P", "20")
reps <- env_int("SIM_REPS", "6"); rep0 <- env_int("SIM_REP0", "1"); G_max <- env_int("SIM_GMAX", "5"); delta <- as.numeric(Sys.getenv("SIM_DELTA", "1"))
warm <- env_int("SIM_WARM", "0"); cores <- env_int("SIM_CORES", "2")
base_model <- Sys.getenv("SIM_BASE", "none")
dev_type <- Sys.getenv("SIM_DEV", "dense")
run_standard <- env_int("SIM_STANDARD", "1") == 1L
split_rule <- Sys.getenv("SIM_SPLIT", "kmeans")
scens <- strsplit(Sys.getenv("SIM_SCEN", "null,s3,s4"), ",")[[1]]
out_file <- Sys.getenv("SIM_OUT", "sim_imbalanced.csv")
K <- 4L; cls <- letters[1:K]

ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); c2 <- function(z) sum(z * (z - 1) / 2)
  s <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); e <- sa * sb / c2(n)
  if (sa + sb == 2 * e) return(1)
  (s - e) / ((sa + sb) / 2 - e)
}
# Contiguous regimes grown from seeds by a random BFS with target sizes.
make_regimes <- function(shares, seed) {
  set.seed(seed)
  B <- side^2; G <- length(shares)
  if (G == 1L) return(rep(1L, B))
  target <- round(shares * B); target[1] <- B - sum(target[-1])
  gx <- rep(seq_len(side), side); gy <- rep(seq_len(side), each = side)
  nb <- lapply(seq_len(B), function(b) which(abs(gx - gx[b]) + abs(gy - gy[b]) == 1))
  lab <- rep(0L, B)
  for (g in G:2) {                       # small regimes first, as compact blobs
    free <- which(lab == 0L)
    s <- free[sample.int(length(free), 1)]
    member <- s; frontier <- s
    while (length(member) < target[g]) {
      cand <- setdiff(unique(unlist(nb[member])), member)
      cand <- cand[lab[cand] == 0L]
      if (!length(cand)) break
      d <- (gx[cand] - gx[s])^2 + (gy[cand] - gy[s])^2
      member <- c(member, cand[order(d + stats::runif(length(d), 0, 2))][1])
    }
    lab[member] <- g
  }
  lab[lab == 0L] <- 1L
  lab
}
make_data <- function(scen, seed) {
  shares <- switch(scen, null = 1, s3 = c(.60, .30, .10), s4 = c(.55, .30, .10, .05),
                   s4t = c(.60, .28, .10, .02))
  reg <- make_regimes(shares, seed)
  set.seed(seed + 1L)
  G <- length(shares); B <- side^2
  base <- matrix(0, K, p + 1); base[-1, 2:7] <- matrix(stats::rnorm((K - 1) * 6, 0, 0.8), K - 1)
  base[-1, 1] <- stats::rnorm(K - 1, 0, 0.3)
  betas <- lapply(seq_len(G), function(g) {
    if (g == 1L) return(base)
    b <- base; k <- sample(2:K, 2)
    if (dev_type == "sparse") {
      b[k, 1] <- b[k, 1] + delta * sample(c(-1, 1), 2, TRUE)
      j <- sample(2:7, 1)
      b[k[1], j] <- b[k[1], j] + delta * sample(c(-1, 1), 1)
    } else {
      j <- sample(2:7, 3)
      b[k, j] <- b[k, j] + delta * matrix(sample(c(-1, 1), length(k) * 3, TRUE), length(k))
      b[k, 1] <- b[k, 1] + stats::rnorm(length(k), 0, 0.5 * delta)
    }
    b
  })
  x <- matrix(stats::rnorm(B * n_b * p), ncol = p, dimnames = list(NULL, paste0("x", 1:p)))
  r <- rep(reg, each = n_b)
  eta <- matrix(0, nrow(x), K)
  for (g in seq_len(G)) { i <- which(r == g); eta[i, ] <- cbind(1, x[i, ]) %*% t(betas[[g]]) }
  pr <- exp(eta - apply(eta, 1, max)); pr <- pr / rowSums(pr)
  u <- stats::runif(nrow(pr)); y <- factor(cls[1 + rowSums(u > t(apply(pr, 1, cumsum)))], levels = cls)
  gx <- rep(seq_len(side), side); gy <- rep(seq_len(side), each = side)
  list(x = x, y = y, id = rep(sprintf("B%04d", seq_len(B)), each = n_b),
       coords = cbind(rep(gx, each = n_b), rep(gy, each = n_b)), reg = reg, G = G)
}
ctrl_of <- function(n_starts) scmr_control(lambda_scale = "sum", phi_update = "pl", min_units = 5,
  min_per_class = 2, k_neighbors = 8, max_iter = 25, tiny_movement_max_units = 0, tiny_movement_rate_tol = 0,
  tiny_movement_revert = FALSE, n_starts = n_starts, start_methods = c("kmeans", "coefficient", "kmeans"),
  coef_init_k = 30, coef_init_anchors = 200, phi_warmup_sweeps = warm, alpha_grid = 0.5, base_model = base_model)

one <- function(job) {
  scen <- job$scen; rep <- job$rep
  d <- make_data(scen, 1000L * rep + match(scen, c("null", "s3", "s4", "s4t")))
  units <- sort(unique(d$id)); truth <- d$reg[match(units, sprintf("B%04d", seq_len(side^2)))]
  gl <- fit_scmr(d$x, d$y, model = "global", penalty = "elastic_net", alpha = 0.5, nfolds = 5, seed = rep)
  lam <- gl$lambda_used[1]
  rows <- list()
  add <- function(method, G, f, sec) {
    cr <- f$criteria; st <- f$diagnostics$starts
    obj <- if (is.null(st) || !nrow(st)) cr$PenalizedObjective else st$Objective[which(st$Selected)[1]]
    rows[[length(rows) + 1L]] <<- data.frame(Scenario = scen, Dev = dev_type, Warm = warm, Base = base_model, Split = split_rule,
      B = side^2, Rep = rep, TrueG = d$G, Method = method, G = G,
      Objective = obj, PLIC_AIC = cr$CriterionPLIC_AIC, PLIC_BIC = cr$CriterionPLIC_BIC,
      PLIC_BIC_rows = cr$CriterionPLIC_BIC_rows, PLIC_BIC_units = cr$CriterionPLIC_BIC_units,
      DF = cr$df_effective_total, Phi = f$phi %||% NA_real_, ARI = ari(truth, f$group_unit[units]),
      MinUnits = min(table(f$group_unit)), Seconds = sec, Lambda = lam)
  }
  if (run_standard) for (G in seq_len(G_max)) {
    t0 <- proc.time()[["elapsed"]]
    f <- tryCatch(fit_scmr(d$x, d$y, model = "scmr", G = G, unit_id = d$id, coords = d$coords, alpha = 0.5,
                           lambda = lam, control = ctrl_of(3L), seed = rep), error = function(e) NULL)
    if (!is.null(f)) add("standard", G, f, proc.time()[["elapsed"]] - t0)
  }
  t0 <- proc.time()[["elapsed"]]
  pth <- scmr_split_merge_path(d$x, d$y, d$id, d$coords, G_max = G_max, alpha = 0.5, lambda = lam,
                               control = ctrl_of(1L), seed = rep, n_split = 2L, n_merge = 2L,
                               split_rule = split_rule, scan_n = 2L)
  tot <- proc.time()[["elapsed"]] - t0
  for (G in seq_along(pth$fits)) add("path", G, pth$fits[[G]], tot / length(pth$fits))
  do.call(rbind, rows)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
jobs <- unlist(lapply(rep0 - 1L + seq_len(reps), function(r) lapply(scens, function(s) list(scen = s, rep = r))), recursive = FALSE)
res <- parallel::mclapply(jobs, function(j) tryCatch(one(j), error = function(e) {
  message("job ", j$scen, " rep ", j$rep, " failed: ", conditionMessage(e)); NULL }), mc.cores = cores)
res <- do.call(rbind, res)
utils::write.csv(res, out_file, row.names = FALSE)

# Summary: chosen G per criterion, and ARI at the chosen G.
crit <- c("PLIC_AIC", "PLIC_BIC", "PLIC_BIC_rows", "PLIC_BIC_units")
sel <- do.call(rbind, lapply(split(res, list(res$Scenario, res$Rep, res$Method), drop = TRUE), function(z) {
  do.call(rbind, lapply(crit, function(k) { i <- which.min(z[[k]]); data.frame(Scenario = z$Scenario[1],
    Rep = z$Rep[1], Method = z$Method[1], Criterion = k, TrueG = z$TrueG[1], ChosenG = z$G[i], ARI = z$ARI[i]) }))
}))
agg <- aggregate(cbind(Correct = ChosenG == TrueG, Over = ChosenG > TrueG, ARI) ~ Scenario + Method + Criterion,
                 data = sel, FUN = mean)
print(agg[order(agg$Scenario, agg$Criterion, agg$Method), ], row.names = FALSE, digits = 3)
tg <- res[res$G == res$TrueG, ]
if (length(unique(tg$Method)) == 2L) {
  cmp <- reshape(tg[, c("Scenario", "Rep", "Method", "Objective", "ARI")], idvar = c("Scenario", "Rep"),
                 timevar = "Method", direction = "wide")
  cat("\nAt the true G: path minus standard objective and ARI\n")
  print(aggregate(cbind(dObj = Objective.path - Objective.standard, ARI.path, ARI.standard) ~ Scenario, cmp, mean), digits = 3)
}
cat("\nSeconds per replicate (all G):\n")
print(aggregate(Seconds ~ Scenario + Method, res, sum), digits = 3)
