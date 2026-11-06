# Feasibility of the phase-transition split rule.
# (1) G = 1 fit: how many spikes leave the Marchenko-Pastur bulk? Expect 0 under
#     the null and about G - 1 with G true regimes (if detectable).
# (2) Fit started at the true partition: each regime should show no spike.
# Variants: rows per block (SPEC_NB), within-block random effect sd (SPEC_TAU).
suppressMessages(library(scmr))
here <- Sys.getenv("SPEC_DIR", "inst/simulation/spectrum")
source(file.path(here, "spectrum_funs.R"))
env_num <- function(k, d) as.numeric(Sys.getenv(k, d))
side <- env_num("SPEC_SIDE", "20"); n_b <- env_num("SPEC_NB", "12"); p <- env_num("SPEC_P", "20")
tau <- env_num("SPEC_TAU", "0"); delta <- env_num("SPEC_DELTA", "1.5"); reps <- env_num("SPEC_REPS", "10")
cores <- env_num("SPEC_CORES", "2"); lam <- env_num("SPEC_LAMBDA", "0.005")
scens <- strsplit(Sys.getenv("SPEC_SCEN", "null,s3,s4,s4t"), ",")[[1]]
out <- Sys.getenv("SPEC_OUT", "spectrum.csv")
nboot <- env_num("SPEC_BOOT", "200")
K <- 4L; cls <- letters[1:K]

make_regimes <- function(shares, seed) {
  set.seed(seed); B <- side^2; G <- length(shares)
  if (G == 1L) return(rep(1L, B))
  target <- round(shares * B); target[1] <- B - sum(target[-1])
  gx <- rep(seq_len(side), side); gy <- rep(seq_len(side), each = side)
  nb <- lapply(seq_len(B), function(b) which(abs(gx - gx[b]) + abs(gy - gy[b]) == 1))
  lab <- rep(0L, B)
  for (g in G:2) {
    free <- which(lab == 0L); s <- free[sample.int(length(free), 1)]; member <- s
    while (length(member) < target[g]) {
      cand <- setdiff(unique(unlist(nb[member])), member); cand <- cand[lab[cand] == 0L]
      if (!length(cand)) break
      dd <- (gx[cand] - gx[s])^2 + (gy[cand] - gy[s])^2
      member <- c(member, cand[order(dd + stats::runif(length(dd), 0, 2))][1])
    }
    lab[member] <- g
  }
  lab[lab == 0L] <- 1L; lab
}
make_data <- function(scen, seed) {
  shares <- switch(scen, null = 1, s3 = c(.60, .30, .10), s4 = c(.55, .30, .10, .05), s4t = c(.60, .28, .10, .02))
  reg <- make_regimes(shares, seed); set.seed(seed + 1L)
  G <- length(shares); B <- side^2
  base <- matrix(0, K, p + 1); base[-1, 2:7] <- matrix(stats::rnorm((K - 1) * 6, 0, 0.8), K - 1)
  base[-1, 1] <- stats::rnorm(K - 1, 0, 0.3)
  betas <- lapply(seq_len(G), function(g) {
    if (g == 1L) return(base)
    b <- base; j <- sample(2:7, 3); k <- sample(2:K, 2)
    b[k, j] <- b[k, j] + delta * matrix(sample(c(-1, 1), length(k) * 3, TRUE), length(k))
    b[k, 1] <- b[k, 1] + stats::rnorm(length(k), 0, 0.5 * delta); b })
  x <- matrix(stats::rnorm(B * n_b * p), ncol = p, dimnames = list(NULL, paste0("x", 1:p)))
  r <- rep(reg, each = n_b)
  re <- matrix(stats::rnorm(B * K, 0, tau), B, K); re[, 1] <- 0     # within-block dependence
  eta <- matrix(0, nrow(x), K)
  for (g in seq_len(G)) { i <- which(r == g); eta[i, ] <- cbind(1, x[i, ]) %*% t(betas[[g]]) }
  eta <- eta + re[rep(seq_len(B), each = n_b), ]
  pr <- exp(eta - apply(eta, 1, max)); pr <- pr / rowSums(pr)
  u <- stats::runif(nrow(pr)); y <- factor(cls[1 + rowSums(u > t(apply(pr, 1, cumsum)))], levels = cls)
  gx <- rep(seq_len(side), side); gy <- rep(seq_len(side), each = side)
  list(x = x, y = y, id = rep(sprintf("B%04d", seq_len(B)), each = n_b),
       coords = cbind(rep(gx, each = n_b), rep(gy, each = n_b)), reg = reg, G = G)
}
ctrl <- scmr_control(lambda_scale = "sum", phi_update = "pl", min_units = 5, min_per_class = 2,
                     k_neighbors = 8, max_iter = 25, tiny_movement_max_units = 0, tiny_movement_rate_tol = 0,
                     tiny_movement_revert = FALSE, n_starts = 1)
adjacency <- function(units_xy) {
  D <- as.matrix(stats::dist(units_xy)); A <- (D > 0 & D <= 1.5) * 1; A
}
one <- function(job) {
  d <- make_data(job$scen, 1000L * job$rep + match(job$scen, c("null", "s3", "s4", "s4t")))
  units <- sort(unique(d$id)); unit_num <- match(d$id, units)
  truth <- d$reg[match(units, sprintf("B%04d", seq_len(side^2)))]
  f1 <- suppressWarnings(fit_scmr(d$x, d$y, model = "scmr", G = 1L, unit_id = d$id, coords = d$coords,
                                  alpha = 0.5, lambda = lam, control = ctrl, seed = job$rep))
  sp1 <- score_spectrum(unit_scores_info(d$x, d$y, f1$prob_train, unit_num))
  bt1 <- spectrum_boot(d$x, d$y, f1$prob_train, unit_num, B = nboot, seed = job$rep)
  uxy <- d$coords[match(units, d$id), , drop = FALSE]
  si1 <- unit_scores_info(d$x, d$y, f1$prob_train, unit_num)
  mo1 <- spike_moran(si1$U, si1$n_b, si1$J, adjacency(uxy), n_perm = 199L, seed = job$rep)
  rows <- data.frame(Scenario = job$scen, Rep = job$rep, TrueG = d$G, Fit = "G1", Regime = 1L,
    Units = length(units), Pure = d$G == 1L, d = sp1$d, gamma = sp1$gamma, sigma2 = sp1$sigma2,
    lambda1 = sp1$lambda[1], lambda2 = sp1$lambda[2], lambda3 = sp1$lambda[3], crit = sp1$crit, Spikes = sp1$spikes,
    R1 = bt1$R[1], Rcrit = bt1$crit, SpikesBoot = bt1$spikes,
    Moran1 = mo1$MoranI[1], MoranP1 = mo1$p[1], Moran2 = mo1$MoranI[2], MoranP2 = mo1$p[2])
  if (d$G > 1L) {
    ft <- suppressWarnings(fit_scmr(d$x, d$y, model = "scmr", G = d$G, unit_id = d$id, coords = d$coords,
                                    alpha = 0.5, lambda = lam, control = ctrl, seed = job$rep, initial_cluster = truth))
    ug <- as.integer(ft$group_unit[units])
    for (g in seq_len(d$G)) {
      ug_rows <- which(ug[unit_num] == g)
      if (length(unique(unit_num[ug_rows])) < 5L) next
      sp <- score_spectrum(unit_scores_info(d$x, d$y, ft$prob_train, unit_num, ug_rows))
      bt <- spectrum_boot(d$x, d$y, ft$prob_train, unit_num, ug_rows, B = nboot, seed = job$rep + g)
      sig <- unit_scores_info(d$x, d$y, ft$prob_train, unit_num, ug_rows)
      ug_units <- sort(unique(unit_num[ug_rows]))
      mo <- spike_moran(sig$U, sig$n_b, sig$J, adjacency(uxy[ug_units, , drop = FALSE]), n_perm = 199L, seed = job$rep + g)
      pure <- mean(truth[ug == g] == as.integer(names(which.max(table(truth[ug == g])))))
      rows <- rbind(rows, data.frame(Scenario = job$scen, Rep = job$rep, TrueG = d$G, Fit = "truthstart", Regime = g,
        Units = sp$m, Pure = pure, d = sp$d, gamma = sp$gamma, sigma2 = sp$sigma2, lambda1 = sp$lambda[1],
        lambda2 = sp$lambda[2], lambda3 = sp$lambda[3], crit = sp$crit, Spikes = sp$spikes,
        R1 = bt$R[1], Rcrit = bt$crit, SpikesBoot = bt$spikes,
        Moran1 = mo$MoranI[1], MoranP1 = mo$p[1], Moran2 = mo$MoranI[2], MoranP2 = mo$p[2]))
    }
  }
  rows$NB <- n_b; rows$Tau <- tau; rows
}
jobs <- unlist(lapply(seq_len(reps), function(r) lapply(scens, function(s) list(scen = s, rep = r))), recursive = FALSE)
res <- do.call(rbind, parallel::mclapply(jobs, function(j) tryCatch(one(j), error = function(e) {
  message("fail ", j$scen, " ", j$rep, ": ", conditionMessage(e)); NULL }), mc.cores = cores))
utils::write.csv(res, out, row.names = FALSE)
g1 <- res[res$Fit == "G1", ]
cat(sprintf("\nn_b = %g, tau = %g: G = 1 fit, spikes above the TW edge (target TrueG - 1)\n", n_b, tau))
print(table(Scenario = g1$Scenario, Spikes = g1$Spikes))
cat("Bootstrap-calibrated:\n"); print(table(Scenario = g1$Scenario, SpikesBoot = g1$SpikesBoot))
cat("Leading eigenvector spatially ordered (Moran p < 0.05):\n"); print(table(Scenario = g1$Scenario, Ordered = g1$MoranP1 < 0.05))
print(aggregate(cbind(sigma2, lambda1, crit, gamma) ~ Scenario, g1, mean), digits = 3)
ts <- res[res$Fit == "truthstart", ]
if (nrow(ts)) { cat("\nTruth-start regimes: spikes (target 0)\n"); print(table(Scenario = ts$Scenario, Spikes = ts$Spikes))
  cat("Bootstrap-calibrated:\n"); print(table(Scenario = ts$Scenario, SpikesBoot = ts$SpikesBoot))
  cat("Spike AND spatially ordered:\n"); print(table(Scenario = ts$Scenario, Split = ts$SpikesBoot > 0 & ts$MoranP1 < 0.05)) }
