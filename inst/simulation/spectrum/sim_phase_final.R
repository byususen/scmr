# Final simulation of the two-stage phase test (scmr_phase_test).
# (1) G = 1 fit: the single regime should be "discrete" when true regimes exist
#     (s3, s4, s4t) and not discrete under the null or smooth variation
#     (grad: linear gradient; wave: sinusoidal surface).
# (2) Fit started at the true partition: no regime should be "discrete".
# Variants: rows per block (SPEC_NB), within-block random effect sd (SPEC_TAU),
# predictors (SPEC_P). Other env: SPEC_REPS, SPEC_SCEN, SPEC_BOOT, SPEC_PERM,
# SPEC_CORES, SPEC_OUT.
suppressMessages(library(scmr))
here <- Sys.getenv("SPEC_DIR", "inst/simulation/spectrum")
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
  smooth <- scen %in% c("grad", "wave")
  shares <- switch(scen, null = 1, grad = 1, wave = 1, s3 = c(.60, .30, .10), s4 = c(.55, .30, .10, .05),
                   s4t = c(.60, .28, .10, .02))
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
  if (smooth) {
    # Continuous spatial variation, no discrete regimes: beta(s) = base + f(s) D.
    gx0 <- rep(seq_len(side), side); gy0 <- rep(seq_len(side), each = side)
    f <- if (scen == "grad") (gx0 + gy0 - side - 1) / (side - 1) else
      sin(pi * gx0 / side * 2) * cos(pi * gy0 / side * 2)
    D <- matrix(0, K, p + 1); kk <- sample(2:K, 2); jj <- sample(2:7, 3)
    D[kk, jj] <- delta * matrix(sample(c(-1, 1), 6, TRUE), 2); D[kk, 1] <- delta * sample(c(-1, 1), 2, TRUE)
    fr <- f[rep(seq_len(B), each = n_b)]
    eta <- eta + fr * (cbind(1, x) %*% t(D))
  }
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
nperm <- env_num("SPEC_PERM", "9999")
one <- function(job) {
  d <- make_data(job$scen, 1000L * job$rep + match(job$scen, c("null", "s3", "s4", "s4t", "grad", "wave")))
  units <- sort(unique(d$id))
  truth <- d$reg[match(units, sprintf("B%04d", seq_len(side^2)))]
  f1 <- suppressWarnings(fit_scmr(d$x, d$y, model = "scmr", G = 1L, unit_id = d$id, coords = d$coords,
                                  alpha = 0.5, lambda = lam, control = ctrl, seed = job$rep))
  t1 <- scmr_phase_test(f1, d$x, d$y, d$id, d$coords, B = nboot, n_perm = nperm, seed = job$rep)
  out <- data.frame(Scenario = job$scen, Rep = job$rep, TrueG = d$G, Fit = "G1", t1)
  if (d$G > 1L) {
    ft <- suppressWarnings(fit_scmr(d$x, d$y, model = "scmr", G = d$G, unit_id = d$id, coords = d$coords,
                                    alpha = 0.5, lambda = lam, control = ctrl, seed = job$rep, initial_cluster = truth))
    tt <- scmr_phase_test(ft, d$x, d$y, d$id, d$coords, B = nboot, n_perm = nperm, seed = job$rep)
    out <- rbind(out, data.frame(Scenario = job$scen, Rep = job$rep, TrueG = d$G, Fit = "truthstart", tt))
  }
  out$NB <- n_b; out$Tau <- tau; out$P <- p
  out
}
jobs <- unlist(lapply(seq_len(reps), function(r) lapply(scens, function(s) list(scen = s, rep = r))), recursive = FALSE)
res <- do.call(rbind, parallel::mclapply(jobs, function(j) tryCatch(one(j), error = function(e) {
  message("fail ", j$scen, " ", j$rep, ": ", conditionMessage(e)); NULL }), mc.cores = cores))
utils::write.csv(res, out, row.names = FALSE)
g1 <- res[res$Fit == "G1", ]
cat(sprintf("\nn_b = %g, tau = %g, p = %g: G = 1 fit, phase of the single regime\n", n_b, tau, p))
print(table(Scenario = g1$Scenario, Phase = g1$Phase))
ts <- res[res$Fit == "truthstart", ]
if (nrow(ts)) { cat("Regimes of the true partition (target: none discrete)\n"); print(table(Scenario = ts$Scenario, Phase = ts$Phase)) }
