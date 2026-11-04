# design_case() copied unchanged from Naskah SCMR-EN/diagnostics/stage26_initialization_benchmark.R (8x8 grid, 6 classes).
design_case <- function(p, shape, size, signal, rep) {
  B <- 64L; C <- 6L; per <- 12L; n <- B * per
  cu <- as.matrix(expand.grid(East = 1:8, North = 1:8))
  truth <- if (shape == "rectangle") {
    which(cu[, 1] %in% 4:5 & cu[, 2] %in% if (size == 6L) 3:5 else 2:6)
  } else if (size == 6L) {
    which((cu[, 1] %in% 4:6 & cu[, 2] == 3) |
          (cu[, 1] == 4 & cu[, 2] == 4) |
          (cu[, 1] %in% 4:5 & cu[, 2] == 5))
  } else {
    which((cu[, 1] %in% 4:7 & cu[, 2] == 3) |
          (cu[, 1] == 4 & cu[, 2] %in% 4:6) |
          (cu[, 1] %in% 5:7 & cu[, 2] == 6))
  }
  stopifnot(length(truth) == size)
  block <- rep(seq_len(B), each = per)
  coords <- cu[block, , drop = FALSE]
  seed <- 8600000L + 100000L * match(p, c(24L, 105L)) +
    10000L * match(shape, c("rectangle", "L")) +
    1000L * match(size, c(6L, 10L)) +
    100L * match(signal, c("null", "intercept", "slope")) + rep
  set.seed(seed)
  shared <- rnorm(n)
  x <- matrix(sqrt(.6) * shared, n, p) +
    matrix(rnorm(n * p, sd = sqrt(.4)), n, p)
  x <- unclass(scale(x)); colnames(x) <- paste0("X", seq_len(p))
  base <- c(.25, .20, .17, .15, .13, .10)
  eta <- matrix(rep(log(base / base[C]), each = n), n, C)
  eta[, 1] <- eta[, 1] + .20 * x[, 1] - .15 * x[, 2]
  eta[, 2] <- eta[, 2] - .20 * x[, 1] + .15 * x[, 3]
  eta[, 3] <- eta[, 3] + .15 * x[, 4]
  eta[, C] <- 0
  inside <- block %in% truth
  if (signal == "intercept")
    eta[inside, ] <- eta[inside, ] - 2 *
      matrix(rep(log(base / base[C]), each = sum(inside)), sum(inside), C)
  if (signal == "slope") {
    v <- 1.3 * x[inside, 1]
    eta[inside, 1] <- eta[inside, 1] + v
    eta[inside, 2] <- eta[inside, 2] - v
  }
  e <- exp(eta - apply(eta, 1L, max)); prob <- e / rowSums(e)
  y <- factor(vapply(seq_len(n), function(i)
    sample.int(C, 1L, prob = prob[i, ]), integer(1L)),
    levels = seq_len(C), labels = paste0("C", seq_len(C)))
  control <- scmr_control(phi = .5, phi_update = "pl", phi_max = 2,
    k_neighbors = 4L, min_units = 5L, min_per_class = 2L,
    max_iter = 8L, n_starts = 1L, init_method = "kmeans",
    start_methods = "kmeans", init_tries = 100L,
    update_order = "fixed", lambda_scale = "sum",
    tune_lambda = FALSE, tune_alpha_after_convergence = FALSE,
    type_multinomial = "ungrouped", verbose = FALSE)
  list(x=x,y=y,block=block,coords=coords,truth=truth,control=control,
       seed=seed,p=p,shape=shape,size=size,signal=signal,rep=rep)
}

