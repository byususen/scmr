# Panel Markov generator for dynamic SCMR studies (unit x time).
#
#   P(y_it = c | y_i,t-1 = k, x_it, g_i = g) proportional to
#     exp(gamma_g[k, c] + x_it' beta_g[, c]),
# with a crop-calendar-like cycle 1 -> 2 -> ... -> (C-1) -> 1 among the first
# C - 1 classes and an "other" class C. Regimes differ in calendar speed and in
# the response to the features; delta scales the feature differences.

panel_regimes <- function(ucoords, pattern, G, seed) {
  m <- nrow(ucoords)
  if (pattern %in% c("global", "smooth")) return(rep(1L, m))
  if (pattern == "blocks") {
    km <- with_scmr_seed(seed, stats::kmeans(scale(ucoords), centers = G, nstart = 20))
    return(as.integer(factor(km$cluster, levels = order(km$centers[, 1]))))
  }
  # "irregular": G x G Latin-square tiling of the bounding box; every regime
  # consists of G disconnected tiles.
  r1 <- (ucoords[, 1] - min(ucoords[, 1])) / max(diff(range(ucoords[, 1])), 1e-12)
  r2 <- (ucoords[, 2] - min(ucoords[, 2])) / max(diff(range(ucoords[, 2])), 1e-12)
  col <- pmin(floor(r1 * G), G - 1)
  row <- pmin(floor(r2 * G), G - 1)
  as.integer((row + col) %% G) + 1L
}

panel_transition_logits <- function(C, stay, advance) {
  cyc <- C - 1L
  gam <- matrix(-4, C, C)
  for (k in seq_len(cyc)) {
    gam[k, k] <- stay
    gam[k, (k %% cyc) + 1L] <- advance
    gam[k, C] <- -2.5
  }
  gam[C, C] <- 2
  gam[C, 1L] <- 0
  gam
}

#' Simulate panel data from a spatially clustered multinomial Markov model
#'
#' @param n_units Number of units when `coords` is not supplied (uniform on the
#'   unit square).
#' @param coords Optional unit coordinates (one row per unit), for example the
#'   real survey locations.
#' @param n_time Number of time points per unit.
#' @param pattern Regime pattern: `"global"`, `"blocks"` (contiguous k-means
#'   regions), `"irregular"` (Latin-square tiles, disconnected regimes) or
#'   `"smooth"` (coefficients varying linearly in space, no regimes).
#' @param G Number of regimes for `"blocks"` and `"irregular"`.
#' @param n_classes Number of classes C (C - 1 cyclic phases plus "other").
#' @param p_active,p_inactive Numbers of active and inactive features.
#' @param rho AR(1) correlation between features.
#' @param delta Scale of the regime differences in the feature effects.
#' @param stay_range Range of the regime-specific stay logits (calendar speed).
#' @param seed Random seed; the caller's random state is restored.
#' @return A list with `x` (features), `lag_x` (one-hot previous class), `y`,
#'   `lag`, `unit_id`, `time`, `coords` (per row), `regime` (per row, `NA` for
#'   `"smooth"`), `unit_coords`, `unit_regime`, `beta` and `gamma`.
#' @export
simulate_scmr_panel <- function(n_units = 300, coords = NULL, n_time = 12,
                                pattern = c("irregular", "blocks", "global", "smooth"),
                                G = 3, n_classes = 6, p_active = 5, p_inactive = 10,
                                rho = 0.5, delta = 1, stay_range = c(0.2, 1.4), seed = 1) {
  pattern <- match.arg(pattern)
  if (pattern == "global") G <- 1L
  C <- as.integer(n_classes)
  if (C < 3L) stop("n_classes must be at least 3.", call. = FALSE)
  p <- p_active + p_inactive
  with_scmr_seed(seed, {
    ucoords <- if (is.null(coords)) cbind(stats::runif(n_units), stats::runif(n_units)) else as.matrix(coords)
    m <- nrow(ucoords)
    reg <- panel_regimes(ucoords, pattern, G, seed + 1L)
    Gr <- max(reg)
    stay <- if (Gr == 1L) mean(stay_range) else seq(stay_range[1], stay_range[2], length.out = Gr)
    gamma <- lapply(seq_len(Gr), function(g) panel_transition_logits(C, stay[g], 1))
    # Feature effects: common part + regime-specific part (centred over classes).
    centre <- function(b) sweep(b, 1, rowMeans(b), "-")
    base <- centre(matrix(stats::rnorm(p_active * C, 0, 0.6), p_active, C))
    dev <- lapply(seq_len(max(Gr, 2L)), function(g) centre(matrix(stats::rnorm(p_active * C), p_active, C)))
    beta_active <- lapply(seq_len(Gr), function(g) base + if (Gr > 1L) delta * dev[[g]] else 0)
    S <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
    R <- chol(S)
    n <- m * n_time
    x <- matrix(stats::rnorm(n * p), n, p) %*% R
    colnames(x) <- paste0("x", seq_len(p))
    unit_row <- rep(seq_len(m), each = n_time)
    time <- rep(seq_len(n_time), m)
    beta_row <- function(i) {
      if (pattern != "smooth") return(beta_active[[reg[i]]])
      s <- (ucoords[i, ] - colMeans(ucoords)) / apply(ucoords, 2, stats::sd)
      base + delta * (s[1] * dev[[1]] + s[2] * dev[[2]]) / sqrt(2)
    }
    y <- integer(n)
    lag <- integer(n)
    init <- c(rep(1, C - 1L), 0.3)
    for (i in seq_len(m)) {
      b <- beta_row(i)
      gam <- gamma[[reg[i]]]
      prev <- sample.int(C, 1L, prob = init)
      for (t in seq_len(n_time)) {
        r <- (i - 1L) * n_time + t
        eta <- gam[prev, ] + as.numeric(x[r, seq_len(p_active)] %*% b)
        pr <- exp(eta - max(eta))
        lag[r] <- prev
        y[r] <- sample.int(C, 1L, prob = pr / sum(pr))
        prev <- y[r]
      }
    }
  })
  classes <- as.character(seq_len(C))
  lag_x <- matrix(0, n, C, dimnames = list(NULL, paste0("lag_X", classes)))
  lag_x[cbind(seq_len(n), lag)] <- 1
  ids <- sprintf("U%04d", seq_len(m))
  list(x = x, lag_x = lag_x, y = factor(classes[y], levels = classes),
       lag = factor(classes[lag], levels = classes), unit_id = ids[unit_row], time = time,
       coords = ucoords[unit_row, , drop = FALSE],
       regime = if (pattern == "smooth") rep(NA_integer_, n) else reg[unit_row],
       unit_coords = ucoords, unit_regime = if (pattern == "smooth") rep(NA_integer_, m) else reg,
       beta = beta_active, gamma = gamma, pattern = pattern)
}
