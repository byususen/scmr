# Scan-type split proposals for the split-merge path.
#
# Given a fitted regime g, look for a connected set W of its units whose rows
# deviate from the regime's fit. For unit b and class c < C, with standardized
# predictors z_i = (1, x_i) and the fitted probabilities p_ic of the current fit,
#   U_b[c, j] = sum_{i in b} z_ij (y_ic - p_ic),
#   V_b[c, j] = sum_{i in b} z_ij^2 p_ic (1 - p_ic),
# and for a region W the diagonal score statistics
#   intercept(W) = sum_c U_W[c, 0]^2 / (V_W[c, 0] + ridge),
#   slope(W)     = max_j sum_c U_W[c, j]^2 / (V_W[c, j] + ridge),
# with U_W = sum_{b in W} U_b (Jung, Kulldorff & Richard, 2010; Lee, Gangnon &
# Zhu, 2017, for the intercept / slope distinction). Connected regions are grown
# by a beam search over the spatial graph restricted to the regime (cf. Tango &
# Takahashi, 2005, for flexible connected scan windows). The statistics only rank
# split proposals; every proposal is refitted by the full SCMR alternation.

# Unit-level score sums U and information sums V (units x (C - 1) * (p + 1)).
scan_unit_scores <- function(Z, Y, P, unit_num, n_units) {
  C <- ncol(Y)
  d <- ncol(Z)
  U <- matrix(0, n_units, (C - 1L) * d)
  V <- matrix(0, n_units, (C - 1L) * d)
  for (c in seq_len(C - 1L)) {
    jj <- (c - 1L) * d + seq_len(d)
    U[, jj] <- rowsum(Z * (Y[, c] - P[, c]), unit_num, reorder = TRUE)
    V[, jj] <- rowsum(Z^2 * (P[, c] * (1 - P[, c])), unit_num, reorder = TRUE)
  }
  list(U = U, V = V, d = d, K = C - 1L)
}

# Channel statistics for many regions at once: rows of SU, SV are region sums.
scan_channel_stats <- function(SU, SV, d, K, ridge) {
  fs <- SU^2 / (SV + ridge)
  if (is.null(dim(fs))) fs <- matrix(fs, nrow = 1L)
  icpt <- (seq_len(K) - 1L) * d + 1L
  intercept <- rowSums(fs[, icpt, drop = FALSE])
  slope <- if (d > 1L) {
    per_pred <- 0
    for (c in seq_len(K)) per_pred <- per_pred + fs[, (c - 1L) * d + 2:d, drop = FALSE]
    apply(per_pred, 1L, max)
  } else rep(0, nrow(fs))
  cbind(intercept = intercept, slope = slope)
}

# Beam search for connected high-score regions inside one regime.
#   members: unit indices of the regime; nb: list of neighbour indices (all units)
#   sizes: admissible region sizes; returns a data frame of candidates with a
#   list column of members.
scan_regime_regions <- function(sc, members, nb, sizes, anchors, beam_width = 4L,
                                ridge = 0.01, channels = c("intercept", "slope"),
                                keep_top = 200L) {
  inside <- logical(length(nb))
  inside[members] <- TRUE
  max_size <- max(sizes)
  # Order-free fingerprints of a unit set (two random sums) replace string
  # keys; a collision of both sums has probability ~ 0.
  h1 <- stats::runif(length(nb)); h2 <- stats::runif(length(nb))
  fp <- function(mem) sum(h1[mem]) + 1e3 * sum(h2[mem])
  found <- list(); found_fp <- numeric(); found_st <- matrix(0, 0, 2)
  prune <- function() {
    # Keep the keep_top best regions of each channel.
    keep <- unique(c(utils::head(order(-found_st[, 1]), keep_top), utils::head(order(-found_st[, 2]), keep_top)))
    found <<- found[keep]; found_fp <<- found_fp[keep]; found_st <<- found_st[keep, , drop = FALSE]
  }
  record <- function(mem, stats) {
    f <- fp(mem)
    if (length(found_fp) && any(abs(found_fp - f) < 1e-9)) return(invisible(NULL))
    found[[length(found) + 1L]] <<- sort(mem)
    found_fp <<- c(found_fp, f)
    found_st <<- rbind(found_st, stats)
    if (length(found) > 4L * keep_top) prune()
    invisible(NULL)
  }
  for (a in anchors) for (channel in channels) {
    beam <- list(list(mem = a, su = sc$U[a, ], sv = sc$V[a, ]))
    for (size in 2:max_size) {
      cand_mem <- list(); cand_par <- integer(); cand_v <- integer(); cand_fp <- numeric()
      for (bi in seq_along(beam)) {
        b <- beam[[bi]]
        fr <- unique(unlist(nb[b$mem], use.names = FALSE))
        fr <- fr[inside[fr] & !(fr %in% b$mem)]
        if (!length(fr)) next
        base_fp <- fp(b$mem)
        cand_par <- c(cand_par, rep(bi, length(fr)))
        cand_v <- c(cand_v, fr)
        cand_fp <- c(cand_fp, base_fp + h1[fr] + 1e3 * h2[fr])
      }
      if (!length(cand_v)) break
      uq <- !duplicated(round(cand_fp, 9))
      cand_par <- cand_par[uq]; cand_v <- cand_v[uq]
      SU <- do.call(rbind, lapply(seq_along(beam), function(bi) beam[[bi]]$su))[cand_par, , drop = FALSE] +
        sc$U[cand_v, , drop = FALSE]
      SV <- do.call(rbind, lapply(seq_along(beam), function(bi) beam[[bi]]$sv))[cand_par, , drop = FALSE] +
        sc$V[cand_v, , drop = FALSE]
      st <- scan_channel_stats(SU, SV, sc$d, sc$K, ridge)
      keep <- utils::head(order(-st[, channel]), beam_width)
      beam <- lapply(keep, function(i) list(mem = c(beam[[cand_par[i]]]$mem, cand_v[i]), su = SU[i, ], sv = SV[i, ]))
      # Record the best region of this size for this anchor and channel.
      if (size %in% sizes) record(beam[[1L]]$mem, st[keep[1L], ])
    }
  }
  if (length(found)) prune()
  lapply(seq_along(found), function(i) list(members = found[[i]], intercept = found_st[i, 1],
                                            slope = found_st[i, 2]))
}

# Choose up to n regions alternating between channels, with a Jaccard limit.
scan_select_regions <- function(found, n, max_jaccard = 0.6, channels = c("intercept", "slope")) {
  if (!length(found)) return(list())
  sc <- sapply(channels, function(ch) vapply(found, function(z) z[[ch]], numeric(1)))
  if (is.null(dim(sc))) sc <- matrix(sc, nrow = 1L, dimnames = list(NULL, channels))
  orders <- lapply(channels, function(ch) order(-sc[, ch]))
  pos <- rep(1L, length(channels))
  chosen <- integer(); chosen_ch <- character()
  jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  while (length(chosen) < n && any(pos <= length(found))) {
    for (k in seq_along(channels)) {
      if (length(chosen) >= n) break
      while (pos[k] <= length(found)) {
        i <- orders[[k]][pos[k]]
        pos[k] <- pos[k] + 1L
        if (i %in% chosen) next
        if (length(chosen) && any(vapply(chosen, function(j)
          jacc(found[[i]]$members, found[[j]]$members) > max_jaccard, logical(1)))) next
        chosen <- c(chosen, i); chosen_ch <- c(chosen_ch, channels[k])
        break
      }
    }
  }
  lapply(seq_along(chosen), function(r) c(found[[chosen[r]]], list(channel = chosen_ch[r],
                                                                   score = sc[chosen[r], chosen_ch[r]])))
}
