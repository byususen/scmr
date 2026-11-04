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
                                ridge = 0.01, channels = c("intercept", "slope")) {
  inside <- logical(length(nb))
  inside[members] <- TRUE
  max_size <- max(sizes)
  found <- list()
  keys <- new.env(hash = TRUE, parent = emptyenv())
  record <- function(mem, stats, origin) {
    key <- paste(sort(mem), collapse = ",")
    if (exists(key, envir = keys, inherits = FALSE)) return(invisible(NULL))
    assign(key, TRUE, envir = keys)
    found[[length(found) + 1L]] <<- list(members = sort(mem), intercept = stats[1], slope = stats[2],
                                         origin = origin)
  }
  for (a in anchors) for (channel in channels) {
    beam <- list(list(mem = a, su = sc$U[a, ], sv = sc$V[a, ]))
    for (size in 2:max_size) {
      cand_mem <- list(); cand_su <- list(); cand_sv <- list()
      seen <- new.env(hash = TRUE, parent = emptyenv())
      for (b in beam) {
        fr <- unique(unlist(nb[b$mem], use.names = FALSE))
        fr <- fr[inside[fr] & !(fr %in% b$mem)]
        if (!length(fr)) next
        for (v in fr) {
          key <- paste(sort(c(b$mem, v)), collapse = ",")
          if (exists(key, envir = seen, inherits = FALSE)) next
          assign(key, TRUE, envir = seen)
          cand_mem[[length(cand_mem) + 1L]] <- c(b$mem, v)
          cand_su[[length(cand_su) + 1L]] <- b$su + sc$U[v, ]
          cand_sv[[length(cand_sv) + 1L]] <- b$sv + sc$V[v, ]
        }
      }
      if (!length(cand_mem)) break
      SU <- do.call(rbind, cand_su)
      SV <- do.call(rbind, cand_sv)
      st <- scan_channel_stats(SU, SV, sc$d, sc$K, ridge)
      keep <- utils::head(order(-st[, channel]), beam_width)
      beam <- lapply(keep, function(i) list(mem = cand_mem[[i]], su = SU[i, ], sv = SV[i, ]))
      if (size %in% sizes) for (i in keep) record(cand_mem[[i]], st[i, ], channel)
    }
  }
  found
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
