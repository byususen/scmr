# Copied unchanged from the Stage 26 working package (Bahan/Syntax/scmr_pkg/R) for the scan comparison.
# Connected, flexible candidate starts. The scores rank starts only; they are
# neither calibrated spatial tests nor criteria for selecting G.

#' Flexible connected starts for multinomial SCMR
#'
#' Combines BFS prefixes, full-grid rectangles, and bounded score-guided graph
#' growth. The candidate set and all screening choices must be repeated during
#' any later null calibration of the final model-selection rule.
#'
#' @inheritParams scmr_score_scan_starts
#' @param beam_width Partial regions retained per anchor and score channel.
#' @param max_jaccard Maximum Jaccard overlap between selected starts.
#' @param compute_whitened Compute optional ridge-whitened slope scores.
#' @param compute_full_intercept Compute optional full-class intercept scores.
#' @param include_motifs Enumerate connected subsets of local nearest units.
#' @param motif_neighbors Number of local units considered around each anchor.
#' @param motif_sizes Sizes enumerated by the local subset search.
#' @return List with named starts, selected candidates, all feasible regions,
#'   and candidate counts. All regions contain unit IDs, not row indices.
#' @export
scmr_flexible_scan_starts <- function(x, y, unit_id, coords, global_fit,
                                      component = NULL, k_neighbors = 4L,
                                      sizes = 5:12, top_k = 5L,
                                      min_per_class = 5L, ridge = 0.01,
                                      beam_width = 4L, max_jaccard = 0.6,
                                      compute_whitened = FALSE,
                                      compute_full_intercept = FALSE,
                                      include_motifs = FALSE,
                                      motif_neighbors = 9L,
                                      motif_sizes = 5:8) {
  data <- prepare_scmr_data(x, y, unit_id, coords, component)
  if (!inherits(global_fit, "scmr_fit") || global_fit$G != 1L ||
      nrow(global_fit$prob_train) != nrow(data$x) ||
      !identical(as.character(global_fit$y_train), as.character(data$y)) ||
      !identical(global_fit$x_colnames, colnames(data$x)) ||
      !identical(as.character(global_fit$unit_id), data$unit_id) ||
      !identical(colnames(global_fit$prob_train), levels(data$y)))
    stop("global_fit must be a one-cluster fit on these ordered rows and classes.",
         call. = FALSE)
  positive_integer <- function(v) length(v) == 1L && is.finite(v) &&
    v >= 1L && v == floor(v)
  if (!positive_integer(k_neighbors) || !positive_integer(top_k) ||
      !positive_integer(beam_width) || !positive_integer(motif_neighbors))
    stop("k_neighbors, top_k, beam_width and motif_neighbors must be positive integers.",
         call. = FALSE)
  if (length(ridge) != 1L || !is.finite(ridge) || ridge <= 0 ||
      length(max_jaccard) != 1L || !is.finite(max_jaccard) ||
      max_jaccard < 0 || max_jaccard > 1)
    stop("ridge must be positive and max_jaccard in [0, 1].", call. = FALSE)
  if (!is.logical(compute_whitened) || length(compute_whitened) != 1L ||
      is.na(compute_whitened) || !is.logical(compute_full_intercept) ||
      length(compute_full_intercept) != 1L || is.na(compute_full_intercept))
    stop("Optional score switches must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(include_motifs) || length(include_motifs) != 1L ||
      is.na(include_motifs))
    stop("include_motifs must be TRUE or FALSE.", call. = FALSE)
  units <- sort(unique(data$unit_id))
  B <- length(units)
  sizes <- sort(unique(as.integer(sizes)))
  sizes <- sizes[!is.na(sizes) & sizes >= 2L & sizes <= B - 2L]
  if (!length(sizes)) stop("No feasible candidate sizes.", call. = FALSE)
  motif_sizes <- intersect(sizes, unique(as.integer(motif_sizes)))
  motif_sizes <- motif_sizes[!is.na(motif_sizes) & motif_sizes <= motif_neighbors]
  min_class <- resolve_min_class(min_per_class, levels(data$y))
  coords_u <- get_unit_coords(data$unit_id, data$coords)
  unit_rows <- split(seq_len(nrow(data$x)),
                     factor(match(data$unit_id, units), levels = seq_along(units)))
  component_u <- if (is.null(data$component)) rep("all", B) else
    vapply(unit_rows, function(idx) unique(data$component[idx]), character(1L))
  A <- as.matrix(build_unit_weights(coords_u, k = k_neighbors,
                                     component = component_u)$W) != 0
  diag(A) <- FALSE
  distance <- as.matrix(stats::dist(coords_u))
  K <- nlevels(data$y) - 1L
  d <- ncol(data$x) + 1L
  U <- matrix(0, B, K * d)
  V <- matrix(0, B, K * d)
  # A ridge-whitened predictor basis is computed only from X. Its scores are
  # diagnostic alternatives to the raw slope ranking, not inferential tests.
  if (compute_whitened) {
    S <- stats::cov(data$x)
    eig <- eigen((S + t(S)) / 2, symmetric = TRUE)
    whiten <- eig$vectors %*%
      (diag(1 / sqrt(pmax(eig$values, 0) + 0.1), nrow = ncol(data$x))) %*%
      t(eig$vectors)
    x_white <- data$x %*% whiten
    Uw <- matrix(0, B, K * d)
    Vw <- matrix(0, B, K * d)
  }
  if (compute_full_intercept) H <- array(0, dim = c(B, K, K))
  counts <- matrix(0L, B, nlevels(data$y))
  p0 <- global_fit$prob_train
  for (b in seq_len(B)) {
    idx <- unit_rows[[b]]
    z <- cbind(1, data$x[idx, , drop = FALSE])
    if (compute_whitened) zw <- cbind(1, x_white[idx, , drop = FALSE])
    yy <- as.integer(data$y[idx])
    counts[b, ] <- tabulate(yy, nbins = nlevels(data$y))
    if (compute_full_intercept) {
      pp <- p0[idx, seq_len(K), drop = FALSE]
      H[b, , ] <- diag(colSums(pp), nrow = K) - crossprod(pp)
    }
    for (c in seq_len(K)) {
      jj <- ((c - 1L) * d + 1L):(c * d)
      r <- as.numeric(yy == c) - p0[idx, c]
      v <- p0[idx, c] * (1 - p0[idx, c])
      U[b, jj] <- colSums(z * r)
      V[b, jj] <- colSums(z^2 * v)
      if (compute_whitened) {
        Uw[b, jj] <- colSums(zw * r)
        Vw[b, jj] <- colSums(zw^2 * v)
      }
    }
  }
  total <- colSums(counts)
  intercept_idx <- seq.int(1L, K * d, by = d)
  if (compute_full_intercept) {
    Htotal <- apply(H, c(2L, 3L), sum)
    Htotal_inv <- solve(Htotal + diag(ridge, K))
  }
  score <- function(members, full = FALSE) {
    u <- colSums(U[members, , drop = FALSE])
    v <- colSums(V[members, , drop = FALSE])
    fs <- u^2 / (v + ridge)
    slopes <- matrix(fs[-intercept_idx], nrow = d - 1L)
    slope_white <- NA_real_
    if (compute_whitened) {
      uw <- colSums(Uw[members, , drop = FALSE])
      vw <- colSums(Vw[members, , drop = FALSE])
      fw <- uw^2 / (vw + ridge)
      slopes_white <- matrix(fw[-intercept_idx], nrow = d - 1L)
      slope_white <- max(rowSums(slopes_white))
    }
    out <- c(intercept = sum(fs[intercept_idx]),
      slope = max(rowSums(slopes)), joint = sum(fs),
      slope_white = slope_white, intercept_full = NA_real_)
    if (full && compute_full_intercept) {
      Hwin <- apply(H[members, , , drop = FALSE], c(2L, 3L), sum)
      # Profile out the common intercept. The global slopes and penalty are
      # still treated as fixed, so this remains a ranking heuristic.
      Heff <- Hwin - Hwin %*% Htotal_inv %*% Hwin
      zint <- u[intercept_idx]
      out[["intercept_full"]] <- as.numeric(crossprod(
        zint, solve(Heff + diag(ridge, K), zint)))
    }
    out
  }
  seen <- new.env(hash = TRUE, parent = emptyenv())
  regions <- list()
  scores <- list()
  anchors <- character()
  origins <- character()
  add <- function(members, anchor, origin) {
    members <- sort(unique(members))
    key <- paste(members, collapse = ",")
    if (exists(key, envir = seen, inherits = FALSE)) return(invisible(NULL))
    assign(key, TRUE, envir = seen)
    inside <- colSums(counts[members, , drop = FALSE])
    if (any(inside < min_class) || any(total - inside < min_class))
      return(invisible(NULL))
    regions[[length(regions) + 1L]] <<- members
    scores[[length(scores) + 1L]] <<- score(members, full = TRUE)
    anchors <<- c(anchors, units[anchor])
    origins <<- c(origins, origin)
    invisible(NULL)
  }
  connected <- function(members) {
    reached <- members[1L]
    repeat {
      fresh <- setdiff(members[vapply(members, function(i)
        any(A[i, reached]), logical(1L))], reached)
      if (!length(fresh)) break
      reached <- c(reached, fresh)
    }
    length(reached) == length(members)
  }
  # The historical BFS prefixes are retained as an exact ablation baseline.
  for (start in seq_len(B)) {
    hops <- rep(Inf, B)
    hops[start] <- 0
    queue <- start
    head <- 1L
    while (head <= length(queue)) {
      u <- queue[head]
      head <- head + 1L
      fresh <- which(A[u, ] & is.infinite(hops))
      if (length(fresh)) {
        hops[fresh] <- hops[u] + 1
        queue <- c(queue, fresh)
      }
    }
    ord <- order(hops, distance[start, ], seq_len(B))
    ord <- ord[is.finite(hops[ord])]
    for (size in sizes[sizes <= length(ord)])
      add(ord[seq_len(size)], start, "bfs")
  }
  # Rectangle enumeration is restricted to complete, equally spaced grids.
  # This prevents accidental rectangles on irregular application geometry.
  for (cc in unique(component_u)) {
    idx <- which(component_u == cc)
    xx <- sort(unique(coords_u[idx, 1L]))
    yy <- sort(unique(coords_u[idx, 2L]))
    regular <- function(v) length(v) < 3L ||
      max(abs(diff(v) - median(diff(v)))) < 1e-8
    if (length(xx) * length(yy) != length(idx) ||
        !regular(xx) || !regular(yy) ||
        anyDuplicated(paste(coords_u[idx, 1L], coords_u[idx, 2L]))) next
    for (x1 in seq_along(xx)) for (x2 in x1:length(xx))
      for (y1 in seq_along(yy)) for (y2 in y1:length(yy)) {
        size <- (x2 - x1 + 1L) * (y2 - y1 + 1L)
        if (!size %in% sizes) next
        members <- idx[coords_u[idx, 1L] >= xx[x1] &
                       coords_u[idx, 1L] <= xx[x2] &
                       coords_u[idx, 2L] >= yy[y1] &
                       coords_u[idx, 2L] <= yy[y2]]
        if (length(members) == size && connected(members))
          add(members, members[1L], "rectangle")
      }
  }
  # Search beyond rectangles, retaining a small beam for each seed and channel.
  for (start in seq_len(B)) for (channel in c("intercept", "slope")) {
    beam <- list(start)
    for (size in 2:max(sizes)) {
      grown <- list()
      local_seen <- new.env(hash = TRUE, parent = emptyenv())
      for (path in beam) {
        frontier <- setdiff(which(colSums(A[path, , drop = FALSE]) > 0), path)
        for (v in frontier) {
          members <- sort(c(path, v))
          key <- paste(members, collapse = ",")
          if (exists(key, envir = local_seen, inherits = FALSE)) next
          assign(key, TRUE, envir = local_seen)
          grown[[length(grown) + 1L]] <- members
        }
      }
      if (!length(grown)) break
      value <- vapply(grown, function(z) score(z)[[channel]], numeric(1L))
      beam <- grown[head(order(-value, seq_along(value)), beam_width)]
      if (size %in% sizes)
        for (members in beam) add(members, start, "beam")
    }
  }
  if (include_motifs && length(motif_sizes)) {
    # Enumerate local connected sets independent of the response score.
    # The fixed geometric neighborhood prevents a single greedy path from
    # eliminating a compact irregular shape before it reaches target size.
    for (start in seq_len(B)) {
      same <- which(component_u == component_u[start])
      local <- head(same[order(distance[start, same], same)], motif_neighbors)
      for (size in motif_sizes[motif_sizes <= length(local)]) {
        utils::combn(local, size, FUN = function(members) {
          if (connected(members)) add(members, start, "motif")
          NULL
        }, simplify = FALSE)
      }
    }
  }
  if (!length(regions)) stop("No feasible flexible-scan region.", call. = FALSE)
  sm <- do.call(rbind, scores)
  pool <- list(regions = lapply(regions, function(i) units[i]),
               origins = origins, anchors = anchors,
               all_scores = sm, unit_ids = units,
               n_candidates = length(regions), n_by_origin = table(origins))
  c(pool, scmr_select_flexible_starts(pool, top_k = top_k,
                                     max_jaccard = max_jaccard))
}

#' Rank a fixed pool of flexible connected starts
#'
#' This permits ranking ablations without changing the candidate pool or
#' refitting the global model. Within-size ranks are empirical screening
#' normalizations, not p-values or calibrated model-selection criteria.
#'
#' @param scan Output from `scmr_flexible_scan_starts()`.
#' @param top_k Number of starts to select.
#' @param max_jaccard Maximum Jaccard overlap among selected regions.
#' @param slope_score Raw or ridge-whitened predictor score.
#' @param intercept_score Diagonal, full covariance, or local offset gain.
#' @param size_calibration Raw score or within-size empirical rank.
#' @return A list with `starts` and `candidates`.
#' @export
scmr_select_flexible_starts <- function(scan, top_k = 5L,
                                        max_jaccard = 0.6,
                                        slope_score = c("raw", "whitened"),
                                        intercept_score = c("diagonal", "full",
                                                            "local_likelihood"),
                                        size_calibration = c("none", "within_size")) {
  slope_score <- match.arg(slope_score)
  intercept_score <- match.arg(intercept_score)
  size_calibration <- match.arg(size_calibration)
  if (length(top_k) != 1L || !is.finite(top_k) || top_k < 1L ||
      top_k != floor(top_k) || length(max_jaccard) != 1L ||
      !is.finite(max_jaccard) || max_jaccard < 0 || max_jaccard > 1)
    stop("Invalid top_k or max_jaccard.", call. = FALSE)
  required <- c("regions", "all_scores", "unit_ids", "anchors", "origins")
  if (!all(required %in% names(scan)))
    stop("scan must contain a flexible candidate pool.", call. = FALSE)
  regions <- scan$regions
  sm <- scan$all_scores
  sizes <- lengths(regions)
  slope_column <- if (slope_score == "raw") "slope" else "slope_white"
  intercept_column <- switch(intercept_score, diagonal = "intercept",
                             full = "intercept_full",
                             local_likelihood = "intercept_local")
  if (!intercept_column %in% colnames(sm))
    stop("scan lacks the requested intercept score.", call. = FALSE)
  if (any(!is.finite(sm[, intercept_column])) ||
      any(!is.finite(sm[, slope_column])))
    stop("Requested optional scores were not computed for this pool.",
         call. = FALSE)
  rank_scores <- cbind(intercept = sm[, intercept_column],
                       slope = sm[, slope_column], joint = sm[, "joint"])
  if (size_calibration == "within_size") {
    for (size in unique(sizes)) {
      idx <- which(sizes == size)
      for (channel in colnames(rank_scores)) {
        percentile <- (rank(rank_scores[idx, channel], ties.method = "average") -
                       0.5) / length(idx)
        rank_scores[idx, channel] <- stats::qnorm(percentile)
      }
    }
  }
  diverse <- function(i, selected) {
    if (!length(selected)) return(TRUE)
    all(vapply(selected, function(j)
      length(intersect(regions[[i]], regions[[j]])) /
        length(union(regions[[i]], regions[[j]])) <= max_jaccard,
      logical(1L)))
  }
  ranked <- integer()
  channels <- character()
  append_rank <- function(ordering, channel, quota) {
    for (i in ordering) {
      if (length(ranked) >= quota) break
      if (i %in% ranked || !diverse(i, ranked)) next
      ranked <<- c(ranked, i)
      channels <<- c(channels, channel)
    }
  }
  append_rank(order(-rank_scores[, "intercept"], seq_len(nrow(sm))), "intercept",
              max(1L, floor(top_k / 2L)))
  append_rank(order(-rank_scores[, "slope"], seq_len(nrow(sm))), "slope", top_k)
  append_rank(order(-rank_scores[, "joint"], seq_len(nrow(sm))), "joint", top_k)
  if (length(ranked) < top_k)
    warning("Fewer diverse feasible starts than top_k.", call. = FALSE)
  starts <- lapply(ranked, function(i) {
    labels <- rep(1L, length(scan$unit_ids))
    labels[match(regions[[i]], scan$unit_ids)] <- 2L
    stats::setNames(labels, scan$unit_ids)
  })
  selected_score <- vapply(seq_along(ranked), function(j)
    rank_scores[ranked[j], channels[j]], numeric(1L))
  candidates <- data.frame(Rank = seq_along(ranked),
    Anchor = scan$anchors[ranked], Units = lengths(regions[ranked]),
    Channel = channels, Origin = scan$origins[ranked], Score = selected_score,
    ScoreIntercept = sm[ranked, intercept_column],
    ScoreSlope = sm[ranked, slope_column], ScoreJoint = sm[ranked, "joint"],
    InterceptScore = intercept_score, SlopeScore = slope_score,
    SizeCalibration = size_calibration,
    MaxJaccard = max_jaccard)
  list(starts = starts, candidates = candidates)
}
