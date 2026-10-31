# Split-merge path initialization over G = 1, ..., G_max.
#
# The fit at G starts from the fit at G - 1 with one regime split in two. The
# split is proposed from block score vectors: for unit b in regime g,
#   s_b = (1 / n_b) sum_{i in b} (e_{y_i} - p_g(x_i)) (x) (1, x_i),
# the gradient of regime g's log-likelihood contributed by b, randomly projected
# to a few dimensions and smoothed over spatial neighbours in the same regime.
# At the fitted coefficients the scores of a regime sum to (nearly) zero, so a
# 2-means split of the scores separates the sub-regions that pull the regime's
# coefficients in opposite directions. After the forward pass an optional
# backward pass merges adjacent regimes whose cross-prediction loss is smallest
# and keeps the merged fit at G - 1 when its penalized objective is larger.
# Every candidate is a complete SCMR fit started from the proposed partition, so
# each accepted fit is a fixed point of the same monotone alternation.

#' Split-merge path of SCMR fits over the number of regimes
#'
#' Fits SCMR for G = 1, ..., `G_max`, each G warm-started from the best fit at
#' G - 1 with one regime split along its block score vectors (see Details).
#' Optionally the usual multi-start fit is also run at every G and the better of
#' the two is kept, and a backward pass tries merges of adjacent regimes.
#'
#' The penalized objectives of different fits at the same G are compared, so
#' they must use the same penalty: supply `lambda` and use
#' `lambda_scale = "sum"` in `control` (the common sum-scale strength is then
#' kappa = lambda * N / G for every fit at that G).
#'
#' @inheritParams fit_scmr
#' @param G_max Largest number of regimes.
#' @param n_split Number of regimes (ranked by the between-group sum of squares
#'   of their 2-means score split) tried as split candidates at each G.
#' @param n_merge Number of adjacent regime pairs tried in the backward pass.
#' @param proj_dim Dimension of the random projection of the score vectors.
#' @param smooth If `TRUE`, average each unit's score with its spatial
#'   neighbours in the same regime before the 2-means split.
#' @param merge_back If `TRUE`, run the backward merge pass.
#' @param standard_starts If `TRUE`, also run the multi-start fit of `control`
#'   at every G and keep it when its objective is larger.
#' @param init_features Optional coefficient features for the standard starts.
#' @param keep_fits If `FALSE`, keep only memberships and criteria per G.
#' @param verbose Print progress.
#' @return An object of class `scmr_path`: `fits` (list by G), `table` (one row
#'   per G with the objective, the information criteria and the origin of the
#'   kept fit), `log` (every candidate fit), `selected` (G chosen by each
#'   criterion) and `seconds`.
#' @export
scmr_split_merge_path <- function(x, y, unit_id, coords, G_max,
                                  penalty = c("elastic_net", "lasso", "ridge", "none"),
                                  alpha = NULL, lambda = NULL, control = scmr_control(),
                                  seed = 123, n_split = 2L, n_merge = 2L, proj_dim = 24L,
                                  smooth = TRUE, merge_back = TRUE, standard_starts = FALSE,
                                  init_features = NULL, keep_fits = TRUE, verbose = FALSE) {
  penalty <- match.arg(penalty)
  t_start <- proc.time()[["elapsed"]]
  if (penalty != "none" && (is.null(lambda) || !identical(control$lambda_scale, "sum"))) {
    warning("Objectives are comparable across candidate fits only with a supplied lambda ",
            "and lambda_scale = \"sum\".", call. = FALSE)
  }
  data <- prepare_scmr_data(x, y, unit_id, coords)
  units <- sort(unique(data$unit_id))
  m <- length(units)
  unit_num <- match(data$unit_id, units)
  rows <- split(seq_len(nrow(data$x)), factor(unit_num, levels = seq_len(m)))
  contrib <- precompute_unit_contributions(rows, data$y, levels(data$y))
  min_class <- resolve_min_class(control$min_per_class, levels(data$y))
  feasible <- function(ug, G) length(unique(ug)) == G &&
    is_feasible_partition(ug, G, contrib$unit_total, contrib$unit_class, control$min_units, min_class)
  G_max <- as.integer(min(G_max, floor(m / control$min_units)))
  if (G_max < 1L) stop("G_max must allow at least one regime of min_units units.", call. = FALSE)

  ctrl1 <- control
  ctrl1$n_starts <- 1L
  ctrl1$start_methods <- NULL
  fit_from <- function(G, init) {
    t0 <- proc.time()[["elapsed"]]
    f <- tryCatch(fit_scmr(x, y, model = "scmr", penalty = penalty, G = G, unit_id = unit_id,
                           coords = coords, alpha = alpha, lambda = lambda, control = ctrl1,
                           initial_cluster = init, seed = seed),
                  error = function(e) e)
    if (inherits(f, "error")) return(list(error = conditionMessage(f), seconds = proc.time()[["elapsed"]] - t0))
    f$path_seconds <- proc.time()[["elapsed"]] - t0
    f
  }
  objective <- function(f) {
    st <- f$diagnostics$starts
    # G = 1 is fitted without the partition search; its penalized objective is
    # the log-likelihood minus the penalty.
    if (is.null(st) || !nrow(st)) return(f$criteria$PenalizedObjective)
    st$Objective[which(st$Selected)[1]]
  }
  unit_groups <- function(f) as.integer(f$group_unit[units])

  # Random projection of (class x predictor) score vectors; fixed for the path.
  xs <- scale(data$x)
  xs[!is.finite(xs)] <- 0
  X1 <- cbind(1, xs)
  K <- nlevels(data$y)
  Y <- matrix(0, nrow(data$x), K)
  Y[cbind(seq_len(nrow(data$x)), as.integer(data$y))] <- 1
  d_full <- ncol(X1) * K
  proj <- d_full > proj_dim
  R <- if (proj) with_scmr_seed(scmr_seed(seed, 5150), array(stats::rnorm(ncol(X1) * K * proj_dim),
                                                             c(ncol(X1), K, proj_dim)) / sqrt(proj_dim)) else NULL
  unit_scores <- function(f) {
    E <- Y - f$prob_train
    Z <- if (proj) {
      Reduce(`+`, lapply(seq_len(K), function(k) E[, k] * (X1 %*% R[, k, ])))
    } else {
      do.call(cbind, lapply(seq_len(K), function(k) E[, k] * X1))
    }
    rowsum(Z, unit_num, reorder = TRUE) / contrib$unit_total
  }
  smooth_scores <- function(S, ug, W) {
    if (!smooth || is.null(W)) return(S)
    e <- Matrix::summary(methods::as(W, "generalMatrix"))
    keep <- ug[e$i] == ug[e$j]
    Wg <- Matrix::sparseMatrix(i = e$i[keep], j = e$j[keep], x = e$x[keep], dims = c(m, m))
    as.matrix((S + Wg %*% S) / (1 + Matrix::rowSums(Wg)))
  }
  split_proposals <- function(f, G_new) {
    ug <- unit_groups(f)
    S <- smooth_scores(unit_scores(f), ug, f$W)
    G_old <- G_new - 1L
    props <- list()
    for (g in seq_len(G_old)) {
      idx <- which(ug == g)
      if (length(idx) < 2L * control$min_units) next
      Sg <- S[idx, , drop = FALSE]
      km <- tryCatch(with_scmr_seed(scmr_seed(seed, 6000 + 37L * g + G_new),
                                    stats::kmeans(Sg, centers = 2L, nstart = 5L, iter.max = 50L)),
                     error = function(e) NULL)
      cand <- list()
      if (!is.null(km)) {
        lab <- ug; lab[idx[km$cluster == 2L]] <- G_new
        cand$kmeans <- list(labels = lab, score = km$betweenss)
      }
      # Fallback when the 2-means split violates min_units/min_per_class:
      # the median cut of the leading principal direction.
      pc <- tryCatch(stats::prcomp(Sg, center = TRUE)$x[, 1], error = function(e) NULL)
      if (!is.null(pc)) {
        lab <- ug; lab[idx[pc > stats::median(pc)]] <- G_new
        cand$pc_median <- list(labels = lab, score = if (is.null(km)) sum(pc^2) else km$betweenss)
      }
      for (nm in names(cand)) if (feasible(cand[[nm]]$labels, G_new)) {
        props[[length(props) + 1L]] <- list(regime = g, rule = nm, labels = cand[[nm]]$labels,
                                            score = cand[[nm]]$score)
        break
      }
    }
    if (!length(props)) return(props)
    props[order(-vapply(props, function(z) z$score, numeric(1)))][seq_len(min(n_split, length(props)))]
  }
  merge_proposals <- function(f) {
    ug <- unit_groups(f)
    G <- f$G
    if (G < 2L || is.null(f$W)) return(list())
    H <- Matrix::sparseMatrix(i = seq_len(m), j = ug, x = 1, dims = c(m, G))
    A <- as.matrix(Matrix::t(H) %*% f$W %*% H)
    lev <- levels(data$y)
    own <- vapply(seq_len(G), function(g) sum(true_class_logp(f$prob_train[f$row_group == g, , drop = FALSE],
                                                              data$y[f$row_group == g], lev)), numeric(1))
    cross <- function(g, h) {           # log-likelihood of regime g's rows under regime h's fit
      idx <- which(f$row_group == g)
      sum(true_class_logp(engine_predict(f$fits[[h]], data$x[idx, , drop = FALSE]), data$y[idx], lev))
    }
    pairs <- which(upper.tri(A) & A > 0, arr.ind = TRUE)
    if (!nrow(pairs)) return(list())
    loss <- apply(pairs, 1L, function(p) {
      g <- p[1]; h <- p[2]
      min(own[g] - cross(g, h), own[h] - cross(h, g))
    })
    ord <- order(loss)[seq_len(min(n_merge, nrow(pairs)))]
    lapply(ord, function(r) {
      g <- pairs[r, 1]; h <- pairs[r, 2]
      lab <- ug; lab[lab == h] <- g
      list(pair = c(g, h), labels = as.integer(factor(lab)), loss = loss[r])
    })
  }

  log <- data.frame()
  add_log <- function(G, step, detail, f) {
    ok <- is.null(f$error)
    log <<- rbind(log, data.frame(G = G, Step = step, Detail = detail,
      Objective = if (ok) objective(f) else NA_real_,
      Seconds = if (ok) f$path_seconds else f$seconds,
      Error = if (ok) "" else f$error, stringsAsFactors = FALSE))
  }
  say <- function(...) if (verbose) message(...)
  best <- vector("list", G_max)
  origin <- character(G_max)

  # Forward pass.
  f1 <- fit_from(1L, rep(1L, m))
  if (!is.null(f1$error)) stop("G = 1 fit failed: ", f1$error, call. = FALSE)
  add_log(1L, "start", "single regime", f1)
  best[[1L]] <- f1
  origin[1L] <- "single"
  if (G_max >= 2L) for (G in 2L:G_max) {
    props <- split_proposals(best[[G - 1L]], G)
    for (p in props) {
      f <- fit_from(G, p$labels)
      add_log(G, "split", sprintf("regime %d (%s)", p$regime, p$rule), f)
      if (is.null(f$error) && (is.null(best[[G]]) || objective(f) > objective(best[[G]]))) {
        best[[G]] <- f
        origin[G] <- sprintf("split %d", p$regime)
      }
    }
    if (standard_starts) {
      t0 <- proc.time()[["elapsed"]]
      f <- tryCatch(fit_scmr(x, y, model = "scmr", penalty = penalty, G = G, unit_id = unit_id,
                             coords = coords, alpha = alpha, lambda = lambda, control = control,
                             seed = seed, init_features = init_features), error = function(e) e)
      f <- if (inherits(f, "error")) list(error = conditionMessage(f), seconds = proc.time()[["elapsed"]] - t0) else
        c(f, list(path_seconds = proc.time()[["elapsed"]] - t0))
      if (is.null(f$error)) class(f) <- class(best[[1L]])
      add_log(G, "standard", "multi-start", f)
      if (is.null(f$error) && (is.null(best[[G]]) || objective(f) > objective(best[[G]]))) {
        best[[G]] <- f
        origin[G] <- "standard"
      }
    }
    if (is.null(best[[G]])) {
      say("No feasible split at G = ", G, "; path stops at G = ", G - 1L)
      G_max <- G - 1L
      best <- best[seq_len(G_max)]
      origin <- origin[seq_len(G_max)]
      break
    }
    say(sprintf("G = %d: objective %.2f (%s)", G, objective(best[[G]]), origin[G]))
  }

  # Backward pass.
  if (merge_back && G_max >= 3L) for (G in G_max:3L) {
    for (p in merge_proposals(best[[G]])) {
      f <- fit_from(G - 1L, p$labels)
      add_log(G - 1L, "merge", sprintf("regimes %d+%d of G = %d", p$pair[1], p$pair[2], G), f)
      if (is.null(f$error) && objective(f) > objective(best[[G - 1L]]) + 1e-8) {
        best[[G - 1L]] <- f
        origin[G - 1L] <- sprintf("merge from G = %d", G)
        say(sprintf("G = %d improved by merging: objective %.2f", G - 1L, objective(f)))
      }
    }
  }

  crit_names <- c("CriterionPLIC_AIC", "CriterionPLIC_BIC", "CriterionPLIC_BIC_rows", "CriterionPLIC_BIC_units")
  table <- do.call(rbind, lapply(seq_len(G_max), function(G) {
    f <- best[[G]]
    cr <- f$criteria
    data.frame(G = G, Origin = origin[G], Objective = objective(f), Phi = f$phi %||% NA_real_,
               MinUnits = min(tabulate(unit_groups(f), G)), DF = cr$df_effective_total,
               PLIC_AIC = cr$CriterionPLIC_AIC, PLIC_BIC = cr$CriterionPLIC_BIC,
               PLIC_BIC_rows = cr$CriterionPLIC_BIC_rows %||% NA_real_,
               PLIC_BIC_units = cr$CriterionPLIC_BIC_units %||% NA_real_,
               stringsAsFactors = FALSE)
  }))
  selected <- vapply(c("PLIC_AIC", "PLIC_BIC", "PLIC_BIC_rows", "PLIC_BIC_units"),
                     function(k) table$G[which.min(table[[k]])], numeric(1))
  if (!keep_fits) best <- lapply(best, function(f) list(G = f$G, group_unit = f$group_unit,
                                                        criteria = f$criteria[setdiff(names(f$criteria), "cluster_components")],
                                                        phi = f$phi))
  structure(list(fits = best, table = table, log = log, selected = selected,
                 seconds = proc.time()[["elapsed"]] - t_start, units = units),
            class = "scmr_path")
}

#' @export
print.scmr_path <- function(x, ...) {
  cat("SCMR split-merge path, G = 1..", nrow(x$table), " (", round(x$seconds, 1), " s, ",
      nrow(x$log), " candidate fits)\n", sep = "")
  print(x$table, row.names = FALSE, digits = 6)
  cat("Selected G:", paste(names(x$selected), x$selected, sep = " = ", collapse = ", "), "\n")
  invisible(x)
}
