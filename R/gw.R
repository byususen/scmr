# Geographically weighted multinomial elastic net (comparator).
#
# At anchor units a multinomial elastic net is fitted with bisquare kernel
# weights over the k nearest units (adaptive bandwidth). A unit is predicted by
# the model of its nearest anchor. k is chosen by held-out log-loss over whole
# units; lambda is the global cross-validated value unless supplied.

gw_bisquare <- function(d, h) ifelse(d < h, (1 - (d / h)^2)^2, 0)

gw_fit_anchors <- function(z, xraw, y, rows, ucoords, anchor_idx, k, alpha, lambda, classes,
                           type_multinomial, pre, control) {
  m <- nrow(ucoords)
  lapply(anchor_idx, function(a) {
    d <- sqrt(rowSums(sweep(ucoords, 2, ucoords[a, ], "-")^2))
    ord <- order(d)
    kk <- min(k, m)
    repeat {
      near <- ord[seq_len(kk)]
      h <- max(d[near]) * 1.0001 + 1e-12
      wu <- gw_bisquare(d[near], h)
      ii <- unlist(rows[near], use.names = FALSE)
      wi <- rep(wu, lengths(rows[near]))
      keep <- wi > 0
      if (all(table(factor(y[ii][keep], levels = classes)) > 0L) || kk >= m) break
      kk <- min(2L * kk, m)
    }
    ii <- ii[keep]
    wi <- wi[keep]
    raw <- tryCatch(fit_fixed_lambda_glmnet(z[ii, , drop = FALSE], factor(y[ii], levels = classes),
                                            alpha, lambda, classes, FALSE, type_multinomial,
                                            weights = wi / mean(wi)),
                    error = function(e) NULL)
    if (is.null(raw)) return(NULL)
    eng <- engine_from_glmnet(raw, z, classes, alpha, lambda, control)
    eng$standardize <- FALSE
    eng$pre <- pre
    eng$bandwidth_units <- kk
    # Discard numerically failed local fits (non-finite probabilities).
    p <- tryCatch(engine_predict(eng, xraw[ii, , drop = FALSE]), error = function(e) NULL)
    if (is.null(p) || any(!is.finite(p))) return(NULL)
    eng
  })
}

#' Geographically weighted multinomial elastic net
#'
#' A local-regression comparator for SCMR: smoothly varying coefficients
#' instead of spatial regimes.
#' @inheritParams fit_scmr
#' @param k_grid Candidate adaptive bandwidths (number of nearest units).
#' @param alpha Elastic-net mixing parameter.
#' @param lambda Fixed glmnet lambda; by default the global `cv.glmnet`
#'   `lambda.min` for `alpha`.
#' @param anchors Maximum number of anchor units with a local fit.
#' @param val_prop Share of units held out to choose the bandwidth.
#' @param type_multinomial Passed to glmnet.
#' @return An object of class `scmr_gw`.
#' @export
fit_gw_multinom_en <- function(x, y, unit_id, coords, k_grid = c(25, 50, 100, 200),
                               alpha = 0.5, lambda = NULL, anchors = 300, val_prop = 0.2,
                               type_multinomial = "ungrouped", seed = 123) {
  data <- prepare_scmr_data(x, y, unit_id, coords)
  classes <- levels(data$y)
  pre <- make_pre_transform(data$x)
  z <- apply_pre_transform(data$x, pre)
  control <- scmr_control(standardize = FALSE, type_multinomial = type_multinomial)
  units <- sort(unique(data$unit_id))
  ucoords <- get_unit_coords(data$unit_id, data$coords)
  rows <- split(seq_len(nrow(z)), factor(data$unit_id, levels = units))
  with_scmr_seed(seed, {
    if (is.null(lambda)) {
      cv <- glmnet::cv.glmnet(glmnet_design(z), data$y, family = "multinomial", alpha = alpha,
                              foldid = make_stratified_foldid(data$y, 5, scmr_seed(seed, 11)),
                              standardize = FALSE, type.multinomial = type_multinomial)
      lambda <- cv$lambda.min
    }
    m <- length(units)
    pick_anchors <- function(pool) if (length(pool) <= anchors) pool else sort(sample(pool, anchors))
    # Bandwidth by held-out log-loss over whole units.
    val_u <- sort(sample.int(m, max(1L, floor(val_prop * m))))
    fit_u <- setdiff(seq_len(m), val_u)
    trace <- data.frame(k = k_grid, LogLoss = NA_real_)
    if (length(k_grid) > 1L) {
      sub_rows <- rows[fit_u]
      sub_anchor <- pick_anchors(seq_along(fit_u))
      for (j in seq_along(k_grid)) {
        engines <- gw_fit_anchors(z, data$x, data$y, sub_rows, ucoords[fit_u, , drop = FALSE], sub_anchor,
                                  k_grid[j], alpha, lambda, classes, type_multinomial, pre, control)
        ok <- !vapply(engines, is.null, logical(1))
        if (!any(ok)) next
        ac <- ucoords[fit_u[sub_anchor[ok]], , drop = FALSE]
        nearest <- apply(ucoords[val_u, , drop = FALSE], 1L, function(s) which.min(colSums((t(ac) - s)^2)))
        ll <- 0
        nn <- 0
        for (v in seq_along(val_u)) {
          ii <- rows[[val_u[v]]]
          p <- engine_predict(engines[ok][[nearest[v]]], data$x[ii, , drop = FALSE])
          ll <- ll + sum(true_class_logp(p, data$y[ii], classes))
          nn <- nn + length(ii)
        }
        trace$LogLoss[j] <- -ll / nn
      }
      k_best <- k_grid[which.min(replace(trace$LogLoss, is.na(trace$LogLoss), Inf))]
    } else {
      k_best <- k_grid
    }
    anchor_idx <- pick_anchors(seq_len(m))
    engines <- gw_fit_anchors(z, data$x, data$y, rows, ucoords, anchor_idx, k_best, alpha, lambda,
                              classes, type_multinomial, pre, control)
  })
  ok <- !vapply(engines, is.null, logical(1))
  if (!any(ok)) stop("All local fits failed; increase k_grid or lambda.", call. = FALSE)
  structure(list(engines = engines[ok], anchor_coords = ucoords[anchor_idx[ok], , drop = FALSE],
                 anchor_units = units[anchor_idx[ok]], k = k_best, k_trace = trace,
                 alpha = alpha, lambda = lambda, class_levels = classes,
                 x_colnames = colnames(data$x)), class = "scmr_gw")
}

#' @rdname fit_gw_multinom_en
#' @param object A fitted `scmr_gw` object.
#' @param newx Predictor matrix.
#' @param new_coords Two-column coordinates, one row per prediction row.
#' @param type `"prob"` or `"class"`.
#' @param ... Unused.
#' @export
predict.scmr_gw <- function(object, newx, new_coords, type = c("prob", "class"), ...) {
  type <- match.arg(type)
  newx <- check_newx(newx, object$x_colnames)
  new_coords <- as.matrix(new_coords)
  if (nrow(new_coords) != nrow(newx) || ncol(new_coords) != 2L) stop("new_coords needs one row per prediction row.", call. = FALSE)
  key <- paste(new_coords[, 1], new_coords[, 2])
  uk <- unique(key)
  uc <- new_coords[match(uk, key), , drop = FALSE]
  nearest <- apply(uc, 1L, function(s) which.min(colSums((t(object$anchor_coords) - s)^2)))
  near_row <- nearest[match(key, uk)]
  prob <- matrix(0, nrow(newx), length(object$class_levels), dimnames = list(NULL, object$class_levels))
  for (a in unique(near_row)) {
    ii <- which(near_row == a)
    prob[ii, ] <- engine_predict(object$engines[[a]], newx[ii, , drop = FALSE])
  }
  if (type == "prob") return(prob)
  factor(object$class_levels[max.col(prob, ties.method = "first")], levels = object$class_levels)
}

#' Block cross-validation folds
#'
#' Assigns whole blocks (for example KSA segments) to folds, balancing the
#' number of blocks per fold within each stratum (for example regency).
#' @param block_id Block identifier per row.
#' @param strata Optional stratum per row (constant within a block).
#' @param nfolds Number of folds.
#' @param seed Random seed; the caller's random state is restored.
#' @return Integer fold per row.
#' @export
scmr_block_folds <- function(block_id, strata = NULL, nfolds = 5, seed = 123) {
  block_id <- as.character(block_id)
  if (anyNA(block_id)) stop("block_id must be nonmissing.", call. = FALSE)
  blocks <- unique(block_id)
  st <- if (is.null(strata)) rep("all", length(blocks)) else {
    s <- tapply(as.character(strata), block_id, function(z) {
      if (length(unique(z)) != 1L) stop("strata must be constant within a block.", call. = FALSE)
      z[1]
    })
    s[blocks]
  }
  fold_of_block <- with_scmr_seed(seed, {
    f <- integer(length(blocks))
    offset <- 0L
    for (s in unique(st)) {
      b <- which(st == s)
      b <- b[sample.int(length(b))]
      f[b] <- ((seq_along(b) - 1L + offset) %% nfolds) + 1L
      offset <- offset + length(b)
    }
    f
  })
  fold_of_block[match(block_id, blocks)]
}
