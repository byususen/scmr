build_unit_weights <- function(coords, k = 5, weight_type = c("binary", "exp"),
                               bandwidth = NULL, symmetrize = TRUE,
                               row_standardize = FALSE) {
  weight_type <- match.arg(weight_type)
  coords <- as.matrix(coords)
  n <- nrow(coords)
  if (n < 2) stop("Need at least two units to build spatial weights.", call. = FALSE)
  if (k >= n) k <- n - 1L
  if (k < 1) stop("Not enough units to build spatial weights.", call. = FALSE)

  dmat <- as.matrix(stats::dist(coords))
  diag(dmat) <- Inf
  nn_index <- t(apply(dmat, 1L, order))[, seq_len(k), drop = FALSE]
  w <- matrix(0, n, n)

  if (weight_type == "binary") {
    for (i in seq_len(n)) w[i, nn_index[i, ]] <- 1
    bandwidth_used <- NA_real_
  } else {
    if (is.null(bandwidth)) {
      kth_dist <- dmat[cbind(seq_len(n), nn_index[, k])]
      bandwidth <- stats::median(kth_dist, na.rm = TRUE)
    }
    bandwidth_used <- bandwidth
    for (i in seq_len(n)) {
      jj <- nn_index[i, ]
      w[i, jj] <- exp(-(dmat[i, jj]^2) / (bandwidth^2))
    }
  }

  if (symmetrize) w <- pmax(w, t(w))
  if (row_standardize) {
    rs <- rowSums(w)
    rs[rs == 0] <- 1
    w <- w / rs
  }

  list(W = Matrix::Matrix(w, sparse = TRUE), bandwidth_used = bandwidth_used)
}

get_unit_coords <- function(unit_id, coords) {
  coord_df <- unique(data.frame(
    unit = as.character(unit_id),
    xcoord = coords[, 1],
    ycoord = coords[, 2],
    stringsAsFactors = FALSE
  ))
  unit_levels <- levels(factor(unit_id))
  coord_df <- coord_df[match(unit_levels, coord_df$unit), , drop = FALSE]
  as.matrix(coord_df[, c("xcoord", "ycoord"), drop = FALSE])
}

resolve_min_class <- function(min_class, class_levels) {
  k <- length(class_levels)
  if (length(min_class) == 1) {
    out <- rep(as.integer(min_class), k)
    names(out) <- class_levels
    return(out)
  }
  out <- as.integer(min_class)
  names(out) <- class_levels
  out
}

precompute_unit_contributions <- function(unit_index, y, class_levels) {
  m <- length(unit_index)
  k <- length(class_levels)
  unit_total <- integer(m)
  unit_class <- matrix(0L, nrow = m, ncol = k)
  colnames(unit_class) <- class_levels
  for (u in seq_len(m)) {
    idx <- unit_index[[u]]
    yu <- factor(y[idx], levels = class_levels)
    unit_total[u] <- length(idx)
    unit_class[u, ] <- as.integer(table(yu))
  }
  list(unit_total = unit_total, unit_class = unit_class)
}

compute_cluster_counts <- function(group_unit, G, unit_total, unit_class) {
  k <- ncol(unit_class)
  unit_count <- integer(G)
  total_count <- integer(G)
  class_count <- matrix(0L, nrow = G, ncol = k)
  colnames(class_count) <- colnames(unit_class)
  for (g in seq_len(G)) {
    idx <- which(group_unit == g)
    unit_count[g] <- length(idx)
    if (length(idx) > 0) {
      total_count[g] <- sum(unit_total[idx])
      class_count[g, ] <- colSums(unit_class[idx, , drop = FALSE])
    }
  }
  list(unit_count = unit_count, total_count = total_count, class_count = class_count)
}

is_feasible_partition <- function(group_unit, G, unit_total, unit_class,
                                  min_units, min_class_vec) {
  cnt <- compute_cluster_counts(group_unit, G, unit_total, unit_class)
  if (any(cnt$unit_count < min_units)) return(FALSE)
  for (g in seq_len(G)) {
    if (any(cnt$class_count[g, ] < min_class_vec)) return(FALSE)
  }
  TRUE
}

relabel_clusters <- function(group_unit) {
  as.integer(factor(group_unit, levels = unique(group_unit)))
}

initialize_feasible_groups <- function(unit_coords, G, unit_total, unit_class,
                                       min_units, min_class_vec,
                                       tries = 100, seed = NULL,
                                       method = c("kmeans", "random"),
                                       initial_cluster = NULL) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)
  m <- nrow(unit_coords)

  if (!is.null(initial_cluster)) {
    if (length(initial_cluster) != m) {
      stop("initial_cluster must have one value per unique spatial unit.", call. = FALSE)
    }
    grp <- relabel_clusters(initial_cluster)
    if (length(unique(grp)) != G) {
      stop("initial_cluster must contain exactly G distinct clusters.", call. = FALSE)
    }
    if (!is_feasible_partition(grp, G, unit_total, unit_class, min_units, min_class_vec)) {
      stop("initial_cluster is not feasible under min_units/min_per_class.", call. = FALSE)
    }
    return(grp)
  }

  if (G == 1) {
    grp <- rep(1L, m)
    if (!is_feasible_partition(grp, G, unit_total, unit_class, min_units, min_class_vec)) {
      stop("G = 1 is infeasible under the supplied constraints.", call. = FALSE)
    }
    return(grp)
  }

  if (G >= m) stop("G must be smaller than the number of spatial units.", call. = FALSE)
  if (m < G * min_units) {
    stop("Not enough unique spatial units for G = ", G, " with min_units = ", min_units, call. = FALSE)
  }

  coords_scaled <- scale(unit_coords)

  for (tt in seq_len(tries)) {
    if (method == "kmeans") {
      grp <- stats::kmeans(coords_scaled, centers = G, nstart = 1)$cluster
    } else {
      base <- rep(seq_len(G), each = min_units)
      rest <- sample(seq_len(G), size = m - length(base), replace = TRUE)
      grp <- sample(c(base, rest), size = m)
    }
    grp <- relabel_clusters(grp)
    if (length(unique(grp)) == G &&
        is_feasible_partition(grp, G, unit_total, unit_class, min_units, min_class_vec)) {
      return(grp)
    }
  }

  stop("Could not find feasible initialization under unique-unit and class-count constraints.", call. = FALSE)
}

move_is_feasible <- function(u, g_new, group_unit, unit_count, class_count,
                             unit_class, min_units, min_class_vec) {
  g_old <- group_unit[u]
  if (g_new == g_old) return(TRUE)
  unit_old_after <- unit_count[g_old] - 1L
  class_old_after <- class_count[g_old, ] - unit_class[u, ]
  if (unit_old_after < min_units) return(FALSE)
  if (any(class_old_after < min_class_vec)) return(FALSE)
  TRUE
}
