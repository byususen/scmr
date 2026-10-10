multinomial_hessian <- function(z, prob) {
  q <- ncol(z)
  k <- ncol(prob)
  h <- matrix(0, k * q, k * q)
  for (a in seq_len(k)) for (b in seq_len(k)) {
    ia <- (a - 1L) * q + seq_len(q)
    ib <- (b - 1L) * q + seq_len(q)
    weight <- prob[, a] * (as.numeric(a == b) - prob[, b])
    h[ia, ib] <- crossprod(z, z * weight)
  }
  (h + t(h)) / 2
}

engine_effective_df <- function(engine, x, tol = 1e-8) {
  native <- engine_coef(engine)
  slopes <- native[, -1, drop = FALSE]
  selected <- abs(slopes) > tol
  k <- nrow(native)
  p <- ncol(x)
  n <- nrow(x)
  nominal <- (k - 1L) * (p + 1L)
  if (engine$kind == "multinom") {
    return(list(df = (k - 1L) * qr(cbind(1, x))$rank,
                df_nominal_active = (k - 1L) * qr(cbind(1, x))$rank,
                nominal = nominal, active_predictors = p,
                active_coefficients = k * p, method = "unpenalized_design_rank"))
  }
  grouped <- engine$type_multinomial == "grouped"
  active_predictors <- colSums(selected) > 0
  if (grouped) selected <- matrix(rep(active_predictors, each = k), k, p)
  scale <- rep(1, p)
  z <- x
  if (!is.null(engine$pre)) {
    # Fitted on globally standardized predictors with standardize = FALSE.
    z <- apply_pre_transform(x, engine$pre)
    scale <- engine$pre$scale
  } else if (engine$standardize) {
    z <- sweep(x, 2, colMeans(x), "-")
    scale <- sqrt(colMeans(z^2))
    scale[scale <= 1e-12] <- 1
    z <- sweep(z, 2, scale, "/")
  }
  z <- cbind(1, z)
  q <- ncol(z)
  h <- multinomial_hessian(z, engine_predict(engine, x))
  active <- as.vector(t(cbind(TRUE, selected)))
  # glmnet uses mean negative log-likelihood plus an elastic-net penalty.
  pen <- diag(rep(c(0, rep(n * engine$lambda * (1 - engine$alpha), p)), k))
  if (grouped && engine$alpha > 0) {
    # Nonzero grouped L1 terms have vector-norm curvature. This term is absent
    # from the scalar L1 active-set approximation for ungrouped multinomial EN.
    beta_std <- sweep(slopes, 2, scale, "*")
    for (j in which(active_predictors)) {
      b <- beta_std[, j]
      norm_b <- sqrt(sum(b^2))
      if (norm_b > tol) {
        idx <- (seq_len(k) - 1L) * q + j + 1L
        pen[idx, idx] <- pen[idx, idx] + n * engine$lambda * engine$alpha /
          norm_b * (diag(k) - tcrossprod(b / norm_b))
      }
    }
  }
  h <- h[active, active, drop = FALSE]
  pen <- pen[active, active, drop = FALSE]
  rank_h <- as.numeric(Matrix::rankMatrix(h, tol = 1e-8))
  df <- sum(diag(MASS::ginv(h + pen) %*% h))
  if (!is.finite(df)) stop("Nonfinite effective degrees of freedom.", call. = FALSE)
  list(df = max(0, min(df, rank_h)), df_nominal_active = rank_h,
       nominal = nominal, active_predictors = sum(active_predictors),
       active_coefficients = sum(selected),
       method = paste0("multinomial_active_hessian_", if (grouped) "grouped" else "ungrouped"))
}
