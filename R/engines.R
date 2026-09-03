# Local numerical engines. No study metadata or file output is used here.

prepare_scmr_data <- function(x, y, unit_id = NULL, coords = NULL) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !nrow(x) || !ncol(x) || any(!is.finite(x))) {
    stop("x must be a nonempty finite numeric matrix.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  if (anyNA(colnames(x)) || anyDuplicated(colnames(x)) || any(!nzchar(colnames(x)))) {
    stop("Predictor names must be nonempty and unique.", call. = FALSE)
  }
  if (length(y) != nrow(x) || anyNA(y)) stop("y must match x and have no missing labels.", call. = FALSE)
  y <- if (is.factor(y)) y else factor(y)
  if (nlevels(y) < 2 || any(table(y) == 0L)) {
    stop("Training data must contain every declared class and at least two classes.", call. = FALSE)
  }
  if (is.null(unit_id)) unit_id <- paste0("row", seq_len(nrow(x)))
  unit_id <- as.character(unit_id)
  if (length(unit_id) != nrow(x) || anyNA(unit_id) || any(!nzchar(unit_id))) {
    stop("unit_id must contain one nonmissing identifier per row.", call. = FALSE)
  }
  if (!is.null(coords)) {
    coords <- as.matrix(coords)
    if (!is.numeric(coords) || nrow(coords) != nrow(x) || ncol(coords) != 2L || any(!is.finite(coords))) {
      stop("coords must be a finite numeric matrix with one row per observation and two columns.", call. = FALSE)
    }
    get_unit_coords(unit_id, coords)
  }
  list(x = x, y = y, unit_id = unit_id, coords = coords)
}

check_newx <- function(newx, feature_names) {
  if (is.null(dim(newx))) newx <- matrix(newx, nrow = 1L, dimnames = list(NULL, names(newx)))
  newx <- as.matrix(newx)
  if (!is.numeric(newx) || ncol(newx) != length(feature_names) || any(!is.finite(newx))) {
    stop("newx must be finite and have the fitted predictor columns.", call. = FALSE)
  }
  if (is.null(colnames(newx))) colnames(newx) <- feature_names
  if (anyDuplicated(colnames(newx)) || !setequal(colnames(newx), feature_names)) {
    stop("newx predictor names do not match the model.", call. = FALSE)
  }
  newx[, feature_names, drop = FALSE]
}

engine_from_glmnet <- function(fit, x, classes, alpha, lambda, control, warnings = character()) {
  raw <- if (inherits(fit, "cv.glmnet")) fit$glmnet.fit else fit
  list(kind = "glmnet", fit = fit, classes = classes, features = colnames(x),
       alpha = alpha, lambda = lambda, standardize = control$standardize,
       type_multinomial = control$type_multinomial,
       converged = is.null(raw$jerr) || identical(as.integer(raw$jerr), 0L),
       warnings = warnings)
}

fit_local_engine <- function(x, y, penalty, alpha, lambda, control) {
  classes <- levels(y)
  if (any(table(factor(y, levels = classes)) == 0L)) {
    stop("Every fitted cluster must contain every response class.", call. = FALSE)
  }
  warnings <- character()
  capture_warning <- function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  if (penalty != "none") {
    raw <- withCallingHandlers(
      fit_fixed_lambda_glmnet(x, y, alpha, lambda, classes,
                             control$standardize, control$type_multinomial,
                             maxit = control$glmnet_maxit), warning = capture_warning)
    return(engine_from_glmnet(raw, x, classes, alpha, lambda, control, warnings))
  }
  center <- if (control$standardize) colMeans(x) else rep(0, ncol(x))
  scale <- if (control$standardize) apply(x, 2, stats::sd) else rep(1, ncol(x))
  scale[!is.finite(scale) | scale == 0] <- 1
  z <- sweep(sweep(x, 2, center, "-"), 2, scale, "/")
  # Stable syntactic engine names permit arbitrary user-facing predictor names.
  colnames(z) <- paste0("V", seq_len(ncol(x)))
  dat <- as.data.frame(z)
  dat$.response <- y
  raw <- withCallingHandlers(nnet::multinom(
    .response ~ ., data = dat, decay = 0, trace = FALSE,
    maxit = control$multinom_maxit, MaxNWts = control$multinom_maxnwts), warning = capture_warning)
  list(kind = "multinom", fit = raw, classes = classes, features = colnames(x),
       center = center, scale = scale, alpha = NA_real_, lambda = NA_real_,
       converged = identical(as.integer(raw$convergence), 0L), warnings = warnings)
}

engine_predict <- function(engine, newx) {
  newx <- check_newx(newx, engine$features)
  if (!nrow(newx)) return(matrix(numeric(), 0L, length(engine$classes),
                                dimnames = list(NULL, engine$classes)))
  if (engine$kind == "glmnet") {
    p <- predict_multinom_prob(engine$fit, newx, engine$lambda)
  } else {
    z <- sweep(sweep(newx, 2, engine$center, "-"), 2, engine$scale, "/")
    colnames(z) <- paste0("V", seq_len(ncol(newx)))
    p <- stats::predict(engine$fit, as.data.frame(z), type = "probs")
    if (length(engine$classes) == 2L && is.null(dim(p))) {
      p <- cbind(1 - as.numeric(p), as.numeric(p))
      colnames(p) <- engine$classes
    } else if (is.null(dim(p))) {
      p <- matrix(p, nrow = 1L, dimnames = list(NULL, names(p)))
    }
  }
  normalize_prob_matrix(p, engine$classes, nrow(newx))
}

engine_coef <- function(engine, parameterization = c("native", "sum_to_zero")) {
  parameterization <- match.arg(parameterization)
  out <- matrix(0, length(engine$classes), length(engine$features) + 1L,
                dimnames = list(engine$classes, c("(Intercept)", engine$features)))
  if (engine$kind == "glmnet") {
    cf <- stats::coef(engine$fit, s = engine$lambda)
    for (cl in engine$classes) out[cl, ] <- as.numeric(cf[[cl]][colnames(out), 1])
  } else {
    cf <- stats::coef(engine$fit)
    if (is.null(dim(cf))) cf <- matrix(cf, nrow = 1L,
                                      dimnames = list(engine$classes[2], names(cf)))
    out[rownames(cf), ] <- cf
    slopes <- out[, -1, drop = FALSE]
    out[, 1] <- out[, 1] - as.numeric(slopes %*% (engine$center / engine$scale))
    out[, -1] <- sweep(slopes, 2, engine$scale, "/")
  }
  if (parameterization == "sum_to_zero") out <- sweep(out, 2, colMeans(out), "-")
  out
}

with_scmr_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  if (length(seed) != 1L || !is.finite(seed) || seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) {
    stop("seed must be NULL or a nonnegative integer within R's range.", call. = FALSE)
  }
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (existed) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)
  force(expr)
}

scmr_seed <- function(seed, offset) {
  if (is.null(seed)) return(NULL)
  as.integer((as.double(seed) + offset) %% .Machine$integer.max)
}
