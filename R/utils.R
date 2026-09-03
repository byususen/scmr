# Internal utility functions ---------------------------------------------------

safe_div <- function(num, den) {
  ifelse(is.finite(den) & den > 0, num / den, NA_real_)
}

softmax_rows <- function(q) {
  q <- as.matrix(q)
  mx <- apply(q, 1L, max)
  z <- exp(q - mx)
  z / rowSums(z)
}

normalize_prob_matrix <- function(out, y_levels, n_test = NULL) {
  y_levels <- as.character(y_levels)
  if (!length(y_levels) || anyNA(y_levels) || anyDuplicated(y_levels)) {
    stop("Class labels must be nonmissing and unique.", call. = FALSE)
  }
  if (is.null(dim(out)) && length(out) == length(y_levels) &&
      (is.null(n_test) || n_test == 1L)) {
    out <- matrix(out, nrow = 1L, dimnames = list(NULL, names(out)))
  }
  out <- as.matrix(out)
  # The public convention is observations by classes. A square matrix already
  # satisfies that convention and must never be transposed by shape alone.
  if (!is.null(n_test) && nrow(out) != n_test &&
      nrow(out) == length(y_levels) && ncol(out) == n_test) {
    out <- t(out)
  }
  if (!is.null(n_test) && nrow(out) != n_test) {
    stop("Probability rows do not match the number of observations.", call. = FALSE)
  }
  if (!is.numeric(out) || any(!is.finite(out)) || any(out < 0)) {
    stop("Probabilities must be finite, nonnegative numbers.", call. = FALSE)
  }
  full_out <- matrix(0, nrow = nrow(out), ncol = length(y_levels))
  colnames(full_out) <- y_levels
  if (!is.null(colnames(out))) {
    if (anyDuplicated(colnames(out)) || any(!colnames(out) %in% y_levels)) {
      stop("Probability columns contain unknown or duplicate class labels.", call. = FALSE)
    }
    common <- intersect(colnames(out), y_levels)
    full_out[, common] <- out[, common, drop = FALSE]
  } else if (ncol(out) == length(y_levels)) {
    full_out <- out
    colnames(full_out) <- y_levels
  } else {
    stop("Probability matrix has incompatible number of columns.", call. = FALSE)
  }
  out <- full_out
  rs <- rowSums(out)
  if (any(rs <= 0)) {
    stop("Each probability row must have positive total mass.", call. = FALSE)
  }
  out <- out / rs
  rownames(out) <- NULL
  out
}

make_one_hot <- function(y, classes) {
  y <- factor(y, levels = classes)
  mat <- matrix(0, nrow = length(y), ncol = length(classes))
  colnames(mat) <- classes
  idx <- match(as.character(y), classes)
  ok <- !is.na(idx)
  mat[cbind(which(ok), idx[ok])] <- 1
  mat
}

true_class_logp <- function(prob, y, classes, eps = 1e-15) {
  y <- factor(y, levels = classes)
  if (anyNA(y)) stop("Observed labels must belong to the model classes.", call. = FALSE)
  prob <- normalize_prob_matrix(prob, classes, length(y))
  idx <- cbind(seq_along(y), match(as.character(y), colnames(prob)))
  log(pmax(prob[idx], eps))
}

multiclass_logloss <- function(prob, y_true, classes = NULL) {
  if (is.null(classes)) classes <- colnames(prob)
  if (is.null(classes) && is.factor(y_true)) classes <- levels(y_true)
  if (is.null(classes)) {
    stop("Supply model class labels or named probability columns.", call. = FALSE)
  }
  -mean(true_class_logp(prob, y_true, classes))
}

multiclass_brier <- function(prob, y_true, classes) {
  prob <- normalize_prob_matrix(prob, classes, length(y_true))
  oh <- make_one_hot(y_true, classes)
  mean(rowSums((prob - oh)^2), na.rm = TRUE)
}

manual_kappa <- function(y_true, y_pred, classes) {
  y_true <- factor(y_true, levels = classes)
  y_pred <- factor(y_pred, levels = classes)
  tab <- table(y_pred, y_true)
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / (n^2)
  if (!is.finite(pe) || abs(1 - pe) < 1e-12) return(NA_real_)
  (po - pe) / (1 - pe)
}

#' Per-class classification metrics
#'
#' @param y_true True class labels.
#' @param y_pred Predicted class labels.
#' @param classes Character vector of class levels.
#' @return A data frame with one row per class.
#' @export
calc_class_metrics <- function(y_true, y_pred, classes) {
  if (length(y_true) != length(y_pred)) {
    stop("True and predicted labels must have equal lengths.", call. = FALSE)
  }
  y_true <- factor(y_true, levels = classes)
  y_pred <- factor(y_pred, levels = classes)
  if (anyNA(y_true) || anyNA(y_pred)) {
    stop("Labels must be nonmissing and belong to classes.", call. = FALSE)
  }
  out <- lapply(classes, function(cl) {
    tp <- sum(y_true == cl & y_pred == cl, na.rm = TRUE)
    fp <- sum(y_true != cl & y_pred == cl, na.rm = TRUE)
    fn <- sum(y_true == cl & y_pred != cl, na.rm = TRUE)
    tn <- sum(y_true != cl & y_pred != cl, na.rm = TRUE)
    precision <- safe_div(tp, tp + fp)
    recall <- safe_div(tp, tp + fn)
    specificity <- safe_div(tn, tn + fp)
    f1 <- safe_div(2 * tp, 2 * tp + fp + fn)
    data.frame(
      Class = cl,
      Support = tp + fn,
      PredictedSupport = tp + fp,
      TP = tp, FP = fp, FN = fn, TN = tn,
      Precision = precision,
      Recall = recall,
      Specificity = specificity,
      F1 = f1,
      BalancedAccuracy = mean(c(recall, specificity), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' Overall multiclass classification metrics
#'
#' @param y_true True labels.
#' @param y_pred Predicted labels.
#' @param prob Matrix of predicted probabilities.
#' @param classes Character vector of class levels.
#' @return One-row data frame of overall metrics.
#' @export
calc_overall_metrics <- function(y_true, y_pred, prob, classes) {
  y_true <- factor(y_true, levels = classes)
  y_pred <- factor(y_pred, levels = classes)
  prob <- normalize_prob_matrix(prob, classes, length(y_true))
  cls <- calc_class_metrics(y_true, y_pred, classes)
  support_weight <- cls$Support / sum(cls$Support)
  data.frame(
    N = length(y_true),
    Accuracy = mean(y_true == y_pred, na.rm = TRUE),
    Kappa = manual_kappa(y_true, y_pred, classes),
    LogLoss = multiclass_logloss(prob, y_true, classes),
    BrierScore = multiclass_brier(prob, y_true, classes),
    MacroPrecision = mean(cls$Precision, na.rm = TRUE),
    MacroRecall = mean(cls$Recall, na.rm = TRUE),
    MacroF1 = mean(cls$F1, na.rm = TRUE),
    MacroBalancedAccuracy = mean(cls$BalancedAccuracy, na.rm = TRUE),
    WeightedPrecision = sum(cls$Precision * support_weight, na.rm = TRUE),
    WeightedRecall = sum(cls$Recall * support_weight, na.rm = TRUE),
    WeightedF1 = sum(cls$F1 * support_weight, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
