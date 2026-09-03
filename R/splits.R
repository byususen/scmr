#' Stratified train/test split
#'
#' @param y Nonmissing class labels.
#' @param region Optional region labels, one per row. Strata are region by class.
#' @param test_prop Test-set proportion in (0, 1).
#' @param seed Random seed. The caller's random state is restored.
#' @return Sorted integer test-row indices; at least one row per stratum remains
#'   in training. Singleton strata contribute no test rows.
#' @export
make_stratified_split <- function(y, region = NULL, test_prop = 0.30, seed = 123) {
  if (anyNA(y) || !length(y)) stop("y must contain nonmissing class labels.", call. = FALSE)
  if (length(test_prop) != 1L || !is.finite(test_prop) || test_prop <= 0 || test_prop >= 1) {
    stop("test_prop must be in (0, 1).", call. = FALSE)
  }
  if (!is.null(region) && (length(region) != length(y) || anyNA(region))) {
    stop("region must contain one nonmissing label per row.", call. = FALSE)
  }
  strata <- if (is.null(region)) factor(y) else interaction(factor(region), factor(y), drop = TRUE)
  idx <- split(seq_along(y), strata)
  with_scmr_seed(seed, sort(as.integer(unlist(lapply(idx, function(rows) {
    if (length(rows) <= 1L) return(integer())
    size <- min(length(rows) - 1L, max(1L, ceiling(test_prop * length(rows))))
    rows[sample.int(length(rows), size)]
  }), use.names = FALSE))))
}
