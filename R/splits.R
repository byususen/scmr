#' Stratified train/test split
#'
#' @param y Class labels.
#' @param region Optional region labels. If supplied, strata are region by class.
#' @param test_prop Test-set proportion.
#' @param seed Random seed.
#' @return Integer vector of test-row indices.
#' @export
make_stratified_split <- function(y, region = NULL, test_prop = 0.30, seed = 123) {
  set.seed(seed)
  idx_all <- seq_along(y)
  if (is.null(region)) {
    strata <- as.character(y)
  } else {
    strata <- paste(as.character(region), as.character(y), sep = "_")
  }
  split_idx <- split(idx_all, strata)
  test_idx <- unlist(lapply(split_idx, function(idx) {
    if (length(idx) <= 1) return(integer(0))
    sample(idx, size = max(1, ceiling(test_prop * length(idx))))
  }))
  sort(unique(test_idx))
}
