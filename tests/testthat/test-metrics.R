test_that("probabilities preserve row orientation, labels, and batch size", {
  p <- matrix(c(.7,.2,.1, .2,.5,.3, .1,.4,.5), 3, byrow = TRUE,
              dimnames = list(NULL, c("a", "b", "c")))
  norm <- scmr:::normalize_prob_matrix
  expect_equal(norm(p, colnames(p), 3), p)
  expect_equal(norm(p[1, ], colnames(p), 1), p[1, , drop = FALSE])
  expect_equal(norm(p[, c("c", "a", "b")], colnames(p), 3), p)
  expect_equal(norm(t(p[1:2, ]), colnames(p), 2), p[1:2, ])
  expect_error(norm(p, colnames(p), 4), "rows")
  expect_error(norm(matrix(NA_real_, 1, 3), colnames(p), 1), "finite")
})

test_that("absent evaluation classes do not change log-loss", {
  p <- matrix(rep(c(.6,.3,.1), 4), 4, byrow = TRUE,
              dimnames = list(NULL, c("a", "b", "c")))
  y <- factor(rep("a", 4), levels = colnames(p))
  expect_equal(scmr:::multiclass_logloss(p, y), -log(.6))
  expect_equal(scmr:::multiclass_logloss(p, as.character(y)), -log(.6))
  expect_equal(scmr:::multiclass_logloss(p[, c(3,1,2)], y), -log(.6))
})

test_that("missed supported classes contribute zero F1", {
  cls <- c("a", "b", "c")
  m <- calc_class_metrics(cls, c("a", "c", "b"), cls)
  expect_equal(m$F1, c(1, 0, 0))
  p <- matrix(c(.8,.1,.1, .1,.1,.8, .1,.8,.1), 3, byrow = TRUE,
              dimnames = list(NULL, cls))
  overall <- calc_overall_metrics(cls, c("a", "c", "b"), p, cls)
  expect_equal(overall$MacroF1, 1/3)
  expect_equal(overall$WeightedF1, 1/3)
  expect_equal(calc_class_metrics(c("a","b"), c("a","a"), cls)$F1[2], 0)
  expect_true(is.na(calc_class_metrics("a", "a", cls)$F1[3]))
})
