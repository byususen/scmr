load_article_helpers <- function() {
  env <- new.env(parent = globalenv())
  for (nm in c("config", "io", "evaluate", "summarize", "run")) {
    sys.source(system.file("simulation", paste0(nm, ".R"), package = "scmr", mustWork = TRUE), env)
  }
  env
}

test_that("the study runs through the installed API and resumes complete bundles", {
  a <- load_article_helpers()
  path <- tempfile("study-")
  cfg <- a$article_config("smoke", path, scenarios = "global", n_obs = 75L, n_new = 3L,
    G = 1L, n_active = 1L, n_inactive = 1L, nfolds = 3L,
    fit_control = scmr_control(min_units = 2, min_per_class = 2,
      alpha_grid = .5, type_multinomial = "ungrouped"))
  expect_message(first <- a$run_article_simulation(cfg), "Results:")
  expect_equal(nrow(first$results), 4L)
  expect_equal(nrow(first$errors), 0L)
  expect_true(all(first$results$CBBIC_WeightedN_Check == first$results$n_rows))
  paths <- list.files(file.path(first$output_dir, "runs"), full.names = TRUE)
  before <- tools::md5sum(paths)
  expect_message(second <- a$run_article_simulation(cfg), "Already complete")
  expect_identical(tools::md5sum(paths), before)
  expect_equal(second$results, first$results)
  expect_true(file.exists(file.path(first$output_dir, "parameters_region.csv")))
  expect_equal(nrow(read.csv(file.path(first$output_dir, "best_G_CB_BIC.csv"))), 1)
  expect_error(a$article_open(within(cfg, resume <- FALSE)), "already exist")
  different <- a$article_open(within(cfg, fit_control$phi <- 2))
  expect_false(identical(different, first$output_dir))
})

test_that("missing G candidates are visible and cannot select a best G", {
  a <- load_article_helpers()
  out <- tempfile("reports-"); dir.create(out)
  dir.create(file.path(out, "runs")); dir.create(file.path(out, "errors"))
  cfg <- a$article_config(G = c(1L, 2L))
  expected <- data.frame(DatasetID = "absent-data", RunID = c("a", "b"))
  a$article_summarize(out, expected, cfg)
  coverage <- read.csv(file.path(out, "selection_coverage_CB_BIC.csv"))
  expect_identical(coverage$Complete, FALSE)
  expect_equal(coverage$NEligible, 0)
  expect_equal(coverage$NExpected, 2)
})
