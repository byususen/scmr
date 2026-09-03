article_bind <- function(rows) {
  rows <- Filter(function(z) is.data.frame(z) && nrow(z) > 0L, rows)
  if (!length(rows)) return(data.frame())
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  out <- do.call(rbind, lapply(rows, function(z) {
    for (nm in setdiff(cols, names(z))) z[[nm]] <- NA
    z[, cols, drop = FALSE]
  }))
  rownames(out) <- NULL
  out
}

article_meta <- function(table, meta) {
  if (!nrow(table)) return(data.frame())
  cbind(meta[rep(1L, nrow(table)), , drop = FALSE], table)
}

article_hash <- function(object) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(object, path, version = 2)
  unname(tools::md5sum(path))
}

article_save <- function(object, path) {
  tmp <- tempfile(pattern = "pending-", tmpdir = dirname(path))
  on.exit(unlink(tmp))
  saveRDS(object, tmp, version = 2)
  if (!file.rename(tmp, path)) stop("Could not commit checkpoint: ", path)
  invisible(path)
}

article_csv <- function(table, path) {
  # CSVs are derived reports. Completed RDS bundles remain the source of truth.
  tmp <- tempfile(pattern = "report-", tmpdir = dirname(path))
  on.exit(unlink(tmp))
  utils::write.csv(table, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path) && !file.copy(tmp, path, overwrite = TRUE)) stop("Could not write report: ", path)
}

article_open <- function(cfg) {
  settings <- cfg[setdiff(names(cfg), c("output_dir", "resume"))]
  ns <- asNamespace("scmr")
  package_code <- lapply(sort(ls(ns, all.names = TRUE)), function(nm) {
    obj <- get(nm, ns)
    if (is.function(obj)) list(name = nm, formals = paste(deparse(formals(obj)), collapse = "\n"),
      body = paste(deparse(body(obj)), collapse = "\n")) else NULL
  })
  identity <- list(settings = settings, code = package_code,
    engines = vapply(c("scmr", "glmnet", "nnet", "Matrix", "MASS"), function(p) as.character(utils::packageVersion(p)), character(1)),
    R = R.version.string)
  # Include the study implementation too: stale evaluation code cannot share checkpoints.
  identity$study <- lapply(sort(ls(environment(article_open), pattern = "^article_|^run_article_")), function(nm) {
    obj <- get(nm, environment(article_open))
    if (is.function(obj)) list(name = nm, formals = paste(deparse(formals(obj)), collapse = "\n"),
      body = paste(deparse(body(obj)), collapse = "\n")) else NULL
  })
  key <- article_hash(identity)
  out <- file.path(cfg$output_dir, paste0(cfg$profile, "-", substr(key, 1, 12)))
  if (!cfg$resume && dir.exists(out)) stop("Results already exist. Use resume=TRUE or a different output_dir.")
  for (d in c(out, file.path(out, c("runs", "errors", "data")))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  manifest <- file.path(out, "manifest.rds")
  if (!file.exists(manifest)) {
    article_save(list(config = cfg, fingerprint = key, versions = identity$engines,
      created_utc = format(Sys.time(), tz = "UTC"), session = utils::sessionInfo(),
      source_sha256 = "921bd7679b8e659824df38dcae3a57e10abe44f2776222375b5ca0b59acd211d"), manifest)
    writeLines(capture.output(dput(cfg)), file.path(out, "config.R"))
    writeLines(capture.output(utils::sessionInfo()), file.path(out, "session-info.txt"))
  } else if (!identical(readRDS(manifest)$fingerprint, key)) stop("Checkpoint identity mismatch.")
  normalizePath(out)
}
