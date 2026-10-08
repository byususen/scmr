article_summary <- function(table, group_columns) {
  if (!nrow(table)) return(data.frame())
  group_columns <- intersect(group_columns, names(table))
  numeric_columns <- setdiff(names(table)[vapply(table, is.numeric, logical(1))],
    c(group_columns, "Repeat", "Seed"))
  groups <- split(seq_len(nrow(table)), interaction(table[, group_columns, drop = FALSE], drop = TRUE, lex.order = TRUE))
  article_bind(lapply(groups, function(ii) {
    z <- table[ii, , drop = FALSE]
    meta <- z[1, group_columns, drop = FALSE]
    article_bind(lapply(numeric_columns, function(nm) {
      values <- z[[nm]][is.finite(z[[nm]])]
      n <- length(values)
      avg <- if (n) mean(values) else NA_real_
      se <- if (n > 1L) stats::sd(values) / sqrt(n) else NA_real_
      half <- if (n > 1L) stats::qt(.975, df = n - 1L) * se else NA_real_
      article_meta(data.frame(Metric = nm, NFinite = n, Mean = avg, MCSE = se,
        Lower95 = avg - half, Upper95 = avg + half), meta)
    }))
  }))
}

article_summarize <- function(out, expected, cfg) {
  paths <- list.files(file.path(out, "runs"), pattern = "\\.rds$", full.names = TRUE)
  bundles <- lapply(paths, readRDS)
  result <- article_bind(lapply(bundles, `[[`, "result"))
  article_csv(result, file.path(out, "result_by_k.csv"))
  error_paths <- list.files(file.path(out, "errors"), pattern = "\\.rds$", full.names = TRUE)
  errors <- article_bind(lapply(error_paths, readRDS))
  article_csv(errors, file.path(out, "errors.csv"))
  expected$Status <- "missing"
  if (nrow(result)) {
    at <- match(expected$RunID, result$RunID); ok <- !is.na(at)
    expected$Status[ok] <- ifelse(result$Converged[at[ok]], "converged", "not_converged")
  }
  if (nrow(errors)) expected$Status[expected$RunID %in% errors$RunID & expected$Status == "missing"] <- "error"
  article_csv(expected, file.path(out, "run_coverage.csv"))
  grouping <- c("Scenario", "Eta", "ClassBalance", "NActive", "NInactive", "Heterogeneity", "P", "Model", "Mode", "Penalty", "G")
  # Retain all fitted results; only eligible fits enter Monte Carlo summaries and G selection.
  valid <- if (nrow(result)) result[result$Eligible, , drop = FALSE] else result
  article_csv(article_summary(valid, grouping), file.path(out, "result_summary.csv"))
  table_names <- unique(unlist(lapply(bundles, function(z) names(z$tables))))
  for (nm in table_names) {
    tab <- article_bind(lapply(bundles, function(z) z$tables[[nm]]))
    article_csv(tab, file.path(out, paste0(nm, ".csv")))
    if (nm %in% c("overall", "class", "region", "region_class") && nrow(tab)) {
      tab <- tab[tab$RunID %in% valid$RunID, , drop = FALSE]
      article_csv(article_summary(tab, c(grouping, "Dataset", "Region", "Class")),
                  file.path(out, paste0(nm, "_summary.csv")))
    }
  }
  criteria <- c(SCR_BIC_original = "CriterionSCR_BIC_original", SCR_BIC_effective = "CriterionSCR_BIC_effective",
    SCR_AIC_effective = "CriterionSCR_AIC_effective", CB_BIC = "CriterionCB_BIC", CB_AIC = "CriterionCB_AIC",
    PLIC_BIC = "CriterionPLIC_BIC", PLIC_AIC = "CriterionPLIC_AIC")
  criteria <- criteria[unname(criteria) %in% names(valid)]
  # SCMR-EN-MS shares the G = 1 global fit with SCMR-EN.
  selection_models <- intersect(c("SCMR-EN", "SCMR-EN-MS", "TwoStage-EN"), unique(valid$Model))
  for (sel_model in selection_models) {
  candidate <- if (nrow(valid)) valid[valid$Model %in% c(sel_model, if (sel_model != "SCMR-EN") "SCMR-EN") &
                                      (valid$Model == sel_model | valid$G == 1L), , drop = FALSE] else valid
  if (nrow(candidate)) candidate$Model[candidate$G == 1L] <- sel_model
  suffix <- if (sel_model == "SCMR-EN") "" else paste0("_", gsub("-", "", sel_model))
  # Every requested G must have an eligible fit before a replicate can select G.
  # This prevents failed candidates from silently improving selection frequencies.
  for (criterion in names(criteria)) {
    chosen <- list(); coverage <- list()
    for (id in unique(expected$DatasetID)) {
      z <- if (nrow(candidate)) candidate[candidate$DatasetID == id, , drop = FALSE] else candidate
      finite <- if (nrow(z)) is.finite(z[[criteria[[criterion]]]]) else logical()
      complete <- nrow(z) > 0 && setequal(z$G[finite], cfg$G)
      coverage[[length(coverage) + 1L]] <- data.frame(DatasetID = id,
        Criterion = criterion, NEligible = sum(finite), NExpected = length(cfg$G), Complete = complete)
      if (!complete) next
      z <- z[finite, , drop = FALSE]
      chosen[[length(chosen) + 1L]] <- z[order(z[[criteria[[criterion]]]], z$G)[1], , drop = FALSE]
    }
    best <- article_bind(chosen)
    article_csv(best, file.path(out, paste0("best_G_", criterion, suffix, ".csv")))
    article_csv(article_summary(best, setdiff(grouping, "G")), file.path(out, paste0("summary_selected_", criterion, suffix, ".csv")))
    freq <- data.frame()
    if (nrow(best)) {
      keys <- c(setdiff(grouping, c("Model", "Mode", "Penalty")), "G")
      keys <- unique(keys)
      freq <- stats::aggregate(list(Count = best$Repeat), best[, keys, drop = FALSE], length)
      denom_keys <- setdiff(keys, "G")
      denominator <- stats::aggregate(list(NSelected = best$Repeat), best[, denom_keys, drop = FALSE], length)
      freq <- merge(freq, denominator, by = denom_keys, all.x = TRUE)
      freq$Frequency <- freq$Count / freq$NSelected
    }
    article_csv(freq, file.path(out, paste0("best_G_frequency_", criterion, suffix, ".csv")))
    article_csv(article_bind(coverage), file.path(out, paste0("selection_coverage_", criterion, suffix, ".csv")))
  }
  }
  if (nrow(result)) {
    audit <- result[, c("RunID", "Model", "G", "Converged", "CBBIC_WeightedN_Check", "n_rows"), drop = FALSE]
    audit$WeightedNError <- audit$CBBIC_WeightedN_Check - audit$n_rows
    audit$FiniteCriteria <- apply(result[, unname(criteria), drop = FALSE], 1, function(z) all(is.finite(z)))
    article_csv(audit, file.path(out, "criterion_audit.csv"))
  }
  invisible(list(results = result, coverage = expected, errors = errors))
}
