# Ulangan E1 (nowcast) dan E3 (imputasi) dengan beberapa split, dalam 4 proses
# R paralel. Split 1 adalah run v4 yang sudah ada dan dipakai ulang; yang
# difit baru: E1 split 2 (latih bulan 1-15, uji 16-21), E1 split 3 (latih
# 1-21, uji 22-24), dan E3 split 2-5 (penyembunyian acak 15% yang berbeda).
# Jalankan seluruh file ini (Ctrl+Shift+S di RStudio). Proses berjalan di latar
# belakang; panggil status() untuk memantau. Perkiraan waktu: 8-10 jam.
syntax_dir <- "C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax"
script <- file.path(syntax_dir, "SCMR_EN_empirical_article_v4.R")
log_dir <- file.path(dirname(syntax_dir), "Results", "Real data v4")
rscript <- file.path(R.home("bin"), "Rscript.exe")
jobs <- list(
  E1_s2  = c(SCMR_DESIGNS = "E1_nowcast",    SCMR_E1_FOLDS = "2",   SCMR_E3_FOLDS = ""),
  E1_s3  = c(SCMR_DESIGNS = "E1_nowcast",    SCMR_E1_FOLDS = "3",   SCMR_E3_FOLDS = ""),
  E3_s23 = c(SCMR_DESIGNS = "E3_imputation", SCMR_E1_FOLDS = "",    SCMR_E3_FOLDS = "2,3"),
  E3_s45 = c(SCMR_DESIGNS = "E3_imputation", SCMR_E1_FOLDS = "",    SCMR_E3_FOLDS = "4,5")
)
for (nm in names(jobs)) {
  Sys.unsetenv(c("SCMR_E1_FOLDS", "SCMR_E3_FOLDS"))
  v <- jobs[[nm]]; v <- v[nzchar(v)]
  do.call(Sys.setenv, as.list(v))
  system2(rscript, shQuote(script), stdout = file.path(log_dir, paste0("log_", nm, ".txt")),
          stderr = file.path(log_dir, paste0("err_", nm, ".txt")), wait = FALSE)
  cat("Diluncurkan:", nm, "\n")
  Sys.sleep(5)
}
Sys.unsetenv(c("SCMR_DESIGNS", "SCMR_E1_FOLDS", "SCMR_E3_FOLDS"))

status <- function() {
  for (nm in names(jobs)) {
    f <- file.path(log_dir, paste0("log_", nm, ".txt"))
    e <- file.path(log_dir, paste0("err_", nm, ".txt"))
    x <- tryCatch(readLines(f, warn = FALSE), error = function(err) NULL, warning = function(w) NULL)
    st <- if (is.null(x)) "MASIH BERJALAN" else if (any(grepl("^Outputs in:", x))) "SELESAI" else "BERHENTI (lihat err)"
    cat(sprintf("%-7s %-22s log %6.1f KB | err %5.1f KB | diubah %s\n", nm, st,
                file.size(f) / 1024, file.size(e) / 1024, format(file.mtime(f), "%H:%M")))
  }
  cat("Proses R aktif:", length(grep("Rterm|Rscript", system("tasklist", intern = TRUE))), "\n")
}
status()
# Setelah keempatnya SELESAI, jalankan bootstrap ulang (beberapa menit):
# Sys.setenv(SCMR_INTERP_STAGE = "boot", SCMR_BOOT = "1000")
# source(file.path(syntax_dir, "SCMR_EN_interpret_v4.R"))
