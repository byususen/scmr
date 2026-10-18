# Menjalankan tafsir rezim KSA (SCMR_EN_interpret_v4.R) di latar belakang.
# Prasyarat: hasil v4 sudah lengkap di Results/Real data v4 (row_predictions.csv,
# results.csv) dan paket scmr terbaru sudah terinstal.
# Jalankan seluruh file ini (Ctrl+Shift+S di RStudio). Satu proses latar belakang
# menyetel model global dan penalti lokal pada seluruh segmen, lalu memfit
# G = 1, 3, 6, 8, 10, 12, 15 dalam 4 proses pekerja, lalu menyusun laporan dan
# bootstrap. Jendela R ini tetap bisa dipakai; panggil status() untuk memantau.
syntax_dir <- "C:/Users/bayus/OneDrive/S3/Disertasi/Article 3/Bahan/Syntax"
script <- file.path(syntax_dir, "SCMR_EN_interpret_v4.R")
log_dir <- file.path(dirname(syntax_dir), "Results", "Real data v4")
int_dir <- file.path(log_dir, "interpret")
rscript <- file.path(R.home("bin"), "Rscript.exe")

Sys.setenv(SCMR_INTERP_G = "1,3,6,8,10,12,15", SCMR_INTERP_WORKERS = "4", SCMR_BOOT = "1000")
system2(rscript, shQuote(script), stdout = file.path(log_dir, "log_interpret.txt"),
        stderr = file.path(log_dir, "err_interpret.txt"), wait = FALSE)
Sys.unsetenv(c("SCMR_INTERP_G", "SCMR_INTERP_WORKERS", "SCMR_BOOT"))
cat("Diluncurkan. Log:", file.path(log_dir, "log_interpret.txt"), "\n")

status <- function() {
  f <- file.path(log_dir, "log_interpret.txt")
  x <- tryCatch(readLines(f, warn = FALSE), error = function(e) NULL, warning = function(w) NULL)
  st <- if (is.null(x)) "MASIH BERJALAN" else
    if (file.exists(file.path(int_dir, "bootstrap_E2_pooled_kappa.csv"))) "SELESAI" else "BERHENTI (lihat err_interpret.txt)"
  cat("Status:", st, "\n")
  fits <- list.files(file.path(int_dir, "fits"), pattern = "\\.rds$")
  cat("Fit selesai:", if (length(fits)) paste(sub("SCMR_EN_S_full_|\\.rds", "", fits), collapse = " ") else "-", "\n")
  w <- file.path(int_dir, "workers_log.txt")
  if (file.exists(w)) {
    y <- tryCatch(readLines(w, warn = FALSE), error = function(e) character(0), warning = function(e) character(0))
    y <- grep("SCMR-EN-S G =|done in|FAILED", y, value = TRUE)
    if (length(y)) cat("Pekerja:\n", paste(" ", utils::tail(y, 8), collapse = "\n"), "\n")
  }
  cat("Proses R aktif:", length(grep("Rterm|Rscript", system("tasklist", intern = TRUE))), "\n")
}
status()
# Selesai bila status() menulis SELESAI. Hasil ada di Results/Real data v4/interpret.
