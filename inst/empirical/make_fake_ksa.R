# Synthetic stand-in for df_month_rec_final.csv (same column conventions), used
# only to exercise SCMR_EN_empirical_article_v3.R in continuous integration.
# Usage: Rscript make_fake_ksa.R <output.csv>
suppressPackageStartupMessages(library(scmr))
out <- commandArgs(trailingOnly = TRUE)[1]
set.seed(1)
reg <- c("3205", "3212", "6171")
centre <- rbind(c(-7.4, 107.6), c(-6.4, 108.2), c(-0.05, 109.3))
seg <- do.call(rbind, lapply(seq_along(reg), function(r) data.frame(
  id_segmen = sprintf("%s%05d", reg[r], seq_len(12)),
  lat = centre[r, 1] + runif(12, -0.15, 0.15), lon = centre[r, 2] + runif(12, -0.15, 0.15))))
sub <- do.call(rbind, lapply(seq_len(nrow(seg)), function(i) {
  g <- expand.grid(a = 1:3, b = 1:3)
  data.frame(id_segmen = seg$id_segmen[i], id_subsegmen = paste0(seg$id_segmen[i], LETTERS[g$a], g$b),
             lat = seg$lat[i] + 0.001 * g$a, lon = seg$lon[i] + 0.001 * g$b)
}))
sim <- simulate_scmr_panel(coords = cbind(sub$lat, sub$lon), n_time = 24, pattern = "blocks", G = 3,
                           n_classes = 6, p_active = 4, p_inactive = 4, delta = 1.5, seed = 3)
u <- match(sim$unit_id, sprintf("U%04d", seq_len(nrow(sub))))
phase <- as.integer(as.character(sim$y))
other <- phase == 6L
phase[other] <- sample(6:8, sum(other), replace = TRUE)
phase[sample(length(phase), round(0.01 * length(phase)))] <- 12L
year <- 2021L + (sim$time - 1L) %/% 12L
month <- (sim$time - 1L) %% 12L + 1L
d <- data.frame(id_segmen = sub$id_segmen[u], id_subsegmen = sub$id_subsegmen[u],
                id_month = sprintf("%s%d%02d", sub$id_subsegmen[u], year, month),
                phase = sprintf("%.1f", phase), STRATA = 1L,
                lati_mean = sub$lat[u], long_mean = sub$lon[u],
                s2r_med_NDVI_stms = sim$x[, 1], s2r_med_EVI_stms = sim$x[, 2],
                s1_med_vh_stms = sim$x[, 3], s2r_month_NDVI_stms_slope1 = sim$x[, 4],
                s2r_med_NDWI_stms = sim$x[, 5], s2r_med_MSI_stms = sim$x[, 6],
                s1_med_vv_stms = sim$x[, 7], s1_month_vh_slope1 = sim$x[, 8],
                s2r_last_NDVI_stms = rnorm(nrow(sim$x)),
                year = year, month = month, month_end = sprintf("%d-%02d-28", year, month))
d$s1_med_vv_stms[sample(nrow(d), 30)] <- NA
utils::write.csv(d, out, row.names = FALSE)
cat("Wrote", nrow(d), "rows to", out, "\n")
