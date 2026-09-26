# 04 - Spatial autocorrelation of the C1b program ---------------------------
#  * global Moran's I (analytical randomisation test + 999 permutations)
#  * comparison to 1,000 expression-matched random gene sets (empirical null)
#  * macrophage-abundance-adjusted C1b (residual) Moran's I, same null
#  * local Getis-Ord Gi* hotspots (BH-FDR) defining C1b niches
source("R/utils/common.R")
suppressPackageStartupMessages({ library(spdep); library(ggplot2); library(patchwork) })
set.seed(CFG$seed)

meta  <- available_samples()
spots <- readRDS(file.path(P$objects, "spot_scores.rds"))
scfg  <- CFG$spatial
mac_cov <- c("MP_Mac", "MP_Inflammatory.Mac")
feat <- c("C1b", "C1b_noRibo", "C1b_noMPoverlap", "C1b_MHCII", "C1b_nonMHCII", "C1b_AMS", "C1b_Tirosh",
          grep("^MP_", names(spots), value = TRUE))

resid_on <- function(y, covars) residuals(lm(y ~ ., data = as.data.frame(covars)))

moran_rows <- list(); null_rows <- list(); hot <- list()
for (s in meta$sample_id) {
  d <- spots[sample_id == s]
  lw <- visium_listw(d)
  W  <- listw_to_sparse(lw)

  # observed Moran's I for C1b, variants and MPs
  for (f in feat) {
    x <- d[[f]]
    if (is.null(x) || anyNA(x)) next
    mt <- moran.test(x, lw, randomisation = TRUE, zero.policy = TRUE)
    mc <- if (f %in% c("C1b", "C1b_noRibo")) moran.mc(x, lw, nsim = scfg$n_perm_moran, zero.policy = TRUE)$p.value else NA
    moran_rows[[paste(s, f)]] <- data.table(sample_id = s, feature = f,
      I = unname(mt$estimate[1]), E_I = unname(mt$estimate[2]), var_I = unname(mt$estimate[3]),
      z = unname(mt$statistic), p_analytic = mt$p.value, p_perm = mc, n_spots = nrow(d))
  }

  # macrophage-adjusted C1b residual
  X <- d[, ..mac_cov]
  r_obs <- resid_on(d$C1b, X)
  mt <- moran.test(r_obs, lw, randomisation = TRUE, zero.policy = TRUE)
  mc <- moran.mc(r_obs, lw, nsim = scfg$n_perm_moran, zero.policy = TRUE)
  moran_rows[[paste(s, "resid")]] <- data.table(sample_id = s, feature = "C1b_resid_macadj",
      I = unname(mt$estimate[1]), E_I = unname(mt$estimate[2]), var_I = unname(mt$estimate[3]),
      z = unname(mt$statistic), p_analytic = mt$p.value, p_perm = mc$p.value, n_spots = nrow(d))

  # random gene-set null
  rnd <- readRDS(obj_path(s, "random_sets"))$scores
  rnd <- rnd[d$barcode, , drop = FALSE]
  I_null <- moran_matrix(rnd, W)
  rnd_res <- apply(rnd, 2, resid_on, covars = X)
  I_null_res <- moran_matrix(rnd_res, W)
  # C1b vs random: same UCell scale -> also compare score magnitude (mean)
  null_rows[[s]] <- data.table(sample_id = s, set = colnames(rnd), I_random = I_null,
                               I_random_resid = I_null_res, mean_random = colMeans(rnd))

  # Getis-Ord Gi* hotspots
  nb_self <- include.self(lw$neighbours)
  lw_b <- nb2listw(nb_self, style = "B", zero.policy = TRUE)
  gz <- as.numeric(localG(d$C1b, lw_b, zero.policy = TRUE))
  q  <- p.adjust(2 * pnorm(-abs(gz)), "BH")
  hot[[s]] <- data.table(sample_id = s, barcode = d$barcode, Gi_z = gz, Gi_q = q,
                         hotspot = fcase(gz > 0 & q < scfg$hotspot_fdr, "hot",
                                         gz < 0 & q < scfg$hotspot_fdr, "cold", default = "ns"))
  logf("%s: Moran I(C1b)=%.3f; random null median=%.3f; hot spots=%d", s,
       moran_rows[[paste(s, "C1b")]]$I, median(I_null), sum(hot[[s]]$hotspot == "hot"))
}

moran <- rbindlist(moran_rows)
nulls <- rbindlist(null_rows)
hot   <- rbindlist(hot)
saveRDS(hot, file.path(P$objects, "c1b_hotspots.rds"))
fwrite(nulls, file.path(P$objects, "moran_random_null.csv.gz"))

# empirical comparison C1b vs random sets
emp <- merge(moran[feature %in% c("C1b", "C1b_resid_macadj")], nulls[, .(
  null_mean = mean(I_random), null_sd = sd(I_random), null_q95 = quantile(I_random, .95),
  null_mean_res = mean(I_random_resid), null_sd_res = sd(I_random_resid),
  I_rand_vec = list(I_random), I_res_vec = list(I_random_resid)),
  by = sample_id], by = "sample_id")
emp[, `:=`(
  z_vs_random = fifelse(feature == "C1b", (I - null_mean) / null_sd, (I - null_mean_res) / null_sd_res),
  p_vs_random = mapply(function(I, f, a, b) {
    v <- if (f == "C1b") a else b; (1 + sum(v >= I)) / (length(v) + 1) }, I, feature, I_rand_vec, I_res_vec)
)]
emp[, c("I_rand_vec", "I_res_vec") := NULL]
moran <- merge(moran, meta[, .(sample_id, cohort, histology, who_grade, is_LGG, patient)], by = "sample_id")
moran[, q_analytic := p.adjust(p_analytic, "BH"), by = feature]
write_tab(moran, "moran_I_all_features.csv")
emp <- merge(emp, meta[, .(sample_id, cohort, histology, who_grade, is_LGG, patient)], by = "sample_id")
write_tab(emp, "moran_I_C1b_vs_random_sets.csv")
hs <- hot[, .(n_spots = .N, n_hot = sum(hotspot == "hot"), n_cold = sum(hotspot == "cold"),
              frac_hot = mean(hotspot == "hot")), by = sample_id]
write_tab(merge(hs, meta[, .(sample_id, cohort, who_grade)], by = "sample_id"), "c1b_hotspot_summary.csv")

# ---- Figures ---------------------------------------------------------------
ord <- meta[order(cohort == "GBM", who_grade, sample_id)]
lab <- setNames(sprintf("%s\n%s G%d", toupper(ord$sample_id), ord$histology, ord$who_grade), ord$sample_id)

# Fig 1: spatial maps of C1b (IDHm) + hotspots
dd <- merge(spots[, .(sample_id, barcode, array_row, array_col, C1b)], hot, by = c("sample_id", "barcode"))
mk <- function(ids, file, ncol) {
  if (!length(ids)) return(invisible(NULL))          # e.g. replication cohort has no GBM
  pl <- lapply(ids, function(s) {
    x <- dd[sample_id == s][, C1b_z := as.numeric(scale(C1b))]
    a <- spot_plot(x, "C1b_z", lab[s], limits = c(-2.5, 2.5), size = 0.45) +
      labs(colour = "C1b (z)")
    b <- spot_plot(x, "hotspot", "Gi* hotspots", discrete = TRUE, size = 0.45) +
      scale_colour_manual(values = c(hot = "#B2182B", cold = "#2166AC", ns = "grey85"),
                          breaks = c("hot", "cold", "ns"), name = "Gi* (FDR<0.05)") +
      guides(colour = guide_legend(override.aes = list(size = 2.5)))
    a / b
  })
  save_pdf(wrap_plots(pl, ncol = ncol) + plot_layout(guides = "collect"), file, 1.8 * ncol, 3.6 * ceiling(length(ids) / ncol))
}
mk(ord[cohort == "IDHm", sample_id], "Fig1_C1b_spatial_IDHm.pdf", 6)
mk(ord[cohort == "GBM", sample_id], "FigS2_C1b_spatial_GBM.pdf", 7)

# Fig 2a: observed I vs random-set null
nl <- merge(nulls, meta[, .(sample_id, cohort)], by = "sample_id")
nl[, sample_id := factor(sample_id, ord$sample_id)]
ob <- moran[feature == "C1b"][, sample_id := factor(sample_id, ord$sample_id)]
g2 <- ggplot(nl, aes(sample_id, I_random)) +
  geom_violin(fill = "grey80", colour = NA, scale = "width") +
  geom_point(data = ob, aes(sample_id, I, colour = cohort), size = 2) +
  scale_colour_manual(values = c(IDHm = "#D6604D", GBM = "#4393C3")) +
  labs(x = NULL, y = "Global Moran's I", title = "C1b vs 1,000 matched random gene sets") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ob2 <- moran[feature == "C1b_resid_macadj"][, sample_id := factor(sample_id, ord$sample_id)]
g2b <- ggplot(nl, aes(sample_id, I_random_resid)) +
  geom_violin(fill = "grey80", colour = NA, scale = "width") +
  geom_point(data = ob2, aes(sample_id, I, colour = cohort), size = 2) +
  scale_colour_manual(values = c(IDHm = "#D6604D", GBM = "#4393C3")) +
  labs(x = NULL, y = "Moran's I (macrophage-adjusted)", title = "C1b residual after Mac/Inflammatory-Mac MP adjustment") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(g2 / g2b + plot_layout(guides = "collect"), "Fig2_MoranI_vs_random.pdf", 7.5, 5.5)
