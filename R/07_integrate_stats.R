# 07 - Cross-sample integration: random-effects meta-analysis, cohort
#      comparison, seed robustness, summary figures, session info
source("R/utils/common.R")
suppressPackageStartupMessages({ library(metafor); library(spdep); library(ggplot2); library(patchwork) })

meta  <- available_samples()
moran <- fread(file.path(P$tables, "moran_I_all_features.csv"))
emp   <- fread(file.path(P$tables, "moran_I_C1b_vs_random_sets.csv"))
co    <- fread(file.path(P$tables, "coloc_C1b_vs_MP_per_sample.csv"))
nh    <- fread(file.path(P$tables, "nhood_enrichment_C1b_hotspots.csv"))
cohorts <- list(IDHm = meta[cohort == "IDHm", sample_id], LGG = meta[is_LGG == TRUE, sample_id],
                GBM = meta[cohort == "GBM", sample_id])

rma_ml <- function(dt, yi, vi) {
  # multilevel RE model (sections nested in patients) when patients repeat
  if (uniqueN(dt$patient) < nrow(dt)) {
    rma.mv(dt[[yi]], dt[[vi]], random = ~ 1 | patient / sample_id, data = dt, method = "REML")
  } else rma(dt[[yi]], dt[[vi]], method = "REML")
}
tidy_rma <- function(f) data.table(estimate = as.numeric(f$b), se = f$se, ci_lb = f$ci.lb, ci_ub = f$ci.ub,
                                   p = f$pval, k = f$k, I2 = if (!is.null(f$I2)) f$I2 else NA_real_)

# ---- 1. Moran's I meta-analysis ---------------------------------------------
mi <- rbindlist(lapply(c("C1b", "C1b_noRibo", "C1b_noMPoverlap", "C1b_resid_macadj", "MP_Mac", "MP_MES.Hyp"), function(f) {
  rbindlist(lapply(names(cohorts), function(ch) {
    dt <- moran[feature == f & sample_id %in% cohorts[[ch]]]
    if (nrow(dt) < 2) return(NULL)
    cbind(data.table(feature = f, cohort = ch), tidy_rma(rma_ml(dt, "I", "var_I")))
  }))
}))
write_tab(mi, "meta_moran_I.csv")

# cohort difference (IDHm vs GBM) for C1b Moran's I, multilevel meta-regression
dm <- moran[feature == "C1b"]
dm[, cohort := factor(cohort, c("GBM", "IDHm"))]
mr <- rma.mv(I, var_I, mods = ~ cohort, random = ~ 1 | patient / sample_id, data = dm, method = "REML")
write_tab(data.table(term = rownames(mr$b), estimate = as.numeric(mr$b), se = mr$se, p = mr$pval),
          "meta_regression_moran_C1b_IDHm_vs_GBM.csv")

# random-set z-scores combined (Stouffer) per cohort
st <- rbindlist(lapply(names(cohorts), function(ch) {
  emp[sample_id %in% cohorts[[ch]], .(cohort = ch, n = .N, median_z = median(z_vs_random),
     n_p_lt_0.05 = sum(p_vs_random < 0.05), stouffer_z = sum(qnorm(1 - p_vs_random)) / sqrt(.N),
     stouffer_p = pnorm(-sum(qnorm(1 - p_vs_random)) / sqrt(.N))), by = feature]
}))
write_tab(st, "moran_vs_random_stouffer.csv")

# C1b rank among all signatures in spatial coherence
rk <- moran[grepl("^MP_|^C1b$", feature), .(feature, I, rank = frank(-I)), by = sample_id][feature == "C1b"]
write_tab(merge(rk, meta[, .(sample_id, cohort)], by = "sample_id"), "moran_C1b_rank_among_MPs.csv")

# ---- 2. Co-localisation meta-analysis (Fisher z, ESS-based variance) -------
co[, `:=`(zr = atanh(rho), vz = 1 / pmax(ess - 3, 1),
          zr_p = atanh(rho_partial_macadj))]
cm <- rbindlist(lapply(names(cohorts), function(ch) {
  rbindlist(lapply(unique(co$MP), function(m) {
    dt <- co[MP == m & sample_id %in% cohorts[[ch]] & is.finite(vz)]
    if (nrow(dt) < 2) return(NULL)
    a <- tidy_rma(rma_ml(dt, "zr", "vz"))
    b <- if (all(is.na(dt$zr_p))) NULL else tidy_rma(rma_ml(dt[!is.na(zr_p)], "zr_p", "vz"))
    out <- cbind(data.table(cohort = ch, MP = m, type = "rho"), a)
    if (!is.null(b)) out <- rbind(out, cbind(data.table(cohort = ch, MP = m, type = "partial_rho_macadj"), b))
    out
  }))
}))
cm[, `:=`(rho = tanh(estimate), rho_lb = tanh(ci_lb), rho_ub = tanh(ci_ub))]
cm[, q := p.adjust(p, "BH"), by = .(cohort, type)]
write_tab(cm, "meta_coloc_C1b_vs_MP.csv")

# specificity vs random gene sets: in how many samples does C1b-MP rho exceed
# (or fall below) 95% of the matched random sets? exact binomial test vs 5%
sp <- rbindlist(lapply(names(cohorts), function(ch) {
  co[sample_id %in% cohorts[[ch]], .(cohort = ch, k = .N,
     n_above95 = sum(pct_random_below > 0.95), n_below5 = sum(pct_random_below < 0.05),
     mean_rho_C1b = mean(rho), mean_rho_random = mean(rho_random_mean)), by = MP]
}))
sp[, `:=`(p_above = mapply(function(x, n) binom.test(x, n, 0.05, "greater")$p.value, n_above95, k),
          p_below = mapply(function(x, n) binom.test(x, n, 0.05, "greater")$p.value, n_below5, k))]
sp[, `:=`(q_above = p.adjust(p_above, "BH"), q_below = p.adjust(p_below, "BH")), by = cohort]
write_tab(sp, "coloc_specificity_vs_random_sets.csv")

# ---- 3. Neighbourhood enrichment combined (Stouffer on z) --------------------
nm <- rbindlist(lapply(names(cohorts), function(ch) {
  nh[sample_id %in% cohorts[[ch]], .(cohort = ch, k = .N, mean_z = mean(z),
     stouffer_z = sum(qnorm(1 - p_enrich)) / sqrt(.N)), by = MP]
}))
nm[, p := pnorm(-stouffer_z)][, q := p.adjust(p, "BH"), by = cohort]
write_tab(nm, "meta_nhood_enrichment.csv")

# ---- 4. Seed robustness of permutation p-values -----------------------------
spots <- readRDS(file.path(P$objects, "spot_scores.rds"))
sr <- rbindlist(lapply(meta[cohort == "IDHm", sample_id], function(s) {
  d <- spots[sample_id == s]; lw <- visium_listw(d)
  rbindlist(lapply(c(1, 2, 3), function(k) {
    set.seed(CFG$seed * k)
    data.table(sample_id = s, seed = CFG$seed * k,
               p_perm = moran.mc(d$C1b, lw, nsim = 999, zero.policy = TRUE)$p.value)
  }))
}))
write_tab(sr, "seed_robustness_moran_perm.csv")

# ---- Figures ---------------------------------------------------------------
fp <- cm[type == "rho" & cohort %in% c("IDHm", "GBM")]
fp[, MP := factor(MP, fp[cohort == "IDHm"][order(rho), MP])]
g1 <- ggplot(fp, aes(rho, MP, colour = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .3) +
  geom_pointrange(aes(xmin = rho_lb, xmax = rho_ub), position = position_dodge(width = .6), size = .2) +
  scale_colour_manual(values = c(IDHm = "#D6604D", GBM = "#4393C3")) +
  labs(x = "Pooled Spearman rho (95% CI, RE meta-analysis)", y = NULL, title = "C1b co-localisation with MPs") +
  theme_pub()
fp2 <- cm[type == "partial_rho_macadj" & cohort %in% c("IDHm", "GBM")]
fp2[, MP := factor(MP, levels(fp$MP))]
g2 <- ggplot(fp2, aes(rho, MP, colour = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .3) +
  geom_pointrange(aes(xmin = rho_lb, xmax = rho_ub), position = position_dodge(width = .6), size = .2) +
  scale_colour_manual(values = c(IDHm = "#D6604D", GBM = "#4393C3")) +
  labs(x = "Pooled partial rho | Mac, Inflammatory-Mac", y = NULL, title = "Macrophage-adjusted") + theme_pub()
save_pdf(g1 + g2 + plot_layout(guides = "collect"), "Fig3c_coloc_meta_forest.pdf", 8, 3.8)

mm <- moran[feature %in% c("C1b", "C1b_resid_macadj")]
mm[, sample_id := factor(sample_id, meta[order(cohort == "GBM", who_grade, sample_id), sample_id])]
g3 <- ggplot(mm, aes(I, sample_id, colour = cohort, shape = feature)) +
  geom_errorbarh(aes(xmin = I - 1.96 * sqrt(var_I), xmax = I + 1.96 * sqrt(var_I)), height = 0, linewidth = .3) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = c(IDHm = "#D6604D", GBM = "#4393C3")) +
  labs(x = "Global Moran's I (95% CI)", y = NULL) + theme_pub()
save_pdf(g3, "Fig2b_MoranI_forest.pdf", 4.5, 4.2)

writeLines(capture.output(sessionInfo()), file.path(P$logs, "sessionInfo.txt"))
logf("integration done")
