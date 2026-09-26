# 05 - Co-localisation of C1b with Greenwald metaprograms / cell composition
#  (a) spot-level Spearman rho, tested against spatially-aware nulls:
#      Moran spectral randomisation (MSR; adespatial::msr) and Dutilleul's
#      modified t-test with effective sample size (SpatialPack)
#  (b) specificity: rho of 1,000 matched random gene sets with each MP
#  (c) macrophage-adjusted partial correlations
#  (d) neighbourhood enrichment of MP labels around C1b hotspots
source("R/utils/common.R")
suppressPackageStartupMessages({
  library(spdep); library(adespatial); library(SpatialPack)
  library(ggplot2); library(patchwork); library(ComplexHeatmap)
})
set.seed(CFG$seed)

meta  <- available_samples()
spots <- readRDS(file.path(P$objects, "spot_scores.rds"))
hot   <- readRDS(file.path(P$objects, "c1b_hotspots.rds"))
scfg  <- CFG$spatial
mp_cols <- grep("^MP_", names(spots), value = TRUE)
mac_cov <- c("MP_Mac", "MP_Inflammatory.Mac")

rk <- function(x) rank(x, ties.method = "average")
resid_rank <- function(y, X) residuals(lm(rk(y) ~ ., data = as.data.frame(apply(X, 2, rk))))

co_rows <- list(); nh_rows <- list(); or_rows <- list()
for (s in meta$sample_id) {
  d <- merge(spots[sample_id == s], hot[sample_id == s, .(barcode, hotspot, Gi_z)], by = "barcode", sort = FALSE)
  lw <- visium_listw(d)
  coords <- cbind(d$array_col, d$array_row * sqrt(3))
  n <- nrow(d)

  # MSR surrogates of rank(C1b) and of rank(C1b | macrophage MPs)
  x_r  <- rk(d$C1b)
  x_rr <- resid_rank(d$C1b, d[, ..mac_cov])
  sur  <- adespatial::msr(x_r, lw, nrepet = scfg$n_msr, method = "pair")
  sur_r <- adespatial::msr(x_rr, lw, nrepet = scfg$n_msr, method = "pair")

  rnd <- readRDS(obj_path(s, "random_sets"))$scores[d$barcode, , drop = FALSE]
  rnd_r <- apply(rnd, 2, rk)

  for (m in mp_cols) {
    y_r <- rk(d[[m]])
    rho <- cor(x_r, y_r)
    null <- as.numeric(cor(sur, y_r))
    p_msr <- (1 + sum(abs(null) >= abs(rho))) / (length(null) + 1)
    mt <- tryCatch(SpatialPack::modified.ttest(x_r, y_r, coords, nclass = 13),
                   error = function(e) NULL)
    # partial (macrophage-adjusted); not defined for the covariates themselves
    if (m %in% mac_cov) {
      prho <- NA; p_msr_partial <- NA
    } else {
      y_rr <- resid_rank(d[[m]], d[, ..mac_cov])
      prho <- cor(x_rr, y_rr)
      nullp <- as.numeric(cor(sur_r, y_rr))
      p_msr_partial <- (1 + sum(abs(nullp) >= abs(prho))) / (length(nullp) + 1)
    }
    rho_rand <- as.numeric(cor(rnd_r, y_r))
    co_rows[[paste(s, m)]] <- data.table(
      sample_id = s, MP = sub("^MP_", "", m), n_spots = n, rho = rho,
      p_msr = p_msr, msr_null_sd = sd(null),
      ess = if (is.null(mt)) NA_real_ else mt$ESS,
      p_dutilleul = if (is.null(mt)) NA_real_ else mt$p.value,
      rho_partial_macadj = prho, p_msr_partial = p_msr_partial,
      rho_random_mean = mean(rho_rand), rho_random_sd = sd(rho_rand),
      pct_random_below = mean(rho_rand < rho))
  }

  # ---- neighbourhood enrichment around C1b hotspots -------------------------
  lab <- d$mp_label
  H <- which(d$hotspot == "hot")
  if (length(H) >= 5) {
    nbh <- sort(unique(c(H, unlist(lw$neighbours[H]))))
    nbh <- nbh[nbh > 0]
    levs <- sort(unique(lab))
    obs <- table(factor(lab[nbh], levs))
    perm <- replicate(scfg$n_perm_nhood, table(factor(sample(lab)[nbh], levs)))
    z <- (as.numeric(obs) - rowMeans(perm)) / pmax(apply(perm, 1, sd), 1e-9)
    p <- (1 + rowSums(perm >= as.numeric(obs))) / (ncol(perm) + 1)
    p_dep <- (1 + rowSums(perm <= as.numeric(obs))) / (ncol(perm) + 1)
    nh_rows[[s]] <- data.table(sample_id = s, MP = levs, n_niche = length(nbh),
                               obs = as.numeric(obs), exp = rowMeans(perm), z = z,
                               p_enrich = p, p_deplete = p_dep)
    # odds ratio: label within hotspots vs rest (Fisher)
    or_rows[[s]] <- rbindlist(lapply(levs, function(l) {
      tb <- table(factor(d$hotspot == "hot", c(TRUE, FALSE)), factor(lab == l, c(TRUE, FALSE)))
      ft <- fisher.test(tb)
      data.table(sample_id = s, MP = l, OR = unname(ft$estimate), p = ft$p.value,
                 frac_in_hot = tb[1, 1] / sum(tb[1, ]), frac_in_rest = tb[2, 1] / sum(tb[2, ]))
    }))
  }
  logf("%s: co-localisation done (%d MPs, %d hotspot spots)", s, length(mp_cols), length(H))
}

co <- merge(rbindlist(co_rows), meta[, .(sample_id, cohort, histology, who_grade, is_LGG, patient)], by = "sample_id")
co[, `:=`(q_msr = p.adjust(p_msr, "BH"), q_dutilleul = p.adjust(p_dutilleul, "BH"),
          q_msr_partial = p.adjust(p_msr_partial, "BH")), by = sample_id]
write_tab(co, "coloc_C1b_vs_MP_per_sample.csv")
nh <- merge(rbindlist(nh_rows), meta[, .(sample_id, cohort, who_grade)], by = "sample_id")
nh[, q_enrich := p.adjust(p_enrich, "BH"), by = sample_id]
write_tab(nh, "nhood_enrichment_C1b_hotspots.csv")
orr <- merge(rbindlist(or_rows), meta[, .(sample_id, cohort)], by = "sample_id")
orr[, q := p.adjust(p, "BH"), by = sample_id]
write_tab(orr, "hotspot_MP_composition_fisher.csv")

# ---- Figures ---------------------------------------------------------------
ord <- meta[order(cohort == "GBM", who_grade, sample_id)]
mk_heat <- function(dt, val, qcol, title, file, lim = 0.8) {
  M <- dcast(dt, MP ~ sample_id, value.var = val); mm <- as.matrix(M[, -1]); rownames(mm) <- M$MP
  Q <- dcast(dt, MP ~ sample_id, value.var = qcol); qq <- as.matrix(Q[, -1]); rownames(qq) <- Q$MP
  cols <- intersect(ord$sample_id, colnames(mm)); mm <- mm[, cols]; qq <- qq[rownames(mm), cols]
  ann <- HeatmapAnnotation(Cohort = ord[match(cols, sample_id), cohort],
                           Grade = factor(ord[match(cols, sample_id), who_grade]),
                           col = list(Cohort = c(IDHm = "#D6604D", GBM = "#4393C3"),
                                      Grade = c(`2` = "#FEE0B6", `3` = "#F1A340", `4` = "#B35806")),
                           annotation_name_gp = grid::gpar(fontsize = 7), simple_anno_size = grid::unit(2.5, "mm"))
  hm <- Heatmap(mm, name = val, col = circlize::colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B")),
                cluster_columns = FALSE, top_annotation = ann, column_title = title,
                column_split = ord[match(cols, sample_id), cohort],
                row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  q <- qq[i, j]
                  if (!is.na(q) && q < 0.05) grid::grid.text(ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", "*")),
                                                             x, y, gp = grid::gpar(fontsize = 6))
                })
  f <- file.path(P$figures, file)
  grDevices::cairo_pdf(f, width = 7.5, height = 3.8); draw(hm); dev.off()
  png(sub("pdf$", "png", f), width = 7.5, height = 3.8, units = "in", res = 200); draw(hm); dev.off()
  logf("saved figure %s", file)
}
mk_heat(co, "rho", "q_msr", "C1b vs MP: Spearman rho (MSR-tested, BH)", "Fig3a_coloc_rho_heatmap.pdf")
mk_heat(co[!MP %in% sub("MP_", "", mac_cov)], "rho_partial_macadj", "q_msr_partial",
        "Partial rho | Mac + Inflammatory-Mac", "Fig3b_coloc_partial_rho_heatmap.pdf")
mk_heat(nh, "z", "q_enrich", "MP label enrichment in C1b hotspot niches (z)", "Fig4_nhood_enrichment_heatmap.pdf", lim = 10)
