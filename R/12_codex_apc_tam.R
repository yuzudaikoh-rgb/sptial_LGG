# 12 - Single-cell protein validation on CODEX (Hoefflin et al., Cancer Cell 2026)
#
#  30 IDH-mutant sections (WHO 2-4) + 7 GBM sections, 3.0 M segmented cells,
#  66-plex CODEX with author cell-type calls; Nimbus scores (per-cell
#  probability of marker positivity, 0-1) are used as marker readouts.
#  Coordinates are in micrometres (median nearest-neighbour distance ~12 um).
#
#  Definitions (within the author-annotated TAM compartment, cell_type == "Mac"):
#    C1b-protein score = mean within-sample z of MHC-II, CD163, CD206, CD44, VIM
#                        (the five C1b members on the panel; validated as a
#                        transcript-level proxy of the 90-gene program in step 10)
#    C1b-high / C1b-low TAM = top / bottom tertile of that score within section
#    MHC-II-high / -low TAM = top / bottom tertile of MHC-II within section
#    (absolute Nimbus positivity varies 5-99% between sections - staining
#    batch - so only within-section relative levels are compared)
#  Questions
#    1. Are C1b-high (and MHC-II-high) TAMs closer to CD4+ / CD8+ T cells than other
#       TAMs of the same section?  (label permutation among TAMs = null that
#       keeps TAM positions; plus a per-section logistic model adjusting for
#       local TAM density, local cell density and distance to vessels)
#    2. Are they perivascular and away from hypoxia (CA9+ cells)?
#    3. Are T cells next to MHC-II-high TAMs in a different state (PD-1+, CD69+)?
#    4. Does this depend on WHO grade (IDHm 2 / 3 / 4) and differ from GBM?
#  Per-section effects are pooled with random-effects models (sections nested
#  in patients) per grade group.
source("R/utils/common.R")
suppressPackageStartupMessages({
  library(arrow); library(hdf5r); library(RANN); library(metafor)
  library(ggplot2); library(patchwork)
})
set.seed(CFG$seed)
CX  <- file.path(ROOT, "data", "raw", "codex")
RADII <- c(15, 27.5, 55)
R0 <- 27.5
N_PERM <- 1000
MARK <- c("MHCII", "CD163", "CD206", "CD44", "VIM", "CD279", "CD69", "CA9", "CD31", "CD68" )
P5 <- c("MHCII", "CD163", "CD206", "CD44", "VIM")

read_nimbus <- function(h5, ds, markers) {
  f <- H5File$new(h5, "r"); on.exit(f$close_all())
  rn <- f[["rownames"]][]; cn <- f[["colnames"]][]
  idx <- match(intersect(markers, rn), rn)
  m <- sapply(idx, function(i) f[[ds]][i, ])
  colnames(m) <- rn[idx]
  data.table(cell_name = cn, m)
}

# ---- load ----------------------------------------------------------------------
idh <- as.data.table(read_parquet(file.path(CX, "cells_df.parquet"),
        col_select = c("cell_name", "sample", "PatientID", "GradeTumor", "GradeSection", "GliomaType",
                       "centroid_x", "centroid_y", "cell_type", "cell_type2", "ivygap")))
setnames(idh, "PatientID", "patient")
idh <- merge(idh, read_nimbus(file.path(CX, "CODEX_IDHm_nimbus_scores.h5"), "CODEX_IDHm_nimbus_scores", MARK),
             by = "cell_name")
gbm <- as.data.table(read_parquet(file.path(CX, "gbm_cells_df.parquet"),
        col_select = c("cell_name", "sample", "patient", "centroid_x", "centroid_y", "cell_type")))
gbm <- merge(gbm, read_nimbus(file.path(CX, "CODEX_IDHm_nimbus_scores_gbm.h5"), "CODEX_IDHm_nimbus_scores_gbm", MARK),
             by = "cell_name")
gbm[, `:=`(GradeTumor = 4L, GradeSection = 4L, GliomaType = "GBM", cell_type2 = cell_type, ivygap = NA_character_)]
cells <- rbind(idh, gbm, fill = TRUE)
cells <- cells[!cell_type %in% c("excluded", "low")]
cells[, group := fcase(GliomaType == "GBM", "GBM",
                       GradeTumor == 2, "IDHm G2", GradeTumor == 3, "IDHm G3", default = "IDHm G4")]
cells[, is_T := cell_type %in% c("TcellCD4", "TcellCD8")]
logf("CODEX: %d cells, %d sections (%d IDHm), %d TAMs, %d CD4 T, %d CD8 T", nrow(cells), uniqueN(cells$sample),
     uniqueN(cells[GliomaType != "GBM", sample]), sum(cells$cell_type == "Mac"),
     sum(cells$cell_type == "TcellCD4"), sum(cells$cell_type == "TcellCD8"))

count_within <- function(ref_xy, qry_xy, r) {
  if (!nrow(ref_xy)) return(rep(0L, nrow(qry_xy)))
  k <- min(nrow(ref_xy), 200)
  nn <- nn2(ref_xy, qry_xy, k = k, searchtype = "radius", radius = r)
  rowSums(nn$nn.idx > 0)
}
nearest <- function(ref_xy, qry_xy) {
  if (!nrow(ref_xy)) return(rep(NA_real_, nrow(qry_xy)))
  as.numeric(nn2(ref_xy, qry_xy, k = 1)$nn.dists)
}

# ---- per-section computation -------------------------------------------------------
tam_rows <- list(); sec_rows <- list(); tstate_rows <- list(); comp_rows <- list()
for (s in unique(cells$sample)) {
  x <- cells[sample == s]; xy <- as.matrix(x[, .(centroid_x, centroid_y)])
  mac <- which(x$cell_type == "Mac")
  if (length(mac) < 100) next
  z <- sapply(P5, function(m) as.numeric(scale(x[[m]][mac])))
  z[is.na(z)] <- 0
  score <- rowMeans(z)
  terc <- cut(score, quantile(score, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low", "mid", "high"))
  mt <- cut(rank(x$MHCII[mac], ties.method = "random"), 3, labels = c("low", "mid", "high"))
  apc <- mt == "high"
  sel <- function(v) xy[which(v), , drop = FALSE]
  T4 <- sel(x$cell_type == "TcellCD4"); T8 <- sel(x$cell_type == "TcellCD8"); TT <- sel(x$is_T)
  V  <- sel(x$cell_type == "Vasc"); H <- sel(x$CA9 > 0.5)
  q  <- xy[mac, , drop = FALSE]
  tm <- data.table(sample = s, cell_name = x$cell_name[mac], group = x$group[1], patient = x$patient[1],
                   C1b_score = score, C1b_tertile = terc, MHCII_tertile = mt, subtype = x$cell_type2[mac],
                   n_all55 = count_within(xy, q, 55) - 1L, n_mac55 = count_within(xy[mac, , drop = FALSE], q, 55) - 1L,
                   d_T = nearest(TT, q), d_CD4 = nearest(T4, q), d_CD8 = nearest(T8, q),
                   d_vessel = nearest(V, q), d_hypoxia = nearest(H, q))
  for (r in RADII) {
    set(tm, j = paste0("nT_", r), value = count_within(TT, q, r))
    set(tm, j = paste0("nCD4_", r), value = count_within(T4, q, r))
    set(tm, j = paste0("nCD8_", r), value = count_within(T8, q, r))
  }
  tam_rows[[s]] <- tm

  # composition of C1b-high TAMs by author TAM subtype
  comp_rows[[s]] <- tm[, .N, by = .(sample, group, C1b_tertile, subtype)]

  # T-cell state next to MHC-II-high TAMs
  if (nrow(TT) >= 10) {
    ti <- which(x$is_T)
    nA <- count_within(q[which(apc), , drop = FALSE], xy[ti, , drop = FALSE], R0)
    nN <- count_within(q[which(mt == "low"), , drop = FALSE], xy[ti, , drop = FALSE], R0)
    tstate_rows[[s]] <- data.table(sample = s, group = x$group[1], patient = x$patient[1], cell_type = x$cell_type[ti],
                                   near_APCpos = nA > 0, near_APCneg_only = nA == 0 & nN > 0, no_TAM = nA == 0 & nN == 0,
                                   PD1 = x$CD279[ti] > 0.5, CD69 = x$CD69[ti] > 0.5)
  }

  # ---- section-level statistics ---------------------------------------------------
  st <- list()
  for (contrast in c("C1b_high_vs_low", "MHCII_high_vs_low")) {
    tert <- if (contrast == "C1b_high_vs_low") tm$C1b_tertile else tm$MHCII_tertile
    g <- ifelse(tert == "high", 1L, ifelse(tert == "low", 0L, NA_integer_))
    ok <- !is.na(g)
    if (sum(g[ok] == 1) < 20 || sum(g[ok] == 0) < 20) next
    for (tt in c("T", "CD4", "CD8")) {
      n <- tm[[paste0("n", tt, "_", R0)]][ok]; gg <- g[ok]
      if (sum(n > 0) < 5) next
      obs <- mean(n[gg == 1]) - mean(n[gg == 0])
      perm <- replicate(N_PERM, { p <- sample(gg); mean(n[p == 1]) - mean(n[p == 0]) })
      dfm <- data.table(y = as.integer(n > 0), g = gg, lmac = log1p(tm$n_mac55[ok]), lall = log1p(tm$n_all55[ok]),
                        lves = log1p(pmin(tm$d_vessel[ok], 1000)))
      fit <- tryCatch(glm(y ~ g + lmac + lall + lves, family = binomial, data = dfm), error = function(e) NULL)
      co <- if (!is.null(fit) && "g" %in% rownames(summary(fit)$coefficients)) summary(fit)$coefficients["g", ] else rep(NA, 4)
      st[[paste(contrast, tt)]] <- data.table(sample = s, group = x$group[1], patient = x$patient[1], contrast = contrast,
        target = tt, n_hi = sum(gg == 1), n_lo = sum(gg == 0),
        frac_hi = mean(n[gg == 1] > 0), frac_lo = mean(n[gg == 0] > 0),
        mean_diff = obs, p_perm = (1 + sum(abs(perm) >= abs(obs))) / (N_PERM + 1),
        logOR_adj = co[1], se_adj = co[2])
    }
    for (dd in c("d_vessel", "d_hypoxia")) {
      v <- log10(tm[[dd]][ok] + 1); gg <- g[ok]
      if (all(is.na(v))) next
      st[[paste(contrast, dd)]] <- data.table(sample = s, group = x$group[1], patient = x$patient[1], contrast = contrast,
        target = dd, n_hi = sum(gg == 1), n_lo = sum(gg == 0), frac_hi = NA_real_, frac_lo = NA_real_,
        mean_diff = mean(v[gg == 1], na.rm = TRUE) - mean(v[gg == 0], na.rm = TRUE),
        p_perm = wilcox.test(v[gg == 1], v[gg == 0])$p.value,
        logOR_adj = mean(v[gg == 1], na.rm = TRUE) - mean(v[gg == 0], na.rm = TRUE),
        se_adj = sqrt(var(v[gg == 1], na.rm = TRUE) / sum(gg == 1) + var(v[gg == 0], na.rm = TRUE) / sum(gg == 0)))
    }
  }
  sec_rows[[s]] <- rbindlist(st)
  logf("CODEX %s (%s): %d TAMs (MHC-II Nimbus>0.5: %.0f%%), %d T cells; C1b-high vs low T-contact %.3f vs %.3f",
       s, x$group[1], length(mac), 100 * mean(x$MHCII[mac] > 0.5), nrow(TT),
       mean(tm$nT_27.5[tm$C1b_tertile == "high"] > 0), mean(tm$nT_27.5[tm$C1b_tertile == "low"] > 0))
}

tams <- rbindlist(tam_rows); sec <- rbindlist(sec_rows); tst <- rbindlist(tstate_rows); comp <- rbindlist(comp_rows)
fwrite(tams, file.path(P$objects, "codex_tams.csv.gz"))
write_tab(sec, "codex_section_level_stats.csv")

# sensitivity: radius
sens <- tams[C1b_tertile != "mid", lapply(.SD, function(v) mean(v > 0)), by = .(sample, group, C1b_tertile),
             .SDcols = patterns("^nT_")]
write_tab(sens, "codex_Tcontact_by_radius.csv")

# ---- pooled estimates per grade group ---------------------------------------------------
meta_fit <- function(dt) {
  dt <- dt[is.finite(logOR_adj) & is.finite(se_adj) & se_adj > 0]
  if (nrow(dt) < 2) return(NULL)
  f <- if (uniqueN(dt$patient) < nrow(dt))
         rma.mv(logOR_adj, se_adj^2, random = ~ 1 | patient / sample, data = dt, method = "REML")
       else rma(logOR_adj, sei = se_adj, data = dt, method = "REML")
  data.table(k = nrow(dt), estimate = as.numeric(f$b), ci_lb = f$ci.lb, ci_ub = f$ci.ub, p = f$pval,
             n_sections_pperm_lt_0.05 = sum(dt$p_perm < 0.05 & dt$mean_diff > 0),
             median_frac_hi = median(dt$frac_hi), median_frac_lo = median(dt$frac_lo))
}
groups <- list(`IDHm all` = c("IDHm G2", "IDHm G3", "IDHm G4"), `IDHm G2-3` = c("IDHm G2", "IDHm G3"),
               `IDHm G2` = "IDHm G2", `IDHm G3` = "IDHm G3", `IDHm G4` = "IDHm G4", GBM = "GBM")
pm <- rbindlist(lapply(names(groups), function(gn) sec[group %in% groups[[gn]],
        meta_fit(.SD), by = .(contrast, target)][, grp := gn]), fill = TRUE)
pm[, q := p.adjust(p, "BH"), by = grp]
write_tab(pm, "codex_meta_by_grade.csv")

# T-cell state next to APC+ TAMs (per-section logistic, pooled)
ts <- tst[, {
  out <- list()
  for (y in c("PD1", "CD69")) {
    d <- data.table(y = as.integer(get(y)), a = as.integer(near_APCpos))
    if (sum(d$y) >= 5 && sum(d$a) >= 5 && sum(1 - d$a) >= 5) {
      co <- summary(glm(y ~ a, family = binomial, data = d))$coefficients["a", ]
      out[[y]] <- data.table(state = y, logOR = co[1], se = co[2], frac_near = mean(d$y[d$a == 1]), frac_far = mean(d$y[d$a == 0]),
                             n_near = sum(d$a))
    }
  }
  rbindlist(out)
}, by = .(sample, group, patient)]
write_tab(ts, "codex_Tstate_near_APC_TAM.csv")
pool_state <- function(dt, gn) {
  if (nrow(dt) < 2) return(NULL)
  f <- rma(yi = dt$logOR, sei = dt$se, method = "REML")
  data.table(grp = gn, k = nrow(dt), estimate = as.numeric(f$b), ci_lb = f$ci.lb, ci_ub = f$ci.ub, p = f$pval,
             median_frac_near = median(dt$frac_near), median_frac_far = median(dt$frac_far))
}
tsm <- rbindlist(lapply(names(groups), function(gn) rbindlist(lapply(c("PD1", "CD69"), function(st) {
  r <- pool_state(ts[group %in% groups[[gn]] & state == st & is.finite(se)], gn)
  if (!is.null(r)) r[, state := st]
  r
}))))
write_tab(tsm, "codex_Tstate_meta.csv")

# TAM subtype composition of C1b-high vs C1b-low TAMs
cp <- comp[C1b_tertile != "mid", .(N = sum(N)), by = .(group = fifelse(group == "GBM", "GBM", "IDHm"), C1b_tertile, subtype)]
cp[, frac := N / sum(N), by = .(group, C1b_tertile)]
write_tab(cp, "codex_C1b_TAM_subtype_composition.csv")

# ---- figures ---------------------------------------------------------------------------
lev <- c("IDHm all", "IDHm G2-3", "IDHm G2", "IDHm G3", "IDHm G4", "GBM")
fp <- pm[target %in% c("T", "CD4", "CD8")][, grp := factor(grp, rev(lev))]
g1 <- ggplot(fp, aes(estimate, grp, colour = contrast)) + geom_vline(xintercept = 0, linewidth = .3) +
  geom_pointrange(aes(xmin = ci_lb, xmax = ci_ub), position = position_dodge(width = .6), size = .25) +
  facet_wrap(~target, labeller = labeller(target = c(T = "Any T cell", CD4 = "CD4+ T", CD8 = "CD8+ T"))) +
  scale_colour_manual(values = c(C1b_high_vs_low = "#762A83", MHCII_high_vs_low = "#E08214"),
                      labels = c(C1b_high_vs_low = "C1b-protein high vs low TAM", MHCII_high_vs_low = "MHC-II high vs low TAM")) +
  labs(x = "Adjusted log OR: T cell within 27.5 um (RE meta)", y = NULL, colour = NULL) +
  theme_pub() + theme(legend.position = "top")
fd <- pm[target %in% c("d_vessel", "d_hypoxia")][, grp := factor(grp, rev(lev))]
g2 <- ggplot(fd, aes(estimate, grp, colour = contrast)) + geom_vline(xintercept = 0, linewidth = .3) +
  geom_pointrange(aes(xmin = ci_lb, xmax = ci_ub), position = position_dodge(width = .6), size = .25) +
  facet_wrap(~target, labeller = labeller(target = c(d_vessel = "Distance to vessel", d_hypoxia = "Distance to CA9+ cell"))) +
  scale_colour_manual(values = c(C1b_high_vs_low = "#762A83", MHCII_high_vs_low = "#E08214"), guide = "none") +
  labs(x = "Difference in log10 distance (high - low)", y = NULL) + theme_pub()
save_pdf(g1 / g2 + plot_layout(heights = c(1.2, 1)), "Fig10a_CODEX_TAM_Tcell_vessel.pdf", 9, 6.5)

rs <- melt(sens, id.vars = c("sample", "group", "C1b_tertile"))
rs[, radius := as.numeric(sub("nT_", "", variable))]
g3 <- ggplot(rs, aes(factor(radius), value, fill = C1b_tertile)) +
  geom_boxplot(outlier.size = .4) + facet_wrap(~group, nrow = 1) +
  scale_fill_manual(values = c(high = "#762A83", low = "grey75")) +
  labs(x = "Radius (um)", y = "Fraction of TAMs with >=1 T cell", fill = "C1b-protein") + theme_pub()
g4 <- ggplot(cp[subtype %in% c("InfMg", "GAM", "MacScav", "MacBorder") | group == "GBM"],
             aes(C1b_tertile, frac, fill = subtype)) + geom_col() + facet_wrap(~group) +
  labs(x = "C1b-protein tertile", y = "Fraction of TAMs", fill = "Author TAM subtype") + theme_pub()
save_pdf(g3 / g4, "Fig10b_CODEX_radius_subtype.pdf", 9, 6)
logf("12 done")
