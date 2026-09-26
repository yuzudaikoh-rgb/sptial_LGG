# 08 - Full Consensus1b vs 3-gene proxy (CD68 / CD14 / HLA-DPB1) and the
#      neighbourhood of C1b-high (TAM) niches
#
#  A. Can a 3-gene panel stand in for the 90-gene C1b TAM program?
#     - spot-level concordance, Moran's I, hotspot / top-decile overlap
#     - null 1: 1,000 random 3-gene subsets of C1b
#     - null 2: 1,000 expression-matched random genome-wide trios
#     - does the 3-gene score reproduce the C1b proximity profile (part B)?
#  B. Are C1b-high regions adjacent to malignant cells, T cells, endothelium,
#     pericytes, IFN / hypoxia / inflammatory programs and other TAM states?
#     - target score inside hotspots and in the 1-2 spot ring around them
#     - distance-decay from the nearest hotspot
#     - significance: Moran spectral randomisation (MSR) of the target
#       (preserves the target's own autocorrelation), 499 surrogates
source("R/utils/common.R")
source("R/utils/scoring.R")
suppressPackageStartupMessages({
  library(Seurat); library(spdep); library(adespatial)
  library(ggplot2); library(patchwork); library(ComplexHeatmap)
})
set.seed(CFG$seed)
N_MSR  <- 499
N_TRIO <- 1000
PANEL  <- c("CD68", "CD14", "HLA-DPB1")

# ---- target signatures ------------------------------------------------------
# Curated marker sets (canonical markers; TAM states follow Pombo Antunes et al.
# Nat Neurosci 2021, Abdelfattah et al. Nat Commun 2022, Miller et al. 2025
# glioma myeloid programs). Genes shared with C1b or the 3-gene panel are
# removed below to avoid circularity.
TARGETS <- list(
  T_cell        = c("CD3D", "CD3E", "CD3G", "CD2", "CD247", "TRAC", "CD8A", "CD8B", "GZMK",
                    "GZMA", "CCL5", "NKG7", "IL7R", "LCK", "CD96"),
  Endothelial   = c("VWF", "CLDN5", "CDH5", "ESAM", "FLT1", "KDR", "EGFL7", "EMCN", "PLVAP",
                    "TIE1", "ENG", "ERG", "PECAM1"),
  Pericyte      = c("RGS5", "PDGFRB", "KCNJ8", "ABCC9", "NOTCH3", "CSPG4", "MCAM", "HIGD1B",
                    "ACTA2", "MYL9", "TAGLN", "DES"),
  IFN_response  = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI44L", "MX1", "MX2", "OAS1",
                    "OAS2", "OAS3", "RSAD2", "IFI27", "STAT1", "IRF7", "XAF1", "HERC5", "CMPK2", "IFITM3"),
  Hypoxia       = c("VEGFA", "ADM", "NDRG1", "BNIP3", "CA9", "PGK1", "SLC2A1", "P4HA1", "ERO1A",
                    "ANKRD37", "EGLN3", "BHLHE40", "ANGPTL4", "LOX", "PDK1", "ENO2", "HILPDA"),
  Inflammation  = c("IL1B", "TNF", "CXCL8", "CCL3", "CCL4", "CCL2", "NFKBIA", "PTGS2", "IL6",
                    "CXCL2", "CXCL3", "IER3", "SOD2", "ICAM1", "NFKBIZ", "IL1RN", "TNFAIP3"),
  TAM_microglia = c("P2RY12", "TMEM119", "CX3CR1", "SALL1", "GPR34", "SELPLG", "OLFML3",
                    "P2RY13", "CSF1R", "SLC2A5", "SIGLEC8"),
  TAM_monocyte  = c("FCN1", "VCAN", "S100A8", "S100A9", "S100A12", "CCR2", "SELL", "EREG",
                    "CD36", "CLEC12A", "LILRA5", "CD14", "LYZ"),
  TAM_SPP1_lipid = c("SPP1", "APOE", "APOC1", "GPNMB", "LPL", "FABP5", "TREM2", "CTSB", "ACP5",
                     "CD9", "PLD3", "LGALS3"),
  TAM_C1Q       = c("C1QA", "C1QB", "C1QC", "SELENOP", "FOLR2", "LYVE1", "F13A1", "MAF", "CD163")
)
excl <- c(C1B_GENES, map_genes(C1B_GENES, C1B_GENES), PANEL, "MARCHF1")
TARGETS <- lapply(TARGETS, setdiff, excl)
write_tab(data.table(signature = names(TARGETS), n_genes = lengths(TARGETS),
                     genes = vapply(TARGETS, paste, "", collapse = ";")), "target_signatures_used.csv")

meta  <- available_samples()
only  <- Sys.getenv("C1B_SAMPLES")          # optional: comma-separated subset for testing
if (nzchar(only)) meta <- meta[sample_id %in% strsplit(only, ",")[[1]]]
spots <- readRDS(file.path(P$objects, "spot_scores.rds"))
hot   <- readRDS(file.path(P$objects, "c1b_hotspots.rds"))
cna_f <- file.path(P$objects, "cna_spot.rds")
cna   <- if (file.exists(cna_f)) readRDS(cna_f) else NULL

zmean <- function(lognorm, genes) {
  g <- intersect(genes, rownames(lognorm))
  x <- as.matrix(lognorm[g, , drop = FALSE])
  z <- t(scale(t(x))); z[is.na(z)] <- 0
  colMeans(z)
}
gi_hot <- function(x, lw_b, fdr = CFG$spatial$hotspot_fdr) {
  gz <- as.numeric(localG(x, lw_b, zero.policy = TRUE))
  gz > 0 & p.adjust(2 * pnorm(-abs(gz)), "BH") < fdr
}
dice <- function(a, b) 2 * sum(a & b) / max(sum(a) + sum(b), 1)
# hex-grid hop distance from every spot to the nearest spot in `H` (BFS)
hop_distance <- function(nb, H, max_d = 8) {
  d <- rep(NA_integer_, length(nb)); d[H] <- 0L; front <- H
  for (k in seq_len(max_d)) {
    nxt <- setdiff(unique(unlist(nb[front])), 0)
    nxt <- nxt[is.na(d[nxt])]
    if (!length(nxt)) break
    d[nxt] <- k; front <- nxt
  }
  d[is.na(d)] <- max_d + 1L
  d
}

tgenes_rows <- list(); conc_rows <- list(); prox_rows <- list(); decay_rows <- list(); trio_rows <- list()
for (i in seq_len(nrow(meta))) {
  s <- meta$sample_id[i]
  so <- readRDS(obj_path(s))
  counts  <- GetAssayData(so, assay = "Spatial", layer = "counts")
  lognorm <- GetAssayData(so, assay = "Spatial", layer = "data")
  d <- merge(spots[sample_id == s, .(barcode, array_row, array_col, C1b, MP_Mac)],
             hot[sample_id == s, .(barcode, hotspot)], by = "barcode", sort = FALSE)
  d <- d[match(colnames(so), barcode)]
  lw <- visium_listw(d); W <- listw_to_sparse(lw)
  lw_b <- nb2listw(include.self(lw$neighbours), style = "B", zero.policy = TRUE)
  me <- adespatial::scores.listw(lw, MEM.autocor = "all")

  # ---- A. 3-gene proxy -------------------------------------------------------
  panel <- intersect(PANEL, rownames(counts))
  d[, P3_z := zmean(lognorm, panel)]
  d[, P3_UCell := UCell::ScoreSignatures_UCell(counts, features = list(P3 = panel),
                                               maxRank = CFG$scoring$ucell_maxrank, name = "")[, "P3"]]
  c1b <- unname(na.omit(map_genes(C1B_GENES, rownames(counts))))
  hot_c1b <- d$hotspot == "hot"
  hot_p3  <- gi_hot(d$P3_z, lw_b)
  top_c1b <- d$C1b >= quantile(d$C1b, .9); top_p3 <- d$P3_z >= quantile(d$P3_z, .9)

  # null 1: random trios drawn from C1b; null 2: expression-matched genome trios
  set.seed(CFG$seed + i)
  trios_in <- replicate(N_TRIO, sample(setdiff(c1b, panel), 3), simplify = FALSE)
  rs <- matched_random_sets(lognorm, panel, n_sets = N_TRIO, seed = CFG$seed + i)
  trio_rho <- function(tr) vapply(tr, function(g) cor(zmean(lognorm, g), d$C1b, method = "spearman"), 0)
  rho_in <- trio_rho(trios_in); rho_gen <- trio_rho(rs)
  rho_p3 <- cor(d$P3_z, d$C1b, method = "spearman")
  trio_rows[[s]] <- data.table(sample_id = s, null = rep(c("C1b_subset", "genome_matched"), each = N_TRIO),
                               rho = c(rho_in, rho_gen))
  conc_rows[[s]] <- data.table(
    sample_id = s, panel_genes_found = length(panel),
    rho_P3z_C1b = rho_p3, rho_P3UCell_C1b = cor(d$P3_UCell, d$C1b, method = "spearman"),
    pct_C1b_trios_below = mean(rho_in < rho_p3), pct_genome_trios_below = mean(rho_gen < rho_p3),
    median_rho_C1b_trios = median(rho_in), median_rho_genome_trios = median(rho_gen),
    moran_C1b = moran_matrix(cbind(d$C1b), W), moran_P3 = moran_matrix(cbind(d$P3_z), W),
    n_hot_C1b = sum(hot_c1b), n_hot_P3 = sum(hot_p3),
    dice_hotspots = dice(hot_c1b, hot_p3), jaccard_top10 = sum(top_c1b & top_p3) / sum(top_c1b | top_p3))

  # ---- B. target scores --------------------------------------------------------
  tg <- ucell_scores(counts, TARGETS, maxRank = CFG$scoring$ucell_maxrank, ncores = CFG$n_cores,
                     min_genes = 3)
  tgenes_rows[[s]] <- data.table(sample_id = s, target = names(TARGETS),
                                 n_genes_detected = vapply(TARGETS, function(g) sum(g %in% rownames(counts)), 0L),
                                 n_genes_total = lengths(TARGETS))
  tg <- as.data.table(tg)
  if (!is.null(cna)) tg[, Malignant_CNA := cna[sample_id == s][match(d$barcode, barcode), CNAtot]]
  targets <- names(tg)[vapply(tg, function(x) !anyNA(x) && sd(x) > 0, TRUE)]
  # MSR surrogates of each target, computed once and shared by C1b and 3-gene
  surs <- lapply(setNames(targets, targets), function(t)
    scale(adespatial::msr(tg[[t]], me, nrepet = N_MSR, method = "pair")))

  for (src in c("C1b", "P3")) {
    H <- which(if (src == "C1b") hot_c1b else hot_p3)
    if (length(H) < 5) next
    hd <- hop_distance(lw$neighbours, H)
    ring <- hd %in% 1:2; inside <- hd == 0
    for (t in targets) {
      y <- tg[[t]]
      ys <- scale(y)[, 1]
      sur <- surs[[t]]
      stat <- function(v, idx) mean(v[idx])
      for (zone in c("inside", "ring")) {
        idx <- if (zone == "inside") inside else ring
        obs <- stat(ys, idx); nul <- colMeans(sur[idx, , drop = FALSE])
        prox_rows[[paste(s, src, t, zone)]] <- data.table(
          sample_id = s, source = src, target = t, zone = zone, n_spots = sum(idx),
          mean_z = obs, z_msr = (obs - mean(nul)) / sd(nul),
          p_enrich = (1 + sum(nul >= obs)) / (N_MSR + 1), p_deplete = (1 + sum(nul <= obs)) / (N_MSR + 1))
      }
      if (src == "C1b")
        decay_rows[[paste(s, t)]] <- data.table(sample_id = s, target = t,
                                                dist = pmin(hd, 6L), z = ys)[, .(mean_z = mean(z), n = .N), by = .(sample_id, target, dist)]
    }
  }
  logf("%s: 3-gene rho=%.2f (C1b-trio pct %.2f), hotspot Dice=%.2f", s, rho_p3,
       mean(rho_in < rho_p3), conc_rows[[s]]$dice_hotspots)
  rm(so, counts, lognorm, me, surs); gc(verbose = FALSE)
}

conc <- merge(rbindlist(conc_rows), meta[, .(sample_id, cohort, histology, who_grade, is_LGG)], by = "sample_id")
write_tab(conc, "C1b_vs_3gene_concordance.csv")
write_tab(rbindlist(tgenes_rows), "target_signature_gene_detection.csv")
trio <- rbindlist(trio_rows); fwrite(trio, file.path(P$objects, "trio_null.csv.gz"))
prox <- merge(rbindlist(prox_rows), meta[, .(sample_id, cohort, who_grade, is_LGG)], by = "sample_id")
prox[, `:=`(q_enrich = p.adjust(p_enrich, "BH"), q_deplete = p.adjust(p_deplete, "BH")), by = .(sample_id, source, zone)]
write_tab(prox, "niche_proximity_per_sample.csv")
decay <- merge(rbindlist(decay_rows), meta[, .(sample_id, cohort)], by = "sample_id")
write_tab(decay, "niche_distance_decay.csv")

# ---- cohort-level combination (Stouffer on MSR z) -------------------------------
cohorts <- list(IDHm = meta[cohort == "IDHm", sample_id], LGG = meta[is_LGG == TRUE, sample_id],
                GBM = meta[cohort == "GBM", sample_id])
pm <- rbindlist(lapply(names(cohorts), function(ch) {
  prox[sample_id %in% cohorts[[ch]] & is.finite(z_msr),
       .(cohort = ch, k = .N, mean_z_msr = mean(z_msr), stouffer_z = sum(z_msr) / sqrt(.N),
         n_enriched = sum(q_enrich < 0.05), n_depleted = sum(q_deplete < 0.05)),
       by = .(source, target, zone)]
}))
pm[, p := 2 * pnorm(-abs(stouffer_z))][, q := p.adjust(p, "BH"), by = .(cohort, source, zone)]
write_tab(pm, "niche_proximity_meta.csv")

# agreement of C1b vs 3-gene proximity profiles
ag <- dcast(prox[zone == "ring"], sample_id + target ~ source, value.var = "z_msr")
ag <- ag[!is.na(C1b) & !is.na(P3)]
cls <- function(z) cut(z, c(-Inf, -1.96, 1.96, Inf), labels = c("depleted", "ns", "enriched"))
agr <- ag[, .(profile_spearman = cor(C1b, P3, method = "spearman"),
              call_agreement = mean(cls(C1b) == cls(P3)), n_targets = .N), by = sample_id]
kap <- function(a, b) { a <- factor(a, c("depleted", "ns", "enriched")); b <- factor(b, levels(a))
  tb <- table(a, b); po <- sum(diag(tb)) / sum(tb); pe <- sum(rowSums(tb) * colSums(tb)) / sum(tb)^2
  (po - pe) / (1 - pe) }
agr <- rbind(agr, data.table(sample_id = "ALL", profile_spearman = cor(ag$C1b, ag$P3, method = "spearman"),
                             call_agreement = mean(cls(ag$C1b) == cls(ag$P3)), n_targets = nrow(ag)),
             fill = TRUE)
agr[, cohen_kappa := c(ag[, kap(cls(C1b), cls(P3)), by = sample_id]$V1, kap(cls(ag$C1b), cls(ag$P3)))]
write_tab(agr, "C1b_vs_3gene_proximity_agreement.csv")

# ---- Figures -------------------------------------------------------------------
ord <- meta[order(cohort == "GBM", who_grade, sample_id)]
col_c <- c(IDHm = "#D6604D", GBM = "#4393C3")

# Fig 6a: 3-gene rho vs trio nulls
tr <- merge(trio, meta[, .(sample_id, cohort)], by = "sample_id")
tr[, sample_id := factor(sample_id, ord$sample_id)]
g1 <- ggplot(tr, aes(sample_id, rho, fill = null)) +
  geom_violin(scale = "width", colour = NA, position = position_dodge(.8)) +
  geom_point(data = conc[, .(sample_id = factor(sample_id, ord$sample_id), rho = rho_P3z_C1b)],
             aes(sample_id, rho), inherit.aes = FALSE, shape = 23, fill = "black", size = 1.8) +
  scale_fill_manual(values = c(C1b_subset = "#FDB863", genome_matched = "grey75"),
                    labels = c("Random C1b trios", "Expression-matched genome trios")) +
  labs(x = NULL, y = "Spearman rho with full C1b", fill = NULL,
       title = "CD68/CD14/HLA-DPB1 (diamond) vs random 3-gene panels") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
cl <- melt(conc[, .(sample_id, cohort, `Moran I (C1b)` = moran_C1b, `Moran I (3-gene)` = moran_P3,
                    `Hotspot Dice` = dice_hotspots, `Top-10% Jaccard` = jaccard_top10)],
           id.vars = c("sample_id", "cohort"))
g2 <- ggplot(cl, aes(cohort, value, colour = cohort)) +
  geom_boxplot(outlier.shape = NA, width = .5) + geom_jitter(width = .12, size = 1) +
  facet_wrap(~variable, nrow = 1, scales = "free_y") + scale_colour_manual(values = col_c) +
  labs(x = NULL, y = NULL) + theme_pub() + theme(legend.position = "none")
save_pdf(g1 / g2 + plot_layout(heights = c(1.2, 1)), "Fig6a_C1b_vs_3gene_concordance.pdf", 8, 6)

# Fig 6b: example maps (IDHm) C1b vs 3-gene
mp <- lapply(ord[cohort == "IDHm", sample_id], function(s) {
  x <- spots[sample_id == s, .(barcode, array_row, array_col, C1b)]
  so <- readRDS(obj_path(s)); ln <- GetAssayData(so, assay = "Spatial", layer = "data")
  x[, P3 := zmean(ln, intersect(PANEL, rownames(ln)))[barcode]]
  lim1 <- quantile(x$C1b, c(.01, .99)); lim2 <- quantile(x$P3, c(.01, .99))
  spot_plot(x, "C1b", paste(toupper(s), "C1b (90 genes)"), limits = lim1, size = .35) /
    spot_plot(x, "P3", "CD68/CD14/HLA-DPB1", limits = lim2, size = .35)
})
save_pdf(wrap_plots(mp, nrow = 1), "Fig6b_C1b_vs_3gene_maps_IDHm.pdf", 11, 4.2)

# Fig 7a: proximity heatmap (ring zone), C1b vs 3-gene, per sample
for (src in c("C1b", "P3")) {
  M <- dcast(prox[zone == "ring" & source == src], target ~ sample_id, value.var = "z_msr")
  mm <- as.matrix(M[, -1]); rownames(mm) <- M$target
  cols <- intersect(ord$sample_id, colnames(mm)); mm <- mm[, cols, drop = FALSE]
  Q <- dcast(prox[zone == "ring" & source == src], target ~ sample_id,
             value.var = "q_enrich")
  Qd <- dcast(prox[zone == "ring" & source == src], target ~ sample_id, value.var = "q_deplete")
  qq <- pmin(as.matrix(Q[, -1])[, cols, drop = FALSE], as.matrix(Qd[, -1])[, cols, drop = FALSE])
  hm <- Heatmap(mm, name = "MSR z", col = circlize::colorRamp2(c(-4, 0, 4), c("#2166AC", "white", "#B2182B")),
                cluster_columns = FALSE, column_split = ord[match(cols, sample_id), cohort],
                column_title = sprintf("%s-high niche: target enrichment in 1-2 spot ring (MSR z)",
                                       ifelse(src == "C1b", "C1b", "3-gene")),
                row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if (!is.na(qq[i, j]) && qq[i, j] < 0.05) grid::grid.text("*", x, y, gp = grid::gpar(fontsize = 7))
                })
  f <- file.path(P$figures, sprintf("Fig7a_niche_proximity_%s.pdf", src))
  grDevices::cairo_pdf(f, 7.5, 3.6); draw(hm); dev.off()
  png(sub("pdf$", "png", f), 7.5, 3.6, units = "in", res = 200); draw(hm); dev.off()
}

# Fig 7b: cohort-level forest, C1b vs 3-gene side by side
fp <- pm[zone == "ring" & cohort %in% c("IDHm", "GBM")]
lev <- fp[cohort == "IDHm" & source == "C1b"][order(stouffer_z), target]
fp[, target := factor(target, lev)]
g3 <- ggplot(fp, aes(stouffer_z, target, colour = source, shape = source)) +
  geom_vline(xintercept = c(-1.96, 0, 1.96), linetype = c(3, 1, 3), linewidth = .3) +
  geom_point(size = 2, position = position_dodge(width = .5)) + facet_wrap(~cohort) +
  scale_colour_manual(values = c(C1b = "#762A83", P3 = "#1B7837"), labels = c("C1b (90 genes)", "CD68/CD14/HLA-DPB1")) +
  scale_shape_manual(values = c(C1b = 16, P3 = 17), labels = c("C1b (90 genes)", "CD68/CD14/HLA-DPB1")) +
  labs(x = "Stouffer z (MSR), target in 1-2 spot ring of hotspots", y = NULL, colour = NULL, shape = NULL) +
  theme_pub() + theme(legend.position = "top")
# Fig 7c: distance decay (IDHm), key targets
dk <- decay[cohort == "IDHm", .(mean_z = weighted.mean(mean_z, n)), by = .(target, dist)]
g4 <- ggplot(dk, aes(dist, mean_z, colour = target)) + geom_hline(yintercept = 0, linewidth = .3) +
  geom_line() + geom_point(size = .8) +
  scale_x_continuous(breaks = 0:6, labels = c("hot", 1:5, ">=6")) +
  labs(x = "Distance to nearest C1b hotspot (spots)", y = "Target score (z, within sample)",
       title = "IDHm: distance decay from C1b niches", colour = NULL) + theme_pub()
save_pdf(g3 | g4, "Fig7b_niche_proximity_meta_decay.pdf", 10, 4.2)
logf("08 done")
