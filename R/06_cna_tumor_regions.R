# 06 - CNA inference and tumour / non-tumour regions ------------------------
#  * CNA per spot (Greenwald et al. approach; external normal-brain Visium
#    references UKF256_C / UKF265_C + iterative in-sample reference refinement)
#  * malignancy calls against a per-sample genomic-order permutation null
#  * positive controls: 1p/19q co-deletion (oligodendroglioma), +7/-10 (GBM)
#  * validation against the authors' GBM malignancy levels
#  * C1b by region (linear mixed models), C1b vs CNA, hotspot-region overlap
source("R/utils/common.R")
source("R/utils/cna.R")
suppressPackageStartupMessages({
  library(Seurat); library(lme4); library(lmerTest); library(SpatialPack)
  library(ggplot2); library(patchwork); library(ComplexHeatmap)
})
set.seed(CFG$seed)
ccfg <- CFG$cna
N_PERM <- 10

meta  <- available_samples()
spots <- readRDS(file.path(P$objects, "spot_scores.rds"))
hot   <- readRDS(file.path(P$objects, "c1b_hotspots.rds"))
genome <- readRDS(file.path(P$inputs_dir, "CNA", "hg38.rds"))
ref_files <- list.files(file.path(P$inputs_dir, "CNA", "normal_brain_ref"), full.names = TRUE)
refs <- lapply(ref_files, function(f) { r <- as(as.matrix(readRDS(f)), "CsparseMatrix"); r })
names(refs) <- sub("\\.rds$", "", basename(ref_files))
for (k in names(refs)) colnames(refs[[k]]) <- paste0(k, "_", colnames(refs[[k]]))

cna_rows <- list(); arm_rows <- list()
done_f <- file.path(P$objects, "cna_spot.rds")
resume <- file.exists(done_f) && Sys.getenv("FORCE_CNA") != "1" &&
  all(file.exists(obj_path(meta$sample_id, "cna_matrix")))
if (resume) logf("CNA per-sample results found; reusing (set FORCE_CNA=1 to recompute)")
for (s in if (resume) character() else meta$sample_id) {
  so <- readRDS(obj_path(s))
  cnt <- GetAssayData(so, assay = "Spatial", layer = "counts")
  genes <- Reduce(intersect, c(list(rownames(cnt)), lapply(refs, rownames)))
  M <- do.call(cbind, c(list(cnt[genes, ]), lapply(refs, function(r) r[genes, ])))
  m <- cna_logcpm(M, ccfg$min_mean)
  query <- colnames(cnt)
  ref0 <- lapply(refs, colnames)

  # pass 1: external reference only
  c1 <- calc_cna(m, query, ref0, genome, ccfg$window, unlist(ccfg$range), ccfg$noise)
  sc1 <- cna_sig_cor(c1, query)
  # refine reference with in-sample low-CNA spots (author thresholds)
  low <- query[sc1$cor[query] <= ccfg$ref_cor_cut & sc1$sig[query] <= ccfg$ref_sig_cut]
  ref1 <- c(ref0, if (length(low) >= 10) list(insample = low))
  c2 <- calc_cna(m, query, ref1, genome, ccfg$window, unlist(ccfg$range), ccfg$noise)
  sc2 <- cna_sig_cor(c2, query)
  tot <- sc2$cor + sc2$sig * max(sc2$cor[query], na.rm = TRUE) / max(sc2$sig[query], na.rm = TRUE)

  # genomic-order permutation null (same references, gene order shuffled)
  null <- unlist(lapply(seq_len(N_PERM), function(k) {
    cp <- calc_cna(m, query, ref1, genome, ccfg$window, unlist(ccfg$range), ccfg$noise,
                   permute = TRUE, seed = CFG$seed + k)
    sp <- cna_sig_cor(cp, query)
    (sp$cor + sp$sig * max(sc2$cor[query]) / max(sc2$sig[query]))[query]
  }))
  thr_hi <- quantile(null, 0.99); thr_lo <- quantile(null, 0.95)
  refspots <- unlist(ref0)
  cls <- ifelse(tot[query] > thr_hi, "malignant", ifelse(tot[query] <= thr_lo, "non-malignant", "intermediate"))

  cna_rows[[s]] <- data.table(sample_id = s, barcode = query, CNAsig = sc2$sig[query],
                              CNAcor = sc2$cor[query], CNAtot = tot[query], cna_class = cls,
                              thr_hi = thr_hi, thr_lo = thr_lo,
                              ref_fpr = mean(tot[refspots] > thr_hi))
  # arm-level CNA (malignant spots; all spots if < 20 malignant)
  cq <- c2[, query]                       # subsetting drops attributes; restore them
  attr(cq, "chr") <- attr(c2, "chr"); attr(cq, "arm") <- attr(c2, "arm")
  am <- arm_means(cq)
  use <- if (sum(cls == "malignant") >= 20) query[cls == "malignant"] else query
  arm_rows[[s]] <- data.table(sample_id = s, arm = rownames(am), mean_cna = rowMeans(am[, use, drop = FALSE]),
                              n_spots_used = length(use))
  saveRDS(cq, obj_path(s, "cna_matrix"))
  logf("%s: malignant %.0f%%, intermediate %.0f%%, non-malignant %.0f%% (ref FPR %.3f)", s,
       100 * mean(cls == "malignant"), 100 * mean(cls == "intermediate"),
       100 * mean(cls == "non-malignant"), mean(tot[refspots] > thr_hi))
  rm(so, cnt, M, m, c1, c2, cq); gc(verbose = FALSE)
}

if (resume) {
  cna <- readRDS(done_f)
  arm <- fread(file.path(P$tables, "cna_arm_means_malignant_spots.csv"))
} else {
  cna <- rbindlist(cna_rows)
  arm <- merge(rbindlist(arm_rows), meta[, .(sample_id, cohort, histology, who_grade)], by = "sample_id")
  saveRDS(cna, done_f)
}
write_tab(cna[, .(n = .N, malignant = mean(cna_class == "malignant"),
                  intermediate = mean(cna_class == "intermediate"),
                  non_malignant = mean(cna_class == "non-malignant"),
                  thr_hi = thr_hi[1], ref_false_positive_rate = ref_fpr[1]), by = sample_id],
          "cna_class_summary.csv")
write_tab(arm, "cna_arm_means_malignant_spots.csv")

# ---- positive controls -----------------------------------------------------
pc <- dcast(arm[arm %in% c("1p", "19q", "7p", "7q", "10p", "10q")], sample_id + cohort + histology + who_grade ~ arm,
            value.var = "mean_cna")
# call an arm event when the mean CNA of CNA-positive spots passes +-0.05
pc[, codel_1p19q := `1p` < -0.05 & `19q` < -0.05]
pc[, gain7_loss10 := (`7p` + `7q`) / 2 > 0.05 & (`10p` + `10q`) / 2 < -0.05]
write_tab(pc, "cna_positive_controls_1p19q_chr7_10.csv")

# ---- validation against author GBM malignancy levels -----------------------
vm <- fread(file.path(P$inputs_dir, "general", "visium_metadata.csv"))
mal_lev <- readRDS(file.path(P$inputs_dir, "CNA", "mal_lev.rds"))
val <- rbindlist(lapply(meta[!is.na(author_id), sample_id], function(s) {
  a <- meta[sample_id == s, author_id]
  x <- copy(cna[sample_id == s])
  ml <- mal_lev[[a]]
  bins <- vm[sample == a]
  x$author_mal_lev <- unname(ml[x$barcode])
  x$author_bin <- bins$cna_bin[match(x$barcode, bins$spot_id)]
  data.table(sample_id = s, n = sum(!is.na(x$author_mal_lev)),
             spearman_CNAtot_vs_author = cor(x$CNAtot, x$author_mal_lev, method = "spearman", use = "complete.obs"),
             agreement_malignant_vs_nonmal = {
               y <- x[!is.na(author_bin) & author_bin %in% c("malignant", "non_malignant") & cna_class != "intermediate"]
               mean((y$author_bin == "malignant") == (y$cna_class == "malignant")) })
}))
write_tab(val, "cna_validation_vs_authors_GBM.csv")

# ---- C1b by region -----------------------------------------------------------
d <- merge(spots[, .(sample_id, barcode, C1b, C1b_noMPoverlap, MP_Mac, MP_Inflammatory.Mac, array_row, array_col)],
           cna, by = c("sample_id", "barcode"))
d <- merge(d, hot[, .(sample_id, barcode, hotspot)], by = c("sample_id", "barcode"))
d <- merge(d, meta[, .(sample_id, cohort, who_grade, is_LGG, patient)], by = "sample_id")
d[, cna_class := factor(cna_class, c("non-malignant", "intermediate", "malignant"))]
d[, C1b_z := as.numeric(scale(C1b)), by = sample_id]
saveRDS(d, file.path(P$objects, "spot_integrated.rds"))

fit_lmm <- function(dd, label) {
  dd <- droplevels(dd)
  if (nlevels(dd$cna_class) < 2) return(NULL)
  re <- if (uniqueN(dd$patient) < uniqueN(dd$sample_id)) "(1 | patient/sample_id)" else "(1 | sample_id)"
  f <- lmer(as.formula(paste("C1b_z ~ cna_class +", re)), data = dd, REML = TRUE)
  co <- as.data.table(summary(f)$coefficients, keep.rownames = "term")
  co[, cohort := label]
  f2 <- lmer(as.formula(paste("C1b_z ~ cna_class + scale(MP_Mac) + scale(MP_Inflammatory.Mac) +", re)), data = dd)
  co2 <- as.data.table(summary(f2)$coefficients, keep.rownames = "term"); co2[, cohort := paste0(label, "_macadj")]
  rbind(co, co2)
}
lmm <- rbindlist(list(fit_lmm(d[cohort == "IDHm"], "IDHm"), fit_lmm(d[is_LGG == TRUE], "LGG"),
                      fit_lmm(d[cohort == "GBM"], "GBM")))
write_tab(lmm, "lmm_C1b_by_cna_class.csv")

# per-sample C1b vs CNAtot (Dutilleul modified t-test)
cc <- rbindlist(lapply(unique(d$sample_id), function(s) {
  x <- d[sample_id == s]
  mt <- modified.ttest(rank(x$C1b), rank(x$CNAtot), cbind(x$array_col, x$array_row * sqrt(3)), nclass = 13)
  data.table(sample_id = s, rho = cor(x$C1b, x$CNAtot, method = "spearman"), ess = mt$ESS, p = mt$p.value)
}))
cc <- merge(cc, meta[, .(sample_id, cohort, who_grade)], by = "sample_id")
cc[, q := p.adjust(p, "BH")]
write_tab(cc, "C1b_vs_CNAtot_per_sample.csv")

# hotspot vs region overlap
ov <- d[, {
  tb <- table(factor(hotspot == "hot", c(TRUE, FALSE)), factor(cna_class == "malignant", c(TRUE, FALSE)))
  ft <- fisher.test(tb)
  .(frac_hot_in_malignant = tb[1, 1] / sum(tb[1, ]), frac_rest_in_malignant = tb[2, 1] / sum(tb[2, ]),
    OR_malignant = unname(ft$estimate), p = ft$p.value, n_hot = sum(tb[1, ]))
}, by = .(sample_id, cohort)]
ov[, q := p.adjust(p, "BH")]
write_tab(ov, "hotspot_vs_malignant_region_fisher.csv")

# ---- Figures ---------------------------------------------------------------
ord <- meta[order(cohort == "GBM", who_grade, sample_id)]
# CNA heatmap per sample (random 300 spots each, ordered by CNAtot)
set.seed(CFG$seed)
mats <- lapply(ord$sample_id, function(s) {
  cm <- readRDS(obj_path(s, "cna_matrix"))
  x <- cna[sample_id == s]
  sp <- x[sample(.N, min(300, .N))][order(CNAtot), barcode]
  list(m = cm[, sp], chr = attr(cm, "chr"), s = s)
})
common <- Reduce(intersect, lapply(mats, function(z) rownames(z$m)))
H <- do.call(cbind, lapply(mats, function(z) z$m[common, ]))
chr <- mats[[1]]$chr[match(common, rownames(mats[[1]]$m))]
col_s <- rep(ord$sample_id, sapply(mats, function(z) ncol(z$m)))
cls_v <- cna[match(paste(col_s, colnames(H)), paste(sample_id, barcode)), cna_class]
hm <- Heatmap(t(H), name = "CNA", cluster_rows = FALSE, cluster_columns = FALSE,
              col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("#2166AC", "white", "#B2182B")),
              column_split = factor(chr, levels = 1:22), row_split = factor(col_s, ord$sample_id),
              show_row_names = FALSE, show_column_names = FALSE, row_title_rot = 0,
              row_title_gp = grid::gpar(fontsize = 6), column_title_gp = grid::gpar(fontsize = 6),
              left_annotation = rowAnnotation(Class = cls_v,
                col = list(Class = c(malignant = "#B2182B", intermediate = "#F4A582", `non-malignant` = "#4D4D4D")),
                simple_anno_size = grid::unit(2, "mm")),
              use_raster = TRUE)
f <- file.path(P$figures, "Fig5a_CNA_heatmap.pdf")
grDevices::cairo_pdf(f, 9, 10); draw(hm); dev.off()
png(sub("pdf$", "png", f), 9, 10, units = "in", res = 150); draw(hm); dev.off()

# arm-level positive-control plot
pa <- arm[arm %in% c("1p", "19q", "7p", "7q", "10p", "10q")]
pa[, sample_id := factor(sample_id, ord$sample_id)]
g1 <- ggplot(pa, aes(sample_id, factor(arm, c("1p", "19q", "7p", "7q", "10p", "10q")), fill = mean_cna)) +
  geom_tile() + scale_fill_gradient2(low = "#2166AC", high = "#B2182B", limits = c(-0.3, 0.3), oob = scales::squish) +
  labs(x = NULL, y = "Arm", fill = "Mean CNA") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# C1b by region
g2 <- ggplot(d[, .(C1b_z = mean(C1b_z)), by = .(sample_id, cohort, cna_class)],
             aes(cna_class, C1b_z)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_line(aes(group = sample_id, colour = cohort), alpha = .5) + geom_point(aes(colour = cohort), size = 1) +
  facet_wrap(~cohort) + scale_colour_manual(values = c(IDHm = "#D6604D", GBM = "#4393C3")) +
  labs(x = NULL, y = "Mean C1b (z within sample)") + theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
maps <- lapply(ord[cohort == "IDHm", sample_id], function(s) {
  spot_plot(d[sample_id == s], "cna_class", toupper(s), discrete = TRUE, size = 0.4) +
    scale_colour_manual(values = c(malignant = "#B2182B", intermediate = "#F4A582", `non-malignant` = "#4D4D4D"))
})
save_pdf((g1 / g2) | (wrap_plots(maps, ncol = 2) + plot_layout(guides = "collect") &
                      guides(colour = guide_legend(override.aes = list(size = 2.5)))),
         "Fig5b_CNA_regions_C1b.pdf", 10, 6.5)
