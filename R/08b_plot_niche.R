# 08b - Figures for step 08 (reads the tables written by R/08_c1b_vs_3gene_niche.R)
source("R/utils/common.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork); library(ComplexHeatmap) })
PANEL <- c("CD68", "CD14", "HLA-DPB1")
zmean <- function(lognorm, genes) {
  g <- intersect(genes, rownames(lognorm))
  x <- as.matrix(lognorm[g, , drop = FALSE])
  z <- t(scale(t(x))); z[is.na(z)] <- 0
  colMeans(z)
}
meta  <- available_samples()
spots <- readRDS(file.path(P$objects, "spot_scores.rds"))
conc  <- fread(file.path(P$tables, "C1b_vs_3gene_concordance.csv"))
prox  <- fread(file.path(P$tables, "niche_proximity_per_sample.csv"))
decay <- fread(file.path(P$tables, "niche_distance_decay.csv"))
pm    <- fread(file.path(P$tables, "niche_proximity_meta.csv"))
trio  <- fread(file.path(P$objects, "trio_null.csv.gz"))

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
  geom_vline(xintercept = 0, linewidth = .3) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 3, linewidth = .3) +
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
logf("08b done")
