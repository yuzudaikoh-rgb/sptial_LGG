# 09b - Figures for step 09 (APC-TAM analysis); reads tables written by 09
source("R/utils/common.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

meta  <- available_samples()
ord   <- meta[order(cohort == "GBM", who_grade, sample_id)]
col_c <- c(IDHm = "#D6604D", GBM = "#4393C3")
col_s <- c(C1b = "#762A83", APC = "#E08214", nonMHC = "#1B7837")
lab_s <- c(C1b = "C1b (full)", APC = "APC index (MHC-II | TAM)", nonMHC = "non-MHC-II C1b | TAM")

sp    <- readRDS(file.path(P$objects, "apc09_spots.rds"))
moran <- fread(file.path(P$tables, "apc_index_moran.csv"))
pm    <- fread(file.path(P$tables, "meta_apc_vs_nonMHC_proximity.csv"))
pb    <- fread(file.path(P$tables, "niche_pseudobulk_Tcell.csv"))
lr    <- fread(file.path(P$tables, "lr_colocalization.csv"))
lrm   <- fread(file.path(P$tables, "meta_lr_colocalization.csv"))
lra   <- fread(file.path(P$tables, "meta_lr_colocalization_TAMadjusted.csv"))

# ---- Fig 8a: APC index maps and hotspot overlap (IDHm) ---------------------------
ids <- ord[cohort == "IDHm", sample_id]
pl <- lapply(ids, function(s) {
  x <- sp[sample_id == s][, APC_z := as.numeric(scale(APC_index))]
  x[, overlap := fcase(hot_APC & hot_C1b, "both", hot_APC, "APC-high only",
                       hot_C1b, "C1b-high only", default = "neither")]
  m <- moran[sample_id == s]
  a <- spot_plot(x, "APC_z", sprintf("%s\nI=%.2f, z=%.1f", toupper(s), m$I_APC_index, m$z_APC_vs_random),
                 limits = c(-2.5, 2.5), size = 0.45) + labs(colour = "APC index (z)")
  b <- spot_plot(x, "overlap", sprintf("Dice=%.2f", m$dice_APC_C1b), discrete = TRUE, size = 0.45) +
    scale_colour_manual(values = c(both = "#000000", `APC-high only` = "#E08214",
                                   `C1b-high only` = "#762A83", neither = "grey88"),
                        breaks = c("APC-high only", "C1b-high only", "both", "neither"), name = "Gi* hotspots") +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
  a / b
})
save_pdf(wrap_plots(pl, nrow = 1) + plot_layout(guides = "collect"), "Fig8a_APC_index_maps_IDHm.pdf", 11, 4.4)

# ---- Fig 8b: proximity profile of APC vs nonMHC vs C1b hotspots -----------------
fp <- pm[zone == "inside" & cohort %in% c("IDHm", "GBM")]
lev <- fp[cohort == "IDHm" & source == "C1b"][order(stouffer_z), target]
fp[, target := factor(target, unique(c(lev, fp$target)))]
g <- ggplot(fp, aes(stouffer_z, target, colour = source)) +
  geom_vline(xintercept = 0, linewidth = .3) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 3, linewidth = .3) +
  geom_point(size = 1.8, position = position_dodge(width = .6)) + facet_wrap(~cohort) +
  scale_colour_manual(values = col_s, labels = lab_s) +
  labs(x = "Stouffer z (MSR), target inside hotspots", y = NULL, colour = NULL) +
  theme_pub() + theme(legend.position = "top")
save_pdf(g, "Fig8b_APC_vs_nonMHC_proximity.pdf", 7.5, 4.5)

# ---- Fig 8c: niche pseudobulk T cells ---------------------------------------------
x <- merge(pb[tset != "IFNG" & umi_set_total > 0], meta[, .(sample_id)], by = "sample_id")
x[, tset := factor(tset, c("T_all", "CD4_T", "CD8_T", "Treg"))]
x[, sig := p_enrich < 0.05]
g <- ggplot(x, aes(tset, log2_ratio, colour = cohort)) +
  geom_hline(yintercept = 0, linewidth = .3) +
  geom_boxplot(outlier.shape = NA, fill = NA, position = position_dodge(width = .8), width = .6) +
  geom_point(aes(size = umi_set_in_niche, shape = sig), position = position_jitterdodge(.15, dodge.width = .8), alpha = .8) +
  facet_wrap(~niche, labeller = labeller(niche = c(C1b = "C1b niche (hotspot + 2-spot ring)",
                                                   APC = "APC-high niche (hotspot + 2-spot ring)"))) +
  scale_colour_manual(values = col_c) + scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), name = "p<0.05") +
  scale_size_area(max_size = 3.5, name = "UMI in niche") +
  labs(x = NULL, y = "log2 (niche / shape-matched random regions)") + theme_pub()
save_pdf(g, "Fig8c_niche_pseudobulk_Tcells.pdf", 8, 3.8)

# ---- Fig 8d: ligand-receptor co-localisation -----------------------------------
y <- merge(lr[is.finite(z)], meta[, .(sample_id)], by = "sample_id")
y[, sample_id := factor(sample_id, ord$sample_id)]
g1 <- ggplot(y, aes(sample_id, pair, fill = z)) + geom_tile(colour = "white") +
  geom_text(data = y[q_pos < 0.05], label = "*", size = 3) +
  scale_fill_gradient2(low = "#2166AC", high = "#B2182B", limits = c(-5, 5), oob = scales::squish) +
  facet_grid(~cohort, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, fill = "MSR z") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
lrb <- rbind(lrm[, .(pair, cohort, k, stouffer_z, adj = "depth-adjusted")],
             lra[, .(pair, cohort, k, stouffer_z, adj = "depth + TAM-adjusted")])
g2 <- ggplot(lrb[cohort %in% c("IDHm", "GBM")], aes(stouffer_z, pair, colour = cohort, shape = adj)) +
  geom_vline(xintercept = 0, linewidth = .3) + geom_vline(xintercept = c(-1.96, 1.96), linetype = 3, linewidth = .3) +
  geom_point(aes(size = k), position = position_dodge(width = .5)) + scale_colour_manual(values = col_c) +
  scale_shape_manual(values = c(`depth-adjusted` = 16, `depth + TAM-adjusted` = 1), name = NULL) +
  scale_size_area(max_size = 3, name = "n samples") +
  labs(x = "Stouffer z (MSR)", y = NULL) + theme_pub()
save_pdf(g1 + g2 + plot_layout(widths = c(2.6, 1.4)), "Fig8d_ligand_receptor_colocalisation.pdf", 12, 3.8)
logf("09b done")
