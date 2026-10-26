# 14 - Figure 4 (spatial-transcriptomics-only version)
#  Same result tables as step 13, but keeps only the Visium conclusions that
#  replicate in both cohorts: CODEX panels, the protein-panel proxy and the
#  T-cell rows are left out. Cell Press layout: 174 mm wide, Arial, >= 6 pt
#  text, upper-case panel tags. Writes to manuscript/figure4_ST/.
source("R/utils/common.R")
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))  # non-ASCII labels (\u2265, \u00b1) are dropped in a C locale
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(viridis) })

OUT <- file.path(ROOT, "manuscript", "figure4_ST"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
TD <- file.path(ROOT, "results", "tables"); TR <- file.path(ROOT, "results", "replication", "tables")
OD <- file.path(ROOT, "data", "objects"); OR <- file.path(ROOT, "data", "objects_replication")
rd <- function(dir, f) fread(file.path(dir, f))
C1 <- "Cohort 1 (GSE237183)"; C2 <- "Cohort 2 (Hoefflin)"
COH <- setNames(c("#E08214", "#542788"), c(C1, C2))
FAM <- "Arial"   # falls back to the metric-compatible Liberation Sans where Arial is not installed
FS <- 7; FT <- 6 # body text / smallest text (ticks, legends, annotations)
update_geom_defaults("text", list(family = FAM, size = FT / .pt))
th <- theme_pub(FS) + theme(text = element_text(family = FAM),
                            plot.title = element_text(size = FS, face = "bold"),
                            axis.text = element_text(size = FT, colour = "black"), strip.text = element_text(size = FT, face = "bold"),
                            legend.text = element_text(size = FT), legend.title = element_text(size = FT),
                            plot.tag = element_text(size = 9, face = "bold"))
mm <- function(x) x / 25.4
savep <- function(p, name, w, h) {
  ggsave(file.path(OUT, name), p, width = mm(w), height = mm(h), units = "in", device = grDevices::cairo_pdf)
}

# ---- A: representative sections --------------------------------------------------
sec_map <- function(objdir, s, title) {
  sp <- readRDS(file.path(objdir, "spot_scores.rds"))[sample_id == s]
  ho <- readRDS(file.path(objdir, "c1b_hotspots.rds"))[sample_id == s]
  x <- merge(sp[, .(barcode, array_row, array_col, C1b)], ho[, .(barcode, hotspot)], by = "barcode")
  x[, C1b_z := as.numeric(scale(C1b))]
  base <- function(v) ggplot(x, aes(array_col, -array_row, colour = .data[[v]])) +
    geom_point(size = 0.26, shape = 16) + coord_fixed(sqrt(3)) + theme_void(base_size = FS, base_family = FAM) +
    theme(plot.title = element_text(size = FT, hjust = 0.5, margin = margin(b = 2)), plot.margin = margin(4, 2, 2, 8),
          legend.key.size = grid::unit(2.5, "mm"),
          legend.text = element_text(size = FT), legend.title = element_text(size = FT))
  p1 <- base("C1b_z") + scale_colour_viridis(option = "magma", limits = c(-2.5, 2.5), oob = scales::squish,
                                             name = "C1b (z)") + ggtitle(title)
  p2 <- base("hotspot") + scale_colour_manual(values = c(hot = "#B2182B", cold = "#2166AC", ns = "grey88"),
                                              breaks = c("hot", "cold"), labels = c("Hotspot", "Coldspot"),
                                              name = "Gi* (FDR<0.05)") +
    ggtitle(" ") + guides(colour = guide_legend(override.aes = list(size = 1.8)))
  p1 | p2
}
pa <- (sec_map(OD, "mgh259", "MGH259, oligodendroglioma G3\ncohort 1") |
       sec_map(OR, "idha17", "IDHA17, astrocytoma G2\ncohort 2")) + plot_layout(guides = "collect")

# ---- B: coherence vs matched random gene sets --------------------------------------
mb <- rbind(rd(TD, "moran_I_C1b_vs_random_sets.csv")[cohort == "IDHm"][, cohort2 := C1],
            rd(TR, "moran_I_C1b_vs_random_sets.csv")[, cohort2 := C2], fill = TRUE)
mb[, feature := factor(feature, c("C1b", "C1b_resid_macadj"), c("C1b", "C1b | TAM"))]
st <- rbind(rd(TD, "moran_vs_random_stouffer.csv")[cohort == "IDHm"][, cohort2 := C1],
            rd(TR, "moran_vs_random_stouffer.csv")[cohort == "IDHm"][, cohort2 := C2])
st[, feature := factor(feature, c("C1b", "C1b_resid_macadj"), c("C1b", "C1b | TAM"))]
st[, lab := sprintf("Z=%.1f", stouffer_z)]
st[, ylab := fifelse(cohort2 == C1, 11.2, 10.2)]
mb[, cohort2 := factor(cohort2, c(C1, C2))]; st[, cohort2 := factor(cohort2, c(C1, C2))]
pb <- ggplot(mb, aes(feature, z_vs_random, colour = cohort2)) +
  geom_hline(yintercept = 0, linewidth = .3) + geom_hline(yintercept = 1.645, linetype = 3, linewidth = .3) +
  geom_point(position = position_jitterdodge(jitter.width = .15, dodge.width = .7), size = .9, alpha = .85) +
  geom_text(data = st, aes(feature, ylab, label = lab), show.legend = FALSE) +
  scale_colour_manual(values = COH, name = NULL) + coord_cartesian(ylim = c(-6.5, 11.5)) +
  labs(x = NULL, y = "Moran's I vs 1,000 matched\nrandom gene sets (z)") + th +
  theme(legend.position = "top", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0))

# ---- C: specific co-localisation (random-gene-set baseline), IDHm only -------------
sc <- rbind(rd(TD, "coloc_specificity_vs_random_sets.csv")[cohort == "IDHm"][, col := "Cohort 1"],
            rd(TR, "coloc_specificity_vs_random_sets.csv")[cohort == "IDHm"][, col := "Cohort 2"])
sc[, signed := (n_above95 - n_below5) / k]
sc[, size := pmax(n_above95, n_below5) / k]
ordr <- sc[col == "Cohort 2"][order(signed), MP]
sc[, MP := factor(gsub("\\.", "-", MP), gsub("\\.", "-", ordr))]
pc <- ggplot(sc, aes(col, MP, fill = signed, size = size)) + geom_point(shape = 21, colour = "grey30", stroke = .2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1),
                       name = "Specific\nco-loc. -\navoidance") +
  scale_size_area(max_size = 3.6, limits = c(0, 1), name = "Fraction of\nsections") +
  labs(x = NULL, y = NULL) + th + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- D: niche composition inside C1b hotspots (MSR null, depth-adjusted) -----------
tgt <- c("Endothelial", "Pericyte", "TAM_C1Q", "TAM_microglia", "TAM_SPP1_lipid", "TAM_monocyte",
         "Inflammation", "IFN_response", "Hypoxia")
labs_t <- c(Endothelial = "Endothelium", Pericyte = "Pericytes", TAM_C1Q = "C1Q+ TAM", TAM_microglia = "Microglia-like TAM",
            TAM_SPP1_lipid = "SPP1+/lipid TAM", TAM_monocyte = "Monocyte-like TAM", Inflammation = "Inflammation",
            IFN_response = "Type I IFN", IFNG_response = "IFN-gamma", Hypoxia = "Hypoxia")
nd <- rbind(rd(TD, "niche_proximity_meta.csv")[cohort == "IDHm"][, cohort2 := C1],
            rd(TR, "niche_proximity_meta.csv")[cohort == "IDHm"][, cohort2 := C2])
nd <- nd[source == "C1b" & zone == "inside" & target %in% tgt, .(cohort2, target, stouffer_z, q, k)]
na <- rbind(rd(TD, "meta_apc_vs_nonMHC_proximity.csv")[cohort == "IDHm"][, cohort2 := C1],
            rd(TR, "meta_apc_vs_nonMHC_proximity.csv")[cohort == "IDHm"][, cohort2 := C2])
nd <- rbind(nd, na[source == "C1b" & zone == "inside" & target == "IFNG_response", .(cohort2, target, stouffer_z, q, k)])
ord_d <- nd[cohort2 == C2][order(stouffer_z), target]
nd[, target := factor(target, ord_d, labs_t[ord_d])]
nd[, cohort2 := factor(cohort2, c(C1, C2))]
nd[, sig := factor(q < 0.05, c(TRUE, FALSE), c("q < 0.05", "n.s."))]
pd <- ggplot(nd, aes(stouffer_z, target, colour = cohort2, shape = sig)) +
  geom_vline(xintercept = 0, linewidth = .3) + geom_vline(xintercept = c(-1.96, 1.96), linetype = 3, linewidth = .3) +
  geom_point(size = 1.6, position = position_dodge(width = .6)) +
  scale_colour_manual(values = COH, guide = "none") +
  scale_shape_manual(values = c(`q < 0.05` = 16, n.s. = 1), name = NULL) +
  labs(x = "Enrichment inside C1b hotspots\n(Stouffer z; MSR null, depth-adjusted)", y = NULL) + th +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))

# ---- E: distance decay from the nearest C1b hotspot --------------------------------
tg_e <- c(Endothelial = "Endothelium", Pericyte = "Pericytes", TAM_C1Q = "C1Q+ TAM", Hypoxia = "Hypoxia")
de <- rbind(rd(TD, "niche_distance_decay.csv")[cohort == "IDHm"][, cohort2 := C1],
            rd(TR, "niche_distance_decay.csv")[cohort == "IDHm"][, cohort2 := C2])[target %in% names(tg_e)]
de <- de[, .(m = mean(mean_z), se = sd(mean_z) / sqrt(.N), n = .N), by = .(cohort2, target, dist)]
de[, target := factor(target, names(tg_e), tg_e)]
de[, cohort2 := factor(cohort2, c(C1, C2))]
pe <- ggplot(de, aes(dist, m, colour = cohort2, fill = cohort2)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
  geom_ribbon(aes(ymin = m - se, ymax = m + se), colour = NA, alpha = .2) +
  geom_line(linewidth = .5) + geom_point(size = .8) + facet_wrap(~target, nrow = 1) +
  scale_x_continuous(breaks = 0:6, labels = c("0", "1", "2", "3", "4", "5", "≥6")) +
  scale_colour_manual(values = COH, guide = "none") + scale_fill_manual(values = COH, guide = "none") +
  labs(x = "Distance to nearest C1b hotspot (spots; 0 = inside)", y = "Target score (z, mean ± SE)") + th

# ---- F: antigen-presenting vs mesenchymal-like arm ----------------------------------
ah <- rbind(rd(TR, "meta_apc_vs_nonMHC_proximity.csv")[cohort == "IDHm"][, cohort2 := "Cohort 2 IDHm"],
            rd(TD, "meta_apc_vs_nonMHC_proximity.csv")[cohort == "GBM"][, cohort2 := "Cohort 1 GBM"])
tg_f <- c(Pericyte = "Pericytes", Endothelial = "Endothelium", IFNG_response = "IFN-gamma", IFN_response = "Type I IFN",
          Inflammation = "Inflammation", Hypoxia = "Hypoxia")
ah <- ah[zone == "inside" & target %in% names(tg_f)]
ah[, target := factor(target, rev(names(tg_f)), rev(tg_f))]
ah[, source := factor(source, c("C1b", "APC", "nonMHC"), c("C1b (full)", "MHC-II arm | TAM", "Mesenchymal arm | TAM"))]
ah[, cohort2 := factor(cohort2, c("Cohort 2 IDHm", "Cohort 1 GBM"))]
pf <- ggplot(ah, aes(stouffer_z, target, colour = source)) +
  geom_vline(xintercept = 0, linewidth = .3) + geom_vline(xintercept = c(-1.96, 1.96), linetype = 3, linewidth = .3) +
  geom_point(size = 1.4, position = position_dodge(width = .65)) + facet_wrap(~cohort2, nrow = 1) +
  scale_colour_manual(values = c("#762A83", "#E08214", "#1B7837"), name = NULL) +
  labs(x = "Enrichment inside hotspots (Stouffer z)", y = NULL) + th +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0)) + guides(colour = guide_legend(nrow = 2))

# ---- assemble ---------------------------------------------------------------------------
W <- function(p) wrap_elements(full = p)
fig <- W(pa) / (W(pb) + W(pc) + W(pd) + plot_layout(widths = c(0.85, 0.95, 1.3))) /
       (W(pe) + W(pf) + plot_layout(widths = c(1.35, 1))) +
  plot_layout(heights = c(0.85, 1.25, 0.95)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 9, face = "bold", family = FAM))
savep(fig, "Figure4.pdf", 174, 190)
ggsave(file.path(OUT, "Figure4_preview.png"), fig, width = mm(174), height = mm(190), units = "in", dpi = 250)

# ---- individual panels ----------------------------------------------------------------------
savep(pa, "Figure4A.pdf", 174, 48); savep(pb, "Figure4B.pdf", 50, 68); savep(pc, "Figure4C.pdf", 55, 68)
savep(pd, "Figure4D.pdf", 70, 68); savep(pe, "Figure4E.pdf", 100, 55); savep(pf, "Figure4F.pdf", 75, 55)
logf("ST-only Figure 4 written to %s", OUT)
