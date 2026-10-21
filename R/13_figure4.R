# 13 - Figure 4: spatial organisation of the C1b antigen-presenting TAM program
#  Re-draws all panels from the result tables of steps 04-12 (discovery,
#  replication, CODEX) in one style; writes the assembled figure and every
#  panel as a separate PDF to manuscript/figure4/, plus renamed supplementary
#  PDFs to manuscript/figure4/supplementary/.
source("R/utils/common.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(arrow); library(viridis) })

OUT <- file.path(ROOT, "manuscript", "figure4"); dir.create(file.path(OUT, "supplementary"), recursive = TRUE, showWarnings = FALSE)
TD <- file.path(ROOT, "results", "tables"); TR <- file.path(ROOT, "results", "replication", "tables")
OD <- file.path(ROOT, "data", "objects"); OR <- file.path(ROOT, "data", "objects_replication")
rd <- function(dir, f) fread(file.path(dir, f))
COH <- c(Discovery = "#E08214", Replication = "#542788", `Discovery GBM` = "#4393C3")
FS <- 7
th <- theme_pub(FS) + theme(plot.title = element_text(size = FS, face = "bold"),
                            legend.text = element_text(size = FS - 1), legend.title = element_text(size = FS - 1),
                            plot.tag = element_text(size = 9, face = "bold"))
mm <- function(x) x / 25.4
savep <- function(p, name, w, h) {
  ggsave(file.path(OUT, name), p, width = mm(w), height = mm(h), units = "in", device = grDevices::cairo_pdf)
}

# ---- a: representative sections --------------------------------------------------
sec_map <- function(objdir, s, title) {
  sp <- readRDS(file.path(objdir, "spot_scores.rds"))[sample_id == s]
  ho <- readRDS(file.path(objdir, "c1b_hotspots.rds"))[sample_id == s]
  x <- merge(sp[, .(barcode, array_row, array_col, C1b)], ho[, .(barcode, hotspot)], by = "barcode")
  x[, C1b_z := as.numeric(scale(C1b))]
  base <- function(v) ggplot(x, aes(array_col, -array_row, colour = .data[[v]])) +
    geom_point(size = 0.28, shape = 16) + coord_fixed(sqrt(3)) + theme_void(base_size = FS) +
    theme(plot.title = element_text(size = FS, hjust = 0.5, margin = margin(b = 2)), plot.margin = margin(4, 2, 2, 8),
          legend.key.size = grid::unit(2.5, "mm"),
          legend.text = element_text(size = FS - 1), legend.title = element_text(size = FS - 1))
  p1 <- base("C1b_z") + scale_colour_viridis(option = "magma", limits = c(-2.5, 2.5), oob = scales::squish,
                                             name = "C1b (z)") + ggtitle(title)
  p2 <- base("hotspot") + scale_colour_manual(values = c(hot = "#B2182B", cold = "#2166AC", ns = "grey88"),
                                              breaks = c("hot", "cold"), name = "Gi* hotspot\n(FDR<0.05)") +
    ggtitle(" ") + guides(colour = guide_legend(override.aes = list(size = 1.8)))
  p1 | p2
}
meta_d <- fread(file.path(TD, "sample_metadata.csv")); meta_r <- fread(file.path(TR, "sample_metadata.csv"))
pa <- (sec_map(OD, "mgh259", "MGH259\noligodendroglioma G3, discovery") |
       sec_map(OR, "idha17", "IDHA17\nastrocytoma G2, replication")) + plot_layout(guides = "collect")

# ---- b: coherence vs matched random gene sets --------------------------------------
mb <- rbind(rd(TD, "moran_I_C1b_vs_random_sets.csv")[cohort == "IDHm"][, cohort2 := "Discovery"],
            rd(TR, "moran_I_C1b_vs_random_sets.csv")[, cohort2 := "Replication"], fill = TRUE)
mb[, feature := factor(feature, c("C1b", "C1b_resid_macadj"), c("C1b", "C1b | TAM"))]
st <- rbind(rd(TD, "moran_vs_random_stouffer.csv")[cohort == "IDHm"][, cohort2 := "Discovery"],
            rd(TR, "moran_vs_random_stouffer.csv")[cohort == "IDHm"][, cohort2 := "Replication"])
st[, feature := factor(feature, c("C1b", "C1b_resid_macadj"), c("C1b", "C1b | TAM"))]
st[, lab := sprintf("Z=%.1f", stouffer_z)]
st[, ylab := fifelse(cohort2 == "Discovery", 8.0, 7.1)]
pb <- ggplot(mb, aes(feature, z_vs_random, colour = cohort2)) +
  geom_hline(yintercept = c(0, 1.645), linetype = c(1, 3), linewidth = .3) +
  geom_point(position = position_jitterdodge(jitter.width = .15, dodge.width = .7), size = .9, alpha = .85) +
  geom_text(data = st, aes(feature, ylab, label = lab), size = 2.1, show.legend = FALSE) +
  scale_colour_manual(values = COH, name = NULL) + coord_cartesian(ylim = c(-6.5, 8.5)) +
  labs(x = NULL, y = "Moran's I vs 1,000 matched\nrandom gene sets (z)") + th + theme(legend.position = "top")

# ---- c: specific co-localisation (random-gene-set baseline) ------------------------
sc <- rbind(rd(TD, "coloc_specificity_vs_random_sets.csv")[cohort %in% c("IDHm", "GBM")],
            rd(TR, "coloc_specificity_vs_random_sets.csv")[cohort == "IDHm"][, cohort := "IDHm_rep"])
sc[, col := factor(cohort, c("IDHm", "IDHm_rep", "GBM"), c("Disc.\nIDHm", "Rep.\nIDHm", "Disc.\nGBM"))]
sc[, signed := (n_above95 - n_below5) / k]
sc[, size := pmax(n_above95, n_below5) / k]
ordr <- sc[cohort == "IDHm_rep"][order(signed), MP]
sc[, MP := factor(gsub("\\.", "-", MP), gsub("\\.", "-", ordr))]
pc <- ggplot(sc, aes(col, MP, fill = signed, size = size)) + geom_point(shape = 21, colour = "grey30", stroke = .2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1),
                       name = "Specific\nco-loc. -\navoidance") +
  scale_size_area(max_size = 3.6, limits = c(0, 1), name = "Fraction of\nsections") +
  labs(x = NULL, y = NULL) + th + theme(axis.text.x = element_text(size = FS - 1))

# ---- d: niche composition (inside hotspots, depth-adjusted, MSR) ------------------
tgt <- c("Endothelial", "Pericyte", "TAM_C1Q", "TAM_microglia", "TAM_SPP1_lipid", "TAM_monocyte",
         "Inflammation", "IFN_response", "Hypoxia")
nd <- rbind(rd(TD, "niche_proximity_meta.csv")[cohort == "IDHm"][, cohort2 := "Discovery"],
            rd(TR, "niche_proximity_meta.csv")[cohort == "IDHm"][, cohort2 := "Replication"])
nd <- nd[source == "C1b" & zone == "inside" & target %in% tgt, .(cohort2, target, stouffer_z, q, k)]
na <- rbind(rd(TD, "meta_apc_vs_nonMHC_proximity.csv")[cohort == "IDHm"][, cohort2 := "Discovery"],
            rd(TR, "meta_apc_vs_nonMHC_proximity.csv")[cohort == "IDHm"][, cohort2 := "Replication"])
na <- na[source == "C1b" & zone == "inside" & target %in% c("T_cell", "IFNG_response"), .(cohort2, target, stouffer_z, q, k)]
nd <- rbind(nd, na)
labs_t <- c(Endothelial = "Endothelium", Pericyte = "Pericytes", TAM_C1Q = "C1Q+ TAM", TAM_microglia = "Microglia-like TAM",
            TAM_SPP1_lipid = "SPP1+/lipid TAM", TAM_monocyte = "Monocyte-like TAM", Inflammation = "Inflammation",
            IFN_response = "Type I IFN", IFNG_response = "IFN-gamma", Hypoxia = "Hypoxia", T_cell = "T cells")
ord_d <- nd[cohort2 == "Replication"][order(stouffer_z), target]
nd[, target := factor(target, ord_d, labs_t[ord_d])]
nd[, sig := q < 0.05]
pd <- ggplot(nd, aes(stouffer_z, target, colour = cohort2, shape = sig)) +
  geom_vline(xintercept = 0, linewidth = .3) + geom_vline(xintercept = c(-1.96, 1.96), linetype = 3, linewidth = .3) +
  geom_point(size = 1.6, position = position_dodge(width = .6)) +
  scale_colour_manual(values = COH, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "q<0.05") +
  labs(x = "Enrichment inside C1b hotspots\n(Stouffer z; MSR null, depth-adjusted)", y = NULL) + th

# ---- e: protein proxy -------------------------------------------------------------------
pp <- rbind(rd(TD, "protein_panel_proxy_per_sample.csv")[cohort == "IDHm"][, cohort2 := "Discovery"],
            rd(TR, "protein_panel_proxy_per_sample.csv")[, cohort2 := "Replication"], fill = TRUE)
pp[, panel := factor(panel, c("P5", "P3m", "P3"), c("MHCII/CD163/\nCD206/CD44/VIM", "MHCII/CD163/\nCD206", "CD68/CD14/\nHLA-DPB1"))]
beat <- pp[, .(lab = sprintf("%d/%d", sum(pct_C1b_subsets_below > .95), .N)), by = .(cohort2, panel)]
pe <- ggplot(pp, aes(panel, rho, colour = cohort2)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = .6, linewidth = .3, position = position_dodge(width = .75)) +
  geom_point(position = position_jitterdodge(jitter.width = .12, dodge.width = .75), size = .7, alpha = .8) +
  geom_text(data = beat, aes(panel, 0.82, label = lab, group = cohort2), position = position_dodge(width = .75),
            size = 2, show.legend = FALSE) +
  annotate("text", x = 0.45, y = 0.9, hjust = 0, size = 2, label = "sections beating 95% of random C1b subsets") +
  scale_colour_manual(values = COH, name = NULL) + coord_cartesian(ylim = c(0, 0.93)) +
  labs(x = NULL, y = "Spearman rho with\n90-gene C1b") + th + theme(legend.position = "none")

# ---- f: CODEX field of view --------------------------------------------------------------
CS <- "IDHA14"
cells <- as.data.table(read_parquet(file.path(ROOT, "data", "raw", "codex", "cells_df.parquet"),
                                    col_select = c("cell_name", "sample", "centroid_x", "centroid_y", "cell_type")))[sample == CS]
tams <- fread(file.path(OD, "codex_tams.csv.gz"))[sample == CS, .(cell_name, C1b_tertile)]
cells <- merge(cells, tams, by = "cell_name", all.x = TRUE)
win <- 600
Tc <- cells[cell_type %in% c("TcellCD4", "TcellCD8")]
# window centred on the 200-um bin with the most T cells in this section
gx <- floor(Tc$centroid_x / 200); gy <- floor(Tc$centroid_y / 200)
best <- as.data.table(table(paste(gx, gy)))[order(-N)][1, V1]
bx <- (as.numeric(strsplit(best, " ")[[1]]) + 0.5) * 200
x0 <- bx[1] - win / 2; y0 <- bx[2] - win / 2
fv <- cells[centroid_x >= x0 & centroid_x < x0 + win & centroid_y >= y0 & centroid_y < y0 + win]
fv[, cls := fcase(cell_type == "TcellCD4", "CD4+ T", cell_type == "TcellCD8", "CD8+ T", cell_type == "Vasc", "Vascular",
                  C1b_tertile == "high", "C1b-high TAM", C1b_tertile == "low", "C1b-low TAM",
                  cell_type == "Mac", "Other TAM", default = "Other")]
fv[, cls := factor(cls, c("Other", "Other TAM", "C1b-low TAM", "Vascular", "C1b-high TAM", "CD8+ T", "CD4+ T"))]
setorder(fv, cls)
pf <- ggplot(fv, aes(centroid_x, centroid_y, colour = cls, size = cls)) + geom_point(shape = 16) +
  scale_colour_manual(values = c(Other = "grey90", `Other TAM` = "grey70", `C1b-low TAM` = "#92C5DE",
                                 Vascular = "#1B7837", `C1b-high TAM` = "#762A83", `CD8+ T` = "#F4A582", `CD4+ T` = "#B2182B"),
                      name = NULL) +
  scale_size_manual(values = c(Other = .3, `Other TAM` = .6, `C1b-low TAM` = .9, Vascular = .7,
                               `C1b-high TAM` = 1.1, `CD8+ T` = 1.8, `CD4+ T` = 1.8), guide = "none") +
  annotate("segment", x = x0 + 30, xend = x0 + 130, y = y0 + 25, yend = y0 + 25, linewidth = .6) +
  annotate("text", x = x0 + 80, y = y0 + 55, label = "100 \u00b5m", size = 2) +
  scale_y_reverse() + coord_fixed() + theme_void(base_size = FS) +
  ggtitle(sprintf("CODEX %s (astrocytoma G3), 600 x 600 \u00b5m", CS)) +
  theme(plot.title = element_text(size = FS, hjust = .5), legend.text = element_text(size = FS - 1),
        legend.key.size = grid::unit(2.5, "mm")) +
  guides(colour = guide_legend(override.aes = list(size = 1.8), ncol = 1))

# ---- g: CODEX meta-analysis by grade -----------------------------------------------------------
cg <- rd(TD, "codex_meta_by_grade.csv")[target %in% c("T", "CD4", "CD8") & grp %in% c("IDHm G2", "IDHm G3", "IDHm G2-3", "IDHm G4", "GBM")]
cg[, grp := factor(grp, rev(c("IDHm G2", "IDHm G3", "IDHm G2-3", "IDHm G4", "GBM")),
                   rev(c("IDHm G2", "IDHm G3", "IDHm G2-3", "IDHm G4", "GBM")))]
cg[, target := factor(target, c("T", "CD4", "CD8"), c("Any T cell", "CD4+ T", "CD8+ T"))]
cg[, contrast := factor(contrast, c("C1b_high_vs_low", "MHCII_high_vs_low"), c("C1b 5-protein high vs low", "MHC-II high vs low"))]
pg <- ggplot(cg, aes(exp(estimate), grp, colour = contrast)) + geom_vline(xintercept = 1, linewidth = .3) +
  geom_errorbarh(aes(xmin = exp(ci_lb), xmax = exp(ci_ub)), height = 0, linewidth = .4, position = position_dodge(width = .6)) +
  geom_point(size = 1.2, position = position_dodge(width = .6)) + facet_wrap(~target, nrow = 1) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8)) +
  scale_colour_manual(values = c("#762A83", "#E08214"), name = NULL) +
  labs(x = "Adjusted OR: T cell within 27.5 \u00b5m of TAM\n(random-effects meta-analysis)", y = NULL) +
  th + theme(legend.position = "top")

# ---- h: antigen-presenting vs mesenchymal-like arm ----------------------------------------------
ah <- rbind(rd(TR, "meta_apc_vs_nonMHC_proximity.csv")[cohort == "IDHm"][, cohort2 := "Replication IDHm"],
            rd(TD, "meta_apc_vs_nonMHC_proximity.csv")[cohort == "GBM"][, cohort2 := "Discovery GBM"])
ah <- ah[zone == "inside" & target %in% c("T_cell", "Endothelial", "IFN_response", "Hypoxia", "Inflammation")]
ah[, target := factor(target, rev(c("T_cell", "Endothelial", "IFN_response", "Inflammation", "Hypoxia")),
                      rev(c("T cells", "Endothelium", "Type I IFN", "Inflammation", "Hypoxia")))]
ah[, source := factor(source, c("C1b", "APC", "nonMHC"), c("C1b (full)", "MHC-II arm | TAM", "Mesenchymal arm | TAM"))]
ah[, cohort2 := factor(cohort2, c("Replication IDHm", "Discovery GBM"))]
ph <- ggplot(ah, aes(stouffer_z, target, colour = source)) +
  geom_vline(xintercept = 0, linewidth = .3) + geom_vline(xintercept = c(-1.96, 1.96), linetype = 3, linewidth = .3) +
  geom_point(size = 1.4, position = position_dodge(width = .65)) + facet_wrap(~cohort2, nrow = 1) +
  scale_colour_manual(values = c("#762A83", "#E08214", "#1B7837"), name = NULL) +
  labs(x = "Enrichment inside hotspots (Stouffer z)", y = NULL) + th +
  theme(legend.position = "top") + guides(colour = guide_legend(nrow = 2))

# ---- assemble ---------------------------------------------------------------------------
W <- function(p) wrap_elements(full = p)
fig <- (W(pa)) / (W(pb) + W(pc) + W(pd) + plot_layout(widths = c(0.8, 1, 1.3))) /
       (W(pe) + W(pf) + plot_layout(widths = c(1, 1))) /
       (W(pg) + W(ph) + plot_layout(widths = c(1.15, 1))) +
  plot_layout(heights = c(0.8, 1.2, 1.15, 1.05)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 9, face = "bold"))
savep(fig, "Figure4.pdf", 180, 245)
ggsave(file.path(OUT, "Figure4_preview.png"), fig, width = mm(180), height = mm(245), units = "in", dpi = 220)

# ---- individual panels ----------------------------------------------------------------------
savep(pa, "Figure4a.pdf", 180, 50); savep(pb, "Figure4b.pdf", 55, 65); savep(pc, "Figure4c.pdf", 70, 65)
savep(pd, "Figure4d.pdf", 90, 65); savep(pe, "Figure4e.pdf", 90, 65); savep(pf, "Figure4f.pdf", 90, 70)
savep(pg, "Figure4g.pdf", 100, 60); savep(ph, "Figure4h.pdf", 90, 60)
fwrite(data.table(panel = "f", sample = CS, x0 = x0, y0 = y0, window_um = win,
                  n_cells = nrow(fv), n_T = sum(fv$cls %in% c("CD4+ T", "CD8+ T"))),
       file.path(OUT, "Figure4f_field_of_view.csv"))

# ---- supplementary: renamed copies of existing PDFs ------------------------------------------
FIG <- file.path(ROOT, "results", "figures"); FIGR <- file.path(ROOT, "results", "replication", "figures")
supp <- list(
  S4a_QC_discovery = file.path(FIG, "FigS1_QC.pdf"), S4a_QC_replication = file.path(FIGR, "FigS1_QC.pdf"),
  S4b_C1b_maps_IDHm_discovery = file.path(FIG, "Fig1_C1b_spatial_IDHm.pdf"),
  S4b_C1b_maps_IDHm_replication = file.path(FIGR, "Fig1_C1b_spatial_IDHm.pdf"),
  S4b_C1b_maps_GBM_discovery = file.path(FIG, "FigS2_C1b_spatial_GBM.pdf"),
  S4c_MoranI_vs_random_discovery = file.path(FIG, "Fig2_MoranI_vs_random.pdf"),
  S4c_MoranI_vs_random_replication = file.path(FIGR, "Fig2_MoranI_vs_random.pdf"),
  S4d_coloc_heatmap_discovery = file.path(FIG, "Fig3a_coloc_rho_heatmap.pdf"),
  S4d_coloc_heatmap_replication = file.path(FIGR, "Fig3a_coloc_rho_heatmap.pdf"),
  S4d_coloc_partial_discovery = file.path(FIG, "Fig3b_coloc_partial_rho_heatmap.pdf"),
  S4e_nhood_enrichment_discovery = file.path(FIG, "Fig4_nhood_enrichment_heatmap.pdf"),
  S4e_nhood_enrichment_replication = file.path(FIGR, "Fig4_nhood_enrichment_heatmap.pdf"),
  S4f_CNA_heatmap_discovery = file.path(FIG, "Fig5a_CNA_heatmap.pdf"),
  S4f_CNA_heatmap_replication = file.path(FIGR, "Fig5a_CNA_heatmap.pdf"),
  S4f_CNA_regions_discovery = file.path(FIG, "Fig5b_CNA_regions_C1b.pdf"),
  S4f_CNA_regions_replication = file.path(FIGR, "Fig5b_CNA_regions_C1b.pdf"),
  S4g_3gene_concordance_discovery = file.path(FIG, "Fig6a_C1b_vs_3gene_concordance.pdf"),
  S4g_3gene_concordance_replication = file.path(FIGR, "Fig6a_C1b_vs_3gene_concordance.pdf"),
  S4g_3gene_maps_discovery = file.path(FIG, "Fig6b_C1b_vs_3gene_maps_IDHm.pdf"),
  S4h_proximity_per_section_discovery = file.path(FIG, "Fig7a_niche_proximity_C1b.pdf"),
  S4h_proximity_per_section_replication = file.path(FIGR, "Fig7a_niche_proximity_C1b.pdf"),
  S4h_proximity_decay_discovery = file.path(FIG, "Fig7b_niche_proximity_meta_decay.pdf"),
  S4h_proximity_decay_replication = file.path(FIGR, "Fig7b_niche_proximity_meta_decay.pdf"),
  S4i_APC_index_maps_discovery = file.path(FIG, "Fig8a_APC_index_maps_IDHm.pdf"),
  S4i_APC_index_maps_replication = file.path(FIGR, "Fig8a_APC_index_maps_IDHm.pdf"),
  S4j_pseudobulk_Tcells_discovery = file.path(FIG, "Fig8c_niche_pseudobulk_Tcells.pdf"),
  S4j_pseudobulk_Tcells_replication = file.path(FIGR, "Fig8c_niche_pseudobulk_Tcells.pdf"),
  S4k_ligand_receptor_discovery = file.path(FIG, "Fig8d_ligand_receptor_colocalisation.pdf"),
  S4k_ligand_receptor_replication = file.path(FIGR, "Fig8d_ligand_receptor_colocalisation.pdf"),
  S4l_CODEX_vessel_hypoxia = file.path(FIG, "Fig10a_CODEX_TAM_Tcell_vessel.pdf"),
  S4l_CODEX_radius_subtype = file.path(FIG, "Fig10b_CODEX_radius_subtype.pdf"))
ok <- vapply(names(supp), function(n) file.copy(supp[[n]], file.path(OUT, "supplementary", paste0("Figure", n, ".pdf")),
                                                  overwrite = TRUE), TRUE)
if (!all(ok)) logf("missing supplementary source(s): %s", paste(names(supp)[!ok], collapse = ", "))
logf("Figure 4 written to %s (%d supplementary PDFs)", OUT, sum(ok))
