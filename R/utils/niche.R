# Shared helpers for the niche-proximity analyses (steps 08 and 09) ----------
suppressPackageStartupMessages({ library(spdep); library(adespatial); library(Matrix) })

PANEL <- c("CD68", "CD14", "HLA-DPB1")

# ---- target signatures ------------------------------------------------------
# Curated marker sets (canonical markers; TAM states follow Pombo Antunes et al.
# Nat Neurosci 2021, Abdelfattah et al. Nat Commun 2022, Miller et al. 2025
# glioma myeloid programs). Genes shared with C1b or the 3-gene panel are
# removed (see clean_targets) to avoid circularity.
TARGETS_RAW <- list(
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

# APC-TAM specific additions (step 09): IFN-gamma response kept separate from
# the type-I ISG set; T-cell subsets for pseudobulk niche analysis. CD4 itself
# is left out of the CD4 T set because TAMs also express CD4.
TARGETS_APC <- list(
  IFNG_response = c("CIITA", "CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1", "GBP2", "GBP4", "GBP5",
                    "STAT1", "IRF1", "IFNG")
)
TCELL_SETS <- list(
  T_all  = c("CD3D", "CD3E", "CD3G", "CD2", "CD247", "TRAC"),
  CD4_T  = c("CD40LG", "IL7R", "TRAT1", "ICOS"),
  CD8_T  = c("CD8A", "CD8B", "GZMK", "GZMA", "NKG7"),
  Treg   = c("FOXP3", "IL2RA", "CTLA4"),
  IFNG   = c("IFNG")
)

clean_targets <- function(tl) {
  excl <- c(C1B_GENES, map_genes(C1B_GENES, C1B_GENES), PANEL, "MARCHF1")
  lapply(tl, setdiff, excl)
}

# ---- small helpers ------------------------------------------------------------
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

# Full (unfiltered) GEO count matrix restricted to the QC-passing spots; used
# where low-abundance genes removed by the >=10-spot gene filter are needed
read_raw_counts <- function(sid, barcodes) {
  m <- read_sample_meta()
  pref <- m$prefix[m$sample_id == sid]
  x <- Seurat::Read10X_h5(file.path(P$geo_dir, paste0(pref, "_filtered_feature_bc_matrix.h5")))
  if (is.list(x)) x <- x[["Gene Expression"]]
  rownames(x) <- make.unique(rownames(x))
  x[, barcodes]
}
