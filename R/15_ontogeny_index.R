# 15 - Monocyte-derived vs microglial TAM identity inside C1b hotspots
#
#  Visium spots mix cells, so a TAM-rich hotspot scores high on every TAM
#  signature. To ask whether the TAMs in C1b hotspots lean towards a
#  monocyte-derived (MDM) or a microglial (MG) identity, each identity score is
#  residualised on a pan-TAM score (genes shared by MG and MDM) and on
#  log10(UMI depth): the residual is "identity per unit TAM".
#   - MG identity: homeostatic microglia genes (Butovsky 2014; Bennett 2016;
#     Gosselin 2017)
#   - MDM identity: genes reported higher in bone-marrow-derived than in
#     microglia-derived brain/glioma macrophages (ITGA4/CD49d and ITGAL/CD11a:
#     Bowman 2016, Friebel 2020; CLEC12A, EMILIN2, GDA, MS4A7: Bowman 2016 and
#     follow-up mouse brain-macrophage studies) - treat as provisional in human
#   - Mo-vs-MG index = z(MDM residual) - z(MG residual)
#  Genes shared with C1b or the 3-gene panel are removed. Enrichment inside C1b
#  hotspots is tested against 499 MSR surrogates of each score (as step 08) and
#  combined per cohort with Stouffer's method.
source("R/utils/common.R")
source("R/utils/scoring.R")
suppressPackageStartupMessages({ library(Seurat); library(spdep); library(adespatial) })
set.seed(CFG$seed)
N_MSR <- 499
source("R/utils/niche.R")

ONTO_RAW <- list(
  MG_identity  = c("P2RY12", "P2RY13", "TMEM119", "CX3CR1", "SALL1", "GPR34", "SIGLEC8", "SLC2A5",
                   "OLFML3", "SELPLG"),
  MDM_identity = c("ITGA4", "ITGAL", "CLEC12A", "EMILIN2", "GDA", "MS4A7"),
  panTAM       = c("AIF1", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "ITGAM", "LAPTM5", "SPI1",
                   "CD53", "CSF1R")
)
ONTO <- clean_targets(ONTO_RAW)
write_tab(data.table(signature = names(ONTO), n_genes = lengths(ONTO),
                     genes = vapply(ONTO, paste, "", collapse = ";")), "ontogeny_signatures_used.csv")

meta  <- available_samples()
spots <- readRDS(file.path(P$objects, "spot_scores.rds"))
hot   <- readRDS(file.path(P$objects, "c1b_hotspots.rds"))
resid_on <- function(y, X) as.numeric(scale(resid(lm(y ~ ., data = as.data.frame(X)))))

rows <- list(); decay <- list(); det <- list()
for (i in seq_len(nrow(meta))) {
  s <- meta$sample_id[i]
  d <- merge(spots[sample_id == s, .(barcode, array_row, array_col, nCount)],
             hot[sample_id == s, .(barcode, hotspot)], by = "barcode", sort = FALSE)
  H <- which(d$hotspot == "hot")
  if (length(H) < 5) next
  so <- readRDS(obj_path(s))
  counts <- GetAssayData(so, assay = "Spatial", layer = "counts")
  d <- d[match(colnames(so), barcode)]
  H <- which(d$hotspot == "hot")
  det[[s]] <- data.table(sample_id = s, signature = names(ONTO),
                         n_detected = vapply(ONTO, function(g) sum(g %in% rownames(counts)), 0L),
                         n_total = lengths(ONTO))
  sc <- ucell_scores(counts, ONTO, maxRank = CFG$scoring$ucell_maxrank, ncores = CFG$n_cores, min_genes = 3)
  if (!all(c("MG_identity", "panTAM") %in% colnames(sc))) { logf("%s: MG/panTAM not scorable, skipped", s); next }
  depth <- log10(d$nCount)
  tg <- list(MG_raw = resid_on(sc[, "MG_identity"], data.frame(depth)),
             MG_perTAM = resid_on(sc[, "MG_identity"], data.frame(pan = sc[, "panTAM"], depth)),
             panTAM = resid_on(sc[, "panTAM"], data.frame(depth)))
  if ("MDM_identity" %in% colnames(sc)) {
    tg$MDM_raw    <- resid_on(sc[, "MDM_identity"], data.frame(depth))
    tg$MDM_perTAM <- resid_on(sc[, "MDM_identity"], data.frame(pan = sc[, "panTAM"], depth))
    tg$Mo_vs_MG   <- as.numeric(scale(tg$MDM_perTAM - tg$MG_perTAM))
  }
  lw <- visium_listw(d)
  me <- adespatial::scores.listw(lw, MEM.autocor = "all")
  hd <- hop_distance(lw$neighbours, H)
  inside <- hd == 0
  for (t in names(tg)) {
    y <- tg[[t]]
    sur <- scale(adespatial::msr(y, me, nrepet = N_MSR, method = "pair"))
    obs <- mean(y[inside]); nul <- colMeans(sur[inside, , drop = FALSE])
    rows[[paste(s, t)]] <- data.table(sample_id = s, target = t, n_hot = sum(inside), mean_z = obs,
                                      z_msr = (obs - mean(nul)) / sd(nul),
                                      p_enrich = (1 + sum(nul >= obs)) / (N_MSR + 1),
                                      p_deplete = (1 + sum(nul <= obs)) / (N_MSR + 1))
    decay[[paste(s, t)]] <- data.table(sample_id = s, target = t, dist = pmin(hd, 6L), z = y)[
      , .(mean_z = mean(z), n = .N), by = .(sample_id, target, dist)]
  }
  logf("%s: %d hotspot spots; Mo-vs-MG inside z_msr = %s", s, sum(inside),
       if ("Mo_vs_MG" %in% names(tg)) sprintf("%.2f", rows[[paste(s, "Mo_vs_MG")]]$z_msr) else "NA")
  rm(so, counts, me); gc(verbose = FALSE)
}

per <- merge(rbindlist(rows), meta[, .(sample_id, cohort, who_grade, is_LGG)], by = "sample_id")
write_tab(per, "ontogeny_hotspot_per_sample.csv")
write_tab(rbindlist(det), "ontogeny_gene_detection.csv")
write_tab(merge(rbindlist(decay), meta[, .(sample_id, cohort)], by = "sample_id"), "ontogeny_decay.csv")
pm <- per[is.finite(z_msr), .(k = .N, mean_z_msr = mean(z_msr), stouffer_z = sum(z_msr) / sqrt(.N),
                              n_enriched = sum(p_enrich < 0.05), n_depleted = sum(p_deplete < 0.05)),
          by = .(cohort, target)]
pm[, p := 2 * pnorm(-abs(stouffer_z))][, q := p.adjust(p, "BH"), by = cohort]
write_tab(pm, "ontogeny_hotspot_meta.csv")
print(pm[order(cohort, target)])
