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
#     - targets residualised on log10(UMI depth); depth reported as a target
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
source("R/utils/niche.R")                    # PANEL, TARGETS_RAW, helpers
TARGETS <- clean_targets(TARGETS_RAW)
write_tab(data.table(signature = names(TARGETS), n_genes = lengths(TARGETS),
                     genes = vapply(TARGETS, paste, "", collapse = ";")), "target_signatures_used.csv")

meta  <- available_samples()
only  <- Sys.getenv("C1B_SAMPLES")          # optional: comma-separated subset for testing
if (nzchar(only)) meta <- meta[sample_id %in% strsplit(only, ",")[[1]]]
spots <- readRDS(file.path(P$objects, "spot_scores.rds"))
hot   <- readRDS(file.path(P$objects, "c1b_hotspots.rds"))
cna_f <- file.path(P$objects, "cna_spot.rds")
cna   <- if (file.exists(cna_f)) readRDS(cna_f) else NULL

tgenes_rows <- list(); conc_rows <- list(); prox_rows <- list(); decay_rows <- list(); trio_rows <- list()
for (i in seq_len(nrow(meta))) {
  s <- meta$sample_id[i]
  so <- readRDS(obj_path(s))
  counts  <- GetAssayData(so, assay = "Spatial", layer = "counts")
  lognorm <- GetAssayData(so, assay = "Spatial", layer = "data")
  d <- merge(spots[sample_id == s, .(barcode, array_row, array_col, C1b, MP_Mac, nCount)],
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
  # Depth control: every target is residualised on log10(UMI) so that
  # co-occurrence of high-depth regions cannot create spurious proximity;
  # depth itself is carried as a reference "target".
  ldepth <- log10(d$nCount)
  for (t in targets) set(tg, j = t, value = residuals(lm(tg[[t]] ~ ldepth)))
  tg[, UMI_depth := ldepth]
  targets <- c(targets, "UMI_depth")
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

logf("08 computations done; figures are drawn by R/08b_plot_niche.R")
