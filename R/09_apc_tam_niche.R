# 09 - Consensus1b as an antigen-presenting TAM (APC-TAM) program ------------
#
#  1. Separate antigen-presentation intensity from TAM abundance:
#       APC index    = C1b MHC-II module | Mac MP + log10 UMI   (residual)
#       nonMHC index = rest of C1b (non-MHC-II, non-ribosomal) | same covariates
#     Moran's I vs 1,000 expression-matched random MHC-sized gene sets
#     (residualised identically); Gi* hotspots; overlap with C1b hotspots.
#  2. Proximity of APC-high vs nonMHC-high vs C1b hotspots to cell types /
#     programs incl. a separate IFN-gamma response set (MSR null, as step 08).
#  3. T cells around C1b / APC niches by niche pseudobulk (spots of the niche
#     and its 1-2 spot ring pooled), against 1,000 random regions of identical
#     shape (lattice translations / reflections kept inside the tissue).
#  4. Ligand-receptor spatial co-localisation (spot + first-order neighbours):
#     statistic mean(zL * W_self zR), null from 499 MSR surrogates of zR.
#
#  5. TAM-adjusted LR: ligand and receptor residualised on log10 UMI + Mac MP,
#     separating APC-T cell co-localisation from myeloid co-expression (TAMs
#     themselves express CD4 / HAVCR2 / CD44).
#
# Per-sample results are cached (data/objects/<sample>_apc09.rds) so the
# script can be resumed or split across processes (C1B_SAMPLES=a,b,...);
# every run re-aggregates all cached samples. FORCE_APC=1 recomputes.
source("R/utils/common.R")
source("R/utils/scoring.R")
source("R/utils/niche.R")
suppressPackageStartupMessages({ library(Seurat); library(spdep); library(adespatial) })

N_MSR    <- 499
N_REGION <- 1000
N_NULL   <- 1000
TARGETS  <- clean_targets(c(TARGETS_RAW, TARGETS_APC))
MHC_MOD  <- c1b_variants()$C1b_MHCII
MHCII_LIGAND <- c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1")
LR <- data.table(
  pair     = c("MHCII-CD4", "CD80-CD28", "CD86-CD28", "CD80-CTLA4", "CD86-CTLA4", "CD274-PDCD1",
               "PDCD1LG2-PDCD1", "LGALS9-HAVCR2", "CD48-CD2", "CD48-CD244", "SPP1-CD44"),
  ligand   = c("MHCII", "CD80", "CD86", "CD80", "CD86", "CD274", "PDCD1LG2", "LGALS9", "CD48", "CD48", "SPP1"),
  receptor = c("CD4", "CD28", "CD28", "CTLA4", "CTLA4", "PDCD1", "PDCD1", "HAVCR2", "CD2", "CD244", "CD44"))

meta  <- available_samples()
only  <- Sys.getenv("C1B_SAMPLES")
run_ids <- if (nzchar(only)) intersect(meta$sample_id, strsplit(only, ",")[[1]]) else meta$sample_id
spots <- readRDS(file.path(P$objects, "spot_scores.rds"))
hot   <- readRDS(file.path(P$objects, "c1b_hotspots.rds"))
cna   <- readRDS(file.path(P$objects, "cna_spot.rds"))

# random regions with the same shape as `idx` (translations + reflections on
# the hex lattice; parity of row + col preserved), >= 80% inside the tissue
random_regions <- function(d, idx, n, min_in = 0.8, max_try = 40 * n) {
  key <- paste(d$array_row, d$array_col)
  r0 <- d$array_row[idx]; c0 <- d$array_col[idx]
  out <- vector("list", n); k <- 0L; tries <- 0L
  while (k < n && tries < max_try) {
    tries <- tries + 1L
    fr <- sample(c(-1L, 1L), 1); fc <- sample(c(-1L, 1L), 1)
    r1 <- fr * r0; c1 <- fc * c0
    dr <- sample(seq(min(d$array_row) - min(r1), max(d$array_row) - max(r1) + 1L), 1)
    dc <- sample(seq(min(d$array_col) - min(c1), max(d$array_col) - max(c1) + 1L), 1)
    if ((dr + dc) %% 2 != 0) dc <- dc + 1L
    m <- match(paste(r1 + dr, c1 + dc), key)
    if (mean(!is.na(m)) < min_in) next
    k <- k + 1L; out[[k]] <- unique(m[!is.na(m)])
  }
  out[seq_len(k)]
}

# mid-p permutation p-value: ties (common with sparse counts, e.g. 0 vs 0)
# count one half, so that absent signal is not scored as depletion
midp <- function(obs, null, side = c("greater", "less")) {
  side <- match.arg(side)
  gt <- if (side == "greater") sum(null > obs) else sum(null < obs)
  (1 + gt + 0.5 * sum(null == obs)) / (length(null) + 1)
}

msr_z <- function(obs, null) c(z = (obs - mean(null)) / sd(null),
                               p_pos = (1 + sum(null >= obs)) / (length(null) + 1),
                               p_neg = (1 + sum(null <= obs)) / (length(null) + 1))

run_sample <- function(i) {
  s <- meta$sample_id[i]; seed <- CFG$seed + 900 + i
  set.seed(seed)
  so <- readRDS(obj_path(s))
  counts  <- GetAssayData(so, assay = "Spatial", layer = "counts")
  lognorm <- GetAssayData(so, assay = "Spatial", layer = "data")
  raw <- read_raw_counts(s, colnames(so))
  raw_ln <- raw %*% Diagonal(x = 1e4 / Matrix::colSums(raw)); raw_ln@x <- log1p(raw_ln@x)
  dimnames(raw_ln) <- dimnames(raw)

  d <- spots[sample_id == s][match(colnames(so), barcode)]
  d[, hotspot := hot[sample_id == s][match(d$barcode, barcode), hotspot]]
  lw <- visium_listw(d); W <- listw_to_sparse(lw)
  lw_b <- nb2listw(include.self(lw$neighbours), style = "B", zero.policy = TRUE)
  Ws <- listw_to_sparse(nb2listw(include.self(lw$neighbours), style = "W", zero.policy = TRUE))
  me <- adespatial::scores.listw(lw, MEM.autocor = "all")
  ldepth <- log10(d$nCount)

  # ---- 1. APC index -----------------------------------------------------------
  apc    <- residuals(lm(d$C1b_MHCII ~ d$MP_Mac + ldepth))
  nonmhc <- residuals(lm(d$C1b_nonMHCII ~ d$MP_Mac + ldepth))
  mhc <- unname(na.omit(map_genes(MHC_MOD, rownames(counts))))
  rs  <- matched_random_sets(lognorm, mhc, n_sets = N_NULL, seed = seed)
  rsc <- ucell_scores(counts, rs, maxRank = CFG$scoring$ucell_maxrank, ncores = CFG$n_cores, min_genes = 3)
  rres <- apply(rsc, 2, function(y) residuals(lm(y ~ d$MP_Mac + ldepth)))
  I_apc <- moran_matrix(cbind(apc), W); I_null <- moran_matrix(rres, W)
  I_raw <- moran_matrix(cbind(d$C1b_MHCII), W); I_raw_null <- moran_matrix(rsc, W)
  I_non <- moran_matrix(cbind(nonmhc), W)
  h_c1b <- d$hotspot == "hot"; h_apc <- gi_hot(apc, lw_b); h_non <- gi_hot(nonmhc, lw_b)
  moran <- data.table(sample_id = s, n_mhc_genes = length(mhc),
    I_MHCII_raw = I_raw, p_MHCII_raw_vs_random = (1 + sum(I_raw_null >= I_raw)) / (N_NULL + 1),
    I_APC_index = I_apc, I_APC_null_mean = mean(I_null), I_APC_null_sd = sd(I_null),
    z_APC_vs_random = (I_apc - mean(I_null)) / sd(I_null),
    p_APC_vs_random = (1 + sum(I_null >= I_apc)) / (N_NULL + 1),
    I_nonMHC_index = I_non, rho_APC_nonMHC = cor(apc, nonmhc, method = "spearman"),
    n_hot_C1b = sum(h_c1b), n_hot_APC = sum(h_apc), n_hot_nonMHC = sum(h_non),
    dice_APC_C1b = dice(h_apc, h_c1b), dice_nonMHC_C1b = dice(h_non, h_c1b), dice_APC_nonMHC = dice(h_apc, h_non))

  # ---- 2. proximity (targets from full raw counts, depth-residualised) -------
  tg <- as.data.table(ucell_scores(raw, TARGETS, maxRank = CFG$scoring$ucell_maxrank,
                                   ncores = CFG$n_cores, min_genes = 3))
  tg[, Malignant_CNA := cna[sample_id == s][match(d$barcode, barcode), CNAtot]]
  targets <- names(tg)[vapply(tg, function(x) !anyNA(x) && sd(x) > 0, TRUE)]
  for (t in targets) set(tg, j = t, value = residuals(lm(tg[[t]] ~ ldepth)))
  surs <- lapply(setNames(targets, targets), function(t)
    scale(adespatial::msr(tg[[t]], me, nrepet = N_MSR, method = "pair")))
  prox <- list()
  for (src in c("C1b", "APC", "nonMHC")) {
    H <- which(switch(src, C1b = h_c1b, APC = h_apc, nonMHC = h_non))
    if (length(H) < 5) next
    hd <- hop_distance(lw$neighbours, H)
    for (t in targets) {
      ys <- scale(tg[[t]])[, 1]
      for (zone in c("inside", "ring")) {
        idx <- if (zone == "inside") hd == 0 else hd %in% 1:2
        st <- msr_z(mean(ys[idx]), colMeans(surs[[t]][idx, , drop = FALSE]))
        prox[[paste(src, t, zone)]] <- data.table(sample_id = s, source = src, target = t, zone = zone,
          n_spots = sum(idx), mean_z = mean(ys[idx]), z_msr = st[["z"]], p_enrich = st[["p_pos"]],
          p_deplete = st[["p_neg"]])
      }
    }
  }

  # ---- 3. niche pseudobulk T cells ---------------------------------------------
  tot <- Matrix::colSums(raw)
  gs <- lapply(TCELL_SETS, function(g) Matrix::colSums(raw[intersect(g, rownames(raw)), , drop = FALSE]))
  pb <- list()
  for (src in c("C1b", "APC")) {
    H <- which(if (src == "C1b") h_c1b else h_apc)
    if (length(H) < 5) next
    niche <- which(hop_distance(lw$neighbours, H) <= 2)
    regs <- random_regions(d, niche, N_REGION)
    for (ts in names(gs)) {
      obs <- sum(gs[[ts]][niche]) / sum(tot[niche])
      nul <- vapply(regs, function(r) sum(gs[[ts]][r]) / sum(tot[r]), 0)
      pb[[paste(src, ts)]] <- data.table(sample_id = s, niche = src, tset = ts, n_niche_spots = length(niche),
        umi_set_in_niche = sum(gs[[ts]][niche]), umi_set_total = sum(gs[[ts]]),
        cpm_niche = 1e6 * obs, cpm_null_mean = 1e6 * mean(nul), n_regions = length(regs),
        log2_ratio = log2((obs + 1e-7) / (mean(nul) + 1e-7)),
        z = if (sd(nul) > 0) (obs - mean(nul)) / sd(nul) else NA_real_,
        p_enrich = midp(obs, nul, "greater"), p_deplete = midp(obs, nul, "less"))
    }
  }

  # ---- 4. ligand-receptor co-localisation ---------------------------------------
  gexp <- function(g) {
    if (g == "MHCII") return(Matrix::colMeans(raw_ln[intersect(MHCII_LIGAND, rownames(raw_ln)), , drop = FALSE]))
    if (!g %in% rownames(raw_ln)) return(NULL)
    as.numeric(raw_ln[g, ])
  }
  prep <- function(x) { if (is.null(x) || sd(x) == 0) return(NULL); as.numeric(scale(residuals(lm(x ~ ldepth)))) }
  rec_sur <- list(); lr <- list()
  for (k in seq_len(nrow(LR))) {
    L <- prep(gexp(LR$ligand[k])); R <- prep(gexp(LR$receptor[k]))
    nL <- sum(raw[intersect(if (LR$ligand[k] == "MHCII") MHCII_LIGAND else LR$ligand[k], rownames(raw)), ])
    nR <- if (LR$receptor[k] %in% rownames(raw)) sum(raw[LR$receptor[k], ]) else 0
    if (is.null(L) || is.null(R) || nR < 20) {
      lr[[k]] <- data.table(sample_id = s, pair = LR$pair[k], umi_ligand = nL, umi_receptor = nR,
                            S = NA_real_, z = NA_real_, p_pos = NA_real_, p_neg = NA_real_)
      next
    }
    rc <- LR$receptor[k]
    if (is.null(rec_sur[[rc]])) rec_sur[[rc]] <- scale(adespatial::msr(R, me, nrepet = N_MSR, method = "pair"))
    obs <- mean(L * as.numeric(Ws %*% R))
    nul <- colMeans(L * as.matrix(Ws %*% rec_sur[[rc]]))
    st <- msr_z(obs, nul)
    lr[[k]] <- data.table(sample_id = s, pair = LR$pair[k], umi_ligand = nL, umi_receptor = nR,
                          S = obs, z = st[["z"]], p_pos = st[["p_pos"]], p_neg = st[["p_neg"]])
  }

  res <- list(moran = moran, prox = rbindlist(prox), pb = rbindlist(pb), lr = rbindlist(lr),
              spots = data.table(sample_id = s, barcode = d$barcode, array_row = d$array_row,
                                 array_col = d$array_col, APC_index = apc, nonMHC_index = nonmhc,
                                 hot_C1b = h_c1b, hot_APC = h_apc, hot_nonMHC = h_non))
  saveRDS(res, obj_path(s, "apc09"))
  logf("%s: APC index I=%.3f (z vs random %.2f), APC hot=%d, Dice(APC,C1b)=%.2f, T_all niche log2FC=%.2f",
       s, I_apc, moran$z_APC_vs_random, sum(h_apc), moran$dice_APC_C1b,
       if (length(pb)) pb[[paste("C1b", "T_all")]]$log2_ratio %||% NA else NA)
  invisible(res)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- patch pass (PATCH09=1): recompute pseudobulk with mid-p and add
# TAM-adjusted ligand-receptor statistics, reusing the cached hotspots -------
patch_sample <- function(i) {
  s <- meta$sample_id[i]; seed <- CFG$seed + 900 + i
  set.seed(seed)
  res <- readRDS(obj_path(s, "apc09"))
  so <- readRDS(obj_path(s))
  raw <- read_raw_counts(s, colnames(so))
  raw_ln <- raw %*% Diagonal(x = 1e4 / Matrix::colSums(raw)); raw_ln@x <- log1p(raw_ln@x)
  dimnames(raw_ln) <- dimnames(raw)
  d <- spots[sample_id == s][match(colnames(so), barcode)]
  hs <- res$spots[match(d$barcode, barcode)]
  lw <- visium_listw(d)
  Ws <- listw_to_sparse(nb2listw(include.self(lw$neighbours), style = "W", zero.policy = TRUE))
  ldepth <- log10(d$nCount)

  tot <- Matrix::colSums(raw)
  gs <- lapply(TCELL_SETS, function(g) Matrix::colSums(raw[intersect(g, rownames(raw)), , drop = FALSE]))
  pb <- list()
  for (src in c("C1b", "APC")) {
    H <- which(if (src == "C1b") hs$hot_C1b else hs$hot_APC)
    if (length(H) < 5) next
    niche <- which(hop_distance(lw$neighbours, H) <= 2)
    regs <- random_regions(d, niche, N_REGION)
    for (ts in names(gs)) {
      obs <- sum(gs[[ts]][niche]) / sum(tot[niche])
      nul <- vapply(regs, function(r) sum(gs[[ts]][r]) / sum(tot[r]), 0)
      pb[[paste(src, ts)]] <- data.table(sample_id = s, niche = src, tset = ts, n_niche_spots = length(niche),
        umi_set_in_niche = sum(gs[[ts]][niche]), umi_set_total = sum(gs[[ts]]),
        cpm_niche = 1e6 * obs, cpm_null_mean = 1e6 * mean(nul), n_regions = length(regs),
        log2_ratio = log2((obs + 1e-7) / (mean(nul) + 1e-7)),
        z = if (sd(nul) > 0) (obs - mean(nul)) / sd(nul) else NA_real_,
        p_enrich = midp(obs, nul, "greater"), p_deplete = midp(obs, nul, "less"))
    }
  }
  res$pb <- rbindlist(pb)

  # TAM-adjusted LR: ligand and receptor residualised on log10 UMI + Mac MP
  need <- res$lr[is.finite(z)]
  lr_adj <- list()
  if (nrow(need)) {
    me <- adespatial::scores.listw(lw, MEM.autocor = "all")
    gexp <- function(g) {
      if (g == "MHCII") return(Matrix::colMeans(raw_ln[intersect(MHCII_LIGAND, rownames(raw_ln)), , drop = FALSE]))
      as.numeric(raw_ln[g, ])
    }
    prep <- function(x) as.numeric(scale(residuals(lm(x ~ ldepth + d$MP_Mac))))
    rec_sur <- list()
    for (k in seq_len(nrow(need))) {
      row <- LR[pair == need$pair[k]]
      L <- prep(gexp(row$ligand)); R <- prep(gexp(row$receptor))
      if (is.null(rec_sur[[row$receptor]]))
        rec_sur[[row$receptor]] <- scale(adespatial::msr(R, me, nrepet = N_MSR, method = "pair"))
      obs <- mean(L * as.numeric(Ws %*% R))
      nul <- colMeans(L * as.matrix(Ws %*% rec_sur[[row$receptor]]))
      st <- msr_z(obs, nul)
      lr_adj[[k]] <- data.table(sample_id = s, pair = row$pair, S = obs, z = st[["z"]],
                                p_pos = st[["p_pos"]], p_neg = st[["p_neg"]])
    }
  }
  res$lr_adj <- rbindlist(lr_adj)
  res$patched <- TRUE
  saveRDS(res, obj_path(s, "apc09"))
  logf("%s: 09 patch done (pseudobulk mid-p, %d TAM-adjusted LR pairs)", s, length(lr_adj))
}

for (i in seq_len(nrow(meta))) {
  s <- meta$sample_id[i]
  if (!s %in% run_ids) next
  if (Sys.getenv("PATCH09") == "1") {
    if (file.exists(obj_path(s, "apc09")) && !isTRUE(readRDS(obj_path(s, "apc09"))$patched)) patch_sample(i)
    next
  }
  if (file.exists(obj_path(s, "apc09")) && Sys.getenv("FORCE_APC") != "1") { logf("%s: cached", s); next }
  run_sample(i)
  patch_sample(i)          # mid-p pseudobulk + TAM-adjusted LR (same code path as PATCH09=1)
  gc(verbose = FALSE)
}

# ---- aggregate all cached samples ------------------------------------------------
fs <- obj_path(meta$sample_id, "apc09"); fs <- fs[file.exists(fs)]
if (length(fs) < nrow(meta)) logf("09: %d/%d samples cached; aggregating what is available", length(fs), nrow(meta))
R <- lapply(fs, readRDS)
mcol <- meta[, .(sample_id, cohort, histology, who_grade, is_LGG, patient)]
grab <- function(k) merge(rbindlist(lapply(R, `[[`, k), fill = TRUE), mcol, by = "sample_id")
moran <- grab("moran"); prox <- grab("prox"); pb <- grab("pb"); lr <- grab("lr")
lr_adj <- if (all(vapply(R, function(r) !is.null(r$lr_adj), TRUE))) grab("lr_adj") else NULL
saveRDS(rbindlist(lapply(R, `[[`, "spots")), file.path(P$objects, "apc09_spots.rds"))

prox[, `:=`(q_enrich = p.adjust(p_enrich, "BH"), q_deplete = p.adjust(p_deplete, "BH")), by = .(sample_id, source, zone)]
pb[, q_enrich := p.adjust(p_enrich, "BH"), by = .(sample_id, niche)]
lr[, q_pos := p.adjust(p_pos, "BH"), by = sample_id]
write_tab(moran, "apc_index_moran.csv")
write_tab(prox, "apc_vs_nonMHC_proximity.csv")
write_tab(pb, "niche_pseudobulk_Tcell.csv")
write_tab(lr, "lr_colocalization.csv")
if (!is.null(lr_adj)) {
  lr_adj[, q_pos := p.adjust(p_pos, "BH"), by = sample_id]
  write_tab(lr_adj, "lr_colocalization_TAMadjusted.csv")
}

cohorts <- list(IDHm = meta[cohort == "IDHm", sample_id], LGG = meta[is_LGG == TRUE, sample_id],
                GBM = meta[cohort == "GBM", sample_id])
stouffer <- function(z) sum(z) / sqrt(length(z))
pz <- function(p_enr) qnorm(pmin(pmax(1 - p_enr, 1e-12), 1 - 1e-12))   # one-sided p -> z

mm <- rbindlist(lapply(names(cohorts), function(ch) moran[sample_id %in% cohorts[[ch]], .(
  cohort = ch, k = .N, median_I_APC = median(I_APC_index), median_z_vs_random = median(z_APC_vs_random),
  n_p_lt_0.05 = sum(p_APC_vs_random < 0.05), stouffer_z = stouffer(pz(p_APC_vs_random)),
  median_rho_APC_nonMHC = median(rho_APC_nonMHC), median_dice_APC_C1b = median(dice_APC_C1b),
  n_no_APC_hotspot = sum(n_hot_APC < 5))]))
mm[, p := pnorm(-stouffer_z)]
write_tab(mm, "meta_apc_index_moran.csv")

pm <- rbindlist(lapply(names(cohorts), function(ch) prox[sample_id %in% cohorts[[ch]] & is.finite(z_msr),
  .(cohort = ch, k = .N, stouffer_z = stouffer(z_msr), n_enriched = sum(q_enrich < 0.05),
    n_depleted = sum(q_deplete < 0.05)), by = .(source, target, zone)]))
pm[, p := 2 * pnorm(-abs(stouffer_z))][, q := p.adjust(p, "BH"), by = .(cohort, source, zone)]
write_tab(pm, "meta_apc_vs_nonMHC_proximity.csv")

# only samples with >= 10 UMI of the gene set anywhere in the section are combined
pbm <- rbindlist(lapply(names(cohorts), function(ch) pb[sample_id %in% cohorts[[ch]] & umi_set_total >= 10,
  .(cohort = ch, k = .N, median_log2_ratio = median(log2_ratio), n_enriched = sum(p_enrich < 0.05),
    stouffer_z = stouffer(pz(p_enrich)), total_umi = sum(umi_set_total)), by = .(niche, tset)]))
pbm[, p := pnorm(-stouffer_z)][, q := p.adjust(p, "BH"), by = .(cohort, niche)]
write_tab(pbm, "meta_niche_pseudobulk_Tcell.csv")

lrm <- rbindlist(lapply(names(cohorts), function(ch) lr[sample_id %in% cohorts[[ch]] & is.finite(z),
  .(cohort = ch, k = .N, mean_z = mean(z), stouffer_z = stouffer(z), n_pos = sum(q_pos < 0.05)), by = pair]))
lrm[, p := 2 * pnorm(-abs(stouffer_z))][, q := p.adjust(p, "BH"), by = cohort]
write_tab(lrm, "meta_lr_colocalization.csv")
if (!is.null(lr_adj)) {
  lra <- rbindlist(lapply(names(cohorts), function(ch) lr_adj[sample_id %in% cohorts[[ch]] & is.finite(z),
    .(cohort = ch, k = .N, mean_z = mean(z), stouffer_z = stouffer(z), n_pos = sum(q_pos < 0.05)), by = pair]))
  lra[, p := 2 * pnorm(-abs(stouffer_z))][, q := p.adjust(p, "BH"), by = cohort]
  write_tab(lra, "meta_lr_colocalization_TAMadjusted.csv")
}
logf("09 aggregation done (%d samples)", length(fs))
