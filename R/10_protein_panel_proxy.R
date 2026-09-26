# 10 - Which CODEX-measurable protein panel can stand in for Consensus1b? ----
#  CODEX (Hoefflin et al., Cancer Cell 2026) measures five C1b members as
#  proteins: MHC-II (pan HLA-DR/DP/DQ), CD163, CD206 (MRC1), CD44, VIM.
#  Before using them to define "C1b-like TAMs" at single-cell resolution we
#  test, in the Visium data, whether the transcript-level equivalents of these
#  panels reproduce the full 90-gene program:
#     P5  = MHCII + CD163 + MRC1 + CD44 + VIM
#     P3m = MHCII + CD163 + MRC1              (myeloid-restricted)
#     P3  = CD68 + CD14 + HLA-DPB1           (original 3-gene panel, reference)
#  Metrics: Spearman rho with C1b, Moran's I, Gi* hotspot Dice, top-10% Jaccard.
#  Nulls: 1,000 random same-size C1b subsets and 1,000 expression-matched sets.
source("R/utils/common.R")
source("R/utils/scoring.R")
source("R/utils/niche.R")
suppressPackageStartupMessages({ library(Seurat); library(spdep); library(ggplot2); library(patchwork) })

N_NULL <- 1000
MHCII  <- c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1")
PANELS <- list(P5  = list(MHCII = MHCII, CD163 = "CD163", MRC1 = "MRC1", CD44 = "CD44", VIM = "VIM"),
               P3m = list(MHCII = MHCII, CD163 = "CD163", MRC1 = "MRC1"),
               P3  = list(CD68 = "CD68", CD14 = "CD14", HLADPB1 = "HLA-DPB1"))

# panel score: each marker (MHC-II = mean of its genes) z-scored, then averaged,
# mirroring how per-marker protein intensities are combined on CODEX
panel_score <- function(lognorm, panel) {
  f <- vapply(panel, function(g) {
    g <- intersect(g, rownames(lognorm)); length(g) > 0 }, TRUE)
  m <- sapply(panel[f], function(g) {
    g <- intersect(g, rownames(lognorm))
    as.numeric(Matrix::colMeans(lognorm[g, , drop = FALSE]))
  })
  m <- scale(m); m[is.na(m)] <- 0
  rowMeans(m)
}

meta  <- available_samples()
spots <- readRDS(file.path(P$objects, "spot_scores.rds"))
hot   <- readRDS(file.path(P$objects, "c1b_hotspots.rds"))
rows <- list(); nulls <- list()
done_f <- file.path(P$tables, "protein_panel_proxy_per_sample.csv")
resume <- file.exists(done_f) && Sys.getenv("FORCE10") != "1"
for (i in if (resume) integer() else seq_len(nrow(meta))) {
  s <- meta$sample_id[i]; set.seed(CFG$seed + 1000 + i)
  so <- readRDS(obj_path(s)); lognorm <- GetAssayData(so, assay = "Spatial", layer = "data")
  d <- spots[sample_id == s][match(colnames(so), barcode)]
  hc <- hot[sample_id == s][match(d$barcode, barcode), hotspot] == "hot"
  lw <- visium_listw(d); W <- listw_to_sparse(lw)
  lw_b <- nb2listw(include.self(lw$neighbours), style = "B", zero.policy = TRUE)
  top_c <- d$C1b >= quantile(d$C1b, .9)
  c1b <- unname(na.omit(map_genes(C1B_GENES, rownames(lognorm))))
  metr <- function(x) {
    h <- gi_hot(x, lw_b); tp <- x >= quantile(x, .9)
    c(rho = cor(x, d$C1b, method = "spearman"), moran = unname(moran_matrix(cbind(x), W)),
      moran_C1b = unname(moran_matrix(cbind(d$C1b), W)), n_hot = sum(h), dice = dice(h, hc), jaccard10 = sum(tp & top_c) / sum(tp | top_c))
  }
  for (pn in names(PANELS)) {
    pan <- PANELS[[pn]]; k <- length(pan)
    x <- panel_score(lognorm, pan)
    m <- metr(x)
    # nulls (rho only; cheap): random C1b subsets of k genes, expression-matched k-gene sets
    anchor <- vapply(pan, function(g) intersect(g, rownames(lognorm))[1], "")
    rin  <- vapply(seq_len(N_NULL), function(b) cor(zmean(lognorm, sample(setdiff(c1b, unlist(pan)), k)),
                                                    d$C1b, method = "spearman"), 0)
    rs   <- matched_random_sets(lognorm, anchor[!is.na(anchor)], n_sets = N_NULL, seed = CFG$seed + i)
    rgen <- vapply(rs, function(g) cor(zmean(lognorm, g), d$C1b, method = "spearman"), 0)
    rows[[paste(s, pn)]] <- data.table(sample_id = s, panel = pn, n_markers = k, t(m),
      pct_C1b_subsets_below = mean(rin < m[["rho"]]), pct_genome_below = mean(rgen < m[["rho"]]),
      median_rho_C1b_subsets = median(rin))
    nulls[[paste(s, pn)]] <- data.table(sample_id = s, panel = pn, null = rep(c("C1b_subset", "genome"), each = N_NULL),
                                        rho = c(rin, rgen))
  }
  logf("%s: rho(C1b) P5=%.2f P3m=%.2f P3=%.2f", s, rows[[paste(s, "P5")]]$rho,
       rows[[paste(s, "P3m")]]$rho, rows[[paste(s, "P3")]]$rho)
}
if (resume) {
  res <- fread(done_f)
  if ("moran.x" %in% names(res)) setnames(res, "moran.x", "moran")
  if (!"moran_C1b" %in% names(res)) {               # older runs: add the reference Moran's I
    mc <- rbindlist(lapply(unique(res$sample_id), function(s) {
      d <- spots[sample_id == s]
      data.table(sample_id = s, moran_C1b = moran_matrix(cbind(d$C1b), listw_to_sparse(visium_listw(d))))
    }))
    res <- merge(res, mc, by = "sample_id")
  }
} else {
  res <- merge(rbindlist(rows), meta[, .(sample_id, cohort, who_grade, is_LGG)], by = "sample_id")
  nl <- rbindlist(nulls); fwrite(nl, file.path(P$objects, "protein_panel_null.csv.gz"))
}
write_tab(res, "protein_panel_proxy_per_sample.csv")
sm <- res[, .(k = .N, median_rho = median(rho), min_rho = min(rho), median_pct_C1b_subsets = median(pct_C1b_subsets_below),
              n_beats_95pct_subsets = sum(pct_C1b_subsets_below > .95), median_moran_ratio = median(moran / moran_C1b),
              n_no_hotspot = sum(n_hot < 5), median_dice = median(dice), median_jaccard10 = median(jaccard10)),
          by = .(cohort, panel)]
write_tab(sm, "protein_panel_proxy_summary.csv")
# ---- figure -----------------------------------------------------------------------
ord <- meta[order(cohort == "GBM", who_grade, sample_id), sample_id]
res[, sample_id := factor(sample_id, ord)]
lab <- c(P5 = "MHCII/CD163/CD206/CD44/VIM", P3m = "MHCII/CD163/CD206", P3 = "CD68/CD14/HLA-DPB1")
g1 <- ggplot(res, aes(sample_id, rho, colour = panel, group = panel)) +
  geom_line(alpha = .5) + geom_point(size = 1.8) +
  scale_colour_manual(values = c(P5 = "#B2182B", P3m = "#E08214", P3 = "grey50"), labels = lab) +
  labs(x = NULL, y = "Spearman rho with 90-gene C1b", colour = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
m2 <- melt(res[, .(sample_id, cohort, panel, `Moran I / C1b Moran I` = moran / moran_C1b,
                   `Hotspot Dice` = dice, `Top-10% Jaccard` = jaccard10,
                   `Pct random C1b subsets beaten` = pct_C1b_subsets_below)], id.vars = c("sample_id", "cohort", "panel"))
g2 <- ggplot(m2, aes(panel, value, fill = panel)) + geom_boxplot(outlier.shape = NA, alpha = .6) +
  geom_jitter(aes(colour = cohort), width = .15, size = .8) + facet_wrap(~variable, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(P5 = "#B2182B", P3m = "#E08214", P3 = "grey50"), guide = "none") +
  scale_colour_manual(values = c(IDHm = "#D6604D", GBM = "#4393C3")) + labs(x = NULL, y = NULL) + theme_pub()
save_pdf(g1 / g2, "Fig9a_protein_panel_proxy.pdf", 10, 6)
