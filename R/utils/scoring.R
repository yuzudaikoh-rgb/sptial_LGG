# Signature scoring helpers --------------------------------------------------
suppressPackageStartupMessages({ library(UCell); library(Matrix) })

# UCell on a genes x spots count matrix; returns spots x signatures
ucell_scores <- function(counts, sigs, maxRank = 1500, ncores = 1, min_genes = 5) {
  sigs <- lapply(sigs, function(g) intersect(g, rownames(counts)))
  sigs <- sigs[lengths(sigs) >= min_genes]
  sc <- UCell::ScoreSignatures_UCell(counts, features = sigs, maxRank = maxRank,
                                     ncores = ncores, name = "")
  sc[, names(sigs), drop = FALSE]
}

# Author-style expression matrix (Greenwald et al. Module 1/5):
# log2(1 + CPM/10), MT/ribosomal genes removed, genes with mean > 0.4 kept
author_logcpm <- function(counts, drop_mt_ribo = TRUE, min_mean = 0.4) {
  m <- counts
  if (drop_mt_ribo) m <- m[!grepl("^MT-|^RPL|^RPS", rownames(m)), ]
  cpm <- t(t(m) / Matrix::colSums(m)) * 1e6
  lg <- log2(1 + as.matrix(cpm) / 10)
  lg <- lg[apply(lg, 1, var) > 0, ]
  lg[rowMeans(lg) > min_mean, ]
}

# Tirosh et al. 2016 signature score (as in scalop::sigScores with
# expr.center = TRUE): per-gene centred expression, signature mean minus the
# mean of expression-bin-matched control genes (30 bins, 100 ctrl / gene).
tirosh_scores <- function(lg, sigs, nbin = 30, nctrl = 100, conserved = 0.5, seed = 1) {
  set.seed(seed)
  avg <- rowMeans(lg)
  bins <- cut(rank(avg, ties.method = "first"), nbin, labels = FALSE)
  names(bins) <- rownames(lg)
  cen <- lg - avg
  out <- sapply(sigs, function(g) {
    g0 <- g; g <- intersect(g, rownames(lg))
    if (length(g) < conserved * length(g0)) return(rep(NA_real_, ncol(lg)))
    ctrl <- unique(unlist(lapply(g, function(x) {
      pool <- names(bins)[bins == bins[x]]
      sample(pool, min(nctrl, length(pool)))
    })))
    colMeans(cen[g, , drop = FALSE]) - colMeans(cen[ctrl, , drop = FALSE])
  })
  rownames(out) <- colnames(lg)
  out
}

# Expression/detection-matched random gene sets. For every signature gene a
# non-signature gene is drawn from the same bin of average expression x
# detection rate (n_bins quantile bins of each, crossed).
matched_random_sets <- function(expr, sig, n_sets = 1000, n_bins = 25, seed = 1) {
  set.seed(seed)
  sig <- intersect(sig, rownames(expr))
  mu  <- Matrix::rowMeans(expr)
  det <- Matrix::rowMeans(expr > 0)
  b_mu  <- cut(rank(mu, ties.method = "first"), n_bins, labels = FALSE)
  b_det <- cut(rank(det, ties.method = "first"), 5, labels = FALSE)
  key <- paste(b_mu, b_det)
  names(key) <- rownames(expr)
  pool <- split(setdiff(rownames(expr), sig), key[setdiff(rownames(expr), sig)])
  # fall back to expression bin only if the crossed bin is empty
  pool_mu <- split(setdiff(rownames(expr), sig), b_mu[!rownames(expr) %in% sig])
  names(b_mu) <- rownames(expr)
  sets <- lapply(seq_len(n_sets), function(k) {
    vapply(sig, function(g) {
      cand <- pool[[key[g]]]
      if (length(cand) < 5) cand <- pool_mu[[as.character(b_mu[g])]]
      sample(cand, 1)
    }, character(1))
  })
  names(sets) <- sprintf("rand%04d", seq_len(n_sets))
  lapply(sets, unique)
}
