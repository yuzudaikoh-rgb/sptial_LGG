# CNA inference helpers ------------------------------------------------------
# Port of the spatial CNA approach of Greenwald et al. (Cell 2024;
# tiroshlab/Spatial_Glioma, Module5_CNA.R: calc_cna / cna_sig_cor), rewritten
# for sparse input and with an optional genomic-order permutation used to
# build a per-sample null distribution for malignancy calls.
suppressPackageStartupMessages({ library(Matrix) })

# Column-wise centred running mean with shrinking windows at the ends
# (identical to caTools::runmean(x, k, endrule = "mean", align = "center"),
# but vectorised over all columns via cumulative sums).
runmean_mat <- function(x, k) {
  n <- nrow(x); k <- min(k, n); k2 <- k %/% 2
  i <- seq_len(n)
  lo <- pmax(1L, i - k2 + (k %% 2 == 0))   # even k: window [i - k/2 + 1, i + k/2]
  hi <- pmin(n, i + k2)
  cs <- rbind(0, apply(x, 2, cumsum))
  (cs[hi + 1, , drop = FALSE] - cs[lo, , drop = FALSE]) / (hi - lo + 1)
}

log_norm <- function(x) log2(1 + x / 10)
un_log   <- function(x) 10 * (2^x - 1)

# genes x spots counts (sparse) -> log2(1 + CPM/10), genes with mean > min_mean
cna_logcpm <- function(counts, min_mean = 0.4) {
  cpm <- counts %*% Diagonal(x = 1e6 / Matrix::colSums(counts))
  cpm@x <- log2(1 + cpm@x / 10)
  mu <- Matrix::rowMeans(cpm)
  keep <- mu > min_mean
  m <- as.matrix(cpm[keep, ])
  colnames(m) <- colnames(counts)
  m[apply(m, 1, var) > 0, ]
}

genome_order <- function(genes, genome) {
  g <- as.data.frame(genome)
  g <- g[g$symbol %in% genes & !duplicated(g$symbol), ]
  g$chr <- as.integer(as.character(factor(g$chromosome_name, levels = c(1:22, "X", "Y"),
                                          labels = 1:24)))
  g <- g[!is.na(g$chr) & g$chr <= 22, ]               # autosomes only
  g <- g[order(g$chr, g$start_position), ]
  g
}

# Core CNA score matrix (genes ordered by genome x spots)
calc_cna <- function(m, query, ref, genome, window = 150, range = c(-3, 3), noise = 0.15,
                     permute = FALSE, seed = 1) {
  go <- genome_order(rownames(m), genome)
  x <- m[go$symbol, c(query, unlist(ref)), drop = FALSE]
  chr <- go$chr
  if (permute) { set.seed(seed); x <- x[sample(nrow(x)), ]; rownames(x) <- go$symbol }
  x <- x - rowMeans(x)
  x[x > range[2]] <- range[2]; x[x < range[1]] <- range[1]
  x <- un_log(x)
  out <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
  for (ch in unique(chr)) {
    idx <- which(chr == ch)
    k <- min(window, length(idx))
    out[idx, ] <- runmean_mat(x[idx, , drop = FALSE], k)
  }
  out <- log_norm(out)
  out <- sweep(out, 2, apply(out, 2, median))
  out <- un_log(out)
  if (!is.list(ref)) ref <- list(ref)
  ref_mat <- sapply(ref, function(r) rowMeans(out[, r, drop = FALSE]))
  ref_mat <- log_norm(ref_mat)
  rmax <- apply(ref_mat, 1, max) + noise
  rmin <- apply(ref_mat, 1, min) - noise
  out <- log_norm(out)
  sc <- pmax(out - rmax, 0) + pmin(out - rmin, 0)   # rmax > rmin, so at most one term is non-zero
  dimnames(sc) <- dimnames(out)
  attr(sc, "chr") <- chr
  attr(sc, "arm") <- as.character(go$arm)
  sc
}

# CNA signal (mean |CNA| over top-1/3 CNA regions) and CNA correlation
# (mean correlation with the top-1/3 CNA-signal query spots)
cna_sig_cor <- function(cna, query, top_region = 1/3, top_cells = 1/3) {
  q <- cna[, query, drop = FALSE]
  absvals <- rowMeans(abs(q))
  top <- absvals > quantile(absvals, 1 - top_region)
  sig <- colMeans(abs(cna[top, , drop = FALSE]))
  cut <- quantile(sig[query], 1 - top_cells)
  # as in the original: mean correlation with each of the top-signal query spots
  cr <- rowMeans(cor(cna[top, , drop = FALSE], q[top, sig[query] > cut, drop = FALSE]))
  names(cr) <- colnames(cna)
  list(sig = sig, cor = cr)
}

# Mean CNA per chromosome arm for each spot
arm_means <- function(cna) {
  arm <- attr(cna, "arm")
  a <- rowsum(cna, arm) / as.numeric(table(arm)[rownames(rowsum(cna, arm))])
  a
}
