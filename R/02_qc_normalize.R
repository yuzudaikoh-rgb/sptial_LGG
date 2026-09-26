# 02 - Per-sample QC and normalisation --------------------------------------
source("R/utils/common.R")
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})
set.seed(CFG$seed)
meta <- read_sample_meta()
qcfg <- CFG$qc

read_positions <- function(prefix) {
  f <- file.path(P$geo_dir, paste0(prefix, "_tissue_positions_list.csv.gz"))
  pos <- fread(f, header = FALSE)
  if (is.character(pos$V2)) pos <- fread(f, header = TRUE)  # spaceranger >= 2 format
  setnames(pos, c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  pos
}

mad_bounds <- function(x, n) {
  m <- median(x); d <- mad(x)
  c(lo = m - n * d, hi = m + n * d)
}

qc_rows <- list(); spot_qc <- list()
for (i in seq_len(nrow(meta))) {
  s <- meta$sample_id[i]; pref <- meta$prefix[i]
  counts <- Read10X_h5(file.path(P$geo_dir, paste0(pref, "_filtered_feature_bc_matrix.h5")))
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  rownames(counts) <- make.unique(rownames(counts))
  pos <- read_positions(pref)
  pos <- pos[in_tissue == 1 & barcode %in% colnames(counts)]
  counts <- counts[, pos$barcode]

  n_umi  <- Matrix::colSums(counts)
  n_gene <- Matrix::colSums(counts > 0)
  mt     <- grep("^MT-", rownames(counts))
  pct_mt <- 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / pmax(n_umi, 1)

  b_umi  <- mad_bounds(log10(n_umi + 1), qcfg$n_mad)
  b_gene <- mad_bounds(log10(n_gene + 1), qcfg$n_mad)
  b_mt   <- mad_bounds(pct_mt, qcfg$n_mad)
  thr_umi  <- max(qcfg$min_umi, 10^b_umi[["lo"]])
  thr_gene <- max(qcfg$min_genes, 10^b_gene[["lo"]])
  thr_mt   <- min(qcfg$max_mt_pct, max(10, b_mt[["hi"]]))

  keep <- n_umi >= thr_umi & n_gene >= thr_gene & pct_mt <= thr_mt
  # remove spots left without any hexagonal neighbour (spatial statistics need >= 1)
  pos_k <- pos[keep]
  nb <- visium_neighbours(pos_k)
  iso <- pos_k$barcode[lengths(lapply(nb, function(z) z[z > 0])) == 0]
  keep[pos$barcode %in% iso] <- FALSE

  spot_qc[[s]] <- data.table(sample_id = s, barcode = pos$barcode, n_umi, n_gene, pct_mt,
                             pass = keep, array_row = pos$array_row, array_col = pos$array_col)

  counts_f <- counts[, keep]
  gkeep <- Matrix::rowSums(counts_f > 0) >= qcfg$min_spots_per_gene
  counts_f <- counts_f[gkeep, ]

  so <- CreateSeuratObject(counts_f, project = s, assay = "Spatial")
  md <- as.data.frame(pos[match(colnames(so), barcode)])
  rownames(md) <- md$barcode
  so <- AddMetaData(so, md[, c("array_row", "array_col", "pxl_row", "pxl_col")])
  so$pct_mt <- pct_mt[colnames(so)]
  so$sample_id <- s
  so <- NormalizeData(so, verbose = FALSE)                       # log(1 + CP10K)
  so <- SCTransform(so, assay = "Spatial", vst.flavor = "v2", verbose = FALSE)
  so <- RunPCA(so, verbose = FALSE, npcs = 30)
  so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
  so <- FindClusters(so, resolution = 0.8, verbose = FALSE, random.seed = CFG$seed)
  DefaultAssay(so) <- "Spatial"
  saveRDS(so, obj_path(s))

  qc_rows[[s]] <- data.table(sample_id = s, spots_in_tissue = length(keep),
                             spots_pass = sum(keep), spots_isolated_removed = length(iso),
                             genes_kept = sum(gkeep),
                             thr_umi = round(thr_umi), thr_gene = round(thr_gene),
                             thr_mt = round(thr_mt, 1),
                             median_umi = median(n_umi[keep]), median_genes = median(n_gene[keep]),
                             median_pct_mt = round(median(pct_mt[keep]), 2))
  logf("%s: %d/%d spots pass, %d genes", s, sum(keep), length(keep), sum(gkeep))
  rm(so, counts, counts_f); gc(verbose = FALSE)
}

qc <- merge(meta[, .(sample_id, cohort, histology, who_grade)], rbindlist(qc_rows), by = "sample_id")
write_tab(qc, "qc_summary.csv")
sq <- merge(rbindlist(spot_qc), meta[, .(sample_id, cohort)], by = "sample_id")
fwrite(sq, file.path(P$objects, "spot_qc_all.csv.gz"))

# ---- QC figures -----------------------------------------------------------
ord <- meta[order(cohort, who_grade, sample_id), sample_id]
sq[, sample_id := factor(sample_id, ord)]
p1 <- ggplot(sq, aes(sample_id, log10(n_umi), fill = cohort)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_point(data = qc[, .(sample_id = factor(sample_id, ord), y = log10(thr_umi))],
             aes(sample_id, y), inherit.aes = FALSE, shape = 95, size = 5, colour = "red") +
  labs(x = NULL, y = "log10 UMI / spot") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- ggplot(sq, aes(sample_id, pct_mt, fill = cohort)) +
  geom_violin(scale = "width", linewidth = 0.2) + labs(x = NULL, y = "% mitochondrial") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- ggplot(qc, aes(factor(sample_id, ord), spots_pass, fill = cohort)) + geom_col() +
  labs(x = NULL, y = "Spots passing QC") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p1 / p2 / p3 + plot_layout(guides = "collect"), "FigS1_QC.pdf", 7.5, 7)

maps <- lapply(ord, function(s) {
  d <- sq[sample_id == s]; d[, log_umi := log10(n_umi)]
  d[pass == FALSE, log_umi := NA]
  spot_plot(d, "log_umi", s, pal = "viridis", size = 0.3)
})
save_pdf(wrap_plots(maps, ncol = 5), "FigS1b_QC_spatial_logUMI.pdf", 10, 8)
