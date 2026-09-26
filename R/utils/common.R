# Common helpers: configuration, logging, gene sets, sample I/O -------------

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(Matrix)
})

find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "config", "config.yaml"))) {
    p <- dirname(d)
    if (identical(p, d)) stop("Cannot find project root (config/config.yaml)")
    d <- p
  }
  d
}

ROOT <- find_root()
CFG  <- yaml::read_yaml(file.path(ROOT, "config", "config.yaml"))
P    <- lapply(CFG$paths, function(x) file.path(ROOT, x))
invisible(lapply(P, dir.create, recursive = TRUE, showWarnings = FALSE))

logf <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...))
  message(msg)
  cat(msg, "\n", file = file.path(P$logs, "pipeline.log"), append = TRUE)
}

write_tab <- function(x, name) {
  f <- file.path(P$tables, name)
  data.table::fwrite(as.data.table(x), f)
  logf("wrote %s (%d rows)", name, nrow(x))
  invisible(f)
}

# ---- Consensus1b signature (90 genes, as supplied) ------------------------
C1B_GENES <- c(
  "LYZ", "VIM", "AREG", "HLA-DPB1", "CD163", "LGALS3", "HLA-DQA1",
  "S100A4", "HLA-DPA1", "S100A6", "HLA-DRB1", "ANXA1", "HLA-DQB1",
  "HLA-DRA", "SLC40A1", "MS4A6A", "TGFBI", "LGALS1", "CD74", "MRC1",
  "CST3", "ARL4C", "ANXA2", "CTSD", "S100A10", "MS4A4A", "FGL2",
  "CD44", "AHR", "AHNAK", "CSTA", "CD48", "IQGAP1", "NAMPT", "FCHO2",
  "MPEG1", "EMP3", "HLA-DMA", "FILIP1L", "SH3BP5", "FXYD5", "GRN",
  "LITAF", "CTSS", "CFD", "STAB1", "IFITM3", "ARHGAP18", "NPC2", "DAB2",
  "TXNIP", "TXN", "BTG1", "FNBP1", "AP1S2", "PSAP", "CD68", "BLVRB",
  "RNF144B", "YBX1", "C1orf162", "CHPT1", "TMSB10", "GPR183", "BRI3",
  "RPS18", "HEXA", "FCGRT", "RNASE6", "TNFAIP3", "SDCBP", "HLA-DMB",
  "PPDPF", "RPS2", "RPS13", "HLA-DOA", "GSTP1", "PTMA", "PHLDA1",
  "CARD16", "PPA1", "FTH1", "MARCH1", "PECAM1", "RPL12", "PCBD1",
  "RPS17", "SNX6", "RBPJ", "DRAM2"
)
stopifnot(length(C1B_GENES) == 90, !anyDuplicated(C1B_GENES))

# HGNC renamings that differ between annotation versions (old -> new)
GENE_ALIASES <- c("MARCH1" = "MARCHF1")

# Map signature symbols onto the features of a matrix (tries alias both ways)
map_genes <- function(genes, features) {
  out <- vapply(genes, function(g) {
    if (g %in% features) return(g)
    a <- GENE_ALIASES[g]
    if (!is.na(a) && a %in% features) return(unname(a))
    rev <- names(GENE_ALIASES)[GENE_ALIASES == g]
    if (length(rev) && rev[1] %in% features) return(rev[1])
    NA_character_
  }, character(1))
  out
}

# Signature variants used for sensitivity analyses
c1b_variants <- function() {
  ribo <- grep("^RP[SL]", C1B_GENES, value = TRUE)
  mhc2 <- grep("^HLA-D|^CD74$", C1B_GENES, value = TRUE)
  list(
    C1b          = C1B_GENES,
    C1b_noRibo   = setdiff(C1B_GENES, ribo),
    C1b_MHCII    = mhc2,
    C1b_nonMHCII = setdiff(C1B_GENES, c(mhc2, ribo))
  )
}

# ---- Sample metadata ------------------------------------------------------
read_sample_meta <- function() {
  f <- file.path(P$tables, "sample_metadata.csv")
  if (!file.exists(f)) stop("Run R/01_metadata.R first")
  fread(f)
}

obj_path <- function(sample_id, what = "seurat") {
  file.path(P$objects, sprintf("%s_%s.rds", sample_id, what))
}

# Seurat object per sample only for the samples that passed through 02_
available_samples <- function(meta = read_sample_meta()) {
  meta[file.exists(obj_path(meta$sample_id)), ]
}

# ---- Visium grid neighbour graph ------------------------------------------
# Visium spots sit on a hexagonal lattice; in array coordinates the 6
# nearest neighbours are (row +-0, col +-2) and (row +-1, col +-1).
visium_neighbours <- function(coords) {
  suppressPackageStartupMessages(library(spdep))
  xy <- cbind(coords$array_col, coords$array_row * sqrt(3))
  nb <- spdep::dnearneigh(xy, 0, 2.05)
  nb
}

visium_listw <- function(coords, style = "W") {
  nb <- visium_neighbours(coords)
  spdep::nb2listw(nb, style = style, zero.policy = TRUE)
}

# sparse row-standardised W for fast vectorised Moran's I
listw_to_sparse <- function(lw) {
  n <- length(lw$neighbours)
  i <- rep(seq_len(n), lengths(lw$neighbours))
  j <- unlist(lw$neighbours)
  keep <- j > 0
  x <- unlist(lw$weights)
  sparseMatrix(i = i[keep], j = j[keep], x = x[keep], dims = c(n, n))
}

# Moran's I for every column of X (spots x features), W row-standardised
moran_matrix <- function(X, W) {
  X <- as.matrix(X)
  Z <- sweep(X, 2, colMeans(X))
  n <- nrow(Z)
  S0 <- sum(W)
  num <- colSums(Z * as.matrix(W %*% Z))
  den <- colSums(Z^2)
  (n / S0) * num / den
}

# ---- Plot theme -----------------------------------------------------------
theme_pub <- function(base = 8) {
  ggplot2::theme_classic(base_size = base) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(colour = "black"),
      legend.key.size = grid::unit(0.35, "cm")
    )
}

save_pdf <- function(p, name, w, h) {
  f <- file.path(P$figures, name)
  ggplot2::ggsave(f, p, width = w, height = h, units = "in", device = grDevices::cairo_pdf)
  png_f <- sub("\\.pdf$", ".png", f)
  ggplot2::ggsave(png_f, p, width = w, height = h, units = "in", dpi = 200)
  logf("saved figure %s", name)
}

# Spatial spot plot on array coordinates (hex grid)
spot_plot <- function(df, value, title = NULL, limits = NULL, pal = "magma",
                      discrete = FALSE, size = 0.55) {
  g <- ggplot2::ggplot(df, ggplot2::aes(x = array_col, y = -array_row,
                                        colour = .data[[value]])) +
    ggplot2::geom_point(size = size, shape = 16) +
    ggplot2::coord_fixed(ratio = sqrt(3)) +
    ggplot2::theme_void(base_size = 7) +
    ggplot2::labs(title = title, colour = NULL) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 7))
  if (discrete) {
    g
  } else {
    g + viridis::scale_colour_viridis(option = pal, limits = limits, oob = scales::squish)
  }
}
