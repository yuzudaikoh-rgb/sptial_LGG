# 00 - Download GEO GSE237183 supplementary files + author inputs -----------
# Run from the repository root:  Rscript R/00_download.R
# Idempotent: files already present are not downloaded again.
source("R/utils/common.R")

fetch <- function(url, dest, tries = 5) {
  if (file.exists(dest) && file.size(dest) > 0) {
    logf("present: %s", basename(dest)); return(invisible(dest))
  }
  for (k in seq_len(tries)) {
    st <- system2("curl", c("-sSL", "--retry", "3", "-C", "-", "-o", shQuote(dest), shQuote(url)))
    if (st == 0 && file.exists(dest) && file.size(dest) > 0) break
    Sys.sleep(2^k)
  }
  if (!file.exists(dest)) stop("download failed: ", url)
  logf("downloaded %s (%.1f MB)", basename(dest), file.size(dest) / 1e6)
  invisible(dest)
}

raw <- P$raw
tar_f <- fetch(CFG$geo$raw_tar_url, file.path(raw, "GSE237183_RAW.tar"))
sm_f  <- fetch(CFG$geo$series_matrix_url, file.path(raw, "GSE237183_series_matrix.txt.gz"))
inp_f <- fetch(sprintf("https://drive.usercontent.google.com/download?id=%s&export=download&confirm=t",
                       CFG$author_inputs$gdrive_id), file.path(raw, "Inputs.zip"))

# md5 log for provenance
md5 <- tools::md5sum(c(tar_f, sm_f, inp_f))
write_tab(data.table(file = basename(names(md5)), md5 = unname(md5),
                     bytes = file.size(names(md5))), "download_md5.csv")

# Untar GEO per-sample files
if (length(list.files(P$geo_dir, pattern = "h5$")) < 19) {
  utils::untar(tar_f, exdir = P$geo_dir)
}
logf("GEO files: %d", length(list.files(P$geo_dir)))

# Extract only the author inputs needed (the zip is ~6 GB, mostly NMF/CODEX)
need <- c("Inputs/CNA/*", "Inputs/MP/clean_spatial_gbm_metaprograms_124.rds",
          "Inputs/MP/mp_assign_124/*", "Inputs/general/GBM_samples.txt",
          "Inputs/general/visium_metadata.csv",
          "Inputs/Spatial_coh_zones/spatial_zonesv3.rds")
if (!file.exists(file.path(P$inputs_dir, "MP", "clean_spatial_gbm_metaprograms_124.rds"))) {
  system2("unzip", c("-q", "-o", shQuote(inp_f), shQuote(need), "-d", shQuote(raw)))
}
logf("author inputs ready: %s", P$inputs_dir)
