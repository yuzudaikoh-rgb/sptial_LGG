# Run the full Consensus1b spatial pipeline (from the repository root):
#   Rscript run_all.R            # all steps
#   Rscript run_all.R 04 05      # selected steps
steps <- c("00_download", "01_metadata", "02_qc_normalize", "03_score_c1b",
           "04_spatial_autocorr", "05_colocalization", "06_cna_tumor_regions", "07_integrate_stats", "08_c1b_vs_3gene_niche", "08b_plot_niche",
           "09_apc_tam_niche", "09b_plot_apc")
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) steps <- steps[substr(steps, 1, 2) %in% args | substr(steps, 1, 3) %in% args]
for (s in steps) {
  message("==== ", s, " ====")
  t0 <- Sys.time()
  status <- system2("Rscript", file.path("R", paste0(s, ".R")))
  if (status != 0) stop("step failed: ", s)
  message(sprintf("==== %s done in %.1f min ====", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
