# 11 - Independent replication cohort (Hoefflin et al., Cancer Cell 2026) -----
# Run with:  C1B_CONFIG=config/config_replication.yaml Rscript R/11_replication_setup.R
# then steps 02-09 with the same C1B_CONFIG.
#
#  * downloads the IDH-mutant Visium archive (Zenodo 18380571) and the sample
#    metadata table (Zenodo 21335411, Table S1)
#  * removes sections already contained in GSE237183 (identical in-tissue
#    barcode set and identical total UMI) and the low-quality "LQ" section
#  * lays the remaining sections out in the file naming used by step 02 and
#    writes results/replication/tables/sample_metadata.csv
source("R/utils/common.R")
suppressPackageStartupMessages({ library(Seurat); library(readxl) })
stopifnot(grepl("replication", CFG_FILE))

raw <- P$raw
tgz <- file.path(raw, "Visium_IDHmut_samples.tar.gz")
if (!file.exists(tgz))
  system2("curl", c("-sL", "--retry", "5", "-o", shQuote(tgz),
                    shQuote("https://zenodo.org/records/18380571/files/Visium_IDHmut_samples.tar.gz?download=1")))
xl <- file.path(raw, "TableS1_SampleMetadata_v13.xlsx")
if (!file.exists(xl))
  system2("curl", c("-sL", "--retry", "5", "-o", shQuote(xl),
                    shQuote("https://zenodo.org/records/21335411/files/TableS1_SampleMetadata_v13.xlsx?download=1")))
src <- file.path(raw, "IDHmut")
if (!dir.exists(src)) {
  system2("tar", c("-xzf", shQuote(tgz), "-C", shQuote(raw), "--wildcards",
                   "'*/outs/filtered_feature_bc_matrix.h5'", "'*/outs/spatial/tissue_positions*'",
                   "'*/outs/spatial/scalefactors_json.json'"))
  inner <- list.files(raw, pattern = "^IDHmut$", recursive = TRUE, include.dirs = TRUE, full.names = TRUE)
  file.rename(inner[1], src); unlink(file.path(raw, "home"), recursive = TRUE)
}

# ---- de-duplicate against GSE237183 ---------------------------------------------
gse_meta <- fread(file.path(ROOT, "results", "tables", "sample_metadata.csv"))[cohort == "IDHm"]
gse_geo  <- file.path(ROOT, "data", "raw", "geo")
sig <- function(h5) { x <- Read10X_h5(h5); if (is.list(x)) x <- x[[1]]; list(bc = colnames(x), tot = sum(x)) }
G <- lapply(gse_meta$prefix, function(p) sig(file.path(gse_geo, paste0(p, "_filtered_feature_bc_matrix.h5"))))
names(G) <- gse_meta$sample_id
ids <- setdiff(basename(list.dirs(src, recursive = FALSE)), "LQ")
dup <- rbindlist(lapply(ids, function(s) {
  q <- sig(file.path(src, s, "outs", "filtered_feature_bc_matrix.h5"))
  hit <- names(G)[vapply(G, function(g) setequal(g$bc, q$bc) && isTRUE(all.equal(g$tot, q$tot)), TRUE)]
  data.table(sample_id = s, duplicate_of_GSE237183 = if (length(hit)) hit else NA_character_)
}))
write_tab(dup, "replication_dedup_vs_GSE237183.csv")
keep <- dup[is.na(duplicate_of_GSE237183), sample_id]
logf("replication: %d sections, %d duplicates of GSE237183 removed, %d kept", nrow(dup),
     sum(!is.na(dup$duplicate_of_GSE237183)), length(keep))

# ---- file layout expected by step 02 --------------------------------------------
dir.create(P$geo_dir, showWarnings = FALSE, recursive = TRUE)
for (s in keep) {
  pref <- paste0("HOEF_", tolower(s))
  h5 <- file.path(P$geo_dir, paste0(pref, "_filtered_feature_bc_matrix.h5"))
  if (!file.exists(h5)) file.symlink(normalizePath(file.path(src, s, "outs", "filtered_feature_bc_matrix.h5")), h5)
  pos <- list.files(file.path(src, s, "outs", "spatial"), pattern = "^tissue_positions", full.names = TRUE)[1]
  x <- fread(pos, header = FALSE)
  if (!is.numeric(x[[2]])) x <- fread(pos, header = TRUE)      # Space Ranger >= 2 has a header
  fwrite(x, file.path(P$geo_dir, paste0(pref, "_tissue_positions_list.csv.gz")), col.names = FALSE)
}

# ---- metadata -----------------------------------------------------------------------
tm <- as.data.table(read_excel(xl, sheet = "This_cohort"))
tm <- tm[SampleID %in% keep]
meta <- data.table(sample_id = tolower(tm$SampleID), gsm = NA_character_, patient = tm$PatientID,
                   cohort = "IDHm", idh_status = "IDH-mutant",
                   histology = ifelse(tm$GliomaType == "IDH_O", "Oligodendroglioma", "Astrocytoma"),
                   who_grade = as.integer(tm$GradeTumor), section_grade = as.integer(tm$GradeSection),
                   is_LGG = as.integer(tm$GradeTumor) %in% 2:3, tumor_region = tm$Location,
                   mgmt_status = NA_character_, author_id = NA_character_,
                   title = tm$SampleID, prefix = paste0("HOEF_", tolower(tm$SampleID)),
                   sex = tm$Sex, age = tm$Age)
write_tab(meta, "sample_metadata.csv")
write_tab(meta[, .N, by = .(histology, who_grade)][order(histology, who_grade)], "sample_counts.csv")
logf("replication cohort: %d sections from %d patients (grade 2: %d, 3: %d, 4: %d)", nrow(meta),
     uniqueN(meta$patient), sum(meta$who_grade == 2), sum(meta$who_grade == 3), sum(meta$who_grade == 4))
