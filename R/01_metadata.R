# 01 - Sample metadata from the GEO series matrix ---------------------------
source("R/utils/common.R")

sm <- readLines(gzfile(file.path(P$raw, "GSE237183_series_matrix.txt.gz")))
grab <- function(key) {
  rows <- grep(paste0("^!", key, "\t"), sm, value = TRUE)
  lapply(rows, function(r) gsub('"', "", strsplit(r, "\t")[[1]][-1]))
}
title <- grab("Sample_title")[[1]]
gsm   <- grab("Sample_geo_accession")[[1]]
ch    <- grab("Sample_characteristics_ch1")
kv <- lapply(ch, function(v) {
  key <- sub(":.*", "", v[1]); list(key = key, val = trimws(sub("^[^:]*:", "", v)))
})
meta <- data.table(gsm = gsm, title = title)
for (x in kv) meta[[gsub(" ", "_", x$key)]] <- x$val

# Map to the file prefix in the GEO tarball
files <- list.files(P$geo_dir, pattern = "_filtered_feature_bc_matrix.h5$")
pref  <- sub("_filtered_feature_bc_matrix.h5$", "", files)
meta[, prefix := pref[match(gsm, sub("_.*", "", pref))]]
meta[, sample_id := sub("^GSM[0-9]+_", "", prefix)]

meta[, idh_status := fifelse(grepl("IDH mutant", tissue), "IDH-mutant", "IDH-wildtype")]
meta[, who_grade := as.integer(tumor_stage)]
meta[, histology := fcase(grepl("oligodendroglioma", title), "Oligodendroglioma",
                          grepl("astrocytoma", title), "Astrocytoma",
                          default = "Glioblastoma")]
meta[, patient := toupper(sub("^(GBM|IDHm) ([A-Z]+[0-9]+).*", "\\2", title))]
meta[, cohort := fifelse(idh_status == "IDH-mutant", "IDHm", "GBM")]
meta[, is_LGG := idh_status == "IDH-mutant" & who_grade %in% 2:3]

# Sample name used by Greenwald et al. (author annotation files)
author_names <- if (file.exists(file.path(P$inputs_dir, "general", "GBM_samples.txt")))
  readLines(file.path(P$inputs_dir, "general", "GBM_samples.txt")) else character()
meta[, author_id := vapply(sample_id, function(s) {
  hit <- author_names[tolower(author_names) == s]
  if (!length(hit)) hit <- author_names[tolower(author_names) == paste0(s, "bulk")]
  if (length(hit)) hit[1] else NA_character_
}, character(1))]

setcolorder(meta, c("sample_id", "gsm", "patient", "cohort", "idh_status", "histology",
                    "who_grade", "is_LGG", "tumor_region", "mgmt_status", "author_id"))
meta[, c("tissue", "tumor_stage") := NULL]
write_tab(meta, "sample_metadata.csv")

cnt <- meta[, .N, by = .(cohort, histology, who_grade)][order(cohort, histology, who_grade)]
write_tab(cnt, "sample_counts.csv")
logf("Total samples: %d | IDH-mutant: %d | LGG (IDHm WHO 2-3): %d | GBM: %d",
     nrow(meta), sum(meta$cohort == "IDHm"), sum(meta$is_LGG), sum(meta$cohort == "GBM"))
print(meta[, .(sample_id, gsm, cohort, histology, who_grade, is_LGG, tumor_region, author_id)])
