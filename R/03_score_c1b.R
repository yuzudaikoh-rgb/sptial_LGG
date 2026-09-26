# 03 - Consensus1b scoring, sensitivity variants, matched random controls,
#      Greenwald metaprogram (MP) scores and spot MP annotation
source("R/utils/common.R")
source("R/utils/scoring.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })

meta <- available_samples()
scfg <- CFG$scoring
MPS  <- readRDS(file.path(P$inputs_dir, "MP", "clean_spatial_gbm_metaprograms_124.rds"))
author_dir <- file.path(P$inputs_dir, "MP", "mp_assign_124")

cov_rows <- list(); loo_rows <- list(); agree_rows <- list(); spots <- list()

for (i in seq_len(nrow(meta))) {
  s <- meta$sample_id[i]
  so <- readRDS(obj_path(s))
  counts <- GetAssayData(so, assay = "Spatial", layer = "counts")
  lognorm <- GetAssayData(so, assay = "Spatial", layer = "data")
  seed_s <- CFG$seed + i

  # ---- gene coverage --------------------------------------------------------
  mapped <- map_genes(C1B_GENES, rownames(counts))
  ok_g <- !is.na(mapped)
  cov_rows[[s]] <- data.table(sample_id = s, gene = C1B_GENES, feature = mapped, detected = ok_g,
                              mean_lognorm = NA_real_, frac_spots = NA_real_)
  cov_rows[[s]][ok_g, `:=`(mean_lognorm = Matrix::rowMeans(lognorm[mapped[ok_g], , drop = FALSE]),
                           frac_spots = Matrix::rowMeans(counts[mapped[ok_g], , drop = FALSE] > 0))]
  c1b <- unname(mapped[!is.na(mapped)])

  vars <- lapply(c1b_variants(), function(g) { m <- map_genes(g, rownames(counts)); unname(m[!is.na(m)]) })
  # C1b without any gene shared with a Greenwald MP (guards against circularity)
  vars$C1b_noMPoverlap <- setdiff(vars$C1b, unlist(MPS))

  # ---- UCell: C1b + variants + Greenwald MPs --------------------------------
  sc_c1b <- ucell_scores(counts, vars, maxRank = scfg$ucell_maxrank, ncores = CFG$n_cores)
  sc_mp  <- ucell_scores(counts, MPS, maxRank = scfg$ucell_maxrank, ncores = CFG$n_cores)
  colnames(sc_mp) <- paste0("MP_", colnames(sc_mp))

  # ---- AddModuleScore (Seurat) as secondary scoring method -------------------
  so <- AddModuleScore(so, features = list(c1b), name = "C1b_AMS", seed = seed_s,
                       assay = "Spatial", nbin = 24, ctrl = 100)

  # ---- Author-style Tirosh MP scores and spot annotation ---------------------
  lg <- author_logcpm(counts)
  tsc <- tirosh_scores(lg, MPS, seed = seed_s)
  mp_label <- colnames(tsc)[max.col(replace(tsc, is.na(tsc), -Inf), ties.method = "first")]
  names(mp_label) <- rownames(tsc)
  tsc_c1b <- tirosh_scores(lg, list(C1b = setdiff(c1b, grep("^RP[SL]", c1b, value = TRUE))), seed = seed_s)

  author_label <- rep(NA_character_, ncol(so))
  if (!is.na(meta$author_id[i]) && file.exists(file.path(author_dir, paste0(meta$author_id[i], ".rds")))) {
    a <- readRDS(file.path(author_dir, paste0(meta$author_id[i], ".rds")))
    author_label <- as.character(a$spot_type_meta_new[match(colnames(so), a$SpotID)])
    ok <- !is.na(author_label)
    agree_rows[[s]] <- data.table(sample_id = s, n_shared = sum(ok),
                                  agreement = mean(author_label[ok] == mp_label[colnames(so)][ok]))
  }

  # ---- leave-one-gene-out stability -----------------------------------------
  loo_sets <- setNames(lapply(c1b, function(g) setdiff(c1b, g)), paste0("LOO_", c1b))
  loo <- ucell_scores(counts, loo_sets, maxRank = scfg$ucell_maxrank, ncores = CFG$n_cores)
  loo_rows[[s]] <- data.table(sample_id = s, gene_removed = c1b,
                              spearman_vs_full = apply(loo, 2, cor, y = sc_c1b[, "C1b"], method = "spearman"))

  # ---- matched random gene sets ---------------------------------------------
  rs <- matched_random_sets(lognorm, c1b, n_sets = scfg$n_random_sets,
                            n_bins = scfg$n_expr_bins, seed = seed_s)
  rnd <- ucell_scores(counts, rs, maxRank = scfg$ucell_maxrank, ncores = CFG$n_cores)
  saveRDS(list(sets = rs, scores = rnd), obj_path(s, "random_sets"))

  df <- data.table(sample_id = s, barcode = colnames(so),
                   array_row = so$array_row, array_col = so$array_col,
                   nCount = so$nCount_Spatial, cluster = as.character(so$seurat_clusters))
  df <- cbind(df, as.data.table(sc_c1b), C1b_AMS = so$C1b_AMS1,
              C1b_Tirosh = tsc_c1b[colnames(so), 1],
              as.data.table(sc_mp),
              as.data.table(setNames(as.data.frame(tsc[colnames(so), ]), paste0("T_", colnames(tsc)))),
              mp_label = mp_label[colnames(so)], author_mp = author_label)
  spots[[s]] <- df
  logf("%s: C1b genes %d/90, random sets done, UCell-AMS rho=%.2f", s, length(c1b),
       cor(df$C1b, df$C1b_AMS, method = "spearman"))
  rm(so, counts, lognorm, lg, rnd); gc(verbose = FALSE)
}

spots <- rbindlist(spots, fill = TRUE)
saveRDS(spots, file.path(P$objects, "spot_scores.rds"))
cov <- rbindlist(cov_rows); write_tab(cov, "c1b_gene_coverage_per_sample.csv")
covs <- cov[, .(n_detected = sum(detected), missing = paste(gene[!detected], collapse = ";")), by = sample_id]
write_tab(covs, "c1b_gene_coverage_summary.csv")
loo <- rbindlist(loo_rows); write_tab(loo, "c1b_leave_one_out.csv")
if (length(agree_rows)) write_tab(rbindlist(agree_rows), "mp_annotation_agreement_vs_authors.csv")

# method concordance per sample
conc <- spots[, .(rho_UCell_AMS = cor(C1b, C1b_AMS, method = "spearman"),
                  rho_UCell_Tirosh = cor(C1b, C1b_Tirosh, method = "spearman"),
                  rho_full_noRibo = cor(C1b, C1b_noRibo, method = "spearman"),
                  rho_MHCII_nonMHCII = cor(C1b_MHCII, C1b_nonMHCII, method = "spearman")),
              by = sample_id]
write_tab(conc, "c1b_scoring_method_concordance.csv")
