source("scripts/R/offline-pipeline.R")
snapshot <- readRDS("data/offline-g2-inputs.rds")
pairs <- canonical_pairs(snapshot)
stopifnot(nrow(pairs) == 33L, length(unique(pairs$mediator_id)) == 26L,
  length(unique(pairs$target_snp)) == 6L,
  length(unique(snapshot$phewas_result$rsid)) == 8L,
  all(unique(pairs$target_snp) %in% snapshot$outgwasf$SNP))
out <- file.path(tempdir(), "gbc-offline-validation")
analysis <- run_offline_analysis("data/offline-g2-inputs.rds", out, file.path(out, "figures"))
stopifnot(nrow(analysis$comparisons) == 33L, nrow(analysis$pooled) == 26L,
  nrow(analysis$decomposition) == 33L,
  all(analysis$comparisons$target_nsnp == 1L),
  nrow(analysis$tables$window_sensitivity) == 66L,
  sum(analysis$tables$lead_inventory$retained_pairs == 0L) == 2L,
  all(abs(analysis$decomposition$t - analysis$decomposition$i - analysis$decomposition$d) < 1e-12),
  nrow(read.csv(file.path(out, "comparisons.csv"))) == 33L)
# API calls fail at entry, before credentials or network access can be used.
stopifnot(inherits(try(TwoSampleMR::extract_instruments("fake-id"), silent = TRUE), "try-error"),
  inherits(try(ieugwasr::api_query("gwasinfo"), silent = TRUE), "try-error"))
cat("Portable offline pipeline, canonical inventory and network guards passed\n")
