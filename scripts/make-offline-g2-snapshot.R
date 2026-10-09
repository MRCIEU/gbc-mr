#!/usr/bin/env Rscript
# Derive a portable association subset from the already supplied local data.
source("scripts/R/offline-core.R")
original <- "data/gbc_phewas_results-g2.Rdata"
snapshot <- "data/offline-g2-inputs.rds"
inputs <- load_saved_inputs(original)
saveRDS(inputs, snapshot, compress = "xz")
stopifnot(identical(inputs, readRDS(snapshot)))
dir.create("results/g2", recursive = TRUE, showWarnings = FALSE)
provenance <- cbind(data.frame(source = original, source_md5 = unname(tools::md5sum(original)),
  snapshot = snapshot, snapshot_md5 = unname(tools::md5sum(snapshot)),
  selection = "Full saved instruments/int_chd/PheWAS; outcome rows matching their SNPs and all eight original leads"),
  inputs$inventory_counts)
write.csv(provenance, "results/g2/input-snapshot-provenance.csv", row.names = FALSE)
cat("Lossless saved-association snapshot verified\n")
