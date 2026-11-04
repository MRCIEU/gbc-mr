#!/usr/bin/env Rscript
# Run from the repository root. No remote extraction or cached MR results.
source("scripts/R/offline-pipeline.R")
args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1] else "data/gbc_phewas_results-g2.Rdata"
analysis <- run_offline_analysis(path)
print(analysis$summary$families)
