#!/usr/bin/env Rscript
files <- list.files("tests", pattern = "^test-.*\\.R$", full.names = TRUE)
for (file in files) {
  cat("Running", file, "\n")
  source(file, local = new.env(parent = globalenv()))
}
cat("All", length(files), "test files passed\n")
