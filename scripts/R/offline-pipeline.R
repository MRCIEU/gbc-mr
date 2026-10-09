# Run from the repository root or a Qmd which sets the root explicitly.
for (module in c("core", "pooled", "decomposition", "inference", "plots"))
  source(file.path("scripts", "R", paste0("offline-", module, ".R")))

disable_remote_access <- function() {
  Sys.setenv(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
    ALL_PROXY = "http://127.0.0.1:9")
  block <- quote(stop("Remote access disabled: use the saved local G2 inputs"))
  for (package in c("TwoSampleMR", "ieugwasr")) {
    if (requireNamespace(package, quietly = TRUE)) {
      functions <- if (package == "TwoSampleMR") c("extract_instruments", "extract_outcome_data") else "api_query"
      for (fun in functions) invisible(capture.output(trace(fun, tracer = block,
        where = asNamespace(package), print = FALSE)))
    }
  }
  for (fun in c("download.file", "download.packages", "install.packages"))
    invisible(capture.output(trace(fun, tracer = block, where = asNamespace("utils"), print = FALSE)))
  invisible(TRUE)
}

run_offline_analysis <- function(input_path = "data/gbc_phewas_results-g2.Rdata",
                                 output_dir = "results/g2", figures_dir = "figures") {
  disable_remote_access()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  inputs <- if (grepl("\\.rds$", input_path, ignore.case = TRUE)) readRDS(input_path) else load_saved_inputs(input_path)
  pairs <- canonical_pairs(inputs)
  comparisons <- annotate_comparisons(compute_comparisons(inputs))
  pooled <- compute_pooled(inputs)
  decomposition <- annotate_decomposition(compute_decomposition(comparisons))
  sensitivity <- do.call(rbind, lapply(c(FALSE, TRUE), function(mhc) {
    x <- annotate_comparisons(compute_comparisons(inputs, window_bp = 2e6, mhc = mhc))
    x$analysis <- if (mhc) "+/-2 Mb plus extended MHC chr6 25-34 Mb" else "+/-2 Mb"
    x$q_family <- if (mhc) "sensitivity_pair_heterogeneity_2mb_mhc" else "sensitivity_pair_heterogeneity_2mb"
    x
  }))
  covariance <- decomposition_covariance_sensitivity(comparisons)
  summary <- summarise_inference(comparisons, pooled, decomposition)
  for (name in unique(sensitivity$q_family)) {
    x <- sensitivity[sensitivity$q_family == name, ]
    summary$families <- rbind(summary$families, data.frame(family = name, inventory = nrow(x),
      tested = sum(is.finite(x$q_pval)), unavailable = sum(!is.finite(x$q_pval)),
      nominal_below_005 = sum(x$q_pval < 0.05, na.rm = TRUE), bh_below_005 = sum(x$q_padj < 0.05, na.rm = TRUE)))
  }
  original_snps <- unique(inputs$phewas_result$rsid)
  lead_inventory <- data.frame(target_snp = original_snps,
    retained_pairs = vapply(original_snps, function(snp) sum(pairs$target_snp == snp), integer(1)))
  lead_inventory$followup_status <- ifelse(lead_inventory$retained_pairs > 0,
    "Retained saved follow-up", "No retained saved SNP-mediator association; not analysed")
  oi <- match(lead_inventory$target_snp, inputs$outgwasf$SNP)
  lead_inventory$gbc_pval <- inputs$outgwasf$pval.outcome[oi]
  # Export the exact saved mediator associations alongside the harmonised results.
  source_associations <- inputs$int_chd[!duplicated(inputs$int_chd[c("id.outcome", "SNP")]), ]
  stopifnot(nrow(comparisons) == nrow(pairs), nrow(decomposition) == nrow(pairs),
    nrow(pooled) == length(unique(pairs$mediator_id)),
    all(abs(decomposition$t - decomposition$i - decomposition$d) < 1e-10, na.rm = TRUE))
  representative <- comparisons[c(1, which(comparisons$precision == "Comparator CI includes zero")[1],
    which(comparisons$q_padj < 0.05)[1], 1), ]
  representative$example <- c("Available", "Imprecise comparator", "Incompatible", "Illustrative unavailable comparator")
  representative$source_kind <- c(rep("Saved empirical comparison", 3), "Validation scenario: remove comparator; not an empirical result")
  representative[4, c("comparator_b", "comparator_se", "comparator_lo", "comparator_hi", "q", "q_pval", "q_padj")] <- NA_real_
  representative$interpretation[4] <- "Unavailable comparison"
  representative$precision[4] <- "Comparator unavailable"
  representative$comparator_status[4] <- "unavailable"
  representative$comparator_nsnp[4] <- 0L
  representative$reason[4] <- "Illustrative validation only: comparator association deliberately omitted"
  representative$status[4] <- "unavailable"
  tables <- list(comparisons = comparisons, pooled_mr = pooled, decomposition = decomposition,
    window_sensitivity = sensitivity, covariance_sensitivity = covariance,
    inference_families = summary$families, per_snp_summary = summary$per_snp,
    lead_inventory = lead_inventory, exact_saved_associations = source_associations,
    representative_relationships = representative)
  legacy_path <- file.path(output_dir, "legacy-pooled-mr.csv")
  if (file.exists(legacy_path)) tables$pooled_legacy_diagnostic <- compare_pooled_legacy_diagnostic(pooled, read.csv(legacy_path))
  for (name in names(tables)) {
    path <- file.path(output_dir, paste0(name, ".csv"))
    write.csv(tables[[name]], path, row.names = FALSE, na = "")
    stopifnot(nrow(read.csv(path)) == nrow(tables[[name]]))
  }
  provenance <- data.frame(input = input_path, md5 = unname(tools::md5sum(input_path)),
    retained_pairs = nrow(pairs), mediators = nrow(pooled), targets = length(unique(pairs$target_snp)),
    original_lead_snps = length(original_snps),
    policy = "Exact SNP; +/-1 Mb geographic exclusion; all palindromes removed; no remote access or result cache")
  write.csv(provenance, file.path(output_dir, "provenance.csv"), row.names = FALSE)
  versions <- data.frame(package = c("R", "TwoSampleMR", "ieugwasr", "knitr", "rmarkdown"),
    version = c(as.character(getRversion()), vapply(c("TwoSampleMR", "ieugwasr", "knitr", "rmarkdown"),
      function(p) as.character(packageVersion(p)), character(1))))
  write.csv(versions, file.path(output_dir, "package_versions.csv"), row.names = FALSE)
  writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())), file.path(output_dir, "session-info.txt"))
  plots <- list(comparisons = export_comparison_plot(comparisons,
      file.path(figures_dir, "phewas_followup_results-g2.pdf")),
    decomposition = export_decomposition_plot(decomposition,
      file.path(figures_dir, "phewas_decomposition-g2.pdf")),
    pooled = export_pooled_plot(pooled, file.path(figures_dir, "phewas_pooled_mr-g2.pdf")))
  list(inputs = inputs, pairs = pairs, comparisons = comparisons, pooled = pooled,
    decomposition = decomposition, summary = summary, tables = tables, plots = plots,
    provenance = provenance, versions = versions)
}
