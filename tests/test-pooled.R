# Regression check against the saved multi-target associations. The former
# notebook unioned the complements for each target, which lets one target's
# regional instruments return through another target's complement.
source("scripts/R/offline-core.R")
source("scripts/R/offline-pooled.R")
if (file.exists("data/offline-g2-inputs.rds")) {
  inputs <- readRDS("data/offline-g2-inputs.rds")
} else {
  inputs <- load_saved_inputs()
}
pairs <- canonical_pairs(inputs)
multi_target <- names(which(table(pairs$mediator_id) > 1))
stopifnot(length(multi_target) > 0L)

demonstrated <- FALSE
for (mediator_id in multi_target) {
  targets <- unique(pairs$target_snp[pairs$mediator_id == mediator_id])
  per_target_complements <- lapply(targets, function(target) {
    select_instruments(inputs, mediator_id, target, window_bp = 1e6)$instruments$SNP
  })
  old_complement_union <- unique(unlist(per_target_complements, use.names = FALSE))
  corrected <- select_instruments(inputs, mediator_id, targets, window_bp = 1e6)
  reintroduced <- setdiff(old_complement_union, corrected$instruments$SNP)
  if (length(reintroduced)) {
    stopifnot(all(reintroduced %in% corrected$excluded_instruments$SNP))
    demonstrated <- TRUE
    break
  }
}
stopifnot(demonstrated)

pooled <- compute_pooled(inputs)
stopifnot(nrow(pooled) == length(unique(pairs$mediator_id)),
  all(pooled$family == "pooled_mediator_mr"),
  all(is.na(pooled$pval_bh) | (pooled$pval_bh >= 0 & pooled$pval_bh <= 1)),
  all(pooled$n_saved_before_exclusions >= pooled$n_retained_before_outcome_matching, na.rm = TRUE))
cat("Pooled multi-target union exclusion and result-family checks passed\n")
