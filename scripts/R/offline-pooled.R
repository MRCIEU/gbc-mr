# Pooled mediator-to-GBC MR, rebuilt from saved instruments and outcome data.
# Requires the shared offline helpers in offline-core.R.

compute_pooled <- function(inputs, window_bp = 1e6, mhc = FALSE) {
  pairs <- canonical_pairs(inputs)
  loci <- known_gbc_loci(inputs)
  mediators <- unique(pairs[, c("mediator_id", "mediator", "mediator_scale")])
  rows <- lapply(seq_len(nrow(mediators)), function(i) {
    med <- mediators[i, , drop = FALSE]
    targets <- unique(pairs$target_snp[pairs$mediator_id == med$mediator_id])
    selected <- tryCatch(select_instruments(
      inputs, mediator_id = med$mediator_id, target_snps = targets,
      window_bp = window_bp, mhc = mhc
    ), error = identity)
    selection_error <- inherits(selected, "error")
    selection_error_message <- if (selection_error) conditionMessage(selected) else ""
    if (selection_error) {
      selected <- list(
        instruments = inputs$int_inst[0, , drop = FALSE],
        counts = data.frame(n_saved = NA_real_, n_excluded_target = NA_real_,
          n_excluded_region = NA_real_, n_excluded_mhc = NA_real_,
          n_excluded_coordinate_missing = NA_real_, n_excluded_coordinate_mismatch = NA_real_,
          n_retained = NA_real_),
        excluded_instruments = NULL
      )
    }
    exposure <- selected$instruments
    outcome <- inputs$outgwasf[inputs$outgwasf$SNP %in% exposure$SNP, , drop = FALSE]
    harmonised <- harmonise_local(exposure, outcome)
    hcounts <- attr(harmonised, "counts")
    estimate <- estimate_mr(harmonised)
    if (selection_error) {
      estimate$status <- "unavailable"
      estimate$reason <- selection_error_message
    }

    # Count exclusions that would have had a saved GBC association where the
    # selector exposes row-level exclusions. Older selectors still provide
    # the aggregate counts; their post-match count is explicitly unavailable.
    excluded <- selected$excluded_instruments
    post_match_region <- NA_integer_
    post_match_target <- NA_integer_
    post_match_mhc <- NA_integer_
    post_match_missing <- NA_integer_
    post_match_mismatch <- NA_integer_
    if (!is.null(excluded)) {
      matched <- excluded$SNP %in% inputs$outgwasf$SNP
      role <- if ("exclusion_reason" %in% names(excluded)) excluded$exclusion_reason else rep(NA_character_, nrow(excluded))
      post_match_region <- sum(matched & role == "region", na.rm = TRUE)
      post_match_target <- sum(matched & role == "target", na.rm = TRUE)
      post_match_mhc <- sum(matched & role == "mhc", na.rm = TRUE)
      post_match_missing <- sum(matched & role == "coordinate_missing", na.rm = TRUE)
      post_match_mismatch <- sum(matched & role == "coordinate_mismatch", na.rm = TRUE)
    }

    count <- selected$counts
    count_value <- function(name) if (name %in% names(count)) as.numeric(count[[name]][1]) else NA_real_
    cbind(
      data.frame(
        mediator_id = med$mediator_id,
        mediator = med$mediator,
        mediator_scale = med$mediator_scale,
        n_exclusion_loci = nrow(loci),
        excluded_locus_snps = paste(loci$target_snp, collapse = ";"),
        n_targets = length(targets),
        target_snps = paste(targets, collapse = ";"),
        window_bp = window_bp,
        mhc_exclusion = mhc,
        exclusion_family = "pooled_mediator_mr",
        n_outcome_matched_before_harmonisation = length(intersect(exposure$SNP, inputs$outgwasf$SNP)),
        n_excluded_region_after_outcome_matching = post_match_region,
        n_excluded_target_after_outcome_matching = post_match_target,
        n_excluded_mhc_after_outcome_matching = post_match_mhc,
        n_excluded_coordinate_missing_after_outcome_matching = post_match_missing,
        n_excluded_coordinate_mismatch_after_outcome_matching = post_match_mismatch,
        n_excluded_after_outcome_matching = if (!is.null(excluded)) sum(excluded$SNP %in% inputs$outgwasf$SNP) else NA_integer_,
        n_harmonised = if (!is.null(hcounts) && !is.null(hcounts$n_usable)) hcounts$n_usable else NA_integer_,
        n_removed_during_harmonisation = if (!is.null(hcounts) && !is.null(hcounts$n_usable))
          max(0L, length(intersect(exposure$SNP, inputs$outgwasf$SNP)) - hcounts$n_usable) else NA_integer_,
        n_harmonisation_missing_outcome = if (!is.null(hcounts)) hcounts$n_missing_outcome else NA_integer_,
        n_harmonisation_palindromic = if (!is.null(hcounts)) hcounts$n_palindromic else NA_integer_,
        n_harmonisation_incompatible = if (!is.null(hcounts)) hcounts$n_incompatible else NA_integer_,
        n_saved_before_exclusions = count_value("n_saved"),
        n_excluded_target = count_value("n_excluded_target"),
        n_excluded_region = count_value("n_excluded_region"),
        n_excluded_mhc = count_value("n_excluded_mhc"),
        n_excluded_coordinate_missing = count_value("n_excluded_coordinate_missing"),
        n_excluded_coordinate_mismatch = count_value("n_excluded_coordinate_mismatch"),
        n_retained_before_outcome_matching = count_value("n_retained"),
        stringsAsFactors = FALSE
      ),
      estimate
    )
  })
  result <- do.call(rbind, rows)
  result$family <- "pooled_mediator_mr"
  result$pval_bh <- NA_real_
  available <- !is.na(result$pval) & is.finite(result$pval)
  result$pval_bh[available] <- stats::p.adjust(result$pval[available], method = "BH")
  result
}

# Diagnostic-only comparison with the pooled MR values printed in the original
# stored HTML report. The historical exposure label sometimes contains its
# OpenGWAS ID; a plain display-label fallback is explicitly marked when used.
compare_pooled_legacy_diagnostic <- function(pooled, legacy) {
  stopifnot(all(c("mediator_id", "b", "se", "pval") %in% names(pooled)),
            all(c("exposure", "b", "se", "pval") %in% names(legacy)))
  legacy_id <- if ("mediator_id" %in% names(legacy)) as.character(legacy$mediator_id) else rep(NA_character_, nrow(legacy))
  legacy_id[!nzchar(legacy_id)] <- NA_character_
  mediator_name <- if ("mediator" %in% names(pooled)) as.character(pooled$mediator) else rep(NA_character_, nrow(pooled))
  id_match <- match(as.character(pooled$mediator_id), legacy_id)
  label_match <- match(mediator_name, as.character(legacy$exposure))
  index <- id_match
  method <- rep("unmatched", nrow(pooled))
  method[!is.na(id_match)] <- "mediator_id"
  fallback <- is.na(index) & !is.na(label_match)
  index[fallback] <- label_match[fallback]
  method[fallback] <- "display_label_only_no_legacy_id"
  data.frame(
    mediator_id = pooled$mediator_id,
    mediator = mediator_name,
    legacy_exposure = ifelse(is.na(index), NA_character_, as.character(legacy$exposure[index])),
    matching_method = method,
    corrected_b = pooled$b,
    corrected_se = pooled$se,
    corrected_pval = pooled$pval,
    corrected_nsnp = pooled$nsnp,
    legacy_b = ifelse(is.na(index), NA_real_, as.numeric(legacy$b[index])),
    legacy_se = ifelse(is.na(index), NA_real_, as.numeric(legacy$se[index])),
    legacy_pval = ifelse(is.na(index), NA_real_, as.numeric(legacy$pval[index])),
    legacy_nsnp = ifelse(is.na(index), NA_integer_, as.integer(legacy$nsnp[index])),
    delta_b = ifelse(is.na(index), NA_real_, pooled$b - as.numeric(legacy$b[index])),
    stringsAsFactors = FALSE
  )
}
