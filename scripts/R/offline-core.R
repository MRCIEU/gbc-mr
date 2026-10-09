# Offline G2 analysis. SNP and GWAS IDs are keys; labels never identify groups.
# Positions are saved coordinates, not a verified LD reference or build conversion.
load_saved_inputs <- function(path = "data/gbc_phewas_results-g2.Rdata") {
  saved <- new.env(parent = emptyenv())
  load(path, envir = saved)
  required <- c("int_chd", "int_inst", "outgwasf", "phewas_result")
  stopifnot(all(required %in% ls(saved)))
  wanted <- unique(c(saved$int_inst$SNP, saved$int_chd$SNP,
                     saved$phewas_result$rsid))
  # Retain exact targets from the full outcome, never from reduced outdat.
  outcome <- as.data.frame(saved$outgwasf[saved$outgwasf$SNP %in% wanted, ])
  list(int_chd = as.data.frame(saved$int_chd),
       int_inst = as.data.frame(saved$int_inst), outgwasf = outcome,
       phewas_result = as.data.frame(saved$phewas_result),
       inventory_counts = data.frame(n_saved_pair_rows = nrow(saved$int_chd),
                                     n_saved_instrument_rows = nrow(saved$int_inst),
                                     n_full_outcome_rows = nrow(saved$outgwasf),
                                     n_retained_outcome_rows = nrow(outcome),
                                     n_original_phewas_rows = nrow(saved$phewas_result)))
}

offline_deduplicate <- function(dat, role) {
  dat <- as.data.frame(dat)
  if (!nrow(dat)) return(list(data = dat, removed = 0L, conflicts = character()))
  key <- paste(dat[[paste0("id.", role)]], dat$SNP, sep = "\r")
  important <- intersect(c(paste0(c("beta", "se", "effect_allele", "other_allele",
                                    "chr", "pos", "mr_keep"), ".", role),
                            "chr", "pos", "proxy.outcome", "target_snp.outcome"),
                         names(dat))
  conflict_keys <- unique(key[duplicated(key)])
  conflict_keys <- conflict_keys[vapply(conflict_keys, function(k) {
    nrow(unique(dat[key == k, important, drop = FALSE])) > 1L
  }, logical(1))]
  keep <- !duplicated(key) & !key %in% conflict_keys
  list(data = dat[keep, , drop = FALSE],
       removed = sum(duplicated(key)),
       conflicts = unique(dat$SNP[key %in% conflict_keys]))
}

offline_coord <- function(dat, role) {
  chr_col <- if (paste0("chr.", role) %in% names(dat)) paste0("chr.", role) else "chr"
  pos_col <- if (paste0("pos.", role) %in% names(dat)) paste0("pos.", role) else "pos"
  chr <- if (chr_col %in% names(dat)) as.character(dat[[chr_col]]) else rep(NA_character_, nrow(dat))
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  pos <- if (pos_col %in% names(dat)) suppressWarnings(as.numeric(dat[[pos_col]])) else rep(NA_real_, nrow(dat))
  data.frame(chr = chr, pos = pos, stringsAsFactors = FALSE)
}

canonical_pairs <- function(inputs) {
  target <- inputs$int_chd
  pairs <- unique(data.frame(mediator_id = as.character(target$id.outcome),
                             target_snp = as.character(target$SNP)))
  ix <- match(paste(pairs$mediator_id, pairs$target_snp),
              paste(target$id.outcome, target$SNP))
  pairs$mediator <- if ("originalname.outcome" %in% names(target)) {
    as.character(target$originalname.outcome[ix])
  } else as.character(target$outcome[ix])
  pairs$mediator_scale <- "Saved GWAS units; specific transformation/unit unavailable"
  out_ix <- match(pairs$target_snp, inputs$outgwasf$SNP)
  coord <- offline_coord(inputs$outgwasf, "outcome")
  fallback <- offline_coord(target, "outcome")[ix, ]
  pairs$target_chr <- coord$chr[out_ix]
  pairs$target_pos <- coord$pos[out_ix]
  absent <- is.na(out_ix)
  pairs$target_chr[absent] <- fallback$chr[absent]
  pairs$target_pos[absent] <- fallback$pos[absent]
  pairs[order(pairs$mediator_id, pairs$target_snp), , drop = FALSE]
}

# The original lead inventory includes hits without retained mediator pairs.
# Geographic exclusions use this complete inventory for every mediator.
known_gbc_loci <- function(inputs) {
  snps <- sort(unique(c(as.character(inputs$phewas_result$rsid),
                        as.character(inputs$int_chd$SNP))))
  snps <- snps[!is.na(snps) & nzchar(snps)]
  loci <- data.frame(target_snp = snps, target_chr = NA_character_, target_pos = NA_real_)
  for (i in seq_along(snps)) {
    saved <- inputs$outgwasf[inputs$outgwasf$SNP == snps[i], , drop = FALSE]
    exact <- inputs$int_chd[inputs$int_chd$SNP == snps[i], , drop = FALSE]
    coords <- rbind(offline_coord(saved, "outcome"), offline_coord(exact, "outcome"))
    # Any usable saved coordinate may establish the locus, but conflicting
    # coordinates cannot certify the geographic exclusion without a build map.
    coords <- unique(coords[!is.na(coords$chr) & nzchar(coords$chr) & is.finite(coords$pos), ])
    if (nrow(coords) == 1L) {
      loci$target_chr[i] <- coords$chr
      loci$target_pos[i] <- coords$pos
    }
  }
  loci
}

select_instruments <- function(inputs, mediator_id, target_snps,
                               window_bp = 1e6, mhc = FALSE) {
  stopifnot(length(mediator_id) == 1L, is.finite(window_bp), window_bp >= 0)
  raw <- inputs$int_inst[inputs$int_inst$id.exposure == mediator_id, , drop = FALSE]
  dd <- offline_deduplicate(raw, "exposure")
  inst <- dd$data
  coords <- offline_coord(inst, "exposure")
  loci <- known_gbc_loci(inputs)
  if (!all(target_snps %in% loci$target_snp) ||
      !nrow(loci) || any(is.na(loci$target_chr) | loci$target_chr == "" | !is.finite(loci$target_pos))) {
    stop("Known GBC hit coordinates unavailable or conflicting: cannot establish geographic exclusions")
  }
  reason <- rep(NA_character_, nrow(inst))
  mark <- function(which, value) {
    reason[which & is.na(reason)] <<- value
  }
  mark(inst$SNP %in% loci$target_snp, "target")
  mark(is.na(coords$chr) | coords$chr == "" | !is.finite(coords$pos), "coordinate_missing")
  out <- offline_deduplicate(inputs$outgwasf, "outcome")$data
  j <- match(inst$SNP, out$SNP)
  oc <- offline_coord(out, "outcome")
  mismatch <- !is.na(j) & !is.na(oc$chr[j]) & is.finite(oc$pos[j]) &
    (coords$chr != oc$chr[j] | coords$pos != oc$pos[j])
  mismatch[is.na(mismatch)] <- FALSE
  mark(mismatch, "coordinate_mismatch")
  for (i in seq_len(nrow(loci))) {
    regional <- coords$chr == loci$target_chr[i] &
      abs(coords$pos - loci$target_pos[i]) <= window_bp
    regional[is.na(regional)] <- FALSE
    mark(regional, "region")
  }
  if (mhc && any(loci$target_chr == "6")) {
    in_mhc <- coords$chr == "6" & coords$pos >= 25e6 & coords$pos <= 34e6
    in_mhc[is.na(in_mhc)] <- FALSE
    mark(in_mhc, "mhc")
  }
  counts <- data.frame(n_saved = nrow(raw), n_duplicates_removed = dd$removed,
                       n_duplicate_conflicts = length(dd$conflicts),
                       n_excluded_target = sum(reason == "target", na.rm = TRUE),
                       n_excluded_region = sum(reason == "region", na.rm = TRUE),
                       n_excluded_mhc = sum(reason == "mhc", na.rm = TRUE),
                       n_excluded_coordinate_missing = sum(reason == "coordinate_missing", na.rm = TRUE),
                       n_excluded_coordinate_mismatch = sum(reason == "coordinate_mismatch", na.rm = TRUE),
                       n_retained = sum(is.na(reason)))
  excluded <- inst[!is.na(reason), , drop = FALSE]
  excluded$exclusion_reason <- reason[!is.na(reason)]
  kept <- inst[is.na(reason), , drop = FALSE]
  counts$retained_SNPs <- paste(kept$SNP, collapse = ";")
  counts$n_exclusion_loci <- nrow(loci)
  counts$excluded_locus_snps <- paste(loci$target_snp, collapse = ";")
  list(instruments = kept, excluded_instruments = excluded, counts = counts,
       duplicate_conflict_snps = dd$conflicts)
}

offline_as_exposure <- function(dat) {
  names(dat) <- sub("\\.outcome$", ".exposure", names(dat))
  names(dat)[names(dat) == "outcome"] <- "exposure"
  if ("chr" %in% names(dat)) names(dat)[names(dat) == "chr"] <- "chr.exposure"
  if ("pos" %in% names(dat)) names(dat)[names(dat) == "pos"] <- "pos.exposure"
  dat
}

harmonise_local <- function(exposure, outcome) {
  ex <- offline_deduplicate(exposure, "exposure")
  ou <- offline_deduplicate(outcome, "outcome")
  e <- ex$data
  o <- ou$data
  match_ix <- match(e$SNP, o$SNP)
  counts <- data.frame(n_missing_outcome = sum(is.na(match_ix)),
                       n_duplicate_conflicts_exposure = length(ex$conflicts),
                       n_duplicate_conflicts_outcome = length(ou$conflicts),
                       n_coordinate_mismatch = 0L, n_palindromic = 0L,
                       n_incompatible = 0L, n_invalid_association = 0L,
                       n_proxy = 0L, n_usable = 0L)
  e <- e[!is.na(match_ix), , drop = FALSE]
  if (!nrow(e)) {
    result <- data.frame()
    attr(result, "counts") <- counts
    attr(result, "reason") <- if (length(ex$conflicts)) "Conflicting duplicate target/exposure associations" else "No matching outcome associations"
    return(result)
  }
  o <- o[o$SNP %in% e$SNP, , drop = FALSE]
  ec <- offline_coord(e, "exposure")
  oc <- offline_coord(o, "outcome")[match(e$SNP, o$SNP), ]
  mismatch <- !is.na(ec$chr) & !is.na(oc$chr) & is.finite(ec$pos) & is.finite(oc$pos) &
    (ec$chr != oc$chr | ec$pos != oc$pos)
  counts$n_coordinate_mismatch <- sum(mismatch)
  proxy <- if ("proxy.exposure" %in% names(e)) e$proxy.exposure %in% TRUE else rep(FALSE, nrow(e))
  if ("target_snp.exposure" %in% names(e)) {
    proxy <- proxy | (!is.na(e$target_snp.exposure) & e$target_snp.exposure != e$SNP)
  }
  counts$n_proxy <- sum(proxy)
  valid <- is.finite(e$beta.exposure) & is.finite(e$se.exposure) & e$se.exposure > 0 &
    e$beta.exposure != 0
  oi <- match(e$SNP, o$SNP)
  valid <- valid & is.finite(o$beta.outcome[oi]) & is.finite(o$se.outcome[oi]) & o$se.outcome[oi] > 0
  if ("mr_keep.exposure" %in% names(e)) valid <- valid & e$mr_keep.exposure %in% TRUE
  if ("mr_keep.outcome" %in% names(o)) valid <- valid & o$mr_keep.outcome[oi] %in% TRUE
  counts$n_invalid_association <- sum(!valid)
  e <- e[!mismatch & !proxy & valid, , drop = FALSE]
  if (nrow(e)) {
    # Action 3 discards every palindromic A/T or C/G variant, even with EAF.
    result <- suppressMessages(TwoSampleMR::harmonise_data(e, o, action = 3))
    counts$n_palindromic <- sum(result$palindromic %in% TRUE)
    counts$n_incompatible <- sum(result$remove %in% TRUE & !result$palindromic %in% TRUE)
    counts$n_usable <- sum(result$mr_keep %in% TRUE)
  } else result <- data.frame()
  attr(result, "counts") <- counts
  attr(result, "reason") <- if (counts$n_usable) "" else paste(
    "No usable harmonised variants;", paste(names(counts)[counts[1, ] > 0], collapse = ", "))
  result
}

estimate_mr <- function(dat) {
  result <- data.frame(b = NA_real_, se = NA_real_, pval = NA_real_, lo = NA_real_,
                       hi = NA_real_, nsnp = 0L, method = NA_character_,
                       status = "unavailable", reason = "No usable instruments",
                       F_min = NA_real_, F_mean = NA_real_, n_weak_F = NA_integer_)
  if (!nrow(dat)) {
    if (!is.null(attr(dat, "reason"))) result$reason <- attr(dat, "reason")
    return(result)
  }
  d <- dat[dat$mr_keep %in% TRUE, , drop = FALSE]
  if (!nrow(d)) {
    result$reason <- attr(dat, "reason")
    return(result)
  }
  stopifnot(!anyDuplicated(d$SNP))
  strength <- (d$beta.exposure / d$se.exposure)^2
  result$F_min <- min(strength)
  result$F_mean <- mean(strength)
  result$n_weak_F <- sum(strength < 10)
  fit <- if (nrow(d) == 1L) {
    result$method <- "Wald ratio"
    TwoSampleMR::mr_wald_ratio(d$beta.exposure, d$beta.outcome,
                               d$se.exposure, d$se.outcome, TwoSampleMR::default_parameters())
  } else {
    result$method <- "IVW (multiplicative random effects)"
    TwoSampleMR::mr_ivw(d$beta.exposure, d$beta.outcome,
                        d$se.exposure, d$se.outcome)
  }
  result$nsnp <- nrow(d)
  if (!is.finite(fit$b) || !is.finite(fit$se) || fit$se <= 0) {
    result$reason <- "MR estimate or standard error is non-finite"
    return(result)
  }
  result$b <- fit$b
  result$se <- fit$se
  result$pval <- fit$pval
  result$lo <- fit$b - qnorm(0.975) * fit$se
  result$hi <- fit$b + qnorm(0.975) * fit$se
  result$status <- "available"
  result$reason <- ""
  result
}

compute_comparisons <- function(inputs, window_bp = 1e6, mhc = FALSE) {
  pairs <- canonical_pairs(inputs)
  mediator_ids <- unique(pairs$mediator_id)
  loci <- known_gbc_loci(inputs)
  comparators <- setNames(lapply(mediator_ids, function(mediator_id) {
    targets <- unique(pairs$target_snp[pairs$mediator_id == mediator_id])
    selected <- tryCatch(select_instruments(inputs, mediator_id, targets,
                                            window_bp, mhc), error = identity)
    if (inherits(selected, "error")) {
      comparator <- estimate_mr(data.frame())
      comparator$reason <- conditionMessage(selected)
      counts <- data.frame(n_saved = sum(inputs$int_inst$id.exposure == mediator_id),
                           n_duplicates_removed = NA_integer_, n_duplicate_conflicts = NA_integer_,
                           n_excluded_target = NA_integer_, n_excluded_region = NA_integer_,
                           n_excluded_mhc = NA_integer_, n_excluded_coordinate_missing = NA_integer_,
                           n_excluded_coordinate_mismatch = NA_integer_, n_retained = NA_integer_,
                           retained_SNPs = "", n_exclusion_loci = nrow(loci),
                           excluded_locus_snps = paste(loci$target_snp, collapse = ";"))
      hc <- attr(harmonise_local(inputs$int_inst[FALSE, ], inputs$outgwasf), "counts")
    } else {
      comp_h <- harmonise_local(selected$instruments, inputs$outgwasf)
      comparator <- estimate_mr(comp_h)
      counts <- selected$counts
      hc <- attr(comp_h, "counts")
    }
    list(estimate = comparator, counts = counts, harmonisation_counts = hc)
  }), mediator_ids)
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    pair <- pairs[i, , drop = FALSE]
    exact <- inputs$int_chd[inputs$int_chd$id.outcome == pair$mediator_id &
                            inputs$int_chd$SNP == pair$target_snp, , drop = FALSE]
    h <- harmonise_local(offline_as_exposure(exact), inputs$outgwasf)
    target <- estimate_mr(h)
    # A total G->Y association remains interpretable even when G->M is ambiguous.
    raw_out <- inputs$outgwasf[inputs$outgwasf$SNP == pair$target_snp, , drop = FALSE]
    raw_out <- offline_deduplicate(raw_out, "outcome")$data
    raw_valid <- nrow(raw_out) == 1L && is.finite(raw_out$beta.outcome) &&
      is.finite(raw_out$se.outcome) && raw_out$se.outcome > 0
    pair$effect_allele <- if (raw_valid) raw_out$effect_allele.outcome else NA_character_
    pair$a <- pair$a_se <- NA_real_
    pair$t <- if (raw_valid) raw_out$beta.outcome else NA_real_
    pair$t_se <- if (raw_valid) raw_out$se.outcome else NA_real_
    if (target$nsnp == 1L && target$status == "available") {
      stopifnot(nrow(h[h$mr_keep %in% TRUE, ]) == 1L,
                h$SNP[h$mr_keep %in% TRUE] == pair$target_snp)
      hh <- h[h$mr_keep %in% TRUE, ]
      pair$effect_allele <- hh$effect_allele.exposure
      pair$a <- hh$beta.exposure
      pair$a_se <- hh$se.exposure
      pair$t <- hh$beta.outcome
      pair$t_se <- hh$se.outcome
    }
    pair$target_status <- target$status
    pair$target_reason <- target$reason
    for (field in c("b", "se", "pval", "lo", "hi", "nsnp", "method", "F_min", "F_mean", "n_weak_F")) {
      pair[[paste0("target_", field)]] <- target[[field]]
    }
    shared <- comparators[[pair$mediator_id]]
    comparator <- shared$estimate
    counts <- shared$counts
    hc <- shared$harmonisation_counts
    for (field in c("b", "se", "pval", "lo", "hi", "nsnp", "method", "F_min", "F_mean", "n_weak_F", "status", "reason")) {
      pair[[paste0("comparator_", field)]] <- comparator[[field]]
    }
    pair$q <- pair$q_pval <- NA_real_
    pair$status <- if (target$status == "available" && comparator$status == "available") "available" else "unavailable"
    pair$reason <- paste(c(if (target$status != "available") paste("Target:", target$reason),
                          if (comparator$status != "available") paste("Comparator:", comparator$reason)), collapse = "; ")
    if (pair$status == "available") {
      pair$q <- (target$b - comparator$b)^2 / (target$se^2 + comparator$se^2)
      pair$q_pval <- pchisq(pair$q, df = 1, lower.tail = FALSE)
    }
    pair$window_bp <- window_bp
    pair$mhc_exclusion <- mhc
    cbind(pair, counts, hc)
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

comparison_sensitivity <- function(inputs) {
  wider <- compute_comparisons(inputs, window_bp = 2e6, mhc = FALSE)
  wider$sensitivity <- "+/-2 Mb"
  mhc <- compute_comparisons(inputs, window_bp = 2e6, mhc = TRUE)
  mhc$sensitivity <- "+/-2 Mb plus extended MHC (chr6:25-34 Mb)"
  rbind(wider, mhc)
}
