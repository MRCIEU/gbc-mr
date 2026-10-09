# Run from repository root: Rscript tests/test-exact-hits.R
source("scripts/R/offline-core.R")

exposure <- data.frame(
  id.exposure = "mediator-1", exposure = "Example mediator",
  SNP = c("target", "boundary", "outside", "other", "mismatch", "unknown"),
  chr.exposure = c(1, 1, 1, 2, 2, NA),
  pos.exposure = c(10e6, 11e6, 12e6, 4e6, 5e6, NA),
  beta.exposure = c(.2, .1, .3, .4, .5, .6), se.exposure = .02,
  effect_allele.exposure = "A", other_allele.exposure = "G",
  eaf.exposure = .2, pval.exposure = 1e-9, mr_keep.exposure = TRUE)
outcome <- data.frame(
  id.outcome = "gbc", outcome = "GBC", SNP = exposure$SNP,
  chr.outcome = c(1, 1, 1, 2, 2, 2),
  pos.outcome = c(10e6, 11e6, 12e6, 4e6, 6e6, 7e6),
  beta.outcome = c(-.1, .05, .15, .19, .25, .3), se.outcome = .04,
  effect_allele.outcome = c("G", rep("A", 5)),
  other_allele.outcome = c("A", rep("G", 5)),
  eaf.outcome = c(.8, rep(.2, 5)), pval.outcome = .01, mr_keep.outcome = TRUE)
target <- exposure[1, ]
names(target) <- sub("\\.exposure$", ".outcome", names(target))
target$originalname.outcome <- "Example mediator"
inputs <- list(int_chd = rbind(target, target), int_inst = rbind(exposure, exposure[3, ]),
               outgwasf = outcome, phewas_result = data.frame(rsid = "target"))
stopifnot(nrow(canonical_pairs(inputs)) == 1L)

# Geographic boundary is inclusive; exact target is never a comparator.
selected <- select_instruments(inputs, "mediator-1", "target")
stopifnot(setequal(selected$instruments$SNP, c("outside", "other")),
          selected$counts$n_duplicates_removed == 1L,
          selected$counts$n_excluded_target == 1L,
          selected$counts$n_excluded_region == 1L,
          selected$counts$n_excluded_coordinate_mismatch == 1L,
          selected$counts$n_excluded_coordinate_missing == 1L)
wider <- select_instruments(inputs, "mediator-1", "target", window_bp = 2e6)
stopifnot(identical(as.character(wider$instruments$SNP), "other"))

# An allele reversal must reverse the GBC beta and retain the exact SNP only.
h <- harmonise_local(offline_as_exposure(target), outcome)
stopifnot(nrow(h) == 1L, h$SNP == "target", h$mr_keep,
          abs(h$beta.outcome - .1) < 1e-12,
          h$effect_allele.exposure == "A")
r <- compute_comparisons(inputs)
stopifnot(r$target_nsnp == 1L, abs(r$target_b - .5) < 1e-12,
          abs(r$t - .1) < 1e-12, r$comparator_nsnp == 2L,
          abs(r$comparator_b - (.3 * .15 + .4 * .19) / (.3^2 + .4^2)) < 1e-12,
          is.finite(r$q), r$status == "available")

# Contradictory duplicates are removed, rather than choosing a convenient row.
conflict <- inputs
conflict$int_inst$beta.exposure[nrow(conflict$int_inst)] <- .8
dd <- select_instruments(conflict, "mediator-1", "target")
stopifnot(!"outside" %in% dd$instruments$SNP,
          identical(dd$duplicate_conflict_snps, "outside"))
conflict$int_chd$beta.outcome[2] <- .4
bad <- compute_comparisons(conflict)
stopifnot(bad$target_nsnp == 0L, is.na(bad$a), is.finite(bad$t),
          grepl("Conflicting duplicate", bad$target_reason))

# Frequencies do not rescue palindrome ambiguity under action=3.
palindrome <- inputs
palindrome$int_chd$other_allele.outcome <- "T"
palindrome$outgwasf$effect_allele.outcome[1] <- "A"
palindrome$outgwasf$other_allele.outcome[1] <- "T"
pal <- compute_comparisons(palindrome)
stopifnot(is.na(pal$a), is.na(pal$target_b), is.finite(pal$t),
          pal$target_status == "unavailable", pal$effect_allele == "A")

# Coordinate disagreements and proxy records cannot count as exact estimates.
mismatch <- inputs
mismatch$int_chd$pos.outcome <- 10e6 + 1
bad <- compute_comparisons(mismatch)
stopifnot(is.na(bad$a), is.finite(bad$t), grepl("coordinate_mismatch", bad$target_reason))
proxy <- inputs
proxy$int_chd$proxy.outcome <- TRUE
bad <- compute_comparisons(proxy)
stopifnot(is.na(bad$a), grepl("n_proxy", bad$target_reason))

# Extended MHC removes far-away chr6 instruments as a separate sensitivity.
mhc <- inputs
mhc$int_chd$chr.outcome <- 6
mhc$int_chd$pos.outcome <- 31e6
mhc$int_inst$chr.exposure <- c(6, 6, 6, 2, 2, NA, 6)
mhc$int_inst$pos.exposure <- c(31e6, 32e6, 27e6, 4e6, 5e6, NA, 27e6)
mhc$outgwasf$chr.outcome[1:3] <- 6
mhc$outgwasf$pos.outcome[1:3] <- c(31e6, 32e6, 27e6)
primary <- select_instruments(mhc, "mediator-1", "target")
extended <- select_instruments(mhc, "mediator-1", "target", mhc = TRUE)
stopifnot("outside" %in% primary$instruments$SNP,
          !"outside" %in% extended$instruments$SNP,
          extended$counts$n_excluded_mhc == 1L)

cat("Exact-hit membership, harmonisation, duplicate, coordinate and exclusion checks passed.\n")
