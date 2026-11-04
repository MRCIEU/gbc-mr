# Every mediator uses one comparator free of all known GBC lead loci, including
# original leads without a retained association with that mediator.
source("scripts/R/offline-core.R")
source("scripts/R/offline-pooled.R")
inputs <- readRDS("data/offline-g2-inputs.rds")
comparisons <- compute_comparisons(inputs)
known_hits <- unique(c(inputs$phewas_result$rsid, inputs$int_chd$SNP))
for (id in unique(comparisons$mediator_id)) {
  x <- comparisons[comparisons$mediator_id == id, ]
  stopifnot(length(unique(x$comparator_nsnp)) == 1L,
    length(unique(x$comparator_b)) == 1L,
    length(unique(x$comparator_se)) == 1L,
    length(unique(x$retained_SNPs)) == 1L,
    !any(known_hits %in% strsplit(x$retained_SNPs[1], ";", fixed = TRUE)[[1]]))
}
pooled <- compute_pooled(inputs)
i <- match(comparisons$mediator_id, pooled$mediator_id)
stopifnot(identical(comparisons$comparator_nsnp, pooled$nsnp[i]),
  isTRUE(all.equal(comparisons$comparator_b, pooled$b[i])),
  isTRUE(all.equal(comparisons$comparator_se, pooled$se[i])))
cat("One shared comparator per mediator, excluding all known GBC hits, passed.\n")

# A known lead without any retained mediator pair and a SNP near that lead must
# both be excluded; absence from the plotted targets is no exemption.
id <- "synthetic-mediator"
ex <- inputs$int_inst[1:4, ]
ex$id.exposure <- id
ex$SNP <- c("pair-hit", "unpaired-hit", "near-unpaired-hit", "background")
ex$chr.exposure <- c(1, 2, 2, 3)
ex$pos.exposure <- c(10e6, 20e6, 20.5e6, 30e6)
ou <- inputs$outgwasf[1:4, ]
ou$SNP <- ex$SNP
ou$chr.outcome <- ex$chr.exposure
ou$pos.outcome <- ex$pos.exposure
pair <- offline_as_exposure(inputs$int_chd[1, ])
pair$id.exposure <- id
pair$SNP <- "pair-hit"
pair$chr.exposure <- 1
pair$pos.exposure <- 10e6
names(pair) <- sub("\\.exposure$", ".outcome", names(pair))
fixture <- list(int_inst = ex, int_chd = pair, outgwasf = ou,
  phewas_result = data.frame(rsid = c("pair-hit", "unpaired-hit")))
selected <- select_instruments(fixture, id, "pair-hit")
stopifnot(identical(selected$instruments$SNP, "background"),
  selected$counts$n_excluded_target == 2L,
  selected$counts$n_excluded_region == 1L,
  selected$counts$n_exclusion_loci == 2L)
# Fail closed for an unpaired lead whose exclusion coordinates are unavailable.
fixture$outgwasf$pos.outcome[2] <- NA_real_
failed <- compute_comparisons(fixture)
stopifnot(failed$comparator_status == "unavailable",
  grepl("Known GBC hit coordinates", failed$comparator_reason))
cat("Unpaired known leads, surrounding regions and missing coordinates passed.\n")
# Coordinate disagreements between saved GBC and exact mediator records also
# prevent certification that every known-hit region was excluded.
fixture$outgwasf$pos.outcome[2] <- 20e6
fixture$int_chd$pos.outcome <- 11e6
failed <- compute_comparisons(fixture)
stopifnot(failed$comparator_status == "unavailable",
  grepl("conflicting", failed$comparator_reason))
cat("Conflicting lead coordinates fail closed.\n")
