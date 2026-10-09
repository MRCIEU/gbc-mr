source("scripts/R/offline-inference.R")
x <- data.frame(mediator_id = c("m1", "m2", "m3"), target_snp = c("s1", "s1", "s2"),
  q_pval = c(0.001, 0.04, NA), comparator_lo = c(0.1, -2, NA), comparator_hi = c(1, 2, NA))
y <- annotate_comparisons(x)
stopifnot(identical(y$q_padj, c(0.002, 0.04, NA_real_)),
  y$interpretation[3] == "Unavailable comparison", y$precision[2] == "Comparator CI includes zero")
s <- summarise_inference(y, data.frame(pval = c(0.01, NA)),
  data.frame(i_pval = c(0.01, 0.1, NA), d_pval = c(0.04, 0.5, NA)))
stopifnot(s$families$tested[1] == 2, s$families$unavailable[1] == 1,
  s$per_snp$pairs[s$per_snp$target_snp == "s1"] == 2)
duplicate <- rbind(x, x[1, ])
stopifnot(inherits(try(annotate_comparisons(duplicate), silent = TRUE), "try-error"))
cat("Inference checks passed\n")
