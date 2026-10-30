#!/usr/bin/env Rscript
# How much of each GBC lead's total effect runs through its candidate mediators.
#
# Mediators: pooled mediator MR P < 0.05 (results/g2/pooled_mr.csv), with
# "Disorders of gallbladder, biliary tract and pancreas" folded into
# cholelithiasis as the same liability. A lead with one such mediator keeps its
# univariable decomposition (results/g2/decomposition.csv). A lead with several
# gets multivariable MR (MVMR) of those mediators on GBC, so that
#   t = sum_j a_j * b_j + d
# where b_j is mediator j's effect conditional on the others. The sum is the
# effect through the set; it does not double-count correlated mediators.
#
# MVMR instruments come from OpenGWAS (ieugwasr token from ~/.Renviron) and are
# saved to data/mvmr-inputs.rds on first run; later runs use the saved copy.
#
# Usage: Rscript scripts/mvmr-mediation.R [gbc_phewas_results-g2.Rdata]

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ggplot2)
})

rdata <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(rdata)) rdata <- "data/gbc_phewas_results-g2.Rdata"
cache_file <- "data/mvmr-inputs.rds"
results_dir <- "results/g2"
figures_dir <- "figures"

cholelithiasis <- "ebi-a-GCST90018819"
duplicates <- c("finn-b-K11_GALLBILPANC" = cholelithiasis)

pooled <- fread(file.path(results_dir, "pooled_mr.csv"))
decomposition <- fread(file.path(results_dir, "decomposition.csv"))
inventory <- fread(file.path(results_dir, "lead_inventory.csv"))

selected <- setdiff(pooled[pval < 0.05, mediator_id], names(duplicates))
pairs <- decomposition[mediator_id %in% selected]
sets <- split(pairs$mediator_id, pairs$target_snp)
mv_sets <- Filter(function(ids) length(ids) > 1, sets)
set_key <- function(ids) paste(sort(ids), collapse = "+")

# ---- OpenGWAS extraction (first run only) ---------------------------------

if (!file.exists(cache_file)) {
  message("Extracting MVMR instruments from OpenGWAS for ", length(mv_sets), " mediator set(s)")
  exposures <- lapply(mv_sets, function(ids)
    mv_extract_exposures(sort(ids), find_proxies = FALSE))
  names(exposures) <- vapply(mv_sets, set_key, character(1))
  saveRDS(list(exposures = exposures, extracted = format(Sys.time(), tz = "UTC", usetz = TRUE),
               versions = c(TwoSampleMR = as.character(packageVersion("TwoSampleMR")),
                            ieugwasr = as.character(packageVersion("ieugwasr")))),
          cache_file)
}
cache <- readRDS(cache_file)
missing_sets <- setdiff(vapply(mv_sets, set_key, character(1)), names(cache$exposures))
if (length(missing_sets)) stop("Delete ", cache_file, " to re-extract; missing: ", paste(missing_sets, collapse = ", "))

# ---- GBC outcome and lead regions -----------------------------------------

saved <- new.env()
load(rdata, envir = saved)
outcome <- as.data.table(saved$outgwasf)
outcome <- outcome[!SNP %in% SNP[duplicated(SNP)]]
leads <- outcome[SNP %in% inventory$target_snp, .(SNP, chr = chr.outcome, pos = pos.outcome)]
stopifnot(nrow(leads) == nrow(inventory))

in_lead_region <- function(chr, pos, window_bp, mhc) {
  hit <- is.na(chr) | is.na(pos)
  for (k in seq_len(nrow(leads)))
    hit <- hit | (chr == leads$chr[k] & abs(pos - leads$pos[k]) <= window_bp)
  if (mhc) hit <- hit | (chr == "6" & pos >= 25e6 & pos <= 34e6)
  hit
}

palindromic <- function(a1, a2) paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")

# ---- MVMR -----------------------------------------------------------------

# IVW with multiplicative random effects; returns the full covariance matrix so
# the indirect-effect variance can include Cov(b_j, b_k).
mvmr_ivw <- function(X, y, se_y) {
  w <- 1 / se_y^2
  fit <- lm(y ~ 0 + X, weights = w)
  phi <- sum(w * resid(fit)^2) / (nrow(X) - ncol(X))
  V <- solve(crossprod(X * sqrt(w))) * max(1, phi)
  b <- unname(coef(fit))
  list(b = b, V = V, se = sqrt(diag(V)), Q = sum(w * resid(fit)^2), Q_df = nrow(X) - ncol(X))
}

# Sanderson-Windmeijer conditional F (as in MVMR::strength_mvmr), with delta
# from OLS of exposure j on the others and no covariance between the exposures'
# SNP associations (sample overlap unknown).
conditional_f <- function(X, S) {
  vapply(seq_len(ncol(X)), function(j) {
    others <- X[, -j, drop = FALSE]
    delta <- coef(lm(X[, j] ~ 0 + others))
    sigma2 <- S[, j]^2 + drop(S[, -j, drop = FALSE]^2 %*% delta^2)
    sum((X[, j] - others %*% delta)^2 / sigma2) / (nrow(X) - (ncol(X) - 1))
  }, numeric(1))
}

run_mvmr <- function(exposure_dat, window_bp, mhc, label) {
  exposure_dat <- exposure_dat[!palindromic(exposure_dat$effect_allele.exposure, exposure_dat$other_allele.exposure), ]
  out <- as.data.frame(outcome[SNP %in% exposure_dat$SNP])
  out <- out[!in_lead_region(out$chr.outcome, out$pos.outcome, window_bp, mhc), ]
  out <- out[!palindromic(out$effect_allele.outcome, out$other_allele.outcome), ]
  h <- mv_harmonise_data(exposure_dat, out, harmonise_strictness = 2)
  ids <- colnames(h$exposure_beta)
  fit <- mvmr_ivw(h$exposure_beta, h$outcome_beta, h$outcome_se)
  # Univariable IVW on each mediator's own P < 5e-8 SNPs within this SNP set,
  # to separate the effect of conditioning from the effect of losing instruments.
  own <- lapply(ids, function(j) {
    k <- h$exposure_pval[, j] < 5e-8
    if (sum(k) < 2) return(c(NA_real_, NA_real_, sum(k)))
    u <- mr_ivw(h$exposure_beta[k, j], h$outcome_beta[k], h$exposure_se[k, j], h$outcome_se[k])
    c(u$b, u$se, sum(k))
  })
  fit$table <- data.table(
    analysis = label, mediator_id = ids,
    mediator = pooled$mediator[match(ids, pooled$mediator_id)],
    b = fit$b, se = fit$se, pval = 2 * pnorm(-abs(fit$b / fit$se)),
    nsnp = nrow(h$exposure_beta), conditional_F = conditional_f(h$exposure_beta, h$exposure_se),
    Q = fit$Q, Q_df = fit$Q_df, Q_pval = pchisq(fit$Q, fit$Q_df, lower.tail = FALSE),
    univariable_b = pooled$b[match(ids, pooled$mediator_id)],
    univariable_nsnp = pooled$nsnp[match(ids, pooled$mediator_id)],
    same_snps_univariable_b = vapply(own, `[`, numeric(1), 1),
    same_snps_univariable_se = vapply(own, `[`, numeric(1), 2),
    same_snps_own_instruments = vapply(own, `[`, numeric(1), 3))
  colnames(fit$V) <- rownames(fit$V) <- ids
  fit
}

mvmr <- list()
for (snp in names(mv_sets)) {
  dat <- cache$exposures[[set_key(mv_sets[[snp]])]]
  mvmr[[snp]] <- list(
    primary = run_mvmr(dat, 1e6, FALSE, "Primary: +/-1 Mb around all GBC leads"),
    sensitivity = run_mvmr(dat, 2e6, TRUE, "Sensitivity: +/-2 Mb plus extended MHC"))
}
mvmr_table <- rbindlist(lapply(names(mvmr), function(snp)
  rbind(cbind(target_snp = snp, mvmr[[snp]]$primary$table),
        cbind(target_snp = snp, mvmr[[snp]]$sensitivity$table))))

# ---- Decomposition per lead -----------------------------------------------

z <- qnorm(0.975)
decompose <- function(snp, fit_name = "primary") {
  rows <- pairs[target_snp == snp]
  if (!nrow(rows)) {
    o <- outcome[SNP == snp]
    return(data.table(target_snp = snp, effect_allele = o$effect_allele.outcome, method = "No mediator",
                      t = o$beta.outcome, t_se = o$se.outcome, component = "Residual direct",
                      estimate = o$beta.outcome, se = o$se.outcome))
  }
  stopifnot(uniqueN(rows$effect_allele) == 1, uniqueN(rows$t) == 1)
  if (nrow(rows) == 1) {
    b <- rows$comparator_b; V <- matrix(rows$comparator_se^2)
    method <- "Univariable MR"
  } else {
    fit <- mvmr[[snp]][[fit_name]]
    b <- fit$b[match(rows$mediator_id, colnames(fit$V))]
    V <- fit$V[rows$mediator_id, rows$mediator_id]
    method <- "MVMR"
  }
  a <- rows$a
  i <- a * b
  # Delta method: Cov(a_j, a_k) and Cov(a, b) taken as zero; Cov(b) from MVMR.
  i_var <- sum(b^2 * rows$a_se^2) + drop(t(a) %*% V %*% a)
  i_se <- sqrt(b^2 * rows$a_se^2 + a^2 * diag(V))
  d <- rows$t[1] - sum(i)
  d_se <- sqrt(rows$t_se[1]^2 + i_var)
  data.table(target_snp = snp, effect_allele = rows$effect_allele[1], method = method,
             t = rows$t[1], t_se = rows$t_se[1],
             component = c(rows$mediator, "Residual direct"),
             estimate = c(i, d), se = c(i_se, d_se),
             mediated = sum(i), mediated_se = sqrt(i_var))
}

share <- rbindlist(lapply(inventory$target_snp, decompose), fill = TRUE)
share[, `:=`(lo = estimate - z * se, hi = estimate + z * se,
             t_lo = t - z * t_se, t_hi = t + z * t_se,
             share_mediated = mediated / t)]
sensitivity_share <- rbindlist(lapply(names(mv_sets), decompose, fit_name = "sensitivity"))

fwrite(mvmr_table, file.path(results_dir, "mvmr_estimates.csv"))
fwrite(share, file.path(results_dir, "mediated_share.csv"))
fwrite(sensitivity_share, file.path(results_dir, "mediated_share_sensitivity.csv"))

# ---- Plot -----------------------------------------------------------------

totals <- unique(share[, .(target_snp, method, t, t_lo, t_hi, share_mediated)])
totals <- totals[order(-abs(t))]
totals[, label := ifelse(is.na(share_mediated), "no mediator",
  sprintf("%s%% mediated", formatC(100 * share_mediated, format = "f", digits = 0)))]
totals[, snp_label := paste0(target_snp, "\n", method)]
plot_dat <- merge(share, totals[, .(target_snp, snp_label)], by = "target_snp")
plot_dat[, snp_label := factor(snp_label, levels = totals$snp_label)]
totals[, snp_label := factor(snp_label, levels = totals$snp_label)]
plot_dat[, component := paste0(toupper(substr(component, 1, 1)), substring(component, 2))]
mediators <- setdiff(unique(plot_dat$component), "Residual direct")
plot_dat[, component := factor(component, levels = c(mediators, "Residual direct"))]
fills <- c(setNames(c("#1b7837", "#762a83", "#e08214", "#2166ac", "#b2182b")[seq_along(mediators)], mediators),
           "Residual direct" = "grey80")

span <- diff(range(c(0, totals$t_lo, totals$t_hi)))
totals[, label_y := ifelse(t >= 0, t_hi, t_lo) + sign(t) * 0.06 * span]

p <- ggplot(plot_dat, aes(x = snp_label, y = estimate, fill = component)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_errorbar(data = totals, aes(x = snp_label, ymin = t_lo, ymax = t_hi), width = 0.2,
                inherit.aes = FALSE) +
  geom_point(data = totals, aes(x = snp_label, y = t), inherit.aes = FALSE, size = 2) +
  geom_text(data = totals, aes(x = snp_label, y = label_y, label = label), inherit.aes = FALSE, size = 3) +
  scale_fill_manual(values = fills, name = NULL) +
  labs(x = NULL, y = "GBC log odds per target effect allele",
       caption = paste("Bars stack the indirect effect through each mediator and the residual direct effect;",
                       "together they sum to the total (point, 95% CI).",
                       "\nSegments of opposite sign to the total are negative mediation.",
                       "rs35604230 uses MVMR; others use univariable MR.")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank(),
        plot.caption = element_text(hjust = 0, size = 8))

ggsave(file.path(figures_dir, "phewas_mediated_share-g2.pdf"), p, width = 10, height = 6, device = cairo_pdf)
ggsave(file.path(figures_dir, "phewas_mediated_share-g2.png"), p, width = 10, height = 6, dpi = 200)

print(mvmr_table[, .(target_snp, analysis, mediator, b = signif(b, 3), se = signif(se, 3),
                     pval = signif(pval, 2), univariable_b = signif(univariable_b, 3),
                     conditional_F = round(conditional_F, 1), nsnp, Q_pval = signif(Q_pval, 2))])
print(totals[, .(target_snp, method, t = round(t, 3), share_mediated = round(share_mediated, 3))])
