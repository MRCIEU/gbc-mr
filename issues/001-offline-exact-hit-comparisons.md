# Correct exact-SNP PheWAS comparisons using saved GBC data

Type: AFK

## Parent PRD

`issues/prd-phewas-mediation.md`, sections A, B and F.

## What to build

Deliver a reproducible offline exact-hit comparison table and plot from `scripts/phewas-followup-g2.qmd`. Replace regional target estimates with Wald ratios using the exact GBC SNP association saved in `int_chd` and its outcome association from `outgwasf`. Rebuild each comparator from saved `int_inst`/outcome associations; deduplicate before analysis and retain explicit group metadata. Extract reusable analytical helpers where useful so labels no longer drive computation.

The current saved analysis has 53 comparisons but only 33 unique pairs, and none of its target-labelled regional instrument groups contains its labelled SNP. The full saved outcome contains all six required exact target SNPs; reduced `outdat` contains none.

## Acceptance criteria

- [ ] Recompute from `data/gbc_phewas_results-g2.Rdata` and/or other existing project empirical inputs, without association queries, new data downloads or dependence on old result caches.
- [ ] Preserve an inventory of the 33 unique saved SNP–mediator pairs, 26 mediators and six targets; reconcile the original eight-hit PheWAS inventory and explicitly account for hits/pairs not followed up.
- [ ] Every usable target estimate uses exactly the named SNP, harmonised against `outgwasf`, with a recorded effect allele, association SEs and mediator scale. Missing/ambiguous targets are recorded as unavailable rather than replaced silently.
- [ ] Comparator instruments exclude the target and the documented primary target-region window before MR; record excluded/retained counts and methods. Add locally feasible wider-window/MHC sensitivity and document coordinate/LD limitations.
- [ ] Use explicit mediator ID, target SNP and group-role metadata; never infer identity from display strings. Validate harmonisation, including the stated strand/palindrome policy.
- [ ] Produce a tidy per-pair table with target/comparator effects, CIs, instrument counts, Q and failure reasons, plus an updated comparison plot available for inspection.
- [ ] Meaningful checks cover duplicate removal, exact target membership, exclusion behavior and a known allele-reversal case; running the analysis with remote extraction disabled succeeds.

## Blocked by

None - can start immediately.

## User stories addressed

- User stories 1, 2 and 6.
