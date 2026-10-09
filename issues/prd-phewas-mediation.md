# Offline PheWAS mediation follow-up for GBC

## Objective

Correct `scripts/phewas-followup-g2.qmd` so each GBC SNP–mediator pair is evaluated using the exact SNP and independent mediator instruments. Produce an updated comparison figure, a summary of total/indirect/residual direct effects, and a reproducible report.

## User stories

1. As an analyst, I can rerun the G2 analysis using only empirical data already in the project directory, without OpenGWAS or any replacement data service.
2. As an analyst, I can compare the exact GBC hit's mediator-to-GBC Wald estimate with MR using other mediator instruments, once per unique SNP–mediator pair.
3. As an analyst, I can view pooled mediator-to-GBC MR with every relevant GBC target locus excluded.
4. As an analyst, I can view total, indirect and residual direct effects and their uncertainty for every SNP–mediator–GBC relationship.
5. As a reader, I can distinguish evidence against a proposed mediation model from uncertain results, with appropriate multiple-testing summaries and stated assumptions.
6. As a reader, I can inspect an updated, readable figure and HTML report, with downloadable publication figures and result tables.

## A. Offline input contract

- Empirical inputs must already exist inside this repository. No OpenGWAS extraction, replacement association queries, new LD-reference downloads, or other newly fetched empirical datasets.
- Primary input: `data/gbc_phewas_results-g2.Rdata`, containing `outgwasf`, `int_inst`, `int_chd`, `phewas_result`, `outdat` and old `exp_sets`. The full outcome object and local raw G2 GWAS are available; the reduced `outdat` has no exact target SNPs.
- Validation of the saved Rdata found 33 unique SNP–mediator pairs across 26 mediators and six target SNPs, all six present with usable GBC effects in `outgwasf`; `int_inst` contains 8,549 rows. Old 53-row comparisons include duplicates.
- The supplied PheWAS CSV originally covers eight lead SNPs. Account explicitly for the two without retained saved follow-up pairs, rather than claiming all eight were analysed.
- Rebuild scientific outputs from saved association inputs; do not treat old `exp_sets`, `res_sets`, figures or knitr result caches as corrected results.
- Record local input provenance and package versions. Use currently available packages; do not make remote installation a prerequisite.

## B. Exact-hit comparisons

For each unique pair, convert the saved exact target association in `int_chd` to exposure format and harmonise it against its GBC association from `outgwasf`. Estimate its Wald ratio. Estimate mediator-to-GBC MR from other usable, saved mediator instruments, excluding the target and its prespecified local region. Retain explicit SNP, mediator and group identifiers; labels are presentation only.

Keep the original ±1 Mb exclusion as a documented primary geographic rule, check coordinate consistency using local data, and provide wider-window/MHC sensitivity where local data permit. Do not claim a geographic window proves LD independence, especially around rs35604230 in the MHC. Do not replace an exact hit with a regional instrument set or an unvalidated proxy.

Report instrument counts, methods, effect sizes, uncertainty and Q for available comparisons. Unavailable or unusable pairs remain in the result inventory with explicit reasons. Harmonisation must document its forward-strand assumption or use conservative palindrome handling; missing frequency information must not silently resolve ambiguous alleles. Preserve an explicit, consistent target effect allele for every estimate and display it in the output.

## C. Pooled mediator MR

For each mediator, exclude the union of all its target loci directly from the saved instrument list, then harmonise and estimate MR once. Never form this instrument set by unioning per-target complements. Show a revised table with exclusion counts, estimates, confidence intervals, nominal p-values and adjusted p-values.

## D. Effect decomposition

For a target SNP G, mediator M and GBC Y, let a = the harmonised G→M association, b = the M→Y MR estimate using other instruments, and t = the harmonised G→Y association. Report:

- total effect: t;
- indirect effect through M: i = a × b;
- residual direct effect relative to M: d = t − i.

Use the per-pair comparator from section B, not the pooled estimate from section C. Effects are per target effect allele on the GBC log-odds scale, and a/b must use the same mediator scale. Each mediator is evaluated separately; indirect effects across correlated mediators are not additive.

Provide SEs and 95% CIs for total, indirect and residual direct effects. Propagate uncertainty in a, b and t using a documented delta method or reproducible parametric simulation. For independent estimates, first-order Var(i) = b² Var(a) + a² Var(b), and Var(d) = Var(t) + Var(i); general formulas must include covariance terms when known. Saved summary statistics do not identify every covariance: flag any zero-covariance working assumption and the resulting approximate intervals. Interval estimation must not silently treat a or b as fixed. Where covariance is unknown, use a prespecified valid covariance sensitivity range if feasible and report it as sensitivity, without inventing empirical covariance.

Describe d as an estimated residual direct component under the stated model. Binary mediators/outcomes, non-collapsibility, instrument validity, population differences and nonlinearity limit a causal direct/indirect interpretation. Do not present these as identified natural effects, or calculate a mediation percentage by default when total effects are near zero/opposite in sign. Missing associations or comparators yield explicit unavailable components, never zero.

## E. Interpretation and multiplicity

Use unique tested pairs as the Q/residual-test family, and one result per mediator as the pooled-MR family. Report nominal and BH-adjusted p-values within each explicitly named family, with unavailable tests excluded and their count stated. Do not imply the families jointly control one global error rate.

Heterogeneity indicates incompatibility with the common-effect/full-mediation model under its assumptions; it does not uniquely locate a direct pathway. No detected heterogeneity does not establish mediation or equivalence. Report uncertainty and instrument validity/precision diagnostics rather than a binary 'mediated/not mediated' verdict. Update stale prose that calls these GBC associations very weak: the existing report has seven of eight lead variants at P < 5e-8, while rs2131242 is below genome-wide significance (P approximately 4.81e-7). Any general statement must identify its SNP/pair denominator and distinguish candidate explanations from causal proof.

## F. Updated figure and report

Update `figures/phewas_followup_results-g2.pdf` and add a PNG counterpart. Include an exact-hit versus other-instruments comparison and a separate total/indirect/residual-direct forest summary, with CIs, SNP/mediator labels and effect-allele units. Avoid mixing mediator-specific MR scales in the decomposition panel: its components share per-allele GBC log-odds units.

Render `scripts/phewas-followup-g2.html` from the corrected offline pipeline. Export tidy comparison, pooled-MR and decomposition CSVs with provenance, counts and missingness reasons. Figures must remain readable for all pairs through a suitable multipage layout or additional supplementary export. Provide a local rendering command and meaningful checks of the exclusions, uncertainty propagation and result counts. Verify the report without dependence on old knitr caches or network access.

## Scope

No new data acquisition, no new PheWAS, no multivariable mediation model requiring unsaved associations, no automatic mediator selection, and no summing separately estimated indirect effects across mediators. Existing G1 notebook behavior should be preserved if genuinely shared helpers are extracted.

## Methodological references

- User-provided CHD example: https://explodecomputer.github.io/lab-book/posts/2026-05-22-phewas-example/chd_example.html
- Carter et al., Mendelian randomisation for mediation analysis: current methods and challenges for implementation: https://doi.org/10.1007/s10654-021-00757-1
- de Leeuw et al., Understanding the assumptions underlying Mendelian randomization: https://pubmed.ncbi.nlm.nih.gov/35082398/

These references justify methodology; they are not new empirical inputs or required network calls during execution.
