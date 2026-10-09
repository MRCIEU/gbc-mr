# PheWAS + MR analysis of gall bladder cancer GWAS

The G2 follow-up asks whether each GBC lead SNP's mediator association is compatible
with MR from other instruments for the same mediator. The corrected analysis uses
only saved local associations and explicitly separates exact-hit comparisons,
pooled mediator MR and model-based effect decomposition.

[Read the offline G2 report](scripts/phewas-followup-g2.html). There are **33 unique
SNP–mediator pairs, 26 mediators and six retained targets**, rather than 53 repeated
comparison rows. The [eight-hit inventory](results/g2/lead_inventory.csv) accounts
for rs2312959 and rs77155094, which have no retained saved follow-up pairs. Seven
of eight original GBC leads meet P < 5e-8; the saved outcome P for rs2131242 is
9.862e-7 (the older report quoted approximately 4.81e-7).

![GBC GWAS](figures/gwas_manhattan.png)

The source PheWAS is [the supplied eight-lead-SNP CSV](data/phewas/10_phewas_results_of_8_lead_snps_GBC.csv),
with candidate GWAS selected in [the saved trait list](data/gbc-traits-from-phewas.csv).
No new empirical associations or LD references are fetched.

The exact target is harmonised by its named SNP against the full saved GBC outcome.
Other-instruments MR excludes that SNP and its inclusive ±1 Mb region. Pooled MR
instead excludes the union of **all** of a mediator's target loci before
harmonisation. The [regression diagnostic](results/g2/pooled-complement-union-diagnostic.csv)
shows that the former union-of-complements approach reintroduced 32 excluded rows
across five mediators. All palindromic variants are conservatively removed;
coordinate mismatches, missingness and instrument-strength diagnostics are recorded.
Geographic exclusions do not establish LD independence, especially in the MHC.

**32 of 33 unique pairs** have BH-adjusted heterogeneity P < 0.05. This indicates
incompatibility with a common-effect/full-mediation model under its assumptions;
it does not locate a direct pathway. The remaining pair has no detected
heterogeneity, with mediation unresolved. Invalid instruments, LD and population
or scale differences remain alternative explanations.

**None of the 26 pooled mediator estimates passes BH at 0.05.** Platelet count has
estimate 0.211 (SE 0.0894), nominal P = 0.0182 and BH P = 0.0948 after eight regional
instruments are excluded. Cholelithiasis and gallbladder/biliary/pancreatic disorders
have pooled BH P approximately 0.0512. These are candidate explanations with
uncertainty, not demonstrated mediating pathways. The named testing families are
corrected separately, without a claim of global control across families.

- [Exact-hit comparison PDF](figures/phewas_followup_results-g2.pdf), [PNG first page](figures/phewas_followup_results-g2.png), [tidy CSV](results/g2/comparisons.csv).
- [Total/indirect/residual-direct PDF](figures/phewas_decomposition-g2.pdf), [PNG first page](figures/phewas_decomposition-g2.png), [tidy CSV](results/g2/decomposition.csv).
- [Pooled MR PDF](figures/phewas_pooled_mr-g2.pdf), [tidy CSV](results/g2/pooled_mr.csv), [legacy comparison](results/g2/pooled_legacy_diagnostic.csv).
- [Window/MHC sensitivity](results/g2/window_sensitivity.csv), [covariance sensitivity](results/g2/covariance_sensitivity.csv), [family denominators](results/g2/inference_families.csv), [per-SNP summary](results/g2/per_snp_summary.csv).

The report embeds every PNG page. Comparison panels use mediator-specific saved
GWAS units and separate scales; exact units/transformations were not preserved.
Decomposition components share GBC log-odds units per recorded target effect allele:
total t, indirect i = a × b, and model-based residual direct d = t − i. Its
approximate delta-method CIs propagate uncertainty in a, b and t, assuming zero
covariances, with hypothetical positive-definite covariance sensitivities supplied.
Binary-trait non-collapsibility and model assumptions limit causal interpretation.
Separate mediators' indirect effects must not be summed, and no mediation percentage
is reported.

Run from the repository root with the existing local `data/gbc_phewas_results-g2.Rdata`:

```bash
Rscript scripts/run-offline-g2.R
quarto render scripts/phewas-followup-g2.qmd --no-cache
Rscript tests/run-tests.R
```

For a checkout without the original 215 MB local Rdata, a compact lossless subset
of its required saved associations is included. It retains 8,549 saved instrument
rows, all saved exact associations and 7,281 matching outcome rows including all
eight original leads. [Snapshot provenance](results/g2/input-snapshot-provenance.csv)
records source/snapshot checksums and counts. Regenerate it from the original using
`Rscript scripts/make-offline-g2-snapshot.R`.

```bash
Rscript scripts/run-offline-g2.R data/offline-g2-inputs.rds
GBC_G2_INPUT=data/offline-g2-inputs.rds quarto render scripts/phewas-followup-g2.qmd --no-cache
```

Required installed tools: R, TwoSampleMR and its dependencies, Quarto, knitr and
rmarkdown. There is no remote installation step. [Package versions](results/g2/package_versions.csv)
and [session information](results/g2/session-info.txt) document the render environment.
All results are recomputed with caches disabled and remote APIs/downloads blocked.
Tests validate exact membership, duplicate/conflicting rows, allele reversal,
exclusions, pooled union regression, unavailable inputs and decomposition uncertainty.
The older [G1 workflow](scripts/phewas-followup.qmd) remains unchanged.
