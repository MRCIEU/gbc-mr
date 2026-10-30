# Methods

## Gallbladder cancer GWAS

Summary statistics for gallbladder cancer (GBC) were taken from the [G1 + G2 meta-analysis; sample size and ancestry to be added], with effects on the log-odds scale. Eight lead variants were carried forward: rs67274221, rs35604230, rs60541816, rs6733452, rs2312959, rs15734, rs2131242 and rs77155094. Seven of these reached P < 5 × 10⁻⁸; rs2131242 had P = 9.9 × 10⁻⁷.

## PheWAS

We looked up each lead variant in the IEU OpenGWAS database (Elsworth et al., 2020) for associations with other traits at P < 1 × 10⁻⁵. This returned 423 associations across 370 GWAS datasets. From these we chose 26 traits as candidate mediators of the GBC association, covering haematological, lipid, liver-function, anthropometric, blood-pressure, endocrine and gallbladder/biliary traits. Each trait was represented by one OpenGWAS dataset.

## Follow-up Mendelian randomization

For each candidate mediator we used OpenGWAS to extract genome-wide significant instruments (P < 5 × 10⁻⁸). These were LD-clumped at r² < 0.001 within a 10 Mb window, using the 1000 Genomes European reference panel. We also extracted the association of each GBC lead variant with each mediator. A lead-variant–mediator pair was kept if that exact variant reached P < 5 × 10⁻⁸ for the mediator; proxies were not used. This gave 33 unique pairs, covering 26 mediators and six lead variants. Neither rs2312959 nor rs77155094 reached P < 5 × 10⁻⁸ for any selected mediator.

We harmonised mediator and GBC associations to the same effect allele. All palindromic (A/T, C/G) variants were removed, whatever their allele frequency. For each pair we estimated two things:

- **Exact-hit estimate.** The Wald ratio of the lead variant's GBC association to its mediator association.
- **Other-instrument estimate.** An MR estimate of the mediator's effect on GBC, using all other instruments for that mediator. We removed every instrument within ±1 Mb of any of the eight GBC lead variants. This estimate was shared by all pairs for that mediator. We used inverse-variance weighted (IVW) MR with multiplicative random effects, or a Wald ratio when only one instrument remained.

If the lead variant affects GBC only through the mediator, the two estimates should agree. We tested this with a one-degree-of-freedom Cochran's Q statistic, Q = (β_exact − β_other)² / (SE²_exact + SE²_other). We recorded instrument strength (F = β²/SE²) for every instrument. As a sensitivity analysis, we widened the exclusion to ±2 Mb, with and without also excluding the extended MHC (chr6:25–34 Mb). The pooled mediator→GBC IVW estimates are reported for each mediator.

## Direct and indirect effect estimation

For each lead variant G, mediator M and GBC outcome Y (Figure 1), let:

- *a* = the G→M association;
- *b* = the M→Y other-instrument MR estimate above;
- *t* = the G→Y association.

All three use the same effect allele. We estimated:

- the total effect, *t*;
- the indirect effect through M, *i* = *a* × *b*;
- the residual direct effect, *d* = *t* − *i*.

All three are expressed in GBC log odds per effect allele. Standard errors came from the first-order delta method: Var(*i*) = *b*²Var(*a*) + *a*²Var(*b*) and Var(*d*) = Var(*t*) + Var(*i*). The covariances between *a*, *b* and *t* can't be identified from summary data, so the primary analysis set them to zero. As a sensitivity analysis, we recomputed intervals assuming equal correlations of −0.3, 0 and 0.3. Each mediator was analysed on its own, and indirect effects were not summed across mediators. Because GBC is binary and odds ratios are non-collapsible, we treat *d* as a model-based residual. It should not be read as an identified natural direct effect.

## Effect captured through mediators

To summarise how much of each lead variant's effect runs through its mediators, we took as candidates the mediators whose pooled MR estimate had P < 0.05. We treated cholelithiasis and FinnGen disorders of gallbladder, biliary tract and pancreas as one liability, represented by cholelithiasis. A lead variant with one candidate mediator kept its decomposition from above. For a lead variant with several candidate mediators (rs35604230: platelet count, rheumatoid arthritis and diastolic blood pressure), we estimated the mediators' effects on GBC jointly by multivariable MR (MVMR).

MVMR instruments were extracted from OpenGWAS with TwoSampleMR's `mv_extract_exposures()`. This took the union of each mediator's instruments (P < 5 × 10⁻⁸), clumped them jointly (r² < 0.001, 10 Mb) and looked up each SNP in every mediator GWAS, without proxies. As before, we removed palindromic variants and all variants within ±1 Mb of any GBC lead. Conditional effects were estimated by IVW multivariable regression with multiplicative random effects. Instrument strength was assessed with the Sanderson–Windmeijer conditional F statistic, assuming no covariance between the mediator GWAS (Sanderson et al., 2021). As a sensitivity analysis, we also excluded variants within ±2 Mb of a lead and the extended MHC.

The effect through the set of mediators was Σ *a*ⱼ*b*ⱼ, where *b*ⱼ are the conditional MVMR estimates, so correlated mediators are not double-counted. The residual direct effect was *d* = *t* − Σ *a*ⱼ*b*ⱼ, and the mediated share was Σ *a*ⱼ*b*ⱼ / *t*. Standard errors used the delta method with the full MVMR covariance matrix of *b*, treating the *a*ⱼ as independent of each other and of *b*. Because mediators were selected on their association with GBC in the same data, these shares are exploratory.

## Multiple testing

P-values were adjusted with the Benjamini–Hochberg method, separately within each family of tests: the 33 pairwise heterogeneity tests, the 26 pooled mediator MR estimates, the 33 indirect effects and the 33 residual direct effects.

## Software

Data were extracted from OpenGWAS with ieugwasr and TwoSampleMR in August–September 2026, and MVMR instruments in October 2026. All analyses were rerun offline from these saved extracts, using R 4.6.1, TwoSampleMR 0.7.12 (Hemani et al., 2018), ieugwasr 1.2.0, data.table 1.18.4, knitr 1.51, rmarkdown 2.31 and Quarto 1.9.37. Code and results are available at [repository URL].

## Figure 1

![Analysis design](../figures/analysis-design.png)

**Figure 1. Analysis design.** Each GBC GWAS hit is paired with a mediator identified in the PheWAS. *a* is the hit's association with the mediator. *b* is the mediator's effect on gallbladder cancer, estimated by MR using other instruments for the mediator, with all regions within ±1 Mb of a GBC lead variant excluded. *d* is the residual direct effect of the hit that does not act through the mediator. The total association of the hit with gallbladder cancer, *t*, decomposes as *t* = *d* + *a* × *b*.

## References

- Elsworth B, Lyon M, Alexander T, et al. The MRC IEU OpenGWAS data infrastructure. *bioRxiv* 2020. doi:10.1101/2020.08.10.244293
- Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports systematic causal inference across the human phenome. *eLife* 2018;7:e34408. doi:10.7554/eLife.34408
- Sanderson E, Spiller W, Bowden J. Testing and correcting for weak and pleiotropic instruments in two-sample multivariable Mendelian randomization. *Stat Med* 2021;40:5434–5452. doi:10.1002/sim.9133
- Carter AR, Sanderson E, Hammerton G, et al. Mendelian randomisation for mediation analysis: current methods and challenges for implementation. *Eur J Epidemiol* 2021;36:465–478. doi:10.1007/s10654-021-00757-1
