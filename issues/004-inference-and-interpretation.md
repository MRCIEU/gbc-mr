# Report unique comparisons, multiplicity and uncertainty without overclaiming mediation

Type: AFK

## Parent PRD

`issues/prd-phewas-mediation.md`, sections B–F.

## What to build

Deliver a coherent inferential summary in the notebook, README and analytical figures, covering exact-hit heterogeneity, pooled mediator effects and effect decomposition. Replace statements that identify a direct pathway from Q alone or equate lack of heterogeneity with demonstrated mediation. Base summaries on unique pairs and report candidate explanations alongside uncertainty and instrument limitations.

## Acceptance criteria

- [ ] Define and name the unique-pair heterogeneity family, pooled-mediator family and any additional residual/indirect testing families; apply BH correction within each family and show nominal and adjusted p-values. Report tested/unavailable counts and avoid implying global control across separate families.
- [ ] Compute denominators from unique IDs, retaining both per-pair and per-SNP summaries; repeated exposure rows do not increase counts. Summaries reconcile to the retained input inventory.
- [ ] Present Q as incompatibility with the common-effect/full-mediation model under assumptions, with possible target/comparator instrument invalidity and LD/population differences acknowledged.
- [ ] Report non-significant comparisons as no detected heterogeneity or unresolved, rather than proven mediation/equivalence. Include precision and available instrument-strength diagnostics for candidate pathways.
- [ ] Interpret total/indirect/residual direct estimates with their CIs and covariance/model assumptions. Distinguish approximate components from identified causal natural effects.
- [ ] Update stale weak-GBC-association prose and qualify any 'mostly not' statement with its explicit denominator and evidence definition. Keep scientific narrative and plot legends consistent with computed results.
- [ ] Update README and figure annotations so platelet-count evidence and other candidate pathways reflect revised exclusions and multiplicity rather than the original unadjusted table.
- [ ] A representative available, unavailable, imprecise and incompatible relationship can each be traced from its local inputs to the corresponding table, interpretation and visual annotation.

## Blocked by

- Blocked by `issues/001-offline-exact-hit-comparisons.md`.
- Blocked by `issues/002-pooled-mediator-exclusions.md`.
- Blocked by `issues/003-direct-indirect-effects.md`.

## User stories addressed

- User stories 2, 3, 4 and 5.
