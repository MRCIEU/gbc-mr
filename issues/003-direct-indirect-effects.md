# Estimate total, indirect and residual direct effects for each GBC SNP–mediator pair

Type: AFK

## Parent PRD

`issues/prd-phewas-mediation.md`, sections A, B, D and F.

## What to build

Deliver an effect-decomposition table and forest plot for every saved SNP–mediator–GBC relationship, using the exact target association and its per-target other-instruments MR estimate. On the same target effect-allele orientation, estimate total t = beta(G→GBC), indirect i = beta(G→mediator) × beta(mediator→GBC), and residual direct d = t − i. Follow PRD section D for units, uncertainty and interpretation.

## Acceptance criteria

- [ ] Join by explicit SNP and mediator IDs to the exact-hit results from ticket 001, using the per-target comparator rather than pooled MR from ticket 002.
- [ ] Report t, i and d on a common per-effect-allele GBC log-odds scale, recording the effect allele and mediator units/scale. Confirm t = i + d for every complete pair.
- [ ] Provide SEs and 95% CIs with uncertainty in both multiplicands and the total effect propagated using a documented delta method or seeded parametric simulation.
- [ ] Document available covariance information and flag zero-covariance working assumptions; implement a clearly labelled covariance sensitivity where feasible without inventing empirical data.
- [ ] Preserve each relationship in the inventory; missing/unusable components have explicit reasons and are never imputed as zero.
- [ ] Label the direct component as a model-based residual direct estimate; explain binary-mediator/outcome and non-collapsibility limitations. Do not sum indirect effects across mediators or report unstable mediation percentages by default.
- [ ] Export a tidy decomposition table and readable component forest plot; retain signed opposing-path estimates rather than clipping them.
- [ ] Validate an interpretable complete-mediation case, a residual-direct case, an opposing-path case and allele reversal; demonstrate that increasing input uncertainty broadens the resulting intervals. Use local or synthetic validation inputs, with no newly fetched empirical data.

## Blocked by

- Blocked by `issues/001-offline-exact-hit-comparisons.md`.

## User stories addressed

- User stories 1, 4 and 6.
