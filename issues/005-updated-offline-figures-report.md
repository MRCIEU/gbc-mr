# Export the updated GBC mediation figures and reproducible offline report

Type: AFK

## Parent PRD

`issues/prd-phewas-mediation.md`, sections A, E and F.

## What to build

Deliver the final updated figure and report from the corrected results: exact-hit versus other-instruments MR comparisons, total/indirect/residual-direct component summaries, and the revised pooled mediator MR summary. Integrate the individually verifiable outputs from tickets 001–004 into readable publication exports and an HTML report that runs entirely from saved local empirical data.

## Acceptance criteria

- [ ] Update `figures/phewas_followup_results-g2.pdf` and produce a PNG counterpart; export an additional decomposition PDF/PNG with total, indirect and residual direct estimates and 95% CIs.
- [ ] All unique retained pairs appear once in analytical summaries. Unavailable estimates are accounted for explicitly; supplemental or multipage layouts keep SNP/mediator labels and CIs readable.
- [ ] Clearly distinguish mediator-specific MR units in comparison panels from the common per-allele GBC log-odds units in decomposition panels. Effect orientation, instrument counts, inference labels and adjustment families are understandable from the figure/caption.
- [ ] Render `scripts/phewas-followup-g2.html` from the corrected local pipeline and link/export tidy comparison, pooled-MR and decomposition CSVs alongside figures.
- [ ] Document the local rendering command, required existing inputs and installed-package versions. Remove executable remote extraction from the G2 workflow and ensure outdated cached outputs cannot mask changed associations/functions.
- [ ] Validate a fresh render with networking/extraction disabled and old knitr caches bypassed; check exported counts against the canonical tables.
- [ ] Inspect the resulting HTML and PDF/PNG pages visually for clipped labels, overlapping points, misleading scales, missing CIs and duplicate rows; repair layout before delivery.
- [ ] Update README links/captions to the final exports and revised findings. Preserve existing G1 analytical behavior if shared helper extraction affects the older notebook.

## Blocked by

- Blocked by `issues/004-inference-and-interpretation.md`.

## User stories addressed

- User stories 1, 2, 3, 4, 5 and 6.
