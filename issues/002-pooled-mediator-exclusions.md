# Correct pooled mediator MR by excluding every relevant GBC target locus

Type: AFK

## Parent PRD

`issues/prd-phewas-mediation.md`, sections A, C, E and F.

## What to build

Deliver a corrected pooled mediator-to-GBC MR table and supporting figure directly from saved associations. For each mediator, exclude its union of target loci before harmonisation and MR. The existing union of per-target complements reintroduces 32 excluded regional instrument rows across five traits, including the principal proposed mediators.

This pooled result is a distinct summary from the per-target comparator: do not substitute it into the exact-hit comparisons or decomposition.

## Acceptance criteria

- [ ] Start from each mediator's saved instruments and apply the complete target-locus union using the exclusion policy/metadata established in ticket 001.
- [ ] An instrument excluded for target A cannot return through target B. Report exclusions before and after outcome matching/harmonisation.
- [ ] Return one pooled result or explicit unavailable status per mediator, with SNP count, MR method, estimate, SE, CI, nominal p-value and BH-adjusted p-value across the named pooled-MR family.
- [ ] Update the notebook's pooled MR table and supporting forest figure, showing target exclusion counts and any changed candidate-mediator conclusions.
- [ ] Include a meaningful regression check using an existing mediator with multiple targets that demonstrates the old complement-union bug and validates the new set.
- [ ] Recompute without external data access or old knitr result caches. Compare the revised results with the saved original results as a diagnostic, not an expected-result oracle.

## Blocked by

- Blocked by `issues/001-offline-exact-hit-comparisons.md`.

## User stories addressed

- User stories 1, 3, 5 and 6.
