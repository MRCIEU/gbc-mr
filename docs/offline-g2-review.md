# Integrated review of issues 1–5

Fixed base: `67fe24519222a1d6228d9d722df35569d33319cf` (original `main`).
Diff: `git diff 67fe24519222a1d6228d9d722df35569d33319cf...HEAD`.
Specifications: `issues/001-*.md` through `issues/005-*.md`, plus the approved PRD.

## Standards

The independent standards reviewer found no actionable documented-standard
violations, Fowler-baseline smells or clear bugs. The analytical helpers share
one offline association/harmonisation domain; the pipeline is an orchestration
boundary. Names, recorded identifiers, scales and assumptions are explicit.

## Spec

The independent specification reviewer found no actionable requirements or
scientific correctness defects. Exact targets, complete pooled locus exclusion,
consistent allele orientation, per-target decomposition, uncertainty propagation,
PSD covariance sensitivity, separate BH families and original-hit reconciliation
match the approved specifications. No unrequested scientific scope was found.

Total findings: Standards **0**; Spec **0**. Neither axis has an outstanding issue.

## Verification

- All five test files passed with `Rscript tests/run-tests.R`.
- The portable pipeline test passed again after final export/metadata refinements.
- R source parsing and Git whitespace checks passed.
- A fresh Quarto render succeeded using the compact saved-association snapshot,
  `--no-cache`, blocked remote APIs/downloads and unusable network proxies.
- Exported inventories reconcile: 33 comparison/decomposition pairs, 26 pooled
  mediators, six retained targets, eight original leads, 66 window-sensitivity
  rows and 99 covariance-sensitivity rows.
- The old pooled complement-union regression reintroduces 32 excluded instrument
  rows across five mediators; the corrected sets exclude them.
- All 14 figure pages were inspected; PDF fonts were embedded to preserve spacing,
  and label/CIs, signed estimates, scales and page counts were checked.
- HTML source checks confirm 14 embedded PNG pages, two 33-row pair tables,
  working local download links and no remote script sources.

Browser security policy refused a `file:` URL for the HTML preview. Interactive
browser appearance was therefore not verified; the generated report, embedded
figure pages and local assets were checked without bypassing that policy.

## Scientific interpretation

All 33 empirical pairs are usable under the primary policy. 32/33 pair Q tests
are BH-significant, indicating incompatibility with the common-effect model
under its assumptions. This does not uniquely identify direct pathways.
No pooled mediator passes BH at 0.05 (0/26); platelet count has nominal
P=0.0182 and BH P=0.0948. Approximate residual/indirect components remain
model-dependent; saved summary data do not identify empirical covariances.
