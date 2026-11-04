# Shared comparator follow-up

The user's correction on 2026-10-09 supersedes the original per-target comparator
policy in issues 1–5:

> for any mediator, i want all known GBC hits to be removed, so that means that
> the 'other snps' estimate is always the same.

> in the phewas follow up results plot there can happily be one plot per mediator,
> with one 'other snps' estimate and one 'exact hit' estimate for each SNP that hit
> that mediator.

Exclude the complete saved GBC lead inventory (eight original leads, plus any
retained exact targets) for every mediator. Preserve the existing inclusive
±1 Mb geographic exclusions around every known lead, ±2 Mb sensitivity and MHC
sensitivity policy. Include the two original leads without retained follow-up
pairs. Missing or conflicting lead coordinates make the comparator unavailable.

Compute one comparator per mediator and reuse it across exact-hit Q tests and
model-based decomposition. Pooled MR must agree with this comparator. Keep the
33 unique pair tests and 26 pooled tests in their existing separate BH families.
Redraw comparisons with one panel per mediator, one other-SNP estimate, and a
labelled exact estimate and Q BH P value for each retained target. Regenerate
all downstream tables, figures, interpretation and the offline HTML report.
