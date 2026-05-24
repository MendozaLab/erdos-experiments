# EHP114 n=14 CELL-02-06 Slab Source Sharpening Diagnostic

Experiment: `EXP-MATH-EHP114-N14-CELL-02-06-SLAB-SOURCE-SHARPENING-20260506-01`

## Verdict

- Status: `SLAB_SOURCE_SHARPENING_DIAGNOSTIC_FAIL_OVER_BUDGET`
- Source cell: `CELL-02-06`
- Source total upper: `21.11393239317112`
- Exact cap: `20.672796062619668`
- Margin to cap: `-0.4411363305514513`
- Unresolved branches: `4401`
- Worst branch length: `3.6906218429575897`
- Worst branch slope: `4330.329513605012`
- Top-10 branch length sum: `8.60602095967236`
- Source-chain recommendation: `targeted_high_slope_slab_repair_before_downstream_pipeline`

## Interpretation

This diagnostic audits the first non-hard-cell slab source before downstream local-pipeline promotion. The key question is whether the over-budget bound is localized to a few high-slope/low-denominator slabs or spread globally. It does not certify a cell and does not promote length.

## Claim Ceiling

CELL-02-06 slab source-sharpening diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
