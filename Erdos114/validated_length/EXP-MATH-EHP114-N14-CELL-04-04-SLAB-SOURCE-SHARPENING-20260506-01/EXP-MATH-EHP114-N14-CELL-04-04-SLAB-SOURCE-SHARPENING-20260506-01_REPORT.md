# EHP114 n=14 CELL-04-04 Slab Source Sharpening Diagnostic

Experiment: `EXP-MATH-EHP114-N14-CELL-04-04-SLAB-SOURCE-SHARPENING-20260506-01`

## Verdict

- Status: `SLAB_SOURCE_SHARPENING_DIAGNOSTIC_FAIL_OVER_BUDGET`
- Source cell: `CELL-04-04`
- Source total upper: `30.986477958908345`
- Exact cap: `20.672796062619668`
- Margin to cap: `-10.313681896288678`
- Unresolved branches: `4400`
- Worst branch length: `5.636765226771244`
- Worst branch slope: `6613.804457148043`
- Top-10 branch length sum: `17.503032997028853`
- Source-chain recommendation: `targeted_high_slope_slab_repair_before_downstream_pipeline`

## Interpretation

This diagnostic audits the first non-hard-cell slab source before downstream local-pipeline promotion. The key question is whether the over-budget bound is localized to a few high-slope/low-denominator slabs or spread globally. It does not certify a cell and does not promote length.

## Claim Ceiling

CELL-04-04 slab source-sharpening diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
