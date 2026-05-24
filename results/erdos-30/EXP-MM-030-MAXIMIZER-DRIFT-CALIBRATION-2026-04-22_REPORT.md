# EXP-MM-030-MAXIMIZER-DRIFT-CALIBRATION-2026-04-22 — Prefix Drift Calibration on Exact Maximizers

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-MAXIMIZER-DRIFT-CALIBRATION-2026-04-22 |
| Source packet | EXP-MM-030-MAXIMIZER-RIGIDITY-SMALLN-2026-04-22_RESULTS.json |
| Data integrity | REAL_COMPUTATION_DERIVED_FROM_EXACT_PACKET |
| Scan window | n = 10 through n = 50 |
| Exact maximizers covered | 76,368 |

## Probe Question

> On the exact maximizer packet, how much deterministic prefix drift is actually forced by the small-n data?

## Answer

The small-n data does not force the current deterministic drift nearly as hard as the safe theorem wrapper suggests. The present theorem-aligned drift improves the observed constant, but even the raw prefix discrepancy with no drift subtracted at all stays below the same literature-scale normalization on this window.

## Drift Comparison

| Drift ansatz subtracted from raw prefix discrepancy | Worst residual in n^(7/8) units | Worst n |
|---|---|---|
| 0 | 0.8927 | 13 |
| sqrt(n) | 0.5875 | 43 |
| ||A|-sqrt(n)| * sqrt(n) | 0.5334 | 10 |
| max(||A|-sqrt(n)|, 1) * sqrt(n) | 0.5303 | 16 |

## Interpretation

The important comparison is between the first and last rows. Subtracting the current theorem-aligned drift lowers the worst small-n constant from about 0.8928 to about 0.5304, so the wrapper is not useless. But the raw prefix discrepancy is already O(n^(7/8)) on all scanned exact maximizers. That shifts the next target: instead of asking for a bigger deterministic endpoint term, the right question is whether a sharper theorem can absorb some or all of the current drift back into the Balasubramanian-Dutta error scale.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-MAXIMIZER-DRIFT-CALIBRATION-2026-04-22_RESULTS.json | Structured drift-calibration results |
| EXP-MM-030-MAXIMIZER-DRIFT-CALIBRATION-2026-04-22_REPORT.md | Human-readable report |
| EXP-MM-030-MAXIMIZER-DRIFT-CALIBRATION-2026-04-22_RESULTS.sha256 | Integrity checksum |
