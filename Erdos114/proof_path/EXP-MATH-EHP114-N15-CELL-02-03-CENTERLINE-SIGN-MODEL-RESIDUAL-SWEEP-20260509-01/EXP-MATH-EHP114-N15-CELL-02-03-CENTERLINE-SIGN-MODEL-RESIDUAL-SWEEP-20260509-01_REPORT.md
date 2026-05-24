# EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-RESIDUAL-SWEEP-20260509-01

## Scope

Internal sweep artifact only. Track B1 residual-owner-family sweep of EHP114
bridge program. Not a proof of Erdos #114, not an n=15 certificate, not a
CELL-02-03 closure. Diagnostic of WS-01-CENTER-STRIP-CANCELLATION (Branch-
Centered Moving-Frame Collar rewrite) on the 5 owner families of CELL-02-03's
108 wall-separation failures that were not already covered by the 4488:3 /
2484:4 / 4571:2 runs. Python numerical demonstration at the same claim level
as the per-family runs; NOT Rust interval-certified.

## Pattern choice

Pattern B (per-family artifacts, mirroring the existing 2484-4 / 4571-2 layout)
for the actual WS-01 runs, plus this Pattern A umbrella aggregate to give a
single ledger across all 108 CELL-02-03 wall failures. The 5 residual families
were small enough that per-family granularity is more useful than a monolithic
sweep file, and the umbrella here pulls all of it together.

## Status

`RESIDUAL_SWEEP_PARTIAL`

## Residual-sweep headline numbers (5 owner families, 65 boxes)

- Closed by WS-01 (conservative 2R bound): `27`
- Closed by WS-01 (tight R only): `37`
- Still failing: `38`
- Branch point not found: `0`
- Median required factor (RHS/LHS) before WS-01 (overall residual): `0.005520414975142834`
- Median required factor (RHS/LHS) after WS-01 at 2R (overall residual): `0.8439888666232092`

## Per-family table (residual sweep)

| owner_key | failures | closed (2R) | closed (tight R only) | still_failing | branch_not_found | outcome |
|-----------|----------|-------------|------------------------|----------------|-------------------|---------|
| `3104:0` | 16 | 0 | 16 | 16 | 0 | `BLOCKED` |
| `3546:0` | 16 | 6 | 10 | 10 | 0 | `PARTIAL` |
| `3982:0` | 16 | 4 | 11 | 12 | 0 | `PARTIAL` |
| `4404:0` | 16 | 16 | 0 | 0 | 0 | `FULL_PASS` |
| `2932:1` | 1 | 1 | 0 | 0 | 0 | `FULL_PASS` |

## Per-family table (prior runs, for context only — read-only)

| owner_key | failures | closed (2R) | closed (tight R only) | still_failing | branch_not_found | outcome |
|-----------|----------|-------------|------------------------|----------------|-------------------|---------|
| `4488:3` | 16 | 16 | 0 | 0 | 0 | `FULL_PASS` |
| `2484:4` | 16 | 16 | 0 | 0 | 0 | `FULL_PASS` |
| `4571:2` | 11 | 4 | 1 | 7 | 0 | `PARTIAL` |

## Aggregate after today (all 108 CELL-02-03 wall failures)

- Total failures: `108`
- Closed by WS-01 at conservative 2R: `63`
- Closed by WS-01 at tight R only: `38`
- Still failing: `45`
- Branch point not found: `0`
- Closure rate at conservative 2R: `58.33%`
- Owner families full-pass: `['4404:0', '2932:1', '4488:3', '2484:4']`
- Owner families partial: `['3546:0', '3982:0', '4571:2']`
- Owner families blocked: `['3104:0']`

## Interpretation

Across the 5 residual owner families (65 boxes), WS-01 closes 27 boxes by the conservative 2R bound and an additional 37 boxes by the tight-R bound; 38 remain failing and 0 have no validated branch point. Two of the five (4404:0, 2932:1; 17 boxes total) are uniform full-pass — those rows have R*|F_n|/|F_t|-style geometry where the rewrite gives a large factor of safety. The remaining three (3104:0, 3546:0, 3982:0; 48 boxes total) are dominated by R^2*|F_tt| where the conservative 2R-radius penalty overshoots the wall RHS by ~30-50%, but the tight-R sensitivity figure clears comfortably (~99% of the 48 close at tight R). Combined with prior runs (4488:3=16/16, 2484:4=16/16, 4571:2=4/11 at 2R), the aggregate across all 108 CELL-02-03 wall failures is 63/108 = 58.3% closed at conservative 2R. This is a Python numerical demonstration, NOT a Rust interval-certified result.

## Next dependency (residual sweep)

Partial pass on residual sweep. WS-01 generalizes but is not uniform at the conservative 2R bound. Three owner families (3104:0, 3546:0, 3982:0) need either a tighter R bound (interval-certified z* via interval Newton, which would replace 2R with a true off-center radius) or a refined T3/F_tt bound. Recommended next: produce a tight-R (z*-anchored) variant of the WS-01 closure for these three families and verify whether the tight-R closure (which already clears) survives interval-Newton certification.

## Next dependency (aggregate, all 108 boxes)

Aggregate closure rate < 70%. WS-01 alone is insufficient for CELL-02-03. A different rewrite or substantive cell subdivision is needed before Track B1 can be claimed as a CELL-02-03 closure path.

## Source provenance

Per-family residual artifacts (each with its own SHA-256 sidecar in its own folder):
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3104-0-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3104-0-20260509-01_RESULTS.json`
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3546-0-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3546-0-20260509-01_RESULTS.json`
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3982-0-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3982-0-20260509-01_RESULTS.json`
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4404-0-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4404-0-20260509-01_RESULTS.json`
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2932-1-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2932-1-20260509-01_RESULTS.json`

Prior (read-only) artifacts:
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01_RESULTS.json`
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2484-4-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2484-4-20260509-01_RESULTS.json`
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4571-2-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4571-2-20260509-01_RESULTS.json`

Math reference (consume, do not re-derive):
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_2026-05-06.md`
  lines 122-146 (Branch-Centered Moving-Frame Collar lemma).

## Safety

This artifact is review-only. Forbidden writes (D1, curated morphisms, public
pages, proof registries, existing EHP114 finite packets, any of the per-family
centerline-sign-model artifacts) were not attempted. The aggregate script
`build_n15_cell0203_residual_sweep_aggregate.py` writes only this experiment
folder. No prior script or artifact was modified.
