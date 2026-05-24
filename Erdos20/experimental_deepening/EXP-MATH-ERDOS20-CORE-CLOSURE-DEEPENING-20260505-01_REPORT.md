# EXP-MATH-ERDOS20-CORE-CLOSURE-DEEPENING-20260505-01 Report

## Status

Secondary deterministic diagnostic from saved Erdos #20 artifacts only. No new family enumeration, no D1, no scorecard, no public docs.

## Verdict

- Classification: `DEEPENING_GEOMETRY_DRIFT_PRESENT`
- Claim ceiling: `A-axis A0; shadow signature, not universal law; no theorem, no lower-bound progress, no Leg-4 pass`
- Meaning: existing data support a measurable per-core channel, but not a Leg-4 pass or theorem progress.

## Drift Screen

| Series | points | slope vs N | late CV | class |
|---|---:|---:|---:|---|
| w=2 strongest per-core | 5 | 0.625853 | 0.123378 | GEOMETRY_DRIFT |
| w=3 strongest per-core | 4 | 0.061889 | 0.106989 | STABLE_BY_SLOPE_ONLY |
| w=4 strongest per-core | 2 | None | None | INSUFFICIENT_POINTS |
| w=2, fixed core size s=0 | 5 | 0.98886 | 0.123378 | GEOMETRY_DRIFT |
| w=2, fixed core size s=1 | 5 | -0.676349 | 0.120621 | GEOMETRY_DRIFT |
| w=3, fixed core size s=2 | 4 | -0.076613 | 0.182257 | STABLE_BY_SLOPE_ONLY |

## Interpretation

The w=3 fixed-core-size series is the useful next lane because it now spans n=5..8 after the Rust W3N8 artifact. The screen is still only geometry-drift analysis: it does not define a Mendoza-floor numerator, compare Abbott-Hansen-Sauer constructions, or show a theorem beyond encoding.

## Next Measurable

Pre-register fixed core-size local closure cost `I_core_local(s,m*)` with a defect/hysteresis window around `m*`, then run a symmetry-reduced or Rust exact sweep where feasible. The immediate target is not a sunflower bound; it is whether the per-core channel stays stable once ordinary ambient geometry is accounted for.

## Source Boundary

Inputs are listed with SHA-256 hashes in the companion results JSON.
