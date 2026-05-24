# Track B Pilot Readiness — Observable Freeze

**Date:** 2026-05-05  
**Protocol:** `TRACK_B_PREREGISTRATION.json`  
**Scope:** QC ↔ #30 generic control and QC ↔ #20 sunflower scoped-positive-control pilot  
**Status:** Ready for run-script implementation; not yet executed.

## Readiness Verdict

Track B is ready to move from preregistration to implementation only if the run script obeys the frozen choices below. This file does **not** change the preregistration thresholds. It records the observable freeze required by the preregistration before code is written.

The pilot remains a Leg-4 experiment only after both exp10 and exp11 run under the frozen protocol:

- `exp10_qc_sidon_generic_control`: expected to confirm #30 as generic / geometry-dominated.
- `exp11_qc_sunflower_scoped_positive_control`: expected to test #20 as the scoped / geometry-exhausted positive control.

## Existing #20 Artifacts

The #20 side has two checked-in experiment families:

- `Erdos20/EXP-MATH-ERDOS20-SUNFLOWER-001_RESULTS.json`
  - w = 3, k = 3, n = 4..8.
  - Has full partition function and max-family-size sequence for n = 4..8.
  - Has growth-rate profiles only for n = 4..7.
- `Erdos20/EXP-MATH-ERDOS20-SUNFLOWER-002_RESULTS.json`
  - Cross-w sweep w = 2, 3, 4 at k = 3.
  - Confirms w = 2 saturation at M = 6.
  - w = 3 saturation at M = 12 remains tentative.
  - w = 4 is incomplete beyond n = 7.

## Frozen #20 Verdict Observable

**Primary verdict observable:** max-family-size sequence

```text
M(n) = largest sunflower-free family size for w = 3, k = 3, n ∈ {4,5,6,7,8}
```

**Reason for choosing M(n):** the preregistration domain includes n ∈ {4,5,6,7,8}. Growth-rate profiles `g_n(m)` are preferred in the preregistration, but the existing #20 artifact lacks growth rates for n = 8. Choosing `g_n(m)` now would either drop n = 8 or require new computation before the first Leg-4 run, both of which would violate the "freeze before code" discipline. `M(n)` is the complete existing observable over the registered domain.

**Secondary diagnostic only:** growth-rate profile `g_n(m)` for n ∈ {4,5,6,7}. This can diagnose whether the max-family signal is a real jamming transition, but it cannot decide the Leg-4 verdict in exp11 unless a dated preregistration amendment explicitly changes the primary observable.

## Registered Geometry-Killed Nulls for exp11

The run script must implement all three nulls and must not add or drop nulls after seeing results.

1. **Occupancy-scale null**
   - Preserve the ambient number of w-subsets and family-size m.
   - Destroy sunflower core-closure geometry by assigning random forbidden triples matched only on count.

2. **Low-order marginal null**
   - Preserve pairwise intersection-size distribution approximately.
   - Randomize triple-level core alignment so actual sunflower closure structure is broken.

3. **State-space-size nuisance null**
   - Match `C(n,w)` and the observed maximum feasible family size scale.
   - Replace sunflower closure with a nuisance exclusion rule sampled from the same ambient state count.

## exp11 Verdict Discipline

Use the preregistered Track A v2 criteria without tuning:

- **C1 strength ratio:** M(n)-derived signal must separate from all three nulls by the frozen Track A threshold.
- **C2 scaling:** M(n) must follow the scoped/geometry-exhausted scaling behavior, not the generic geometric class.
- **C3 fixed-knot model selection:** the M_L fixed-knot model must be strongly preferred on #20 and not equally preferred on at least one null.

The only allowed verdict labels remain:

- `LEG4_PASS`
- `LEG4_FAIL`

No "near pass", "suggestive pass", or scoped-bucket promotion is allowed from narrative interpretation alone.

## exp10 Coupling

Track B still requires the paired #30 generic-control run. The #30 side should use the already-compiled Sidon strict-upper / sharp-difference observable from `Erdos30_SharpDiff.lean` as the generic control anchor, with the expectation that the fixed-knot M_L model fails to distinguish #30 from at least one geometry-killed null.

## Immediate Implementation Target

Create `exp11_qc_sunflower_leg4_run.py` with:

- deterministic seed,
- JSON output,
- report output,
- SHA-256 checksum,
- frozen observable extraction from the existing #20 result JSONs,
- all three nulls above,
- Track A v2 model-selection code reused rather than rewritten ad hoc.

The implementation may read the existing #20 JSON artifacts. It may not recompute or replace them during the Leg-4 run.
