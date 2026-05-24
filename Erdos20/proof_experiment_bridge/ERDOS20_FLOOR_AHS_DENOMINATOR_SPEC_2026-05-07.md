# Erdos #20 Floor / AHS Denominator Specification

**Date:** 2026-05-07
**Worker:** #20-C
**Lane:** sunflower core-closure proof/experiment bridge — denominator specification
**Scope:** specification only. No compute, no enumerator change, no Lean edit, no scorecard touch, no D1 write, no public-doc update, no theorem statement.
**Claim ceiling:** A0.
**Required preserved phrase:** shadow signature, not universal law.

This packet pins down two co-primary normalized closure ratios — `floor_ratio` and `ahs_ratio` — precisely enough that a downstream analysis script can convert any per-core summary row produced by Worker A's enumerator into both columns, without writing that script today and without reissuing existing JSON files. The specification is additive: every key it introduces is optional in `erdos20.safe_candidates.v1` and required in `erdos20.safe_candidates.v1.1`.

This packet does not advance, improve, resolve, prove, or establish anything about the sunflower conjecture. It specifies how to read the channel measurement, not what the channel says.

## Inputs Read

- `erdos-experiments/Erdos20/SUNFLOWER_CORE_CLOSURE_LEG4_PACKET_2026-05-05.md` — defines `I_core_local`, the Mendoza floor `E_L = I * k_B * T * ln 2`, PASS/FAIL/INCONCLUSIVE classifier, and the rule that the floor numerator must be defined before the run.
- `erdos-experiments/Erdos20/third_party_review/ERDOS20_THIRD_PARTY_REVIEW_BASELINE_SAFE_CANDIDATES_2026-05-06.md` — promotes Abbott-Hansen-Sauer (AHS) normalization from secondary control to co-primary reporting and requires `floor_ratio` and `ahs_ratio` side-by-side.
- `erdos-experiments/Erdos20/proof_experiment_bridge/ERDOS20_SAFE_CANDIDATES_BRIDGE_PACKET_2026-05-06.md` — defines the `erdos20.safe_candidates.v1` enumerator schema (per-family, per-core, per-candidate rows) and the `safeCandidates` Lean predicate target.
- `erdos-experiments/Erdos20/Q1_LITERATURE_GATE_2026-04-17.md` — records that the AHS lower bound is `f(w, 3) >= 10^(w/2 - O(log w))`, hence `c_3 >= sqrt(10) ~= 3.162`, and that our small-n c_3 measurements are strictly weaker than AHS.
- `erdos-experiments/Erdos20/proof_experiment_bridge/EXP-MATH-ERDOS20-SAFE-CANDIDATES-BRIDGE-20260506-01_RESULTS.json` — sampled to confirm the actual layout of `samples_by_m[*]` rows including the existing fields `I_core_local_bits`, `candidate_count`, `I_AHS_bits`, `I_floor_bits`, `ahs_ratio`. The May 6 run already shipped a *preliminary* `I_floor_bits = 1.0` and a global `I_AHS_bits = log2(sqrt(10))`. This spec replaces both with sharper definitions and disambiguates per-core-size scaling.
- `erdos-experiments/Erdos20/EXP-MATH-ERDOS20-PER-CORE-CLOSURE-20260505-01_RESULTS.json` — sampled to confirm that the older per-core JSON does not carry these denominator fields, which is precisely why the analysis script must accept missing fields as `null`.

## What Each Denominator Means

`I_floor` and `I_AHS` answer different questions, and conflating them was the mistake the third-party review caught.

`I_floor` asks: **is the measured per-core closure channel above the Landauer-style information floor at the chosen scratch temperature**. The floor bit budget is the unconstrained-channel maximum — how many bits would be needed to record an arbitrary admissible extension through the same core if no sunflower constraint were imposed at all. The ratio `floor_ratio = I_core_local / I_floor` tells us how far the constraint has compressed the unconstrained channel. A low `floor_ratio` near zero means the sunflower constraint is barely biting at this `(s, m)`; a `floor_ratio` near 1 means the channel is operating near saturation of the floor; values >1 are unphysical for an information-theoretic ratio of this construction and would flag a definitional error.

`I_AHS` asks: **is the measured per-core closure channel large compared to the bit-cost of placing a new safe extension in a known sunflower-free construction at the same parameters**. Abbott-Hansen-Sauer gives `f(w, 3) >= 10^(w/2 - O(log w))` and thus `c_3 >= sqrt(10) ~= 3.162` — a construction-side lower bound on the extremal family size. The AHS-implied per-core information budget is the bits needed to record which AHS petal slot a new safe extension occupies, given the AHS template at parameters `(n, w)`. The ratio `ahs_ratio = I_core_local / I_AHS` tells us whether our channel measurement is rebadged trivia — if `ahs_ratio < 1`, the measured channel is below the construction baseline and any "shadow" we are seeing is no stronger than what AHS already supplies. If `ahs_ratio` is stable and `>= 1` across a late window, the channel is at least carrying as much information as the AHS construction does at the same `(n, w, s, m)`.

The two denominators exist for different reasons. `I_floor` is the Leg-4 thermodynamic-floor classifier (the Mendoza-floor probe). `I_AHS` is the literature-baseline calibration that prevents us from declaring a discovery that AHS already accounts for. Reporting them side-by-side is the point.

## Formal Definition: I_floor(s, m)

Inputs (per-core summary row at fixed `(n, w, s, m)`):

- `C(s, m)` := `candidate_count` — number of candidate extensions through any core of size `s` at family size `m`. In the v1 enumerator schema this is the `candidate_count` field on each `core_rows[*]` row, summed or averaged over cores of that size as the run-design specifies.
- `L(s, m)` := `local_safe_count` — number of locally-safe candidate extensions (no 3-sunflower closed through that exact core).
- `I_core_local(s, m)` := `I_core_local_bits` — `-log2( L(s, m) / C(s, m) )` when `0 < L < C`. This field already exists in v1.

Closed-form spec for the floor numerator:

```
I_floor(s, m) := log2( C(s, m) )
T_landauer    := 1.0   (dimensionless scratch; declared with the run)
floor_ratio   := I_core_local(s, m) / I_floor(s, m)
```

Justification: `log2(C(s, m))` is the bit count of an unconstrained random extension through a fixed core — the maximum information that a saturated boolean channel through `C(s, m)` candidate slots can carry. It is the right denominator for "how much of the unconstrained channel did the sunflower constraint actually consume" because `I_core_local` is itself measured in those same bits (it is `-log2` of an empirical probability over the same `C(s, m)` candidate slots). The `T_landauer = 1.0` choice fixes the dimensionless scratch regime; the Mendoza floor `E_L = I * k_B * T * ln 2` reduces to `I * ln 2` of energy per measured bit when `k_B = 1` and `T = 1`, so the ratio `I_core_local / I_floor` equals the energy ratio in those units. Higher-T regimes (if a future run wants to test temperature scaling explicitly) are handled by overriding `T_landauer` in the run config and reporting the new value with the row.

The May 6 bridge run shipped `I_floor_bits = 1.0` as a placeholder ("one predeclared bit of local closure information; preliminary Leg-4 denominator"). That is *not* the spec defined here. The May 6 placeholder must be re-derived in the analysis script using `log2( C(s, m) )` from the row's `candidate_count`. If the analysis script encounters an old row carrying `I_floor_bits = 1.0` from the May 6 placeholder definition, it should overwrite it using the spec here, not trust the cached value.

Edge cases (mandatory `null` returns, no silent coercion):

- `C(s, m) <= 1` (zero or one candidate, so `log2(C) <= 0`): emit `I_floor_bits = null` and `floor_ratio = null`. Tag basis as `null_too_few_candidates`.
- `L(s, m) = C(s, m)` (all candidates safe, so `I_core_local = 0`): emit `floor_ratio = 0.0` (the channel is genuinely carrying zero closure information at this `(s, m)`), keep `I_floor_bits` defined. Do not emit `null`.
- `L(s, m) = 0` (no candidate is locally safe, so `I_core_local` formally `+inf`): emit `I_core_local_bits = null` upstream and `floor_ratio = null`. Tag basis as `null_no_safe_candidates`.
- `floor_ratio` is never coerced to 1 or 0 to fill a missing value. `null` is the only allowed missing-value marker.

## Formal Definition: I_AHS(s, m)

The AHS construction is asymptotic. At small `(n, w)` it does not yet supply a meaningful per-core safe-extension count. The spec must therefore be `null`-safe.

Inputs:

- The same per-core summary row at fixed `(n, w, s, m)`.
- A function `M_AHS(n, w, s, m)` := count of safe extensions through a core of size `s`, evaluated inside the AHS-template family at parameters `(n, w)` at family size `m`. This count is *not* taken from our enumerator; it is computed from the AHS construction template directly. The construction is the standard `X(f) = {(x, f(x)) : x in [w]}` blow-up family used in the Q1 gate (and refined by AHS to the `10^(w/2 - O(log w))` bound). The exact reference enumeration of safe extensions through a fixed core inside that template is a small auxiliary computation owned by the analysis script; it is not in scope for this spec to compute it.

Closed-form spec for the AHS numerator:

```
I_AHS(s, m) := max( 0.0, log2( max( 1, M_AHS(n, w, s, m) ) ) )
ahs_ratio   := I_core_local(s, m) / I_AHS(s, m)
```

Small-`(n, w)` carve-out (mandatory):

The AHS construction is asymptotic; for `n < n_AHS_min(w)`, the AHS template either does not embed in `[1..n]` at all, or the safe-extension count through a fixed `s`-core is degenerate and not informative. In that regime, the analysis script must emit:

```
I_AHS_bits = null
ahs_ratio  = null
```

and tag `denominator_basis_ahs = "null_below_AHS_minimum"`.

Concrete `n_AHS_min(w)` defaults (subject to revision when the analysis script is written):

| w | `n_AHS_min(w)` | Reason |
|---|----------------|--------|
| 2 | infinity | AHS asymptotic, w=2 not used for AHS comparison; classical M(infinity, 2, 3) = 6 already covers w=2 sanity. |
| 3 | 8 | `n=5,6,7` calibration runs populate `floor_ratio` only; AHS comparison begins at `n >= 8`. |
| 4 | 7 | AHS template embeds at `n >= 7` for w=4; matches existing run reach. |
| >=5 | (define when the run reaches it) | Out of current compute scope. |

Calibration runs at `w=3, n=5` (and any other regime under `n_AHS_min(w)`) populate `floor_ratio` only. They are calibration, not Leg-4 AHS evidence. The downstream analysis script must enforce that an entire run is tagged `ahs_below_minimum_calibration_only` if every row has `ahs_ratio = null`.

Coarse-grained fallback:

If AHS template enumeration of `M_AHS(n, w, s, m)` is too expensive at chosen `(n, w)` (the per-core safe extension count inside the AHS template is itself a non-trivial enumeration), the analysis script may use the construction lower bound:

```
M_AHS_LB(w)              := floor( c_3^w )    where c_3 := sqrt(10) ~= 3.162
core_count(s, n, w)      := number of size-s cores in the run (already in the v1 schema)
M_AHS_LB_per_core(s, n, w) := max( 1, floor( M_AHS_LB(w) / core_count(s, n, w) ) )
I_AHS_coarse(w, s, n)    := max( 0.0, log2( M_AHS_LB_per_core(s, n, w) ) )
ahs_ratio_coarse         := I_core_local(s, m) / I_AHS_coarse(w, s, n)
```

The coarse-grained ratio is reported as `ahs_ratio_coarse` and tagged `denominator_basis_ahs = "ahs_lower_bound_coarse"`. It is *not* a substitute for `ahs_ratio`. A run that reports only `ahs_ratio_coarse` cannot PASS Leg-4. It can be INCONCLUSIVE.

## JSON Schema Patch (additive, schema_version bump to `erdos20.safe_candidates.v1.1`)

The patch is purely additive at the per-core summary row level. Old `v1` JSON files remain valid: the analysis script treats missing fields as `null`. No candidate-level rows are touched. No family-row or witness-pair fields are touched.

Per-core summary row, exact keys to add:

```json
{
  "I_floor_bits":               <float|null>,
  "floor_ratio":                <float|null>,
  "I_AHS_bits":                 <float|null>,
  "ahs_ratio":                  <float|null>,
  "ahs_ratio_coarse":           <float|null>,
  "denominator_basis_floor":    "log2_candidate_count_unconstrained",
  "denominator_basis_ahs":      "ahs_template_per_core_safe_count" | "ahs_lower_bound_coarse" | "null_below_AHS_minimum",
  "T_landauer":                 1.0
}
```

Run-level header additions:

```json
{
  "schema_version": "erdos20.safe_candidates.v1.1",
  "denominator_spec_ref": "proof_experiment_bridge/ERDOS20_FLOOR_AHS_DENOMINATOR_SPEC_2026-05-07.md"
}
```

Worker A's enumerator does not need to be modified to emit these fields. The analysis script consumes a `v1` per-core JSON, computes the seven new fields from existing `candidate_count`, `local_safe_count`, `I_core_local_bits`, and `core_count` (or from a passed-in AHS template enumeration), and writes a sibling `*_FLOOR_AHS_OVERLAY.json` carrying the `v1.1` header and the per-core rows extended with the new keys. Worker A's original JSON is left untouched.

## PASS / FAIL / INCONCLUSIVE Reading

The Leg-4 packet's PASS criteria, restated in terms of the two co-primary ratios:

PASS requires:

- `floor_ratio` is defined (not `null`) across at least three consecutive `(n, m)` values in the late window for at least one fixed core size `s`.
- `ahs_ratio` is defined (not `null`) and `>= 1` across the same late window — the construction-normalized ratio does not collapse below the literature baseline.
- Both ratios are stable: late-window coefficient of variation `<= 0.25` for `floor_ratio` and `<= 0.25` for `ahs_ratio`, separately.
- Log-log slope of either ratio against `N(n, w) = C(n, w)` has absolute value `<= 0.15` across the accepted window.
- The same qualitative pattern (floor and AHS together) appears in at least two core sizes or two `w` values, unless a pre-registered reason restricts to one core size.

FAIL is triggered if any of:

- `ahs_ratio < 1` consistently in the late window (channel is below construction baseline; signal is rebadged trivia).
- Either ratio drifts monotonically with absolute log-log slope `> 0.25`.
- The signal disappears under the AHS template-true `M_AHS` (only `ahs_ratio_coarse` looks favorable, true `ahs_ratio` does not).
- Per-core data contradicts the aggregate reading.

INCONCLUSIVE is triggered if any of:

- Run is entirely below `n_AHS_min(w)` and `ahs_ratio` is `null` for every row.
- Only `ahs_ratio_coarse` is available (template enumeration not yet done).
- `floor_ratio` is defined but `ahs_ratio` is not, or vice versa, across the late window.
- Window has fewer than three consecutive `(n, m)` rows.

The May 6 bridge run, evaluated under this spec, is INCONCLUSIVE: its `I_floor_bits = 1.0` placeholder does not match the `log2(C)` definition, and its global `I_AHS_bits = log2(sqrt(10))` is not a per-core, per-`(s, m)` AHS template count. The May 6 run is reclassified as preliminary calibration once the analysis script overlays `v1.1` fields onto it.

## What This Spec Does Not Do

- Does not advance, prove, resolve, improve, or establish the sunflower conjecture.
- Does not improve any Erdos-Rado or Abbott-Hansen-Sauer lower bound.
- Does not classify the May 5 per-core data, the May 6 bridge data, or any prior #20 run as Leg-4 PASS evidence.
- Does not authorize any public claim, social post, scorecard upgrade, A-axis change, or NotebookLM upload.
- Does not define a Lean theorem. The Lean theorem track (the `safeCandidates` predicate and the filtered-subset theorem) is owned by Worker B and tracked by the May 6 bridge packet, not by this spec.
- Does not modify Worker A's enumerator code, the Rust runner, or any existing per-core JSON file. The schema patch is additive and produced by a downstream overlay script.
- Does not commit any data to git or D1.

## Files Changed

- Created: `erdos-experiments/Erdos20/proof_experiment_bridge/ERDOS20_FLOOR_AHS_DENOMINATOR_SPEC_2026-05-07.md`

No other files written. No existing files edited.

## Next Concrete Step

A separate downstream analysis script — *not written by this spec packet* — will read a per-core JSON in either the existing `erdos20.safe_candidates.v1` format (or any earlier per-core RESULTS shape) and emit a sibling `*_FLOOR_AHS_OVERLAY.json` carrying the `v1.1` header and the seven new keys per per-core summary row defined above. Suggested location and name: `erdos-experiments/Erdos20/proof_experiment_bridge/run_floor_ahs_overlay.py`. Inputs: a `v1`/`v1.1` per-core JSON path. Outputs: a sibling JSON adding `I_floor_bits`, `floor_ratio`, `I_AHS_bits`, `ahs_ratio`, `ahs_ratio_coarse`, `denominator_basis_floor`, `denominator_basis_ahs`, `T_landauer` to every summary row, and a `denominator_spec_ref` pointer in the run header. The script must enforce the `null`-safe edge cases above and tag every row with the chosen AHS basis. It must not edit Worker A's JSON in place. It must not infer Leg-4 PASS/FAIL/INCONCLUSIVE — that classification is a separate consumer of the overlay file. The status remains A0: shadow signature, not universal law.
