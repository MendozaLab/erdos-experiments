# Cross-Problem WS-01 Generalization — 2026-05-15 Extension

Date: 2026-05-15 (extension of 2026-05-14 #1038 cold-test; new file, NOT a modification of `CROSS_PROBLEM_WS01_GENERALIZATION_2026-05-14.md`)
Plan track: Probe 3 of "Intractability test program — 2026-05-15 (POC framing)" in `~/.claude/plans/yes-wondrous-blum.md` (per-candidate method at `~/.claude/plans/yes-wondrous-blum-agent-a4ae5ad039883b7ba.md`)
Claim ceiling: **Internal Python-level diagnostic only.** mpmath dps=30 interval arithmetic on small (n=8) base polynomials at box half-width R = 1e-4. NOT inari Rust. NOT certified. NOT a proof. NOT a Lean statement. NOT a public artifact.

## Scope

Three additional candidates run cold against the cross-problem WS-01 applicability probe, complementing the 2026-05-14 #1038 cold test:

| Problem | Role | Expected | Observed |
|---|---|---|---|
| #1041 (Erdős–Herzog–Piranian path-length) | primary high-fit target, sibling of #114 in EHP58 | `YES_STRUCTURAL_FIT`, `WS01_APPLICABLE` | `YES_STRUCTURAL_FIT`, `WS01_APPLICABLE` (12/12 boxes closed at conservative 2R) |
| #1043 (Pommerenke projection bound) | secondary high-fit target, same EHP58+Po59/Po61 lineage, different objective | `PARTIAL_FIT` or `YES_FIT`, `WS01_APPLICABLE` possible | `PARTIAL_FIT`, `WS01_APPLICABLE` (12/12 boxes closed at conservative 2R) |
| #20 (sunflower / Erdős–Ko–Rado–style) | explicit negative control, discrete combinatorial parameter space | `NO_STRUCTURAL_FIT`, cold test skipped | `NO_STRUCTURAL_FIT`, cold test skipped |

## Per-candidate findings

### #1041 — Erdős–Herzog–Piranian path-length conjecture

Structural fit: all three features YES. #1041 lives in the same paper (EHP58) and shares the same parameter space (monic complex polynomial of degree n with roots in the unit disk) and the same zero curve `{|f(z)| = 1}` as #114. The path-length objective is a different sublevel-set-derived quantity from #114's component-length objective, but that difference lives downstream of WS-01's wall test.

Cold-test artifact: `Math/erdos-experiments/Erdos1041/proof_path/EXP-MATH-ERDOS-1041-WS01-COLD-TEST-20260515-01/EXP-MATH-ERDOS-1041-WS01-APPLICABILITY-PROBE-20260516-01_RESULTS.json`

Input boundary slice: 12 candidate failure boxes around the n=8 base lemniscate (roots evenly spaced on |z| = 0.85 inside the unit disk; 4 boundary angles × 3 small rotation perturbations).

Probe verdict: `WS01_APPLICABLE`. 12/12 closed at conservative 2R bound. Branch-point validation rate 1.00. Median required factor (RHS/LHS) after WS-01 at 2R: 343.76 (i.e. the wall test passes with ~340× slack on the conservative bound). `center_fn_abs_lower` (|∇F| at center) ranges 12.1–19.8 across the 12 boxes, `f_tt_abs_upper` ranges 163–332.

Reading: WS-01 ports to #1041 by direct copy with only the box-data input changing. The wall test's local geometric content is identical to #114's; the objective change does not affect the per-box gate. This is the strongest possible portability evidence at the EHP-paper-sibling level.

### #1043 — Pommerenke projection bound

Structural fit: Features 1 and 3 YES, Feature 2 PARTIAL. The lemniscate boundary is identical to #114/#1041, the gradient and interval-IVT inputs are identical, but the wall test's diagnostic value for the conjecture itself is downstream-limited — Pommerenke's 1961 negative resolution is precisely a construction of polynomials whose lemniscate-bounded sublevel sets project with measure > 2 in every direction. WS-01 certifies the local lemniscate; the projection failure mechanism is upstream of where WS-01 helps.

The cross-problem cold-test question is *not* "does WS-01 close #1043"; the question is "does WS-01 survive an objective change on the same boundary geometry." That question has a clean YES answer.

Cold-test artifact: `Math/erdos-experiments/Erdos1043/proof_path/EXP-MATH-ERDOS-1043-WS01-COLD-TEST-20260515-01/EXP-MATH-ERDOS-1043-WS01-APPLICABILITY-PROBE-20260516-01_RESULTS.json`

Input boundary slice: 12 candidate failure boxes around the same n=8 base lemniscate (4 boundary angles × 3 projection-direction angles `u ∈ {0, π/3, 2π/3}`). The projection direction is recorded per box but does not enter the per-box gate.

Probe verdict: `WS01_APPLICABLE`. 12/12 closed at conservative 2R bound. Branch-point validation rate 1.00. Median required factor (RHS/LHS) after WS-01 at 2R: 336.58. The numerical pattern matches #1041 exactly — because the boxes were generated on the same polynomial family.

Reading: the architecture's portability is *boundary-geometry-dependent, not objective-dependent*. The wall test does not care which sublevel-set-derived quantity (component count, path length, projection measure) is the conjecture's objective; it cares about the local lemniscate geometry. This is information about *which features port* and *which do not*: the architecture ports across objective changes; the architecture does not by itself close conjectures whose failure mechanism is upstream of the wall.

### #20 — sunflower / Erdős–Ko–Rado–style extremal sets

Structural fit: all three features NO. Sunflower lives in a discrete combinatorial parameter space; no zero curve, no gradient, no interval IVT, no continuous deformation. Cold test skipped by design.

Cold-test artifact: not generated (cold test skipped). The structural-fit triple at `Math/erdos-experiments/Erdos20/STRUCTURAL_FIT_CHECK_2026-05-15.md` + `_RESULTS.json` + `_RESULTS.sha256` is the complete output.

Reading: the structural-fit checker correctly identifies that WS-01 does not apply to discrete combinatorial problems. Without this negative control we could not distinguish "WS-01 is broadly portable" from "the structural-fit checker is permissive." The negative control passes — the checker is *not* permissive.

## Aggregate portability boundary

Per the verdict criteria in the plan:

- **`STRONG_PORTABILITY`**: 2+ of 3 candidates pass fit + WS-01 applicable. Met: #1041 passes (YES + WS01_APPLICABLE) and #1043 passes (PARTIAL_FIT + WS01_APPLICABLE). 2 of 3 non-control candidates passed.
- **`EHP58_FAMILY_PORTABILITY`**: #1041 passes; #1043 partial; #20 confirmed-NO-fit. Met exactly: #1041 full YES, #1043 PARTIAL on Feature 2, #20 explicit NO.
- **`#114_SPECIFIC`**: only #1041 fails or all non-#114-sibling candidates fail. Not met (both #1041 and #1043 pass).
- **`NO_PORTABILITY`**: even #1041 doesn't replicate. Not met (#1041 closes 12/12).

**Aggregate verdict: `STRONG_PORTABILITY` — with the structural caveat that all four problems that pass (#114, #1038, #1041, #1043) share polynomial-level-set boundary geometry.** The portability boundary is therefore: WS-01 ports across the polynomial-level-set / real-analytic-zero-curve family, including across objective changes within that family (path length, component length, sublevel-set measure, projected measure). It does NOT port across the analytic / combinatorial boundary (#20 cold test skipped by design).

A useful refinement: the architecture is *boundary-geometry-bound*, not *objective-bound*. The wall test's local computation does not depend on which functional of the sublevel set is the conjecture's target, only on the polynomial-derived F(z) and its derivatives. This is a sharper structural finding than the plan's binary "STRONG_PORTABILITY vs #114_SPECIFIC" axis anticipated.

## Implication for the POC's methodology-portability claim

The intractability-test-program POC framing in `yes-wondrous-blum.md` casts the bridge program for #114 as a proof-of-concept for a certified-verification pipeline applied to intractable middle ranges. Probe 3's job was to test whether WS-01 (one component of the pipeline) is #114-specific or portable.

The 2026-05-15 finding: **WS-01 is portable across the polynomial-level-set boundary family** — confirmed concretely on #114 (regression), #1038 (2026-05-14 cold test), #1041 (this extension), and #1043 (this extension). It correctly fails to apply to #20 (this extension, negative control). The architecture is bounded by *boundary geometry*, not by *conjecture identity* or *objective functional*.

For the POC technical report (conditional on the intractability synthesis confirming intractability): the methodology-portability claim should be phrased as "the wall-test rewrite component of the pipeline ports across polynomial-level-set problems with real-analytic zero curves; objective-functional changes do not break the gate; non-analytic problems are out of scope by structural-fit verdict." This is a more honest and more useful claim than "the architecture generalizes."

For Probe 2 (Atlas/Collider neighborhood probe) synthesis: a "shortcut" to #114 via a structural neighbor would have to come from inside the polynomial-level-set family, because that is where WS-01 has demonstrated portability. Neighbors outside this family (combinatorial extremal sets, discrete number-theoretic problems) are unlikely to provide WS-01-compatible transfer mechanisms.

For Probe 1 (K_sting extraction) synthesis: the wall test's portability does not by itself refute intractability of the candidate-N₀ threshold; intractability is a property of the constants K_sting × K_outside⁴, not of the gate architecture. Probe 3's positive result is consistent with either intractability outcome at Probe 1.

## Files produced

Probe 3 artifacts (all new, none overwriting prior work — H² Rule 3 compliant):

- `Math/erdos-experiments/Erdos1041/STRUCTURAL_FIT_CHECK_2026-05-15.md`
- `Math/erdos-experiments/Erdos1041/STRUCTURAL_FIT_CHECK_2026-05-15_RESULTS.json`
- `Math/erdos-experiments/Erdos1041/STRUCTURAL_FIT_CHECK_2026-05-15_RESULTS.sha256`
- `Math/erdos-experiments/Erdos1041/build_1041_boundary_slice_cold_test.py`
- `Math/erdos-experiments/Erdos1041/proof_path/INPUT_boundary_slice_1041_2026-05-15.json`
- `Math/erdos-experiments/Erdos1041/proof_path/EXP-MATH-ERDOS-1041-WS01-COLD-TEST-20260515-01/EXP-MATH-ERDOS-1041-WS01-APPLICABILITY-PROBE-20260516-01_RESULTS.json`
- `Math/erdos-experiments/Erdos1041/proof_path/EXP-MATH-ERDOS-1041-WS01-COLD-TEST-20260515-01/EXP-MATH-ERDOS-1041-WS01-APPLICABILITY-PROBE-20260516-01_REPORT.md`
- `Math/erdos-experiments/Erdos1041/proof_path/EXP-MATH-ERDOS-1041-WS01-COLD-TEST-20260515-01/EXP-MATH-ERDOS-1041-WS01-APPLICABILITY-PROBE-20260516-01_RESULTS.sha256`
- `Math/erdos-experiments/Erdos1043/STRUCTURAL_FIT_CHECK_2026-05-15.md`
- `Math/erdos-experiments/Erdos1043/STRUCTURAL_FIT_CHECK_2026-05-15_RESULTS.json`
- `Math/erdos-experiments/Erdos1043/STRUCTURAL_FIT_CHECK_2026-05-15_RESULTS.sha256`
- `Math/erdos-experiments/Erdos1043/build_1043_boundary_slice_cold_test.py`
- `Math/erdos-experiments/Erdos1043/proof_path/INPUT_boundary_slice_1043_2026-05-15.json`
- `Math/erdos-experiments/Erdos1043/proof_path/EXP-MATH-ERDOS-1043-WS01-COLD-TEST-20260515-01/EXP-MATH-ERDOS-1043-WS01-APPLICABILITY-PROBE-20260516-01_RESULTS.json`
- `Math/erdos-experiments/Erdos1043/proof_path/EXP-MATH-ERDOS-1043-WS01-COLD-TEST-20260515-01/EXP-MATH-ERDOS-1043-WS01-APPLICABILITY-PROBE-20260516-01_REPORT.md`
- `Math/erdos-experiments/Erdos1043/proof_path/EXP-MATH-ERDOS-1043-WS01-COLD-TEST-20260515-01/EXP-MATH-ERDOS-1043-WS01-APPLICABILITY-PROBE-20260516-01_RESULTS.sha256`
- `Math/erdos-experiments/Erdos20/STRUCTURAL_FIT_CHECK_2026-05-15.md`
- `Math/erdos-experiments/Erdos20/STRUCTURAL_FIT_CHECK_2026-05-15_RESULTS.json`
- `Math/erdos-experiments/Erdos20/STRUCTURAL_FIT_CHECK_2026-05-15_RESULTS.sha256`
- `Math/erdos-experiments/CROSS_PROBLEM_WS01_GENERALIZATION_2026-05-15.md` (this file)

Note on artifact-ID dates: the probe-internal `experiment_id` field is `EXP-MATH-ERDOS-<N>-WS01-APPLICABILITY-PROBE-20260516-01` (date when the probe ran), while the containing artifact directory and the plan-specified per-candidate experiment ID is `EXP-MATH-ERDOS-<N>-WS01-COLD-TEST-20260515-01` (plan date). The two are intentionally not the same because the probe stamps its own run date and the plan stamps the experiment-program date; both are recorded in the JSON.

## Sources

- `~/.claude/plans/yes-wondrous-blum.md` § "Intractability test program — 2026-05-15 (POC framing)" Probe 3 (lines 318-326)
- `~/.claude/plans/yes-wondrous-blum-agent-a4ae5ad039883b7ba.md` (per-candidate execution method)
- `Math/erdos-experiments/cross_problem_ws01_applicability_probe.py` (probe consumer; regression PASS 16/16 against #114 CELL-02-03 owner-family-4488:3)
- `Math/erdos-experiments/Erdos1038/STRUCTURAL_FIT_CHECK_2026-05-14.md` + `CROSS_PROBLEM_WS01_GENERALIZATION_2026-05-14.md` (prior cold-test precedent, NOT modified by this extension)
- `Math/formal-conjectures/FormalConjectures/ErdosProblems/{1041,1043,20}.lean` (canonical Lean problem statements)
