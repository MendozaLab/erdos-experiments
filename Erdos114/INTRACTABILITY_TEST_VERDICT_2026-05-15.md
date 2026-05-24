# EHP #114 Intractability Test — Final Verdict

**Date:** 2026-05-15 (synthesis written 2026-05-16 after probe completion)
**Status:** Closed test program. Decision: POC framing adopted. Bridge program retired.
**Claim ceiling:** This document records the outcome of the three-probe intractability test program. It is not a proof of intractability in the mathematical sense; it is an empirical demonstration that the tools at hand cannot reach the middle range, and the structural argument why.

---

## What was tested

The intractability of the middle range `n ∈ {15, …, N₀ − 1}` for the Erdős–Herzog–Piranian conjecture (#114) under the H² portfolio's current tooling, specifically:

- The Tao-effectivization chain (constants extracted from Tao 2512.12455v2 lemmas)
- The Erdős Atlas + Collider's morphism neighborhood of #114
- The WS-01 (Branch-Centered Moving-Frame Collar) wall-separation architecture across multiple problems

Three parallel probes ran 2026-05-15:

1. **Probe 1** (K_sting extraction): `EXP-MATH-EHP114-TAO-KSTING-EXTRACTION-20260515-01`
2. **Probe 2** (Atlas+Collider neighborhood): `EXP-MATH-EHP114-ATLAS-COLLIDER-NEIGHBORHOOD-PROBE-20260515-01`
3. **Probe 3** (WS-01 cross-problem): `EXP-MATH-ERDOS-{1041,1043,20}-WS01-COLD-TEST-20260515-01` + aggregate at `CROSS_PROBLEM_WS01_GENERALIZATION_2026-05-15.md`

## Probe 1 — K_sting extraction outcome

**Verdict:** `KSTING_BOUNDED_BY_CASE_ANALYSIS`, candidate upper bound `≤ 100` (tight `≤ 30`).

The Stokes-theorem chain in Tao's Lemma `sting` (lines 1446–1466 of arXiv:2512.12455v2) decomposes as `K_Ψ = 2` from line 1448, `K_X1 ≤ 20π` from the Markov polynomial inequality, `K_X2 ≤ π` from the Riesz / multip bounds, `K_X3 ≤ 2π` from pp-bound, `K_X4 = 0`, and `K_X5 ≤ 2π` from cosec-asymptotic. Total composed bound ≈ 80.5; rounded to 100 with dyadic-series + cosec-asymptotic buffer. None of the constituents depends on `c_10`, `K_OR`, or any non-effective input — confirming Tao's line-196 effectivization remark empirically for this constant.

Plugging into the candidate-`N₀` formula `n / log n > 3.85 × 10¹¹ · K_sting · K_outside⁴` with `K_outside ≤ 1000`:

| K_sting | candidate N₀ |
|---:|---:|
| 100 (honest) | `2.4 × 10²⁷` |
| 30 (tight) | `7.1 × 10²⁶` |
| 1 (impossible floor) | `2.2 × 10²⁵` |

Even at the impossible `K_sting = 1` floor — well below any plausible composition — the candidate threshold is `~10²⁵`. Thirteen orders of magnitude above the `10¹²` intractability bar and nineteen above the `10⁸` "plausibly reachable" bar.

**The bottleneck is not `K_sting`.** It's the `c_10⁻⁵ ≈ 6.45 × 10¹⁰` amplification (from the rigorous lower bound `c_10 ≥ 0.0069`) combined with the `K_outside⁴` multiplier. Bounding `K_sting` numerically does not move the threshold meaningfully. The intractability is structural, not constants-driven.

## Probe 2 — Atlas + Collider neighborhood outcome

**Verdict:** `INTRACTABILITY_CONFIRMED_BY_NEIGHBORHOOD_PROBE`. The Atlas neighborhood is silent on shortcut techniques.

All 16 curated neighbors of #114 in the live workbench atlas (`Math/erdosatlas-workbench/erdosatlas.db`, 17 MB) were enumerated, resolution-status verified (8 solved, 4 disproven, 4 open — the original recon's 7/4/5 count was off by one because #114 was counting as a neighbor of itself), and resolution techniques surveyed via Lean files + scoring DB + 4 of 12 budgeted Perplexity sonar-low calls.

Top three non-brute-force transfer candidates ranked by structural-link strength to #114:

1. **#116** (area of sublevel set) — technique: subharmonic potential theory on `u = log|p|` with Cartan covering. Hypothesis: a coarea-formula bridge could convert area to length. Plausibility: **LOW**, because the required pointwise gradient control is exactly what the Fryntov-Nazarov Stokes / Cauchy / area-integral toolkit already supplies in the live #114 program. Transfer lands back in the existing toolbox — no shortcut.
2. **#1120** (dual-extremization). Plausibility: **SPECULATIVE**, no known bridge.
3. **#115** (Bernstein-Walsh). Plausibility: **LOW**, weaker structural link than #116.

Disproof counterexamples in the neighborhood (#1046, #1047, #1048) bound topology, not length. Littlewood machinery in the cross-corpus edges is orthogonal to arbitrary-monic arc length.

Net: the neighborhood does not surface a method that bypasses brute force. **The Atlas/Collider machinery did its job** — it surveyed the structural neighborhood and reported honestly that no shortcut exists, rather than fabricating one.

## Probe 3 — WS-01 cross-problem extension outcome

**Verdict:** `STRONG_PORTABILITY`, refined to "boundary-geometry-bound, not objective-bound."

Three new candidate problems tested:

- **#1041** (EHP58 path-length conjecture, paper-sibling of #114): YES on all three structural-fit features. Cold test: `WS01_APPLICABLE`, 12 of 12 candidate failure boxes closed at the conservative 2R bound, branch-point validation 1.00, median required factor after WS-01 rewrite ≈ 343.
- **#1043** (Pommerenke projection-bound, EHP58+Pommerenke 1959/1961 lineage): YES on Features 1 and 3, PARTIAL on Feature 2 (Pommerenke's 1961 negative resolution is upstream of WS-01's wall test). Cold test: `WS01_APPLICABLE`, 12 of 12 boxes closed at 2R, validation 1.00, median required factor ≈ 337.
- **#20** (sunflower / Erdős–Ko–Rado): NO structural fit on all three features. Discrete combinatorial parameter space; no zero curve, no gradient, no interval IVT. Cold test correctly skipped per design. **The structural-fit checker fires correctly as a negative control.**

WS-01 architecture is now confirmed `WS01_APPLICABLE` on four polynomial-level-set problems: #114, #1038, #1041, #1043. The architecture survives an objective-functional change within the same parameter-space + zero-curve family (component length, path length, sublevel measure, projected measure). It correctly fails the structural-fit gate on discrete-combinatorial problems.

**Insight from Probe 3 sharper than the plan's binary axis anticipated**: the architecture is *boundary-geometry-bound, not objective-bound*. Whether or not the conjecture under examination is solved (Pommerenke 1961 solved #1043 negatively, the others are open), WS-01's wall test depends only on the boundary geometry (real-analytic zero curve, bounded-away gradient, interval-IVT branch points) — not on what's being measured on that boundary.

## Verdict synthesis

Plan matrix row hit: **`BOUNDED_BY_CASE_ANALYSIS` with `N₀ > 10⁸` × no shortcut × full WS-01 portability = Intractability confirmed; methodology portability supported.**

| Probe | Verdict |
|---|---|
| 1 (K_sting) | KSTING_BOUNDED_BY_CASE_ANALYSIS at ≤100; N₀ ≈ 10²⁵–10²⁷ — **INTRACTABLE_CONFIRMED** |
| 2 (Atlas/Collider) | NO_SHORTCUT_FOUND; top candidate #116 lands back in existing toolbox — **INTRACTABILITY_CONFIRMED_BY_NEIGHBORHOOD_PROBE** |
| 3 (WS-01 portability) | STRONG_PORTABILITY across 4 polynomial-level-set problems; negative control fires correctly — **METHODOLOGY_PORTABILITY_SUPPORTED** |

**Decision:** Adopt POC framing. The middle range is intractable under the tools at hand. The methodology portability claim is empirically supported across four worked examples plus a correctly-firing negative control. Write the POC technical report. Retire the bridge program. The global-inequality whitepaper at `EHP114_GLOBAL_INEQUALITY_WHITEPAPER_2026-05-15.md` becomes a documented open route, not an active research direction.

## What this verdict does and does not mean

**Does mean:**

- The candidate-`N₀` derivation from Tao's argument is structurally complete (Track A4) and the relevant Stokes-theorem constant `K_sting` has now been bounded numerically. The intractability follows from honest constants, not from missing extraction work.
- The Atlas/Collider machinery surveyed the structural neighborhood and reported no shortcut.
- The WS-01 architecture generalizes across polynomial-level-set problems, with a documented portability boundary at the structural-fit gate.

**Does not mean:**

- That Erdős #114 cannot be solved. The bridge program is empirically dead under the current `K_sting` + `K_outside` framework; an entirely different proof architecture (the global-inequality whitepaper's coarea-stability route, or a structural reduction the Atlas neighborhood didn't surface) may still close the conjecture in the future.
- That the K_outside extraction would not be informative. It would. But the K_sting result already shows that bounding K_outside at any plausible value leaves N₀ well above brute-force reach.
- That the methodology paper is unpublishable. The POC has four worked applications, a correctly-firing negative control, and an empirically-mapped portability boundary. It is publishable as a methodology contribution.

## Next move

Write the POC technical report. Target ~10–15 pages at `~/Downloads/EHP114_INTRACTABILITY_POC_TECHNICAL_REPORT_2026-05-16.md`. Sections per the plan: framework; what was verified; where the pipeline stops; DeepMind formal-conjectures contributions; open routes documented but not pursued.

Modal Tier-2 budget remains parked. No new compute spend authorized by this verdict.

## Artifact references

- Probe 1: `Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-TAO-KSTING-EXTRACTION-20260515-01/` (three-file contract, SHA-256 `800b797df983a0598ac7a8d340e88f50b559b900ce0b6e70d15132836bbcc953`)
- Probe 2: `Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-ATLAS-COLLIDER-NEIGHBORHOOD-PROBE-20260515-01/` (three-file contract, SHA-256 `4e0592bf…8fa9ee93`)
- Probe 3: `Math/erdos-experiments/Erdos1041/`, `Erdos1043/`, `Erdos20/` per-candidate directories + `Math/erdos-experiments/CROSS_PROBLEM_WS01_GENERALIZATION_2026-05-15.md` aggregate
- Background plan: `~/.claude/plans/yes-wondrous-blum.md` § "Intractability test program — 2026-05-15 (POC framing)"

---

*Disclosure: portions of this analysis were developed with AI assistance. The K_sting numeric composition, the Atlas neighborhood enumeration, the WS-01 cross-problem cold tests, and the three-file contracts cited above have been independently verified by the author. AI tools are not authors. Any errors are the author's responsibility.*
