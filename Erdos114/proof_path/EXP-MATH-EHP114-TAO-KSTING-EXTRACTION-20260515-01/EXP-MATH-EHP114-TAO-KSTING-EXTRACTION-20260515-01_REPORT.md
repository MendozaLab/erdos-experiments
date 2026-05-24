# EXP-MATH-EHP114-TAO-KSTING-EXTRACTION-20260515-01

**Track A5 — K_sting extraction from Tao 2512.12455v2 Lemma `sting`**

**Date:** 2026-05-15 (POC framing; experiment-ID date preserved from plan)
**Status:** `KSTING_BOUNDED_BY_CASE_ANALYSIS`
**Claim ceiling:** Internal extraction packet only.

## Plain-language summary

K_sting was the suspected universal-Stokes-theorem bottleneck shared between Tao's `inside-2` and `annulus-2` regions (per Track A4). Walking the proof at lines 1446-1466 with the TeX source on disk shows the bottleneck is structural-only: every absolute constant on the right-hand side of Lemma sting traces to a named Tao label (`markov` line 242, `riesz-bound` line 739, `multip` line 743, `pp-bound` line 1397, `cosec-asym` line 585) verified directly in the source. None depends on c_10, K_OR, or any non-effective input.

Honest upper bound: `K_sting <= 100`. Tight upper bound (if the Cauchy constant K_pp in pp-bound is taken at 1): `K_sting <= 30`. Either qualifies as `BOUNDED_BY_CASE_ANALYSIS` for the candidate-N₀ formula.

The candidate-N₀ formula `n / log n > 3.85e11 · K_sting · K_outside⁴` then yields, at K_outside ≤ 1000 (the brief's generous placeholder):

| K_sting | RHS | Candidate N₀ |
|---|---|---|
| 1 (floor, impossibly tight) | 3.85e23 | ≈ 2.2 × 10²⁵ |
| 30 (tight) | 1.16e25 | ≈ 7.1 × 10²⁶ |
| 100 (honest UB) | 3.85e25 | ≈ 2.4 × 10²⁷ |

All scenarios land **17+ orders of magnitude above** the 10⁸ middle-range bar and **13+ orders above** the 10¹² INTRACTABLE_CONFIRMED bar.

**Verdict:** `INTRACTABLE_CONFIRMED`. K_sting was not the bottleneck. The c_10⁻⁵ ≈ 6.45 × 10¹⁰ amplification (from Tao's a-ineq trade-off at line 1868, with c_10 ≥ 0.0069 from Track A2) dominates; the K_outside⁴ = 10¹² multiplier dominates further. Refining K_sting numerically does not move the threshold.

## Step-by-step extraction

Tao's proof of Lemma `sting` (lines 1446-1466) applies Theorem `stokes` (lines 357-370) with `Ω = D(z_0, r)` and `λ = n/|z_0|`, yielding

```
ℓ(∂E_1 ∩ D(z_0, r)) ≤ Ψ(D(z_0, r)) + O(X_1 + X_2 + X_3 + X_4 + X_5)
```

with the X_i defined at lines 364-368 (X1-def through X5-def). Direct calculation produces:

- **X_4 = 0** exactly (line 1451). Tao's λ = n/|z_0| is radially constant; |λ'|/|λ| = 0 in D(z_0, r) since z_0 ≠ 0 and r < |z_0|/2 keeps the disk away from origin. *K_X4 = 0.*

- **X_2 ≪ nr²/|z_0|** (line 1451). By definition X_2 = ∫_{E_2 ∩ Ω, |ψ| ≤ λ} |ψ| dA ≤ λ · area(Ω) = (n/|z_0|) · πr². Quote: *"X_2 ≪ n r² / |z_0|"* (line 1451). *K_X2 ≤ π.*

- **X_5 = ℓ(∂Ω) = 2πr** exactly (line 1451, line 368). *K_X5 = 2π.*

- **X_1 ≪ nr²** via `pp-bound` (line 1397/1452). Quote: *"|p'(z)| ≪ n for all z ∈ E_4"* (line 1398). X_1 = ∫_{E_2 ∩ Ω} |p'| dA ≤ K_pp · n · πr². K_pp is the Cauchy-bound absolute constant in pp-bound; ≤ 10 generously, ≤ 1 if the upstream Proposition pvar chain is tightened. Using r < |z_0|/2 to land back in the n r²/|z_0| family adds a factor 2. *K_X1 ≤ 20π (honest) or 2π (tight).*

- **X_3 ≪ nr²/|z_0| + r + ‖p‖_1 log n / n** via `markov` + `multip` + log-integral (lines 1454-1459). Tao decomposes X_3 into three pieces: a minor contribution near critical points bounded by O(n · n · n⁻²r) = O(r) using `multip`'s factor 2π; a main contribution `|z_0|/n · n · r²/|z_0|²` for critical points in D(0, |z_0|/4) (which the improvement `r²/|z_0|²` is available for); and a markov-controlled contribution `|z_0|/n · ‖p‖_1/|z_0| · log n` for critical points outside D(0, |z_0|/4). Quote: *"Each integral is O(log n), but with the improved bound of O(r²/|z_0|²) if ζ ∈ D(0, |z_0|/4)"* (line 1458). *K_X3 ≤ 2π.*

- **Ψ(D(z_0, r)) ≤ 2(nr²/|z_0| + ‖p‖_1 r/|z_0|)** via triangle + `multip` + markov (lines 1460-1465). Tao quotes: *"The integral here is O(r) for all ζ (by `riesz-bound`), but the bound can be improved to r²/|z_0| if ζ ∈ D(0, |z_0|/4)"* (line 1464). Combined with Ψ-def's 1/π prefactor (line 345), the riesz 2π collapses to 2. *K_Ψ = 2 exactly.*

Sum: `K_sting ≤ K_Ψ + K_X1 + K_X2 + K_X3 + K_X5 = 2 + 20π + π + 2π + 2π ≈ 80.5`. Round up to **100** as honest UB with buffer for the cosec-asym (line 585) leading-order Taylor coefficient (= 1 exactly) and the dyadic geometric series sum 1 + 1/4 + 1/16 + … = 4/3 in the X_3 minor-contribution step. Tight UB: **30** if K_pp = 1.

## N₀ formula evaluation

Using A4's structural formula `n / log n > 3.85e11 · K_sting · K_outside⁴` (which encodes the c_10⁻⁵ amplification, the K_logn_universal · C_0⁴ scaling, and Tao's trade-off optimization at C_0⁵ ~ K_outside · n/(4·K_logn·log n)):

```
At K_outside = 1000, K_sting = 100:
RHS = 3.85e11 × 100 × 10^12 = 3.85e25
Solve n/log n = 3.85e25 by fixed-point: n ≈ 2.43e27
```

Threshold matrix from brief:

- `N₀ ≤ 10⁴` → INTRACTABILITY_REFUTED. **Not approached** (off by 23 orders).
- `N₀ ≤ 10⁸` → MIDDLE_RANGE_PLAUSIBLY_REACHABLE. **Not approached** (off by 19 orders).
- `N₀ ≥ 10¹²` → INTRACTABLE_CONFIRMED. **Vastly exceeded** (off by 15 orders in the safe direction).

## Honest scope

This packet extracts K_sting symbolically. It does NOT:

- Prove Erdős #114.
- Upgrade any public claim above `n ≤ 14` finite certificate (Zenodo v3.1.0).
- Extract K_outside (Track A6, queued).
- Refute Tao's effectivization remark at line 196 (extraction is demonstrated; the threshold is still astronomical).

It DOES:

- Confirm that the symbolic-only verdict on K_sting can be lifted to a numeric bound by direct case analysis, validating Tao's line-196 claim that all implied constants are effectively computable.
- Localize the structural intractability source: c_10⁻⁵ amplification × K_outside⁴ blow-up. K_sting is a sideshow.
- Close one of the three remaining symbolic constants in Track A's chain (c_10 numeric A2, K_OR bounded A3, K_logn_universal decomposed A4, K_sting now bounded A5). One symbolic constant remains: K_outside.

## Asymmetry note

The candidate N₀ scales as `(1/c_10)^5` from Tao's optimization. With c_10 ≥ 0.0069 (A2 rigorous interval bound), 1/c_10⁵ ≈ 6.45 × 10¹⁰, which is the dominant universal coefficient inside the 3.85 × 10¹¹. K_sting bounded does not shift this. Refining c_10 to 0.069 (10× larger, well beyond any plausible improvement on a rigorous interval-arithmetic disk integral) would only drop N₀ by 10⁵, leaving N₀ ≈ 10²² — still INTRACTABLE_CONFIRMED.

## Next dependency

Track A6: K_outside extraction from Lemma `outside-again` (pvar + pform + out chain, Tao lines ~1832-1864). Free local symbolic work. Expected verdict: BOUNDED_BY_CASE_ANALYSIS at K_outside ~ 30-300. Even at K_outside = 30, candidate N₀ drops only to ≈ 10²¹ — INTRACTABLE_CONFIRMED holds.

## Cross-references

- A1: `EXP-MATH-EHP114-TAO-INSIDE-2-CONSTANT-EXTRACTION-20260509-01` (7 symbolic-only constants mapped, K_sting flagged as the shared bottleneck for inside-2)
- A2: `EXP-MATH-EHP114-TAO-DEFECT-PSI-COR-C10-LOWER-BOUND-20260509-01` (c_10 ≥ 0.006900495524140533 rigorous lower bound)
- A3: `EXP-MATH-EHP114-TAO-ORIGIN-REPULSION-CONSTANT-20260509-01` (K_OR bounded by Taylor-remainder case analysis)
- A4: `EXP-MATH-EHP114-TAO-ANNULUS-2-LOG-N-CONSTANT-20260509-01` (K_logn = K_logn_universal · C_0⁴ decomposition; K_sting identified as shared bottleneck)
- Tao 2512.12455v2 final-section lines 1632-1873 (Tao endorses extraction at line 196)
