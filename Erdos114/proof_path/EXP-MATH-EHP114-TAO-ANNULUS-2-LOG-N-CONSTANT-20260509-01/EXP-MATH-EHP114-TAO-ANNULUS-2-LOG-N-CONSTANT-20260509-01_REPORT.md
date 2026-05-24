# Track A4 — annulus-2 log-n Constant Extraction

**Experiment ID:** EXP-MATH-EHP114-TAO-ANNULUS-2-LOG-N-CONSTANT-20260509-01
**Status:** ANNULUS_2_LOG_N_BOUNDED_BY_CASE_ANALYSIS
**Date:** 2026-05-09

---

## Claim ceiling

Internal extraction packet only. Not a proof of Erdos #114, not a numeric N0 candidate, and not a public claim. One step in the Tao-dependency chain.

## Tao lemma under examination

**annulus-2** (Tao 2512.12455v2, line 1781, "Bound in intermediate region"):

> l(dE_1 cap Ann(0, r_-, r_+)) <= 2n(r_+ - r_-) - c sum_{zeta notin D(0, 10 r_-)} |zeta| + O_{C_0}(log n / n * ||p||) for an absolute constant c>0 (independent of C_0).

The `c` here is the same absolute constant that A2 lower-bounded numerically (c_10 >= 0.006900495524140533). The target of A4 is the SEPARATE constant: the implied multiplier inside `O_{C_0}(log n / n * ||p||)`.

## Upstream inputs from prior tracks

- **A2:** c_10 >= 0.006900495524140533 (rigorous interval-arithmetic lower bound on the disk integral over D(0, 1/21)).
- **A3:** K_OR (origin-repulsion's implied constant) is BOUNDED_BY_CASE_ANALYSIS, extractable from poz/pocl/Taylor-remainder chain.

A3's recommendation to queue A4 first turned out correct: this is the first chain in which c_10 feeds in as a numeric input.

---

## Chain walk through annulus-2's proof (lines 1786–1827)

The `O_{C_0}(log n / n * ||p||)` remainder in annulus-2's statement decomposes into TWO additive contributions (and three structural feeds that absorb into them):

### Step 1 — sting-disk-contribution (line 1787)

> "By Lemma sting, the combined contribution of these disks D(zeta, 2 log^{1/2} n / n * r_+) to the lemniscate length is O_{C_0}(log n ||p||/n)."

**Verdict: BOUNDED_BY_CASE_ANALYSIS.**

Number of bad critical points (those whose ε-neighborhood intersects the annulus) is at most `||p||_1 / (r_-/2) = 2 ||p||_1 / r_- <= 2 C_0^2` by markov (eq. \eqref{markov}, line 242). Each disk has radius `rho = 2 log^{1/2} n / n * r_+`. Lemma sting (line 1442) gives per-disk bound `K_sting * n rho^2 / |z_0|` with `|z_0| asymp_{C_0} r_+`, yielding per-disk contribution `4 K_sting C_0^2 * (log n / n) * ||p||`. Total:

```
sting-disk-total <= 8 * K_sting * C_0^4 * (log n / n) * ||p||
```

`K_sting` is the absolute constant in Lemma sting's `<<`. Its proof at lines 1446–1466 chains Theorem `stokes` (the X1..X5 estimates), markov, log-integral, and triangle inequality. This is structurally identical to A3's K_OR chain — constructive, traceable, ~200–400 lines of real analysis to numerify, no Roth-style non-effective input.

### Step 2 — ute-2-correction (line 1802)

> "ell(dE_1 cap Omega) = int_{r_-}^{r_+} sum_{z in dE_1 cap Omega cap dD(0, r)} 1 + O_{C_0}(1/(n^2 delta(z)^2)) dr"

The `1/(n^2 delta^2)` correction is bounded by ute-2 (line 1807):

> int_{r_-}^{r_+} sum 1/delta(z)^2 dr <<_{C_0} n ||p|| log n

Tao establishes ute-2 (lines 1810–1818) by reducing to a per-critical-point bound `int 1/|z-zeta|^2 dr <<_{C_0} n log n / ||p||`, dyadic-decomposed and bounded again via Lemma sting.

**Verdict: BOUNDED_BY_CASE_ANALYSIS.**

The 1/n² prefactor on the delta-correction (from cosec-asym, line 585: a Taylor expansion `1/|sin| = 1 + O(dist²)` with absolute numeric coefficient ~1) cancels one of the n's, giving a contribution

```
ute-2-correction <= K_arc2 * K_ute2 * (log n / n) * ||p||
```

where `K_arc2` is the cosec-asym Taylor coefficient (numeric absolute, ~1) and `K_ute2` is a product of `(markov-count) * K_sting * (dyadic geometric series sum)`. Same status as Step 1: BOUNDED_BY_CASE_ANALYSIS, K_sting unfilled.

### Step 3 — rouche-step (line 1826)

> "But this follows from Rouche's theorem by repeating the rest of the proof of Proposition ets (with implied constants now depending on C_0 instead of eps)."

**Verdict: NUMERIC. Contribution to log-n constant: 0.**

Rouche's theorem itself is exact. This step contributes ONLY to the lower-bound count of zeroes (used in the gain term `-c sum |zeta|`), not to the log-n loss. So it does not enter the constant we are tracking.

### Step 4 — pform-error-feed (line 1793)

Proposition pform (line 889) gives `p(z) = -1 + (z p'(z) / n) (1 + O(||p||² / (n |z|² delta(z))))`. For z in Omega cap Ann(0, r_-, r_+), |z| asymp_{C_0} r_- asymp_{C_0} ||p||, so the relative error is O_{C_0}(1/n) under Omega's delta cutoff `delta >= log^{1/2} n / n` (line 1791).

**Verdict: BOUNDED_BY_CASE_ANALYSIS, but fully absorbed into K_ute2's chain.**

No new bottleneck. pform's error is dominated by the 1/(n delta) term from arcl-2 + cosec-asym (line 1799).

### Step 5 — stokes-X1-X5-feed (line 1448, internal to Lemma sting)

Theorem stokes is applied in Lemma sting's proof with Omega = D(z_0, r) and lambda = n / |z_0|. The five error integrals X1..X5 (defined at lines 1051–1058 region) are bounded as:

- `X_2 << n r² / |z_0|` (direct calculation, line 1451)
- `X_4 = 0`
- `X_5 << r`
- `X_1 << n r²` via pp-bound (line 1397, `|p'| << n` on E_4)
- `X_3` chains markov + log-integral (lines 1454–1459)

**Verdict: BOUNDED_BY_CASE_ANALYSIS.** All constituents are constructive Vinogradov bounds. No dependence on c_10. Same case-analysis structure as A3's K_OR.

---

## Verdict on the log-n constant

**`K_logn = (8 K_sting + K_arc2 * K_ute2) * C_0^4`**

with:
- `K_sting`: BOUNDED_BY_CASE_ANALYSIS (Stokes X1..X5 chain, absolute, extractable)
- `K_arc2`: numeric absolute (cosec-asym Taylor coefficient ~1)
- `K_ute2`: numeric * K_sting (dyadic geometric series + markov + sting)

Symbolically:

```
implied multiplier in O_{C_0}(log n / n * ||p||) <= K_logn_universal * C_0^4
```

where `K_logn_universal` is an absolute constant fully determined by Lemma sting's proof. The C_0^4 polynomial scaling comes from chasing `r_+ = C_0^2 ||p||` through the bad-disk count, the per-disk sting bound, and the ute-2 dyadic decomposition.

We do NOT pin C_0 here — C_0 is the free parameter that the next-stage trade-off (line 1868: `c ||p|| (gain) vs O(||p||/C_0) (loss from outside-again) + O_{C_0}(log n / n * ||p||) (loss from annulus-2)`) optimizes. Pinning C_0 prematurely throws away the tradeoff's degree of freedom.

**Status: ANNULUS_2_LOG_N_BOUNDED_BY_CASE_ANALYSIS.**

---

## Candidate N0 sketch

With c_10 numeric (`c >= 0.0069` from A2), the gain in annulus-2 is at least `0.0069 * sum_{zeta notin D(0, 10 r_-)} |zeta|`. Combining with inside-2's gain (which uses defect-psi-cor at C=10 too, so c >= c_10 carries through there) and outside-again's `O(||p||/C_0)` loss, Tao's overall a-ineq trade-off (line 1868) reduces to:

```
loss(C_0, n) = K_logn_universal * C_0^4 * (log n / n) + K_outside / C_0  <  c = 0.0069
```

Optimal C_0 from `d/dC_0 = 0`:

```
C_0_opt^5 ~ K_outside * n / (4 * K_logn_universal * log n)
=> C_0_opt ~ (n / log n)^{1/5} * (K_outside / (4 K_logn_universal))^{1/5}
```

Loss at optimum:

```
loss_opt ~ 5 * (K_logn_universal * K_outside^4)^{1/5} * (log n / n)^{1/5} / 4^{4/5}
```

Setting `loss_opt < 0.0069`:

```
n / log n  >  (5/4^{4/5})^5 * (K_logn_universal * K_outside^4) / 0.0069^5
```

Numerically `(5/4^{4/5})^5 ~ 5.96`, `1/0.0069^5 ~ 6.45e10`. So:

```
n / log n  >  ~3.85e11 * K_logn_universal * K_outside^4
```

This formula is **structurally complete** but **not numeric** — pending K_logn_universal (this session's bottleneck) and K_outside (a future session's bottleneck). With both extracted, a concrete N0 lands.

**Feasibility status:**

- Today, after A4: not feasible (two universal constants still symbolic).
- After A5 (K_sting extraction → K_logn_universal numeric): one constant remaining (K_outside).
- After A5 + A5-followup (K_outside extraction): fully numeric N0 candidate.

---

## Next dependency

**Queue A5: extract K_sting** (Lemma sting's absolute Vinogradov constant, Tao lines 1446–1466).

Subroutines required:
1. Absolute constants in X1..X5 of Theorem stokes (line 1448).
2. Cauchy-bound K_pp from `|p'| << n` on E_4 (line 1397).
3. Markov-bound numerics (line 242, exact, no extraction needed).

All constructive Vinogradov; no non-effective input. Estimated complexity: 300–400 lines of focused real analysis, no compute. Once K_sting is in hand, K_logn_universal reduces to a polynomial expression in K_sting and a handful of small absolute numerics (cosec-asym Taylor, dyadic geometric series), all extractable in the same session.

A5-followup (K_outside): independent, feeds outside-again's O(||p||/C_0). Together they close out the loss-side of Tao's a-ineq trade-off.

---

## Honest assessment

- **Did c_10 feed into A4?** Yes — directly, as the gain coefficient. A4 is the first place where A2's numeric input is consumed.
- **Did A4 reveal a brand-new bottleneck or a known one?** A KNOWN one. K_sting was already flagged by A1 as the symbolic core for inside-2. A4 reveals it ALSO feeds annulus-2's log-n constant. **A single dedicated extraction unlocks both.**
- **Is Track A still alive?** Yes. The parts inventory is now: c_10 (numeric), K_OR (BOUNDED, extractable), K_sting (BOUNDED, extractable, shared between inside-2 and annulus-2), K_outside (BOUNDED, extractable). Three focused real-analysis sessions, no compute, no novel theory.
- **Track A is one shared extraction (K_sting) away from a structurally complete N0 formula, and three extractions away from a fully numeric N0 candidate.** The shared-extraction observation is the structural-cartographic value of A4.

---

## Cross-references

- A1 dependency map: `EXP-MATH-EHP114-TAO-INSIDE-2-CONSTANT-EXTRACTION-20260509-01`
- A2 c_10 numeric: `EXP-MATH-EHP114-TAO-DEFECT-PSI-COR-C10-LOWER-BOUND-20260509-01`
- A3 origin-repulsion: `EXP-MATH-EHP114-TAO-ORIGIN-REPULSION-CONSTANT-20260509-01`
- Tao 2512.12455v2 source: `/tmp/tao114_src/lemniscate.tex`
- Key lines: 1781 (annulus-2 statement), 1786–1827 (proof), 1442 (sting), 889 (pform), 563–582 (arclength), 1640 (ets), 242 (markov)
