# Track A3 — Origin-Repulsion Constant Extraction

**Experiment ID:** EXP-MATH-EHP114-TAO-ORIGIN-REPULSION-CONSTANT-20260509-01
**Date:** 2026-05-09
**Lemma target:** Tao 2512.12455v2, Lemma `origin-repulsion`, line 869
**Upstream artifact consumed:** EXP-MATH-EHP114-TAO-DEFECT-PSI-COR-C10-LOWER-BOUND-20260509-01 (`c_10 ≥ 0.006900495524140533`)
**Claim ceiling:** Internal extraction packet only. Not a proof of Erdős #114, not an N₀ candidate, not a public claim.

## Statement of the lemma

> `‖p‖_0 ≪ ‖p‖_1 + n · dist(0, ∂E_1(p))`, hence `‖p‖ ≪ ‖p‖_1 + n · dist(0, ∂E_1(p))`.

The `≪` hides a single positive absolute constant. Call it `K_OR`. The question A3 asks: with Track A2's numeric `c_10` in hand, can we now write `K_OR` as a positive rational?

## Headline answer

**No, but the reason is more useful than A1 suggested.** `c_10` does not feed into origin-repulsion at all — they live in disjoint chains. A1's classification was directionally right (origin-repulsion is `BOUNDED_BY_CASE_ANALYSIS` once one does the upstream work) but identified the wrong upstream. The actual upstream is a Taylor remainder coefficient inside Proposition `pvar`, not `defect-psi-cor`'s disk integral.

The good news: this puts `K_OR` firmly in `BOUNDED_BY_CASE_ANALYSIS` with a concrete extraction recipe. The bad news: this session does not perform that extraction; another focused pass is required.

## The chain, line by line

Tao's proof at lines 877–881 splits into two cases on `|z|` where `z` is the closest point of `∂E_1(p)` to the origin.

### npz — line 253

> `|p(z) - p(0)| ≥ n^{-n} ‖p‖_0^n` for any `z ∈ ∂E_1(p)`.

This is **exact**, not Vinogradov. Derivation: `p(0)` is normalized to a non-positive real, so `dist(p(0), ∂D(0,1)) = |1 + p(0)| = n^{-n} ‖p‖_0^n` by definition of `‖p‖_0` (lines 230–231 and 250–251). On `∂E_1(p)`, `|p(z)| = 1`. Triangle inequality finishes.

**Verdict: NUMERIC, constant = 1.** No opacity.

### poz — line 803

> `p(z) = p(0) + O(|z|^n · e^{O(min(‖p‖_1/|z|, ‖p‖_1²/|z|²))})`

Two implied constants stack. The outer envelope is a triangle-inequality artifact after Weierstrass-style factorization (Tao line 822):
`|p'(z)| = n|z|^{n-1} · ∏_ζ |(1 - ζ/z) · e^{ζ/z}|`.
The inner `O(min(...))` comes from the Taylor remainder of `(1 - w) · e^w` at `w = ζ/z` (line 826). This is `(1-w)e^w = 1 - w²/2 + O(w³)`, so for `|w| ≤ 1/2`, `|(1-w)e^w| ≤ e^{2|w|²}` is a clean explicit bound; the constant 2 (or smaller, with sharper analysis) is numeric.

The proof of `poz` itself (lines 838–850) introduces "a sufficiently large constant `C`" and then "another large constant `C'`" — each chosen so that the radial derivative of `|p(z) - p(0)| - C|z|^n e^{C‖p‖_1/|z|}` (and its second variant) stays negative. Closing the inequality requires `C·e^{C·u} - K·e^{K·u} - C²·u·e^{C·u} > 0` where `K` is the implied constant from `p'-bound` and `u = ‖p‖_1/|z|`. This is a one-variable inequality in `C`; solving for the smallest `C` is mechanical real analysis but unfilled in the paper.

**Verdict: BOUNDED_BY_CASE_ANALYSIS.** Roughly 50 lines of work to numericize.

### pocl — line 813

> `p(z) = p(0) + O(e^{O(n)} · n^{-n} · ‖p‖_1^n)` for `z = O(‖p‖_1/n)`.

Tao's proof at line 837: "For `|z| ≍ ‖p‖_1/n` the claim `ppcl` comes from `p'-bound`; the case `|z| ≪ ‖p‖_1/n` then follows from the maximum principle. To obtain `pocl`, we integrate `ppcl`."

The inner `O(n)` in the exponent is `p'-bound`'s `e^{O(min(‖p‖_1/|z|, ‖p‖_1²/|z|²))}` evaluated at `|z| ≍ ‖p‖_1/n`, where the min equals `‖p‖_1/|z| ≍ n`. So the `O(n)` exponent constant is the same Taylor coefficient as in `poz`. The outer envelope is a radial integration, lifting that constant by a numeric factor (≤ 2).

**Verdict: BOUNDED_BY_CASE_ANALYSIS.** Same machinery as `poz`, no new opacity.

## Composing the cases into K_OR

**Case 1** (`|z| ≥ ‖p‖_1/n`): from `npz` ≤ `poz`,
`n^{-n} ‖p‖_0^n ≤ K_poz · |z|^n · e^{K_poz_exp · n}`
(since `min` reduces to `‖p‖_1/|z| ≤ n`). Take n-th roots:
`‖p‖_0 ≤ n · |z| · (K_poz)^{1/n} · e^{K_poz_exp}`,
giving `K_OR_case1 = e^{K_poz_exp} · max(1, K_poz)`.

**Case 2** (`|z| < ‖p‖_1/n`): from `npz` ≤ `pocl`,
`n^{-n} ‖p‖_0^n ≤ K_pocl · e^{K_pocl_exp · n} · n^{-n} · ‖p‖_1^n`.
Take n-th roots:
`‖p‖_0 ≤ ‖p‖_1 · (K_pocl)^{1/n} · e^{K_pocl_exp}`,
giving `K_OR_case2 = e^{K_pocl_exp} · max(1, K_pocl)`.

**Final:** `K_OR = max(K_OR_case1, K_OR_case2)`.

Each of `K_poz`, `K_poz_exp`, `K_pocl`, `K_pocl_exp` is in `BOUNDED_BY_CASE_ANALYSIS` via Taylor remainder of `(1 - w)e^w` plus the derivative-closure inequality in `poz`'s proof. With ~100–200 lines of focused real analysis, a numeric `K_OR` could be produced.

## Where c_10 actually feeds in

The Track A2 result `c_10 ≥ 0.006900495524140533` is the value of `c_C` at `C = 10` from `defect-psi-cor` (line 1025). Tao cites `defect-psi-cor` at lines 1234, 1665, and 1765 — none of which are in origin-repulsion's chain. `defect-psi-cor` is upstream of `inside-2` (line 1739) and `annulus-2` (line 1781); origin-repulsion (line 869) is upstream of those.

The two chains share no symbols. `c_10` does not unblock A3.

This is not a discovery against A1; it is a tightening of A1's verdict. A1 said "possibly upgradable to `BOUNDED_BY_CASE_ANALYSIS` after upstream extraction." The line-by-line read for A3 confirms `BOUNDED_BY_CASE_ANALYSIS` and identifies the actual upstream as the Taylor coefficient in `p'-bound`, with no involvement of `c_10`.

## Status of the broader extraction program

After Tracks A1, A2, A3:

- `c_10` — numeric (A2)
- `K_OR` — `BOUNDED_BY_CASE_ANALYSIS` with concrete extraction recipe (A3)
- `pp0`'s `C` — was `SYMBOLIC_ONLY` per A1; **now `BOUNDED_BY_CASE_ANALYSIS`** since per Tao line 1770, `pp0_C ≥ 4 · K_OR + buffer`; opacity inherited from `K_OR`'s `BOUNDED_BY_CASE_ANALYSIS` status, not its own
- `inside-2` absolute `c` — still `SYMBOLIC_ONLY` per A1; depends on `c_10` (now numeric) plus Stokes-theorem `X1..X5` estimates (still opaque)
- `annulus-2` absolute `c` — still `SYMBOLIC_ONLY` per A1; depends on `inside-2` plus Lemma `sting` plus Proposition `pform`
- `outside-again` `O(‖p‖/C_0)` envelope — still `SYMBOLIC_ONLY` per A1

So one constant became numeric (A2's `c_10`), one moved from `SYMBOLIC_ONLY` to `BOUNDED_BY_CASE_ANALYSIS` with recipe (A3's `K_OR`), and one moved indirectly (`pp0_C` via inheritance from `K_OR`). The remaining bottlenecks are the Stokes-theorem `X1..X5` estimates inside `inside-2`, which both Tracks A2 and A3 leave untouched.

## Honest closing verdict

`origin-repulsion`'s implied constant **is not numeric after A3**. The verdict moves from `SYMBOLIC_ONLY` (A1) to `BOUNDED_BY_CASE_ANALYSIS` (A3) — meaningful progress in the parts inventory, no progress in the numeric ledger. Track A is alive. The next two parallel queues are:

1. **A3-followup**: actually run the `K_OR` extraction (Taylor remainder + derivative-closure inequality, ~100–200 lines of real analysis, no compute).
2. **A4**: annulus-2's `log n / n` constant, which is the first place A2's `c_10` numeric input actually does work.

Recommend: queue A4 first since it consumes A2's output directly; defer A3-followup to a separate focused session. Track A's clean N₀ remains gated by the `inside-2` Stokes-theorem `X1..X5` estimates, which neither this session nor A2 touched.
