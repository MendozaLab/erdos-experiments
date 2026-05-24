# EXP-MATH-EHP114-TAO-INSIDE-2-CONSTANT-EXTRACTION-20260509-01

## Scope

Track A1 of the Erdős #114 (EHP) bridge program. Target: extract numeric values, or honest verdicts, for the existential constants invoked in the final-section proof of Tao's "The maximal length of the Erdős–Herzog–Piranian lemniscate in high degree" (arXiv:2512.12455v2, Dec 2025). This is internal extraction infrastructure. Not a proof of #114, not an N0 candidate, not a public claim.

## Fetch

Paper retrieved successfully via two channels: `WebFetch` on `https://arxiv.org/abs/2512.12455` for metadata confirmation, then `curl` on `https://arxiv.org/e-print/2512.12455v2` for the full TeX source (`lemniscate.tex`, 1965 lines, 142 KB). All five requested labels (`inside-2`, `annulus-2`, `outside-again`, `defect-psi-cor`, `origin-repulsion`) and the `pp0` condition resolved to exact line numbers. The line-number guesses supplied in the task brief were close: `inside-2` at line 1739, `pp0` at 1744, `annulus-2` at 1781, `outside-again` at 1831, `defect-psi-cor` at 1025, `origin-repulsion` at 869.

## inside-2 (Inner bound, Proposition at line 1739)

Tao states two inequalities, both governed by the same symbolic absolute constant `c > 0` independent of `C_0`. The first (`in2-first`, line 1741) bounds inner-region length by `2n r_- − c · Disp{ζ in D(0, 10 r_-)} + O(‖p‖/C_0)`. The second (`in2-second`, line 1749), gated by condition `pp0`, gains an additional `−c‖p‖` term.

Verdict: **SYMBOLIC_ONLY** (both inequalities). The proof at lines 1753–1777 chains a Stokes-theorem application with five Vinogradov error terms `X_1...X_5`, each bounded by `«` constants, and routes the displacement gain through `defect-psi-cor` at line 1765 with `C = 10`. The `c` here IS the `c_C` of `defect-psi-cor` evaluated at `C = 10`. It could in principle be extracted, but Tao does not extract it; the chain is opaque at every link.

## annulus-2 (Intermediate region, Proposition at line 1781)

Tao bounds intermediate length by `2n(r_+ − r_-) − c · Σ_{ζ ∉ D(0,10r_-)}|ζ| + O_{C_0}(log n / n · ‖p‖)`. The `c` is again declared absolute and independent of `C_0`.

Verdict: **SYMBOLIC_ONLY**. The proof (lines 1786–1827) decomposes via the arclength formula `arcl-2`, uses Proposition `pform` for fine pointwise behavior in the annulus, and reduces to Rouché's theorem applied to a polynomial `z^n p(z) p̃(r²/z) − z^n`. The trade-off here is between the absolute `c` (gain) and the `O_{C_0}(log n / n · ‖p‖)` remainder (loss). The latter's implicit `C_0`-dependent constant is not quantified anywhere in the proof. The user's note that this lemma "invokes Corollary `defect-psi-cor` (line ~1767), whose integral is C-dependent" is correct in spirit: the `c` carried over from inside-2's defect-psi-cor application at line 1765–1767 is reused here under the same absolute label. The `C_0`-dependence does not enter the gain term but instead lives in the additive log-n term.

## outside-again (Outer bound, Proposition at line 1831)

Tao bounds outer length by `ℓ(∂E_1(p_0)) − 2n r_+ + O(‖p‖/C_0)`. There is no named constant `c > 0` in the statement; the result is stated entirely in terms of remainders.

Verdict: **SYMBOLIC_ONLY**. Proof (lines 1835–1865) runs through Proposition `pform` (line 889, with a `O(‖p‖²/(n|z|²))` error term proven by a fundamental-theorem-of-calculus argument plus a substitution `t = 1 − s/n` controlled by `exp(−s/2)`), then through Proposition `pvar` and Lemma `out`, with a final residual integral evaluating to `O(n / C_0² ‖p‖)` at line 1864. There is one numerically-explicit exponent in this chain (`C_0²` versus `C_0`), but every Vinogradov `«` carries an unextracted absolute constant.

## defect-psi-cor (Corollary at line 1025) — the bottleneck

This is the source of the `c` that propagates into both `inside-2` and `annulus-2`. Statement: `Ψ(E) ≤ 2(n−1) r̄[E] − c_C · Disp{ζ : ζ in D(z_0, Cr)}` with `c_C > 0` "depending only on `C`."

Verdict: **BOUNDED_BY_CASE_ANALYSIS**. This is the only constant in the chain that is genuinely extractable in principle. The proof at lines 1032–1038 reduces (via Lemma `tridef` line 1012) to a CONCRETE integral over `D(0, 1/(2C+1))` of the function `f(z) = 1/|z| + 1/|z−1| − |1/z + 1/(z−1)|`. This function is universal, geometrically explicit, non-negative, absolutely integrable, and non-vanishing off the real axis. For the value `C = 10` (which is what `inside-2` uses), the disk is `D(0, 1/21)`, and a positive rational lower bound on `c_10` could be obtained by a focused computer-algebra evaluation or a piecewise-constant inner approximation on a sufficiently fine grid. A few hundred lines of focused real-analysis + interval-arithmetic work would settle it. Tao does not perform this calculation in the paper.

## origin-repulsion (Lemma at line 869)

Statement: `‖p‖_0 « ‖p‖_1 + n · dist(0, ∂E_1(p))`, hence `‖p‖ « ‖p‖_1 + n · dist(0, ∂E_1(p))`.

Verdict: **SYMBOLIC_ONLY**. Proof is a clean two-case argument (lines 877–881) that chains earlier bounds `poz`, `npz`, `pocl`, each carrying its own implicit constant. The argument is constructive — no non-effective Roth-style input — so the `«` could be quantified after upstream extraction, but Tao does not do so.

## pp0 (Condition at line 1744)

Statement: `‖p‖ ≥ C ‖p‖_1` for "a large absolute constant `C` (independent of `C_0, c`)."

Verdict: **SYMBOLIC_ONLY**. The threshold `C` is chosen at line 1770 so that `origin-repulsion` forces `∂E_1` to avoid `D(0, ‖p‖/4n)`; thus `C ≈ 4 · (origin-repulsion's implied constant) + buffer`. Inherits opacity from the upstream lemma.

## Verdict

Six of seven named constants are **SYMBOLIC_ONLY**. Exactly one — `defect-psi-cor` `c_C` at `C = 10` — is **BOUNDED_BY_CASE_ANALYSIS**. Zero are **NUMERIC**. Zero are **NON_EFFECTIVE** (the chain is structurally constructive end-to-end; no Roth-style step is invoked, no Cauchy-Schwarz over an unbounded set, etc.).

Track A1 alone does not produce a numeric `N_0`. The trade-off at line 1868, where `c‖p‖` must dominate `O(‖p‖/C_0)`, requires a numeric `c` to fix `C_0`; the `n`-threshold then propagates from `r_+ = C_0² ‖p‖ = o(1)` (line 1730), but without numeric `c` no numeric chain closes.

## Path forward

Track A is **not dead, but heavily blocked**. The single live extraction target is **Track A2**: compute a positive rational lower bound on `c_10`, the universal disk integral over `D(0, 1/21)` of `f(z) = 1/|z| + 1/|z−1| − |1/z + 1/(z−1)|`. This is a one-shot, low-risk computation on a known geometric object. If `c_10` lands cleanly, the secondary blockers — `origin-repulsion`'s implied constant (Track A3) and the `O_{C_0}(log n / n · ‖p‖)` constant in annulus-2 (Track A4) — both become tractable. If any one of A2/A3/A4 falls through, Track A is dead and Track B (direct small-n verification at `n ≤ 14` plus per-`n` upper-bound certificates extending toward Tao's high-`n` regime) carries the load.

**Recommendation:** Queue A2 immediately. Run Track B in parallel; do not gate B on A2.
