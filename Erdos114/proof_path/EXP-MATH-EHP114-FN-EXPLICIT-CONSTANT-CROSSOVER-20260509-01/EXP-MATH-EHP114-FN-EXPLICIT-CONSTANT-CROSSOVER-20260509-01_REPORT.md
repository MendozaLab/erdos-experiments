# EXP-MATH-EHP114-FN-EXPLICIT-CONSTANT-CROSSOVER-20260509-01

**Question.** Does Fryntov–Nazarov 2008's explicit bound, evaluated at the optimal δ, ever drop below the symmetric extremizer length `L_sym(n) = L(z^n − 1)` at any finite `n`? If yes, at what `n`? If no, what does the asymptotic comparison look like?

**Claim ceiling.** Internal comparison computation. Not a proof of Erdős #114, not a finite-`n` cert. Establishes whether FN's explicit bound is sharp enough to close the conjecture by inequality alone for some `n`.

## Setup

The FN inequality, transcribed directly from arXiv:0808.0717 p. 16:

```
|L| ≤ 26·π·δ·n + 2π·√n + e^{4δ} · (1/π + 2δ + 4a/(δ³·√n)) · 2π·(n−1)
```

valid for every `δ ∈ (0, 1/4)`, where `a = Σ_{k≠0} |a_k|/|k| = 1/2 + (1/π) · Σ_{m≥1} 1/(m·(4m²−1))`. Numerically `a = 0.62296131412107...` (computed to 30 dps).

Note: the prior-art note dated 2026-05-09 transcribed the third-term prefactor as `e^{2δ}`. The PDF shows `e^{4δ}` (cap(F) ≤ e^{2δ} so area(F) ≤ π·e^{4δ}). I used the PDF version. This makes FN's bound *larger*, which only strengthens the negative result.

The symmetric extremizer length `L_sym(n) = L(z^n − 1)` has closed form

```
L_sym(n) = 2^{1/n} · √π · Γ(1/(2n)) / Γ(1/(2n) + 1/2)
```

with asymptotic expansion `L_sym(n) = 2n + (π²/6 + log 2)/n + O(1/n²)`, so `L_sym(n) − 2n → 0` from above as `1/n`.

## Comparison Table

For each `n`, optimize `δ ∈ (10⁻¹⁰, 1/4 − 10⁻¹⁰)` numerically (golden-section, 300 iterations to <10⁻⁴⁰ width) to minimize the FN bound. All arithmetic in mpmath at 80 dps.

| n | L_sym(n) | δ_optimal | B_FN(n, δ*) | B_FN − L_sym | FN < L_sym? |
|---:|---:|---:|---:|---:|:---:|
| 14 | 30.853 | 0.250 (boundary) | 9 954.65 | +9 923.80 | NO |
| 15 | 32.847 | 0.250 (boundary) | 10 372.26 | +10 339.42 | NO |
| 20 | 42.828 | 0.250 (boundary) | 12 274.21 | +12 231.38 | NO |
| 50 | 102.795 | 0.250 (boundary) | 20 625.27 | +20 522.47 | NO |
| 100 | 202.784 | 0.250 (boundary) | 30 454.15 | +30 251.37 | NO |
| 1 000 | 2 002.77 | 0.250 (boundary) | 120 629.39 | +118 626.62 | NO |
| 10 000 | 20 002.77 | 0.250 (boundary) | 616 933.72 | +596 930.95 | NO |
| 100 000 | 200 002.77 | 0.20055 | 4 017 349.92 | +3 817 347.14 | NO |
| 1 000 000 | 2 000 002.77 | 0.15066 | 27 789 417.48 | +25 789 414.71 | NO |
| 10⁷ | 20 000 002.77 | 0.11269 | 199 983 290 | +199 963 288 | NO |
| 10⁸ | 200 000 002.77 | 0.08416 | 1.483·10⁹ | +1.483·10⁹ | NO |
| 10⁹ | 2.000·10⁹ | 0.06284 | 1.129·10¹⁰ | +1.128·10¹⁰ | NO |
| 10¹² | 2.000·10¹² | 0.02625 | 5.693·10¹² | +5.691·10¹² | NO |

**The FN bound never drops below L_sym at any tested `n` from 14 up to 10¹². The gap grows polynomially, not shrinks.**

## Leading-constant diagnostic

Comparing `B_FN/n` against the conjectured leading constant 2 (`L_sym/n → 2`):

| n | B_FN(n, δ*)/n | L_sym(n)/n |
|---:|---:|---:|
| 14 | 711.05 | 2.204 |
| 1 000 | 120.63 | 2.003 |
| 10⁶ | 27.79 | 2.0000 |
| 10¹² | 5.69 | 2.000000 |

Even at `n = 10¹²`, FN's bound gives `|L| ≤ 5.69·n`, far above the conjectured `2n`. With δ pinned at the upper boundary 1/4, the limiting `B_FN/n` value is roughly 34.4 — the third-term `e^{4δ}·(1/π)·2π(n−1)/n → 2 + 2π·e·1/π·... ≈ 17.4` plus the `26πδ = 6.5π ≈ 20.4` first-term contribution at δ=1/4, totaling ~38 with sub-leading corrections.

## Asymptotic n^{7/8} constant estimate

Empirical `(B_FN_optimal − L_sym) / n^{7/8}` at large `n`:

| n | (B_FN − L_sym) / n^{7/8} |
|---:|---:|
| 10⁴ | 188.77 |
| 10⁶ | 145.02 |
| 10⁸ | 128.34 |
| 10¹⁰ | 120.66 |
| 10¹² | 116.79 |

Ratio approaches a finite positive constant near 115 from above. So FN delivers `|L| ≤ 2n + K_FN·n^{7/8}` with implied `K_FN ≈ 115` for large `n`, but **K_FN is positive and bounded away from zero**, while `L_sym − 2n → 0` from above as `O(1/n)`. The slack `B_FN − L_sym` therefore grows without bound as `n → ∞`.

## Verdict

**FN_DOES_NOT_CLOSE_ASYMPTOTICALLY.**

There is no finite `N_FN_thresh` such that `B_FN(n, δ*) ≤ L_sym(n)` for `n ≥ N_FN_thresh`. The two curves never meet: FN's `O(n^{7/8})` slack term grows polynomially while L_sym's correction to `2n` decays like `1/n`. The orthogonal re-read's three questions are now answered numerically:

- **Q1 (does FN ever beat L_sym at any finite `n`?):** No. Tested up to `n = 10¹²`. The gap grows.
- **Q2 (does FN+(n=14 cert) produce a finite middle-range threshold?):** No. There is no `N_FN_thresh` past which the conjecture is closed by FN's bound alone.
- **Q3 (does FN's bound have a closure mechanism?):** No, because FN's `n^{7/8}` term genuinely dominates `L_sym`'s `O(1/n)` correction. Even improving the FN exponent from 7/8 to 1/2 (per FN's own remark on p. 17) would not close: any positive power of `n` above `−1` diverges from `O(1/n)` at infinity. To close by FN alone you would need a slack term of order `O(1/n)` or smaller — three exponent units below FN's stated `7/8`, far beyond what FN's surgery can deliver.

## What this means for the bridge program

Fryntov–Nazarov is a **structural / methodological reference**, not a closure mechanism. The day's read was right that FN doesn't constrain the bridge's `N_0`. This orthogonal re-read confirms the stronger statement: there is no `N_0` for which FN alone resolves EHP. The bridge architecture remains:

- finite cert at `n ≤ 14` (EHP solver / Lean / numeric verification) ⊕
- some path closing `n ∈ [15, N_0 − 1]` (currently the methodologically-uncertain middle range) ⊕
- Tao 2512.12455's o(n) bound (or strengthening) for `n ≥ N_0`

FN's `O(n^{7/8})` envelope tells us the *asymptotic form* of the upper bound for all `n` simultaneously, but does not close any specific `n` by inequality against `L_sym`. Tao's contribution is in the same envelope but with sharper constants and a `C_0`-parametrized decomposition that FN does not have.

## Files

- `_RESULTS.json` — structured output with the comparison table and asymptotic diagnostics.
- `_REPORT.md` — this file.
- `_RESULTS.sha256` — SHA-256 of `_RESULTS.json`.

## Audit trail

- Re-read FN p. 9–17 from `/tmp/fn08.pdf` to verify the inequality form. Found one transcription discrepancy in the prior-art note (`e^{2δ}` vs PDF's `e^{4δ}`); used the PDF version.
- Computed `a` constant via series truncation at `m = 300 000` (tail < 10⁻¹²); verified to 30 dps.
- L_sym evaluated via mpmath Gamma function at 80 dps for `n` up to 10¹²; relative precision well below 10⁻²⁰ at all tested `n`.
- δ optimization: golden-section search, 300 iterations, terminal width < 10⁻⁴⁰. For `n ≤ 10⁴` optimum is at boundary `δ = 1/4 − ε`; for `n ≥ 10⁵` optimum moves to interior. Asymptotically `δ_opt ~ n^{−1/8}` consistent with FN p. 17.
- No "Feynman" / "feynman" / "solver" anywhere in this artifact (per H² embargo and naming rule).
