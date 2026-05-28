# Round 20 — Reproducibility Check

## Method

Copied `genus_growth_sweep.py` and `transform_diagnostic.py` from the PC bundle into a clean temp directory (`/tmp/round20-repro/`). Ran each script with local Python 3.9.6 / numpy 2.0.2. Compared resulting JSON reports against PC's submitted versions.

## G1 — Genus-growth conditioning probe

Structural identity:

- `verdict`: `FAIL_CANONICAL_MONOMIAL_AT_OR_BEFORE_G24` (match)
- `trials` count: 147 (match)
- `sweep_axes`, `quadrature`, `threshold_for_first_crossing`: all identical
- `first_crossings_in_genus_per_eps_jitter`: all 21 entries report `first_genus_crossing_1e10_or_singular: 24` (match)

Numerical values — per-genus max cond(M):

| genus | PC value         | local rerun       | relative diff |
|-------|------------------|-------------------|---------------|
| 2     | 2.747573e+00     | 2.747573e+00      | 4.85e-16      |
| 4     | 1.867468e+01     | 1.867468e+01      | 1.90e-16      |
| 8     | 1.440372e+03     | 1.440372e+03      | 3.79e-15      |
| 12    | 1.443703e+05     | 1.443703e+05      | 3.19e-14      |
| 16    | 1.580721e+07     | 1.580721e+07      | 2.82e-12      |
| 20    | 1.806932e+09     | 1.806932e+09      | 2.21e-11      |
| 24    | 2.118886e+11     | 2.118886e+11      | 9.67e-11      |

Relative agreement degrades from ~1e-16 at g=2 to ~1e-10 at g=24. This is the expected f64 ULP-amplification regime: ULP error in the SVD computation scales as roughly `eps_mach × cond(M)`, so at g=24 where cond ~ 2e11, agreement to 10 sig digits is exactly what f64 supports. **No anomaly.**

**Verdict: REPRODUCIBLE at f64 precision.**

## G2 — Non-tautological transform diagnostic

Structural identity:

- `verdict`: `DEFENSIVE_PIVOT_TO_CHEBYSHEV_RESCALED_IS_ACTIONABLE_UNDER_EQUILIBRATION` (match)
- `option_chosen`: `(b) basis-change matrix conditioning, monomial -> Chebyshev-rescaled` (match)
- `rejected_alternatives`: same reasoning for deferring (a) and (c)
- 21 (eps, jitter) configurations × 7 genera = 147 result rows (match)

Numerical values — g=24 aggregates:

| field                                          | PC value     | local rerun  | match |
|-----------------------------------------------|--------------|--------------|-------|
| `g24_max_block_cond_equilibrated`              | 1.1762e+07   | 1.1762e+07   | yes (to displayed precision) |
| `g24_max_block_cond_raw`                       | 4.6305e+118  | 4.6305e+118  | yes (to displayed precision) |
| `weighted_qr_failure_cond` (reference)         | 4.7e+16      | 4.7e+16      | yes (constant)               |

The **equilibrated** cond reproduces exactly to displayed precision because PC's basis-change matrix is built in `Fraction` arithmetic — there is no f64 round-off during construction. Only the final SVD on the f64-cast matrix carries any numerical noise, which appears below the displayed precision in this case.

The **raw** cond at e+118 is a deterministic property of the rational coefficients of `(mid + half·t)^n` at `g=24, eps=1e-4`. Building in f64 directly would produce `inf` here (catastrophic cancellation); Fraction arithmetic preserves the exact value through to the SVD. Reproducibility confirms this.

**Verdict: REPRODUCIBLE at the precision Fraction arithmetic affords.**

## Performance

| script                   | wall time     | dominant cost            |
|--------------------------|---------------|--------------------------|
| `genus_growth_sweep.py`  | ~5 seconds    | 147 × O(g³) SVDs in f64 |
| `transform_diagnostic.py`| ~64 seconds   | 441 × O(g²) Fraction matrix construction + 441 × O(g³) SVDs |

Both fast enough for routine reruns. The G2 cost is dominated by `Fraction` arithmetic at high genus (the binomial expansion of `(mid + half·t)^n` for `n` up to 23 produces rationals with very large numerator/denominator at small `half`); replacing `Fraction` with a fixed-precision library (`mpmath`, ~50 digits) would reduce wall time by perhaps 5-10× without changing the diagnostic.

## What reproducibility does NOT establish

- The scripts faithfully implement what PC claims; this does not validate the underlying mathematical model.
- Both sweeps remain `F64_SAMPLED_ONLY`. Re-running yields the same f64 numbers but does not produce an interval certificate.
- Reproducibility of a FAIL conditioning result (G1 at g=24) is consistent with both "the unrescaled monomial basis is genuinely unsuitable at g=24" and "the f64 arithmetic itself is failing." The Vandermonde-theory grounding in `CONDITIONING_LITERATURE.md` (Gautschi 1990, Pan 2016) is what licenses the first reading over the second.
- Reproducibility of the G2 actionable verdict (equilibrated cond ~1.18e7) does not certify the Chebyshev-rescaled period matrix's own conditioning — only the basis change between the two families. Period-matrix conditioning of `M_T` is a separate computation deferred to G3 / future local-only work.

## Local rerun environment

```text
python3 --version       → Python 3.9.6
numpy.__version__       → 2.0.2
platform                → macOS / Apple Silicon (Accelerate BLAS)
working directory       → /tmp/round20-repro/
input files             → genus_growth_sweep.py (sha256 c94dd90...90f8)
                          transform_diagnostic.py (sha256 e41a50d0...2afb6)
output paths            → /tmp/round20-repro/GENUS_GROWTH_CONDITIONING_REPORT.json
                          /tmp/round20-repro/TRANSFORM_DIAGNOSTIC_REPORT.json
```
