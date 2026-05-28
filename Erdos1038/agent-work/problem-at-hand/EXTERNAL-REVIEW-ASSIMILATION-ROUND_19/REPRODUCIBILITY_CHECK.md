# Round 19 — Reproducibility Check

## Method

Copied `round19_canonical_falsifier_sweep.py` from the PC bundle into a clean temp directory (`/tmp/round19-repro/`). Ran with local Python 3.9.6 / numpy 2.0.2. Compared resulting `NULL_FALSIFIER_REPORT.json` against PC's submitted version.

## Result

Structural identity:
- `status`: NO_FALSIFIER_FOUND_IN_TOY_RANGE (match)
- `trial_count`: 63 (match)
- `coverage`: identical (match)
- `claim_level`: 0 (match)
- `claim_scope`: identical text (match)

Numerical values:

| Field | PC value | Local rerun value | Delta |
|---|---|---|---|
| `max_condition` | 20.81295453370498 | 20.812954533704993 | ~1e-14 |
| `max_transform_condition` | 20.81295453370497 | 20.812954533704982 | ~1e-14 |
| `min_singular_min` | 1.7326932764839196 | 1.7326932764839196 | 0 |

The two condition values differ from PC's only in the last digit. This is consistent with f64 last-bit roundoff arising from different BLAS backend implementations of QR/SVD (Accelerate on macOS vs. whichever backend PC's sandbox used). Structural and quantitative agreement to ~14 significant digits.

**Verdict: REPRODUCIBLE at f64 precision.** No discrepancy beyond float roundoff noise.

## What reproducibility does NOT establish

- Reproducibility verifies the script faithfully implements what PC claims; it does not validate the scientific claim about #1038.
- The script is f64-sampled only. Re-running yields the same f64 numbers but does not produce an interval certificate.
- The script tests genus 2-4 only. Re-running cannot extrapolate to genus 24.
- Reproducibility of a null falsification result does not constitute proof of anything; null result remains scope-limited evidence per PC's own claim_scope.
