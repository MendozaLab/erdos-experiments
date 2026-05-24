# Track A2: Rigorous Lower Bound on Tao's c_10 Constant

**Experiment ID:** `EXP-MATH-EHP114-TAO-DEFECT-PSI-COR-C10-LOWER-BOUND-20260509-01`
**Status:** `TAO_C10_LOWER_BOUND_RIGOROUS_PASS`
**Date (UTC):** 2026-05-09 18:56:45

## Quantity

c_10 is the universal disk integral

  c_10 = ∫∫_{D(0, 1/21)} f(z) dA,    f(z) = 1/|z| + 1/|z-1| - |1/z + 1/(z-1)|

referenced as c_C at C = 10 in the corollary `defect-psi-cor` of arXiv:2512.12455v2 (line ~1025).

## Math preconditions (verified before computing)

- (P1) f(z) >= 0 by triangle inequality |a+b| <= |a|+|b|.
- (P2) Only singular point in D(0, 1/21) is z = 0. We have |z-1| > 20/21 for z in D(0, 1/21).
- (P3) The Jacobian r in polar coordinates absorbs the 1/|z| singularity. Algebraic simplification:

      f(z) = (|z-1| + |z| - |2z-1|) / (|z| |z-1|),
      g(r, theta) := f * r = (|z-1| + |z| - |2z-1|) / |z-1|

  is continuous and bounded on the closed disk, with g(0) = 0. Hence c_10 = ∫∫ g dr dtheta over the polar rectangle [0, 1/21] x [0, 2 pi] of a CONTINUOUS BOUNDED integrand. No singular-strip handling required.

- Spot-check (50 x 50 sample grid): min g = 9.378679e-07, max g = 9.882729e-02. PASS.

## Method

- Working precision: mpmath dps = 60.
- Rigorous: mpmath.iv interval arithmetic on uniform 200 x 720 polar grid. Each cell evaluates the interval enclosure [g_lo, g_hi]; cell contribution to lower bound is g_lo * dr * dtheta (mpmath.iv handles outward rounding internally).
- Non-rigorous: midpoint Riemann sum at 2000 x 7200 for sanity comparison.

## Result

- **Rigorous lower bound:** c_10 >= 0.00690049552414053277372859578815
- **Non-rigorous midpoint estimate:** c_10 ~ 0.007069155243213715
- **Positive rational lower bound (truncated to 30 decimals):**
  num/den = [6900495524140532773728595788, 1000000000000000000000000000000]
- **Consistency check:** rigorous_LB <= non_rigorous_estimate
- **Cells with negative interval-min:** 13288 / 144000

## Honest scope

This is **internal infrastructure**. The result is one number. It is

- NOT a proof of Erdos problem #114.
- NOT an N0 candidate by itself.
- NOT a Tao threshold extraction.

It produces one numeric input that propagates into Track A3 (origin-repulsion implied constant) and Track A4. If the lower bound is strictly positive, Track A3 is unblocked; otherwise we report obstacle and refine.

## Reproduce

```
python3 /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/build_tao_c10_disk_integral_lower_bound.py
```

Walltime: 501.7s total (15.1s rigorous + 486.5s non-rigorous).

## Tao reference

Tao (arXiv:2512.12455v2), corollary `defect-psi-cor`, ~line 1025 (snippet stored in JSON `tao_paper_reference.wording`).

## Files

- `EXP-MATH-EHP114-TAO-DEFECT-PSI-COR-C10-LOWER-BOUND-20260509-01_RESULTS.json` (machine readable)
- `EXP-MATH-EHP114-TAO-DEFECT-PSI-COR-C10-LOWER-BOUND-20260509-01_RESULTS.sha256`
- `EXP-MATH-EHP114-TAO-DEFECT-PSI-COR-C10-LOWER-BOUND-20260509-01_REPORT.md` (this file)
