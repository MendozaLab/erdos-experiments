#!/usr/bin/env python3
"""
Track A2 — Rigorous Positive Rational Lower Bound on Tao's c_10 Constant
========================================================================

Experiment:  EXP-MATH-EHP114-TAO-DEFECT-PSI-COR-C10-LOWER-BOUND-20260509-01

Purpose
-------
Compute a certified (interval-arithmetic) lower bound on the universal disk
integral

    c_10 = ∫∫_{D(0, 1/21)} f(z) dA,
    f(z) = 1/|z| + 1/|z-1| - |1/z + 1/(z-1)|,

where D(0, 1/21) is the OPEN complex disk of radius 1/21 about the origin.
This is the constant Tao calls c_C evaluated at C = 10 in the corollary
"defect-psi-cor" of arXiv:2512.12455v2 (line ~1025).

Math preconditions verified before computing:

(P1) f(z) >= 0 for all z in C \ {0, 1}. By the triangle inequality,
        |1/z + 1/(z-1)| <= |1/z| + |1/(z-1)| = 1/|z| + 1/|z-1|.

(P2) The only singularity of f inside D(0, 1/21) is at z = 0. For
     |z| < 1/21,  |z - 1| >= 1 - |z| > 20/21,  so 1/|z-1| < 21/20.

(P3) Singularity is integrable in polar coordinates: with z = r e^{i theta},
     dA = r dr dtheta. Algebraic simplification:

        1/z + 1/(z-1) = (2z-1) / (z (z-1))    =>    |1/z + 1/(z-1)| = |2z-1| / (|z| |z-1|)

     So
        f(z) = 1/|z| + 1/|z-1| - |2z-1| / (|z| |z-1|)
             = (|z-1| + |z| - |2z-1|) / (|z| |z-1|).

     Multiplying by r = |z|:
        g(r, theta) := f(re^{i theta}) * r
                     = (|z-1| + |z| - |2z-1|) / |z-1|     (where z = r e^{i theta}).

     g is bounded and continuous on the CLOSED disk: at z=0, g = (1+0-1)/1 = 0.
     So c_10 = ∫_0^{2 pi} ∫_0^{1/21} g(r, theta) dr dtheta on a continuous
     bounded integrand --- no singular handling needed.

Procedure (rigorous lower bound)
--------------------------------
1. Use mpmath.iv (interval / ball arithmetic) at high precision (dps = 60).
2. Partition the rectangle [0, 1/21] x [0, 2 pi] into N_r * N_theta cells.
3. For each cell, build an INTERVAL bounding box [r_lo, r_hi] x [theta_lo, theta_hi]
   and evaluate g on it via mpmath.iv. The result is an interval enclosure
   [g_lo, g_hi] containing the true minimum and maximum of g over the cell.
4. The integral over the cell is bounded below by g_lo * (cell area, lower
   end of width interval). Multiplying interval widths by interval g-values
   yields a sound lower bound on the cell integral (mpmath.iv handles all
   directed-rounding internally).
5. Sum across cells; the final mpf is a certified lower bound of c_10.
6. Also report a non-rigorous fine-grid midpoint estimate for sanity check.

f(z) >= 0 (P1) means g >= 0 too. So if a cell's interval enclosure
[g_lo, g_hi] returns g_lo < 0, that's a numerical-overestimate artifact of
interval widening (the dependency problem); we floor it at 0 only when we
can independently verify g >= 0 there. Default: do NOT clamp; if g_lo is
negative we let it contribute negatively to the lower bound. This is honest:
if the grid is too coarse, the rigorous bound may be < 0, in which case
we report status BLOCKED.

Output
------
- RESULTS.json       (machine-readable, schema_version 1.0)
- REPORT.md          (human-readable)
- RESULTS.sha256     (sha256 over RESULTS.json)
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import mpmath as mp
from mpmath import iv, mpf, mp as mpgblc

# ----------------------------------------------------------------------
# Constants and precision
# ----------------------------------------------------------------------
EXPERIMENT_ID = "EXP-MATH-EHP114-TAO-DEFECT-PSI-COR-C10-LOWER-BOUND-20260509-01"
ARTIFACT_DIR = Path(
    "/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/"
    "Erdos114/proof_path/" + EXPERIMENT_ID
)
SCRIPT_PATH = Path(__file__).resolve()

DPS = 60  # decimal places for both float and interval contexts
mpgblc.dps = DPS
iv.dps = DPS

R_MAX = mpf(1) / mpf(21)        # disk radius
TWO_PI = 2 * mp.pi              # used for grid construction (mpf, then converted)

# Grid resolution. Cell-count tuning: small enough to run in a few seconds,
# large enough that interval-widening doesn't drag the lower bound below
# the non-rigorous estimate by more than a factor of ~2-3.
N_R_RIGOROUS = 200
N_THETA_RIGOROUS = 720

# Non-rigorous midpoint estimate: finer grid for a tight reference.
N_R_NONRIGOROUS = 2000
N_THETA_NONRIGOROUS = 7200


# ----------------------------------------------------------------------
# Interval integrand
# ----------------------------------------------------------------------
def integrand_interval(r_iv, theta_iv):
    """
    g(r, theta) = (|z - 1| + |z| - |2z - 1|) / |z - 1|, with z = r * e^{i theta}.

    All inputs are mpmath.iv intervals. Returns an interval enclosing
    the range of g over (r_iv, theta_iv).

    Identity: |z| = r (real positive), |z-1|^2 = r^2 - 2 r cos(theta) + 1,
              |2z-1|^2 = 4 r^2 - 4 r cos(theta) + 1.
    """
    cos_t = iv.cos(theta_iv)

    # |z| = r (already nonneg interval)
    abs_z = r_iv

    # |z - 1|^2 = r^2 - 2 r cos(theta) + 1
    sq_zm1 = r_iv * r_iv - 2 * r_iv * cos_t + 1
    # Numerical safety: |z-1|^2 must be >= (1 - 1/21)^2 = (20/21)^2; clamp lower
    # by intersecting with a known positive interval to avoid sqrt of weakly
    # negative interval boundary due to widening.
    lo_safe = (mpf(20) / mpf(21)) ** 2
    sq_zm1 = iv.mpf([max(mpf(sq_zm1.a), lo_safe), mpf(sq_zm1.b)])
    abs_zm1 = iv.sqrt(sq_zm1)

    # |2z - 1|^2 = 4 r^2 - 4 r cos(theta) + 1
    sq_2zm1 = 4 * r_iv * r_iv - 4 * r_iv * cos_t + 1
    # |2z - 1| = sqrt(...). At r=1/2, theta=0 this hits zero, but for r <= 1/21
    # we have |2z-1| >= 1 - 2 r >= 1 - 2/21 = 19/21 > 0.
    lo_safe2 = (mpf(19) / mpf(21)) ** 2
    sq_2zm1 = iv.mpf([max(mpf(sq_2zm1.a), lo_safe2), mpf(sq_2zm1.b)])
    abs_2zm1 = iv.sqrt(sq_2zm1)

    numerator = abs_zm1 + abs_z - abs_2zm1
    g = numerator / abs_zm1
    return g


def integrand_float(r, theta):
    """Non-rigorous mpmath float evaluation of g(r, theta). Used for sanity grid."""
    cos_t = mp.cos(theta)
    abs_z = r
    abs_zm1 = mp.sqrt(r * r - 2 * r * cos_t + 1)
    abs_2zm1 = mp.sqrt(4 * r * r - 4 * r * cos_t + 1)
    return (abs_zm1 + abs_z - abs_2zm1) / abs_zm1


# ----------------------------------------------------------------------
# Rigorous lower bound via interval grid
# ----------------------------------------------------------------------
def rigorous_lower_bound(n_r: int, n_theta: int):
    """Returns (lower_bound_mpf, num_cells_with_negative_lower)."""
    print(f"[rigorous] N_r = {n_r}, N_theta = {n_theta}, dps = {DPS}")

    # Build grid edges as mpf
    r_edges = [R_MAX * mpf(i) / mpf(n_r) for i in range(n_r + 1)]
    theta_edges = [TWO_PI * mpf(j) / mpf(n_theta) for j in range(n_theta + 1)]

    total_lo = iv.mpf(0)
    neg_cells = 0

    # Iterate cells. For speed, we keep cell area in interval form too.
    for i in range(n_r):
        r_lo, r_hi = r_edges[i], r_edges[i + 1]
        r_iv = iv.mpf([r_lo, r_hi])
        dr_iv = iv.mpf([r_hi - r_lo, r_hi - r_lo])  # exact width

        for j in range(n_theta):
            t_lo, t_hi = theta_edges[j], theta_edges[j + 1]
            t_iv = iv.mpf([t_lo, t_hi])
            dt_iv = iv.mpf([t_hi - t_lo, t_hi - t_lo])

            g_iv = integrand_interval(r_iv, t_iv)

            # Cell integral lower bound: g.a (minimum over cell) * dr * dtheta.
            g_lower = mpf(g_iv.a)
            cell_lower = iv.mpf([g_lower, g_lower]) * dr_iv * dt_iv
            total_lo += cell_lower

            if g_lower < 0:
                neg_cells += 1

        if (i + 1) % max(1, n_r // 10) == 0:
            mid = (mpf(total_lo.a) + mpf(total_lo.b)) / 2
            print(f"  [{i + 1}/{n_r}] running lower-end estimate ~ {mp.nstr(mid, 8)}, "
                  f"neg-cells so far = {neg_cells}")

    return mpf(total_lo.a), neg_cells


def nonrigorous_midpoint_estimate(n_r: int, n_theta: int):
    """Standard midpoint Riemann sum at mpf precision. Not certified."""
    print(f"[non-rigorous] midpoint, N_r = {n_r}, N_theta = {n_theta}")
    dr = R_MAX / mpf(n_r)
    dtheta = TWO_PI / mpf(n_theta)
    total = mpf(0)
    for i in range(n_r):
        r_mid = (mpf(i) + mpf("0.5")) * dr
        for j in range(n_theta):
            t_mid = (mpf(j) + mpf("0.5")) * dtheta
            total += integrand_float(r_mid, t_mid)
    return total * dr * dtheta


# ----------------------------------------------------------------------
# Spot-check: f(z) >= 0 numerically on a sample grid
# ----------------------------------------------------------------------
def spot_check_nonnegativity(n_samples: int = 50):
    """Evaluate g on a sample grid; min should be >= 0. Confirms P1."""
    print(f"[spot-check] sampling g on {n_samples} x {n_samples} grid")
    dr = R_MAX / mpf(n_samples)
    dtheta = TWO_PI / mpf(n_samples)
    min_g = mpf("1e9")
    max_g = mpf(0)
    for i in range(n_samples):
        r = (mpf(i) + mpf("0.5")) * dr
        for j in range(n_samples):
            theta = (mpf(j) + mpf("0.5")) * dtheta
            g = integrand_float(r, theta)
            if g < min_g:
                min_g = g
            if g > max_g:
                max_g = g
    print(f"  min g = {mp.nstr(min_g, 12)}, max g = {mp.nstr(max_g, 12)}")
    return float(min_g), float(max_g)


# ----------------------------------------------------------------------
# Convert lower bound to a positive rational (for the JSON record)
# ----------------------------------------------------------------------
def to_positive_rational_lower(x_mpf, decimals: int = 30):
    """
    Returns (num, den) such that num/den <= x_mpf and num/den is positive.
    We truncate (round-down) the decimal expansion to 'decimals' places.

    For x = 0.0123456..., decimals=10: returns (123456..., 10**10) such that
    num/den is the round-DOWN truncation, hence <= x.
    """
    if x_mpf <= 0:
        return None
    # Use mpmath.floor on x * 10^decimals to be safe with high-precision rounding.
    scaled = x_mpf * mpf(10) ** decimals
    n_int = int(mp.floor(scaled))
    if n_int <= 0:
        # x is very small; shrink representation
        # Find the smallest decimals d such that floor(x * 10^d) >= 1, then return that.
        d = decimals
        while True:
            d += 5
            s = x_mpf * mpf(10) ** d
            n = int(mp.floor(s))
            if n >= 1:
                return [n, 10 ** d]
            if d > 200:
                return None
    return [n_int, 10 ** decimals]


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    if not ARTIFACT_DIR.exists():
        print(f"FATAL: artifact dir does not exist: {ARTIFACT_DIR}", file=sys.stderr)
        sys.exit(1)
    print(f"Artifact dir: {ARTIFACT_DIR}")

    t0 = time.time()
    print("=" * 70)
    print(f"Track A2: rigorous lower bound on Tao c_10 = ∫∫_{{D(0,1/21)}} f(z) dA")
    print(f"          where f(z) = 1/|z| + 1/|z-1| - |1/z + 1/(z-1)|")
    print(f"Math precondition checks:")
    print(f"  (P1) f >= 0 everywhere by triangle inequality.")
    print(f"  (P2) z=1 outside disk: |z-1| > 20/21 for z in D(0, 1/21).")
    print(f"  (P3) Algebraic simplification g = f*r = (|z-1|+|z|-|2z-1|)/|z-1|")
    print(f"       is continuous and bounded on closed disk; no singular strip.")
    print("=" * 70)

    # Step 0: spot-check P1 numerically
    sc_min, sc_max = spot_check_nonnegativity(50)
    if sc_min < -1e-12:
        print(f"FATAL: spot check found negative g (min = {sc_min})")
        sys.exit(1)
    print(f"[spot-check] PASS: g >= {sc_min:.6e} (essentially 0 from below by P1).")

    # Step 1: non-rigorous midpoint estimate (sanity reference)
    t_nr0 = time.time()
    nr_estimate = nonrigorous_midpoint_estimate(N_R_NONRIGOROUS, N_THETA_NONRIGOROUS)
    t_nr = time.time() - t_nr0
    print(f"[non-rigorous] estimate = {mp.nstr(nr_estimate, 16)}  ({t_nr:.1f}s)")

    # Step 2: rigorous interval lower bound
    t_r0 = time.time()
    rig_lower, neg_cells = rigorous_lower_bound(N_R_RIGOROUS, N_THETA_RIGOROUS)
    t_r = time.time() - t_r0
    print(f"[rigorous] lower bound = {mp.nstr(rig_lower, 16)}  ({t_r:.1f}s)")
    print(f"           cells with negative interval-min: {neg_cells} / "
          f"{N_R_RIGOROUS * N_THETA_RIGOROUS}")

    # Step 3: classify
    if rig_lower > 0:
        status = "TAO_C10_LOWER_BOUND_RIGOROUS_PASS"
    elif rig_lower == 0:
        status = "TAO_C10_LOWER_BOUND_NUMERIC_ONLY"
    else:
        status = "TAO_C10_LOWER_BOUND_BLOCKED"
    print(f"[classify] status = {status}")

    # Step 4: consistency check
    if mpf(rig_lower) <= mpf(nr_estimate):
        consistency = "rigorous_LB <= non_rigorous_estimate"
    else:
        consistency = "RIGOROUS_LB_EXCEEDS_ESTIMATE_INVESTIGATE"
    print(f"[consistency] {consistency}")

    # Step 5: rational
    rational = to_positive_rational_lower(rig_lower, decimals=30) \
        if rig_lower > 0 else None
    print(f"[rational] {rational}")

    # Step 6: build artifact JSON
    elapsed = time.time() - t0
    payload = {
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "claim_ceiling": (
            "Internal rigorous lower bound on a universal disk integral. "
            "Not a proof of Erdos #114, not an N0 candidate by itself, and "
            "not a Tao threshold extraction. Produces one numeric input "
            "that propagates into Track A3 and A4."
        ),
        "constant_name": "c_10",
        "definition": (
            "Disk integral over D(0, 1/21) of f(z) = 1/|z| + 1/|z-1| - "
            "|1/z + 1/(z-1)| dA"
        ),
        "rigorous_lower_bound": float(rig_lower),
        "rigorous_lower_bound_mpf_str": mp.nstr(rig_lower, 30),
        "rigorous_lower_bound_rational": rational,
        "rigorous_lower_bound_method": (
            "mpmath.iv ball/interval arithmetic on uniform polar grid; "
            "g = (|z-1| + |z| - |2z-1|) / |z-1| is bounded and "
            "continuous (singular factor 1/|z| canceled by Jacobian r). "
            "Cell lower bound = interval-min(g) * cell width (rounded "
            "outward by mpmath.iv)."
        ),
        "rigorous_method_notes": (
            f"Working precision dps={DPS}. Polar uniform grid "
            f"{N_R_RIGOROUS} radial x {N_THETA_RIGOROUS} angular cells. "
            "Used a sound algebraic rewrite to remove the 1/|z| singularity "
            "before interval evaluation, so no singular strip is needed. "
            "Negative interval-min cells (if any) NOT clamped — let them "
            "drag the lower bound; an honest negative bound would yield "
            "BLOCKED status."
        ),
        "non_rigorous_midpoint_estimate": float(nr_estimate),
        "non_rigorous_midpoint_estimate_mpf_str": mp.nstr(nr_estimate, 16),
        "non_rigorous_method": (
            f"Midpoint Riemann sum on uniform polar grid "
            f"{N_R_NONRIGOROUS} radial x {N_THETA_NONRIGOROUS} angular at "
            f"dps={DPS}; not certified, used only for sanity-check."
        ),
        "estimate_consistency_check": consistency,
        "grid_resolution": {
            "radial_cells": N_R_RIGOROUS,
            "angular_cells": N_THETA_RIGOROUS,
            "singular_strip_radius": 0.0,
            "singular_strip_method": (
                "closed-form (algebraic cancellation of 1/r)"
            ),
        },
        "negative_interval_min_cells": int(neg_cells),
        "total_cells": int(N_R_RIGOROUS * N_THETA_RIGOROUS),
        "computation_time_seconds": float(elapsed),
        "rigorous_walltime_seconds": float(t_r),
        "nonrigorous_walltime_seconds": float(t_nr),
        "library_used": "mpmath (mpmath.iv for intervals, dps=60)",
        "next_dependency": (
            "If status is RIGOROUS_PASS with c_10 > 0, queue Track A3 "
            "(origin-repulsion implied constant). If BLOCKED or "
            "NUMERIC_ONLY, document obstacle and refine grid or method."
        ),
        "tao_paper_reference": {
            "arxiv_id": "2512.12455v2",
            "lemma": "defect-psi-cor",
            "line": 1025,
            "wording": (
                "c_C := ∫∫_{D(0, 1/(2C+1))} (1/|z| + 1/|z-1| - "
                "|1/z + 1/(z-1)|) dA, evaluated here at C = 10."
            ),
        },
        "math_preconditions_verified": {
            "P1_nonnegativity": "f >= 0 by triangle inequality |a+b| <= |a|+|b|",
            "P2_singularity_inside_disk": (
                "Only z=0 lies in D(0,1/21); z=1 has |z-1| > 20/21"
            ),
            "P3_singularity_is_integrable": (
                "g(r,theta) = f * r = (|z-1|+|z|-|2z-1|)/|z-1| continuous bounded; "
                "g(0) = 0; integration over polar rectangle is over a continuous "
                "bounded function."
            ),
            "spot_check_min_g": float(sc_min),
            "spot_check_max_g": float(sc_max),
        },
        "timestamp_unix": int(time.time()),
        "schema_version": "1.0",
    }

    results_path = ARTIFACT_DIR / f"{EXPERIMENT_ID}_RESULTS.json"
    with open(results_path, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"[write] {results_path}")

    # SHA-256
    with open(results_path, "rb") as f:
        sha = hashlib.sha256(f.read()).hexdigest()
    sha_path = ARTIFACT_DIR / f"{EXPERIMENT_ID}_RESULTS.sha256"
    with open(sha_path, "w") as f:
        f.write(f"{sha}  {results_path.name}\n")
    print(f"[write] {sha_path}")

    # REPORT.md
    report_path = ARTIFACT_DIR / f"{EXPERIMENT_ID}_REPORT.md"
    write_report(report_path, payload, rig_lower, nr_estimate, sc_min, sc_max)
    print(f"[write] {report_path}")

    print("=" * 70)
    print(f"DONE. status = {status}")
    print(f"      c_10 >= {mp.nstr(rig_lower, 16)}  (rigorous)")
    print(f"      c_10 ~  {mp.nstr(nr_estimate, 16)}  (non-rigorous estimate)")
    print(f"      total walltime = {elapsed:.1f}s")
    print("=" * 70)


def write_report(path, payload, rig_lower, nr_estimate, sc_min, sc_max):
    md = f"""# Track A2: Rigorous Lower Bound on Tao's c_10 Constant

**Experiment ID:** `{payload["experiment_id"]}`
**Status:** `{payload["status"]}`
**Date (UTC):** {time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(payload["timestamp_unix"]))}

## Quantity

c_10 is the universal disk integral

  c_10 = ∫∫_{{D(0, 1/21)}} f(z) dA,    f(z) = 1/|z| + 1/|z-1| - |1/z + 1/(z-1)|

referenced as c_C at C = 10 in the corollary `defect-psi-cor` of arXiv:2512.12455v2 (line ~1025).

## Math preconditions (verified before computing)

- (P1) f(z) >= 0 by triangle inequality |a+b| <= |a|+|b|.
- (P2) Only singular point in D(0, 1/21) is z = 0. We have |z-1| > 20/21 for z in D(0, 1/21).
- (P3) The Jacobian r in polar coordinates absorbs the 1/|z| singularity. Algebraic simplification:

      f(z) = (|z-1| + |z| - |2z-1|) / (|z| |z-1|),
      g(r, theta) := f * r = (|z-1| + |z| - |2z-1|) / |z-1|

  is continuous and bounded on the closed disk, with g(0) = 0. Hence c_10 = ∫∫ g dr dtheta over the polar rectangle [0, 1/21] x [0, 2 pi] of a CONTINUOUS BOUNDED integrand. No singular-strip handling required.

- Spot-check (50 x 50 sample grid): min g = {sc_min:.6e}, max g = {sc_max:.6e}. PASS.

## Method

- Working precision: mpmath dps = {DPS}.
- Rigorous: mpmath.iv interval arithmetic on uniform {N_R_RIGOROUS} x {N_THETA_RIGOROUS} polar grid. Each cell evaluates the interval enclosure [g_lo, g_hi]; cell contribution to lower bound is g_lo * dr * dtheta (mpmath.iv handles outward rounding internally).
- Non-rigorous: midpoint Riemann sum at {N_R_NONRIGOROUS} x {N_THETA_NONRIGOROUS} for sanity comparison.

## Result

- **Rigorous lower bound:** c_10 >= {mp.nstr(rig_lower, 30)}
- **Non-rigorous midpoint estimate:** c_10 ~ {mp.nstr(nr_estimate, 16)}
- **Positive rational lower bound (truncated to 30 decimals):**
  num/den = {payload["rigorous_lower_bound_rational"]}
- **Consistency check:** {payload["estimate_consistency_check"]}
- **Cells with negative interval-min:** {payload["negative_interval_min_cells"]} / {payload["total_cells"]}

## Honest scope

This is **internal infrastructure**. The result is one number. It is

- NOT a proof of Erdos problem #114.
- NOT an N0 candidate by itself.
- NOT a Tao threshold extraction.

It produces one numeric input that propagates into Track A3 (origin-repulsion implied constant) and Track A4. If the lower bound is strictly positive, Track A3 is unblocked; otherwise we report obstacle and refine.

## Reproduce

```
python3 {SCRIPT_PATH}
```

Walltime: {payload["computation_time_seconds"]:.1f}s total ({payload["rigorous_walltime_seconds"]:.1f}s rigorous + {payload["nonrigorous_walltime_seconds"]:.1f}s non-rigorous).

## Tao reference

Tao (arXiv:2512.12455v2), corollary `defect-psi-cor`, ~line 1025 (snippet stored in JSON `tao_paper_reference.wording`).

## Files

- `{payload["experiment_id"]}_RESULTS.json` (machine readable)
- `{payload["experiment_id"]}_RESULTS.sha256`
- `{payload["experiment_id"]}_REPORT.md` (this file)
"""
    with open(path, "w") as f:
        f.write(md)


if __name__ == "__main__":
    main()
