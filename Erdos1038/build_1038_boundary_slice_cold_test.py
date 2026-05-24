#!/usr/bin/env python3
"""Build a candidate failure-box JSON for Erdos #1038 in the #114 boundary-slice
schema, then run the cross-problem WS-01 probe on it.

Track B3 of `~/.claude/plans/yes-wondrous-blum.md` (Tao-gap-closing plan).
Internal diagnostic only. NOT a proof. NOT a Lean statement. Python-level
small-N (N=8) test at mpmath dps=30. The B1 verdict's 1e-8 enclosure
pre-condition is NOT met by this script -- this is a smaller-scale validity
check meant to answer:

  "If we run the B2 probe on a cold #1038 input that imitates the #114
  failure-box schema, does the WS-01 architecture even apply?"

Canonical #1038 statement (per B1 STRUCTURAL_FIT_CHECK):
  inf over nonconstant monic polynomials f with all real roots in [-1, 1] of
  |{x in R : |f(x)| < 1}|

Method:
  Base polynomial p_0(x) = prod_{k=1..N} (x - cos(k*pi/(N+1))), N=8 Chebyshev
  type interior roots in (-1, 1).
  Perturbation: shift the left half of roots inward by delta*w (so w controls
  the fraction shifted and delta controls magnitude). This generates the
  "boundary components merge/split" failure mode flagged in B1's Feature 2.
  Boundary points: real zeros of F(x) = log|p(x)|, i.e., zeros of |p(x)|^2 - 1.
  For each (w, delta) pair (a 4x3 grid centered roughly at B1's (0.174, 0.200)),
  we find a real boundary point with a small bracket search at modest precision,
  then box a region of half-width R around it and compute:
    - f_interval = interval enclosure of F over [x_b - R, x_b + R]
    - center_fn_abs_lower = mpmath-interval lower bound on |F'(x_b)|
    - f_tt_abs_upper, third_directional_upper, normal_remainder_bound
    - wall_lhs, wall_rhs = R * lower(|F'|)  (safety factor S = R)
  Then emit the 12 boxes as remaining_unresolved[*] items with
  reason = 'third_order_wall_separation_failed', and run the B2 probe.

Honest scope: the kernel is mpmath dps=30 Python interval arithmetic on a
single base polynomial. It is NOT the certified inari Rust kernel referenced
in B1. Interval blowup, if any, is reported but not catastrophic for this
diagnostic.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from mpmath import iv, mp, mpf

# -----------------------------------------------------------------------------
# Configuration.
# -----------------------------------------------------------------------------

# Modest precision for a Python-level diagnostic. NOT the B1-stretch 1e-8.
iv.dps = 30
mp.dps = 50  # crisper float backbone for bracket searches

N_ROOTS = 8

# 2D parameter sweep around B1's center (w*, delta*) ~ (0.174, 0.200).
W_GRID = [0.10, 0.13, 0.16, 0.19]
DELTA_GRID = [0.10, 0.15, 0.20]

# Box radius around the boundary point (in x).
R_BOX = mpf("1e-4")

# Safety factor for the wall RHS. Matches #114 convention wall_rhs = S * lower(|F_n|).
SAFETY_FACTOR_S = float(R_BOX)

ERDOS_EXP_DIR = Path("/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments")
PROBE_SCRIPT = ERDOS_EXP_DIR / "cross_problem_ws01_applicability_probe.py"
OUT_BASE = ERDOS_EXP_DIR / "Erdos1038" / "proof_path"
ARTIFACT_DIRNAME = "EXP-MATH-ERDOS-1038-WS01-COLD-TEST-20260514-01"
STAGING_JSON = (
    ERDOS_EXP_DIR
    / "Erdos1038"
    / "proof_path"
    / "_staging_1038_boundary_slice_2026-05-14.json"
)


# -----------------------------------------------------------------------------
# Base polynomial (real-roots-only monic in [-1, 1]).
# -----------------------------------------------------------------------------

def chebyshev_like_roots(n: int) -> list:
    """Return mpmath floats: cos(k * pi / (n + 1)) for k = 1..n, all in (-1, 1)."""
    return [
        mp.cos(mp.mpf(k) * mp.pi / mp.mpf(n + 1)) for k in range(1, n + 1)
    ]


BASE_ROOTS = chebyshev_like_roots(N_ROOTS)


def perturbed_roots(w: float, delta: float) -> list:
    """Shift the left k roots inward by delta, where k = ceil(w * N).
    Resulting polynomial is still monic with real roots in [-1, 1] when
    w, delta are modest.
    """
    n = len(BASE_ROOTS)
    k_shift = max(1, math.ceil(w * n))
    # Sort ascending so 'left' means smaller x.
    sorted_roots = sorted(BASE_ROOTS)
    shifted = []
    for i, r in enumerate(sorted_roots):
        if i < k_shift:
            new_r = r + mp.mpf(delta)  # shift inward (left -> right)
            # Clamp to [-1, 1] (still real roots in interval).
            if new_r > mp.mpf("0.99"):
                new_r = mp.mpf("0.99")
            shifted.append(new_r)
        else:
            shifted.append(r)
    return shifted


def log_p_float(x: float, roots: list) -> float:
    """Floating evaluation of log|p(x)| with mpmath precision; returns Python float."""
    x_mp = mp.mpf(x)
    s = mp.mpf("0")
    for r in roots:
        d = x_mp - r
        if d == 0:
            return float("-inf")
        s += mp.log(abs(d))
    return float(s)


def log_p_iv(x_iv, roots: list):
    """Interval evaluation of log|p(x)| over x_iv using mpmath.iv.

    Each factor (x - r) is interval-evaluated; |.| is the interval absolute
    value (lo = max(0, min(|lo|, |hi|)) when interval straddles 0, else
    just abs of the endpoints). Caller is responsible for ensuring the
    interval does not contain a root.
    """
    s = iv.mpf(0)
    for r in roots:
        # r is exact (mpmath mp.mpf); convert to interval mid point.
        r_iv = iv.mpf(str(r))
        d = x_iv - r_iv
        # abs interval
        d_abs = iv.absmin(d)  # may not exist; do manually:
        # mpmath.iv has |.|; use built-in
        abs_d = abs(d)  # iv.mpf supports __abs__
        # log of an interval that contains 0 -> diverges. Reject if so.
        if iv.mpf(0) in abs_d:
            return None
        s = s + iv.log(abs_d)
    return s


def d_log_p_iv(x_iv, roots: list):
    """Interval evaluation of (log|p|)'(x) = sum 1/(x - r_i)."""
    s = iv.mpf(0)
    for r in roots:
        r_iv = iv.mpf(str(r))
        d = x_iv - r_iv
        if iv.mpf(0) in d:
            return None
        s = s + (iv.mpf(1) / d)
    return s


def d2_log_p_iv(x_iv, roots: list):
    """Interval evaluation of (log|p|)''(x) = -sum 1/(x - r_i)^2."""
    s = iv.mpf(0)
    for r in roots:
        r_iv = iv.mpf(str(r))
        d = x_iv - r_iv
        if iv.mpf(0) in d:
            return None
        s = s - (iv.mpf(1) / (d * d))
    return s


def d3_log_p_iv(x_iv, roots: list):
    """Interval evaluation of (log|p|)'''(x) = 2 * sum 1/(x - r_i)^3."""
    s = iv.mpf(0)
    for r in roots:
        r_iv = iv.mpf(str(r))
        d = x_iv - r_iv
        if iv.mpf(0) in d:
            return None
        s = s + (iv.mpf(2) / (d * d * d))
    return s


# -----------------------------------------------------------------------------
# Find a real boundary point x_b such that log|p(x_b)| = 0, i.e. |p(x_b)| = 1.
# Bracket-then-bisect over a known sign-change interval in x outside the root
# support.
# -----------------------------------------------------------------------------

def find_boundary_point(roots: list, x_lo: float, x_hi: float, tol: float = 1e-12) -> float | None:
    """Bisection on the float-evaluation of log|p|. Returns x_b or None."""
    f_lo = log_p_float(x_lo, roots)
    f_hi = log_p_float(x_hi, roots)
    if f_lo * f_hi > 0:
        return None
    a, b = x_lo, x_hi
    fa = f_lo
    for _ in range(200):
        m = 0.5 * (a + b)
        fm = log_p_float(m, roots)
        if abs(fm) < tol or (b - a) < tol:
            return m
        if fa * fm < 0:
            b = m
        else:
            a = m
            fa = fm
    return 0.5 * (a + b)


def search_boundary_outside_support(roots: list) -> float | None:
    """Search for the rightmost real boundary point (where |p| crosses 1
    moving outward from the root support to +infinity)."""
    rmax = float(max(roots))
    # log|p(rmax+eps)| is very negative; at x=2 it's large positive. Bracket.
    x_lo = rmax + 1e-6
    x_hi = 3.0
    return find_boundary_point(roots, x_lo, x_hi)


# -----------------------------------------------------------------------------
# Build one failure-box dict for a given (w, delta).
# -----------------------------------------------------------------------------

def iv_to_pair(x_iv) -> tuple[float, float]:
    """Best-effort conversion of an iv.mpf interval to (lo_float, hi_float).

    mpmath's iv.mpf endpoints (x_iv.a, x_iv.b) are themselves zero-width
    ivmpf objects; float() on them returns the underlying real number.
    """
    lo = float(x_iv.a)
    hi = float(x_iv.b)
    return lo, hi


def iv_abs_lower(x_iv) -> float:
    """Lower bound on |x| for an interval x_iv. If interval straddles 0,
    returns 0; else returns min of |endpoints|."""
    lo, hi = iv_to_pair(x_iv)
    if lo <= 0.0 <= hi:
        return 0.0
    return min(abs(lo), abs(hi))


def iv_abs_upper(x_iv) -> float:
    """Upper bound on |x| for an interval x_iv."""
    lo, hi = iv_to_pair(x_iv)
    return max(abs(lo), abs(hi))


def build_failure_box(
    w: float,
    delta: float,
    src_index: int,
) -> dict:
    """Construct one #114-schema failure box for a (w, delta) parameter pair.

    Honest reporting fields:
      source_ownership_key  = f"w={w:.3f}/d={delta:.3f}"   (#1038 stand-in)
      source_index          = src_index
      split_path            = "1038/cold/w={w}/delta={delta}"
      x_interval            = [x_b - R, x_b + R]
      y_interval            = [0, 0]  (1D problem; no spatial y)
      tangent, normal       = (1, 0) and (0, 1) sentinel (1D, no real frame)
      f_interval            = interval enclosure of F = log|p(x)| over the box
      inequality.center_fn_abs_lower
                            = lower bound on |F'(x_b)| (the 1D 'normal' grad)
      inequality.tangent_radius_upper
                            = R (the box half-width)
      inequality.f_tt_abs_upper
                            = upper bound on |F''| over the box
      inequality.third_directional_upper
                            = upper bound on |F'''| over the box
      inequality.normal_remainder_bound
                            = a small residual bound (mid-box residual after
                              subtracting interval-Taylor terms)
      inequality.wall_lhs   = sup_box |F| + R_n  (matches #114 convention)
      inequality.wall_rhs   = R * lower(|F'|)    (S = R, the box scale)
      inequality.center_strip_f_abs_upper
                            = sup_box |F|  (for reporting)
    """
    roots = perturbed_roots(w, delta)
    x_b_float = search_boundary_outside_support(roots)
    if x_b_float is None or not math.isfinite(x_b_float):
        return _failed_box(w, delta, src_index, "no_boundary_point_found")

    # Build a box of half-width R around x_b.
    x_b = mp.mpf(x_b_float)
    R = R_BOX
    x_iv = iv.mpf([str(x_b - R), str(x_b + R)])

    # Interval evaluations.
    f_iv = log_p_iv(x_iv, roots)
    if f_iv is None:
        return _failed_box(w, delta, src_index, "log_pole_in_box")

    df_iv = d_log_p_iv(x_iv, roots)
    if df_iv is None:
        return _failed_box(w, delta, src_index, "first_deriv_pole_in_box")

    d2f_iv = d2_log_p_iv(x_iv, roots)
    if d2f_iv is None:
        return _failed_box(w, delta, src_index, "second_deriv_pole_in_box")

    d3f_iv = d3_log_p_iv(x_iv, roots)
    if d3f_iv is None:
        return _failed_box(w, delta, src_index, "third_deriv_pole_in_box")

    # Center-point gradient lower bound.
    x_b_iv_thin = iv.mpf([str(x_b - mpf("1e-20")), str(x_b + mpf("1e-20"))])
    df_at_xb = d_log_p_iv(x_b_iv_thin, roots)
    if df_at_xb is None:
        return _failed_box(w, delta, src_index, "first_deriv_pole_at_center")

    center_fn_abs_lower = iv_abs_lower(df_at_xb)
    f_tt_abs_upper = iv_abs_upper(d2f_iv)
    third_directional_upper = iv_abs_upper(d3f_iv)

    # Normal remainder bound (residual after subtracting first/second order
    # Taylor expansion around x_b). For interval Taylor with remainder, the
    # third-derivative bound times R^3 / 6 is a sound enclosure.
    R_f = float(R)
    normal_remainder_bound = (third_directional_upper / 6.0) * (R_f ** 3)

    # Wall LHS = sup_box |F| + R_n  (matches #114 convention).
    f_abs_upper = iv_abs_upper(f_iv)
    wall_lhs = f_abs_upper + normal_remainder_bound

    # Wall RHS = S * lower(|F'|), with S = box half-width R.
    wall_rhs = SAFETY_FACTOR_S * center_fn_abs_lower

    f_lo_float, f_hi_float = iv_to_pair(f_iv)

    return {
        "f_interval": {"lo": f_lo_float, "hi": f_hi_float},
        "failing_inequality": (
            f"sup |F| + R_n = {wall_lhs:.6g} compared to S * lower(|F'|) = "
            f"{wall_rhs:.6g}; w={w} delta={delta}"
        ),
        "fn_full_first_order_interval": {
            "lo": iv_to_pair(df_iv)[0],
            "hi": iv_to_pair(df_iv)[1],
        },
        "ft_full_first_order_interval": {
            "lo": iv_to_pair(df_iv)[0],
            "hi": iv_to_pair(df_iv)[1],
        },
        "inequality": {
            "center_fn_abs_lower": center_fn_abs_lower,
            "tangent_radius_upper": R_f,
            "f_tt_abs_upper": f_tt_abs_upper,
            "third_directional_upper": third_directional_upper,
            "normal_remainder_bound": normal_remainder_bound,
            "wall_lhs": wall_lhs,
            "wall_rhs": wall_rhs,
            "center_strip_f_abs_upper": f_abs_upper,
        },
        "normal": {"x": 0.0, "y": 1.0},
        "reason": "third_order_wall_separation_failed",
        "source_index": src_index,
        "source_ownership_key": f"w={w:.3f}/d={delta:.3f}",
        "source_reason": "third_order_wall_separation_failed",
        "split_path": f"1038/cold/w={w}/delta={delta}",
        "tangent": {"x": 1.0, "y": 0.0},
        "x_interval": [float(x_b - R), float(x_b + R)],
        "y_interval": [0.0, 0.0],
        "_diag": {
            "x_b_float": float(x_b_float),
            "n_roots": len(roots),
            "roots_left_shifted": math.ceil(w * len(roots)),
        },
    }


def _failed_box(w: float, delta: float, src_index: int, reason: str) -> dict:
    """Emit a degenerate box with a clear failure tag. The probe will report
    this as BRANCH_POINT_NOT_FOUND or schema-validation failure depending on
    fields. We provide zeros so the probe's classify_box runs cleanly."""
    return {
        "f_interval": {"lo": 1.0, "hi": 2.0},  # does NOT straddle zero
        "failing_inequality": f"box construction failed: {reason}",
        "fn_full_first_order_interval": {"lo": 0.0, "hi": 0.0},
        "ft_full_first_order_interval": {"lo": 0.0, "hi": 0.0},
        "inequality": {
            "center_fn_abs_lower": 0.0,
            "tangent_radius_upper": float(R_BOX),
            "f_tt_abs_upper": 0.0,
            "third_directional_upper": 0.0,
            "normal_remainder_bound": 0.0,
            "wall_lhs": 0.0,
            "wall_rhs": 0.0,
            "center_strip_f_abs_upper": 0.0,
        },
        "normal": {"x": 0.0, "y": 1.0},
        "reason": "third_order_wall_separation_failed",
        "source_index": src_index,
        "source_ownership_key": f"w={w:.3f}/d={delta:.3f}",
        "source_reason": "construction_failed",
        "split_path": f"1038/cold/w={w}/delta={delta}/FAILED",
        "tangent": {"x": 1.0, "y": 0.0},
        "x_interval": [0.0, 0.0],
        "y_interval": [0.0, 0.0],
        "_diag": {"construction_failure": reason},
    }


# -----------------------------------------------------------------------------
# Main.
# -----------------------------------------------------------------------------

def main() -> int:
    boxes = []
    src_index = 0
    for w in W_GRID:
        for delta in DELTA_GRID:
            print(f"[build] (w={w}, delta={delta}) ...", flush=True)
            box = build_failure_box(w, delta, src_index)
            print(
                f"        x_interval={box['x_interval']}, "
                f"center_fn_abs_lower={box['inequality']['center_fn_abs_lower']:.4g}, "
                f"f_tt_up={box['inequality']['f_tt_abs_upper']:.4g}",
                flush=True,
            )
            boxes.append(box)
            src_index += 1

    boundary_slice = {
        "schema_version": "1.0",
        "problem_id": 1038,
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "generator": "build_1038_boundary_slice_cold_test.py",
        "claim_ceiling": (
            "Internal Python-level diagnostic. NOT a proof. NOT a Lean "
            "statement. mpmath dps=30 interval kernel on N=8 perturbed "
            "Chebyshev-like roots. Box half-width R = 1e-4. Generated for "
            "Track B3 cold-test of the cross-problem WS-01 applicability "
            "probe (B2)."
        ),
        "n_roots": N_ROOTS,
        "w_grid": W_GRID,
        "delta_grid": DELTA_GRID,
        "box_half_width_R": float(R_BOX),
        "safety_factor_S": SAFETY_FACTOR_S,
        "remaining_unresolved": boxes,
    }

    STAGING_JSON.parent.mkdir(parents=True, exist_ok=True)
    STAGING_JSON.write_text(
        json.dumps(boundary_slice, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"[build] wrote staging JSON: {STAGING_JSON}", flush=True)

    # Run the B2 probe.
    out_dir = OUT_BASE / ARTIFACT_DIRNAME
    if out_dir.exists():
        print(
            f"[build] ERROR: artifact dir already exists: {out_dir}",
            file=sys.stderr,
        )
        return 1

    cmd = [
        sys.executable,
        str(PROBE_SCRIPT),
        "--problem-id",
        "1038",
        "--input-json",
        str(STAGING_JSON),
        "--output-dir",
        str(out_dir),
    ]
    print(f"[build] running probe: {' '.join(cmd)}", flush=True)
    rc = subprocess.run(cmd, check=False).returncode
    print(f"[build] probe exit code: {rc}", flush=True)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
