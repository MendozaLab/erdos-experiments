#!/usr/bin/env python3
"""Build a candidate failure-box JSON for Erdos #1041 in the #114 boundary-slice
schema, then run the cross-problem WS-01 probe on it.

Probe 3 of `~/.claude/plans/yes-wondrous-blum.md` (Intractability test program,
2026-05-15). Internal diagnostic only. NOT a proof. NOT a Lean statement.
Python-level small-N (n=8) test at mpmath dps=30. Box half-width R = 1e-4.

Canonical #1041 statement (per STRUCTURAL_FIT_CHECK_2026-05-15.md):
  Let f(z) = prod_{i=1..n} (z - z_i) with |z_i| < 1 for all i. Must there
  always exist a path of length less than 2 in {z in C : |f(z)| < 1}
  which connects two of the roots of f? (open conjecture)

Method:
  Base polynomial p_0(z) = prod_{k=1..n} (z - z_k), with z_k = r * exp(2 pi i
  (k - 0.5) / n), where r = 0.85 (well inside the unit disk). n=8 gives an
  8-petal lemniscate with critical points well inside the unit disk by
  Gauss-Lucas (the critical points are the roots of p_0').

  Boundary points: real-analytic points on the lemniscate {|p(z)| = 1}, found
  along radial rays out to infinity from the origin. For each (theta, k_rot)
  pair on a 4x3 grid (theta varies the angular sweep, k_rot is a small
  Gaussian rotation of the lemniscate), we find a boundary point z_b on the
  lemniscate by 1D radial bisection of log|p(r e^{i theta})| = 0, then box a
  2D region of half-width R around z_b in (Re z, Im z) coordinates and
  compute interval enclosures of F, F_x, F_y, F_xx, F_yy, F_xy, etc using
  mpmath.iv at dps=30.

  F(z) = |p(z)|^2 - 1 (real-valued; zero on the lemniscate)
  grad F = 2 * (Re(pbar p'), Im(pbar p'))
  F_xx + F_yy + F_xy from second derivatives of |p|^2

  At a real-analytic boundary point z_b in the smooth part of the lemniscate
  (away from critical points of p), grad F is bounded away from zero and the
  WS-01 wall test applies in the same form as for #114.

Honest scope: the kernel is mpmath dps=30 Python interval arithmetic on a
single base polynomial family. It is NOT the certified inari Rust kernel.
Interval blowup, if any, is reported but not catastrophic for this diagnostic.
"""
from __future__ import annotations

import cmath
import json
import math
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from mpmath import iv, mp, mpc, mpf

# -----------------------------------------------------------------------------
# Configuration.
# -----------------------------------------------------------------------------

iv.dps = 30
mp.dps = 50

N_ROOTS = 8
R_BASE = mpf("0.85")  # all roots strictly inside the unit disk

# Angular sweep around the lemniscate (4 angles outside the symmetric petals).
THETA_GRID = [
    mpf("0.10"),   # near positive real axis
    mpf("0.50"),
    mpf("0.95"),
    mpf("1.40"),   # ~80 degrees off real axis
]
# Small rotation perturbations to spread the boxes (in radians).
K_ROT_GRID = [mpf("0.00"), mpf("0.05"), mpf("0.10")]

R_BOX = mpf("1e-4")
SAFETY_FACTOR_S = float(R_BOX)

ERDOS_EXP_DIR = Path(
    "/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments"
)
PROBE_SCRIPT = ERDOS_EXP_DIR / "cross_problem_ws01_applicability_probe.py"
OUT_BASE = ERDOS_EXP_DIR / "Erdos1041" / "proof_path"
ARTIFACT_DIRNAME = "EXP-MATH-ERDOS-1041-WS01-COLD-TEST-20260515-01"
STAGING_JSON = OUT_BASE / "INPUT_boundary_slice_1041_2026-05-15.json"


# -----------------------------------------------------------------------------
# Base polynomial: roots evenly spaced on a circle of radius r inside unit disk.
# -----------------------------------------------------------------------------

def base_roots(n: int, r: mpf) -> list:
    """Return mpmath complex roots: r * exp(2 pi i (k - 0.5) / n) for k=1..n."""
    roots = []
    for k in range(1, n + 1):
        angle = mp.mpf(2) * mp.pi * (mp.mpf(k) - mp.mpf("0.5")) / mp.mpf(n)
        z = mp.mpc(r * mp.cos(angle), r * mp.sin(angle))
        roots.append(z)
    return roots


def rotated_roots(roots: list, k_rot: mpf) -> list:
    """Multiply all roots by exp(i * k_rot)."""
    factor = mp.mpc(mp.cos(k_rot), mp.sin(k_rot))
    return [factor * r for r in roots]


# -----------------------------------------------------------------------------
# Float evaluations (for boundary-point bracket searches).
# -----------------------------------------------------------------------------

def log_abs_p_float(x: float, y: float, roots: list) -> float:
    """log|p(z)| at z = x + i*y. Returns Python float."""
    z = mp.mpc(x, y)
    s = mp.mpf("0")
    for r in roots:
        d = z - r
        a = abs(d)
        if a == 0:
            return float("-inf")
        s += mp.log(a)
    return float(s)


# -----------------------------------------------------------------------------
# Interval evaluations.
#
# We represent z = x + i y where x, y are real intervals. p(z) = prod (z - r_k)
# is evaluated as a sequence of complex multiplications with real/imag parts
# tracked as separate intervals. From p(z) = u + i v we compute
#   F(z) = u^2 + v^2 - 1     (real-valued |p|^2 - 1)
#   F_x = 2 (u u_x + v v_x)
#   F_y = 2 (u u_y + v v_y)
# where u, v, u_x, v_x, u_y, v_y are obtained from the chain
#   p(z) = prod_k (z - r_k);  p'(z) = sum_k prod_{j != k} (z - r_j);
# and the Cauchy-Riemann identity u_x = v_y, u_y = -v_x (so we only need
# the complex derivative p'(z) = u_x + i v_x).
#
# For second derivatives:
#   F_xx + F_yy = 4 |p'|^2 + 4 Re(pbar * p'')   (Laplacian)
#   F_xx - F_yy = 4 Re(pbar * p'' - p'^2 * ??)
# but for the WS-01 probe we only need bounds on |F_tt| (tangent) and the
# third directional derivative. A simpler and sound enclosure is
#   |F_xx|, |F_yy|, |F_xy| <= 2 |p'|^2 + 2 |p| |p''|   (componentwise)
# so the upper bound on |F_tt| <= the same scalar for any unit direction.
# We use this as our f_tt_abs_upper.
# Third-directional bound similarly:
#   |F_ttt| <= 6 |p'| |p''| + 2 |p| |p'''|
# again upper-bounded for any unit direction.
# -----------------------------------------------------------------------------

def p_interval(x_iv, y_iv, roots: list) -> tuple:
    """Return (u_iv, v_iv) for p(z) = u + i v over the box (x_iv, y_iv).

    Each factor (z - r_k) = (x - Re r_k) + i (y - Im r_k), and we multiply
    them out as we go.
    """
    u = iv.mpf(1)
    v = iv.mpf(0)
    for r in roots:
        ar = iv.mpf(str(mp.re(r)))
        br = iv.mpf(str(mp.im(r)))
        # current p = u + i v; factor = (x_iv - ar) + i (y_iv - br)
        fr_re = x_iv - ar
        fr_im = y_iv - br
        new_u = u * fr_re - v * fr_im
        new_v = u * fr_im + v * fr_re
        u = new_u
        v = new_v
    return u, v


def p_and_derivs_interval(x_iv, y_iv, roots: list) -> dict:
    """Return interval enclosures of p, p', p'', p''' over the box.

    Computed by direct interval-arithmetic of the polynomial in factored form.
    Returns a dict with keys 'p_u', 'p_v', 'dp_u', 'dp_v', 'ddp_u', 'ddp_v',
    'dddp_u', 'dddp_v' (real and imaginary parts of each).
    """
    n = len(roots)
    # We use the standard recursion: if p_k(z) = prod_{j=1}^k (z - r_j),
    # then p_{k+1}(z) = (z - r_{k+1}) * p_k(z); and
    # p'_{k+1}(z) = (z - r_{k+1}) * p'_k(z) + p_k(z);
    # p''_{k+1}(z) = (z - r_{k+1}) * p''_k(z) + 2 p'_k(z);
    # p'''_{k+1}(z) = (z - r_{k+1}) * p'''_k(z) + 3 p''_k(z);
    p_u, p_v = iv.mpf(1), iv.mpf(0)
    dp_u, dp_v = iv.mpf(0), iv.mpf(0)
    ddp_u, ddp_v = iv.mpf(0), iv.mpf(0)
    dddp_u, dddp_v = iv.mpf(0), iv.mpf(0)
    for r in roots:
        ar = iv.mpf(str(mp.re(r)))
        br = iv.mpf(str(mp.im(r)))
        fr_re = x_iv - ar
        fr_im = y_iv - br

        # new_p'''   = fr * p''' + 3 * p''
        new_dddp_u = fr_re * dddp_u - fr_im * dddp_v + iv.mpf(3) * ddp_u
        new_dddp_v = fr_re * dddp_v + fr_im * dddp_u + iv.mpf(3) * ddp_v
        # new_p''    = fr * p'' + 2 * p'
        new_ddp_u = fr_re * ddp_u - fr_im * ddp_v + iv.mpf(2) * dp_u
        new_ddp_v = fr_re * ddp_v + fr_im * ddp_u + iv.mpf(2) * dp_v
        # new_p'     = fr * p' + p
        new_dp_u = fr_re * dp_u - fr_im * dp_v + p_u
        new_dp_v = fr_re * dp_v + fr_im * dp_u + p_v
        # new_p      = fr * p
        new_p_u = fr_re * p_u - fr_im * p_v
        new_p_v = fr_re * p_v + fr_im * p_u

        p_u, p_v = new_p_u, new_p_v
        dp_u, dp_v = new_dp_u, new_dp_v
        ddp_u, ddp_v = new_ddp_u, new_ddp_v
        dddp_u, dddp_v = new_dddp_u, new_dddp_v
    return {
        "p_u": p_u, "p_v": p_v,
        "dp_u": dp_u, "dp_v": dp_v,
        "ddp_u": ddp_u, "ddp_v": ddp_v,
        "dddp_u": dddp_u, "dddp_v": dddp_v,
    }


def iv_to_pair(x_iv) -> tuple[float, float]:
    lo = float(x_iv.a)
    hi = float(x_iv.b)
    return lo, hi


def iv_abs_lower(x_iv) -> float:
    lo, hi = iv_to_pair(x_iv)
    if lo <= 0.0 <= hi:
        return 0.0
    return min(abs(lo), abs(hi))


def iv_abs_upper(x_iv) -> float:
    lo, hi = iv_to_pair(x_iv)
    return max(abs(lo), abs(hi))


def iv_mod_upper(u_iv, v_iv) -> float:
    """Upper bound on sqrt(u^2 + v^2) over the box."""
    u_hi = iv_abs_upper(u_iv)
    v_hi = iv_abs_upper(v_iv)
    return math.sqrt(u_hi * u_hi + v_hi * v_hi)


def iv_mod_lower(u_iv, v_iv) -> float:
    """Conservative lower bound on sqrt(u^2 + v^2) over the box.

    sqrt(lower(u^2) + lower(v^2)) where lower(u^2) = (iv_abs_lower(u))^2.
    """
    u_lo = iv_abs_lower(u_iv)
    v_lo = iv_abs_lower(v_iv)
    return math.sqrt(u_lo * u_lo + v_lo * v_lo)


# -----------------------------------------------------------------------------
# Find a real-analytic boundary point on the lemniscate.
# -----------------------------------------------------------------------------

def find_boundary_point_on_ray(roots: list, theta: float) -> tuple[float, float] | None:
    """Bisect log|p(rho * e^{i theta})| = 0 for rho in [0.95, 3.0]."""
    f_lo = log_abs_p_float(0.95 * math.cos(theta), 0.95 * math.sin(theta), roots)
    f_hi = log_abs_p_float(3.0 * math.cos(theta), 3.0 * math.sin(theta), roots)
    if f_lo * f_hi > 0:
        return None
    a, b = 0.95, 3.0
    fa = f_lo
    for _ in range(200):
        m = 0.5 * (a + b)
        fm = log_abs_p_float(m * math.cos(theta), m * math.sin(theta), roots)
        if abs(fm) < 1e-14 or (b - a) < 1e-15:
            return (m * math.cos(theta), m * math.sin(theta))
        if fa * fm < 0:
            b = m
        else:
            a = m
            fa = fm
    rho = 0.5 * (a + b)
    return (rho * math.cos(theta), rho * math.sin(theta))


# -----------------------------------------------------------------------------
# Build one failure-box dict for a (theta, k_rot) pair.
# -----------------------------------------------------------------------------

def build_failure_box(theta: mpf, k_rot: mpf, src_index: int) -> dict:
    roots = rotated_roots(base_roots(N_ROOTS, R_BASE), k_rot)
    boundary = find_boundary_point_on_ray(roots, float(theta))
    if boundary is None:
        return _failed_box(theta, k_rot, src_index, "no_boundary_point_found")
    x_b, y_b = boundary

    if not (math.isfinite(x_b) and math.isfinite(y_b)):
        return _failed_box(theta, k_rot, src_index, "boundary_point_not_finite")

    R = R_BOX
    R_f = float(R)
    x_iv = iv.mpf([str(mp.mpf(x_b) - R), str(mp.mpf(x_b) + R)])
    y_iv = iv.mpf([str(mp.mpf(y_b) - R), str(mp.mpf(y_b) + R)])

    # Interval evaluations.
    derivs = p_and_derivs_interval(x_iv, y_iv, roots)
    p_u, p_v = derivs["p_u"], derivs["p_v"]
    dp_u, dp_v = derivs["dp_u"], derivs["dp_v"]
    ddp_u, ddp_v = derivs["ddp_u"], derivs["ddp_v"]
    dddp_u, dddp_v = derivs["dddp_u"], derivs["dddp_v"]

    # F = |p|^2 - 1 over the box.
    # |p|^2 in [|p|_lo^2, |p|_hi^2]
    p_mod_hi = iv_mod_upper(p_u, p_v)
    p_mod_lo = iv_mod_lower(p_u, p_v)
    F_lo = p_mod_lo * p_mod_lo - 1.0
    F_hi = p_mod_hi * p_mod_hi - 1.0

    # grad F = 2 (Re(pbar p'), Im(pbar p')).
    # Re(pbar p') = p_u * dp_u + p_v * dp_v
    # Im(pbar p') = p_u * dp_v - p_v * dp_u
    rebar = p_u * dp_u + p_v * dp_v
    imbar = p_u * dp_v - p_v * dp_u
    Fx = iv.mpf(2) * rebar
    Fy = iv.mpf(2) * imbar

    # Center-point evaluations (thin intervals around (x_b, y_b)).
    x_thin = iv.mpf([str(mp.mpf(x_b) - mpf("1e-20")), str(mp.mpf(x_b) + mpf("1e-20"))])
    y_thin = iv.mpf([str(mp.mpf(y_b) - mpf("1e-20")), str(mp.mpf(y_b) + mpf("1e-20"))])
    deriv_thin = p_and_derivs_interval(x_thin, y_thin, roots)
    rebar_thin = deriv_thin["p_u"] * deriv_thin["dp_u"] + deriv_thin["p_v"] * deriv_thin["dp_v"]
    imbar_thin = deriv_thin["p_u"] * deriv_thin["dp_v"] - deriv_thin["p_v"] * deriv_thin["dp_u"]
    Fx_c = iv.mpf(2) * rebar_thin
    Fy_c = iv.mpf(2) * imbar_thin

    # |grad F| at center: lower bound = sqrt(Fx_c_lo^2 + Fy_c_lo^2)
    grad_lo = math.sqrt(iv_abs_lower(Fx_c) ** 2 + iv_abs_lower(Fy_c) ** 2)
    center_fn_abs_lower = grad_lo

    # |F_tt|_upper bound. We use the sound enclosure
    #   |F_xx|, |F_yy|, |F_xy| <= 2 |p'|^2_box + 2 |p|_box * |p''|_box
    # so for any unit direction t, |F_tt| <= same scalar.
    dp_mod_hi = iv_mod_upper(dp_u, dp_v)
    ddp_mod_hi = iv_mod_upper(ddp_u, ddp_v)
    dddp_mod_hi = iv_mod_upper(dddp_u, dddp_v)
    f_tt_abs_upper = 2.0 * dp_mod_hi * dp_mod_hi + 2.0 * p_mod_hi * ddp_mod_hi

    # Third directional: |F_ttt| <= 6 |p'| |p''| + 2 |p| |p'''|.
    third_directional_upper = (
        6.0 * dp_mod_hi * ddp_mod_hi + 2.0 * p_mod_hi * dddp_mod_hi
    )

    # Normal remainder bound: third_directional_upper * R^3 / 6 (interval Taylor R3).
    normal_remainder_bound = (third_directional_upper / 6.0) * (R_f ** 3)

    # Wall LHS (old convention) = sup_box |F| + R_n.
    F_abs_upper = max(abs(F_lo), abs(F_hi))
    wall_lhs = F_abs_upper + normal_remainder_bound

    # Wall RHS = S * lower(|grad F|).
    wall_rhs = SAFETY_FACTOR_S * center_fn_abs_lower

    return {
        "f_interval": {"lo": float(F_lo), "hi": float(F_hi)},
        "failing_inequality": (
            f"sup |F| + R_n = {wall_lhs:.6g} compared to S * lower(|grad F|) = "
            f"{wall_rhs:.6g}; theta={float(theta):.3f} k_rot={float(k_rot):.3f}"
        ),
        "fn_full_first_order_interval": {
            "lo": iv_to_pair(Fx)[0], "hi": iv_to_pair(Fx)[1],
        },
        "ft_full_first_order_interval": {
            "lo": iv_to_pair(Fy)[0], "hi": iv_to_pair(Fy)[1],
        },
        "inequality": {
            "center_fn_abs_lower": center_fn_abs_lower,
            "tangent_radius_upper": R_f,
            "f_tt_abs_upper": f_tt_abs_upper,
            "third_directional_upper": third_directional_upper,
            "normal_remainder_bound": normal_remainder_bound,
            "wall_lhs": wall_lhs,
            "wall_rhs": wall_rhs,
            "center_strip_f_abs_upper": F_abs_upper,
        },
        "normal": {"x": 0.0, "y": 1.0},
        "reason": "third_order_wall_separation_failed",
        "source_index": src_index,
        "source_ownership_key": f"theta={float(theta):.3f}/krot={float(k_rot):.3f}",
        "source_reason": "third_order_wall_separation_failed",
        "split_path": f"1041/cold/theta={float(theta):.3f}/krot={float(k_rot):.3f}",
        "tangent": {"x": 1.0, "y": 0.0},
        "x_interval": [float(mp.mpf(x_b) - R), float(mp.mpf(x_b) + R)],
        "y_interval": [float(mp.mpf(y_b) - R), float(mp.mpf(y_b) + R)],
        "_diag": {
            "x_b_float": float(x_b),
            "y_b_float": float(y_b),
            "p_mod_at_center_upper": p_mod_hi,
            "p_mod_at_center_lower": p_mod_lo,
            "n_roots": len(roots),
        },
    }


def _failed_box(theta: mpf, k_rot: mpf, src_index: int, reason: str) -> dict:
    return {
        "f_interval": {"lo": 1.0, "hi": 2.0},
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
        "source_ownership_key": f"theta={float(theta):.3f}/krot={float(k_rot):.3f}",
        "source_reason": "construction_failed",
        "split_path": f"1041/cold/theta={float(theta):.3f}/krot={float(k_rot):.3f}/FAILED",
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
    for theta in THETA_GRID:
        for k_rot in K_ROT_GRID:
            print(f"[build] (theta={float(theta):.3f}, k_rot={float(k_rot):.3f}) ...", flush=True)
            box = build_failure_box(theta, k_rot, src_index)
            print(
                f"        x_interval={box['x_interval']}, "
                f"y_interval={box['y_interval']}, "
                f"center_fn_abs_lower={box['inequality']['center_fn_abs_lower']:.4g}, "
                f"f_tt_up={box['inequality']['f_tt_abs_upper']:.4g}",
                flush=True,
            )
            boxes.append(box)
            src_index += 1

    boundary_slice = {
        "schema_version": "1.0",
        "problem_id": 1041,
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "generator": "build_1041_boundary_slice_cold_test.py",
        "claim_ceiling": (
            "Internal Python-level diagnostic. NOT a proof. NOT a Lean "
            "statement. mpmath dps=30 interval kernel on n=8 base "
            "polynomial with roots evenly spaced on a circle of radius "
            "0.85 inside the unit disk. Box half-width R = 1e-4. "
            "Generated for Probe 3 cross-problem cold-test of the WS-01 "
            "applicability probe."
        ),
        "n_roots": N_ROOTS,
        "r_base": float(R_BASE),
        "theta_grid": [float(t) for t in THETA_GRID],
        "k_rot_grid": [float(k) for k in K_ROT_GRID],
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

    out_dir = OUT_BASE / ARTIFACT_DIRNAME
    if out_dir.exists() and any(out_dir.iterdir()):
        print(
            f"[build] ERROR: artifact dir already exists and is non-empty: {out_dir}",
            file=sys.stderr,
        )
        return 1

    cmd = [
        sys.executable,
        str(PROBE_SCRIPT),
        "--problem-id",
        "1041",
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
