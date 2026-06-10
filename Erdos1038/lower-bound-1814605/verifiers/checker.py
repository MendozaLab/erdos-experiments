#!/usr/bin/env python3
"""
Clean required-domain checker for a piecewise five-atom dual-forcing tail
certificate (Erdos Problem 1038, finite-atom lower-bound route).

This MIRRORS the official Hua Xu verifier
(reference-huaxu-cert/1038/verify_piecewise_tail_correct_domain.py)
domain logic EXACTLY, so that running it on the reference 1.8146 / 560-block
JSON reproduces the published worst required margin ~= 9.534343713646365e-06.

SIGN / NORMALIZATION CONVENTION (locked from the reference verifier):
  - supp(mu) subset {-1} union [0,1]  (Tao/natso normalization).
  - Dual measure on a block parameterized by a in [-A, -C]:
        lambda_a = delta_a + sum_j w_j delta_{a + d_j}
  - Potential, in the y = x - a variable:
        V(y) = log(1/|y|) + sum_j w_j * log(1/|y - d_j|)
    (the leading delta_a has implicit weight 1; the d_j are the shifts).
  - Required positivity domains per block (the ONLY domains that matter):
        x = -1      ->  y in [C - 1, A - 1]
        x in [0,1]  ->  y in [C, A + 1]
    The middle interval [A-1, C] (x in (-1,0)) is an OVERCHECK region and
    is NOT part of the requirement.
  - Block margin = min over the required domains of V(y); positive => block OK.

The checker is exact-rational-safe in spirit (it isolates poles, uses analytic
V'(y) for critical points via bisection/brentq, and never evaluates at a pole).
It is float64 numerics matching the reference; we report the worst margin and
compare to the published baseline to ~1e-7.
"""
from __future__ import annotations

import argparse
import json
import math
from typing import List, Sequence, Tuple

import numpy as np

try:
    from scipy.optimize import brentq  # type: ignore
except Exception:  # pragma: no cover
    brentq = None


# ---------------------------------------------------------------------------
# Core potential V(y) and its analytic derivative V'(y)
# ---------------------------------------------------------------------------
def _safe_log1_over_abs(x: float) -> float:
    return math.log(1.0 / abs(x))


def V_value(y: float, weights: Sequence[float], shifts: Sequence[float]) -> float:
    """V(y) = log(1/|y|) + sum_j w_j log(1/|y - d_j|). Leading atom weight = 1."""
    return _safe_log1_over_abs(y) + sum(
        float(w) * math.log(1.0 / abs(y - d)) for w, d in zip(weights, shifts)
    )


def g_derivative(y: float, weights: Sequence[float], shifts: Sequence[float]) -> float:
    """Analytic V'(y) = 1/y + sum_j w_j / (y - d_j)."""
    return 1.0 / y + sum(float(w) / (y - d) for w, d in zip(weights, shifts))


# ---------------------------------------------------------------------------
# Root finding for critical points (mirror of reference _bisection_root)
# ---------------------------------------------------------------------------
def _bisection_root(f, a, b, max_iter: int = 80):
    fa = f(a)
    fb = f(b)
    if not (math.isfinite(fa) and math.isfinite(fb)):
        return None
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0:
        return None
    lo, hi = a, b
    flo = fa
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if not math.isfinite(fm):
            return None
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def _interval_samples(shifts: Sequence[float], lo: float, hi: float):
    """Split (lo,hi) at any interior shift (pole), nudging away from poles."""
    pts = [lo, hi] + [d for d in shifts if lo < d < hi]
    pts = sorted(pts)
    out = []
    for a, b in zip(pts[:-1], pts[1:]):
        eps_a = 1e-10 if any(abs(a - d) < 1e-14 for d in shifts) else 0.0
        eps_b = 1e-10 if any(abs(b - d) < 1e-14 for d in shifts) else 0.0
        aa = a + eps_a
        bb = b - eps_b
        if aa < bb:
            out.append((aa, bb))
    return out


def critical_points(intervals, weights, shifts) -> List[float]:
    pts: List[float] = []
    for lo, hi in intervals:
        pts.extend([lo, hi])
        for a, b in _interval_samples(shifts, lo, hi):
            try:
                ga = g_derivative(a, weights, shifts)
                gb = g_derivative(b, weights, shifts)
                if math.isfinite(ga) and math.isfinite(gb) and ga * gb < 0:
                    if brentq is not None:
                        root = brentq(
                            lambda yy: g_derivative(yy, weights, shifts), a, b
                        )
                    else:
                        root = _bisection_root(
                            lambda yy: g_derivative(yy, weights, shifts), a, b
                        )
                    if root is not None:
                        pts.append(root)
            except Exception:
                pass
    # de-duplicate; drop anything sitting on a pole
    uniq: List[float] = []
    for y in pts:
        if all(abs(y - z) > 1e-8 for z in uniq) and all(
            abs(y - d) > 1e-8 for d in shifts
        ):
            uniq.append(y)
    # robustness fallback: if no critical points found, sample densely
    if not uniq:
        for lo, hi in intervals:
            for y in np.linspace(lo, hi, 3000):
                if all(abs(float(y) - d) > 1e-9 for d in shifts):
                    uniq.append(float(y))
    return sorted(uniq)


# ---------------------------------------------------------------------------
# Block margin under the per-block shift set (supports per-block shifts)
# ---------------------------------------------------------------------------
def block_required_min(
    A: float, C: float, weights: Sequence[float], shifts: Sequence[float]
) -> Tuple[float, float]:
    """Return (min V over required domains, argmin y) for one block.

    Required domains: y in [C-1, A-1] (x=-1) and y in [C, A+1] (x in [0,1]).
    """
    intervals = [(C - 1.0, A - 1.0), (C, A + 1.0)]
    pts = critical_points(intervals, weights, shifts)
    vals = [(V_value(y, weights, shifts), y) for y in pts]
    return min(vals, key=lambda t: t[0])


def verify_required_domain(cert: dict):
    """Worst required margin across all blocks.

    Supports both the reference schema (single global ``shifts`` array) and a
    per-block ``shifts`` override (used by the tiling/assignment certificates,
    where d_{j,i} differs per block).
    """
    global_shifts = cert.get("shifts")
    worst = None  # (value, block_i, y)
    bad_blocks = []
    for block in cert["blocks"]:
        A = float(block["A"])
        C = float(block["C"])
        weights = [float(w) for w in block["weights"]]
        shifts = [float(d) for d in block.get("shifts", global_shifts)]
        block_min, block_y = block_required_min(A, C, weights, shifts)
        if worst is None or block_min < worst[0]:
            worst = (block_min, int(block["i"]), block_y)
        if block_min <= 0:
            bad_blocks.append((int(block["i"]), block_min, block_y))
    return worst, bad_blocks


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("json_path")
    ap.add_argument(
        "--expect-margin",
        type=float,
        default=None,
        help="if given, assert worst margin matches to --tol",
    )
    ap.add_argument("--tol", type=float, default=1e-7)
    args = ap.parse_args()

    with open(args.json_path, "r") as f:
        cert = json.load(f)

    worst, bad = verify_required_domain(cert)
    print("Certificate:", cert.get("name", args.json_path))
    print("M =", cert.get("M"), "K =", cert.get("K"))
    print("Required-domain worst margin: {:.15g}".format(worst[0]))
    print("Required-domain worst block:", worst[1], "y = {:.15g}".format(worst[2]))
    print("Bad required-domain blocks:", len(bad))
    print("VERDICT:", "PASS" if not bad else "FAIL")

    if args.expect_margin is not None:
        diff = abs(worst[0] - args.expect_margin)
        ok = diff <= args.tol
        print(
            "GATE: expected {:.15g} got {:.15g} |diff|={:.3g} tol={:.3g} -> {}".format(
                args.expect_margin, worst[0], diff, args.tol, "MATCH" if ok else "MISMATCH"
            )
        )
        return 0 if ok else 2
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
