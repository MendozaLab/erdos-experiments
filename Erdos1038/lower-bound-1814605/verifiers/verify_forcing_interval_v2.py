#!/usr/bin/env python3
"""Independent verifier for erdos1038.conservative_forcing_interval.v2-fullrange.


This script re-checks, with a DIFFERENT interval backend (mpmath.iv) than the
generator (directed interval arithmetic), an emitted v2 forcing certificate:

  1. STRUCTURE  — metadata consistency (leaf_count, split_boxes binary-tree
     relation, max_depth).
  2. DOMAIN     — every leaf box lies inside [a_min, a_max] x [0,1] x [0,1],
     with the f64 domain endpoints verified to be OUTWARD roundings of the
     decimal endpoints (superset coverage of [-1.708, -sqrt(2)]).
  3. PARTITION  — leaves tile the domain exactly: pairwise positive-volume
     overlap check (exact f64 comparisons, vectorized) + EXACT total-volume
     equality via Fraction arithmetic (f64 endpoints are converted exactly).
  4. INEQUALITY — per leaf, recompute the rigorous lower bound of U_nu over
     the box with mpmath interval arithmetic (decimal constants enclosed from
     strings) and require it strictly positive; report how many clear the
     declared target margin.

Per-check verdict is PASS / FAIL / NEED_DATA with diagnostics; overall PASS
requires every check PASS. Never silently passes on missing fields.

Usage: verify_forcing_interval_v2.py <certificate.json>
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction

import mpmath
from mpmath import iv

iv.prec = 80

K = iv.mpf("1.395")
CSHIFT = iv.mpf("1.071")
EPS0 = iv.mpf("0.0001")
ONE = iv.mpf(1)
SQRT2_NEG_DEC = "-1.414213562373095"  # rational endpoint > -sqrt(2): superset


def ivpt(x: float):
    return iv.mpf(x)  # exact: f64 embeds exactly at 80-bit precision


def neg_log_lower(d_sup: float):
    """Rigorous lower bound of -log(d_max) given d_max <= d_sup, d_max > 0."""
    if not d_sup > 0.0:
        return None
    return (-iv.log(ivpt(d_sup))).a


def product_lower(w_iv, base_lo):
    """Lower bound of w*t with w in w_iv (>0) and t >= base_lo (t may be +inf)."""
    if base_lo >= 0:
        return (iv.mpf([w_iv.a, w_iv.a]) * iv.mpf([base_lo, base_lo])).a
    return (iv.mpf([w_iv.b, w_iv.b]) * iv.mpf([base_lo, base_lo])).a


def lower_bound_on_box(leaf, cap_iv):
    a = iv.mpf([leaf["a_lo"], leaf["a_hi"]])
    s = iv.mpf([leaf["s_lo"], leaf["s_hi"]])

    b = s * (cap_iv + a)
    w = K - b
    if not w.a > 0:
        return None
    c_loc = CSHIFT - b

    m1ma = -ONE - a
    if not m1ma.a > 0:
        return None
    n = -EPS0 - iv.log(m1ma) - w * iv.log(ONE + b)
    d = iv.log(CSHIFT + ONE - b)
    if not (n.a > 0 and d.a > 0):
        return None
    cw = n / d
    if not cw.a > 0:
        return None

    x_lo, x_hi = leaf["x_lo"], leaf["x_hi"]

    d_a = (ivpt(x_hi) - ivpt(leaf["a_lo"])).b
    base_a = neg_log_lower(float(mpmath.mpf(d_a)))
    if base_a is None:
        return None

    e1 = abs(ivpt(x_lo) - iv.mpf([b.b, b.b])).b
    e2 = abs(ivpt(x_hi) - iv.mpf([b.a, b.a])).b
    base_b = neg_log_lower(float(mpmath.mpf(max(e1, e2))))
    if base_b is None:
        return None

    f1 = abs(ivpt(x_lo) - iv.mpf([c_loc.b, c_loc.b])).b
    f2 = abs(ivpt(x_hi) - iv.mpf([c_loc.a, c_loc.a])).b
    base_c = neg_log_lower(float(mpmath.mpf(max(f1, f2))))
    if base_c is None:
        return None

    total = (
        iv.mpf([base_a, base_a])
        + iv.mpf([product_lower(w, base_b)] * 2)
        + iv.mpf([product_lower(cw, base_c)] * 2)
    )
    return total.a


def main(path: str) -> int:
    checks = {}

    payload = json.loads(open(path, encoding="utf-8").read())

    # ---- 1. STRUCTURE ----
    try:
        assert payload["format"] == "erdos1038.conservative_forcing_interval.v2-fullrange"
        leaves = payload["leaves"]
        assert payload["leaf_count"] == len(leaves) > 0
        assert payload["split_boxes"] == len(leaves) - 1, "binary-tree leaf relation"
        assert payload["max_depth"] == max(l["depth"] for l in leaves)
        wlb = min(l["lower_bound"] for l in leaves)
        assert payload["worst_lower_bound"] == wlb
        checks["structure"] = ("PASS", f"{len(leaves)} leaves, max_depth {payload['max_depth']}")
    except (AssertionError, KeyError) as e:
        checks["structure"] = ("FAIL", repr(e))
        return report(checks)

    a_min, a_max = payload["a_min_f64"], payload["a_max_f64"]
    cap_dec = payload["cap_decimal"]
    cap_iv = iv.mpf(cap_dec)
    margin = payload["target_margin_f64"]

    # ---- 2. DOMAIN (incl. outward-rounding of decimal endpoints) ----
    try:
        amin_dec = iv.mpf(payload["a_min_decimal"])
        amax_dec = iv.mpf(payload["a_max_decimal"])
        assert ivpt(a_min).a <= amin_dec.a, "a_min_f64 must round outward (down)"
        assert ivpt(a_max).b >= amax_dec.b, "a_max_f64 must round outward (up)"
        # superset of the required range: a_max must be >= -sqrt(2)
        assert payload["a_max_decimal"] == SQRT2_NEG_DEC
        assert ivpt(a_max).b >= (-iv.sqrt(iv.mpf(2))).a, "a_max must exceed -sqrt(2)"
        assert payload["a_min_decimal"] in ("-1.708", "-1.7")
        for i, l in enumerate(leaves):
            assert a_min <= l["a_lo"] <= l["a_hi"] <= a_max, f"leaf {i} a-range"
            assert 0.0 <= l["s_lo"] <= l["s_hi"] <= 1.0, f"leaf {i} s-range"
            assert 0.0 <= l["x_lo"] <= l["x_hi"] <= 1.0, f"leaf {i} x-range"
            assert l["a_lo"] < l["a_hi"] and l["s_lo"] < l["s_hi"] and l["x_lo"] < l["x_hi"]
        checks["domain"] = ("PASS", f"a in [{a_min!r}, {a_max!r}] superset of [{payload['a_min_decimal']}, -sqrt2]")
    except (AssertionError, KeyError) as e:
        checks["domain"] = ("FAIL", repr(e))
        return report(checks)

    # ---- 3. PARTITION ----
    try:
        try:
            import numpy as np
        except ImportError:
            checks["partition"] = ("NEED_DATA", "numpy unavailable for overlap check")
            return report(checks)
        arr = np.array(
            [[l["a_lo"], l["a_hi"], l["s_lo"], l["s_hi"], l["x_lo"], l["x_hi"]] for l in leaves]
        )
        n = len(arr)
        overlaps = 0
        for i in range(n):
            lo = np.maximum(arr[i, 0::2], arr[i + 1 :, 0::2])
            hi = np.minimum(arr[i, 1::2], arr[i + 1 :, 1::2])
            overlaps += int(np.count_nonzero(np.all(lo < hi, axis=1)))
        assert overlaps == 0, f"{overlaps} pairwise positive-volume overlaps"
        vol = sum(
            (Fraction(l["a_hi"]) - Fraction(l["a_lo"]))
            * (Fraction(l["s_hi"]) - Fraction(l["s_lo"]))
            * (Fraction(l["x_hi"]) - Fraction(l["x_lo"]))
            for l in leaves
        )
        dom = (Fraction(a_max) - Fraction(a_min)) * 1 * 1
        assert vol == dom, f"exact volume mismatch: {float(vol - dom)}"
        checks["partition"] = ("PASS", "disjoint + exact volume coverage (Fraction)")
    except AssertionError as e:
        checks["partition"] = ("FAIL", repr(e))
        return report(checks)

    # ---- 4. INEQUALITY (independent backend recompute) ----
    worst = None
    n_pos = n_margin = 0
    fail = None
    for i, l in enumerate(leaves):
        lb = lower_bound_on_box(l, cap_iv)
        if lb is None or not lb > 0:
            fail = (i, l, None if lb is None else float(mpmath.mpf(lb)))
            break
        lbf = float(mpmath.mpf(lb))
        n_pos += 1
        if lbf > margin:
            n_margin += 1
        worst = lbf if worst is None else min(worst, lbf)
    if fail is not None:
        checks["inequality"] = ("FAIL", f"leaf {fail[0]} not certified positive by mpmath.iv: {fail[2]}, box={fail[1]}")
    else:
        checks["inequality"] = (
            "PASS",
            f"all {n_pos} leaves > 0 (mpmath.iv recompute); {n_margin} clear margin {margin}; worst {worst:.6e}",
        )
    return report(checks)


def report(checks) -> int:
    overall = "PASS"
    for name, (verdict, msg) in checks.items():
        print(f"{name:>10}: {verdict} — {msg}")
        if verdict != "PASS":
            overall = "FAIL" if verdict == "FAIL" else ("NEED_DATA" if overall == "PASS" else overall)
    for name in ("structure", "domain", "partition", "inequality"):
        if name not in checks:
            print(f"{name:>10}: NEED_DATA — check not reached")
            overall = "FAIL" if overall == "FAIL" else "NEED_DATA"
    print(f"overall: {overall}")
    return 0 if overall == "PASS" else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
