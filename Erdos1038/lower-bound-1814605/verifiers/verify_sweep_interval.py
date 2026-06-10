#!/usr/bin/env python3
"""
verify_sweep_interval.py — pure-Python directed-interval re-verifier for the
piecewise five-atom dual-forcing SWEEP certificate (Erdos Problem 1038,
finite-atom lower-bound route).

WHAT THIS PROVES, INDEPENDENTLY OF ANY COMPILED KERNEL
------------------------------------------------------
Given the sweep certificate JSON (e.g. certificates/w2_uniform4_M1814605.json),
this script certifies, block by block, that the dual potential

    V(y) = log(1/|y|) + sum_j w_j * log(1/|y - d_j|)        (leading weight 1)

is strictly positive on BOTH required domains of every block,

    x = -1      ->  y in [C - 1, A - 1]
    x in [0,1]  ->  y in [C,     A + 1]

using ONLY interval arithmetic with outward (directed) rounding, via mpmath's
``iv`` context. Dependencies: Python standard library + mpmath. No numpy, no
scipy, no compiled code. A verdict of CERTIFIED_POSITIVE on all blocks is an
interval-arithmetic certificate of positivity, not a float spot-check.

CONVENTIONS (locked from the bundled float checker, verifiers/checker.py,
which mirrors the official thread verifier)
-------------------------------------------
- supp(mu) subset {-1} union [0,1] (Tao/natso normalization).
- Dual measure on a block parameterized by a in [-A, -C]:
      lambda_a = delta_a + sum_j w_j delta_{a + d_j}
- The middle interval [A-1, C] (x in (-1,0)) is an OVERCHECK region and is
  NOT part of the requirement; it is not checked here.
- CERTIFICATE VALUES ARE EXACT: the IEEE-754 binary64 numbers stored in the
  JSON file ARE the certificate data (every binary64 value is a dyadic
  rational, represented exactly by mpmath). The verifier does not "trust"
  any decimal rendering; json.load yields the exact binary64 values, and
  iv.mpf(x) embeds each as an exact singleton interval. Derived endpoints
  (C-1, A+1, ...) are computed in interval arithmetic and taken OUTWARD, so
  the certified region is a superset of the true required domain.

ALGORITHM (adaptive bisection, fail-closed)
-------------------------------------------
Each required y-interval is bisected adaptively. For a box [a, b], the
certified lower bound is the max of
  (1) the natural interval extension of V over [a, b], and
  (2) the mean-value form V(m) + V'([a,b]) * ([a,b] - m), m = midpoint,
both evaluated with outward rounding. If the lower bound clears 0 the box is
certified and bisection stops there. If a rigorous point enclosure at the
midpoint has sup < 0 the block is REFUTED with a witness. A box that cannot
be resolved within width / depth / box-count budgets fails the whole block
loudly (FAIL_TO_CERTIFY) — never a silent pass. A point or box evaluation
that touches a pole exactly (|y - d| = [0,0]) is treated as unusable
(lower bound -inf / upper bound +inf), again fail-closed.

Certified lower bounds are BOUNDS, not minima: bisection stops the moment a
box clears zero, so per-block values are naturally looser than float-precision
point margins, and a different bisection path (this one vs any other interval
engine) gives different — but equally valid — positive bounds. Verdicts must
agree; bound values need only be positive.

USAGE
-----
    python3 verify_sweep_interval.py certificates/w2_uniform4_M1814605.json
        [--out RECEIPT.json] [--compare certify-output/REFERENCE_RECEIPT.json]
        [--blocks 191,196] [--self-test] [--progress-every 20]
        [--min-width 1e-13] [--max-depth 80] [--max-boxes 200000]
        [--progress-file PROGRESS.log]

Exit code 0 iff every checked block is CERTIFIED_POSITIVE (and, when
--self-test is given, the deliberate-perturbation negative control fires the
REFUTED/FAIL path; and, when --compare is given, per-block verdicts match).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time

from mpmath import iv, mp

# Working precision for interval endpoints. 53 bits = binary64-equivalent
# endpoints; outward rounding at this precision is rigorous (wider intervals
# are always sound).
iv.prec = 53

INF = float("inf")


# ---------------------------------------------------------------------------
# Interval enclosures of V and V'
# ---------------------------------------------------------------------------
def _ln_abs(x_iv):
    """Enclosure of ln|x| for an interval x, or None if |x| is exactly [0,0]
    (singleton on a pole: the value is undefined; fail closed)."""
    a = abs(x_iv)
    if a.b == 0:  # |x| identically zero -> pole hit exactly
        return None
    return iv.log(a)  # a.a == 0 gives [-inf, ln(a.b)] which is a sound enclosure


def v_interval(y_iv, weights_iv, shifts_iv):
    """Natural interval extension of V(y) = -ln|y| - sum_j w_j ln|y - d_j|.

    (Algebraically identical to the log(1/|.|) form; outward rounding
    throughout.) Returns None if any term sits exactly on a pole (fail-closed
    at the caller).
    """
    lead = _ln_abs(y_iv)
    if lead is None:
        return None
    acc = -lead
    for w_iv, d_iv in zip(weights_iv, shifts_iv):
        if w_iv.a == 0 and w_iv.b == 0:
            continue  # exact-zero weight contributes exactly 0
        t = _ln_abs(y_iv - d_iv)
        if t is None:
            return None
        acc = acc - w_iv * t
    return acc


def vprime_interval(y_iv, weights_iv, shifts_iv):
    """Interval enclosure of V'(y) = -1/y - sum_j w_j/(y - d_j).

    Division by an interval straddling 0 yields [-inf, +inf] in mpmath.iv —
    still a rigorous (if useless) enclosure. Returns None only on exceptions.
    """
    one = iv.mpf(1)
    try:
        acc = -(one / y_iv)
        for w_iv, d_iv in zip(weights_iv, shifts_iv):
            if w_iv.a == 0 and w_iv.b == 0:
                continue
            acc = acc - (w_iv / (y_iv - d_iv))
    except Exception:
        return None
    return acc


def box_lower_bound(bl, bh, weights_iv, shifts_iv):
    """Certified lower bound (mpf, may be -inf) of V on the box [bl, bh]:
    max of the natural extension and the mean-value form."""
    y_iv = iv.mpf([bl, bh])
    nat = v_interval(y_iv, weights_iv, shifts_iv)
    nat_lb = nat.a if nat is not None else mp.ninf

    m = 0.5 * (bl + bh)
    m_iv = iv.mpf(m)
    vm = v_interval(m_iv, weights_iv, shifts_iv)
    vp = vprime_interval(y_iv, weights_iv, shifts_iv)
    if vm is not None and vp is not None:
        mv = vm + vp * (y_iv - m_iv)
        mv_lb = mv.a
    else:
        mv_lb = mp.ninf
    return nat_lb if nat_lb > mv_lb else mv_lb


def point_upper_bound(yp, weights_iv, shifts_iv):
    """Rigorous upper bound of V at the point yp (sup < 0 => V(yp) < 0)."""
    v = v_interval(iv.mpf(yp), weights_iv, shifts_iv)
    return v.b if v is not None else mp.inf


# ---------------------------------------------------------------------------
# Required domains, outward-enclosed
# ---------------------------------------------------------------------------
def required_domains_outer(a, c):
    """Outward-rounded float endpoints for [C-1, A-1] and [C, A+1].
    Positivity on these (slightly larger) domains implies positivity on the
    true required domains."""
    ia, ic, one = iv.mpf(a), iv.mpf(c), iv.mpf(1)
    d1_lo, d1_hi = (ic - one), (ia - one)
    d2_lo, d2_hi = ic, (ia + one)
    return [
        (float(d1_lo.a), float(d1_hi.b)),
        (float(d2_lo.a), float(d2_hi.b)),
    ]


# ---------------------------------------------------------------------------
# Per-block certification (adaptive bisection, fail-closed)
# ---------------------------------------------------------------------------
def certify_block(a, c, weights, shifts, i, min_width=1e-13, max_depth=80,
                  max_boxes=200000):
    weights_iv = [iv.mpf(w) for w in weights]
    shifts_iv = [iv.mpf(d) for d in shifts]

    certified_lb = mp.inf
    boxes = 0
    max_depth_reached = 0

    def result(verdict, lb, refute_y=None, refute_ub=None, reason=None):
        return {
            "i": i,
            "verdict": verdict,
            "certified_lower_bound": float(lb),
            "boxes_examined": boxes,
            "max_depth_reached": max_depth_reached,
            "refute_y": refute_y,
            "refute_value_upper": refute_ub,
            "fail_reason": reason,
        }

    for lo, hi in required_domains_outer(a, c):
        if not lo < hi:
            return result("FAIL_TO_CERTIFY", mp.ninf,
                          reason=f"degenerate required domain [{lo}, {hi}]")
        stack = [(lo, hi, 0)]
        while stack:
            bl, bh, depth = stack.pop()
            boxes += 1
            if depth > max_depth_reached:
                max_depth_reached = depth
            if boxes > max_boxes:
                return result("FAIL_TO_CERTIFY", mp.ninf,
                              reason="box budget exceeded")
            lb = box_lower_bound(bl, bh, weights_iv, shifts_iv)
            if lb > 0:
                if lb < certified_lb:
                    certified_lb = lb
                continue  # certified positive on this box
            # try rigorous refutation at the midpoint
            mid = 0.5 * (bl + bh)
            ub = point_upper_bound(mid, weights_iv, shifts_iv)
            if ub < 0:
                return result("REFUTED", mp.ninf,
                              refute_y=mid, refute_ub=float(ub))
            width = bh - bl
            if width < min_width or depth >= max_depth:
                return result(
                    "FAIL_TO_CERTIFY", mp.ninf,
                    reason=(f"unresolvable box [{bl!r}, {bh!r}] width "
                            f"{width:.3e} depth {depth} (lb={float(lb):.3e})"))
            # bisect; midpoint shared by both closed halves covers the parent
            stack.append((bl, mid, depth + 1))
            stack.append((mid, bh, depth + 1))

    return result("CERTIFIED_POSITIVE", certified_lb)


# ---------------------------------------------------------------------------
# Certificate-level run
# ---------------------------------------------------------------------------
def block_shifts(cert, block):
    s = block.get("shifts", cert.get("shifts"))
    if s is None:
        raise ValueError(f"block {block['i']} has no shifts and certificate "
                         "has no global shifts")
    return s


def certify_certificate(cert, block_filter=None, progress=None, **cfg):
    blocks = cert["blocks"]
    if block_filter is not None:
        want = set(block_filter)
        blocks = [b for b in blocks if int(b["i"]) in want]
    if not blocks:
        raise ValueError("no blocks selected")

    t0 = time.time()
    per_block = []
    worst_lb, worst_i = INF, -1
    for n, block in enumerate(blocks):
        r = certify_block(
            float(block["A"]), float(block["C"]),
            [float(w) for w in block["weights"]],
            [float(d) for d in block_shifts(cert, block)],
            int(block["i"]), **cfg)
        per_block.append(r)
        if (r["verdict"] == "CERTIFIED_POSITIVE"
                and r["certified_lower_bound"] < worst_lb):
            worst_lb = r["certified_lower_bound"]
            worst_i = r["i"]
        if progress and ((n + 1) % progress[1] == 0 or n + 1 == len(blocks)):
            progress[0](n + 1, len(blocks), r, worst_lb, worst_i,
                        time.time() - t0)

    n_cert = sum(1 for r in per_block if r["verdict"] == "CERTIFIED_POSITIVE")
    n_ref = sum(1 for r in per_block if r["verdict"] == "REFUTED")
    n_fail = sum(1 for r in per_block if r["verdict"] == "FAIL_TO_CERTIFY")
    verdict = ("REFUTED" if n_ref else
               "FAIL_TO_CERTIFY" if n_fail else "CERTIFIED_POSITIVE")
    return {
        "mode": "certify-pyiv",
        "engine": ("pure Python mpmath.iv directed interval arithmetic "
                   "(outward rounding), stdlib + mpmath only"),
        "certificate": cert.get("name", ""),
        "M": cert.get("M"),
        "K": cert.get("K"),
        "config": {"min_width": cfg.get("min_width", 1e-13),
                   "max_depth": cfg.get("max_depth", 80),
                   "max_boxes": cfg.get("max_boxes", 200000),
                   "iv_prec_bits": iv.prec},
        "n_blocks_total": len(cert["blocks"]),
        "n_blocks_checked": len(per_block),
        "blocks_checked_subset": (sorted(set(int(i) for i in block_filter))
                                  if block_filter is not None else None),
        "n_certified_positive": n_cert,
        "n_refuted": n_ref,
        "n_fail_to_certify": n_fail,
        "worst_certified_lower_bound": worst_lb if worst_i >= 0 else None,
        "worst_certified_block": worst_i,
        "total_boxes": sum(r["boxes_examined"] for r in per_block),
        "wall_seconds": time.time() - t0,
        "verdict": verdict,
        "per_block": per_block,
    }


# ---------------------------------------------------------------------------
# Cross-validation against an independently produced certify receipt
# ---------------------------------------------------------------------------
def cross_validate(py_out, ref_path):
    """Compare per-block verdicts (must match exactly) and certified lower
    bounds (order-of-magnitude agreement expected; exact equality is NOT
    expected — different engines take different bisection paths and stop at
    different, equally valid, positive bounds)."""
    with open(ref_path) as f:
        ref = json.load(f)
    ref_by_i = {int(r["i"]): r for r in ref["per_block"]}
    mismatches, ratios, spot = [], [], {}
    for r in py_out["per_block"]:
        i = r["i"]
        rr = ref_by_i.get(i)
        if rr is None:
            mismatches.append({"i": i, "reason": "missing in reference"})
            continue
        if rr["verdict"] != r["verdict"]:
            mismatches.append({"i": i, "py": r["verdict"],
                               "ref": rr["verdict"]})
            continue
        if (r["verdict"] == "CERTIFIED_POSITIVE"
                and rr["certified_lower_bound"] > 0
                and r["certified_lower_bound"] > 0):
            ratios.append((i, r["certified_lower_bound"]
                           / rr["certified_lower_bound"]))
    if ratios:
        rs = sorted(v for _, v in ratios)
        ratio_stats = {"n": len(rs), "min": rs[0], "max": rs[-1],
                       "median": rs[len(rs) // 2]}
    else:
        ratio_stats = None
    for i in (191, 196):
        if i in ref_by_i:
            mine = next((r for r in py_out["per_block"] if r["i"] == i), None)
            if mine:
                spot[str(i)] = {
                    "py_verdict": mine["verdict"],
                    "ref_verdict": ref_by_i[i]["verdict"],
                    "py_certified_lb": mine["certified_lower_bound"],
                    "ref_certified_lb": ref_by_i[i]["certified_lower_bound"],
                }
    return {
        "reference_file": ref_path,
        "reference_engine": "independent directed-interval certify receipt",
        "n_compared": len(py_out["per_block"]),
        "n_verdict_matches": len(py_out["per_block"]) - len(mismatches),
        "verdict_mismatches": mismatches,
        "certified_lb_ratio_py_over_ref": ratio_stats,
        "spot_checks": spot,
        "note": ("verdicts must match block-for-block; certified lower "
                 "bounds are path-dependent positive bounds and are only "
                 "expected to agree in order of magnitude"),
    }


# ---------------------------------------------------------------------------
# Deliberate-failure negative control
# ---------------------------------------------------------------------------
def self_test(cert, block_i, **cfg):
    """Perturb one weight of one block IN MEMORY (never on disk) and confirm
    the fail-closed path fires (REFUTED or FAIL_TO_CERTIFY). A verifier that
    cannot reject a corrupted certificate certifies nothing."""
    block = next(b for b in cert["blocks"] if int(b["i"]) == block_i)
    weights = [float(w) for w in block["weights"]]
    shifts = [float(d) for d in block_shifts(cert, block)]
    j = max(range(len(weights)), key=lambda k: abs(weights[k]))
    perturbed = list(weights)
    perturbed[j] = -perturbed[j]  # flip the dominant weight's sign

    clean = certify_block(float(block["A"]), float(block["C"]),
                          weights, shifts, block_i, **cfg)
    broken = certify_block(float(block["A"]), float(block["C"]),
                           perturbed, shifts, block_i, **cfg)
    fired = broken["verdict"] in ("REFUTED", "FAIL_TO_CERTIFY")
    return {
        "perturbed_in_memory": True,
        "block": block_i,
        "weight_index": j,
        "weight_before": weights[j],
        "weight_after": perturbed[j],
        "clean_verdict": clean["verdict"],
        "broken_verdict": broken["verdict"],
        "broken_refute_y": broken["refute_y"],
        "broken_refute_value_upper": broken["refute_value_upper"],
        "broken_fail_reason": broken["fail_reason"],
        "negative_control_passed": fired and
        clean["verdict"] == "CERTIFIED_POSITIVE",
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("json_path", help="sweep certificate JSON")
    ap.add_argument("--out", default=None, help="write receipt JSON here")
    ap.add_argument("--compare", default=None,
                    help="cross-validate against an existing certify receipt")
    ap.add_argument("--blocks", default=None,
                    help="comma-separated block ids (subset run)")
    ap.add_argument("--self-test", action="store_true",
                    help="also run the deliberate-perturbation negative control")
    ap.add_argument("--self-test-block", type=int, default=191)
    ap.add_argument("--progress-every", type=int, default=20)
    ap.add_argument("--progress-file", default=None,
                    help="append incremental progress lines here")
    ap.add_argument("--min-width", type=float, default=1e-13)
    ap.add_argument("--max-depth", type=int, default=80)
    ap.add_argument("--max-boxes", type=int, default=200000)
    args = ap.parse_args()

    with open(args.json_path) as f:
        cert = json.load(f)

    pf = open(args.progress_file, "a") if args.progress_file else None

    def report(done, total, last, worst_lb, worst_i, elapsed):
        line = (f"[{elapsed:8.1f}s] {done}/{total} blocks  last=i{last['i']} "
                f"{last['verdict']} lb={last['certified_lower_bound']:.3e}  "
                f"running-worst lb={worst_lb:.3e} @ block {worst_i}")
        print(line, flush=True)
        if pf:
            pf.write(line + "\n")
            pf.flush()

    cfg = {"min_width": args.min_width, "max_depth": args.max_depth,
           "max_boxes": args.max_boxes}
    block_filter = ([int(t) for t in args.blocks.split(",")]
                    if args.blocks else None)

    out = certify_certificate(cert, block_filter=block_filter,
                              progress=(report, args.progress_every), **cfg)
    out["python_version"] = sys.version.split()[0]
    import mpmath as _m
    out["mpmath_version"] = _m.__version__

    ok = out["verdict"] == "CERTIFIED_POSITIVE"

    if args.compare:
        out["cross_validation"] = cross_validate(out, args.compare)
        cv_ok = not out["cross_validation"]["verdict_mismatches"]
        ok = ok and cv_ok

    if args.self_test:
        out["self_test"] = self_test(cert, args.self_test_block, **cfg)
        ok = ok and out["self_test"]["negative_control_passed"]

    print()
    print("Certificate:", out["certificate"], " M =", out["M"], " K =", out["K"])
    print(f"Blocks: {out['n_certified_positive']} CERTIFIED_POSITIVE / "
          f"{out['n_refuted']} REFUTED / {out['n_fail_to_certify']} "
          f"FAIL_TO_CERTIFY of {out['n_blocks_checked']} checked "
          f"({out['n_blocks_total']} total)")
    print(f"Worst certified lower bound: {out['worst_certified_lower_bound']} "
          f"at block {out['worst_certified_block']}")
    print(f"Total boxes: {out['total_boxes']}  wall: {out['wall_seconds']:.1f}s")
    if args.compare:
        cv = out["cross_validation"]
        print(f"Cross-validation vs {cv['reference_file']}: "
              f"{cv['n_verdict_matches']}/{cv['n_compared']} verdicts match; "
              f"lb ratio (py/ref) stats: {cv['certified_lb_ratio_py_over_ref']}")
    if args.self_test:
        st = out["self_test"]
        print(f"Self-test (negative control): clean={st['clean_verdict']}, "
              f"perturbed={st['broken_verdict']} -> "
              f"{'PASS' if st['negative_control_passed'] else 'FAIL'}")
    print("VERDICT:", out["verdict"], "| overall:", "PASS" if ok else "FAIL")

    if args.out:
        with open(args.out, "w") as f:
            json.dump(out, f, indent=1)
        print("Receipt written:", args.out)
    if pf:
        pf.close()
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
