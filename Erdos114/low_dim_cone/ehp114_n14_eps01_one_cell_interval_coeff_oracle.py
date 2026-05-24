#!/usr/bin/env python3
"""Interval-coefficient oracle prototype for one EHP114 n=14 cell.

This is a deliberately narrow certificate attempt for the selected eps=0.1
low-dimensional spectral cell. It encloses the polynomial coefficients induced
by the full coefficient rectangle, then runs a conservative interval
marching-squares upper bound for the oracle length functional.

The bound is for the current grid/oracle functional, not a formal theorem about
the exact lemniscate length.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


sys.dont_write_bytecode = True

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-EPS01-ONE-CELL-INTERVAL-COEFF-ORACLE-20260505-01"
PARENT_ID = "EXP-MATH-EHP114-N14-EPS01-ONE-CELL-VARIATION-PROBE-20260505-01"
DEGREE = 14
EPS = 0.1
RES = 220
EXTENT = 3.0
RADIAL_HALF = 12.0


@dataclass(frozen=True)
class I:
    lo: float
    hi: float

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError((self.lo, self.hi))

    def __add__(self, other: "I") -> "I":
        return I(self.lo + other.lo, self.hi + other.hi)

    def __sub__(self, other: "I") -> "I":
        return I(self.lo - other.hi, self.hi - other.lo)

    def __neg__(self) -> "I":
        return I(-self.hi, -self.lo)

    def __mul__(self, other: "I") -> "I":
        vals = [self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi]
        return I(min(vals), max(vals))

    def scale(self, c: float) -> "I":
        vals = [self.lo * c, self.hi * c]
        return I(min(vals), max(vals))

    def div(self, other: "I") -> "I":
        if other.lo <= 0.0 <= other.hi:
            raise ZeroDivisionError(other)
        vals = [self.lo / other.lo, self.lo / other.hi, self.hi / other.lo, self.hi / other.hi]
        return I(min(vals), max(vals))

    def square(self) -> "I":
        if self.lo <= 0.0 <= self.hi:
            return I(0.0, max(self.lo * self.lo, self.hi * self.hi))
        vals = [self.lo * self.lo, self.hi * self.hi]
        return I(min(vals), max(vals))

    def sqrt(self) -> "I":
        return I(math.sqrt(max(0.0, self.lo)), math.sqrt(max(0.0, self.hi)))

    def sign(self) -> str:
        if self.lo > 0.0:
            return "+"
        if self.hi < 0.0:
            return "-"
        return "?"

    def to_json(self) -> dict[str, float]:
        return {"lo": self.lo, "hi": self.hi, "width": self.hi - self.lo}


@dataclass(frozen=True)
class CI:
    re: I
    im: I

    def __add__(self, other: "CI") -> "CI":
        return CI(self.re + other.re, self.im + other.im)

    def __sub__(self, other: "CI") -> "CI":
        return CI(self.re - other.re, self.im - other.im)

    def __neg__(self) -> "CI":
        return CI(-self.re, -self.im)

    def __mul__(self, other: "CI") -> "CI":
        return CI(self.re * other.re - self.im * other.im, self.re * other.im + self.im * other.re)

    def abs_sq_minus_one(self) -> I:
        return self.re.square() + self.im.square() - I(1.0, 1.0)

    def to_json(self) -> dict[str, Any]:
        return {"re": self.re.to_json(), "im": self.im.to_json()}


def math_root_from_script() -> Path:
    return Path(__file__).resolve().parents[3]


def erdos114_dir_from_script() -> Path:
    return Path(__file__).resolve().parents[1]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def interval_const(x: float) -> I:
    return I(float(x), float(x))


def interval_from_pair(pair: list[float]) -> I:
    return I(float(pair[0]), float(pair[1]))


def eps_shape_scale(eps: float) -> float:
    return float(eps ** (1.0 / 28.0))


def scalar_target(eps: float) -> float:
    return float(RADIAL_HALF * eps ** (1.0 / DEGREE))


def fixed_complex(z: complex) -> CI:
    return CI(interval_const(z.real), interval_const(z.imag))


def roots_of_unity(n: int) -> np.ndarray:
    j = np.arange(n, dtype=np.float64)
    return np.exp(2j * np.pi * j / n)


def complex_interval_linear(
    base: complex,
    scale: float,
    a: I,
    u0: complex,
    b: I,
    u1: complex,
) -> CI:
    re = interval_const(base.real) + (a.scale(u0.real) + b.scale(u1.real)).scale(scale)
    im = interval_const(base.imag) + (a.scale(u0.imag) + b.scale(u1.imag)).scale(scale)
    return CI(re, im)


def coeff_intervals_from_root_intervals(roots: list[CI]) -> list[CI]:
    zero = CI(I(0.0, 0.0), I(0.0, 0.0))
    one = CI(I(1.0, 1.0), I(0.0, 0.0))
    desc = [one]
    for root in roots:
        new = [zero for _ in range(len(desc) + 1)]
        for i, coeff in enumerate(desc):
            new[i] = new[i] + coeff
            new[i + 1] = new[i + 1] + coeff * (-root)
        desc = new
    return list(reversed(desc[1:]))


def eval_poly_interval(z: complex, coeffs_asc: list[CI]) -> CI:
    acc = CI(I(1.0, 1.0), I(0.0, 0.0))
    zc = fixed_complex(z)
    for k in range(len(coeffs_asc) - 1, -1, -1):
        acc = acc * zc + coeffs_asc[k]
    return acc


def interp_interval(fa: I, fb: I) -> I:
    try:
        out = fa.div(fa - fb)
        lo = max(0.0, out.lo)
        hi = min(1.0, out.hi)
        if lo > hi:
            return I(0.0, 1.0)
        return I(lo, hi)
    except ZeroDivisionError:
        return I(0.0, 1.0)


def seg_upper(ax: I, ay: I, bx: I, by: I) -> float:
    dx = ax - bx
    dy = ay - by
    return (dx.square() + dy.square()).sqrt().hi


def cell_upper_from_fixed_case(case: int, vals: tuple[I, I, I, I], x0: float, y0: float, step: float) -> tuple[float, bool]:
    fsw, fse, fne, fnw = vals
    x0i = I(x0, x0)
    y0i = I(y0, y0)
    x1i = I(x0 + step, x0 + step)
    y1i = I(y0 + step, y0 + step)
    stepi = I(step, step)
    s = (x0i + interp_interval(fsw, fse) * stepi, y0i)
    e = (x1i, y0i + interp_interval(fse, fne) * stepi)
    n = (x0i + interp_interval(fnw, fne) * stepi, y1i)
    w = (x0i, y0i + interp_interval(fsw, fnw) * stepi)

    def seg(p: tuple[I, I], q: tuple[I, I]) -> float:
        return seg_upper(p[0], p[1], q[0], q[1])

    if case in (1, 14):
        return seg(s, w), True
    if case in (2, 13):
        return seg(s, e), True
    if case in (3, 12):
        return seg(w, e), True
    if case in (4, 11):
        return seg(e, n), True
    if case in (6, 9):
        return seg(s, n), True
    if case in (7, 8):
        return seg(w, n), True
    if case == 5:
        avg = (fsw + fse + fne + fnw).scale(0.25)
        if avg.sign() == "+":
            return seg(s, w) + seg(e, n), True
        if avg.sign() == "-":
            return seg(s, e) + seg(w, n), True
        return 2.0 * math.sqrt(2.0) * step, False
    if case == 10:
        avg = (fsw + fse + fne + fnw).scale(0.25)
        if avg.sign() == "+":
            return seg(s, e) + seg(w, n), True
        if avg.sign() == "-":
            return seg(s, w) + seg(e, n), True
        return 2.0 * math.sqrt(2.0) * step, False
    return 0.0, True


def interval_marching_upper(coeffs_asc: list[CI]) -> dict[str, Any]:
    step = 2.0 * EXTENT / RES
    xs = np.linspace(-EXTENT, EXTENT, RES + 1)
    ys = np.linspace(-EXTENT, EXTENT, RES + 1)
    total_upper = 0.0
    active_cells = 0
    definite_case_cells = 0
    uncertain_corner_cells = 0
    ambiguous_avg_cells = 0
    crude_penalty_upper = 2.0 * math.sqrt(2.0) * step
    max_cell_upper = 0.0

    for ix in range(RES):
        x0 = float(xs[ix])
        x1 = float(xs[ix + 1])
        for iy in range(RES):
            y0 = float(ys[iy])
            y1 = float(ys[iy + 1])
            fsw = eval_poly_interval(complex(x0, y0), coeffs_asc).abs_sq_minus_one()
            fse = eval_poly_interval(complex(x1, y0), coeffs_asc).abs_sq_minus_one()
            fne = eval_poly_interval(complex(x1, y1), coeffs_asc).abs_sq_minus_one()
            fnw = eval_poly_interval(complex(x0, y1), coeffs_asc).abs_sq_minus_one()
            signs = [fsw.sign(), fse.sign(), fne.sign(), fnw.sign()]
            if all(s == "+" for s in signs) or all(s == "-" for s in signs):
                continue
            active_cells += 1
            if "?" in signs:
                uncertain_corner_cells += 1
                upper = crude_penalty_upper
                total_upper += upper
                max_cell_upper = max(max_cell_upper, upper)
                continue
            case = (1 if signs[0] == "+" else 0) | ((1 if signs[1] == "+" else 0) << 1) | ((1 if signs[2] == "+" else 0) << 2) | ((1 if signs[3] == "+" else 0) << 3)
            upper, avg_known = cell_upper_from_fixed_case(case, (fsw, fse, fne, fnw), x0, y0, step)
            definite_case_cells += 1
            if not avg_known:
                ambiguous_avg_cells += 1
            total_upper += upper
            max_cell_upper = max(max_cell_upper, upper)
    return {
        "length_upper": total_upper,
        "active_cells": active_cells,
        "definite_case_cells": definite_case_cells,
        "uncertain_corner_cells": uncertain_corner_cells,
        "ambiguous_avg_cells": ambiguous_avg_cells,
        "crude_penalty_upper_per_uncertain_cell": crude_penalty_upper,
        "max_cell_upper": max_cell_upper,
    }


def build_result() -> dict[str, Any]:
    root = math_root_from_script()
    out_dir = Path(__file__).resolve().parent
    erdos114_dir = erdos114_dir_from_script()
    variation = load_json(out_dir / f"{PARENT_ID}_RESULTS.json")
    cell = variation["selected_cell"]
    tensor = load_module("ehp114_tensor_cone_scaffold_interval_coeff", root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py")
    low_dim = load_module("ehp114_low_dim_cone_box_probe_interval_coeff", out_dir / "ehp114_n14_low_dim_cone_box_probe.py")
    _, directions = low_dim.build_low_dim_directions()
    u0 = directions[EPS][0]
    u1 = directions[EPS][1]
    base = (1.0 - EPS) ** (1.0 / DEGREE) * roots_of_unity(DEGREE)
    scale = eps_shape_scale(EPS)
    a = interval_from_pair(cell["u0_interval"])
    b = interval_from_pair(cell["u1_interval"])
    roots = [
        complex_interval_linear(base[i], scale, a, u0["direction"][i], b, u1["direction"][i])
        for i in range(DEGREE)
    ]
    coeffs = coeff_intervals_from_root_intervals(roots)
    oracle = interval_marching_upper(coeffs)
    inari = load_json(root / "erdos-experiments" / "results" / "erdos-114" / "EXP-MM-EHP-007-n14-inari_RESULTS.json")
    lstar_lower = float(inari["l_star_lower"])
    target = scalar_target(EPS)
    deficit_lower = lstar_lower - oracle["length_upper"]
    margin_lower = deficit_lower - target
    status = "ONE_CELL_INTERVAL_COEFF_ORACLE_PASS" if margin_lower >= 0.0 else "ONE_CELL_INTERVAL_COEFF_ORACLE_FAIL"
    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "parent_packet": PARENT_ID,
        "status": status,
        "parameters": {
            "degree": DEGREE,
            "eps": EPS,
            "res": RES,
            "extent": EXTENT,
            "scalar_target": "D14(eps,s) >= 12 * eps^(1/14)",
        },
        "selected_cell": cell,
        "root_box_from_parent": variation["root_box"],
        "coefficient_intervals": [c.to_json() for c in coeffs],
        "interval_marching_upper": oracle,
        "lstar_lower": lstar_lower,
        "target": target,
        "deficit_lower": deficit_lower,
        "margin_lower": margin_lower,
        "certificate_scope": (
            "Conservative interval-coefficient upper bound for the current marching-squares oracle functional "
            "over one selected eps=0.1 coefficient rectangle."
        ),
        "not_a_full_ehp_certificate": (
            "This does not certify exact lemniscate length and does not prove Erdős #114. It is a prototype "
            "continuous-cell enclosure for the existing grid oracle."
        ),
        "next_blocker": (
            "Replace the grid-oracle functional by an exact/validated lemniscate-length enclosure, or formalize "
            "why this marching-squares interval oracle is an accepted certificate layer."
        ),
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    oracle = result["interval_marching_upper"]
    lines = [
        "# EHP114 n=14 eps=0.1 One-Cell Interval-Coefficient Oracle Prototype",
        "",
        f"Experiment: `{EXPERIMENT_ID}`",
        "",
        f"Parent packet: `{PARENT_ID}`",
        "",
        "## Meaning",
        "",
        "This run tries the first continuous-cell enclosure for the selected eps=0.1",
        "coefficient rectangle. It propagates the full coefficient box through a",
        "conservative interval marching-squares upper bound.",
        "",
        "## Verdict",
        "",
        f"- Status: `{result['status']}`",
        f"- Length upper bound: `{oracle['length_upper']}`",
        f"- Deficit lower bound: `{result['deficit_lower']}`",
        f"- Target: `{result['target']}`",
        f"- Margin lower bound: `{result['margin_lower']}`",
        f"- Active cells: `{oracle['active_cells']}`",
        f"- Definite-case cells: `{oracle['definite_case_cells']}`",
        f"- Uncertain-corner cells: `{oracle['uncertain_corner_cells']}`",
        f"- Ambiguous-average cells: `{oracle['ambiguous_avg_cells']}`",
        "",
        "## Certificate Scope",
        "",
        result["certificate_scope"],
        "",
        "## Claim Ceiling",
        "",
        result["not_a_full_ehp_certificate"],
        "",
        "## Next Blocker",
        "",
        result["next_blocker"],
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    out_dir = Path(__file__).resolve().parent
    result_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {path}")
    result = build_result()
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(result, report_path)
    sha = sha256_file(result_path)
    sha_path.write_text(f"{sha}  {result_path.name}\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "experiment_id": EXPERIMENT_ID,
                "status": result["status"],
                "length_upper": result["interval_marching_upper"]["length_upper"],
                "deficit_lower": result["deficit_lower"],
                "margin_lower": result["margin_lower"],
                "report": str(report_path),
                "sha256": sha,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
