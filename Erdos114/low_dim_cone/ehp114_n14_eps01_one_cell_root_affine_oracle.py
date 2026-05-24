#!/usr/bin/env python3
"""Root-affine interval oracle prototype for one EHP114 n=14 cell.

This improves on the naive coefficient-box oracle by evaluating

    p(z) = prod_i (z - r_i)

directly with root intervals induced by the selected (u0,u1) coefficient cell.
That avoids the severe dependency blowup from interval elementary-symmetric
coefficient conversion.

The scope remains the current marching-squares grid oracle functional, not the
exact lemniscate length.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


sys.dont_write_bytecode = True

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-ORACLE-20260505-01"
PARENT_ID = "EXP-MATH-EHP114-N14-EPS01-ONE-CELL-VARIATION-PROBE-20260505-01"
DEGREE = 14
EPS = 0.1
RES = 220
EXTENT = 3.0
RADIAL_HALF = 12.0


def math_root_from_script() -> Path:
    return Path(__file__).resolve().parents[3]


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


def roots_of_unity(n: int) -> np.ndarray:
    j = np.arange(n, dtype=np.float64)
    return np.exp(2j * np.pi * j / n)


def eps_shape_scale(eps: float) -> float:
    return float(eps ** (1.0 / 28.0))


def scalar_target(eps: float) -> float:
    return float(RADIAL_HALF * eps ** (1.0 / DEGREE))


def root_intervals_from_cell(prev: Any, tensor: Any, u0: dict[str, Any], u1: dict[str, Any], cell: dict[str, Any]) -> list[Any]:
    base = (1.0 - EPS) ** (1.0 / DEGREE) * roots_of_unity(DEGREE)
    scale = eps_shape_scale(EPS)
    a = prev.interval_from_pair(cell["u0_interval"])
    b = prev.interval_from_pair(cell["u1_interval"])
    return [
        prev.complex_interval_linear(base[i], scale, a, u0["direction"][i], b, u1["direction"][i])
        for i in range(DEGREE)
    ]


def eval_poly_root_affine(prev: Any, z: complex, roots: list[Any]) -> Any:
    acc = prev.CI(prev.I(1.0, 1.0), prev.I(0.0, 0.0))
    zc = prev.fixed_complex(z)
    for root in roots:
        acc = acc * (zc - root)
    return acc


def interval_marching_upper(prev: Any, roots: list[Any]) -> dict[str, Any]:
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
            fsw = eval_poly_root_affine(prev, complex(x0, y0), roots).abs_sq_minus_one()
            fse = eval_poly_root_affine(prev, complex(x1, y0), roots).abs_sq_minus_one()
            fne = eval_poly_root_affine(prev, complex(x1, y1), roots).abs_sq_minus_one()
            fnw = eval_poly_root_affine(prev, complex(x0, y1), roots).abs_sq_minus_one()
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
            upper, avg_known = prev.cell_upper_from_fixed_case(case, (fsw, fse, fne, fnw), x0, y0, step)
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
    prev = load_module("ehp114_one_cell_interval_coeff_oracle_prev", out_dir / "ehp114_n14_eps01_one_cell_interval_coeff_oracle.py")
    low_dim = load_module("ehp114_low_dim_cone_box_probe_root_affine", out_dir / "ehp114_n14_low_dim_cone_box_probe.py")
    variation = load_json(out_dir / f"{PARENT_ID}_RESULTS.json")
    cell = variation["selected_cell"]
    tensor = load_module("ehp114_tensor_cone_scaffold_root_affine", root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py")
    _, directions = low_dim.build_low_dim_directions()
    u0 = directions[EPS][0]
    u1 = directions[EPS][1]
    roots = root_intervals_from_cell(prev, tensor, u0, u1, cell)
    oracle = interval_marching_upper(prev, roots)
    inari = load_json(root / "erdos-experiments" / "results" / "erdos-114" / "EXP-MM-EHP-007-n14-inari_RESULTS.json")
    lstar_lower = float(inari["l_star_lower"])
    target = scalar_target(EPS)
    deficit_lower = lstar_lower - oracle["length_upper"]
    margin_lower = deficit_lower - target
    status = "ONE_CELL_ROOT_AFFINE_ORACLE_PASS" if margin_lower >= 0.0 else "ONE_CELL_ROOT_AFFINE_ORACLE_FAIL"
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
        "interval_marching_upper": oracle,
        "lstar_lower": lstar_lower,
        "target": target,
        "deficit_lower": deficit_lower,
        "margin_lower": margin_lower,
        "certificate_scope": (
            "Root-affine interval upper bound for the current marching-squares oracle functional "
            "over one selected eps=0.1 coefficient rectangle."
        ),
        "not_a_full_ehp_certificate": (
            "This does not certify exact lemniscate length and does not prove Erdős #114. It is a prototype "
            "continuous-cell enclosure for the existing grid oracle."
        ),
        "next_blocker": (
            "If this passes, formalize the root-affine interval oracle and then replace the grid functional with "
            "a validated exact-length enclosure. If it fails, the remaining route is analytic Cauchy/Lipschitz."
        ),
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    oracle = result["interval_marching_upper"]
    lines = [
        "# EHP114 n=14 eps=0.1 One-Cell Root-Affine Oracle Prototype",
        "",
        f"Experiment: `{EXPERIMENT_ID}`",
        "",
        f"Parent packet: `{PARENT_ID}`",
        "",
        "## Meaning",
        "",
        "This run avoids coefficient-box dependency blowup by evaluating the polynomial",
        "directly as a product over root-affine intervals.",
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
