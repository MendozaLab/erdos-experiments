#!/usr/bin/env python3
"""One-cell variation probe for the EHP114 n=14 eps=0.1 low-dimensional lane.

This script selects the strongest cell from the eps=0.1 dense box-candidate
packet, proves a floating affine root-radius enclosure for that coefficient
rectangle, and densely samples the scalar deficit margin inside the cell.

It does not certify the deficit continuously. The missing proof ingredient is a
rigorous Lipschitz/Cauchy variation bound or an interval-coefficient oracle.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


sys.dont_write_bytecode = True

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-EPS01-ONE-CELL-VARIATION-PROBE-20260505-01"
PARENT_ID = "EXP-MATH-EHP114-N14-EPS01-LOW-DIM-BOX-CANDIDATE-20260505-01"
DEGREE = 14
EPS = 0.1
GRID_N = 17
RES = 220
EXTENT = 3.0
RADIAL_HALF = 12.0
ROOT_FLOAT_GUARD = 1e-12


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


def coeffs_json(coeffs: np.ndarray) -> list[list[float]]:
    return [[float(z.real), float(z.imag)] for z in coeffs]


def max_root_radius(roots: np.ndarray) -> float:
    return float(np.max(np.abs(roots)))


def eps_shape_scale(eps: float) -> float:
    return float(eps ** (1.0 / 28.0))


def scalar_target(eps: float) -> float:
    return float(RADIAL_HALF * eps ** (1.0 / DEGREE))


def label_for(a: float, b: float) -> str:
    return f"eps:1e-01:one-cell:u0:{a:+.10f}:u1:{b:+.10f}"


def run_length_oracle(root: Path, input_path: Path, output_path: Path) -> None:
    binary = root / "erdos-experiments" / "scripts" / "erdos-114" / "target" / "release" / "ehp114_batch_interval_lengths"
    if not binary.exists():
        raise SystemExit(f"Missing oracle binary: {binary}")
    subprocess.run(
        [str(binary), "--input", str(input_path), "--output", str(output_path), "--quiet"],
        cwd=root / "erdos-experiments" / "scripts" / "erdos-114",
        check=True,
    )


def select_best_cell(parent: dict[str, Any]) -> dict[str, Any]:
    cells = [cell for cell in parent["cells"] if cell["all_points_pass"] and cell["min_sample_margin"] is not None]
    if not cells:
        raise SystemExit("No passing cells found in parent packet")
    return max(cells, key=lambda cell: cell["min_sample_margin"])


def build_eps01_directions(low_dim: Any) -> tuple[Any, Any]:
    _, directions = low_dim.build_low_dim_directions()
    return directions[EPS][0], directions[EPS][1]


def affine_root_radius_bound(
    tensor: Any,
    u0: dict[str, Any],
    u1: dict[str, Any],
    u0_interval: list[float],
    u1_interval: list[float],
) -> dict[str, Any]:
    base = (1.0 - EPS) ** (1.0 / DEGREE) * tensor.roots_of_unity(DEGREE)
    scale = eps_shape_scale(EPS)
    a_mid = 0.5 * (u0_interval[0] + u0_interval[1])
    b_mid = 0.5 * (u1_interval[0] + u1_interval[1])
    a_rad = 0.5 * (u0_interval[1] - u0_interval[0])
    b_rad = 0.5 * (u1_interval[1] - u1_interval[0])
    center_roots = base + scale * (a_mid * u0["direction"] + b_mid * u1["direction"])
    per_root = []
    for idx, z in enumerate(center_roots):
        variation = scale * (a_rad * abs(u0["direction"][idx]) + b_rad * abs(u1["direction"][idx]))
        upper = abs(z) + variation + ROOT_FLOAT_GUARD
        per_root.append(
            {
                "root_index": idx,
                "center_abs": float(abs(z)),
                "variation_radius": float(variation),
                "upper_radius_bound": float(upper),
            }
        )
    max_upper = max(row["upper_radius_bound"] for row in per_root)
    return {
        "u0_mid": float(a_mid),
        "u1_mid": float(b_mid),
        "u0_radius": float(a_rad),
        "u1_radius": float(b_rad),
        "max_upper_radius_bound": float(max_upper),
        "root_box_admissible_by_affine_bound": bool(max_upper <= 1.0),
        "per_root_bounds": per_root,
    }


def build_points(
    tensor: Any,
    u0: dict[str, Any],
    u1: dict[str, Any],
    cell: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    base = (1.0 - EPS) ** (1.0 / DEGREE) * tensor.roots_of_unity(DEGREE)
    scale = eps_shape_scale(EPS)
    u0_values = np.linspace(cell["u0_interval"][0], cell["u0_interval"][1], GRID_N)
    u1_values = np.linspace(cell["u1_interval"][0], cell["u1_interval"][1], GRID_N)
    points: list[dict[str, Any]] = []
    meta: dict[str, dict[str, Any]] = {}
    for a in u0_values:
        for b in u1_values:
            roots = base + scale * (float(a) * u0["direction"] + float(b) * u1["direction"])
            label = label_for(float(a), float(b))
            points.append({"label": label, "coeffs": coeffs_json(tensor.coeffs_from_roots(roots))})
            meta[label] = {
                "eps": EPS,
                "u0_coefficient": float(a),
                "u1_coefficient": float(b),
                "target": scalar_target(EPS),
                "max_root_radius": max_root_radius(roots),
            }
    return points, meta


def finite_difference_lipschitz(rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_coord = {(row["u0_coefficient"], row["u1_coefficient"]): row for row in rows}
    coords0 = sorted(set(row["u0_coefficient"] for row in rows))
    coords1 = sorted(set(row["u1_coefficient"] for row in rows))
    slopes = []
    for i in range(len(coords0) - 1):
        for b in coords1:
            r0 = by_coord[(coords0[i], b)]
            r1 = by_coord[(coords0[i + 1], b)]
            dist = abs(coords0[i + 1] - coords0[i])
            slopes.append(abs(r1["margin"] - r0["margin"]) / dist)
    for a in coords0:
        for j in range(len(coords1) - 1):
            r0 = by_coord[(a, coords1[j])]
            r1 = by_coord[(a, coords1[j + 1])]
            dist = abs(coords1[j + 1] - coords1[j])
            slopes.append(abs(r1["margin"] - r0["margin"]) / dist)
    return {
        "max_neighbor_margin_slope": max(slopes) if slopes else None,
        "median_neighbor_margin_slope": float(np.median(slopes)) if slopes else None,
        "slope_count": len(slopes),
    }


def build_result() -> dict[str, Any]:
    root = math_root_from_script()
    out_dir = Path(__file__).resolve().parent
    erdos114_dir = erdos114_dir_from_script()
    parent = load_json(out_dir / f"{PARENT_ID}_RESULTS.json")
    cell = select_best_cell(parent)
    packet = load_module("ehp114_n14_interval_taylor_m14_packet_one_cell", erdos114_dir / "ehp114_n14_interval_taylor_m14_packet.py")
    tensor = load_module("ehp114_tensor_cone_scaffold_one_cell", root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py")
    low_dim = load_module("ehp114_low_dim_cone_box_probe_one_cell", out_dir / "ehp114_n14_low_dim_cone_box_probe.py")
    u0, u1 = build_eps01_directions(low_dim)
    root_box = affine_root_radius_bound(tensor, u0, u1, cell["u0_interval"], cell["u1_interval"])
    points, meta = build_points(tensor, u0, u1, cell)

    input_path = out_dir / f"{EXPERIMENT_ID}_ORACLE_INPUT.json"
    output_path = out_dir / f"{EXPERIMENT_ID}_ORACLE_OUTPUT.json"
    for path in (input_path, output_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing oracle artifact: {path}")
    input_path.write_text(
        json.dumps(
            {
                "degree": DEGREE,
                "res": RES,
                "extent": EXTENT,
                "parent_packet": PARENT_ID,
                "selected_cell": {
                    "i": cell["i"],
                    "j": cell["j"],
                    "u0_interval": cell["u0_interval"],
                    "u1_interval": cell["u1_interval"],
                    "parent_min_sample_margin": cell["min_sample_margin"],
                },
                "grid_n": GRID_N,
                "point_count": len(points),
                "points": points,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    run_length_oracle(root, input_path, output_path)

    output = load_json(output_path)
    lengths = {
        row["label"]: packet.Interval(float(row["length_lower"]), float(row["length_upper"]))
        for row in output["points"]
    }
    inari = load_json(root / "erdos-experiments" / "results" / "erdos-114" / "EXP-MM-EHP-007-n14-inari_RESULTS.json")
    lstar = packet.Interval(float(inari["l_star_lower"]), float(inari["l_star_upper"]))
    rows = []
    for label, row_meta in sorted(meta.items()):
        length = lengths[label]
        deficit = packet.Interval(lstar.lo - length.hi, lstar.hi - length.lo)
        margin = deficit.lo - row_meta["target"]
        rows.append(
            {
                "label": label,
                **row_meta,
                "deficit_interval": deficit.to_json(),
                "margin": margin,
                "pass": margin >= 0.0,
            }
        )
    failures = [row for row in rows if not row["pass"]]
    min_margin = min(row["margin"] for row in rows)
    max_margin = max(row["margin"] for row in rows)
    slope = finite_difference_lipschitz(rows)
    cell_radius = float(np.sqrt(2.0) * 0.5 * (cell["u0_interval"][1] - cell["u0_interval"][0]))
    required_lipschitz_for_margin = min_margin / cell_radius
    empirical_margin_after_cell_radius = (
        min_margin - slope["max_neighbor_margin_slope"] * cell_radius
        if slope["max_neighbor_margin_slope"] is not None
        else None
    )
    status = (
        "ONE_CELL_ROOT_BOX_PASS_EMPIRICAL_VARIATION_PASS"
        if root_box["root_box_admissible_by_affine_bound"] and not failures
        else "ONE_CELL_VARIATION_PROBE_FAIL"
    )
    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "parent_packet": PARENT_ID,
        "status": status,
        "parameters": {
            "degree": DEGREE,
            "eps": EPS,
            "grid_n": GRID_N,
            "res": RES,
            "extent": EXTENT,
            "oracle_point_count": len(points),
            "scalar_target": "D14(eps,s) >= 12 * eps^(1/14)",
        },
        "selected_cell": cell,
        "root_box": root_box,
        "failure_count": len(failures),
        "min_margin": min_margin,
        "max_margin": max_margin,
        "cell_radius": cell_radius,
        "required_lipschitz_for_margin": required_lipschitz_for_margin,
        "finite_difference_lipschitz_diagnostic": slope,
        "empirical_margin_after_cell_radius": empirical_margin_after_cell_radius,
        "continuous_deficit_certificate": False,
        "why_not_continuous_deficit_certificate": (
            "The affine root-radius bound encloses root admissibility for the selected coefficient cell. "
            "The deficit lower bound is still sampled at grid points. Certifying it continuously requires "
            "a rigorous Lipschitz/Cauchy variation bound or an interval-coefficient length oracle."
        ),
        "rows": rows,
        "oracle": {
            "input": str(input_path),
            "output": str(output_path),
            "output_point_count": len(output["points"]),
        },
        "claim_ceiling": (
            "One-cell root-box enclosure plus dense deficit variation probe at eps=0.1. "
            "This is not a continuous deficit certificate, not a Lean proof, and not a proof of Erdős #114."
        ),
        "next_blocker": (
            "Turn the empirical margin-variation diagnostic into a rigorous Lipschitz/Cauchy bound over the selected cell."
        ),
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    root_box = result["root_box"]
    slope = result["finite_difference_lipschitz_diagnostic"]
    lines = [
        "# EHP114 n=14 eps=0.1 One-Cell Variation Probe",
        "",
        f"Experiment: `{EXPERIMENT_ID}`",
        "",
        f"Parent packet: `{PARENT_ID}`",
        "",
        "## Meaning",
        "",
        "This run selects the strongest sampled eps=0.1 coefficient cell and tests the",
        "first bridge from point evidence toward a cell certificate.",
        "",
        "## Verdict",
        "",
        f"- Status: `{result['status']}`",
        f"- Root box admissible by affine radius bound: `{root_box['root_box_admissible_by_affine_bound']}`",
        f"- Max affine root-radius upper bound: `{root_box['max_upper_radius_bound']}`",
        f"- Oracle grid points inside selected cell: `{result['parameters']['oracle_point_count']}`",
        f"- Deficit failures on sampled points: `{result['failure_count']}`",
        f"- Minimum sampled margin: `{result['min_margin']}`",
        f"- Maximum sampled margin: `{result['max_margin']}`",
        f"- Max neighbor margin slope: `{slope['max_neighbor_margin_slope']}`",
        f"- Required Lipschitz upper bound for this margin: `{result['required_lipschitz_for_margin']}`",
        f"- Empirical margin after one cell-radius variation: `{result['empirical_margin_after_cell_radius']}`",
        f"- Continuous deficit certificate: `{result['continuous_deficit_certificate']}`",
        "",
        "## Why This Is Not Yet A Deficit Certificate",
        "",
        result["why_not_continuous_deficit_certificate"],
        "",
        "## Next Blocker",
        "",
        result["next_blocker"],
        "",
        "## Claim Ceiling",
        "",
        result["claim_ceiling"],
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
                "root_box_admissible": result["root_box"]["root_box_admissible_by_affine_bound"],
                "min_margin": result["min_margin"],
                "empirical_margin_after_cell_radius": result["empirical_margin_after_cell_radius"],
                "report": str(report_path),
                "sha256": sha,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
