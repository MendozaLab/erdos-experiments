#!/usr/bin/env python3
"""Dense eps=0.1 two-dimensional box-candidate probe for EHP114 n=14.

This is the next step after the low-dimensional grid pass. It densifies the
two-dimensional spectral coefficient disk at eps=0.1 and summarizes grid cells.

Important: the current Rust oracle evaluates concrete points, not coefficient
boxes. This script therefore produces a box-candidate artifact, not a
continuous interval-box proof.
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

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-EPS01-LOW-DIM-BOX-CANDIDATE-20260505-01"
PARENT_ID = "EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01"
DEGREE = 14
EPS = 0.1
ETA = 0.014
GRID_N = 17
RES = 220
EXTENT = 3.0
RADIAL_HALF = 12.0
ROOT_TOL = 1e-12


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


def label_for(kind: str, a: float, b: float) -> str:
    return f"eps:1e-01:{kind}:u0:{a:+.8f}:u1:{b:+.8f}"


def run_length_oracle(root: Path, input_path: Path, output_path: Path) -> None:
    binary = root / "erdos-experiments" / "scripts" / "erdos-114" / "target" / "release" / "ehp114_batch_interval_lengths"
    if not binary.exists():
        raise SystemExit(f"Missing oracle binary: {binary}")
    subprocess.run(
        [str(binary), "--input", str(input_path), "--output", str(output_path), "--quiet"],
        cwd=root / "erdos-experiments" / "scripts" / "erdos-114",
        check=True,
    )


def build_eps01_directions(low_dim: Any) -> tuple[Any, Any]:
    _, directions = low_dim.build_low_dim_directions()
    return directions[EPS][0], directions[EPS][1]


def add_point(
    seen: set[str],
    points: list[dict[str, Any]],
    meta: dict[str, dict[str, Any]],
    tensor: Any,
    label: str,
    roots: np.ndarray,
    payload: dict[str, Any],
) -> None:
    radius = max_root_radius(roots)
    admissible = radius <= 1.0 + ROOT_TOL
    meta[label] = {**payload, "max_root_radius": radius, "admissible": admissible}
    if label in seen or not admissible:
        return
    seen.add(label)
    points.append({"label": label, "coeffs": coeffs_json(tensor.coeffs_from_roots(roots))})


def coefficient_grid() -> list[float]:
    return [float(x) for x in np.linspace(-ETA, ETA, GRID_N)]


def inside_disk(a: float, b: float, slack: float = 1e-15) -> bool:
    return a * a + b * b <= ETA * ETA + slack


def build_points(tensor: Any, u0: dict[str, Any], u1: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]], list[dict[str, Any]]]:
    roots0 = (1.0 - EPS) ** (1.0 / DEGREE) * tensor.roots_of_unity(DEGREE)
    scale = eps_shape_scale(EPS)
    grid = coefficient_grid()
    seen: set[str] = set()
    points: list[dict[str, Any]] = []
    meta: dict[str, dict[str, Any]] = {}
    cells: list[dict[str, Any]] = []

    for a in grid:
        for b in grid:
            if not inside_disk(a, b):
                continue
            label = label_for("node", a, b)
            roots = roots0 + scale * (a * u0["direction"] + b * u1["direction"])
            add_point(
                seen,
                points,
                meta,
                tensor,
                label,
                roots,
                {
                    "kind": "node",
                    "eps": EPS,
                    "u0_coefficient": a,
                    "u1_coefficient": b,
                    "eta_norm": float(np.hypot(a, b)),
                    "target": scalar_target(EPS),
                    "eps_shape_scale": scale,
                },
            )

    for i in range(len(grid) - 1):
        for j in range(len(grid) - 1):
            corners = [
                (grid[i], grid[j]),
                (grid[i + 1], grid[j]),
                (grid[i], grid[j + 1]),
                (grid[i + 1], grid[j + 1]),
            ]
            if not all(inside_disk(a, b) for a, b in corners):
                continue
            center = (0.5 * (grid[i] + grid[i + 1]), 0.5 * (grid[j] + grid[j + 1]))
            center_label = label_for("cell-center", center[0], center[1])
            roots = roots0 + scale * (center[0] * u0["direction"] + center[1] * u1["direction"])
            add_point(
                seen,
                points,
                meta,
                tensor,
                center_label,
                roots,
                {
                    "kind": "cell_center",
                    "eps": EPS,
                    "u0_coefficient": center[0],
                    "u1_coefficient": center[1],
                    "eta_norm": float(np.hypot(center[0], center[1])),
                    "target": scalar_target(EPS),
                    "eps_shape_scale": scale,
                },
            )
            cells.append(
                {
                    "i": i,
                    "j": j,
                    "u0_interval": [grid[i], grid[i + 1]],
                    "u1_interval": [grid[j], grid[j + 1]],
                    "corner_labels": [label_for("node", a, b) for a, b in corners],
                    "center_label": center_label,
                    "cell_radius": float(np.sqrt(2.0) * (grid[i + 1] - grid[i]) / 2.0),
                }
            )
    return points, meta, cells


def build_result() -> dict[str, Any]:
    root = math_root_from_script()
    out_dir = Path(__file__).resolve().parent
    erdos114_dir = erdos114_dir_from_script()
    packet = load_module("ehp114_n14_interval_taylor_m14_packet_box_candidate", erdos114_dir / "ehp114_n14_interval_taylor_m14_packet.py")
    tensor = load_module("ehp114_tensor_cone_scaffold_box_candidate", root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py")
    low_dim = load_module("ehp114_low_dim_cone_box_probe_parent", out_dir / "ehp114_n14_low_dim_cone_box_probe.py")
    u0, u1 = build_eps01_directions(low_dim)
    points, meta, cells = build_points(tensor, u0, u1)

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
                "eps": EPS,
                "eta_cap": ETA,
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
        evaluated = label in lengths
        deficit = None
        margin = None
        passed = False
        if evaluated:
            length = lengths[label]
            deficit = packet.Interval(lstar.lo - length.hi, lstar.hi - length.lo)
            margin = deficit.lo - row_meta["target"]
            passed = margin >= 0.0
        rows.append(
            {
                "label": label,
                **row_meta,
                "evaluated": evaluated,
                "deficit_interval": deficit.to_json() if deficit else None,
                "margin": margin,
                "pass": passed,
            }
        )
    row_by_label = {row["label"]: row for row in rows}
    evaluated_rows = [row for row in rows if row["evaluated"]]
    failures = [row for row in evaluated_rows if not row["pass"]]
    cell_rows = []
    for cell in cells:
        labels = cell["corner_labels"] + [cell["center_label"]]
        eval_rows = [row_by_label[label] for label in labels if label in row_by_label and row_by_label[label]["evaluated"]]
        cell_rows.append(
            {
                **cell,
                "evaluated_points": len(eval_rows),
                "all_points_evaluated": len(eval_rows) == 5,
                "all_points_pass": len(eval_rows) == 5 and all(row["pass"] for row in eval_rows),
                "min_sample_margin": min((row["margin"] for row in eval_rows if row["margin"] is not None), default=None),
                "max_sample_margin": max((row["margin"] for row in eval_rows if row["margin"] is not None), default=None),
            }
        )
    certified_candidate_cells = [row for row in cell_rows if row["all_points_pass"]]
    min_margin = min((row["margin"] for row in evaluated_rows if row["margin"] is not None), default=None)
    min_cell_margin = min((row["min_sample_margin"] for row in certified_candidate_cells if row["min_sample_margin"] is not None), default=None)
    status = "EPS01_LOW_DIM_BOX_CANDIDATE_PASS" if evaluated_rows and not failures else "EPS01_LOW_DIM_BOX_CANDIDATE_FAIL"
    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "parent_packet": PARENT_ID,
        "status": status,
        "parameters": {
            "degree": DEGREE,
            "eps": EPS,
            "eta_cap": ETA,
            "grid_n": GRID_N,
            "res": RES,
            "extent": EXTENT,
            "oracle_point_count": len(points),
            "candidate_cell_count": len(cells),
            "scalar_target": "D14(eps,s) >= 12 * eps^(1/14)",
            "subspace": "span of the two lowest midpoint eigendirections from the eps=0.1 admissible spectral Taylor matrix",
        },
        "u0_direction": {
            "rank": int(u0["rank"]),
            "midpoint_eigenvalue": float(u0["midpoint_eigenvalue"]),
            "dominant_shape_components": u0["dominant_shape_components"],
        },
        "u1_direction": {
            "rank": int(u1["rank"]),
            "midpoint_eigenvalue": float(u1["midpoint_eigenvalue"]),
            "dominant_shape_components": u1["dominant_shape_components"],
        },
        "evaluated_point_count": len(evaluated_rows),
        "failure_count": len(failures),
        "min_point_margin": min_margin,
        "candidate_cells_all_five_points_pass": len(certified_candidate_cells),
        "min_candidate_cell_sample_margin": min_cell_margin,
        "continuous_certificate": False,
        "why_not_continuous_certificate": (
            "The current oracle evaluates concrete coefficient points. It does not enclose the full "
            "(u0,u1) coefficient boxes. A rigorous box certificate still needs an interval lower bound "
            "or a Lipschitz/Cauchy variation bound over each cell."
        ),
        "rows": rows,
        "cells": cell_rows,
        "oracle": {
            "input": str(input_path),
            "output": str(output_path),
            "output_point_count": len(output["points"]),
        },
        "claim_ceiling": (
            "Dense eps=0.1 2D spectral-span box-candidate evidence. This is not a continuous "
            "interval-box certificate, not a Lean proof, and not a proof of Erdős #114."
        ),
        "next_blocker": (
            "Add a rigorous cell-wise variation bound for D14 over (u0,u1) coefficient boxes, "
            "or upgrade the Rust oracle to accept interval coefficient boxes."
        ),
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    lines = [
        "# EHP114 n=14 eps=0.1 Two-Dimensional Box-Candidate Probe",
        "",
        f"Experiment: `{EXPERIMENT_ID}`",
        "",
        f"Parent packet: `{PARENT_ID}`",
        "",
        "## Meaning",
        "",
        "This run densifies the eps=0.1 two-dimensional spectral span where the previous",
        "low-dimensional probe had full root-admissibility. It is a bridge toward an",
        "interval-box proof, but it is still point-evaluation evidence.",
        "",
        "## Verdict",
        "",
        f"- Status: `{result['status']}`",
        f"- Evaluated oracle points: `{result['evaluated_point_count']}`",
        f"- Failure count: `{result['failure_count']}`",
        f"- Minimum point margin: `{result['min_point_margin']}`",
        f"- Candidate cells whose four corners and center all pass: `{result['candidate_cells_all_five_points_pass']}`",
        f"- Minimum sample margin among candidate cells: `{result['min_candidate_cell_sample_margin']}`",
        f"- Continuous certificate: `{result['continuous_certificate']}`",
        "",
        "## Why This Is Not Yet A Box Certificate",
        "",
        result["why_not_continuous_certificate"],
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
                "evaluated_point_count": result["evaluated_point_count"],
                "failure_count": result["failure_count"],
                "min_point_margin": result["min_point_margin"],
                "candidate_cells_all_five_points_pass": result["candidate_cells_all_five_points_pass"],
                "report": str(report_path),
                "sha256": sha,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
