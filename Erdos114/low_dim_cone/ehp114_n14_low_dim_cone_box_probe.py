#!/usr/bin/env python3
"""Finite low-dimensional cone probe for EHP114 n=14.

This script samples the two-dimensional coefficient disk spanned by the two
lowest midpoint eigendirections from the admissible spectral Taylor packet.
It is a finite probe, not a continuous cone proof.
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

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01"
PARENT_ID = "EXP-MATH-EHP114-N14-EPS-SCALED-SPECTRAL-DIRECTION-SEARCH-20260505-01"
DEGREE = 14
EPS_VALUES = [1e-4, 1e-3, 1e-2, 1e-1]
ETA = 0.014
COEFFICIENT_GRID = [-ETA, -0.5 * ETA, 0.0, 0.5 * ETA, ETA]
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


def run_length_oracle(root: Path, input_path: Path, output_path: Path) -> None:
    binary = root / "erdos-experiments" / "scripts" / "erdos-114" / "target" / "release" / "ehp114_batch_interval_lengths"
    if not binary.exists():
        raise SystemExit(f"Missing oracle binary; refusing to build outside low_dim_cone: {binary}")
    subprocess.run(
        [str(binary), "--input", str(input_path), "--output", str(output_path), "--quiet"],
        cwd=root / "erdos-experiments" / "scripts" / "erdos-114",
        check=True,
    )


def build_low_dim_directions() -> tuple[list[Any], dict[float, list[dict[str, Any]]]]:
    root = math_root_from_script()
    erdos114_dir = erdos114_dir_from_script()
    packet = load_module("ehp114_n14_interval_taylor_m14_packet", erdos114_dir / "ehp114_n14_interval_taylor_m14_packet.py")
    tensor = load_module("ehp114_tensor_cone_scaffold", root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py")
    spectral = load_module(
        "ehp114_n14_eps_scaled_spectral_direction_search",
        erdos114_dir / "ehp114_n14_eps_scaled_spectral_direction_search.py",
    )
    shape, directions = spectral.eigen_directions(packet, tensor, root)
    return shape, directions


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
    if not admissible or label in seen:
        return
    seen.add(label)
    points.append({"label": label, "coeffs": coeffs_json(tensor.coeffs_from_roots(roots))})


def build_points(
    tensor: Any,
    directions: dict[float, list[dict[str, Any]]],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]], list[dict[str, Any]]]:
    unit_roots = tensor.roots_of_unity(DEGREE)
    seen: set[str] = set()
    points: list[dict[str, Any]] = []
    meta: dict[str, dict[str, Any]] = {}
    subspaces: list[dict[str, Any]] = []

    for eps in EPS_VALUES:
        base = (1.0 - eps) ** (1.0 / DEGREE) * unit_roots
        eps_tag = f"{eps:.0e}"
        u0 = directions[eps][0]
        u1 = directions[eps][1]
        subspaces.append(
            {
                "eps": eps,
                "basis": [
                    {
                        "rank": int(u0["rank"]),
                        "midpoint_eigenvalue": float(u0["midpoint_eigenvalue"]),
                        "dominant_shape_components": u0["dominant_shape_components"],
                    },
                    {
                        "rank": int(u1["rank"]),
                        "midpoint_eigenvalue": float(u1["midpoint_eigenvalue"]),
                        "dominant_shape_components": u1["dominant_shape_components"],
                    },
                ],
            }
        )
        for a in COEFFICIENT_GRID:
            for b in COEFFICIENT_GRID:
                eta_norm = float(np.sqrt(a * a + b * b))
                if eta_norm > ETA + 1e-15:
                    continue
                scale = eps_shape_scale(eps)
                direction = a * u0["direction"] + b * u1["direction"]
                roots = base + scale * direction
                label = f"eps:{eps_tag}:u0:{a:+.6f}:u1:{b:+.6f}"
                add_point(
                    seen,
                    points,
                    meta,
                    tensor,
                    label,
                    roots,
                    {
                        "kind": "low_dim_coefficient_disk_grid",
                        "eps": eps,
                        "u0_coefficient": float(a),
                        "u1_coefficient": float(b),
                        "eta_norm": eta_norm,
                        "eta_cap": ETA,
                        "eps_shape_scale": scale,
                        "target": scalar_target(eps),
                    },
                )
    return points, meta, subspaces


def build_result() -> dict[str, Any]:
    root = math_root_from_script()
    out_dir = Path(__file__).resolve().parent
    erdos114_dir = erdos114_dir_from_script()
    packet = load_module("ehp114_n14_interval_taylor_m14_packet_for_result", erdos114_dir / "ehp114_n14_interval_taylor_m14_packet.py")
    tensor = load_module("ehp114_tensor_cone_scaffold_for_result", root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py")
    _, directions = build_low_dim_directions()
    points, meta, subspaces = build_points(tensor, directions)

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
                "eps_values": EPS_VALUES,
                "eta_cap": ETA,
                "coefficient_grid": COEFFICIENT_GRID,
                "subspace": "span of lowest two admissible-Taylor midpoint eigendirections per epsilon",
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
    evaluated_rows = [row for row in rows if row["evaluated"]]
    failures = [row for row in evaluated_rows if not row["pass"]]
    eps_summaries = []
    for eps in EPS_VALUES:
        subset = [row for row in rows if row["eps"] == eps]
        evaluated = [row for row in subset if row["evaluated"]]
        eps_summaries.append(
            {
                "eps": eps,
                "grid_points": len(subset),
                "evaluated_admissible_points": len(evaluated),
                "inadmissible_points": len(subset) - len(evaluated),
                "all_evaluated_points_pass": all(row["pass"] for row in evaluated),
                "min_margin": min((row["margin"] for row in evaluated if row["margin"] is not None), default=None),
                "worst_row": min(evaluated, key=lambda row: row["margin"] if row["margin"] is not None else float("inf"))
                if evaluated
                else None,
            }
        )
    status = "LOW_DIM_CONE_GRID_PASS" if evaluated_rows and not failures else "LOW_DIM_CONE_GRID_FAIL"
    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "parent_packet": PARENT_ID,
        "status": status,
        "parameters": {
            "degree": DEGREE,
            "eps_values": EPS_VALUES,
            "eta_cap": ETA,
            "coefficient_grid": COEFFICIENT_GRID,
            "res": RES,
            "extent": EXTENT,
            "oracle_point_count": len(points),
            "grid_points_before_admissibility_filter": len(rows),
            "scalar_target": "D14(eps,s) >= 12 * eps^(1/14)",
            "cone_scaling": "s = eps^(1/28) * (a*u0(eps) + b*u1(eps)), sqrt(a^2+b^2) <= 0.014",
            "subspace": "lowest two midpoint eigendirections from the admissible spectral Taylor matrix for each epsilon",
        },
        "subspaces": subspaces,
        "eps_summaries": eps_summaries,
        "failure_count": len(failures),
        "failures": failures[:20],
        "rows": rows,
        "oracle": {
            "input": str(input_path),
            "output": str(output_path),
            "output_point_count": len(output["points"]),
        },
        "claim_ceiling": (
            "Finite 2D coefficient-disk grid evidence for the epsilon-scaled scalar theorem target at n=14. "
            "This is not a continuous cone certificate, not a Lean proof, and not a proof of Erdős #114."
        ),
        "next_blocker": (
            "Promote from sampled grid evidence to interval boxes over the two eigencoefficients, "
            "then bind those boxes to an admissibility certificate and a scalar-deficit lower bound."
        ),
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    lines = [
        "# EHP114 n=14 Low-Dimensional Cone Box Probe",
        "",
        f"Experiment: `{EXPERIMENT_ID}`",
        "",
        f"Parent packet: `{PARENT_ID}`",
        "",
        "## Meaning",
        "",
        "This run samples the two-dimensional coefficient disk generated by the two lowest",
        "midpoint eigendirections from the admissible spectral Taylor matrix. It tests the",
        "same scalar target as the epsilon-scaled axis and spectral-direction runs, but now",
        "on mixed points inside the low-dimensional span.",
        "",
        "This is sampled grid evidence, not a continuous cone certificate.",
        "",
        "## Verdict",
        "",
        f"- Status: `{result['status']}`",
        f"- Failure count among evaluated admissible grid points: `{result['failure_count']}`",
        f"- Oracle points evaluated: `{result['parameters']['oracle_point_count']}`",
        f"- Grid points before admissibility filter: `{result['parameters']['grid_points_before_admissibility_filter']}`",
        "",
        "## Epsilon Grid Summary",
        "",
        "| eps | admissible/evaluated | inadmissible | all evaluated pass | min margin |",
        "|---:|---:|---:|:---:|---:|",
    ]
    for row in result["eps_summaries"]:
        min_margin = "" if row["min_margin"] is None else f"{row['min_margin']:.12g}"
        lines.append(
            f"| {row['eps']:.0e} | {row['evaluated_admissible_points']}/{row['grid_points']} | "
            f"{row['inadmissible_points']} | {'yes' if row['all_evaluated_points_pass'] else 'no'} | {min_margin} |"
        )
    lines.extend(
        [
            "",
            "## Claim Ceiling",
            "",
            result["claim_ceiling"],
            "",
            "## Next Blocker",
            "",
            result["next_blocker"],
            "",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


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
                "failure_count": result["failure_count"],
                "oracle_point_count": result["parameters"]["oracle_point_count"],
                "report": str(report_path),
                "sha256": sha,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
