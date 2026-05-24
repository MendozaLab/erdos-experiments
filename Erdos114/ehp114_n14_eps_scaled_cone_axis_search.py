#!/usr/bin/env python3
"""Epsilon-scaled cone search for the EHP114 n=14 local theorem target.

The admissible-stencil spectral scan shows that transported uniform shape
positivity is the wrong target. This run tests the weaker scalar target

    D14(eps, s) >= 12 * eps^(1/14)

on axis directions in the quotient shape cone after the normalization

    s = eta * eps^(1/28) * u.

This is finite evidence for a theorem target, not a proof of Erdős #114.
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

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-EPS-SCALED-CONE-AXIS-SEARCH-20260505-01"
PARENT_ID = "EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01"
DEGREE = 14
EPS_VALUES = [1e-4, 1e-3, 1e-2, 1e-1]
ETA_GRID = [0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014, 0.016, 0.018, 0.020]
WORST_MODE_INDEX = 21
WORST_MODE_SIGN = "-"
WORST_ETA_GRID = [round(0.004 + 0.001 * k, 6) for k in range(33)]
RES = 220
EXTENT = 3.0
RADIAL_HALF = 12.0
ROOT_TOL = 1e-12


def repo_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


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


def add_point(
    seen: set[str],
    points: list[dict[str, Any]],
    meta: dict[str, dict[str, Any]],
    tensor: Any,
    label: str,
    roots: np.ndarray,
    payload: dict[str, Any],
) -> bool:
    root_radius = max_root_radius(roots)
    admissible = root_radius <= 1.0 + ROOT_TOL
    meta[label] = {**payload, "max_root_radius": root_radius, "admissible": admissible}
    if not admissible or label in seen:
        return False
    seen.add(label)
    points.append({"label": label, "coeffs": coeffs_json(tensor.coeffs_from_roots(roots))})
    return True


def build_points(tensor: Any) -> tuple[list[Any], list[dict[str, Any]], dict[str, dict[str, Any]]]:
    unit_roots = tensor.roots_of_unity(DEGREE)
    shape = [b for b in tensor.quotient_basis(DEGREE) if b.label != "radial_singular_m0"]
    seen: set[str] = set()
    points: list[dict[str, Any]] = []
    meta: dict[str, dict[str, Any]] = {}

    eta_values = sorted(set(ETA_GRID + WORST_ETA_GRID))
    for eps in EPS_VALUES:
        radius = (1.0 - eps) ** (1.0 / DEGREE)
        base = radius * unit_roots
        eps_tag = f"{eps:.0e}"
        for eta in eta_values:
            t = eta * eps_shape_scale(eps)
            eta_tag = f"{eta:.6f}"
            for i, b in enumerate(shape):
                for sign_label, signed in [("+", b.vector), ("-", -b.vector)]:
                    label = f"eps:{eps_tag}:eta:{eta_tag}:axis:{i}:{sign_label}"
                    add_point(
                        seen,
                        points,
                        meta,
                        tensor,
                        label,
                        base + t * signed,
                        {
                            "kind": "axis",
                            "eps": eps,
                            "eta": eta,
                            "t": t,
                            "shape_index": i,
                            "shape_label": b.label,
                            "sign": sign_label,
                            "target": scalar_target(eps),
                            "eps_shape_scale": eps_shape_scale(eps),
                            "is_worst_mode": i == WORST_MODE_INDEX and sign_label == WORST_MODE_SIGN,
                        },
                    )
    return shape, points, meta


def run_length_oracle(root: Path, input_path: Path, output_path: Path) -> None:
    binary = root / "erdos-experiments" / "scripts" / "erdos-114" / "target" / "release" / "ehp114_batch_interval_lengths"
    if not binary.exists():
        subprocess.run(
            ["cargo", "build", "--release", "--bin", "ehp114_batch_interval_lengths"],
            cwd=root / "erdos-experiments" / "scripts" / "erdos-114",
            check=True,
        )
    subprocess.run(
        [str(binary), "--input", str(input_path), "--output", str(output_path), "--quiet"],
        cwd=root / "erdos-experiments" / "scripts" / "erdos-114",
        check=True,
    )


def summarize_eta(rows: list[dict[str, Any]], eta: float) -> dict[str, Any]:
    eta_rows = [row for row in rows if abs(row["eta"] - eta) < 5e-10]
    evaluated = [row for row in eta_rows if row["evaluated"]]
    missing = len(eta_rows) - len(evaluated)
    passing = [row for row in evaluated if row["pass"]]
    margins = [row["margin"] for row in passing if row["margin"] is not None]
    return {
        "eta": eta,
        "expected_axis_points": len(eta_rows),
        "evaluated_admissible_points": len(evaluated),
        "inadmissible_points": missing,
        "all_axis_points_admissible": missing == 0,
        "all_evaluated_points_pass": len(evaluated) == len(passing),
        "all_axis_certified": missing == 0 and len(evaluated) == len(passing),
        "min_margin": min(margins) if margins else None,
        "worst_row": min(evaluated, key=lambda row: row["margin"] if row["margin"] is not None else float("inf"))
        if evaluated
        else None,
    }


def build_result() -> dict[str, Any]:
    root = repo_root_from_script()
    out_dir = Path(__file__).resolve().parent
    packet = load_module("ehp114_n14_interval_taylor_m14_packet", out_dir / "ehp114_n14_interval_taylor_m14_packet.py")
    tensor = load_module("ehp114_tensor_cone_scaffold", root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py")
    shape, points, meta = build_points(tensor)

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
                "eps_values": EPS_VALUES,
                "eta_grid": ETA_GRID,
                "worst_mode": {"shape_index": WORST_MODE_INDEX, "sign": WORST_MODE_SIGN},
                "worst_eta_grid": WORST_ETA_GRID,
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
    inari_path = root / "erdos-experiments" / "results" / "erdos-114" / "EXP-MM-EHP-007-n14-inari_RESULTS.json"
    inari = load_json(inari_path)
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

    eta_summaries = [summarize_eta(rows, eta) for eta in ETA_GRID]
    certified = [row for row in eta_summaries if row["all_axis_certified"]]
    largest_eta = max((row["eta"] for row in certified), default=None)
    worst_rows = [
        row for row in rows if row["is_worst_mode"] and row["evaluated"]
    ]
    worst_eta_summaries = [
        {
            "eta": eta,
            "rows": [
                row for row in worst_rows if abs(row["eta"] - eta) < 5e-10
            ],
        }
        for eta in WORST_ETA_GRID
    ]
    worst_certified = [
        eta_row["eta"]
        for eta_row in worst_eta_summaries
        if len(eta_row["rows"]) == len(EPS_VALUES) and all(row["pass"] for row in eta_row["rows"])
    ]

    status = "EPS_SCALED_AXIS_PASS" if largest_eta is not None else "EPS_SCALED_AXIS_FAIL"
    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "parent_packet": PARENT_ID,
        "status": status,
        "parameters": {
            "degree": DEGREE,
            "eps_values": EPS_VALUES,
            "eta_grid": ETA_GRID,
            "worst_mode_index": WORST_MODE_INDEX,
            "worst_mode_sign": WORST_MODE_SIGN,
            "worst_mode_label": shape[WORST_MODE_INDEX].label,
            "worst_eta_grid": WORST_ETA_GRID,
            "res": RES,
            "extent": EXTENT,
            "axis_basis_rank": len(shape),
            "oracle_point_count": len(points),
            "scalar_target": "D14(eps,s) >= 12 * eps^(1/14)",
            "cone_scaling": "||s|| = eta * eps^(1/28) on normalized axis directions",
        },
        "largest_eta_all_axis_certified": largest_eta,
        "eta_summaries": eta_summaries,
        "worst_mode_largest_eta_certified": max(worst_certified, default=None),
        "worst_mode_rows": worst_rows,
        "rows": rows,
        "oracle": {
            "input": str(input_path),
            "output": str(output_path),
            "output_point_count": len(output["points"]),
        },
        "claim_ceiling": (
            "Finite axis-direction evidence for the epsilon-scaled scalar theorem target at n=14. "
            "This is not a full cone certificate, not a Lean proof, and not a proof of Erdős #114."
        ),
        "next_blocker": (
            "Move from axis directions to low-dimensional spectral/admissible sub-cones, or prove analytic Cauchy "
            "bounds that dominate the observed radial-base shape softening."
        ),
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    lines = [
        "# EHP114 n=14 Epsilon-Scaled Cone Axis Search",
        "",
        f"Experiment: `{EXPERIMENT_ID}`",
        "",
        f"Parent packet: `{PARENT_ID}`",
        "",
        "## Meaning",
        "",
        "This run tests the weaker scalar target after the normalized change of variables",
        "`s = eta * eps^(1/28) * u`. The exponent is treated as a coordinate",
        "hypothesis, not a theorem.",
        "",
        "## Verdict",
        "",
        f"- Status: `{result['status']}`",
        f"- Largest eta certified on every axis/sign/epsilon in the coarse grid: `{result['largest_eta_all_axis_certified']}`",
        f"- Worst-mode largest eta certified: `{result['worst_mode_largest_eta_certified']}`",
        f"- Oracle points evaluated: `{result['parameters']['oracle_point_count']}`",
        "",
        "## Coarse Eta Grid",
        "",
        "| eta | admissible/evaluated | inadmissible | all pass | min margin |",
        "|---:|---:|---:|:---:|---:|",
    ]
    for row in result["eta_summaries"]:
        min_margin = "" if row["min_margin"] is None else f"{row['min_margin']:.12g}"
        lines.append(
            f"| {row['eta']:.3f} | {row['evaluated_admissible_points']}/{row['expected_axis_points']} | "
            f"{row['inadmissible_points']} | {'yes' if row['all_evaluated_points_pass'] else 'no'} | {min_margin} |"
        )
    lines.extend(
        [
            "",
            "## Worst Mode",
            "",
            f"Mode: index `{result['parameters']['worst_mode_index']}`, label "
            f"`{result['parameters']['worst_mode_label']}`, sign `{result['parameters']['worst_mode_sign']}`.",
            "",
            "This is the `m6_sin_tangent` danger lane identified by the interval Taylor packet.",
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
                "largest_eta_all_axis_certified": result["largest_eta_all_axis_certified"],
                "worst_mode_largest_eta_certified": result["worst_mode_largest_eta_certified"],
                "report": str(report_path),
                "sha256": sha,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
