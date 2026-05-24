#!/usr/bin/env python3
"""Epsilon-scaled scalar search along negative spectral directions.

The axis scan tests signed basis directions. This run extracts the lowest
midpoint eigenvectors from the admissible-stencil Taylor matrices and evaluates
the weaker scalar target

    D14(eps, s) >= 12 * eps^(1/14)

along those mixed directions after the normalization s = eta * eps^(1/28) u.
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

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-EPS-SCALED-SPECTRAL-DIRECTION-SEARCH-20260505-01"
PARENT_ID = "EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01"
DEGREE = 14
EPS_VALUES = [1e-4, 1e-3, 1e-2, 1e-1]
ETA_GRID = [0.001, 0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014]
EIGENVECTOR_COUNT = 3
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


def build_interval_matrix(packet: Any, eps: float, h: float, dim: int, deficits: dict[str, Any]) -> list[list[Any]]:
    eps_tag = f"{eps:.0e}"
    base = deficits[f"eps:{eps_tag}:base"]
    h2 = h * h
    matrix = [[packet.Interval(0.0, 0.0) for _ in range(dim)] for _ in range(dim)]
    for i in range(dim):
        plus = deficits[f"eps:{eps_tag}:diag:{i}:+"]
        minus = deficits[f"eps:{eps_tag}:diag:{i}:-"]
        matrix[i][i] = (plus + minus - base.scale(2.0)).scale(0.5 / h2)
    for i in range(dim):
        for j in range(i + 1, dim):
            matrix[i][j] = (
                deficits[f"eps:{eps_tag}:off:{i}:{j}:++"]
                - deficits[f"eps:{eps_tag}:off:{i}:{j}:+-"]
                - deficits[f"eps:{eps_tag}:off:{i}:{j}:-+"]
                + deficits[f"eps:{eps_tag}:off:{i}:{j}:--"]
            ).scale(1.0 / (4.0 * h2))
            matrix[j][i] = matrix[i][j]
    return matrix


def eigen_directions(packet: Any, tensor: Any, root: Path) -> tuple[list[Any], dict[float, list[dict[str, Any]]]]:
    out_dir = Path(__file__).resolve().parent
    spectral = load_json(out_dir / f"{PARENT_ID}_RESULTS.json")
    oracle = load_json(out_dir / f"{PARENT_ID}_ORACLE_OUTPUT.json")
    shape = [b for b in tensor.quotient_basis(DEGREE) if b.label != "radial_singular_m0"]
    inari = load_json(root / "erdos-experiments" / "results" / "erdos-114" / "EXP-MM-EHP-007-n14-inari_RESULTS.json")
    lstar = packet.Interval(float(inari["l_star_lower"]), float(inari["l_star_upper"]))
    lengths = {
        row["label"]: packet.Interval(float(row["length_lower"]), float(row["length_upper"]))
        for row in oracle["points"]
    }
    deficits = {label: packet.Interval(lstar.lo - length.hi, lstar.hi - length.lo) for label, length in lengths.items()}
    h_by_eps = {float(row["eps"]): float(row["admissible_taylor_step"]) for row in spectral["eps_rows"]}
    directions: dict[float, list[dict[str, Any]]] = {}
    for eps in EPS_VALUES:
        matrix = build_interval_matrix(packet, eps, h_by_eps[eps], len(shape), deficits)
        midpoint = np.array([[(x.lo + x.hi) / 2.0 for x in row] for row in matrix], dtype=float)
        eigvals, eigvecs = np.linalg.eigh(midpoint)
        eps_dirs = []
        for rank in range(EIGENVECTOR_COUNT):
            coeffs = eigvecs[:, rank]
            direction = sum(float(c) * b.vector for c, b in zip(coeffs, shape))
            eps_dirs.append(
                {
                    "rank": rank,
                    "midpoint_eigenvalue": float(eigvals[rank]),
                    "coefficients": [float(c) for c in coeffs],
                    "dominant_shape_components": [
                        {"shape_index": int(i), "shape_label": shape[int(i)].label, "coefficient": float(coeffs[int(i)])}
                        for i in np.argsort(np.abs(coeffs))[-5:][::-1]
                    ],
                    "direction": direction,
                }
            )
        directions[eps] = eps_dirs
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
    root_radius = max_root_radius(roots)
    admissible = root_radius <= 1.0 + ROOT_TOL
    meta[label] = {**payload, "max_root_radius": root_radius, "admissible": admissible}
    if not admissible or label in seen:
        return
    seen.add(label)
    points.append({"label": label, "coeffs": coeffs_json(tensor.coeffs_from_roots(roots))})


def build_points(tensor: Any, directions: dict[float, list[dict[str, Any]]]) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    unit_roots = tensor.roots_of_unity(DEGREE)
    seen: set[str] = set()
    points: list[dict[str, Any]] = []
    meta: dict[str, dict[str, Any]] = {}
    for eps in EPS_VALUES:
        base = (1.0 - eps) ** (1.0 / DEGREE) * unit_roots
        eps_tag = f"{eps:.0e}"
        for d in directions[eps]:
            for eta in ETA_GRID:
                eta_tag = f"{eta:.6f}"
                t = eta * eps_shape_scale(eps)
                for sign_label, sign in [("+", 1.0), ("-", -1.0)]:
                    label = f"eps:{eps_tag}:rank:{d['rank']}:eta:{eta_tag}:spectral:{sign_label}"
                    add_point(
                        seen,
                        points,
                        meta,
                        tensor,
                        label,
                        base + sign * t * d["direction"],
                        {
                            "kind": "spectral_direction",
                            "eps": eps,
                            "rank": d["rank"],
                            "midpoint_eigenvalue": d["midpoint_eigenvalue"],
                            "dominant_shape_components": d["dominant_shape_components"],
                            "eta": eta,
                            "t": t,
                            "sign": sign_label,
                            "target": scalar_target(eps),
                            "eps_shape_scale": eps_shape_scale(eps),
                        },
                    )
    return points, meta


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


def build_result() -> dict[str, Any]:
    root = repo_root_from_script()
    out_dir = Path(__file__).resolve().parent
    packet = load_module("ehp114_n14_interval_taylor_m14_packet", out_dir / "ehp114_n14_interval_taylor_m14_packet.py")
    tensor = load_module("ehp114_tensor_cone_scaffold", root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py")
    shape, directions = eigen_directions(packet, tensor, root)
    points, meta = build_points(tensor, directions)
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
                "eigenvector_count": EIGENVECTOR_COUNT,
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
    eta_rows = []
    for eta in ETA_GRID:
        subset = [row for row in rows if abs(row["eta"] - eta) < 5e-10]
        evaluated = [row for row in subset if row["evaluated"]]
        eta_rows.append(
            {
                "eta": eta,
                "expected_points": len(subset),
                "evaluated_admissible_points": len(evaluated),
                "inadmissible_points": len(subset) - len(evaluated),
                "all_evaluated_points_pass": all(row["pass"] for row in evaluated),
                "min_margin": min((row["margin"] for row in evaluated if row["margin"] is not None), default=None),
                "worst_row": min(evaluated, key=lambda row: row["margin"] if row["margin"] is not None else float("inf"))
                if evaluated
                else None,
            }
        )
    certified_eta = [row["eta"] for row in eta_rows if row["evaluated_admissible_points"] > 0 and row["all_evaluated_points_pass"]]
    status = "EPS_SCALED_SPECTRAL_DIRECTIONS_PASS" if not failures and evaluated_rows else "EPS_SCALED_SPECTRAL_DIRECTIONS_FAIL"
    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "parent_packet": PARENT_ID,
        "status": status,
        "parameters": {
            "degree": DEGREE,
            "eps_values": EPS_VALUES,
            "eta_grid": ETA_GRID,
            "eigenvector_count": EIGENVECTOR_COUNT,
            "res": RES,
            "extent": EXTENT,
            "shape_basis_rank": len(shape),
            "oracle_point_count": len(points),
            "scalar_target": "D14(eps,s) >= 12 * eps^(1/14)",
            "cone_scaling": "||s|| = eta * eps^(1/28) along lowest midpoint eigenvectors",
        },
        "largest_eta_with_all_evaluated_passing": max(certified_eta, default=None),
        "eta_summaries": eta_rows,
        "failure_count": len(failures),
        "failures": failures[:20],
        "rows": rows,
        "oracle": {
            "input": str(input_path),
            "output": str(output_path),
            "output_point_count": len(output["points"]),
        },
        "claim_ceiling": (
            "Finite mixed-direction evidence for the epsilon-scaled scalar theorem target at n=14. "
            "This is not a full cone certificate, not a Lean proof, and not a proof of Erdős #114."
        ),
        "next_blocker": (
            "Promote the passing admissible spectral-direction evidence to a certified low-dimensional cone "
            "or derive analytic Cauchy bounds that remove reliance on sampled eigenvectors."
        ),
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    lines = [
        "# EHP114 n=14 Epsilon-Scaled Spectral-Direction Search",
        "",
        f"Experiment: `{EXPERIMENT_ID}`",
        "",
        f"Parent packet: `{PARENT_ID}`",
        "",
        "## Meaning",
        "",
        "This run follows the negative midpoint eigenvectors from the admissible Taylor matrices.",
        "It tests the scalar target on mixed directions rather than only signed coordinate axes.",
        "",
        "## Verdict",
        "",
        f"- Status: `{result['status']}`",
        f"- Largest eta with all evaluated admissible spectral points passing: `{result['largest_eta_with_all_evaluated_passing']}`",
        f"- Failure count among evaluated admissible points: `{result['failure_count']}`",
        f"- Oracle points evaluated: `{result['parameters']['oracle_point_count']}`",
        "",
        "## Eta Grid",
        "",
        "| eta | admissible/evaluated | inadmissible | all pass | min margin |",
        "|---:|---:|---:|:---:|---:|",
    ]
    for row in result["eta_summaries"]:
        min_margin = "" if row["min_margin"] is None else f"{row['min_margin']:.12g}"
        lines.append(
            f"| {row['eta']:.3f} | {row['evaluated_admissible_points']}/{row['expected_points']} | "
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
                "largest_eta_with_all_evaluated_passing": result["largest_eta_with_all_evaluated_passing"],
                "failure_count": result["failure_count"],
                "report": str(report_path),
                "sha256": sha,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
