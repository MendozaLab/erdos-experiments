#!/usr/bin/env python3
"""Admissible-stencil spectral scan for EHP114 n=14.

The previous interval Taylor packet used an ambient central-difference stencil;
some off-axis stencil points were not root-admissible. This run shrinks the
Taylor step separately for each epsilon so every diagonal and off-axis stencil
point remains inside the unit-disk root-admissible region.

If the spectral lower bound remains negative, the radial-base shape softening is
not merely an inadmissible-stencil artifact.
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

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01"
PARENT_ID = "EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-M14-PACKET-20260505-01"
DEGREE = 14
EPS_VALUES = [1e-4, 1e-3, 1e-2, 1e-1]
MAX_TAYLOR_STEP = 0.004
RES = 220
EXTENT = 3.0
WORKING_LAMBDA14 = 100000.0


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


def repo_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def coeffs_json(coeffs: np.ndarray) -> list[list[float]]:
    return [[float(z.real), float(z.imag)] for z in coeffs]


def max_root_radius(roots: np.ndarray) -> float:
    return float(np.max(np.abs(roots)))


def all_stencil_admissible(base: np.ndarray, shape: list[Any], h: float) -> bool:
    for b in shape:
        for sign in (1.0, -1.0):
            if max_root_radius(base + h * sign * b.vector) > 1.0 + 1e-12:
                return False
    for i, bi in enumerate(shape):
        for j in range(i + 1, len(shape)):
            bj = shape[j]
            for a, c in ((1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)):
                if max_root_radius(base + h * (a * bi.vector + c * bj.vector)) > 1.0 + 1e-12:
                    return False
    return True


def admissible_step(base: np.ndarray, shape: list[Any]) -> float:
    lo = 0.0
    hi = MAX_TAYLOR_STEP
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if all_stencil_admissible(base, shape, mid):
            lo = mid
        else:
            hi = mid
    return lo


def add_point(
    tensor: Any,
    points: list[dict[str, Any]],
    meta: dict[str, dict[str, Any]],
    label: str,
    roots: np.ndarray,
    payload: dict[str, Any],
) -> None:
    points.append({"label": label, "coeffs": coeffs_json(tensor.coeffs_from_roots(roots))})
    meta[label] = {**payload, "max_root_radius": max_root_radius(roots)}


def build_points(tensor: Any) -> tuple[list[Any], dict[float, float], list[dict[str, Any]], dict[str, dict[str, Any]]]:
    unit_roots = tensor.roots_of_unity(DEGREE)
    shape = [b for b in tensor.quotient_basis(DEGREE) if b.label != "radial_singular_m0"]
    steps: dict[float, float] = {}
    points: list[dict[str, Any]] = []
    meta: dict[str, dict[str, Any]] = {}

    for eps in EPS_VALUES:
        base = (1.0 - eps) ** (1.0 / DEGREE) * unit_roots
        h = admissible_step(base, shape)
        steps[eps] = h
        eps_tag = f"{eps:.0e}"
        add_point(tensor, points, meta, f"eps:{eps_tag}:base", base, {"kind": "base", "eps": eps, "h": h})

        for i, b in enumerate(shape):
            for sign_label, signed in [("+", b.vector), ("-", -b.vector)]:
                add_point(
                    tensor,
                    points,
                    meta,
                    f"eps:{eps_tag}:diag:{i}:{sign_label}",
                    base + h * signed,
                    {"kind": "diag", "eps": eps, "shape_index": i, "shape_label": b.label, "h": h},
                )

        for i, bi in enumerate(shape):
            for j in range(i + 1, len(shape)):
                bj = shape[j]
                for signs, direction in [
                    ("++", bi.vector + bj.vector),
                    ("+-", bi.vector - bj.vector),
                    ("-+", -bi.vector + bj.vector),
                    ("--", -bi.vector - bj.vector),
                ]:
                    add_point(
                        tensor,
                        points,
                        meta,
                        f"eps:{eps_tag}:off:{i}:{j}:{signs}",
                        base + h * direction,
                        {
                            "kind": "off",
                            "eps": eps,
                            "shape_i": i,
                            "shape_j": j,
                            "shape_i_label": bi.label,
                            "shape_j_label": bj.label,
                            "h": h,
                        },
                    )
    return shape, steps, points, meta


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


def build_matrix(packet: Any, eps: float, h: float, dim: int, deficits: dict[str, Any]) -> tuple[list[list[Any]], list[dict[str, Any]]]:
    eps_tag = f"{eps:.0e}"
    base = deficits[f"eps:{eps_tag}:base"]
    matrix = [[packet.Interval(0.0, 0.0) for _ in range(dim)] for _ in range(dim)]
    h2 = h * h
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
    rows = []
    for i in range(dim):
        diag = matrix[i][i]
        radius = sum(matrix[i][j].abs_upper() for j in range(dim) if j != i)
        rows.append(
            {
                "index": i,
                "diagonal_interval": diag.to_json(),
                "off_diagonal_abs_radius_upper": radius,
                "diagonal_lower_minus_radius": diag.lo - radius,
            }
        )
    return matrix, rows


def build_result() -> dict[str, Any]:
    root = repo_root_from_script()
    out_dir = Path(__file__).resolve().parent
    packet = load_module("ehp114_n14_interval_taylor_m14_packet", out_dir / "ehp114_n14_interval_taylor_m14_packet.py")
    tensor = load_module("ehp114_tensor_cone_scaffold", root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py")
    shape, steps, points, meta = build_points(tensor)
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
                "admissible_steps": {f"{k:.0e}": v for k, v in steps.items()},
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
    deficits = {label: packet.Interval(lstar.lo - length.hi, lstar.hi - length.lo) for label, length in lengths.items()}

    eps_rows = []
    for eps in EPS_VALUES:
        matrix, gersh_rows = build_matrix(packet, eps, steps[eps], len(shape), deficits)
        mid = np.array([[(x.lo + x.hi) / 2.0 for x in row] for row in matrix], dtype=float)
        rad = np.array([[(x.hi - x.lo) / 2.0 for x in row] for row in matrix], dtype=float)
        eigvals = np.linalg.eigvalsh(mid)
        fro = float(np.linalg.norm(rad, "fro"))
        lower = float(eigvals[0] - fro)
        eps_rows.append(
            {
                "eps": eps,
                "admissible_taylor_step": steps[eps],
                "midpoint_lambda_min": float(eigvals[0]),
                "midpoint_lambda_max": float(eigvals[-1]),
                "frobenius_interval_radius": fro,
                "interval_spectral_lower_bound": lower,
                "gershgorin_lower_bound": min(row["diagonal_lower_minus_radius"] for row in gersh_rows),
                "working_lambda14_pass": lower > WORKING_LAMBDA14,
                "positive_by_interval_spectral": lower > 0.0,
            }
        )
    all_positive = all(row["positive_by_interval_spectral"] for row in eps_rows)
    if all(row["working_lambda14_pass"] for row in eps_rows):
        status = "ADMISSIBLE_SPECTRAL_LAMBDA_PASS"
    elif all_positive:
        status = "ADMISSIBLE_SPECTRAL_POSITIVE_BUT_LAMBDA_FAIL"
    else:
        status = "ADMISSIBLE_SHAPE_SOFTENING_CONFIRMED"
    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "parent_packet": PARENT_ID,
        "status": status,
        "claim_ceiling": "This is a shadow signature, not universal law. It is an admissible-stencil Taylor diagnosis, not a proof of local stability.",
        "parameters": {
            "degree": DEGREE,
            "eps_values": EPS_VALUES,
            "shape_basis_rank": len(shape),
            "max_taylor_step": MAX_TAYLOR_STEP,
            "res": RES,
            "extent": EXTENT,
            "point_count": len(points),
        },
        "oracle": {
            "input_path": str(input_path),
            "input_sha256": sha256_file(input_path),
            "output_path": str(output_path),
            "output_sha256": sha256_file(output_path),
            "elapsed_secs": output["elapsed_secs"],
        },
        "all_stencil_points_admissible": all(row["max_root_radius"] <= 1.0 + 1e-12 for row in meta.values()),
        "max_root_radius": max(row["max_root_radius"] for row in meta.values()),
        "eps_rows": eps_rows,
        "global_interval_spectral_lower_bound": min(row["interval_spectral_lower_bound"] for row in eps_rows),
        "all_positive_by_interval_spectral": all_positive,
        "sources": {
            "parent_results": str(out_dir / f"{PARENT_ID}_RESULTS.json"),
            "parent_results_sha256": sha256_file(out_dir / f"{PARENT_ID}_RESULTS.json"),
            "inari_n14": str(inari_path),
            "inari_n14_sha256": sha256_file(inari_path),
        },
        "next_blocker": "Use epsilon-scaled coordinates and target a direct scalar total-deficit theorem instead of transported positive shape cone.",
        "guardrails": {
            "no_scorecard_update": True,
            "no_d1_update": True,
            "no_public_claim": True,
            "no_lean_status_change": True,
            "no_email": True,
            "no_git": True,
        },
    }


def build_report(result: dict[str, Any]) -> str:
    table = "\n".join(
        "| {eps:.0e} | {h:.6g} | {lower:.12g} | {g:.12g} | {ok} |".format(
            eps=row["eps"],
            h=row["admissible_taylor_step"],
            lower=row["interval_spectral_lower_bound"],
            g=row["gershgorin_lower_bound"],
            ok="yes" if row["positive_by_interval_spectral"] else "no",
        )
        for row in result["eps_rows"]
    )
    return f"""# EHP114 n=14 Admissible-Stencil Spectral Taylor Scan

Experiment: `{EXPERIMENT_ID}`

Parent packet: `{PARENT_ID}`

## Meaning

The previous Taylor matrix used an ambient stencil with some non-admissible
off-axis points. This scan shrinks the Taylor step per epsilon until every
diagonal and off-axis stencil point is root-admissible.

## Verdict

- Status: `{result["status"]}`
- All stencil points admissible: `{result["all_stencil_points_admissible"]}`
- Max root radius: `{result["max_root_radius"]}`
- Global interval spectral lower bound: `{result["global_interval_spectral_lower_bound"]}`
- All positive by interval spectral bound: `{result["all_positive_by_interval_spectral"]}`

## Per-Epsilon Scan

| eps | admissible h | spectral lower | Gershgorin lower | positive |
|---:|---:|---:|---:|---:|
{table}

## Consequence

If the lower bound remains negative here, radial-base shape softening survives
the admissibility correction. The next target is the epsilon-scaled total
deficit theorem, not a transported positive shape-cone theorem.

## Claim Ceiling

This is a diagnostic of a fixed-n local route. It does not settle Erdős #114 and
does not prove local stability.
"""


def main() -> int:
    out_dir = Path(__file__).resolve().parent
    result_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {path}")
    result = build_result()
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_path.write_text(build_report(result), encoding="utf-8")
    sha_path.write_text(f"{sha256_file(result_path)}  {result_path.name}\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "experiment_id": EXPERIMENT_ID,
                "status": result["status"],
                "global_interval_spectral_lower_bound": result["global_interval_spectral_lower_bound"],
                "report": str(report_path),
                "sha256": sha_path.read_text(encoding="utf-8").split()[0],
                "next_blocker": result["next_blocker"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

