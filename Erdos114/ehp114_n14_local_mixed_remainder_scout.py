#!/usr/bin/env python3
"""Local mixed-remainder scout for the n=14 EHP114 cone.

This is the repaired version of the admissible scout. The first admissible
scout correctly showed that the unit disk alone permits shape excursions far
outside the local quadratic regime. Here we enforce both constraints:

  1. roots remain in the unit disk;
  2. shape amplitude stays below a fixed local cap.

The result is still a finite scout, not a theorem.
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

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-LOCAL-MIXED-REMAINDER-SCOUT-20260505-02"
FAILED_PREDECESSOR_ID = "EXP-MATH-EHP114-N14-ADMISSIBLE-MIXED-REMAINDER-SCOUT-20260505-01"
SHAPE_MATRIX_ID = "EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01"
DEGREE = 14
RES = 220
EXTENT = 3.0
LAMBDA_HALF = 50000.0
RADIAL_HALF = 12.0
SAFETY = 0.9
LOCAL_T_CAP = 0.008
EPS_VALUES = [1e-4, 1e-3, 1e-2, 1e-1]


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


def max_admissible_t(base: np.ndarray, direction: np.ndarray) -> float:
    lo = 0.0
    hi = 0.1
    while np.max(np.abs(base + hi * direction)) <= 1.0 and hi < 10.0:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if np.max(np.abs(base + mid * direction)) <= 1.0:
            lo = mid
        else:
            hi = mid
    return lo


def build_points(tensor: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    roots = tensor.roots_of_unity(DEGREE)
    basis = tensor.quotient_basis(DEGREE)
    shape = [b for b in basis if b.label != "radial_singular_m0"]
    rows: list[dict[str, Any]] = []
    points: list[dict[str, Any]] = []

    for eps in EPS_VALUES:
        radius = (1.0 - eps) ** (1.0 / DEGREE)
        base = radius * roots
        for i, b in enumerate(shape):
            for sign, signed in [("+", b.vector), ("-", -b.vector)]:
                tmax = max_admissible_t(base, signed)
                t = min(SAFETY * tmax, LOCAL_T_CAP)
                perturbed = base + t * signed
                label = f"eps:{eps:.0e}:dir:{i}:{sign}"
                rows.append(
                    {
                        "label": label,
                        "eps": eps,
                        "shape_index": i,
                        "shape_label": b.label,
                        "sign": sign,
                        "radius": radius,
                        "tmax": tmax,
                        "t_used": t,
                        "local_cap_active": t < SAFETY * tmax,
                        "max_root_radius": float(np.max(np.abs(perturbed))),
                    }
                )
                points.append({"label": label, "coeffs": coeffs_json(tensor.coeffs_from_roots(perturbed))})
    return rows, points


def run_length_oracle(root: Path, input_path: Path, output_path: Path) -> None:
    binary = root / "erdos-experiments" / "scripts" / "erdos-114" / "target" / "release" / "ehp114_batch_interval_lengths"
    subprocess.run(
        [str(binary), "--input", str(input_path), "--output", str(output_path), "--quiet"],
        cwd=root / "erdos-experiments" / "scripts" / "erdos-114",
        check=True,
    )


def build_result() -> dict[str, Any]:
    root = repo_root_from_script()
    out_dir = Path(__file__).resolve().parent
    tensor = load_module(
        "ehp114_tensor_cone_scaffold",
        root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py",
    )
    row_meta, points = build_points(tensor)
    input_path = out_dir / f"{EXPERIMENT_ID}_ORACLE_INPUT.json"
    output_path = out_dir / f"{EXPERIMENT_ID}_ORACLE_OUTPUT.json"
    for path in (input_path, output_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing oracle artifact: {path}")

    input_path.write_text(
        json.dumps({"degree": DEGREE, "res": RES, "extent": EXTENT, "points": points}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    run_length_oracle(root, input_path, output_path)
    output = load_json(output_path)
    length_by_label = {row["label"]: row for row in output["points"]}
    inari = load_json(root / "erdos-experiments" / "results" / "erdos-114" / "EXP-MM-EHP-007-n14-inari_RESULTS.json")
    lstar_lower = float(inari["l_star_lower"])

    rows = []
    for row in row_meta:
        length = length_by_label[row["label"]]
        deficit_lower = lstar_lower - float(length["length_upper"])
        rhs = RADIAL_HALF * (row["eps"] ** (1.0 / DEGREE)) + LAMBDA_HALF * row["t_used"] * row["t_used"]
        rows.append(
            {
                **row,
                "length_upper": float(length["length_upper"]),
                "deficit_lower": deficit_lower,
                "rhs_mixed_target": rhs,
                "margin": deficit_lower - rhs,
                "pass": deficit_lower > rhs,
            }
        )

    all_pass = all(row["pass"] for row in rows)
    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "scope": "finite local mixed-remainder scout for n=14 local cone",
        "status": "LOCAL_MIXED_SCOUT_PASS" if all_pass else "LOCAL_MIXED_SCOUT_FAIL",
        "claim_ceiling": "This is a shadow signature, not universal law. It is a finite local scout, not a uniform theorem.",
        "predecessor": {
            "experiment_id": FAILED_PREDECESSOR_ID,
            "lesson": "unit-disk admissibility alone permits nonlocal tangential excursions and is too broad for the quadratic shape-cone theorem",
        },
        "parameters": {
            "degree": DEGREE,
            "eps_values": EPS_VALUES,
            "signed_points": len(rows),
            "res": RES,
            "extent": EXTENT,
            "radial_half_constant": RADIAL_HALF,
            "lambda_half": LAMBDA_HALF,
            "admissible_safety": SAFETY,
            "local_t_cap": LOCAL_T_CAP,
        },
        "sources": {
            "shape_matrix": str(out_dir / f"{SHAPE_MATRIX_ID}_RESULTS.json"),
            "failed_predecessor": str(out_dir / f"{FAILED_PREDECESSOR_ID}_RESULTS.json"),
        },
        "oracle": {
            "input_path": str(input_path),
            "input_sha256": sha256_file(input_path),
            "output_path": str(output_path),
            "output_sha256": sha256_file(output_path),
            "elapsed_secs": output["elapsed_secs"],
        },
        "all_points_pass": all_pass,
        "min_margin": min(row["margin"] for row in rows),
        "worst_point": min(rows, key=lambda row: row["margin"]),
        "cap_active_count": sum(1 for row in rows if row["local_cap_active"]),
        "rows": rows,
        "next_blocker": "Prove a uniform local cone radius, then bound mixed radial/shape remainder on that cone.",
        "guardrails": {
            "no_scorecard_update": True,
            "no_d1_update": True,
            "no_public_claim": True,
            "no_lean_file_created": True,
            "no_email": True,
            "no_git": True,
        },
    }


def build_report(result: dict[str, Any]) -> str:
    by_eps = []
    for eps in result["parameters"]["eps_values"]:
        rows = [row for row in result["rows"] if row["eps"] == eps]
        by_eps.append(
            {
                "eps": eps,
                "min_margin": min(row["margin"] for row in rows),
                "pass_count": sum(1 for row in rows if row["pass"]),
                "count": len(rows),
                "cap_active": sum(1 for row in rows if row["local_cap_active"]),
            }
        )
    table = "\n".join(
        f"| {row['eps']:.0e} | {row['pass_count']}/{row['count']} | {row['min_margin']:.12g} | {row['cap_active']} |"
        for row in by_eps
    )
    worst = result["worst_point"]
    return f"""# EHP114 n=14 Local Mixed-Remainder Scout

Experiment: `{EXPERIMENT_ID}`

## Meaning

The unrestricted admissible scout failed because the unit disk allows large
tangential excursions that are not local shape-cone perturbations. This repaired
scout adds the missing local condition: `||s|| <= {LOCAL_T_CAP}`.

The claim ceiling remains narrow: this is a shadow signature, not universal
law. It is a finite local scout, not a uniform theorem.

## Verdict

- Status: `{result["status"]}`
- Signed local points: `{result["parameters"]["signed_points"]}`
- Local shape cap: `{result["parameters"]["local_t_cap"]}`
- All points pass: `{result["all_points_pass"]}`
- Minimum margin: `{result["min_margin"]}`
- Cap active count: `{result["cap_active_count"]}`
- Next blocker: `{result["next_blocker"]}`

## Per-Epsilon Summary

| eps | pass | min margin | cap active |
|---:|---:|---:|---:|
{table}

Worst point: `{worst["label"]}` with margin `{worst["margin"]}`.

## What This Teaches

The local theorem needs two gates:

```text
roots stay in the unit disk
||s|| <= eta14
```

The first gate is admissibility. The second gate is locality. The predecessor
only used the first gate and therefore failed. This run uses both and gives the
next theorem a defensible domain.

## What Remains

Prove a uniform analytic version:

```text
R14(eps,s) <= 12 eps^(1/14) + 50000 ||s||^2
for 0 < eps <= 1e-1 and ||s|| <= min(eta14_boundary(eps), 0.008).
```

No scorecard, D1, public document, Lean file, git, or email state was changed.
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
                "all_points_pass": result["all_points_pass"],
                "min_margin": result["min_margin"],
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
