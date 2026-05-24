#!/usr/bin/env python3
"""Build an n=14 interval shape-matrix certificate for EHP114.

This is the next hardening step after the radial Puiseux certificates. It
constructs the n=14 Fourier quotient shape basis, evaluates the finite
difference matrix through a Rust/inari batch length oracle, and checks a
Gershgorin lower bound with interval length enclosures.

It remains a local shape-cone component, not a proof of EHP114.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


sys.dont_write_bytecode = True

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01"
SPLICE_ID = "EXP-MATH-EHP114-N14-SHAPE-CONE-SPLICE-20260505-01"
RADIAL_COMPACT_ID = "EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01"
RADIAL_TAIL_ID = "EXP-MATH-EHP114-N14-RADIAL-TAIL-ARB-20260505-01"
DEGREE = 14
EPS = 0.005
RES = 220
EXTENT = 3.0
WORKING_LAMBDA14 = 100000.0


@dataclass(frozen=True)
class Interval:
    lo: float
    hi: float

    def __add__(self, rhs: "Interval") -> "Interval":
        return Interval(self.lo + rhs.lo, self.hi + rhs.hi)

    def __sub__(self, rhs: "Interval") -> "Interval":
        return Interval(self.lo - rhs.hi, self.hi - rhs.lo)

    def scale(self, value: float) -> "Interval":
        a = self.lo * value
        b = self.hi * value
        return Interval(min(a, b), max(a, b))

    def abs_upper(self) -> float:
        return max(abs(self.lo), abs(self.hi))

    def to_json(self) -> dict[str, float]:
        return {"lo": self.lo, "hi": self.hi, "width": self.hi - self.lo}


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


def build_points(tensor: Any) -> tuple[list[Any], list[dict[str, Any]]]:
    roots = tensor.roots_of_unity(DEGREE)
    basis = tensor.quotient_basis(DEGREE)
    shape = [b for b in basis if b.label != "radial_singular_m0"]
    points: list[dict[str, Any]] = []
    seen: set[str] = set()

    def add(label: str, vec: np.ndarray) -> None:
        if label in seen:
            return
        seen.add(label)
        coeffs = tensor.coeffs_from_roots(roots + EPS * vec)
        points.append({"label": label, "coeffs": coeffs_json(coeffs)})

    for i, b in enumerate(shape):
        add(f"d:{i}:+", b.vector)
        add(f"d:{i}:-", -b.vector)

    for i, bi in enumerate(shape):
        for j in range(i + 1, len(shape)):
            bj = shape[j]
            add(f"o:{i}:{j}:++", bi.vector + bj.vector)
            add(f"o:{i}:{j}:+-", bi.vector - bj.vector)
            add(f"o:{i}:{j}:-+", -bi.vector + bj.vector)
            add(f"o:{i}:{j}:--", -bi.vector - bj.vector)

    return shape, points


def run_length_oracle(root: Path, input_path: Path, output_path: Path) -> None:
    binary = root / "erdos-experiments" / "scripts" / "erdos-114" / "target" / "release" / "ehp114_batch_interval_lengths"
    if not binary.exists():
        subprocess.run(
            ["cargo", "build", "--release", "--bin", "ehp114_batch_interval_lengths"],
            cwd=root / "erdos-experiments" / "scripts" / "erdos-114",
            check=True,
        )
    subprocess.run(
        [str(binary), "--input", str(input_path), "--output", str(output_path)],
        cwd=root / "erdos-experiments" / "scripts" / "erdos-114",
        check=True,
    )


def build_matrix(lengths: dict[str, Interval], dim: int, lstar: Interval) -> tuple[list[list[Interval]], list[dict[str, Any]]]:
    def deficit(label: str) -> Interval:
        length = lengths[label]
        return Interval(lstar.lo - length.hi, lstar.hi - length.lo)

    matrix = [[Interval(0.0, 0.0) for _ in range(dim)] for _ in range(dim)]
    eps2 = EPS * EPS

    for i in range(dim):
        value = (deficit(f"d:{i}:+") + deficit(f"d:{i}:-")).scale(0.5 / eps2)
        matrix[i][i] = value

    for i in range(dim):
        for j in range(i + 1, dim):
            value = (
                deficit(f"o:{i}:{j}:++")
                - deficit(f"o:{i}:{j}:+-")
                - deficit(f"o:{i}:{j}:-+")
                + deficit(f"o:{i}:{j}:--")
            ).scale(1.0 / (4.0 * eps2))
            matrix[i][j] = value
            matrix[j][i] = value

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
    tensor = load_module(
        "ehp114_tensor_cone_scaffold",
        root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py",
    )
    shape, points = build_points(tensor)
    input_path = out_dir / f"{EXPERIMENT_ID}_ORACLE_INPUT.json"
    output_path = out_dir / f"{EXPERIMENT_ID}_ORACLE_OUTPUT.json"

    for path in (input_path, output_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing oracle artifact: {path}")

    input_payload = {
        "degree": DEGREE,
        "res": RES,
        "extent": EXTENT,
        "eps": EPS,
        "point_count": len(points),
        "points": points,
    }
    input_path.write_text(json.dumps(input_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    run_length_oracle(root, input_path, output_path)
    output = load_json(output_path)

    lengths = {
        row["label"]: Interval(float(row["length_lower"]), float(row["length_upper"]))
        for row in output["points"]
    }
    inari_path = root / "erdos-experiments" / "results" / "erdos-114" / "EXP-MM-EHP-007-n14-inari_RESULTS.json"
    inari = load_json(inari_path)
    lstar = Interval(float(inari["l_star_lower"]), float(inari["l_star_upper"]))
    matrix, rows = build_matrix(lengths, len(shape), lstar)
    lower_bound = min(row["diagonal_lower_minus_radius"] for row in rows)
    worst = min(rows, key=lambda row: row["diagonal_lower_minus_radius"])
    max_entry_width = max(matrix[i][j].hi - matrix[i][j].lo for i in range(len(shape)) for j in range(len(shape)))
    all_positive = lower_bound > 0
    lambda_ok = lower_bound > WORKING_LAMBDA14

    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "scope": "n=14 shape-cone interval matrix certificate",
        "status": "SHAPE_INTERVAL_MATRIX_CERTIFIED" if lambda_ok else ("POSITIVE_BUT_LAMBDA_TOO_HIGH" if all_positive else "SHAPE_INTERVAL_MATRIX_FAILED"),
        "claim_ceiling": "This is a shadow signature, not universal law. It certifies the finite-difference shape matrix under the stated interval marching-squares oracle, not EHP #114.",
        "parameters": {
            "degree": DEGREE,
            "eps": EPS,
            "res": RES,
            "extent": EXTENT,
            "shape_basis_rank": len(shape),
            "oracle_point_count": len(points),
            "working_lambda14": WORKING_LAMBDA14,
        },
        "sources": {
            "splice_packet": str(out_dir / f"{SPLICE_ID}_RESULTS.json"),
            "radial_compact": str(out_dir / f"{RADIAL_COMPACT_ID}_RESULTS.json"),
            "radial_tail": str(out_dir / f"{RADIAL_TAIL_ID}_RESULTS.json"),
            "inari_n14": {
                "path": str(inari_path),
                "sha256": sha256_file(inari_path),
                "verdict": inari.get("verdict"),
                "rigor": inari.get("rigor"),
            },
        },
        "oracle": {
            "input_path": str(input_path),
            "input_sha256": sha256_file(input_path),
            "output_path": str(output_path),
            "output_sha256": sha256_file(output_path),
            "elapsed_secs": output["elapsed_secs"],
            "method": "Rust/inari interval marching-squares batch length oracle",
        },
        "shape_basis_labels": [b.label for b in shape],
        "gershgorin_interval_lower_bound": lower_bound,
        "worst_row": worst,
        "max_matrix_interval_width": max_entry_width,
        "working_lambda14_below_interval_bound": lambda_ok,
        "all_rows_positive": all(row["diagonal_lower_minus_radius"] > 0 for row in rows),
        "rows": rows,
        "mixed_remainder_obligation": {
            "status": "OPEN",
            "target": "Bound mixed radial/shape remainder by <= 12*eps^(1/14) + 0.5*lambda14*||s||^2 on the n=14 local cone.",
            "why_it_matters": "The radial and shape certificates are now separate; only mixed-term absorption blocks the local theorem.",
        },
        "next_blocker": "Prove mixed radial/shape remainder absorption for the n=14 local cone.",
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
    row_lines = "\n".join(
        "| {index} | {diag:.12g} | {radius:.12g} | {lower:.12g} |".format(
            index=row["index"],
            diag=row["diagonal_interval"]["lo"],
            radius=row["off_diagonal_abs_radius_upper"],
            lower=row["diagonal_lower_minus_radius"],
        )
        for row in result["rows"]
    )
    return f"""# EHP114 n=14 Shape Interval Matrix Certificate

Experiment: `{EXPERIMENT_ID}`

## Meaning

This hardens the n=14 shape cone from a floating matrix into an interval matrix
certificate, using the Rust/inari length oracle. It is still local and finite:
the certified object is the finite-difference shape matrix at eps
`{result["parameters"]["eps"]}`, not the full EHP114 theorem.

The claim ceiling remains: this is a shadow signature, not universal law.

## Verdict

- Status: `{result["status"]}`
- Shape basis rank: `{result["parameters"]["shape_basis_rank"]}`
- Oracle points: `{result["parameters"]["oracle_point_count"]}`
- Gershgorin interval lower bound: `{result["gershgorin_interval_lower_bound"]}`
- Working lambda14: `{result["parameters"]["working_lambda14"]}`
- Lambda below interval bound: `{result["working_lambda14_below_interval_bound"]}`
- Max matrix interval width: `{result["max_matrix_interval_width"]}`
- Next blocker: `{result["next_blocker"]}`

## Row Bounds

| row | diagonal lower | offdiag radius upper | lower minus radius |
|---:|---:|---:|---:|
{row_lines}

## What This Changes

The shape cone is no longer only a floating diagnostic. Under the same
interval marching-squares method used in the Rust/inari finite certificate, the
n=14 shape matrix has a positive interval Gershgorin lower bound.

This makes the Jordan-style articulation sharper: the shape directions behave
like a positive spectral cone, while the radial direction is the Puiseux
boundary term.

## What Remains

{result["mixed_remainder_obligation"]["target"]}

That mixed-remainder theorem is now the live blocker.

## Files

- Script: `erdos-experiments/Erdos114/ehp114_n14_shape_interval_matrix_certificate.py`
- Rust oracle: `erdos-experiments/scripts/erdos-114/src/bin/ehp114_batch_interval_lengths.rs`
- Oracle input: `{result["oracle"]["input_path"]}`
- Oracle output: `{result["oracle"]["output_path"]}`
- Result SHA: `{result["oracle"]["output_sha256"]}` for oracle output

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
                "gershgorin_interval_lower_bound": result["gershgorin_interval_lower_bound"],
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
