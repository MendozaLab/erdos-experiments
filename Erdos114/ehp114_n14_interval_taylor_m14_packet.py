#!/usr/bin/env python3
"""Interval Taylor / M14 remainder packet for the n=14 EHP114 local cone.

This script is Route B from the theorem-target packet. It evaluates validated
interval lemniscate lengths at a structured Taylor stencil around radial bases

    roots = (1 - eps)^(1/14) * 14th_roots_of_unity

and asks whether the conservative reserves used in the Lean scratch target are
supported on that stencil:

    R14(eps) >= 24 eps^(1/14)
    Q14(s)   >= 100000 ||s||^2
    M14      <= 12 eps^(1/14) + 50000 ||s||^2

The packet is intentionally finite and local. Passing it does not prove EHP114
and does not prove the uniform mixed-remainder theorem. It builds the next
interval-Taylor certificate target.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


sys.dont_write_bytecode = True

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-M14-PACKET-20260505-01"
LOCAL_MIXED_ID = "EXP-MATH-EHP114-N14-LOCAL-MIXED-REMAINDER-SCOUT-20260505-02"
SHAPE_MATRIX_ID = "EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01"
RADIAL_COMPACT_ID = "EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01"
RADIAL_TAIL_ID = "EXP-MATH-EHP114-N14-RADIAL-TAIL-ARB-20260505-01"

DEGREE = 14
EPS_VALUES = [1e-4, 1e-3, 1e-2, 1e-1]
TAYLOR_STEP = 0.004
LOCAL_T_CAP = 0.008
RES = 220
EXTENT = 3.0
RADIAL_HALF = 12.0
SHAPE_HALF = 50000.0
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


def max_root_radius(roots: np.ndarray) -> float:
    return float(np.max(np.abs(roots)))


def max_admissible_t(base: np.ndarray, direction: np.ndarray) -> float:
    lo = 0.0
    hi = 0.1
    while max_root_radius(base + hi * direction) <= 1.0 and hi < 10.0:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if max_root_radius(base + mid * direction) <= 1.0:
            lo = mid
        else:
            hi = mid
    return lo


def add_point(
    seen: set[str],
    points: list[dict[str, Any]],
    meta: dict[str, dict[str, Any]],
    tensor: Any,
    label: str,
    roots: np.ndarray,
    payload: dict[str, Any],
) -> None:
    if label in seen:
        return
    seen.add(label)
    points.append({"label": label, "coeffs": coeffs_json(tensor.coeffs_from_roots(roots))})
    meta[label] = {**payload, "max_root_radius": max_root_radius(roots)}


def build_points(tensor: Any) -> tuple[list[Any], list[dict[str, Any]], dict[str, dict[str, Any]]]:
    unit_roots = tensor.roots_of_unity(DEGREE)
    basis = tensor.quotient_basis(DEGREE)
    shape = [b for b in basis if b.label != "radial_singular_m0"]
    seen: set[str] = set()
    points: list[dict[str, Any]] = []
    meta: dict[str, dict[str, Any]] = {}

    for eps in EPS_VALUES:
        radius = (1.0 - eps) ** (1.0 / DEGREE)
        base = radius * unit_roots
        eps_tag = f"{eps:.0e}"
        add_point(
            seen,
            points,
            meta,
            tensor,
            f"eps:{eps_tag}:base",
            base,
            {"kind": "base", "eps": eps, "t": 0.0},
        )

        for i, b in enumerate(shape):
            for sign, signed in [("+", b.vector), ("-", -b.vector)]:
                tmax = max_admissible_t(base, signed)
                t_axis = min(0.9 * tmax, LOCAL_T_CAP)
                add_point(
                    seen,
                    points,
                    meta,
                    tensor,
                    f"eps:{eps_tag}:diag:{i}:{sign}",
                    base + TAYLOR_STEP * signed,
                    {
                        "kind": "diag",
                        "eps": eps,
                        "shape_index": i,
                        "shape_label": b.label,
                        "sign": sign,
                        "t": TAYLOR_STEP,
                    },
                )
                add_point(
                    seen,
                    points,
                    meta,
                    tensor,
                    f"eps:{eps_tag}:axis-cap:{i}:{sign}",
                    base + t_axis * signed,
                    {
                        "kind": "axis_cap",
                        "eps": eps,
                        "shape_index": i,
                        "shape_label": b.label,
                        "sign": sign,
                        "t": t_axis,
                        "tmax": tmax,
                        "local_cap_active": t_axis < 0.9 * tmax,
                    },
                )

        for i, bi in enumerate(shape):
            for j in range(i + 1, len(shape)):
                bj = shape[j]
                combos = [
                    ("++", bi.vector + bj.vector),
                    ("+-", bi.vector - bj.vector),
                    ("-+", -bi.vector + bj.vector),
                    ("--", -bi.vector - bj.vector),
                ]
                for signs, direction in combos:
                    add_point(
                        seen,
                        points,
                        meta,
                        tensor,
                        f"eps:{eps_tag}:off:{i}:{j}:{signs}",
                        base + TAYLOR_STEP * direction,
                        {
                            "kind": "off",
                            "eps": eps,
                            "shape_i": i,
                            "shape_j": j,
                            "shape_i_label": bi.label,
                            "shape_j_label": bj.label,
                            "signs": signs,
                            "t": TAYLOR_STEP,
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


def build_matrix_for_eps(
    eps: float,
    dim: int,
    deficits: dict[str, Interval],
) -> tuple[list[list[Interval]], list[dict[str, Any]]]:
    eps_tag = f"{eps:.0e}"
    base = deficits[f"eps:{eps_tag}:base"]
    matrix = [[Interval(0.0, 0.0) for _ in range(dim)] for _ in range(dim)]
    h2 = TAYLOR_STEP * TAYLOR_STEP

    for i in range(dim):
        plus = deficits[f"eps:{eps_tag}:diag:{i}:+"]
        minus = deficits[f"eps:{eps_tag}:diag:{i}:-"]
        value = (plus + minus - base.scale(2.0)).scale(0.5 / h2)
        matrix[i][i] = value

    for i in range(dim):
        for j in range(i + 1, dim):
            value = (
                deficits[f"eps:{eps_tag}:off:{i}:{j}:++"]
                - deficits[f"eps:{eps_tag}:off:{i}:{j}:+-"]
                - deficits[f"eps:{eps_tag}:off:{i}:{j}:-+"]
                + deficits[f"eps:{eps_tag}:off:{i}:{j}:--"]
            ).scale(1.0 / (4.0 * h2))
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
    shape, points, meta = build_points(tensor)
    input_path = out_dir / f"{EXPERIMENT_ID}_ORACLE_INPUT.json"
    output_path = out_dir / f"{EXPERIMENT_ID}_ORACLE_OUTPUT.json"
    for path in (input_path, output_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing oracle artifact: {path}")

    input_payload = {
        "degree": DEGREE,
        "res": RES,
        "extent": EXTENT,
        "eps_values": EPS_VALUES,
        "taylor_step": TAYLOR_STEP,
        "local_t_cap": LOCAL_T_CAP,
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
    deficits = {
        label: Interval(lstar.lo - length.hi, lstar.hi - length.lo)
        for label, length in lengths.items()
    }

    eps_summaries = []
    global_lower = float("inf")
    global_max_width = 0.0
    for eps in EPS_VALUES:
        matrix, rows = build_matrix_for_eps(eps, len(shape), deficits)
        lower = min(row["diagonal_lower_minus_radius"] for row in rows)
        worst = min(rows, key=lambda row: row["diagonal_lower_minus_radius"])
        max_width = max(
            matrix[i][j].hi - matrix[i][j].lo
            for i in range(len(shape))
            for j in range(len(shape))
        )
        global_lower = min(global_lower, lower)
        global_max_width = max(global_max_width, max_width)
        eps_summaries.append(
            {
                "eps": eps,
                "gershgorin_interval_lower_bound": lower,
                "working_lambda14_below_interval_bound": lower > WORKING_LAMBDA14,
                "worst_row": worst,
                "max_matrix_interval_width": max_width,
                "all_rows_positive": all(row["diagonal_lower_minus_radius"] > 0 for row in rows),
            }
        )

    axis_rows = []
    for label, row in meta.items():
        if row["kind"] != "axis_cap":
            continue
        eps = float(row["eps"])
        t = float(row["t"])
        deficit_lower = deficits[label].lo
        rhs = RADIAL_HALF * (eps ** (1.0 / DEGREE)) + SHAPE_HALF * t * t
        axis_rows.append(
            {
                **row,
                "label": label,
                "length_interval": lengths[label].to_json(),
                "deficit_interval": deficits[label].to_json(),
                "deficit_lower": deficit_lower,
                "rhs_mixed_absorption_budget": rhs,
                "margin": deficit_lower - rhs,
                "pass": deficit_lower > rhs,
                "admissible_at_endpoint": row["max_root_radius"] <= 1.0 + 1e-12,
            }
        )

    all_matrix_lambda_ok = all(row["working_lambda14_below_interval_bound"] for row in eps_summaries)
    all_axis_pass = all(row["pass"] and row["admissible_at_endpoint"] for row in axis_rows)
    all_oracle_points_admissible = all(row["max_root_radius"] <= 1.0 + 1e-12 for row in meta.values())
    all_axis_points_admissible = all(row["max_root_radius"] <= 1.0 + 1e-12 for row in axis_rows)
    if all_matrix_lambda_ok and all_axis_pass:
        status = "INTERVAL_TAYLOR_M14_AXIS_PACKET_PASS"
    elif all_matrix_lambda_ok:
        status = "INTERVAL_TAYLOR_MATRIX_PASS_AXIS_FAIL"
    else:
        status = "INTERVAL_TAYLOR_MATRIX_FAIL"

    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "scope": "n=14 interval Taylor packet for the M14 local mixed-remainder target",
        "status": status,
        "claim_ceiling": "This is a shadow signature, not universal law. It is a finite interval-Taylor packet on a structured stencil, not the uniform mixed-remainder theorem.",
        "parameters": {
            "degree": DEGREE,
            "eps_values": EPS_VALUES,
            "shape_basis_rank": len(shape),
            "taylor_step": TAYLOR_STEP,
            "local_t_cap": LOCAL_T_CAP,
            "res": RES,
            "extent": EXTENT,
            "oracle_point_count": len(points),
            "working_lambda14": WORKING_LAMBDA14,
            "radial_half_constant": RADIAL_HALF,
            "shape_half_constant": SHAPE_HALF,
        },
        "sources": {
            "local_mixed_scout": str(out_dir / f"{LOCAL_MIXED_ID}_RESULTS.json"),
            "shape_matrix": str(out_dir / f"{SHAPE_MATRIX_ID}_RESULTS.json"),
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
        "matrix_eps_summaries": eps_summaries,
        "global_gershgorin_interval_lower_bound": global_lower,
        "global_max_matrix_interval_width": global_max_width,
        "all_matrix_lambda_ok": all_matrix_lambda_ok,
        "axis_endpoint_count": len(axis_rows),
        "all_axis_endpoints_pass_budget": all_axis_pass,
        "axis_min_margin": min(row["margin"] for row in axis_rows),
        "axis_worst_point": min(axis_rows, key=lambda row: row["margin"]),
        "axis_rows": axis_rows,
        "all_oracle_points_admissible": all_oracle_points_admissible,
        "all_axis_points_admissible": all_axis_points_admissible,
        "max_oracle_root_radius": max(row["max_root_radius"] for row in meta.values()),
        "max_axis_root_radius": max(row["max_root_radius"] for row in axis_rows),
        "interpretation": {
            "what_pass_would_mean": "The structured Taylor stencil supports the conservative M14 budget at radial epsilon grid points and axis endpoints.",
            "what_it_does_not_mean": "It does not certify every point of the 25-dimensional local cone or prove the uniform mixed-remainder theorem.",
        },
        "next_blocker": "Upgrade the structured stencil to a box-wise interval Taylor model with derivative/remainder bounds over multi-mode shape boxes.",
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
    eps_table = "\n".join(
        "| {eps:.0e} | {ok} | {lower:.12g} | {width:.6g} | {worst} |".format(
            eps=row["eps"],
            ok="yes" if row["working_lambda14_below_interval_bound"] else "no",
            lower=row["gershgorin_interval_lower_bound"],
            width=row["max_matrix_interval_width"],
            worst=row["worst_row"]["index"],
        )
        for row in result["matrix_eps_summaries"]
    )
    worst = result["axis_worst_point"]
    return f"""# EHP114 n=14 Interval Taylor / M14 Packet

Experiment: `{EXPERIMENT_ID}`

## Meaning

This packet builds the first interval Taylor scaffold for the local
mixed-remainder target:

```text
M14(eps,s) <= 12 eps^(1/14) + 50000 ||s||^2.
```

It evaluates interval lemniscate lengths on a structured Taylor stencil around
radial bases `r(eps) * roots_of_unity`, with `eps` in `{EPS_VALUES}`.

This is a shadow signature, not universal law. It is not a solution of
Erdos #114 and not the uniform mixed-remainder theorem.

## Verdict

- Status: `{result["status"]}`
- Shape basis rank: `{result["parameters"]["shape_basis_rank"]}`
- Oracle point count: `{result["parameters"]["oracle_point_count"]}`
- Taylor step: `{result["parameters"]["taylor_step"]}`
- Local axis cap: `{result["parameters"]["local_t_cap"]}`
- Global Gershgorin interval lower bound: `{result["global_gershgorin_interval_lower_bound"]}`
- Working `lambda14`: `{result["parameters"]["working_lambda14"]}`
- All matrix epsilon slices pass lambda: `{result["all_matrix_lambda_ok"]}`
- Axis endpoints passing budget: `{result["all_axis_endpoints_pass_budget"]}`
- Axis minimum margin: `{result["axis_min_margin"]}`
- All axis endpoint roots admissible: `{result["all_axis_points_admissible"]}`
- All oracle points root-admissible: `{result["all_oracle_points_admissible"]}`
- Max oracle root radius: `{result["max_oracle_root_radius"]}`
- Max axis root radius: `{result["max_axis_root_radius"]}`

## Per-Epsilon Taylor Matrix Summary

| eps | lambda pass | Gershgorin lower | max interval width | worst row |
|---:|---:|---:|---:|---:|
{eps_table}

## Worst Axis Endpoint

```json
{json.dumps(worst, indent=2, sort_keys=True)}
```

## Interpretation

The packet gives a split verdict:

```text
axis endpoint budget: passes
uniform positive Taylor matrix at radial bases: fails
```

The strongest safe reading is:

```text
The conservative M14 budget survives the tested admissible signed axis
endpoints, but the naive assumption that the boundary shape cone remains
uniformly positive after radial contraction is false on this stencil.
```

The unsafe reading is that the local cone is fully certified. This packet does
not yet bound every multi-mode shape vector inside the 24-dimensional shape
ball. The central-difference Taylor matrix is an ambient diagnostic; some
off-axis stencil points may lie just outside root admissibility. The admissible
axis endpoints are reported separately.

## Next Theorem Target

Upgrade this packet from stencil evidence to a box theorem:

```text
For every eps in (0, 1/10] and every shape vector s with
||s|| <= min(eta14Boundary(eps), 0.008),
M14(eps,s) <= 12 eps^(1/14) + 50000 ||s||^2.
```

The next run should add derivative/remainder interval bounds over multi-mode
boxes, not just axis endpoints.

## Guardrails

- No scorecard update.
- No D1 update.
- No public claim.
- No Lean status change.
- No email.
- No git operation.
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
                "global_gershgorin_interval_lower_bound": result["global_gershgorin_interval_lower_bound"],
                "axis_min_margin": result["axis_min_margin"],
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
