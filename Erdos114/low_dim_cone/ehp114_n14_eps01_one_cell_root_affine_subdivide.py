#!/usr/bin/env python3
"""Subdivide the selected eps=0.1 cell for root-affine interval certification."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


sys.dont_write_bytecode = True

PARENT_ID = "EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-ORACLE-20260505-01"
SUBDIV = int(os.environ.get("EHP114_SUBDIV", "4"))
MAX_WORKERS = int(os.environ.get("EHP114_JOBS", str(min(os.cpu_count() or 1, 8))))
EXPERIMENT_ID = f"EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV{SUBDIV}-20260505-01"
EPS = 0.1
RADIAL_HALF = 12.0
DEGREE = 14
_WORKER_STATE: dict[str, Any] = {}


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


def math_root_from_script() -> Path:
    return Path(__file__).resolve().parents[3]


def scalar_target(eps: float) -> float:
    return float(RADIAL_HALF * eps ** (1.0 / DEGREE))


def split_interval(pair: list[float], n: int) -> list[list[float]]:
    lo, hi = float(pair[0]), float(pair[1])
    step = (hi - lo) / n
    return [[lo + i * step, lo + (i + 1) * step] for i in range(n)]


def init_worker(out_dir_str: str, root_str: str, u0: dict[str, Any], u1: dict[str, Any], lstar_lower: float, target: float) -> None:
    out_dir = Path(out_dir_str)
    root = Path(root_str)
    root_affine = load_module(
        f"ehp114_root_affine_prev_subdivide_worker_{os.getpid()}",
        out_dir / "ehp114_n14_eps01_one_cell_root_affine_oracle.py",
    )
    interval_prev = load_module(
        f"ehp114_interval_coeff_prev_subdivide_worker_{os.getpid()}",
        out_dir / "ehp114_n14_eps01_one_cell_interval_coeff_oracle.py",
    )
    tensor = load_module(
        f"ehp114_tensor_subdivide_worker_{os.getpid()}",
        root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py",
    )
    _WORKER_STATE.clear()
    _WORKER_STATE.update(
        {
            "root_affine": root_affine,
            "interval_prev": interval_prev,
            "tensor": tensor,
            "u0": u0,
            "u1": u1,
            "lstar_lower": lstar_lower,
            "target": target,
        }
    )


def evaluate_subcell_worker(task: tuple[int, int, list[float], list[float], dict[str, Any]]) -> dict[str, Any]:
    i, j, u0_pair, u1_pair, cell = task
    root_affine = _WORKER_STATE["root_affine"]
    interval_prev = _WORKER_STATE["interval_prev"]
    tensor = _WORKER_STATE["tensor"]
    subcell = {**cell, "u0_interval": u0_pair, "u1_interval": u1_pair, "sub_i": i, "sub_j": j}
    roots = root_affine.root_intervals_from_cell(interval_prev, tensor, _WORKER_STATE["u0"], _WORKER_STATE["u1"], subcell)
    oracle = root_affine.interval_marching_upper(interval_prev, roots)
    deficit_lower = _WORKER_STATE["lstar_lower"] - oracle["length_upper"]
    margin_lower = deficit_lower - _WORKER_STATE["target"]
    return {
        "sub_i": i,
        "sub_j": j,
        "u0_interval": u0_pair,
        "u1_interval": u1_pair,
        "length_upper": oracle["length_upper"],
        "deficit_lower": deficit_lower,
        "margin_lower": margin_lower,
        "pass": margin_lower >= 0.0,
        "active_cells": oracle["active_cells"],
        "definite_case_cells": oracle["definite_case_cells"],
        "uncertain_corner_cells": oracle["uncertain_corner_cells"],
        "ambiguous_avg_cells": oracle["ambiguous_avg_cells"],
    }


def build_result() -> dict[str, Any]:
    out_dir = Path(__file__).resolve().parent
    root = math_root_from_script()
    low_dim = load_module("ehp114_low_dim_subdivide", out_dir / "ehp114_n14_low_dim_cone_box_probe.py")
    parent = load_json(out_dir / f"{PARENT_ID}_RESULTS.json")
    cell = parent["selected_cell"]
    _, directions = low_dim.build_low_dim_directions()
    u0 = directions[EPS][0]
    u1 = directions[EPS][1]
    inari = load_json(root / "erdos-experiments" / "results" / "erdos-114" / "EXP-MM-EHP-007-n14-inari_RESULTS.json")
    lstar_lower = float(inari["l_star_lower"])
    target = scalar_target(EPS)

    tasks = [
        (i, j, u0_pair, u1_pair, cell)
        for i, u0_pair in enumerate(split_interval(cell["u0_interval"], SUBDIV))
        for j, u1_pair in enumerate(split_interval(cell["u1_interval"], SUBDIV))
    ]
    worker_args = (str(out_dir), str(root), u0, u1, lstar_lower, target)
    if MAX_WORKERS <= 1:
        init_worker(*worker_args)
        rows = [evaluate_subcell_worker(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=MAX_WORKERS, initializer=init_worker, initargs=worker_args) as pool:
            rows = list(pool.map(evaluate_subcell_worker, tasks))
    rows.sort(key=lambda row: (row["sub_i"], row["sub_j"]))
    failures = [row for row in rows if not row["pass"]]
    status = (
        f"ONE_CELL_ROOT_AFFINE_SUBDIV{SUBDIV}_PASS"
        if not failures
        else f"ONE_CELL_ROOT_AFFINE_SUBDIV{SUBDIV}_FAIL"
    )
    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "parent_packet": PARENT_ID,
        "status": status,
        "parameters": {
            "eps": EPS,
            "subdivision": SUBDIV,
            "max_workers": MAX_WORKERS,
            "subcell_count": len(rows),
            "scalar_target": "D14(eps,s) >= 12 * eps^(1/14)",
        },
        "selected_cell": cell,
        "lstar_lower": lstar_lower,
        "target": target,
        "failure_count": len(failures),
        "min_margin_lower": min(row["margin_lower"] for row in rows),
        "max_length_upper": max(row["length_upper"] for row in rows),
        "rows": rows,
        "claim_ceiling": (
            "Root-affine interval subdivision certificate for the current marching-squares oracle functional "
            "on one selected eps=0.1 coefficient cell. This is not exact lemniscate certification and not a proof of Erdős #114."
        ),
        "next_blocker": (
            "If pass, lift from marching-squares oracle certification to exact/validated lemniscate length; "
            "if fail, subdivide further or use analytic Cauchy bounds."
        ),
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    lines = [
        "# EHP114 n=14 eps=0.1 One-Cell Root-Affine Subdivision",
        "",
        f"Experiment: `{EXPERIMENT_ID}`",
        "",
        f"Parent packet: `{PARENT_ID}`",
        "",
        "## Verdict",
        "",
        f"- Status: `{result['status']}`",
        f"- Subcells: `{result['parameters']['subcell_count']}`",
        f"- Failure count: `{result['failure_count']}`",
        f"- Minimum margin lower bound: `{result['min_margin_lower']}`",
        f"- Maximum length upper bound: `{result['max_length_upper']}`",
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
                "failure_count": result["failure_count"],
                "min_margin_lower": result["min_margin_lower"],
                "max_length_upper": result["max_length_upper"],
                "report": str(report_path),
                "sha256": sha,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
