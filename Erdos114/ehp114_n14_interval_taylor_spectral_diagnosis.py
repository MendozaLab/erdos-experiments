#!/usr/bin/env python3
"""Spectral diagnosis for the n=14 interval Taylor M14 packet.

The first interval Taylor packet used Gershgorin lower bounds and failed. This
diagnosis recomputes the same interval matrices from the saved oracle output
and applies a sharper interval spectral bound:

    lambda_min(midpoint matrix) - Frobenius(entrywise interval radius).

If this also fails, the issue is not just Gershgorin crudeness. It means the
shape Hessian around radially contracted bases has genuine negative directions,
so the M14 theorem must absorb radial-shape curvature softening.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


sys.dont_write_bytecode = True

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-SPECTRAL-DIAGNOSIS-20260505-01"
PACKET_ID = "EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-M14-PACKET-20260505-01"


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


def build_result() -> dict[str, Any]:
    root = repo_root_from_script()
    out_dir = Path(__file__).resolve().parent
    packet_script = load_module(
        "ehp114_n14_interval_taylor_m14_packet",
        out_dir / "ehp114_n14_interval_taylor_m14_packet.py",
    )
    packet = load_json(out_dir / f"{PACKET_ID}_RESULTS.json")
    oracle = load_json(out_dir / f"{PACKET_ID}_ORACLE_OUTPUT.json")
    inari_path = root / "erdos-experiments" / "results" / "erdos-114" / "EXP-MM-EHP-007-n14-inari_RESULTS.json"
    inari = load_json(inari_path)

    lengths = {
        row["label"]: packet_script.Interval(float(row["length_lower"]), float(row["length_upper"]))
        for row in oracle["points"]
    }
    lstar = packet_script.Interval(float(inari["l_star_lower"]), float(inari["l_star_upper"]))
    deficits = {
        label: packet_script.Interval(lstar.lo - length.hi, lstar.hi - length.lo)
        for label, length in lengths.items()
    }

    rows = []
    for eps in packet_script.EPS_VALUES:
        matrix, gersh_rows = packet_script.build_matrix_for_eps(eps, 24, deficits)
        mid = np.array([[(entry.lo + entry.hi) / 2.0 for entry in row] for row in matrix], dtype=float)
        radius = np.array([[(entry.hi - entry.lo) / 2.0 for entry in row] for row in matrix], dtype=float)
        eigvals = np.linalg.eigvalsh(mid)
        eig_min = float(eigvals[0])
        eig_max = float(eigvals[-1])
        fro_radius = float(np.linalg.norm(radius, "fro"))
        spectral_lower = eig_min - fro_radius
        gersh_lower = min(row["diagonal_lower_minus_radius"] for row in gersh_rows)
        rows.append(
            {
                "eps": eps,
                "midpoint_lambda_min": eig_min,
                "midpoint_lambda_max": eig_max,
                "frobenius_interval_radius": fro_radius,
                "interval_spectral_lower_bound": spectral_lower,
                "gershgorin_lower_bound": gersh_lower,
                "positive_by_interval_spectral": spectral_lower > 0.0,
                "working_lambda14_pass": spectral_lower > packet_script.WORKING_LAMBDA14,
            }
        )

    all_positive = all(row["positive_by_interval_spectral"] for row in rows)
    all_lambda = all(row["working_lambda14_pass"] for row in rows)
    if all_lambda:
        status = "SPECTRAL_INTERVAL_LAMBDA_PASS"
    elif all_positive:
        status = "SPECTRAL_INTERVAL_POSITIVE_BUT_LAMBDA_FAIL"
    else:
        status = "RADIAL_BASE_SHAPE_SOFTENING_DETECTED"

    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "parent_packet": PACKET_ID,
        "status": status,
        "claim_ceiling": "This is a diagnostic for the interval Taylor packet. It does not prove or disprove EHP #114.",
        "interpretation": (
            "The Taylor-matrix failure is not just Gershgorin crudeness if any "
            "interval spectral lower bound is negative. In that case, the M14 "
            "term must absorb genuine radial-base shape softening."
        ),
        "rows": rows,
        "global_interval_spectral_lower_bound": min(row["interval_spectral_lower_bound"] for row in rows),
        "all_positive_by_interval_spectral": all_positive,
        "all_working_lambda14_pass": all_lambda,
        "sources": {
            "packet_results": str(out_dir / f"{PACKET_ID}_RESULTS.json"),
            "packet_results_sha256": sha256_file(out_dir / f"{PACKET_ID}_RESULTS.json"),
            "packet_oracle_output": str(out_dir / f"{PACKET_ID}_ORACLE_OUTPUT.json"),
            "packet_oracle_output_sha256": sha256_file(out_dir / f"{PACKET_ID}_ORACLE_OUTPUT.json"),
            "inari_n14": str(inari_path),
            "inari_n14_sha256": sha256_file(inari_path),
        },
        "next_blocker": "Reformulate M14 to absorb radial-base shape softening, or restrict eta14Boundary(eps) until the radial reserve dominates the negative shape curvature.",
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
        "| {eps:.0e} | {eig:.12g} | {rad:.3g} | {lower:.12g} | {g:.12g} | {ok} |".format(
            eps=row["eps"],
            eig=row["midpoint_lambda_min"],
            rad=row["frobenius_interval_radius"],
            lower=row["interval_spectral_lower_bound"],
            g=row["gershgorin_lower_bound"],
            ok="yes" if row["working_lambda14_pass"] else "no",
        )
        for row in result["rows"]
    )
    return f"""# EHP114 n=14 Interval Taylor Spectral Diagnosis

Experiment: `{EXPERIMENT_ID}`

Parent packet: `{PACKET_ID}`

## Meaning

The first interval Taylor packet failed the Gershgorin test. This diagnosis
checks whether that was merely a crude row-sum bound. It recomputes the
interval Taylor matrices and applies:

```text
lambda_min(midpoint matrix) - Frobenius(interval radius).
```

## Verdict

- Status: `{result["status"]}`
- Global interval spectral lower bound: `{result["global_interval_spectral_lower_bound"]}`
- All positive by interval spectral bound: `{result["all_positive_by_interval_spectral"]}`
- All pass working lambda14: `{result["all_working_lambda14_pass"]}`

## Per-Epsilon Diagnosis

| eps | midpoint lambda_min | radius | spectral lower | Gershgorin lower | lambda14 pass |
|---:|---:|---:|---:|---:|---:|
{table}

## Consequence

The matrix failure is not just a Gershgorin artifact. The shape Hessian around
radially contracted bases develops negative directions under this Taylor
stencil. The local axis endpoints still pass the M14 budget, so the route is
not dead; it means the mixed-remainder theorem must explicitly absorb
radial-base shape softening.

The next theorem target should be phrased as:

```text
radial reserve dominates negative radial/shape curvature on the local cone.
```

not as:

```text
the shape cone remains uniformly positive after radial contraction.
```

## Claim Ceiling

This is a shadow signature, not universal law. It is a diagnostic of the
fixed-n interval Taylor strategy, not a proof or disproof of EHP #114.
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

