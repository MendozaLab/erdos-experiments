#!/usr/bin/env python3
"""
EXP-MATH-RH-BEURLING-NYMAN-MDL-QUANTIZATION-STABILITY-CHECK-20260506-01

Finite audit for the frozen Beurling-Nyman MDL objects. This checks that the
stored probe and sensitivity artifacts expose the quantization fields needed
for the modest stability theorem target. It is not RH evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-RH-BEURLING-NYMAN-MDL-QUANTIZATION-STABILITY-CHECK-20260506-01"
CLAIM_CEILING = (
    "Finite Beurling-Nyman MDL quantization-stability audit only; not RH "
    "evidence, not theorem progress, and not an asymptotic statement."
)
REQUIRED_Q_FIELDS = [
    "quantization_step",
    "quantization_penalty_l2_bound",
    "certified_upper_relative_residual",
    "relative_residual",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_sha(results_path: Path) -> dict[str, Any]:
    sha_path = results_path.with_name(results_path.name.replace("_RESULTS.json", "_RESULTS.sha256"))
    if not sha_path.exists():
        return {"path": str(results_path), "sha_path": str(sha_path), "status": "NO_SHA_SOURCE"}
    text = sha_path.read_text().strip().split()
    expected = text[0] if text else ""
    actual = sha256_file(results_path)
    return {
        "path": str(results_path),
        "sha_path": str(sha_path),
        "status": "PASS" if expected == actual else "FAIL",
        "expected": expected,
        "actual": actual,
    }


def finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def audit_result(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text())
    failures: list[dict[str, Any]] = []
    rows = data.get("rows", [])
    quantized_count = 0
    certified_monotone_fail_count = 0

    if "not RH evidence" not in str(data.get("claim_ceiling", "")):
        failures.append({"scope": "top_level", "reason": "claim_ceiling_missing_not_rh_evidence"})

    if not isinstance(rows, list) or not rows:
        failures.append({"scope": "top_level", "reason": "rows_missing_or_empty"})
        rows = []

    for row_index, row in enumerate(rows):
        if not finite_number(row.get("condition_number")) or float(row["condition_number"]) <= 0:
            failures.append({"row_index": row_index, "reason": "condition_number_missing_or_nonpositive"})
        quantized = row.get("quantized")
        if not isinstance(quantized, list) or not quantized:
            failures.append({"row_index": row_index, "reason": "quantized_array_missing_or_empty"})
            continue
        for q_index, q in enumerate(quantized):
            quantized_count += 1
            for field in REQUIRED_Q_FIELDS:
                if field not in q or not finite_number(q[field]):
                    failures.append(
                        {"row_index": row_index, "quantized_index": q_index, "reason": f"{field}_missing_or_nonfinite"}
                    )
            if all(field in q and finite_number(q[field]) for field in REQUIRED_Q_FIELDS):
                if float(q["quantization_step"]) <= 0:
                    failures.append(
                        {"row_index": row_index, "quantized_index": q_index, "reason": "quantization_step_nonpositive"}
                    )
                if float(q["quantization_penalty_l2_bound"]) < 0:
                    failures.append(
                        {"row_index": row_index, "quantized_index": q_index, "reason": "quantization_penalty_negative"}
                    )
                if float(q["certified_upper_relative_residual"]) + 1e-12 < float(q["relative_residual"]):
                    certified_monotone_fail_count += 1
                    failures.append(
                        {
                            "row_index": row_index,
                            "quantized_index": q_index,
                            "reason": "certified_upper_below_quantized_residual",
                        }
                    )

    return {
        "path": str(path),
        "experiment_id": data.get("experiment_id"),
        "row_count": len(rows),
        "quantized_entry_count": quantized_count,
        "failure_count": len(failures),
        "certified_monotone_fail_count": certified_monotone_fail_count,
        "failures": failures[:20],
        "claim_ceiling": data.get("claim_ceiling"),
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    report = f"""# {EXPERIMENT_ID} Report

## Verdict

- Status: `{result["status"]}`
- Source artifacts audited: `{len(result["source_audits"])}`
- Total rows audited: `{result["total_row_count"]}`
- Total quantized entries audited: `{result["total_quantized_entry_count"]}`
- Total audit failures: `{result["total_failure_count"]}`
- Source SHA status: `{result["source_sha_status"]}`

## Meaning

This is a finite audit for the Beurling-Nyman MDL artifacts. It checks that the
stored rows expose the fields needed for the quantization-stability theorem:
quantization step, quantization penalty, certified upper residual, and condition
number.

The `~3.2` bit bend remains named only as the first observed finite-N MDL
conditioning crossover. It is not treated as a constant, RH evidence, theorem
progress, or an asymptotic statement.

## Claim Ceiling

{CLAIM_CEILING}
"""
    path.write_text(report)


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--probe",
        type=Path,
        default=here / "EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-03_RESULTS.json",
    )
    parser.add_argument(
        "--sensitivity",
        type=Path,
        default=here / "EXP-MATH-RH-BEURLING-NYMAN-MDL-SENSITIVITY-20260506-01_RESULTS.json",
    )
    parser.add_argument("--outdir", type=Path, default=here)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result_path = args.outdir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = args.outdir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = args.outdir / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for path in [result_path, report_path, sha_path]:
        if path.exists():
            raise SystemExit(f"refusing to overwrite existing artifact: {path}")

    source_paths = [args.probe.resolve(), args.sensitivity.resolve()]
    sha_checks = [verify_sha(path) for path in source_paths]
    audits = [audit_result(path) for path in source_paths]
    total_failure_count = sum(audit["failure_count"] for audit in audits)
    source_sha_status = "PASS" if all(check["status"] == "PASS" for check in sha_checks) else "FAIL"
    status = (
        "BN_MDL_QUANTIZATION_STABILITY_CHECK_PASS_NOT_RH_EVIDENCE"
        if source_sha_status == "PASS" and total_failure_count == 0
        else "BN_MDL_QUANTIZATION_STABILITY_CHECK_FAIL"
    )
    result = {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_unix": int(time.time()),
        "status": status,
        "source_artifacts": [str(path) for path in source_paths],
        "source_sha_checks": sha_checks,
        "source_sha_status": source_sha_status,
        "source_audits": audits,
        "total_row_count": sum(audit["row_count"] for audit in audits),
        "total_quantized_entry_count": sum(audit["quantized_entry_count"] for audit in audits),
        "total_failure_count": total_failure_count,
        "theorem_target": "||A q_b(a)-y|| <= ||A a-y|| + sigma_max(A) sqrt(k) Delta_b / 2",
        "bend_name": "first observed finite-N MDL conditioning crossover",
        "claim_ceiling": CLAIM_CEILING,
    }
    result_path.write_text(json.dumps(result, indent=2) + "\n")
    write_report(result, report_path)
    digest = sha256_file(result_path)
    sha_path.write_text(f"{digest}  {result_path.name}\n")
    print(json.dumps({
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "total_row_count": result["total_row_count"],
        "total_quantized_entry_count": result["total_quantized_entry_count"],
        "total_failure_count": total_failure_count,
        "result": str(result_path),
        "report": str(report_path),
        "sha256": digest,
    }, indent=2))


if __name__ == "__main__":
    main()
