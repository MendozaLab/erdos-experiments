#!/usr/bin/env python3
"""Build a machine-readable row contract for the EHP114 n=14 bridge candidate.

This is a certificate-shaping script, not a proof checker for Erdős #114. It
extracts the retained Taylor-collar rows from the Rust bridge diagnostic and
checks the local numerical obligations that were used to get the normal-drift
candidate under budget.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-N14-TAYLOR-COLLAR-ROW-CONTRACT-20260505-01"
SOURCE_RESULT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "erdos-114"
    / "bridge-diagnostic-taylor-collar-local-bound6-output-z16"
    / "EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01_RESULTS.json"
)
DEFAULT_OUTDIR = Path(__file__).resolve().parent / EXPERIMENT_ID


def finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def contract_row(row_id: int, row: dict[str, Any]) -> dict[str, Any]:
    normal_error = row.get("normal_drift_error_candidate")
    return {
        "row_id": row_id,
        "ix": row["ix"],
        "iy": row["iy"],
        "x_interval": row["x_interval"],
        "y_interval": row["y_interval"],
        "regularity_resolved": row["regularity_resolved"],
        "p_abs_upper_local": row["p_abs_upper_local"],
        "p_prime_abs_lower_taylor": row["p_prime_abs_lower_taylor"],
        "p_prime_abs_upper_local": row["p_prime_abs_upper_local"],
        "p_second_abs_upper_local": row["p_second_abs_upper_local"],
        "gradient_lower": row["gradient_lower_on_level_candidate"],
        "hessian_upper": row["hessian_spectral_upper_candidate"],
        "condition_ratio": row["condition_ratio_candidate"],
        "normal_drift_error": normal_error,
        "local_bound_subdivision": row["local_bound_subdivision"],
    }


def build_contract(source_path: Path) -> dict[str, Any]:
    source = load_json(source_path)
    cells = source["diagnostic_cells"]
    budget = float(source["source_margins"]["available_error_budget"])

    rows = [contract_row(i, row) for i, row in enumerate(cells)]
    row_failures: list[dict[str, Any]] = []
    normal_error_sum = 0.0

    for row in rows:
        normal_error = row["normal_drift_error"]
        checks = {
            "regularity_resolved": row["regularity_resolved"] is True,
            "positive_gradient": finite_number(row["gradient_lower"]) and row["gradient_lower"] > 0.0,
            "finite_hessian": finite_number(row["hessian_upper"]) and row["hessian_upper"] >= 0.0,
            "finite_condition_ratio": finite_number(row["condition_ratio"]) and row["condition_ratio"] >= 0.0,
            "finite_normal_error": finite_number(normal_error) and normal_error >= 0.0,
            "local_bound_subdivision": row["local_bound_subdivision"] == 6,
        }
        if not all(checks.values()):
            row_failures.append({"row_id": row["row_id"], "checks": checks})
        if finite_number(normal_error):
            normal_error_sum += float(normal_error)

    row_contract_bytes = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("utf-8")
    row_digest = hashlib.sha256(row_contract_bytes).hexdigest()
    source_sha = sha256_path(source_path)
    slack = budget - normal_error_sum

    accepted = (
        source["status"] == "BRIDGE_CANDIDATE_ERROR_UNDER_BUDGET_NOT_THEOREM"
        and source["regularity_unresolved_cell_count"] == 0
        and not row_failures
        and normal_error_sum <= budget
    )

    return {
        "experiment_id": EXPERIMENT_ID,
        "source_result": str(source_path),
        "source_result_sha256": source_sha,
        "source_experiment_id": source["experiment_id"],
        "status": "ROW_CONTRACT_ACCEPTED_NOT_THEOREM" if accepted else "ROW_CONTRACT_FAILED",
        "claim_ceiling": (
            "Taylor-collar normal-drift row contract only. This is not an exact "
            "lemniscate-length certificate, not a Lean theorem, and not a proof of Erdos #114."
        ),
        "parameters": {
            "degree": source["parameters"]["degree"],
            "eps": source["parameters"]["eps"],
            "sub_i": source["parameters"]["sub_i"],
            "sub_j": source["parameters"]["sub_j"],
            "z_subdivision": source["parameters"]["z_subdivision"],
            "regularity_strategy": source["regularity_strategy"],
            "derivative_mode": source["derivative_mode"],
            "local_bound_subdivision": source["parameters"]["local_bound_subdivision"],
        },
        "budget": {
            "normal_drift_budget": budget,
            "normal_drift_error_sum": normal_error_sum,
            "normal_drift_slack": slack,
            "normal_drift_error_to_budget": normal_error_sum / budget,
            "relative_length_error_to_budget": source["error_to_budget_ratio_relative"],
        },
        "row_summary": {
            "row_count": len(rows),
            "row_failure_count": len(row_failures),
            "row_contract_sha256": row_digest,
            "min_gradient_lower": min(row["gradient_lower"] for row in rows),
            "max_hessian_upper": max(row["hessian_upper"] for row in rows),
            "max_condition_ratio": max(row["condition_ratio"] for row in rows),
            "max_normal_drift_error": max(row["normal_drift_error"] for row in rows),
        },
        "row_failures": row_failures[:50],
        "rows": rows,
        "lean_shaped_target": (
            "For every retained z16 collar row B in the (6,4) root-affine "
            "parameter cell, the local 6x6 sub-bound encloses |p|, |p'|, "
            "and |p''|; Taylor gives |p'| > 0 on B; and the finite sum of "
            "normal-drift row errors is <= 2.5620009612530126."
        ),
        "next_blocker": (
            "Turn this row contract into a theorem-grade checker or Lean-side "
            "finite-sum certificate, then justify that the normal-drift bridge "
            "is the correct exact-length comparison theorem."
        ),
    }


def write_report(contract: dict[str, Any], path: Path) -> None:
    report = f"""# EHP114 n=14 Taylor-Collar Row Contract

Experiment: `{contract["experiment_id"]}`

## Verdict

- Status: `{contract["status"]}`
- Source result SHA-256: `{contract["source_result_sha256"]}`
- Row count: `{contract["row_summary"]["row_count"]}`
- Row failure count: `{contract["row_summary"]["row_failure_count"]}`
- Normal-drift error sum: `{contract["budget"]["normal_drift_error_sum"]}`
- Normal-drift budget: `{contract["budget"]["normal_drift_budget"]}`
- Normal-drift slack: `{contract["budget"]["normal_drift_slack"]}`
- Normal error / budget: `{contract["budget"]["normal_drift_error_to_budget"]}`
- Relative error / budget: `{contract["budget"]["relative_length_error_to_budget"]}`
- Row contract SHA-256: `{contract["row_summary"]["row_contract_sha256"]}`

## Claim Ceiling

This is a Taylor-collar normal-drift row contract only. It is not an exact
lemniscate-length certificate, not a Lean theorem, and not a proof of Erdős
#114.

## Lean-Shaped Target

```text
{contract["lean_shaped_target"]}
```

## Next Blocker

{contract["next_blocker"]}
"""
    path.write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=SOURCE_RESULT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    if not args.source.exists():
        raise SystemExit(f"missing source result: {args.source}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    result_path = args.outdir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = args.outdir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = args.outdir / f"{EXPERIMENT_ID}_RESULTS.sha256"

    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise SystemExit(f"refusing to overwrite existing artifact: {path}")

    contract = build_contract(args.source.resolve())
    result_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(contract, report_path)
    result_sha = sha256_path(result_path)
    sha_path.write_text(f"{result_sha}  {result_path.name}\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "experiment_id": EXPERIMENT_ID,
                "status": contract["status"],
                "row_count": contract["row_summary"]["row_count"],
                "row_failure_count": contract["row_summary"]["row_failure_count"],
                "normal_drift_error_to_budget": contract["budget"]["normal_drift_error_to_budget"],
                "normal_drift_slack": contract["budget"]["normal_drift_slack"],
                "result": str(result_path),
                "report": str(report_path),
                "sha256": result_sha,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

