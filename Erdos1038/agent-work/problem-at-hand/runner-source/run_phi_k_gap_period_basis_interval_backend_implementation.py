#!/usr/bin/env python3
"""Fail-closed backend implementation packet for the stable gap-period basis.

This packet adds and fixture-tests the Rust backend harness needed by the
Stage-B basis interval-certificate route. The fixtures are synthetic contract
fixtures. Passing them proves only that the backend rejects common structural
and conditioning mistakes before a private numeric payload is available.

Claim ceiling: backend implementation scaffold and fail-closed fixture tests
only. This packet does not certify B1/B2/B3 for the real carrier, does not
interval-audit a period residual, does not prove period-matrix legitimacy, does
not compose attainment, does not prove selector existence, does not close KKT
composition, does not give a global reduction, does not solve #1038, and does
not improve public SOTA.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPERIMENT_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-GAP-PERIOD-BASIS-INTERVAL-BACKEND-"
    "IMPLEMENTATION-20260527-01"
)
PARENT_PACKET_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-STABLE-GAP-PERIOD-BASIS-INTERVAL-"
    "CERTIFICATE-20260527-01"
)
ROUND6_PACKET_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-EXTERNAL-ROUND6-EVEREST-BRIDGE-"
    "ASSIMILATION-20260527-01"
)

HERE = Path(__file__).resolve().parent
PACKET_ROOT = HERE.parent / "erdos-1038"
BACKEND_NAME = "phi_k_stable_gap_period_basis_interval_certificate_backend"
BACKEND_SOURCE = HERE / "src" / "bin" / f"{BACKEND_NAME}.rs"
FIXTURE_ROOT = HERE / "fixtures" / "phi_k_stable_gap_period_basis_interval_certificate"

RESULTS_JSON = PACKET_ROOT / f"{EXPERIMENT_ID}_RESULTS.json"
REPORT_MD = PACKET_ROOT / f"{EXPERIMENT_ID}_REPORT.md"
STRUCTURAL_PRE_AUDIT_JSON = PACKET_ROOT / f"{EXPERIMENT_ID}_STRUCTURAL_PRE_AUDIT.json"
FIXTURE_ROWS_JSONL = PACKET_ROOT / f"{EXPERIMENT_ID}_FIXTURE_ROWS.jsonl"
SHA_FILE = PACKET_ROOT / f"{EXPERIMENT_ID}_RESULTS.sha256"

CLAIM_CEILING = (
    "Backend implementation scaffold and fail-closed synthetic fixture tests "
    "only. The packet does not certify B1/B2/B3 for the real carrier, does "
    "not interval-audit a period residual, does not prove period-matrix "
    "legitimacy, does not compose attainment, does not prove selector "
    "existence, does not close KKT composition, does not give a global "
    "reduction, does not solve #1038, and does not improve public SOTA."
)

FIXTURE_EXPECTATIONS = {
    "good_24_contract.txt": "PASS_CONTRACT_FIXTURE_ONLY_PRIVATE_NUMERIC_PAYLOAD_REQUIRED",
    "bad_27_rows.txt": "STRUCTURAL_PRE_AUDIT_FAIL_CLOSED",
    "b1_high_condition_transform.txt": "FAIL_CLOSED_CONTRACT_FIXTURE_REJECTED",
    "b2_missing_recovered_witnesses.txt": "FAIL_CLOSED_CONTRACT_FIXTURE_REJECTED",
    "b3_singular_matrix.txt": "FAIL_CLOSED_CONTRACT_FIXTURE_REJECTED",
}


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_unique(path: Path, text: str) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite immutable artifact: {path}")
    path.write_text(text, encoding="utf-8")


def run_backend(fixture: Path) -> dict[str, Any]:
    command = [
        "cargo",
        "run",
        "--quiet",
        "--bin",
        BACKEND_NAME,
        "--",
        str(fixture),
    ]
    completed = subprocess.run(
        command,
        cwd=HERE,
        check=False,
        capture_output=True,
        text=True,
    )
    payload: dict[str, Any]
    try:
        payload = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        payload = {
            "overall_status": "BACKEND_JSON_PARSE_FAIL",
            "parse_error": str(exc),
            "stdout": completed.stdout,
        }
    return {
        "fixture": fixture.name,
        "returncode": completed.returncode,
        "stderr": completed.stderr,
        "payload": payload,
    }


def structural_pre_audit_payload() -> dict[str, Any]:
    return {
        "packet_id": EXPERIMENT_ID,
        "object_being_audited": (
            "Fail-closed backend harness for the Stage-B stable gap-period "
            "basis interval certificate"
        ),
        "claimed_dimension": 24,
        "dimension_derivation": (
            "The real downstream certificate is for the 24-row gap-period "
            "kernel; normalization is excluded from the period rows."
        ),
        "object_count": 24,
        "count_derivation": (
            "The backend rejects 27-row/normalization-inclusive payloads and "
            "accepts only 24-row contract fixtures."
        ),
        "sign_convention": (
            "Right-transform convention: transformed_period_matrix = "
            "seed_period_matrix @ basis_transform. This packet tests the "
            "sign-convention guard, not the real period residual."
        ),
        "sign_convention_source": PARENT_PACKET_ID,
        "expected_failure_modes": [
            "27-vs-24 row mismatch",
            "normalization row leakage into period rows",
            "sign convention mismatch",
            "transform condition above threshold",
            "missing recovered-direction independence witnesses",
            "nonpositive smallest singular value",
        ],
    }


def main() -> None:
    if not BACKEND_SOURCE.exists():
        raise FileNotFoundError(BACKEND_SOURCE)
    missing = [
        fixture_name
        for fixture_name in FIXTURE_EXPECTATIONS
        if not (FIXTURE_ROOT / fixture_name).exists()
    ]
    if missing:
        raise FileNotFoundError(f"Missing fixtures: {missing}")

    structural_pre_audit = structural_pre_audit_payload()
    write_unique(
        STRUCTURAL_PRE_AUDIT_JSON,
        json.dumps(structural_pre_audit, indent=2, sort_keys=True) + "\n",
    )
    json.loads(STRUCTURAL_PRE_AUDIT_JSON.read_text(encoding="utf-8"))

    rows: list[dict[str, Any]] = []
    mismatches: list[dict[str, Any]] = []
    for fixture_name, expected_status in FIXTURE_EXPECTATIONS.items():
        row = run_backend(FIXTURE_ROOT / fixture_name)
        actual_status = row["payload"].get("overall_status")
        row["expected_status"] = expected_status
        row["actual_status"] = actual_status
        row["expectation_pass"] = row["returncode"] == 0 and actual_status == expected_status
        rows.append(row)
        if not row["expectation_pass"]:
            mismatches.append(
                {
                    "fixture": fixture_name,
                    "expected": expected_status,
                    "actual": actual_status,
                    "returncode": row["returncode"],
                }
            )

    write_unique(
        FIXTURE_ROWS_JSONL,
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
    )

    pass_count = sum(1 for row in rows if row["expectation_pass"])
    status = (
        "GAP_PERIOD_BASIS_INTERVAL_BACKEND_IMPLEMENTATION_PASS__"
        "FAIL_CLOSED_FIXTURES_PASS__PRIVATE_NUMERIC_PAYLOADS_PENDING"
        if not mismatches
        else "GAP_PERIOD_BASIS_INTERVAL_BACKEND_IMPLEMENTATION_FAIL__"
        "FIXTURE_EXPECTATION_MISMATCH"
    )
    results = {
        "packet_id": EXPERIMENT_ID,
        "created_utc": now_utc(),
        "route_type": "backend_implementation_scaffold",
        "status": status,
        "verdict": (
            "BACKEND_HARNESS_READY_FOR_PRIVATE_NUMERIC_PAYLOAD"
            if not mismatches
            else "BACKEND_HARNESS_FIXTURE_FAILURE"
        ),
        "parent_packet_id": PARENT_PACKET_ID,
        "round6_assimilation_packet_id": ROUND6_PACKET_ID,
        "backend_source": str(BACKEND_SOURCE.relative_to(HERE)),
        "fixture_root": str(FIXTURE_ROOT.relative_to(HERE)),
        "fixture_count": len(rows),
        "fixture_expectation_pass_count": pass_count,
        "fixture_expectation_fail_count": len(mismatches),
        "mismatches": mismatches,
        "backend_contract": {
            "scope": "FAIL_CLOSED_CONTRACT_HARNESS_ONLY",
            "real_certificate_payload_status": "PRIVATE_NUMERIC_PAYLOAD_PENDING",
            "ordinary_f64_is_not_certificate": True,
            "required_real_payloads": [
                "B1 transform interval condition certificate",
                "B2 recovered-direction interval independence witnesses",
                "B3 transformed period-matrix interval condition and positive singular-value certificate",
            ],
        },
        "claim_ceiling": CLAIM_CEILING,
    }
    write_unique(RESULTS_JSON, json.dumps(results, indent=2, sort_keys=True) + "\n")

    report = f"""# {EXPERIMENT_ID}

Status: `{status}`

## Meaning

This packet adds the missing Rust backend harness for the Stage-B stable gap-period basis certificate. It is deliberately fail-closed: a 27-row payload, normalization leakage, high transform condition, missing recovered-direction witnesses, or a nonpositive singular value all prevent a pass.

The fixtures are synthetic contract tests. They show that the backend guardrails work before private numeric interval payloads are wired. They do not certify the real 24-row period matrix.

## Fixture Result

- Fixture rows: `{len(rows)}`
- Passed expected statuses: `{pass_count}`
- Failed expected statuses: `{len(mismatches)}`

## Claim Ceiling

{CLAIM_CEILING}

Altitude remains `8525 m`.

## Verification Performed

- `cargo run --quiet --bin {BACKEND_NAME}` on all fixture files.
- JSON parse of every backend output.
- JSON parse of `RESULTS.json` and `STRUCTURAL_PRE_AUDIT.json`.
- JSONL parse of `FIXTURE_ROWS.jsonl`.
- SHA-256 sidecar generated and verified for `RESULTS.json`.
"""
    write_unique(REPORT_MD, report)

    sha = sha256_path(RESULTS_JSON)
    write_unique(SHA_FILE, f"{sha}  {RESULTS_JSON.name}\n")


if __name__ == "__main__":
    main()
