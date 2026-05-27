#!/usr/bin/env python3
"""Run the dependent Vieta image consumer scaffold packet.

This runner is intentionally scaffold-only. It verifies fail-closed behavior
over synthetic fixtures and emits immutable packet artifacts. It does not
consume the private numeric receipts needed for a real dependent Vieta theorem.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPERIMENT_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-DEPENDENT-VIETA-IMAGE-CONSUMER-20260527-01"
)
CLAIM_CEILING = (
    "SCAFFOLD_ONLY__NO_DEPENDENT_VIETA_THEOREM_PASS__"
    "NO_1038_SOTA_ALTITUDE_KKT_GLOBAL_LEAN_CLAIMS"
)

ROOT = Path(__file__).resolve().parents[1]
BACKEND_SOURCE = ROOT / "backend-source"
FIXTURE_DIR = (
    BACKEND_SOURCE / "fixtures" / "phi_k_dependent_vieta_image_consumer"
)
PRIVATE_ROUTE_ARTIFACTS = ROOT / "private-route-artifacts"
RUST_BACKEND_EXPECTED = (
    BACKEND_SOURCE / "src" / "bin" / "phi_k_dependent_vieta_image_consumer.rs"
)

RESULTS_JSON = PRIVATE_ROUTE_ARTIFACTS / f"{EXPERIMENT_ID}_RESULTS.json"
REPORT_MD = PRIVATE_ROUTE_ARTIFACTS / f"{EXPERIMENT_ID}_REPORT.md"
FIXTURE_ROWS_JSONL = PRIVATE_ROUTE_ARTIFACTS / f"{EXPERIMENT_ID}_FIXTURE_ROWS.jsonl"
STRUCTURAL_PRE_AUDIT_JSON = (
    PRIVATE_ROUTE_ARTIFACTS / f"{EXPERIMENT_ID}_STRUCTURAL_PRE_AUDIT.json"
)
SHA_FILE = PRIVATE_ROUTE_ARTIFACTS / f"{EXPERIMENT_ID}_RESULTS.sha256"

ROOT_BOX_DIAMETER_HI_THRESHOLD = 1.0e-1
SCALED_VIETA_IMAGE_DIAMETER_HI_THRESHOLD = 1.0e-2
FIXED_CLOUD_DEVIATION_HI_THRESHOLD = 1.0e-2


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_unique(path: Path, text: str) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite immutable artifact: {path}")
    path.write_text(text, encoding="utf-8")


def parse_input(text: str) -> tuple[dict[str, str], list[str]]:
    fields: dict[str, str] = {}
    errors: list[str] = []
    for line_index, raw_line in enumerate(text.splitlines()):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            errors.append(f"line {line_index + 1} is not key=value")
            continue
        key, value = line.split("=", 1)
        fields[key.strip()] = value.strip()
    return fields, errors


def _as_usize(fields: dict[str, str], key: str, errors: list[str]) -> int | None:
    raw = fields.get(key)
    if raw is None:
        errors.append(f"missing {key}")
        return None
    try:
        value = int(raw)
        if value < 0:
            raise ValueError
        return value
    except ValueError:
        errors.append(f"{key} is not a usize")
        return None


def _as_f64(fields: dict[str, str], key: str, errors: list[str]) -> float | None:
    raw = fields.get(key)
    if raw is None:
        errors.append(f"missing {key}")
        return None
    try:
        value = float(raw)
    except ValueError:
        errors.append(f"{key} is not an f64")
        return None
    if not math.isfinite(value):
        errors.append(f"{key} is not finite")
        return None
    return value


def _as_bool(fields: dict[str, str], key: str, errors: list[str]) -> bool | None:
    raw = fields.get(key)
    if raw is None:
        errors.append(f"missing {key}")
        return None
    if raw in ("true", "1"):
        return True
    if raw in ("false", "0"):
        return False
    errors.append(f"{key} is not a bool")
    return None


def evaluate_fixture(text: str) -> dict[str, Any]:
    """Reproduce the Rust backend fail-closed verdict in pure Python."""
    fields, errors = parse_input(text)

    claimed_root_count = _as_usize(fields, "claimed_root_count", errors)
    multiplicity_sum = _as_usize(fields, "multiplicity_sum", errors)
    scaling_convention_match = _as_bool(fields, "scaling_convention_match", errors)
    dependency_relation_recorded = _as_bool(
        fields, "dependency_relation_recorded", errors
    )

    root_box_diameter_hi = _as_f64(fields, "root_box_diameter_hi", errors)
    scaled_vieta_image_diameter_hi = _as_f64(
        fields, "scaled_vieta_image_diameter_hi", errors
    )
    vieta_image_dimension = _as_usize(fields, "vieta_image_dimension", errors)
    vieta_image_independent_count = _as_usize(
        fields, "vieta_image_independent_count", errors
    )
    fixed_cloud_anchor_count = _as_usize(fields, "fixed_cloud_anchor_count", errors)
    fixed_cloud_deviation_hi = _as_f64(fields, "fixed_cloud_deviation_hi", errors)

    structure_pass = (
        not errors
        and claimed_root_count is not None
        and multiplicity_sum is not None
        and claimed_root_count == multiplicity_sum
        and claimed_root_count > 0
        and scaling_convention_match is True
        and dependency_relation_recorded is True
    )

    if not structure_pass:
        v1_status = "BLOCKED_PRE_AUDIT"
    elif (
        root_box_diameter_hi is not None
        and 0.0 < root_box_diameter_hi < ROOT_BOX_DIAMETER_HI_THRESHOLD
    ):
        v1_status = "PASS_CONTRACT_FIXTURE"
    else:
        v1_status = "FAIL_CONTRACT_FIXTURE"

    if not structure_pass:
        v2_status = "BLOCKED_PRE_AUDIT"
    elif (
        scaled_vieta_image_diameter_hi is not None
        and 0.0
        < scaled_vieta_image_diameter_hi
        < SCALED_VIETA_IMAGE_DIAMETER_HI_THRESHOLD
        and vieta_image_dimension is not None
        and vieta_image_independent_count == vieta_image_dimension
    ):
        v2_status = "PASS_CONTRACT_FIXTURE"
    else:
        v2_status = "FAIL_CONTRACT_FIXTURE"

    v1_pass = v1_status == "PASS_CONTRACT_FIXTURE"
    v2_pass = v2_status == "PASS_CONTRACT_FIXTURE"
    if not structure_pass or not v1_pass or not v2_pass:
        v3_status = "BLOCKED_PRECONDITION"
    elif (
        fixed_cloud_anchor_count is not None
        and fixed_cloud_anchor_count > 0
        and fixed_cloud_deviation_hi is not None
        and 0.0 <= fixed_cloud_deviation_hi < FIXED_CLOUD_DEVIATION_HI_THRESHOLD
    ):
        v3_status = "PASS_CONTRACT_FIXTURE"
    else:
        v3_status = "FAIL_CONTRACT_FIXTURE"

    if not structure_pass:
        overall_status = "STRUCTURAL_PRE_AUDIT_FAIL_CLOSED"
    elif (
        v1_status.startswith("FAIL")
        or v2_status.startswith("FAIL")
        or v3_status.startswith("FAIL")
    ):
        overall_status = "FAIL_CLOSED_CONTRACT_FIXTURE_REJECTED"
    elif v1_pass and v2_pass and v3_status == "PASS_CONTRACT_FIXTURE":
        overall_status = "PASS_CONTRACT_FIXTURE_ONLY_PRIVATE_NUMERIC_PAYLOAD_REQUIRED"
    else:
        overall_status = "BLOCKED_CONTRACT_FIXTURE"

    return {
        "parse_errors": errors,
        "structural_pre_audit": {
            "claimed_root_count": claimed_root_count,
            "multiplicity_sum": multiplicity_sum,
            "scaling_convention_match": scaling_convention_match,
            "dependency_relation_recorded": dependency_relation_recorded,
            "status": "PASS" if structure_pass else "FAIL",
        },
        "stage_statuses": {
            "V1_root_box_enclosure": v1_status,
            "V2_scaled_vieta_image_enclosure": v2_status,
            "V3_fixed_cloud_certificate": v3_status,
        },
        "overall_status": overall_status,
    }


def structural_pre_audit_payload() -> dict[str, Any]:
    return {
        "packet_id": EXPERIMENT_ID,
        "object_being_audited": (
            "Scaffold contract harness for a dependent Vieta image consumer: "
            "validates synthetic root-box, multiplicity, scaled Vieta image, "
            "and fixed-cloud fields before any future numeric certificate is "
            "allowed to claim a dependent-image pass."
        ),
        "claimed_root_count": 24,
        "claim_derivation": (
            "Synthetic scaffold fixture dimension only; real n=10000 route "
            "receipts are not consumed in this packet."
        ),
        "multiplicity_sum_expected": 24,
        "sign_convention": (
            "Vieta convention: elementary symmetric functions of roots are read "
            "with standard sign alternation e_k = (-1)^k * coefficient_{n-k} "
            "/ leading. Dependent Vieta image means coefficients must carry a "
            "root-box witness, not independent coefficient-box membership."
        ),
        "expected_failure_modes": [
            "claimed_root_count != multiplicity_sum",
            "scaling_convention_match is not recorded",
            "dependency_relation_recorded is not recorded",
            "root-box diameter above threshold",
            "scaled Vieta image diameter above threshold",
            "scaled Vieta image independence witness count below dimension",
            "fixed-cloud anchor count is zero",
            "fixed-cloud deviation above threshold",
        ],
        "scaffold_only": True,
        "no_dependent_vieta_theorem_pass_claimed": True,
    }


def validate_structural_pre_audit(payload: dict[str, Any]) -> dict[str, Any]:
    required = [
        "object_being_audited",
        "claimed_root_count",
        "claim_derivation",
        "multiplicity_sum_expected",
        "sign_convention",
        "expected_failure_modes",
        "scaffold_only",
        "no_dependent_vieta_theorem_pass_claimed",
    ]
    missing = [key for key in required if key not in payload]
    type_errors: list[str] = []
    if payload.get("claimed_root_count") != 24:
        type_errors.append("claimed_root_count must be 24")
    if payload.get("multiplicity_sum_expected") != 24:
        type_errors.append("multiplicity_sum_expected must be 24")
    if not isinstance(payload.get("expected_failure_modes"), list):
        type_errors.append("expected_failure_modes must be a list")
    if payload.get("scaffold_only") is not True:
        type_errors.append("scaffold_only must be True")
    if payload.get("no_dependent_vieta_theorem_pass_claimed") is not True:
        type_errors.append("no_dependent_vieta_theorem_pass_claimed must be True")
    if not str(payload.get("sign_convention", "")).strip():
        type_errors.append("sign_convention must be nonempty")
    return {
        "status": "PASS" if not missing and not type_errors else "FAIL",
        "missing_fields": missing,
        "type_or_value_errors": type_errors,
    }


def _expected_class(fixture_name: str) -> str:
    if fixture_name.startswith("good_"):
        return "PASS_CONTRACT_FIXTURE_ONLY_PRIVATE_NUMERIC_PAYLOAD_REQUIRED"
    if fixture_name.startswith("bad_"):
        return "STRUCTURAL_PRE_AUDIT_FAIL_CLOSED"
    if fixture_name.startswith(("v1_", "v2_", "v3_")):
        return "FAIL_CLOSED_CONTRACT_FIXTURE_REJECTED"
    return "UNCLASSIFIED"


def _row_class_match(row: dict[str, Any]) -> bool:
    return row["expected_overall_status_class"] == row["computed_overall_status"]


def run_all_fixtures() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for fixture_path in sorted(FIXTURE_DIR.glob("*.txt")):
        text = fixture_path.read_text(encoding="utf-8")
        verdict = evaluate_fixture(text)
        rows.append(
            {
                "fixture": fixture_path.name,
                "fixture_path": str(fixture_path),
                "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
                "expected_overall_status_class": _expected_class(fixture_path.name),
                "computed_overall_status": verdict["overall_status"],
                "computed_stage_statuses": verdict["stage_statuses"],
                "computed_structural_pre_audit": verdict["structural_pre_audit"],
                "parse_errors": verdict["parse_errors"],
                "match": _expected_class(fixture_path.name)
                == verdict["overall_status"],
            }
        )
    return rows


def build_results(
    *,
    pre_audit_validation: dict[str, Any],
    fixture_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    all_match = all(_row_class_match(row) for row in fixture_rows)
    return {
        "packet_id": EXPERIMENT_ID,
        "version": "20260527-01",
        "created_utc": now_utc(),
        "route_type": "dependent_vieta_image_consumer_scaffold",
        "status": (
            "DEPENDENT_VIETA_IMAGE_CONSUMER_SCAFFOLD__"
            "SYNTHETIC_FIXTURE_CONTRACT_PASS"
            if all_match
            else "DEPENDENT_VIETA_IMAGE_CONSUMER_SCAFFOLD__"
            "SYNTHETIC_FIXTURE_CONTRACT_FAIL"
        ),
        "verdict": (
            "SYNTHETIC_FAIL_CLOSED_CONTRACT_VERIFIED"
            if all_match
            else "SYNTHETIC_FAIL_CLOSED_CONTRACT_VERIFICATION_FAILED"
        ),
        "claim_ceiling": CLAIM_CEILING,
        "altitude": "8525 m -- unchanged",
        "public_sota": "Problem open; no public SOTA improvement.",
        "scaffold_only": True,
        "no_dependent_vieta_theorem_pass_claimed": True,
        "structural_pre_audit": {
            "path": str(STRUCTURAL_PRE_AUDIT_JSON),
            "validation": pre_audit_validation,
        },
        "rust_backend": {
            "expected_path": str(RUST_BACKEND_EXPECTED),
            "exists": RUST_BACKEND_EXPECTED.exists(),
            "scope": "FAIL_CLOSED_CONTRACT_HARNESS_ONLY",
            "compiled_and_executed_by_runner": False,
            "compilation_status_note": (
                "The Python runner verifies synthetic contract behavior; cargo "
                "check is run separately in local verification."
            ),
        },
        "fixtures": {
            "directory": str(FIXTURE_DIR),
            "rows_path": str(FIXTURE_ROWS_JSONL),
            "all_expected_classes_match": all_match,
            "row_count": len(fixture_rows),
        },
        "thresholds": {
            "root_box_diameter_hi": ROOT_BOX_DIAMETER_HI_THRESHOLD,
            "scaled_vieta_image_diameter_hi": SCALED_VIETA_IMAGE_DIAMETER_HI_THRESHOLD,
            "fixed_cloud_deviation_hi": FIXED_CLOUD_DEVIATION_HI_THRESHOLD,
        },
        "sub_obligations": [
            {
                "name": "structural_pre_audit",
                "status": "PASS"
                if pre_audit_validation["status"] == "PASS"
                else "FAIL",
                "evidence": str(STRUCTURAL_PRE_AUDIT_JSON),
            },
            {
                "name": "synthetic_fixture_contract_match",
                "status": "PASS" if all_match else "FAIL",
                "evidence": str(FIXTURE_ROWS_JSONL),
            },
            {
                "name": "real_private_numeric_payload_consumed",
                "status": "NOT_PERFORMED",
                "evidence": (
                    "No real root-box/Vieta/cloud/attainment numeric payload "
                    "was consumed; this packet is a scaffold only."
                ),
            },
            {
                "name": "dependent_vieta_theorem_pass",
                "status": "NOT_PERFORMED",
                "evidence": (
                    "No theorem pass is attempted; the scaffold only verifies "
                    "synthetic fail-closed contract behavior."
                ),
            },
        ],
        "gate_checks": {
            "structural_pre_audit_exists_and_parseable": True,
            "structural_pre_audit_pass": pre_audit_validation["status"] == "PASS",
            "synthetic_contract_verified": all_match,
            "real_certificate_produced": False,
            "directed_interval_certificate_produced": False,
            "dependent_vieta_theorem_pass": False,
            "numerical_pass_verdict_allowed": False,
            "altitude_reconsideration_allowed": False,
        },
        "products": {
            "structural_pre_audit": str(STRUCTURAL_PRE_AUDIT_JSON),
            "fixture_rows": str(FIXTURE_ROWS_JSONL),
            "report": str(REPORT_MD),
        },
        "recommended_next_pitch": (
            "Wire real receipts ROOT_BOX.json, ROOT_MULTIPLICITY_LEDGER.json, "
            "ORDERED_ROOT_INTERVALS.json, SCALED_VIETA_IMAGE_CONTRACT.json, "
            "FIXED_CLOUD_BOUND_CERTIFICATE.json, and "
            "ATTAINED_WITNESS_TYPED_DUAL_MARGIN_RESULTS.json into this consumer."
        ),
        "next_blocker": (
            "No allowed-scope source of real dependent Vieta image numeric "
            "payload is wired yet."
        ),
    }


def build_report(results: dict[str, Any]) -> str:
    return f"""# Dependent Vieta Image Consumer Scaffold

Packet: `{EXPERIMENT_ID}`

## Verdict

`{results["verdict"]}`

## Meaning

This packet adds a fail-closed contract harness for a future dependent Vieta
image consumer. The Rust backend binary `{RUST_BACKEND_EXPECTED.name}` parses a
key=value input and emits a JSON verdict over a structural pre-audit and three
stage statuses:

- `V1_root_box_enclosure`
- `V2_scaled_vieta_image_enclosure`
- `V3_fixed_cloud_certificate`

The Python runner reproduces the same fail-closed logic in pure Python over
synthetic fixtures, hashes each fixture, and records that all expected
fail-closed classes match.

## Synthetic Fixture Status

```text
all_expected_classes_match = {results["fixtures"]["all_expected_classes_match"]}
row_count = {results["fixtures"]["row_count"]}
```

## Structural Pre-Audit

```text
{STRUCTURAL_PRE_AUDIT_JSON}
```

## Rust Backend Status

```text
expected_path = {results["rust_backend"]["expected_path"]}
exists = {results["rust_backend"]["exists"]}
compiled_and_executed_by_runner = {results["rust_backend"]["compiled_and_executed_by_runner"]}
```

`{results["rust_backend"]["compilation_status_note"]}`

## Claim Ceiling

{CLAIM_CEILING}

Altitude remains 8525 m. Public SOTA is unchanged and #1038 remains open. No
dependent Vieta theorem pass, no #1038/SOTA/altitude/KKT/global/Lean claim is
implied by anything in this packet.
"""


def main() -> None:
    PRIVATE_ROUTE_ARTIFACTS.mkdir(parents=True, exist_ok=True)

    pre_audit_payload = structural_pre_audit_payload()
    write_unique(
        STRUCTURAL_PRE_AUDIT_JSON,
        json.dumps(pre_audit_payload, indent=2, sort_keys=True) + "\n",
    )
    pre_audit_parsed = json.loads(STRUCTURAL_PRE_AUDIT_JSON.read_text("utf-8"))
    pre_audit_validation = validate_structural_pre_audit(pre_audit_parsed)
    if pre_audit_validation["status"] != "PASS":
        raise RuntimeError(
            f"STRUCTURAL_PRE_AUDIT validation failed: {pre_audit_validation}"
        )

    fixture_rows = run_all_fixtures()
    write_unique(
        FIXTURE_ROWS_JSONL,
        "\n".join(json.dumps(row, sort_keys=True) for row in fixture_rows) + "\n",
    )

    results = build_results(
        pre_audit_validation=pre_audit_validation,
        fixture_rows=fixture_rows,
    )

    write_unique(REPORT_MD, build_report(results))
    write_unique(RESULTS_JSON, json.dumps(results, indent=2, sort_keys=True) + "\n")
    write_unique(SHA_FILE, f"{sha256_path(RESULTS_JSON)}  {RESULTS_JSON.name}\n")

    print(
        json.dumps(
            {
                "packet_id": EXPERIMENT_ID,
                "status": results["status"],
                "verdict": results["verdict"],
                "all_expected_classes_match": results["fixtures"][
                    "all_expected_classes_match"
                ],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
