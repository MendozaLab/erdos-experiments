#!/usr/bin/env python3
"""Stage-B interval-certificate gate for the stable gap-period basis.

This packet consumes the Stage-A f64 conditioning target and asks whether the
best transformed basis has a real Rust/Inari directed-interval certificate. The
allowed write scope for this packet does not include adding a Rust backend, so
the runner records the missing backend contract instead of pretending that f64
linear algebra is directed interval arithmetic.

Claim ceiling: Stage-B certificate attempt/contract only. This packet does not
provide a Rust/Inari directed interval condition certificate, does not
interval-audit a period residual, does not prove period-matrix legitimacy, does
not compose attainment, does not prove selector existence, does not close KKT
composition, does not give a global reduction, does not solve #1038, and does
not improve public SOTA.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPERIMENT_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-STABLE-GAP-PERIOD-BASIS-INTERVAL-"
    "CERTIFICATE-20260527-01"
)

CONDITIONING_PACKET_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-STABLE-GAP-PERIOD-BASIS-CONDITIONING-GATE-"
    "20260527-01"
)
INTERIOR_PILOT_PACKET_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-INTERIOR-GAP-SOURCE-REGULARIZATION-PILOT-"
    "20260527-02"
)
GAP_VECTOR_PACKET_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-GAP-PERIOD-VECTOR-LEGITIMACY-REWRITE-"
    "20260527-01"
)

HERE = Path(__file__).resolve().parent
PACKET_ROOT = HERE.parent / "erdos-1038"
RUST_BACKEND_EXPECTED = (
    HERE
    / "src"
    / "bin"
    / "phi_k_stable_gap_period_basis_interval_certificate_backend.rs"
)

RESULTS_JSON = PACKET_ROOT / f"{EXPERIMENT_ID}_RESULTS.json"
REPORT_MD = PACKET_ROOT / f"{EXPERIMENT_ID}_REPORT.md"
STRUCTURAL_PRE_AUDIT_JSON = (
    PACKET_ROOT / f"{EXPERIMENT_ID}_STRUCTURAL_PRE_AUDIT.json"
)
BACKEND_REQUIREMENTS_JSON = (
    PACKET_ROOT / f"{EXPERIMENT_ID}_STAGE_B_BACKEND_REQUIREMENTS.json"
)
SUBOBLIGATION_ROWS_JSONL = (
    PACKET_ROOT / f"{EXPERIMENT_ID}_STAGE_B_SUBOBLIGATION_ROWS.jsonl"
)
SHA_FILE = PACKET_ROOT / f"{EXPERIMENT_ID}_RESULTS.sha256"

CLAIM_CEILING = (
    "Stage-B stable gap-period basis interval-certificate attempt only. The "
    "packet verifies the structural pre-audit and records that the needed "
    "Rust/Inari backend is not available inside the assigned write scope. It "
    "does not provide a directed interval condition certificate, does not "
    "interval-audit a period residual, does not prove period-matrix "
    "legitimacy, does not compose attainment, does not prove selector "
    "existence, does not close KKT composition, does not give a global "
    "reduction, does not solve #1038, and does not improve public SOTA."
)

ALLOWED_CERTIFICATION_SCOPES = {
    "F64_SAMPLED_ONLY",
    "F64_INTERVAL_CONSTRAINED",
    "FIXED_PROJECTION_DIRECTED_INTERVAL",
    "FIXED_ROOT_BOX_DIRECTED_INTERVAL",
    "WITNESS_TYPED_REGION_DIRECTED_INTERVAL",
    "INDEPENDENT_COEFFICIENT_BOX_DIRECTED_INTERVAL",
    "LEAN_PROVEN",
}
ALLOWED_CERTIFICATION_ARITHMETIC = {
    "F64",
    "F64_EXTENDED",
    "RUST_INARI_DIRECTED_INTERVAL",
    "LEAN_PROOF",
}


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def packet_path(packet_id: str, suffix: str) -> Path:
    return PACKET_ROOT / f"{packet_id}_{suffix}"


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


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def typed_number(
    value: float | int | None,
    *,
    scope: str,
    arithmetic: str,
    audit_packet: str,
    note: str,
) -> dict[str, Any]:
    if scope not in ALLOWED_CERTIFICATION_SCOPES:
        raise ValueError(f"Unknown certification_scope: {scope}")
    if arithmetic not in ALLOWED_CERTIFICATION_ARITHMETIC:
        raise ValueError(f"Unknown certification_arithmetic: {arithmetic}")
    return {
        "value": value,
        "certification_scope": scope,
        "certification_arithmetic": arithmetic,
        "audit_packet": audit_packet,
        "note": note,
    }


def structural_pre_audit_payload() -> dict[str, Any]:
    return {
        "packet_id": EXPERIMENT_ID,
        "object_being_audited": (
            "Stage-B interval certificate for the transformed 24-row "
            "gap-period matrix basis selected by the stable-basis "
            "conditioning gate"
        ),
        "claimed_dimension": 24,
        "dimension_derivation": (
            "24 gap-period kernel rows from the 25-component PUBLIC cloud: "
            "one normalization row is not a period row, so the period-kernel "
            "dimension is m-1 for m=25 components"
        ),
        "object_count": 24,
        "count_derivation": (
            "The consumed conditioning gate records the best basis matrix "
            "shape as [24, 24], and the gap-period vector legitimacy rewrite "
            "records the 24-row period contract"
        ),
        "sign_convention": (
            "Right-transform convention: transformed_period_matrix = "
            "seed_period_matrix @ basis_transform. Rows are gap-period "
            "functionals; columns are transformed basis functions. The signed "
            "principal-value source convention is inherited but not evaluated "
            "in this packet."
        ),
        "sign_convention_source": (
            f"{CONDITIONING_PACKET_ID} for the right-transform basis convention; "
            f"{INTERIOR_PILOT_PACKET_ID} for signed principal-value source "
            "orientation; this packet performs no potential-sign replay"
        ),
        "expected_failure_modes": [
            "the f64 orthogonalizing transform cannot be interval-certified",
            "the 11 rank-recovered directions are not interval-independent",
            "the transformed period matrix has interval condition number above threshold",
            "a superseded upstream packet is accidentally consumed",
            "f64 condition evidence is mistaken for directed interval evidence",
        ],
    }


def validate_structural_pre_audit(payload: dict[str, Any]) -> dict[str, Any]:
    required = [
        "object_being_audited",
        "claimed_dimension",
        "dimension_derivation",
        "object_count",
        "count_derivation",
        "sign_convention",
        "sign_convention_source",
        "expected_failure_modes",
    ]
    missing = [key for key in required if key not in payload]
    type_errors: list[str] = []
    if payload.get("claimed_dimension") != 24:
        type_errors.append("claimed_dimension must be 24")
    if payload.get("object_count") != 24:
        type_errors.append("object_count must be 24")
    if not isinstance(payload.get("expected_failure_modes"), list):
        type_errors.append("expected_failure_modes must be a list")
    if not str(payload.get("sign_convention", "")).strip():
        type_errors.append("sign_convention must be nonempty")
    return {
        "status": "PASS" if not missing and not type_errors else "FAIL",
        "missing_fields": missing,
        "type_or_value_errors": type_errors,
    }


def emit_and_parse_structural_pre_audit() -> tuple[dict[str, Any], dict[str, Any]]:
    payload = structural_pre_audit_payload()
    write_unique(
        STRUCTURAL_PRE_AUDIT_JSON,
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    parsed = read_json(STRUCTURAL_PRE_AUDIT_JSON)
    validation = validate_structural_pre_audit(parsed)
    if validation["status"] != "PASS":
        raise RuntimeError(f"STRUCTURAL_PRE_AUDIT validation failed: {validation}")
    return parsed, validation


def supersession_index() -> dict[str, Any]:
    rows = []
    for path in sorted(PACKET_ROOT.glob("*_SUPERSEDED_BY.json")):
        try:
            payload = read_json(path)
        except Exception as exc:  # pragma: no cover - defensive audit
            rows.append({"path": str(path), "parse_error": str(exc)})
            continue
        rows.append({"path": str(path), **payload})
    return {
        "sidecar_count": len(rows),
        "rows": rows,
        "superseded_packet_ids": sorted(
            {
                str(row["superseded_packet_id"])
                for row in rows
                if "superseded_packet_id" in row
            }
        ),
    }


def sidecar_status(packet_id: str) -> dict[str, Any]:
    results_path = packet_path(packet_id, "RESULTS.json")
    sha_path = packet_path(packet_id, "RESULTS.sha256")
    status: dict[str, Any] = {
        "packet_id": packet_id,
        "exists": results_path.exists() and sha_path.exists(),
        "status": "MISSING",
        "sha_ok": False,
    }
    if results_path.exists():
        payload = read_json(results_path)
        status["source_status"] = payload.get("status")
        status["source_verdict"] = payload.get("verdict")
    if results_path.exists() and sha_path.exists():
        expected = sha_path.read_text(encoding="utf-8").strip().split()[0]
        actual = sha256_path(results_path)
        status.update(
            {
                "results_path": str(results_path),
                "expected_sha256": expected,
                "actual_sha256": actual,
                "sha_ok": expected == actual,
                "status": "OK" if expected == actual else "MISMATCH",
            }
        )
    return status


def load_sources() -> dict[str, Any]:
    supersession = supersession_index()
    consumed_packet_ids = [
        CONDITIONING_PACKET_ID,
        INTERIOR_PILOT_PACKET_ID,
        GAP_VECTOR_PACKET_ID,
    ]
    stale = [
        packet_id
        for packet_id in consumed_packet_ids
        if packet_id in supersession["superseded_packet_ids"]
    ]
    if stale:
        raise RuntimeError(f"Refusing to consume superseded packets: {stale}")

    sources = {
        "stable_basis_conditioning_gate": sidecar_status(CONDITIONING_PACKET_ID),
        "interior_gap_source_regularization_pilot": sidecar_status(
            INTERIOR_PILOT_PACKET_ID
        ),
        "gap_period_vector_legitimacy_rewrite": sidecar_status(GAP_VECTOR_PACKET_ID),
        "supersession_scan": supersession,
    }
    bad = [
        name
        for name, row in sources.items()
        if name != "supersession_scan" and row.get("status") != "OK"
    ]
    if bad:
        raise RuntimeError(f"Source sidecar check failed for: {bad}")
    return sources


def build_stage_b_backend_requirements(best_basis: dict[str, Any]) -> dict[str, Any]:
    return {
        "packet_id": EXPERIMENT_ID,
        "target_backend": str(RUST_BACKEND_EXPECTED),
        "target_backend_exists": RUST_BACKEND_EXPECTED.exists(),
        "status": "BLOCKED_BACKEND_NOT_IMPLEMENTED"
        if not RUST_BACKEND_EXPECTED.exists()
        else "BACKEND_PRESENT_REQUIRES_EXECUTION",
        "best_f64_target_basis": best_basis["basis_name"],
        "required_stages": [
            {
                "stage": "B1",
                "name": "transform_interval_certificate",
                "required_fields": [
                    "basis_transform_entry_intervals",
                    "transform_condition_interval_upper_bound",
                    "transform_column_action_interval_bounds",
                    "proof_that_transform_matches_right-transform convention",
                ],
                "pass_condition": "interval upper bound on transform condition < 1e10",
            },
            {
                "stage": "B2",
                "name": "rank_recovered_direction_certificate",
                "required_fields": [
                    "robust_seed_rank",
                    "recovered_direction_count",
                    "interval_independence_witnesses_for_recovered_directions",
                    "span-separation lower bounds from the 13 robust directions",
                ],
                "pass_condition": (
                    "all recovered directions have interval-certified "
                    "separation from the robust 13-dimensional seed span"
                ),
            },
            {
                "stage": "B3",
                "name": "transformed_matrix_interval_certificate",
                "required_fields": [
                    "transformed_period_matrix_entry_intervals",
                    "smallest_singular_value_interval_lower_bound",
                    "largest_singular_value_interval_upper_bound",
                    "condition_number_interval_upper_bound",
                ],
                "pass_condition": "interval upper bound on transformed matrix condition < 1e10",
            },
        ],
        "claim_ceiling": (
            "Backend requirements only. These fields are not computed by this "
            "Python packet."
        ),
    }


def build_subobligations(best_basis: dict[str, Any]) -> list[dict[str, Any]]:
    recovered_count = 24 - int(best_basis["seed_eval_rank_f64"])
    return [
        {
            "name": "structural_pre_audit",
            "status": "PASS",
            "evidence": str(STRUCTURAL_PRE_AUDIT_JSON),
        },
        {
            "name": "supersession_scan_for_consumed_packets",
            "status": "PASS",
            "evidence": "No consumed packet id appears in parsed SUPERSEDED_BY sidecars.",
        },
        {
            "name": "stage_a_target_basis_available",
            "status": "PASS_F64_TARGET_ONLY",
            "evidence": f"{CONDITIONING_PACKET_ID} selected {best_basis['basis_name']}",
        },
        {
            "name": "B1_transform_interval_certificate",
            "status": "BLOCKED_BACKEND_NOT_IMPLEMENTED",
            "evidence": "No allowed-scope Rust/Inari backend emits transform interval bounds.",
        },
        {
            "name": "B2_rank_recovered_direction_certificate",
            "status": "BLOCKED_BACKEND_NOT_IMPLEMENTED",
            "evidence": (
                f"Seed rank is {best_basis['seed_eval_rank_f64']}; "
                f"{recovered_count} directions need interval independence witnesses."
            ),
        },
        {
            "name": "B3_transformed_matrix_interval_certificate",
            "status": "BLOCKED_BACKEND_NOT_IMPLEMENTED",
            "evidence": (
                "No allowed-scope Rust/Inari backend emits transformed matrix "
                "entry intervals or singular-value bounds."
            ),
        },
        {
            "name": "period_residual_audit_readiness",
            "status": "NOT_READY",
            "evidence": "Requires B1, B2, and B3 directed interval certificates first.",
        },
    ]


def build_results(
    *,
    sources: dict[str, Any],
    pre_audit_validation: dict[str, Any],
    conditioning_results: dict[str, Any],
    best_audit: dict[str, Any],
    backend_requirements: dict[str, Any],
    sub_obligations: list[dict[str, Any]],
) -> dict[str, Any]:
    best_basis = best_audit["best_basis"]
    recovered_count = 24 - int(best_basis["seed_eval_rank_f64"])
    typed_scope = "F64_SAMPLED_ONLY"
    typed_arithmetic = "F64"
    return {
        "packet_id": EXPERIMENT_ID,
        "version": "20260527-01",
        "created_utc": now_utc(),
        "route_type": "stable_gap_period_basis_interval_certificate",
        "status": (
            "STABLE_GAP_PERIOD_BASIS_INTERVAL_CERTIFICATE_PARTIAL__"
            "STRUCTURAL_PRE_AUDIT_PASS__STAGE_B_BLOCKED_BACKEND_REQUIRED"
        ),
        "verdict": "NO_DIRECTED_INTERVAL_CERTIFICATE_PRODUCED__BACKEND_CONTRACT_EMITTED",
        "claim_ceiling": CLAIM_CEILING,
        "altitude": "8525 m -- unchanged",
        "public_sota": "Problem open; no public SOTA improvement.",
        "source_packets": sources,
        "structural_pre_audit": {
            "path": str(STRUCTURAL_PRE_AUDIT_JSON),
            "validation": pre_audit_validation,
        },
        "structural_fingerprint": {
            "period_row_count": 24,
            "kernel_dimension": 24,
            "normalization_rows_excluded_from_period_kernel": 1,
            "basis_family": "chebyshev_global_weighted_qr_orthogonalized",
            "best_basis_name": best_basis["basis_name"],
            "seed_basis": best_basis["seed_basis"],
            "transform_kind": best_basis["transform_kind"],
            "f64_condition_number": typed_number(
                best_basis["condition_number_f64"],
                scope=typed_scope,
                arithmetic=typed_arithmetic,
                audit_packet=CONDITIONING_PACKET_ID,
                note="Inherited Stage-A target metric; not an interval certificate.",
            ),
            "transform_condition": typed_number(
                best_basis["transform_condition_f64"],
                scope=typed_scope,
                arithmetic=typed_arithmetic,
                audit_packet=CONDITIONING_PACKET_ID,
                note="Inherited f64 transform condition; requires Stage B1 intervalization.",
            ),
            "seed_eval_rank": typed_number(
                best_basis["seed_eval_rank_f64"],
                scope=typed_scope,
                arithmetic=typed_arithmetic,
                audit_packet=CONDITIONING_PACKET_ID,
                note="Inherited f64 seed evaluation rank; implies 11 recovered directions.",
            ),
            "rank_recovered_direction_count": typed_number(
                recovered_count,
                scope=typed_scope,
                arithmetic=typed_arithmetic,
                audit_packet=CONDITIONING_PACKET_ID,
                note="Computed as 24 - f64 seed_eval_rank; not interval-certified.",
            ),
        },
        "typed_numeric_results": {
            "best_f64_condition_number": typed_number(
                best_basis["condition_number_f64"],
                scope=typed_scope,
                arithmetic=typed_arithmetic,
                audit_packet=CONDITIONING_PACKET_ID,
                note="Diagnostic target only.",
            ),
            "best_f64_min_singular_value": typed_number(
                best_basis["min_singular_value_f64"],
                scope=typed_scope,
                arithmetic=typed_arithmetic,
                audit_packet=CONDITIONING_PACKET_ID,
                note="Diagnostic target only; not directed interval lower bound.",
            ),
            "best_f64_max_singular_value": typed_number(
                best_basis["max_singular_value_f64"],
                scope=typed_scope,
                arithmetic=typed_arithmetic,
                audit_packet=CONDITIONING_PACKET_ID,
                note="Diagnostic target only; not directed interval upper bound.",
            ),
            "best_f64_transform_condition": typed_number(
                best_basis["transform_condition_f64"],
                scope=typed_scope,
                arithmetic=typed_arithmetic,
                audit_packet=CONDITIONING_PACKET_ID,
                note="Diagnostic target only; transform interval certificate is blocked.",
            ),
            "condition_threshold_target": typed_number(
                conditioning_results["summary"]["condition_threshold"],
                scope=typed_scope,
                arithmetic=typed_arithmetic,
                audit_packet=CONDITIONING_PACKET_ID,
                note="Stage-A threshold used to choose the target basis.",
            ),
        },
        "sub_obligations": sub_obligations,
        "stage_b": {
            "backend_requirements_path": str(BACKEND_REQUIREMENTS_JSON),
            "backend_requirements_status": backend_requirements["status"],
            "B1_transform_certificate": "BLOCKED_BACKEND_NOT_IMPLEMENTED",
            "B2_rank_recovered_directions": "BLOCKED_BACKEND_NOT_IMPLEMENTED",
            "B3_transformed_matrix_certificate": "BLOCKED_BACKEND_NOT_IMPLEMENTED",
            "directed_interval_certificate_produced": False,
            "period_residual_audit_ready": False,
        },
        "gate_checks": {
            "structural_pre_audit_exists_and_parseable": True,
            "structural_pre_audit_pass": pre_audit_validation["status"] == "PASS",
            "source_result_sidecars_ok": all(
                row.get("status") == "OK"
                for key, row in sources.items()
                if key != "supersession_scan"
            ),
            "superseded_source_consumed": False,
            "stage_a_f64_target_basis_available": True,
            "stage_b_backend_available_in_write_scope": backend_requirements[
                "target_backend_exists"
            ],
            "stage_b_directed_interval_certificate_performed": False,
            "stage_b_directed_interval_certificate_pass": False,
            "numerical_pass_verdict_allowed": False,
            "altitude_reconsideration_allowed": False,
        },
        "products": {
            "structural_pre_audit": str(STRUCTURAL_PRE_AUDIT_JSON),
            "backend_requirements": str(BACKEND_REQUIREMENTS_JSON),
            "subobligation_rows": str(SUBOBLIGATION_ROWS_JSONL),
            "report": str(REPORT_MD),
        },
        "recommended_next_pitch": (
            "EXP-MATH-ERDOS1038-PHI-K-STABLE-GAP-PERIOD-BASIS-RUST-INARI-"
            "BACKEND-20260527-01"
        ),
        "next_blocker": (
            "Implement the Rust/Inari backend that emits B1 transform intervals, "
            "B2 recovered-direction independence witnesses, and B3 transformed "
            "matrix singular-value/condition bounds."
        ),
    }


def build_report(results: dict[str, Any]) -> str:
    fp = results["structural_fingerprint"]
    stage_b = results["stage_b"]
    return f"""# Stable Gap-Period Basis Interval Certificate

Packet: `{EXPERIMENT_ID}`

## Verdict

`{results["verdict"]}`

## Meaning

The structural pre-audit passed, but this packet did not produce the Stage-B
directed interval certificate. The best basis from Stage A is still only an f64
target until a Rust/Inari backend certifies the transform, recovered directions,
and transformed matrix.

```text
period_row_count = {fp["period_row_count"]}
kernel_dimension = {fp["kernel_dimension"]}
basis_family = {fp["basis_family"]}
best_basis_name = {fp["best_basis_name"]}
best_f64_condition_number = {fp["f64_condition_number"]["value"]}
transform_condition_f64 = {fp["transform_condition"]["value"]}
seed_eval_rank_f64 = {fp["seed_eval_rank"]["value"]}
rank_recovered_direction_count = {fp["rank_recovered_direction_count"]["value"]}
```

## Stage B Status

```text
B1_transform_certificate = {stage_b["B1_transform_certificate"]}
B2_rank_recovered_directions = {stage_b["B2_rank_recovered_directions"]}
B3_transformed_matrix_certificate = {stage_b["B3_transformed_matrix_certificate"]}
directed_interval_certificate_produced = {stage_b["directed_interval_certificate_produced"]}
period_residual_audit_ready = {stage_b["period_residual_audit_ready"]}
```

## Structural Pre-Audit

The required pre-audit was written and parsed before source numeric fields were
consumed:

```text
{STRUCTURAL_PRE_AUDIT_JSON}
```

## Required Backend Fields

The next backend must emit:

- interval bounds for the right-transform entries and column actions;
- an interval upper bound for transform condition;
- interval independence witnesses for the 11 f64 rank-recovered directions;
- interval bounds for transformed matrix entries;
- interval lower/upper singular-value bounds or an accepted determinant/inverse-norm substitute;
- an interval upper bound for transformed matrix condition.

## Products

```text
{STRUCTURAL_PRE_AUDIT_JSON}
{BACKEND_REQUIREMENTS_JSON}
{SUBOBLIGATION_ROWS_JSONL}
```

## Claim Ceiling

{CLAIM_CEILING}

Altitude remains 8525 m. Public SOTA is unchanged and #1038 remains open.
"""


def main() -> None:
    _pre_audit, pre_audit_validation = emit_and_parse_structural_pre_audit()
    sources = load_sources()
    conditioning_results = read_json(packet_path(CONDITIONING_PACKET_ID, "RESULTS.json"))
    best_audit = read_json(packet_path(CONDITIONING_PACKET_ID, "BEST_BASIS_AUDIT.json"))
    best_basis = best_audit["best_basis"]
    backend_requirements = build_stage_b_backend_requirements(best_basis)
    sub_obligations = build_subobligations(best_basis)
    results = build_results(
        sources=sources,
        pre_audit_validation=pre_audit_validation,
        conditioning_results=conditioning_results,
        best_audit=best_audit,
        backend_requirements=backend_requirements,
        sub_obligations=sub_obligations,
    )

    write_unique(
        BACKEND_REQUIREMENTS_JSON,
        json.dumps(backend_requirements, indent=2, sort_keys=True) + "\n",
    )
    write_unique(
        SUBOBLIGATION_ROWS_JSONL,
        "\n".join(json.dumps(row, sort_keys=True) for row in sub_obligations)
        + "\n",
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
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
