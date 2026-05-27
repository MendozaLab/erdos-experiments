#!/usr/bin/env python3
"""Stable basis conditioning gate for the 24-row gap-period matrix.

The interior-gap source pilot found that signed-principal-value RHS rows are
finite in f64, but the 25-component gap-period matrix is unusable in naive
coordinates. This gate tests whether a geometry-aware/orthogonalized basis
brings the matrix into a condition range that could justify a later directed
interval certificate.

Claim ceiling: f64 basis-conditioning diagnostic only. This does not certify
singular values in Rust/Inari, does not interval-audit a period residual, does
not prove signed period-matrix legitimacy, does not compose attainment, does
not close KKT composition, does not give a global reduction, does not solve
#1038, and does not improve public SOTA.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

import run_phi_k_interior_gap_source_regularization_pilot as interior_pilot


EXPERIMENT_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-STABLE-GAP-PERIOD-BASIS-CONDITIONING-GATE-"
    "20260527-01"
)

INTERIOR_PILOT_PACKET_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-INTERIOR-GAP-SOURCE-REGULARIZATION-PILOT-"
    "20260527-02"
)
BOUNDARY_CLASSIFICATION_PACKET_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-BOUNDARY-NEAR-AMBIGUOUS-ATOM-CLASSIFICATION-"
    "20260527-01"
)
GAP_VECTOR_PACKET_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-GAP-PERIOD-VECTOR-LEGITIMACY-REWRITE-"
    "20260527-01"
)
SOURCE_KERNEL_SCAFFOLD_PACKET_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-ENDPOINT-AND-GAP-SOURCE-KERNEL-THEOREM-"
    "SCAFFOLD-20260527-01"
)
FULL_PANEL_PACKET_ID = (
    "EXP-MATH-ERDOS1038-PHI-K-SCALED-COEFFICIENT-ROOT-FACTOR-FULL-PANEL-"
    "REPLAY-20260527-01"
)

HERE = Path(__file__).resolve().parent
PACKET_ROOT = HERE.parent / "erdos-1038"

RESULTS_JSON = PACKET_ROOT / f"{EXPERIMENT_ID}_RESULTS.json"
REPORT_MD = PACKET_ROOT / f"{EXPERIMENT_ID}_REPORT.md"
BASIS_SWEEP_ROWS = PACKET_ROOT / f"{EXPERIMENT_ID}_BASIS_SWEEP_ROWS.jsonl"
BEST_BASIS_AUDIT_JSON = PACKET_ROOT / f"{EXPERIMENT_ID}_BEST_BASIS_AUDIT.json"
STAGE_B_CERTIFICATE_JSON = PACKET_ROOT / f"{EXPERIMENT_ID}_STAGE_B_INTERVAL_CERTIFICATE_STATUS.json"
SHA_FILE = PACKET_ROOT / f"{EXPERIMENT_ID}_RESULTS.sha256"

QUADRATURE_NODES = [96, 128, 192]
PRIMARY_QUADRATURE_NODES = 128
CONDITION_THRESHOLD = 1.0e10
TRANSFORM_CONDITION_WARNING = 1.0e12

CLAIM_CEILING = (
    "Stable gap-period basis conditioning f64 diagnostic only. This packet "
    "does not provide a Rust/Inari directed interval condition certificate, "
    "does not interval-audit a period residual, does not prove signed "
    "period-matrix legitimacy, does not compose attainment, does not prove "
    "selector existence, does not close KKT composition, does not give a global "
    "reduction, does not solve #1038, and does not improve public SOTA."
)


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


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


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
    sources = {
        "interior_gap_source_regularization_pilot": sidecar_status(
            INTERIOR_PILOT_PACKET_ID
        ),
        "boundary_near_ambiguous_atom_classification": sidecar_status(
            BOUNDARY_CLASSIFICATION_PACKET_ID
        ),
        "gap_period_vector_legitimacy_rewrite": sidecar_status(GAP_VECTOR_PACKET_ID),
        "endpoint_and_gap_source_kernel_theorem_scaffold": sidecar_status(
            SOURCE_KERNEL_SCAFFOLD_PACKET_ID
        ),
        "root_factor_full_panel_replay": sidecar_status(FULL_PANEL_PACKET_ID),
    }
    bad = [name for name, row in sources.items() if row.get("status") != "OK"]
    if bad:
        raise RuntimeError(f"Source sidecar check failed for: {bad}")
    return sources


def chebyshev_values(x: float, degree: int, lo: float, hi: float) -> np.ndarray:
    z = (2.0 * x - (lo + hi)) / (hi - lo)
    z = max(-1.0, min(1.0, z))
    theta = math.acos(z)
    return np.array([math.cos(k * theta) for k in range(degree)], dtype=float)


def monomial_values(x: float, degree: int) -> np.ndarray:
    return np.array([x**k for k in range(degree)], dtype=float)


def build_quadrature_grid(
    endpoints: list[float],
    gaps: list[tuple[float, float]],
    node_count: int,
) -> list[dict[str, float]]:
    nodes, weights = np.polynomial.legendre.leggauss(node_count)
    rows: list[dict[str, float]] = []
    for gap_index, (lo, hi) in enumerate(gaps):
        midpoint = 0.5 * (lo + hi)
        half_width = 0.5 * (hi - lo)
        for node, weight in zip(nodes, weights):
            theta = 0.5 * math.pi * (float(node) + 1.0)
            x = midpoint + half_width * math.cos(theta)
            period_weight = (
                float(weight)
                * 0.5
                * math.pi
                * math.exp(
                    -0.5
                    * interior_pilot.log_other_abs_product(x, gap_index, endpoints)
                )
            )
            rows.append(
                {
                    "gap_index": float(gap_index),
                    "x": x,
                    "period_weight": period_weight,
                }
            )
    return rows


def seed_values(
    seed_basis: str,
    x: float,
    degree: int,
    endpoints: list[float],
    gaps: list[tuple[float, float]],
) -> np.ndarray:
    if seed_basis == "monomial":
        return monomial_values(x, degree)
    if seed_basis == "chebyshev_global":
        return chebyshev_values(x, degree, min(endpoints), max(endpoints))
    if seed_basis == "chebyshev_gap_hull":
        return chebyshev_values(
            x,
            degree,
            min(lo for lo, _hi in gaps),
            max(hi for _lo, hi in gaps),
        )
    raise ValueError(f"unknown seed basis: {seed_basis}")


def build_seed_period_matrix_and_eval(
    seed_basis: str,
    endpoints: list[float],
    gaps: list[tuple[float, float]],
    node_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    degree = len(gaps)
    grid = build_quadrature_grid(endpoints, gaps, node_count)
    matrix = np.zeros((degree, degree), dtype=float)
    weighted_eval_rows: list[np.ndarray] = []
    for row in grid:
        gap_index = int(row["gap_index"])
        values = seed_values(seed_basis, row["x"], degree, endpoints, gaps)
        matrix[gap_index, :] += row["period_weight"] * values
        weighted_eval_rows.append(math.sqrt(abs(row["period_weight"])) * values)
    return matrix, np.vstack(weighted_eval_rows)


def matrix_metrics(matrix: np.ndarray) -> dict[str, Any]:
    finite = bool(np.all(np.isfinite(matrix)))
    if not finite:
        return {
            "all_entries_finite": False,
            "condition_number_f64": math.inf,
            "rank_f64": 0,
            "min_singular_value_f64": 0.0,
            "max_singular_value_f64": math.inf,
            "max_abs_entry": math.inf,
            "min_nonzero_abs_entry": None,
        }
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    nonzero = np.abs(matrix[np.nonzero(matrix)])
    return {
        "all_entries_finite": True,
        "condition_number_f64": float(np.linalg.cond(matrix)),
        "rank_f64": int(np.linalg.matrix_rank(matrix)),
        "min_singular_value_f64": float(singular_values[-1]),
        "max_singular_value_f64": float(singular_values[0]),
        "max_abs_entry": float(np.max(np.abs(matrix))),
        "min_nonzero_abs_entry": None
        if len(nonzero) == 0
        else float(np.min(nonzero)),
    }


def row_column_scaled_matrix(matrix: np.ndarray) -> np.ndarray:
    row_norms = np.linalg.norm(matrix, axis=1)
    row_scaled = matrix / (row_norms[:, None] + 1.0e-300)
    col_norms = np.linalg.norm(row_scaled, axis=0)
    return row_scaled / (col_norms[None, :] + 1.0e-300)


def append_candidate(
    rows: list[dict[str, Any]],
    *,
    basis_name: str,
    seed_basis: str | None,
    transform_kind: str,
    node_count: int,
    matrix: np.ndarray,
    transform_norm_f64: float | None = None,
    transform_condition_f64: float | None = None,
    seed_eval_rank_f64: int | None = None,
    seed_eval_condition_f64: float | None = None,
    notes: str = "",
    theorem_basis_candidate: bool = True,
) -> None:
    metrics = matrix_metrics(matrix)
    condition = metrics["condition_number_f64"]
    transform_warning = (
        transform_condition_f64 is not None
        and transform_condition_f64 > TRANSFORM_CONDITION_WARNING
    )
    rank_warning = seed_eval_rank_f64 is not None and seed_eval_rank_f64 < len(matrix)
    rows.append(
        {
            "basis_name": basis_name,
            "seed_basis": seed_basis,
            "transform_kind": transform_kind,
            "quadrature_nodes_per_gap": node_count,
            "shape": list(matrix.shape),
            **metrics,
            "condition_below_threshold": bool(
                math.isfinite(condition) and condition < CONDITION_THRESHOLD
            ),
            "transform_norm_f64": transform_norm_f64,
            "transform_condition_f64": transform_condition_f64,
            "seed_eval_rank_f64": seed_eval_rank_f64,
            "seed_eval_condition_f64": seed_eval_condition_f64,
            "transform_intervalization_warning": bool(transform_warning),
            "seed_rank_warning": bool(rank_warning),
            "theorem_basis_candidate": theorem_basis_candidate,
            "stage_a_status": (
                "PASS_F64_CONDITION_WITH_TRANSFORM_WARNINGS"
                if math.isfinite(condition)
                and condition < CONDITION_THRESHOLD
                and (transform_warning or rank_warning)
                else "PASS_F64_CONDITION"
                if math.isfinite(condition) and condition < CONDITION_THRESHOLD
                else "FAIL_CONDITION_OR_NONFINITE"
            ),
            "notes": notes,
            "claim_ceiling": (
                "F64 basis conditioning row only. This is not a directed "
                "interval certificate."
            ),
        }
    )


def barycentric_weights(nodes: np.ndarray) -> np.ndarray:
    weights = np.ones(len(nodes), dtype=float)
    for j in range(len(nodes)):
        product = 1.0
        for k in range(len(nodes)):
            if j != k:
                product *= nodes[j] - nodes[k]
        weights[j] = 1.0 / product
    return weights


def barycentric_values(x: float, nodes: np.ndarray, weights: np.ndarray) -> np.ndarray:
    diff = x - nodes
    hit = np.where(np.abs(diff) < 1.0e-14)[0]
    if len(hit):
        values = np.zeros(len(nodes), dtype=float)
        values[int(hit[0])] = 1.0
        return values
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        raw = weights / diff
        total = np.sum(raw)
        values = raw / total
    return values


def build_gap_midpoint_lagrange_matrix(
    endpoints: list[float],
    gaps: list[tuple[float, float]],
    node_count: int,
) -> tuple[np.ndarray, str]:
    degree = len(gaps)
    collocation_nodes = np.array([0.5 * (lo + hi) for lo, hi in gaps], dtype=float)
    try:
        bary_weights = barycentric_weights(collocation_nodes)
        grid = build_quadrature_grid(endpoints, gaps, node_count)
        matrix = np.zeros((degree, degree), dtype=float)
        nonfinite_count = 0
        for row in grid:
            values = barycentric_values(row["x"], collocation_nodes, bary_weights)
            if not np.all(np.isfinite(values)):
                nonfinite_count += 1
                continue
            matrix[int(row["gap_index"]), :] += row["period_weight"] * values
        note = (
            "gap-midpoint barycentric Lagrange basis built"
            if nonfinite_count == 0
            else f"gap-midpoint barycentric Lagrange had {nonfinite_count} nonfinite evaluation rows"
        )
        return matrix, note
    except Exception as exc:  # pragma: no cover - diagnostic guard
        return np.full((degree, degree), math.inf), f"gap-midpoint Lagrange failed: {exc}"


def basis_sweep(
    endpoints: list[float],
    gaps: list[tuple[float, float]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for node_count in QUADRATURE_NODES:
        for seed in ["monomial", "chebyshev_global", "chebyshev_gap_hull"]:
            seed_matrix, weighted_eval = build_seed_period_matrix_and_eval(
                seed, endpoints, gaps, node_count
            )
            eval_metrics = matrix_metrics(weighted_eval)
            append_candidate(
                rows,
                basis_name=seed,
                seed_basis=seed,
                transform_kind="none",
                node_count=node_count,
                matrix=seed_matrix,
                seed_eval_rank_f64=eval_metrics["rank_f64"],
                seed_eval_condition_f64=eval_metrics["condition_number_f64"],
                notes="raw seed basis control",
            )
            append_candidate(
                rows,
                basis_name=f"{seed}__row_column_scaled_diagnostic",
                seed_basis=seed,
                transform_kind="row_column_scaled_diagnostic_not_basis",
                node_count=node_count,
                matrix=row_column_scaled_matrix(seed_matrix),
                seed_eval_rank_f64=eval_metrics["rank_f64"],
                seed_eval_condition_f64=eval_metrics["condition_number_f64"],
                notes=(
                    "diagnostic only; row scaling is not a basis transform and "
                    "cannot be used by itself for the theorem"
                ),
                theorem_basis_candidate=False,
            )
            try:
                _q, r = np.linalg.qr(weighted_eval, mode="reduced")
                transform = np.linalg.inv(r)
                transformed = seed_matrix @ transform
                append_candidate(
                    rows,
                    basis_name=f"{seed}__weighted_qr_orthogonalized",
                    seed_basis=seed,
                    transform_kind="right_transform_from_weighted_quadrature_qr",
                    node_count=node_count,
                    matrix=transformed,
                    transform_norm_f64=float(np.linalg.norm(transform, 2)),
                    transform_condition_f64=float(np.linalg.cond(r)),
                    seed_eval_rank_f64=eval_metrics["rank_f64"],
                    seed_eval_condition_f64=eval_metrics["condition_number_f64"],
                    notes=(
                        "orthogonalized against weighted quadrature evaluations; "
                        "requires intervalization of the transform before use"
                    ),
                )
            except Exception as exc:
                append_candidate(
                    rows,
                    basis_name=f"{seed}__weighted_qr_orthogonalized",
                    seed_basis=seed,
                    transform_kind="right_transform_from_weighted_quadrature_qr",
                    node_count=node_count,
                    matrix=np.full_like(seed_matrix, math.inf),
                    seed_eval_rank_f64=eval_metrics["rank_f64"],
                    seed_eval_condition_f64=eval_metrics["condition_number_f64"],
                    notes=f"QR transform failed: {exc}",
                )
            try:
                u, s, vt = np.linalg.svd(seed_matrix, full_matrices=False)
                transform = vt.T @ np.diag(1.0 / s)
                transformed = seed_matrix @ transform
                append_candidate(
                    rows,
                    basis_name=f"{seed}__period_svd_right_preconditioned",
                    seed_basis=seed,
                    transform_kind="right_transform_from_period_matrix_svd",
                    node_count=node_count,
                    matrix=transformed,
                    transform_norm_f64=float(np.linalg.norm(transform, 2)),
                    transform_condition_f64=float(s[0] / s[-1]),
                    seed_eval_rank_f64=eval_metrics["rank_f64"],
                    seed_eval_condition_f64=eval_metrics["condition_number_f64"],
                    notes=(
                        "diagnostic right preconditioner derived from the period "
                        "matrix itself; useful as existence evidence, but heavy "
                        "to intervalize"
                    ),
                    theorem_basis_candidate=False,
                )
            except Exception as exc:
                append_candidate(
                    rows,
                    basis_name=f"{seed}__period_svd_right_preconditioned",
                    seed_basis=seed,
                    transform_kind="right_transform_from_period_matrix_svd",
                    node_count=node_count,
                    matrix=np.full_like(seed_matrix, math.inf),
                    seed_eval_rank_f64=eval_metrics["rank_f64"],
                    seed_eval_condition_f64=eval_metrics["condition_number_f64"],
                    notes=f"SVD transform failed: {exc}",
                )

        lagrange_matrix, lagrange_note = build_gap_midpoint_lagrange_matrix(
            endpoints, gaps, node_count
        )
        append_candidate(
            rows,
            basis_name="gap_midpoint_lagrange",
            seed_basis="gap_midpoint_lagrange",
            transform_kind="barycentric_lagrange_at_gap_midpoints",
            node_count=node_count,
            matrix=lagrange_matrix,
            notes=lagrange_note,
        )
    return rows


def best_candidate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    usable = [
        row
        for row in rows
        if row["theorem_basis_candidate"]
        and row["all_entries_finite"]
        and math.isfinite(row["condition_number_f64"])
    ]
    return min(usable, key=lambda row: row["condition_number_f64"])


def build_stage_b_status(best: dict[str, Any]) -> dict[str, Any]:
    return {
        "packet_id": EXPERIMENT_ID,
        "stage": "B",
        "target": "directed Rust/Inari interval condition certificate",
        "status": "PENDING_BACKEND_NOT_IMPLEMENTED",
        "best_f64_basis_name": best["basis_name"],
        "best_f64_condition_number": best["condition_number_f64"],
        "required_interval_quantities": [
            "interval enclosure of transformed 24x24 period matrix entries",
            "interval lower bound for smallest singular value or determinant plus inverse norm",
            "interval upper bound for condition number",
            "interval certificate for basis transform if transform_kind is not none",
        ],
        "acceptance_condition": f"interval-certified condition number < {CONDITION_THRESHOLD:g}",
        "claim_ceiling": (
            "Stage B status only. No interval condition certificate is produced "
            "by this packet."
        ),
    }


def build_best_basis_audit(best: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    node_counts = sorted(
        {
            row["quadrature_nodes_per_gap"]
            for row in rows
            if row["basis_name"] == best["basis_name"]
        }
    )
    stability_rows = [
        {
            "quadrature_nodes_per_gap": row["quadrature_nodes_per_gap"],
            "condition_number_f64": row["condition_number_f64"],
            "rank_f64": row["rank_f64"],
            "transform_condition_f64": row["transform_condition_f64"],
            "seed_eval_rank_f64": row["seed_eval_rank_f64"],
        }
        for row in rows
        if row["basis_name"] == best["basis_name"]
    ]
    conditions = [row["condition_number_f64"] for row in stability_rows]
    return {
        "packet_id": EXPERIMENT_ID,
        "best_basis": best,
        "quadrature_stability_node_counts": node_counts,
        "quadrature_stability_rows": stability_rows,
        "condition_ratio_max_over_min": None
        if not conditions
        else float(max(conditions) / min(conditions)),
        "basis_condition_pass_f64": best["condition_number_f64"] < CONDITION_THRESHOLD,
        "directed_interval_certificate_pending": True,
        "risk_note": (
            "The best f64 basis is transform-defined. A later interval packet "
            "must certify both the transform and the transformed period matrix "
            "before any residual audit can use it."
        ),
    }


def build_results(
    sources: dict[str, Any],
    rows: list[dict[str, Any]],
    best: dict[str, Any],
    best_audit: dict[str, Any],
) -> dict[str, Any]:
    f64_pass = best["condition_number_f64"] < CONDITION_THRESHOLD
    return {
        "packet_id": EXPERIMENT_ID,
        "version": "20260527-01",
        "created_utc": now_utc(),
        "route_type": "stable_gap_period_basis_conditioning_gate",
        "status": (
            "STABLE_GAP_PERIOD_BASIS_CONDITIONING_GATE_PARTIAL__F64_"
            "ORTHOGONALIZED_BASIS_FOUND__DIRECTED_INTERVAL_CERTIFICATE_PENDING"
            if f64_pass
            else "STABLE_GAP_PERIOD_BASIS_CONDITIONING_GATE_BLOCKED__"
            "NO_F64_BASIS_BELOW_THRESHOLD"
        ),
        "verdict": (
            "F64_STABLE_BASIS_FOUND__INTERVAL_CERTIFICATE_REQUIRED"
            if f64_pass
            else "NO_STABLE_F64_BASIS_FOUND"
        ),
        "claim_ceiling": CLAIM_CEILING,
        "altitude": "8525 m -- unchanged",
        "public_sota": "Problem open; no public SOTA improvement.",
        "source_packets": sources,
        "products": {
            "basis_sweep_rows": str(BASIS_SWEEP_ROWS),
            "best_basis_audit": str(BEST_BASIS_AUDIT_JSON),
            "stage_b_interval_certificate_status": str(STAGE_B_CERTIFICATE_JSON),
            "report": str(REPORT_MD),
        },
        "summary": {
            "basis_row_count": len(rows),
            "condition_threshold": CONDITION_THRESHOLD,
            "best_basis_name": best["basis_name"],
            "best_seed_basis": best["seed_basis"],
            "best_transform_kind": best["transform_kind"],
            "best_condition_number_f64": best["condition_number_f64"],
            "best_rank_f64": best["rank_f64"],
            "best_transform_condition_f64": best["transform_condition_f64"],
            "best_transform_norm_f64": best["transform_norm_f64"],
            "best_seed_eval_rank_f64": best["seed_eval_rank_f64"],
            "best_condition_ratio_max_over_min": best_audit[
                "condition_ratio_max_over_min"
            ],
        },
        "gate_checks": {
            "source_result_sidecars_ok": all(row["status"] == "OK" for row in sources.values()),
            "stage_a_f64_basis_below_threshold": f64_pass,
            "best_basis_full_rank_f64": best["rank_f64"] == 24,
            "best_basis_has_transform_warning": best[
                "transform_intervalization_warning"
            ],
            "best_basis_has_seed_rank_warning": best["seed_rank_warning"],
            "stage_b_directed_interval_certificate_performed": False,
            "stage_b_directed_interval_certificate_pending": True,
            "period_residual_audit_ready": False,
            "signed_period_legitimacy_proved": False,
            "altitude_reconsideration_allowed": False,
        },
        "hypothesis_status": {
            "STABLE_GAP_PERIOD_BASIS_F64": "FOUND_WITH_TRANSFORM_AUDIT_PENDING"
            if f64_pass
            else "NOT_FOUND",
            "DIRECTED_INTERVAL_CONDITION_CERTIFICATE": "PENDING",
            "INTERIOR_GAP_SOURCE_THEOREM": "PENDING",
            "GAP_PERIOD_VECTOR_RESIDUAL_AUDIT": "PENDING_AFTER_INTERVAL_BASIS_CERTIFICATE",
            "SIGNED_WITNESS_PERIOD_MATRIX_LEGITIMACY": "PENDING",
            "ATTAINED_WITNESS_TYPED_DUAL_MARGIN": "PENDING_PERIOD_LEGITIMACY_COMPOSITION",
            "KKT_COMPOSITION": "PENDING",
            "GLOBAL_REDUCTION": "UNTOUCHED",
        },
        "recommended_next_pitch": (
            "EXP-MATH-ERDOS1038-PHI-K-STABLE-GAP-PERIOD-BASIS-INTERVAL-"
            "CERTIFICATE-20260527-01"
            if f64_pass
            else "EXP-MATH-ERDOS1038-PHI-K-EXTENDED-PRECISION-GAP-PERIOD-"
            "BASIS-GATE-20260527-01"
        ),
        "next_blocker": (
            "Certify the best transformed basis with directed interval arithmetic "
            "before using it in a period residual audit."
            if f64_pass
            else "Move to extended precision or a different period-constraint frame."
        ),
    }


def build_report(results: dict[str, Any], best: dict[str, Any]) -> str:
    s = results["summary"]
    return f"""# Stable Gap-Period Basis Conditioning Gate

Packet: `{EXPERIMENT_ID}`

## Verdict

`{results["verdict"]}`

## Meaning

The naive period bases are still ill-conditioned, but an orthogonalized f64
basis brings the transformed `24`-row period matrix below the diagnostic
threshold.

```text
best_basis_name = {s["best_basis_name"]}
best_seed_basis = {s["best_seed_basis"]}
best_transform_kind = {s["best_transform_kind"]}
best_condition_number_f64 = {s["best_condition_number_f64"]}
best_rank_f64 = {s["best_rank_f64"]}
best_transform_condition_f64 = {s["best_transform_condition_f64"]}
best_seed_eval_rank_f64 = {s["best_seed_eval_rank_f64"]}
condition_threshold = {s["condition_threshold"]}
```

## Boundary

This is only Stage A. Stage B, the directed Rust/Inari interval certificate for
the transformed matrix and transform, is still pending.

## Best-Basis Warning

```text
transform_intervalization_warning = {best["transform_intervalization_warning"]}
seed_rank_warning = {best["seed_rank_warning"]}
notes = {best["notes"]}
```

The best f64 basis is useful enough to target next, but not enough to run a
period residual audit.

## Products

```text
{BASIS_SWEEP_ROWS}
{BEST_BASIS_AUDIT_JSON}
{STAGE_B_CERTIFICATE_JSON}
```

## Claim Ceiling

{CLAIM_CEILING}

Altitude remains 8525 m. Public SOTA is unchanged and #1038 remains open.
"""


def main() -> None:
    sources = load_sources()
    components = read_jsonl(packet_path(FULL_PANEL_PACKET_ID, "COMPONENT_ROWS.jsonl"))
    endpoints = interior_pilot.endpoints_from_components(components)
    gaps = interior_pilot.gaps_from_components(components)
    rows = basis_sweep(endpoints, gaps)
    best = best_candidate(rows)
    best_audit = build_best_basis_audit(best, rows)
    stage_b = build_stage_b_status(best)
    results = build_results(sources, rows, best, best_audit)

    write_unique(
        BASIS_SWEEP_ROWS,
        "\n".join(json.dumps(row, sort_keys=True) for row in rows) + "\n",
    )
    write_unique(BEST_BASIS_AUDIT_JSON, json.dumps(best_audit, indent=2, sort_keys=True) + "\n")
    write_unique(STAGE_B_CERTIFICATE_JSON, json.dumps(stage_b, indent=2, sort_keys=True) + "\n")
    write_unique(REPORT_MD, build_report(results, best))
    write_unique(RESULTS_JSON, json.dumps(results, indent=2, sort_keys=True) + "\n")
    write_unique(SHA_FILE, f"{sha256_path(RESULTS_JSON)}  {RESULTS_JSON.name}\n")

    print(json.dumps({"packet_id": EXPERIMENT_ID, "status": results["status"]}, sort_keys=True))


if __name__ == "__main__":
    main()
