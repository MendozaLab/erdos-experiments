"""Emit an EHP114 Tao final-section symbolic parameter propagation packet.

L47 reconciles the L46 pvar result with the earlier final-section symbolic
budget. Its narrow purpose is to upgrade K_origin from opaque to symbolic-ready
inside the bridge table while preserving the remaining blockers and refusing to
emit a global threshold.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-TAO-FINAL-SECTION-SYMBOLIC-PARAMETER-PROPAGATION-20260507-01"
L46_ID = "EXP-MATH-EHP114-TAO-PVAR-ELEMENTARY-INEQUALITY-PACKET-20260507-01"
L44_ID = "EXP-MATH-EHP114-TAO-BASE-LEMMA-SEED-PROPAGATION-20260507-02"
SYMBOLIC_ID = "EXP-MATH-EHP114-TAO-FINAL-SECTION-SYMBOLIC-BUDGET-REDUCTION-20260507-01"
EXPECTED_TAO_SOURCE_SHA = "ab792bd8fe806a985e44f692ccb8305cf386677509846979bb8ff35a276814a6"

NUMERIC_READY = "NUMERIC_READY"
SYMBOLIC_READY = "SYMBOLIC_READY"
OPAQUE = "OPAQUE"
SOURCE_BLOCKED = "SOURCE_BLOCKED"


def erdos114_root() -> Path:
    return Path(__file__).resolve().parent


def proof_path_root() -> Path:
    return erdos114_root() / "proof_path"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha_path_for(path: Path) -> Path:
    if not path.name.endswith("_RESULTS.json"):
        raise ValueError(f"not a result JSON path: {path}")
    return path.with_name(path.name.replace("_RESULTS.json", "_RESULTS.sha256"))


def verify_result_sha(role: str, path: Path) -> dict[str, Any]:
    sha_path = sha_path_for(path)
    expected = sha_path.read_text().split()[0]
    actual = sha256_file(path)
    return {
        "role": role,
        "path": str(path),
        "sha_path": str(sha_path),
        "expected_sha256": expected,
        "actual_sha256": actual,
        "sha_status": "PASS" if expected == actual else "FAIL",
    }


def load_l40_module() -> Any:
    path = erdos114_root() / "tao_final_section_constant_extraction.py"
    spec = importlib.util.spec_from_file_location("tao_l40_extraction", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load L40 extraction module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def tao_source_sha() -> str:
    module = load_l40_module()
    tex = module.extract_lemniscate_tex(module.fetch_arxiv_source())
    return hashlib.sha256(tex.encode("utf-8")).hexdigest()


def l44_base_index(l44: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {row["base_constant"]: row for row in l44["propagated_base_rows"]}


def l46_k_origin_ready(l46: dict[str, Any]) -> bool:
    return (
        l46.get("status") == "PVAR_ELEMENTARY_INEQUALITY_READY"
        and l46.get("k_origin_status_after_elementary_packet") == "K_ORIGIN_SYMBOLIC_READY_NUMERIC_PARAMETERS_PENDING"
        and l46.get("opaque_elementary_count") == 0
    )


def propagated_parent_rows(l44: dict[str, Any], l46: dict[str, Any]) -> list[dict[str, Any]]:
    base = l44_base_index(l44)
    rows: list[dict[str, Any]] = []

    k_disp = base.get("K_disp_split", {})
    rows.append(
        {
            "parameter": "K_disp_split",
            "parameter_status": NUMERIC_READY if k_disp.get("base_status") == "NUMERIC_BASE_READY" else SOURCE_BLOCKED,
            "source": L44_ID,
            "expression": k_disp.get("propagated_expression"),
            "remaining_blockers": k_disp.get("remaining_blockers", ["missing K_disp_split row"]),
        }
    )

    c_defect = base.get("c_defect_C", {})
    rows.append(
        {
            "parameter": "c_defect_C",
            "parameter_status": SYMBOLIC_READY if c_defect.get("base_status") == "SYMBOLIC_BASE_READY" else SOURCE_BLOCKED,
            "source": L44_ID,
            "expression": c_defect.get("propagated_expression"),
            "remaining_blockers": c_defect.get("remaining_blockers", ["missing c_defect_C row"]),
        }
    )

    k_ankh = base.get("K_ankh", {})
    rows.append(
        {
            "parameter": "K_ankh",
            "parameter_status": SYMBOLIC_READY if k_ankh.get("base_status") == "SYMBOLIC_BASE_READY" else SOURCE_BLOCKED,
            "source": L44_ID,
            "expression": k_ankh.get("propagated_expression"),
            "remaining_blockers": k_ankh.get("remaining_blockers", ["missing K_ankh row"]),
        }
    )

    if l46_k_origin_ready(l46):
        rows.append(
            {
                "parameter": "K_origin",
                "parameter_status": SYMBOLIC_READY,
                "source": L46_ID,
                "expression": "K_origin = F(A_pocl, B1, B2, C, C_prime)",
                "remaining_blockers": ["numeric selection of A_pocl, B1, B2, C, and C_prime"],
            }
        )
    else:
        rows.append(
            {
                "parameter": "K_origin",
                "parameter_status": OPAQUE,
                "source": L46_ID,
                "expression": None,
                "remaining_blockers": ["K_origin not symbolic-ready in L46"],
            }
        )

    return rows


def final_condition_rows(symbolic: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in symbolic["symbolic_conditions"]:
        condition_id = row["condition_id"]
        if condition_id == "choose_C_for_pp0_refinement":
            status = SYMBOLIC_READY
            blockers = ["numeric C_origin_repulsion from symbolic K_origin parameters"]
        elif condition_id in {"choose_C0_for_pp0_case", "choose_C0_for_not_pp0_case"}:
            status = SYMBOLIC_READY
            blockers = ["numeric K_pp0, c_pp0, K_not_pp0, c_disp still pending"]
        elif condition_id == "choose_C_for_pots_contradiction":
            status = OPAQUE
            blockers = ["K_pots_upper and N_pots(C,eps) remain opaque"]
        elif condition_id == "choose_epsilon_for_ets":
            status = OPAQUE
            blockers = ["eps0 and N_ets(eps) remain opaque"]
        else:
            status = OPAQUE
            blockers = ["unclassified symbolic condition"]
        rows.append(
            {
                "condition_id": condition_id,
                "condition_status": status,
                "symbolic_condition": row["symbolic_condition"],
                "remaining_blockers": blockers,
            }
        )
    return rows


def write_artifact(outdir: Path, result: dict[str, Any], report: str) -> str:
    outdir.mkdir(parents=True, exist_ok=False)
    result_path = outdir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = outdir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = outdir / f"{EXPERIMENT_ID}_RESULTS.sha256"
    result_text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    result_path.write_text(result_text)
    report_path.write_text(report)
    digest = hashlib.sha256(result_text.encode("utf-8")).hexdigest()
    sha_path.write_text(f"{digest}  {result_path.name}\n")
    return digest


def report_for(result: dict[str, Any]) -> str:
    return (
        "# EHP114 Tao Final-Section Symbolic Parameter Propagation\n\n"
        f"Experiment: `{EXPERIMENT_ID}`\n\n"
        "## Verdict\n\n"
        f"Status: `{result['status']}`.\n\n"
        "This packet propagates L46 back into the final-section symbolic budget. "
        "It upgrades K_origin to symbolic-ready but emits no global threshold.\n\n"
        "## Summary\n\n"
        f"- Parent parameters: `{result['parent_parameter_count']}`\n"
        f"- Numeric-ready parent parameters: `{result['numeric_parent_ready_count']}`\n"
        f"- Symbolic-ready parent parameters: `{result['symbolic_parent_ready_count']}`\n"
        f"- Opaque parent parameters: `{result['opaque_parent_count']}`\n"
        f"- Final symbolic conditions ready: `{result['symbolic_condition_ready_count']}`\n"
        f"- Tao source SHA status: `{result['tao_source_sha_status']}`\n\n"
        "## First Failed Condition\n\n"
        f"`{result['first_failed_condition']}`\n\n"
        "## Claim Ceiling\n\n"
        f"{result['claim_ceiling']}\n"
    )


def main() -> None:
    proof_root = proof_path_root()
    l46_path = proof_root / L46_ID / f"{L46_ID}_RESULTS.json"
    l44_path = proof_root / L44_ID / f"{L44_ID}_RESULTS.json"
    symbolic_path = proof_root / SYMBOLIC_ID / f"{SYMBOLIC_ID}_RESULTS.json"
    l46 = read_json(l46_path)
    l44 = read_json(l44_path)
    symbolic = read_json(symbolic_path)
    l46_sha = verify_result_sha("pvar_elementary_inequality_packet", l46_path)
    l44_sha = verify_result_sha("base_lemma_seed_propagation", l44_path)
    symbolic_sha = verify_result_sha("final_section_symbolic_budget", symbolic_path)

    source_sha = tao_source_sha()
    tao_source_sha_status = "PASS_PINNED_SOURCE_SHA_MATCH" if source_sha == EXPECTED_TAO_SOURCE_SHA else "FAIL_PINNED_SOURCE_SHA_MISMATCH"
    parent_rows = propagated_parent_rows(l44, l46)
    condition_rows = final_condition_rows(symbolic)

    numeric_parent_count = sum(1 for row in parent_rows if row["parameter_status"] == NUMERIC_READY)
    symbolic_parent_count = sum(1 for row in parent_rows if row["parameter_status"] == SYMBOLIC_READY)
    opaque_parent_count = sum(1 for row in parent_rows if row["parameter_status"] == OPAQUE)
    source_blocked_count = sum(1 for row in parent_rows if row["parameter_status"] == SOURCE_BLOCKED)
    condition_ready_count = sum(1 for row in condition_rows if row["condition_status"] == SYMBOLIC_READY)
    condition_opaque_count = sum(1 for row in condition_rows if row["condition_status"] == OPAQUE)

    if tao_source_sha_status != "PASS_PINNED_SOURCE_SHA_MATCH" or source_blocked_count > 0:
        status = "FINAL_SECTION_SYMBOLIC_PARAMETER_PROPAGATION_BLOCKED_SOURCE"
    elif opaque_parent_count == 0 and condition_opaque_count == 0:
        status = "FINAL_SECTION_SYMBOLIC_PARAMETER_PROPAGATION_READY"
    elif symbolic_parent_count + numeric_parent_count > 0:
        status = "FINAL_SECTION_SYMBOLIC_PARAMETER_PROPAGATION_PARTIAL"
    else:
        status = "FINAL_SECTION_SYMBOLIC_PARAMETER_PROPAGATION_REMAINS_OPAQUE"

    if source_blocked_count:
        first_failed = "source parent parameter missing"
    elif opaque_parent_count:
        first_failed = "opaque parent parameter remains"
    elif condition_opaque_count:
        first_failed = "pots and ets final-section conditions remain opaque"
    else:
        first_failed = "none"

    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "source_experiment_ids": [L46_ID, L44_ID, SYMBOLIC_ID],
        "source_artifacts": [
            {"artifact_id": L46_ID, "path": str(l46_path), "status": l46["status"], "sha_check": l46_sha},
            {"artifact_id": L44_ID, "path": str(l44_path), "status": l44["status"], "sha_check": l44_sha},
            {"artifact_id": SYMBOLIC_ID, "path": str(symbolic_path), "status": symbolic["status"], "sha_check": symbolic_sha},
        ],
        "tao_source_sha_status": tao_source_sha_status,
        "tao_source_sha256": source_sha,
        "expected_tao_source_sha256": EXPECTED_TAO_SOURCE_SHA,
        "parent_parameter_count": len(parent_rows),
        "parent_parameter_rows": parent_rows,
        "numeric_parent_ready_count": numeric_parent_count,
        "symbolic_parent_ready_count": symbolic_parent_count,
        "opaque_parent_count": opaque_parent_count,
        "source_blocked_count": source_blocked_count,
        "final_condition_count": len(condition_rows),
        "final_condition_rows": condition_rows,
        "symbolic_condition_ready_count": condition_ready_count,
        "opaque_condition_count": condition_opaque_count,
        "k_origin_bridge_status": "K_ORIGIN_SYMBOLIC_READY_NUMERIC_PARAMETERS_PENDING" if l46_k_origin_ready(l46) else "K_ORIGIN_NOT_READY",
        "partial_max_N_final_section": None,
        "global_candidate_N0": None,
        "n15_or_higher_computation_authorized": False,
        "full_synthesis_allowed": False,
        "first_failed_condition": first_failed,
        "claim_ceiling": "Final-section symbolic parameter propagation only. This packet does not provide a numeric Tao threshold, does not authorize n=15 or higher computation, and does not prove the full all-degree statement.",
    }

    outdir = proof_root / EXPERIMENT_ID
    digest = write_artifact(outdir, result, report_for(result))
    print(
        json.dumps(
            {
                "status": "EHP114_FINAL_SECTION_SYMBOLIC_PARAMETER_PROPAGATION_EMITTED",
                "experiment_id": EXPERIMENT_ID,
                "artifact_status": status,
                "sha256": digest,
                "outdir": str(outdir),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
