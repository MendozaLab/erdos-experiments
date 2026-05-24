"""Emit an EHP114 Tao base-lemma seed propagation packet.

L44 propagates the primitive seeds from L43 back into the four L42 base
constants. It upgrades only constants supported by the primitive seed audit and
leaves the pvar exponential constants opaque.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-TAO-BASE-LEMMA-SEED-PROPAGATION-20260507-02"
L42_ID = "EXP-MATH-EHP114-TAO-BASE-LEMMA-CONSTANT-AUDIT-20260507-01"
L43_ID = "EXP-MATH-EHP114-TAO-PRIMITIVE-CONSTANT-SEED-AUDIT-20260507-01"
EXPECTED_TAO_SOURCE_SHA = "ab792bd8fe806a985e44f692ccb8305cf386677509846979bb8ff35a276814a6"

NUMERIC_BASE_READY = "NUMERIC_BASE_READY"
SYMBOLIC_BASE_READY = "SYMBOLIC_BASE_READY"
OPAQUE_BASE = "OPAQUE_BASE"
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


def primitive_index(l43: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {row["label"]: row for row in l43["primitive_seed_rows"]}


def has_seed(seeds: dict[str, dict[str, Any]], label: str, statuses: set[str]) -> bool:
    row = seeds.get(label)
    return bool(row and row["seed_status"] in statuses)


def propagated_rows(seeds: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    ready = {"NUMERIC_SEED_READY", "SYMBOLIC_SEED_READY"}
    rows: list[dict[str, Any]] = []

    if all(has_seed(seeds, label, ready) for label in ["markov-many", "disp-upper", "semi-norm"]):
        rows.append(
            {
                "base_constant": "K_disp_split",
                "base_status": NUMERIC_BASE_READY,
                "source_base_lemma": "disp-split",
                "propagated_expression": "||p||_1 <= 2*outside_mass + local_dispersion <= 2*outside_mass + n*local_dispersion for n>=4",
                "numeric_seed": {
                    "outside_mass_coefficient": 2.0,
                    "n_times_local_dispersion_coefficient": 1.0,
                    "valid_for_n_min": 4,
                },
                "used_primitive_labels": ["markov-many", "disp-upper", "semi-norm"],
                "remaining_blockers": [],
            }
        )
    else:
        rows.append(
            {
                "base_constant": "K_disp_split",
                "base_status": SOURCE_BLOCKED,
                "source_base_lemma": "disp-split",
                "propagated_expression": None,
                "numeric_seed": None,
                "used_primitive_labels": ["markov-many", "disp-upper", "semi-norm"],
                "remaining_blockers": ["primitive seed missing"],
            }
        )

    if all(has_seed(seeds, label, ready) for label in ["tridef", "markov-many", "disp-sum", "multip"]):
        rows.append(
            {
                "base_constant": "c_defect_C",
                "base_status": SYMBOLIC_BASE_READY,
                "source_base_lemma": "defect-psi-cor",
                "propagated_expression": "c_defect_C = F(c_tridef(C), c_markov_many, 2, 2*pi)",
                "numeric_seed": None,
                "used_primitive_labels": ["tridef", "markov-many", "disp-sum", "multip"],
                "remaining_blockers": ["numeric evaluation of c_tridef(C)"],
            }
        )
    else:
        rows.append(
            {
                "base_constant": "c_defect_C",
                "base_status": SOURCE_BLOCKED,
                "source_base_lemma": "defect-psi-cor",
                "propagated_expression": None,
                "numeric_seed": None,
                "used_primitive_labels": ["tridef", "markov-many", "disp-sum", "multip"],
                "remaining_blockers": ["primitive seed missing"],
            }
        )

    if all(has_seed(seeds, label, ready) for label in ["Psi-def", "semi-norm", "riesz-bound"]):
        rows.append(
            {
                "base_constant": "K_ankh",
                "base_status": SYMBOLIC_BASE_READY,
                "source_base_lemma": "ankh",
                "propagated_expression": "K_ankh = F(2*pi, 1/pi, pointwise_annulus_constant)",
                "numeric_seed": None,
                "used_primitive_labels": ["Psi-def", "semi-norm", "riesz-bound"],
                "exact_source_identity_labels": ["psi-expand"],
                "remaining_blockers": ["pointwise annulus estimate constant in the proof of ankh"],
            }
        )
    else:
        rows.append(
            {
                "base_constant": "K_ankh",
                "base_status": SOURCE_BLOCKED,
                "source_base_lemma": "ankh",
                "propagated_expression": None,
                "numeric_seed": None,
                "used_primitive_labels": ["psi-expand", "Psi-def", "semi-norm", "riesz-bound"],
                "remaining_blockers": ["primitive seed missing"],
            }
        )

    rows.append(
        {
            "base_constant": "K_origin",
            "base_status": OPAQUE_BASE,
            "source_base_lemma": "origin-repulsion",
            "propagated_expression": "K_origin depends on pvar constants in poz and pocl before origin repulsion can be quantified",
            "numeric_seed": None,
            "used_primitive_labels": ["poz", "npz", "pocl", "size-def"],
            "remaining_blockers": [
                "poz exponential pvar constant",
                "pocl exponential pvar constant",
            ],
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
        "# EHP114 Tao Base-Lemma Seed Propagation\n\n"
        f"Experiment: `{EXPERIMENT_ID}`\n\n"
        "## Verdict\n\n"
        f"Status: `{result['status']}`.\n\n"
        "This packet propagates primitive constant seeds back into the four "
        "base constants. It does not emit a global threshold and does not "
        "authorize higher-degree computation.\n\n"
        "## Summary\n\n"
        f"- Processed base constants: `{result['processed_base_constant_count']}`\n"
        f"- Numeric-ready base constants: `{result['numeric_base_ready_count']}`\n"
        f"- Symbolic-ready base constants: `{result['symbolic_base_ready_count']}`\n"
        f"- Opaque base constants: `{result['opaque_base_count']}`\n"
        f"- Tao source SHA status: `{result['tao_source_sha_status']}`\n\n"
        "## First Failed Condition\n\n"
        f"`{result['first_failed_condition']}`\n\n"
        "## Claim Ceiling\n\n"
        f"{result['claim_ceiling']}\n"
    )


def main() -> None:
    proof_root = proof_path_root()
    l42_path = proof_root / L42_ID / f"{L42_ID}_RESULTS.json"
    l43_path = proof_root / L43_ID / f"{L43_ID}_RESULTS.json"
    l42 = read_json(l42_path)
    l43 = read_json(l43_path)
    l42_sha = verify_result_sha("base_lemma_constant_audit", l42_path)
    l43_sha = verify_result_sha("primitive_constant_seed_audit", l43_path)

    source_sha = tao_source_sha()
    tao_source_sha_status = "PASS_PINNED_SOURCE_SHA_MATCH" if source_sha == EXPECTED_TAO_SOURCE_SHA else "FAIL_PINNED_SOURCE_SHA_MISMATCH"
    rows = propagated_rows(primitive_index(l43))

    numeric_count = sum(1 for row in rows if row["base_status"] == NUMERIC_BASE_READY)
    symbolic_count = sum(1 for row in rows if row["base_status"] == SYMBOLIC_BASE_READY)
    opaque_count = sum(1 for row in rows if row["base_status"] == OPAQUE_BASE)
    source_blocked_count = sum(1 for row in rows if row["base_status"] == SOURCE_BLOCKED)

    if tao_source_sha_status != "PASS_PINNED_SOURCE_SHA_MATCH" or source_blocked_count > 0:
        status = "BASE_LEMMA_SEED_PROPAGATION_BLOCKED_SOURCE"
    elif opaque_count == 0:
        status = "BASE_LEMMA_SEED_PROPAGATION_READY"
    elif numeric_count + symbolic_count > 0:
        status = "BASE_LEMMA_SEED_PROPAGATION_PARTIAL"
    else:
        status = "BASE_LEMMA_SEED_PROPAGATION_REMAINS_OPAQUE"

    if source_blocked_count:
        first_failed = "primitive seed missing"
    elif opaque_count:
        first_failed = "K_origin remains blocked by poz/pocl exponential pvar constants"
    else:
        first_failed = "none"

    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "source_experiment_ids": [L42_ID, L43_ID],
        "source_artifacts": [
            {"artifact_id": L42_ID, "path": str(l42_path), "status": l42["status"], "sha_check": l42_sha},
            {"artifact_id": L43_ID, "path": str(l43_path), "status": l43["status"], "sha_check": l43_sha},
        ],
        "tao_source_sha_status": tao_source_sha_status,
        "tao_source_sha256": source_sha,
        "expected_tao_source_sha256": EXPECTED_TAO_SOURCE_SHA,
        "processed_base_constant_count": len(rows),
        "propagated_base_rows": rows,
        "numeric_base_ready_count": numeric_count,
        "symbolic_base_ready_count": symbolic_count,
        "opaque_base_count": opaque_count,
        "source_blocked_count": source_blocked_count,
        "first_failed_condition": first_failed,
        "partial_max_N_final_section": None,
        "global_candidate_N0": None,
        "n15_or_higher_computation_authorized": False,
        "full_synthesis_allowed": False,
        "claim_ceiling": "Base-lemma seed propagation only. This packet does not provide a numeric Tao threshold, does not authorize n=15 or higher computation, and does not prove the full all-degree statement.",
    }

    outdir = proof_root / EXPERIMENT_ID
    digest = write_artifact(outdir, result, report_for(result))
    print(
        json.dumps(
            {
                "status": "EHP114_BASE_LEMMA_SEED_PROPAGATION_EMITTED",
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
