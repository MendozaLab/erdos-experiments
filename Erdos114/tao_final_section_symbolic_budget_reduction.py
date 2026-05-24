"""Emit an EHP114 final-section symbolic budget reduction packet.

L40A named the first Tao final-section blockers. This pass takes the next
proof-engineering step: reduce the final glue algebra to symbolic inequalities
that can later receive numeric constants. It deliberately emits no global
threshold and authorizes no higher-degree computation.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-TAO-FINAL-SECTION-SYMBOLIC-BUDGET-REDUCTION-20260507-01"
L40A_ID = "EXP-MATH-EHP114-TAO-FINAL-SECTION-CONSTANT-EXTRACTION-20260507-01"
L40B_ID = "EXP-MATH-EHP114-TAO-FINAL-SECTION-THRESHOLD-SUBCHECK-20260507-01"
EXPECTED_TAO_SOURCE_SHA = "ab792bd8fe806a985e44f692ccb8305cf386677509846979bb8ff35a276814a6"


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


def source_lines() -> tuple[list[str], str]:
    module = load_l40_module()
    tex = module.extract_lemniscate_tex(module.fetch_arxiv_source())
    source_sha = hashlib.sha256(tex.encode("utf-8")).hexdigest()
    return tex.splitlines(), source_sha


def snippet(lines: list[str], start: int, end: int) -> str:
    text = " ".join(line.strip() for line in lines[start - 1 : end])
    return " ".join(text.split())


def final_deficit_rows(l40a: dict[str, Any]) -> list[dict[str, Any]]:
    wanted = {"final-deficit-inside-2", "final-deficit-annulus-2", "final-deficit-outside-again"}
    rows = []
    for row in l40a["extraction_rows"]:
        if row["dependency_id"] in wanted:
            rows.append(row)
    rows.sort(key=lambda row: row["dependency_id"])
    return rows


def symbolic_conditions() -> list[dict[str, Any]]:
    return [
        {
            "condition_id": "choose_C_for_pp0_refinement",
            "source_rows": ["needle-absolute-constant-1747", "needle-large-enough-1770"],
            "symbolic_condition": "C >= C_origin_repulsion(c_in, C0_independent_constants)",
            "meaning": "The extra inner-region gain is available only when the origin-repulsion input makes the small disk empty under pp0.",
            "numeric_status": "PARENT_CONSTANTS_OPAQUE",
            "blocks": ["inside-2-second", "annulus-2-repeat"],
        },
        {
            "condition_id": "choose_C0_for_pp0_case",
            "source_rows": ["final-deficit-inside-2", "final-deficit-annulus-2", "final-deficit-outside-again"],
            "symbolic_condition": "C0 >= 2*K_pp0/c_pp0",
            "meaning": "The combined O(||p||/C0) terms must be less than half of the negative pp0-case linear gain.",
            "numeric_status": "PARENT_CONSTANTS_OPAQUE",
            "blocks": ["a-ineq-pp0-case"],
        },
        {
            "condition_id": "choose_C0_for_not_pp0_case",
            "source_rows": ["final-deficit-inside-2", "final-deficit-annulus-2"],
            "symbolic_condition": "C0 >= 2*K_not_pp0*C/c_disp",
            "meaning": "When pp0 fails, ||p|| < C||p||_1 converts the remainder into a ||p||_1 cost; C0 must make that cost smaller than the dispersion gain.",
            "numeric_status": "PARENT_CONSTANTS_OPAQUE",
            "blocks": ["a-ineq-not-pp0-case"],
        },
        {
            "condition_id": "choose_C_for_pots_contradiction",
            "source_rows": ["needle-sufficiently-large-1711", "needle-large-enough-1717", "needle-large-enough-1719"],
            "symbolic_condition": "C/2 > K_pots_upper and n >= N_pots(C, eps)",
            "meaning": "The lower bound for Psi(D(0,C eps^(1/2)/n)) must exceed the upper bound forced by the comparison with the reference lemniscate.",
            "numeric_status": "PARENT_CONSTANTS_OPAQUE",
            "blocks": ["pots"],
        },
        {
            "condition_id": "choose_epsilon_for_ets",
            "source_rows": ["needle-sufficiently-large-1649", "needle-absolute-constant-1655", "needle-sufficiently-large-1661"],
            "symbolic_condition": "eps <= eps0 and n >= N_ets(eps)",
            "meaning": "The inherited inside/annulus/outside estimates must be quantitative enough to turn ||p||_1 into a prescribed small quantity.",
            "numeric_status": "PARENT_CONSTANTS_OPAQUE",
            "blocks": ["ets"],
        },
    ]


def final_case_rows(lines: list[str]) -> list[dict[str, Any]]:
    return [
        {
            "case_id": "pp0_holds",
            "source_span": [1867, 1869],
            "source_snippet": snippet(lines, 1867, 1869),
            "needed_parent_rows": ["inside-2-second", "annulus-2", "outside-again", "split-2"],
            "symbolic_deficit": "c_pp0*||p|| - K_pp0*||p||/C0",
            "sufficient_condition": "C0 >= 2*K_pp0/c_pp0",
            "case_status": "SYMBOLIC_BUDGET_READY_NUMERIC_CONSTANTS_PENDING",
        },
        {
            "case_id": "pp0_fails",
            "source_span": [1869, 1873],
            "source_snippet": snippet(lines, 1869, 1873),
            "needed_parent_rows": ["inside-2-first", "annulus-2", "disp-split", "failure-of-pp0"],
            "symbolic_deficit": "c_disp*||p||_1 - K_not_pp0*||p||/C0",
            "sufficient_condition": "C0 >= 2*K_not_pp0*C/c_disp",
            "case_status": "SYMBOLIC_BUDGET_READY_NUMERIC_CONSTANTS_PENDING",
        },
    ]


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
        "# EHP114 Tao Final-Section Symbolic Budget Reduction\n\n"
        f"Experiment: `{EXPERIMENT_ID}`\n\n"
        "## Verdict\n\n"
        f"Status: `{result['status']}`.\n\n"
        "This packet reduces the final-section glue algebra to symbolic budget "
        "conditions. It does not assign numeric constants, does not emit a "
        "global threshold, and does not authorize computation above the finite "
        "packet range.\n\n"
        "## Summary\n\n"
        f"- Final deficit rows: `{result['final_deficit_row_count']}`\n"
        f"- Final case rows: `{result['final_case_row_count']}`\n"
        f"- Symbolic conditions: `{result['symbolic_condition_count']}`\n"
        f"- Numeric-ready conditions: `{result['numeric_ready_condition_count']}`\n"
        f"- Source SHA status: `{result['tao_source']['source_sha_status']}`\n\n"
        "## Next Action\n\n"
        f"{result['next_action']}\n\n"
        "## Claim Ceiling\n\n"
        f"{result['claim_ceiling']}\n"
    )


def main() -> None:
    proof_root = proof_path_root()
    l40a_path = proof_root / L40A_ID / f"{L40A_ID}_RESULTS.json"
    l40b_path = proof_root / L40B_ID / f"{L40B_ID}_RESULTS.json"
    l40a = read_json(l40a_path)
    l40b = read_json(l40b_path)
    l40a_sha = verify_result_sha("final_section_constant_extraction", l40a_path)
    l40b_sha = verify_result_sha("final_section_threshold_subcheck", l40b_path)

    lines, tao_sha = source_lines()
    source_sha_status = "PASS_PINNED_SOURCE_SHA_MATCH" if tao_sha == EXPECTED_TAO_SOURCE_SHA else "FAIL_PINNED_SOURCE_SHA_MISMATCH"

    deficit_rows = final_deficit_rows(l40a)
    conditions = symbolic_conditions()
    cases = final_case_rows(lines)
    numeric_ready_count = sum(1 for row in conditions if row["numeric_status"] == "NUMERIC_READY")
    symbolic_ready_count = len(conditions) - numeric_ready_count

    status = (
        "FINAL_SECTION_SYMBOLIC_BUDGET_READY_NUMERIC_CONSTANTS_PENDING"
        if source_sha_status == "PASS_PINNED_SOURCE_SHA_MATCH" and len(deficit_rows) == 3
        else "FINAL_SECTION_SYMBOLIC_BUDGET_BLOCKED_SOURCE_OR_ROWS"
    )

    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "source_experiment_ids": [L40A_ID, L40B_ID],
        "source_artifacts": [
            {"artifact_id": L40A_ID, "path": str(l40a_path), "status": l40a["status"], "sha_check": l40a_sha},
            {"artifact_id": L40B_ID, "path": str(l40b_path), "status": l40b["status"], "sha_check": l40b_sha},
        ],
        "tao_source": {
            "expected_source_sha256": EXPECTED_TAO_SOURCE_SHA,
            "source_sha256": tao_sha,
            "source_sha_status": source_sha_status,
            "source_line_count": len(lines),
        },
        "final_deficit_row_count": len(deficit_rows),
        "final_deficit_rows": deficit_rows,
        "final_case_row_count": len(cases),
        "final_case_rows": cases,
        "symbolic_condition_count": len(conditions),
        "symbolic_ready_condition_count": symbolic_ready_count,
        "numeric_ready_condition_count": numeric_ready_count,
        "symbolic_conditions": conditions,
        "partial_max_N_final_section": None,
        "global_candidate_N0": None,
        "n15_or_higher_computation_authorized": False,
        "full_synthesis_allowed": False,
        "first_failed_condition": "numeric parent constants remain opaque" if numeric_ready_count == 0 else "none",
        "next_action": "Quantify parent constants K_pp0, K_not_pp0, c_pp0, c_disp, K_pots_upper, N_pots, eps0, and N_ets before any threshold integer can be emitted.",
        "claim_ceiling": "Symbolic final-section budget reduction only. This packet is not a full all-degree proof and does not authorize n=15 or higher computation.",
    }

    outdir = proof_root / EXPERIMENT_ID
    digest = write_artifact(outdir, result, report_for(result))
    print(
        json.dumps(
            {
                "status": "EHP114_FINAL_SECTION_SYMBOLIC_BUDGET_EMITTED",
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
