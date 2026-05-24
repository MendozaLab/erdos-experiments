"""Emit an EHP114 Tao parent-constant source map.

This packet follows the symbolic budget reduction one layer backward. It maps
the symbolic final-section constants to exact Tao source lemmas and proof
phrases, while preserving the fact that no numeric threshold has yet been
extracted.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-TAO-PARENT-CONSTANT-SOURCE-MAP-20260507-01"
SYMBOLIC_ID = "EXP-MATH-EHP114-TAO-FINAL-SECTION-SYMBOLIC-BUDGET-REDUCTION-20260507-01"
L40A_ID = "EXP-MATH-EHP114-TAO-FINAL-SECTION-CONSTANT-EXTRACTION-20260507-01"
EXPECTED_TAO_SOURCE_SHA = "ab792bd8fe806a985e44f692ccb8305cf386677509846979bb8ff35a276814a6"

PARENT_LABELS = [
    "origin-repulsion",
    "disp-split",
    "defect-psi-cor",
    "ankh",
    "inside",
    "annulus",
    "out",
    "x2-lem",
    "x3-lem",
    "sting",
    "arclength",
    "stokes",
]


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


def source_lines() -> tuple[list[str], str, Any]:
    module = load_l40_module()
    tex = module.extract_lemniscate_tex(module.fetch_arxiv_source())
    source_sha = hashlib.sha256(tex.encode("utf-8")).hexdigest()
    return tex.splitlines(), source_sha, module


def compact(lines: list[str], start: int, end: int, max_chars: int = 1100) -> str:
    text = " ".join(line.strip() for line in lines[start - 1 : end])
    text = re.sub(r"\s+", " ", text)
    if len(text) > max_chars:
        return text[: max_chars - 3] + "..."
    return text


def source_parent_rows(lines: list[str], module: Any) -> list[dict[str, Any]]:
    rows = []
    for label in PARENT_LABELS:
        start, end, env, title = module.find_env(lines, label)
        rows.append(
            {
                "label": label,
                "source_span": [start, end],
                "environment": env,
                "title": title,
                "source_snippet": compact(lines, start, end),
                "quantification_status": "SOURCE_MAPPED_NUMERIC_CONSTANTS_OPAQUE",
            }
        )
    return rows


def constant_rows() -> list[dict[str, Any]]:
    return [
        {
            "constant_symbol": "K_pp0",
            "role": "coefficient of the combined O(||p||/C0) remainder in the pp0 case",
            "source_parents": ["inside-2", "annulus-2", "outside-again", "stokes", "x2-lem", "x3-lem"],
            "needed_for": "C0 >= 2*K_pp0/c_pp0",
            "quantification_status": "PARENT_CONSTANTS_OPAQUE",
        },
        {
            "constant_symbol": "c_pp0",
            "role": "linear negative gain in the pp0 case",
            "source_parents": ["inside-2", "origin-repulsion", "ankh", "defect-psi-cor"],
            "needed_for": "C0 >= 2*K_pp0/c_pp0",
            "quantification_status": "PARENT_CONSTANTS_OPAQUE",
        },
        {
            "constant_symbol": "K_not_pp0",
            "role": "coefficient of the combined O(||p||/C0) remainder in the not-pp0 case",
            "source_parents": ["inside-2", "annulus-2", "outside-again", "stokes", "x2-lem", "x3-lem"],
            "needed_for": "C0 >= 2*K_not_pp0*C/c_disp",
            "quantification_status": "PARENT_CONSTANTS_OPAQUE",
        },
        {
            "constant_symbol": "c_disp",
            "role": "conversion from dispersion and exterior critical-point mass to ||p||_1",
            "source_parents": ["disp-split", "inside-2", "annulus-2"],
            "needed_for": "C0 >= 2*K_not_pp0*C/c_disp",
            "quantification_status": "PARENT_CONSTANTS_OPAQUE",
        },
        {
            "constant_symbol": "K_pots_upper",
            "role": "upper-bound coefficient in the Psi contradiction for Proposition pots",
            "source_parents": ["inside", "annulus", "out", "ankh", "origin-repulsion"],
            "needed_for": "C/2 > K_pots_upper",
            "quantification_status": "PARENT_CONSTANTS_OPAQUE",
        },
        {
            "constant_symbol": "N_pots(C,eps)",
            "role": "large-n gate after choosing C and eps in Proposition pots",
            "source_parents": ["pots", "ets", "origin-repulsion", "ankh"],
            "needed_for": "n >= N_pots(C,eps)",
            "quantification_status": "PARENT_CONSTANTS_OPAQUE",
        },
        {
            "constant_symbol": "eps0",
            "role": "small epsilon gate for Proposition ets",
            "source_parents": ["inside", "annulus", "out", "defect-psi-cor"],
            "needed_for": "eps <= eps0",
            "quantification_status": "PARENT_CONSTANTS_OPAQUE",
        },
        {
            "constant_symbol": "N_ets(eps)",
            "role": "large-n gate for the quantitative estimates used in Proposition ets",
            "source_parents": ["inside", "annulus", "out", "defect-psi-cor", "x2-lem", "x3-lem"],
            "needed_for": "n >= N_ets(eps)",
            "quantification_status": "PARENT_CONSTANTS_OPAQUE",
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
        "# EHP114 Tao Parent Constant Source Map\n\n"
        f"Experiment: `{EXPERIMENT_ID}`\n\n"
        "## Verdict\n\n"
        f"Status: `{result['status']}`.\n\n"
        "This packet maps the final symbolic constants to exact source lemmas "
        "and proof phrases. It is a source map, not a numeric threshold "
        "checker.\n\n"
        "## Summary\n\n"
        f"- Parent source rows: `{result['parent_source_row_count']}`\n"
        f"- Constant rows: `{result['constant_row_count']}`\n"
        f"- Numeric-ready constants: `{result['numeric_ready_constant_count']}`\n"
        f"- Source SHA status: `{result['tao_source']['source_sha_status']}`\n\n"
        "## Next Action\n\n"
        f"{result['next_action']}\n\n"
        "## Claim Ceiling\n\n"
        f"{result['claim_ceiling']}\n"
    )


def main() -> None:
    proof_root = proof_path_root()
    symbolic_path = proof_root / SYMBOLIC_ID / f"{SYMBOLIC_ID}_RESULTS.json"
    l40a_path = proof_root / L40A_ID / f"{L40A_ID}_RESULTS.json"
    symbolic = read_json(symbolic_path)
    l40a = read_json(l40a_path)
    symbolic_sha = verify_result_sha("symbolic_budget_reduction", symbolic_path)
    l40a_sha = verify_result_sha("final_section_constant_extraction", l40a_path)

    lines, tao_sha, module = source_lines()
    source_sha_status = "PASS_PINNED_SOURCE_SHA_MATCH" if tao_sha == EXPECTED_TAO_SOURCE_SHA else "FAIL_PINNED_SOURCE_SHA_MISMATCH"
    parent_rows = source_parent_rows(lines, module)
    constants = constant_rows()
    numeric_ready = sum(1 for row in constants if row["quantification_status"] == "NUMERIC_READY")

    status = (
        "PARENT_CONSTANT_SOURCE_MAP_READY_NUMERIC_CONSTANTS_PENDING"
        if source_sha_status == "PASS_PINNED_SOURCE_SHA_MATCH" and len(parent_rows) == len(PARENT_LABELS)
        else "PARENT_CONSTANT_SOURCE_MAP_BLOCKED"
    )

    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "source_experiment_ids": [SYMBOLIC_ID, L40A_ID],
        "source_artifacts": [
            {"artifact_id": SYMBOLIC_ID, "path": str(symbolic_path), "status": symbolic["status"], "sha_check": symbolic_sha},
            {"artifact_id": L40A_ID, "path": str(l40a_path), "status": l40a["status"], "sha_check": l40a_sha},
        ],
        "tao_source": {
            "expected_source_sha256": EXPECTED_TAO_SOURCE_SHA,
            "source_sha256": tao_sha,
            "source_sha_status": source_sha_status,
            "source_line_count": len(lines),
        },
        "parent_source_row_count": len(parent_rows),
        "parent_source_rows": parent_rows,
        "constant_row_count": len(constants),
        "constant_rows": constants,
        "numeric_ready_constant_count": numeric_ready,
        "partial_max_N_final_section": None,
        "global_candidate_N0": None,
        "n15_or_higher_computation_authorized": False,
        "full_synthesis_allowed": False,
        "first_failed_condition": "numeric parent constants remain opaque" if numeric_ready == 0 else "none",
        "next_action": "Start quantifying the smallest parent rows first: disp-split, origin-repulsion, ankh, and defect-psi-cor; then propagate into inside-2 and annulus-2.",
        "claim_ceiling": "Parent constant source map only. This packet does not provide a numeric Tao threshold or a full all-degree proof.",
    }

    outdir = proof_root / EXPERIMENT_ID
    digest = write_artifact(outdir, result, report_for(result))
    print(
        json.dumps(
            {
                "status": "EHP114_PARENT_CONSTANT_SOURCE_MAP_EMITTED",
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
