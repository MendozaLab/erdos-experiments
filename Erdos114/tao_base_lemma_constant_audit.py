"""Emit an EHP114 Tao base-lemma constant audit packet.

This is the L42 bridge pass. It audits the four smallest Tao-side parent
sources behind the final-section symbolic budget:

  disp-split, origin-repulsion, ankh, defect-psi-cor

The artifact maps each target constant to exact TeX snippets and named
subdependencies. It does not manufacture numeric constants, extract a global
threshold, or authorize higher-degree computation.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-TAO-BASE-LEMMA-CONSTANT-AUDIT-20260507-01"
PARENT_MAP_ID = "EXP-MATH-EHP114-TAO-PARENT-CONSTANT-SOURCE-MAP-20260507-01"
SYMBOLIC_ID = "EXP-MATH-EHP114-TAO-FINAL-SECTION-SYMBOLIC-BUDGET-REDUCTION-20260507-01"
EXPECTED_TAO_SOURCE_SHA = "ab792bd8fe806a985e44f692ccb8305cf386677509846979bb8ff35a276814a6"

NUMERIC_SEED_READY = "NUMERIC_SEED_READY"
SYMBOLIC_CONSTANT_READY = "SYMBOLIC_CONSTANT_READY"
OPAQUE_SUBDEPENDENCY = "OPAQUE_SUBDEPENDENCY"
SOURCE_LABEL_MISSING = "SOURCE_LABEL_MISSING"

TARGETS = [
    {
        "constant_symbol": "K_disp_split",
        "lemma_label": "disp-split",
        "target_role": "constant converting total critical-point size into exterior mass plus local dispersion",
        "dependencies": ["markov-many", "disp-upper", "semi-norm"],
        "symbolic_seed": "K_disp_split = F(K_markov_many, K_disp_upper, K_semi_norm)",
        "next_subdependency": "Quantify markov-many and disp-upper constants before K_disp_split can become numeric.",
    },
    {
        "constant_symbol": "K_origin",
        "lemma_label": "origin-repulsion",
        "target_role": "constant converting origin displacement into critical-size and distance-to-lemniscate terms",
        "dependencies": ["poz", "npz", "pocl", "size-def"],
        "symbolic_seed": "K_origin = F(K_poz, K_npz, K_pocl, K_size_def)",
        "next_subdependency": "Quantify the exponential constants in poz/pocl and the normalization constant in npz.",
    },
    {
        "constant_symbol": "K_ankh",
        "lemma_label": "ankh",
        "target_role": "constant in the lower bound for Psi on a dyadic annulus",
        "dependencies": ["psi-expand", "Psi-def", "semi-norm", "riesz-bound"],
        "symbolic_seed": "K_ankh = F(K_riesz, K_pointwise_annulus, K_semi_norm)",
        "next_subdependency": "Quantify the Riesz-bound constant and the pointwise annulus estimate.",
    },
    {
        "constant_symbol": "c_defect_C",
        "lemma_label": "defect-psi-cor",
        "target_role": "C-dependent positive defect constant in the upper bound for Psi",
        "dependencies": ["tridef", "markov-many", "disp-sum", "multip"],
        "symbolic_seed": "c_defect_C = F(c_tridef_C, c_markov_many, c_disp_sum, K_multip)",
        "next_subdependency": "Quantify the C-dependent tridef lower constant and markov-many population constant.",
    },
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


def compact(lines: list[str], start: int, end: int, max_chars: int = 1200) -> str:
    text = " ".join(line.strip() for line in lines[start - 1 : end])
    text = re.sub(r"\s+", " ", text)
    if len(text) > max_chars:
        return text[: max_chars - 3] + "..."
    return text


def find_label_span(lines: list[str], module: Any, label: str) -> dict[str, Any]:
    try:
        start, end, env, title = module.find_env(lines, label)
        return {
            "label": label,
            "source_found": True,
            "source_span": [start, end],
            "environment": env,
            "title": title,
            "source_snippet": compact(lines, start, end),
            "named_blocker": None,
        }
    except RuntimeError:
        return {
            "label": label,
            "source_found": False,
            "source_span": None,
            "environment": None,
            "title": None,
            "source_snippet": None,
            "named_blocker": "source_label_missing",
        }


def has_asymptotic_tokens(text: str | None) -> bool:
    if not text:
        return False
    tokens = [r"\\ll", r"\\gg", r"\\asymp", r"\bO\(", r"o\("]
    return any(re.search(token, text) for token in tokens)


def dependency_rows(lines: list[str], module: Any, labels: list[str]) -> list[dict[str, Any]]:
    rows = []
    for label in labels:
        row = find_label_span(lines, module, label)
        if row["source_found"] and has_asymptotic_tokens(row["source_snippet"]):
            row["dependency_status"] = OPAQUE_SUBDEPENDENCY
            row["named_blocker"] = "asymptotic_or_comparison_constant_not_numeric"
        elif row["source_found"]:
            row["dependency_status"] = SYMBOLIC_CONSTANT_READY
        else:
            row["dependency_status"] = SOURCE_LABEL_MISSING
        rows.append(row)
    return rows


def constant_seed_row(lines: list[str], module: Any, target: dict[str, Any]) -> dict[str, Any]:
    lemma = find_label_span(lines, module, target["lemma_label"])
    deps = dependency_rows(lines, module, target["dependencies"])
    if not lemma["source_found"] or any(row["dependency_status"] == SOURCE_LABEL_MISSING for row in deps):
        row_status = SOURCE_LABEL_MISSING
    elif any(row["dependency_status"] == OPAQUE_SUBDEPENDENCY for row in deps) or has_asymptotic_tokens(lemma["source_snippet"]):
        row_status = OPAQUE_SUBDEPENDENCY
    elif deps:
        row_status = SYMBOLIC_CONSTANT_READY
    else:
        row_status = NUMERIC_SEED_READY

    blockers = [row["label"] for row in deps if row["dependency_status"] in {OPAQUE_SUBDEPENDENCY, SOURCE_LABEL_MISSING}]
    if has_asymptotic_tokens(lemma["source_snippet"]):
        blockers.append(target["lemma_label"])
    blockers = sorted(set(blockers))

    return {
        "constant_symbol": target["constant_symbol"],
        "lemma_label": target["lemma_label"],
        "target_role": target["target_role"],
        "constant_status": row_status,
        "symbolic_seed": target["symbolic_seed"],
        "lemma_source": lemma,
        "dependency_rows": deps,
        "opaque_or_missing_subdependencies": blockers,
        "next_subdependency": target["next_subdependency"] if blockers else "none",
    }


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
        "# EHP114 Tao Base-Lemma Constant Audit\n\n"
        f"Experiment: `{EXPERIMENT_ID}`\n\n"
        "## Verdict\n\n"
        f"Status: `{result['status']}`.\n\n"
        "This packet source-maps the first base constants feeding the Tao "
        "threshold bridge. It intentionally does not emit an all-degree "
        "threshold or authorize computation beyond the finite packet range.\n\n"
        "## Summary\n\n"
        f"- Processed base lemmas: `{result['processed_base_lemma_count']}`\n"
        f"- Numeric seed rows: `{result['numeric_seed_ready_count']}`\n"
        f"- Symbolic constant rows: `{result['symbolic_constant_ready_count']}`\n"
        f"- Opaque subdependency rows: `{result['opaque_subdependency_count']}`\n"
        f"- Tao source SHA status: `{result['tao_source_sha_status']}`\n\n"
        "## First Failed Condition\n\n"
        f"`{result['first_failed_condition']}`\n\n"
        "## Claim Ceiling\n\n"
        f"{result['claim_ceiling']}\n"
    )


def main() -> None:
    proof_root = proof_path_root()
    parent_path = proof_root / PARENT_MAP_ID / f"{PARENT_MAP_ID}_RESULTS.json"
    symbolic_path = proof_root / SYMBOLIC_ID / f"{SYMBOLIC_ID}_RESULTS.json"
    parent = read_json(parent_path)
    symbolic = read_json(symbolic_path)
    parent_sha = verify_result_sha("parent_constant_source_map", parent_path)
    symbolic_sha = verify_result_sha("symbolic_budget_reduction", symbolic_path)

    lines, tao_sha, module = source_lines()
    tao_source_sha_status = "PASS_PINNED_SOURCE_SHA_MATCH" if tao_sha == EXPECTED_TAO_SOURCE_SHA else "FAIL_PINNED_SOURCE_SHA_MISMATCH"
    rows = [constant_seed_row(lines, module, target) for target in TARGETS]

    numeric_count = sum(1 for row in rows if row["constant_status"] == NUMERIC_SEED_READY)
    symbolic_count = sum(1 for row in rows if row["constant_status"] == SYMBOLIC_CONSTANT_READY)
    opaque_count = sum(1 for row in rows if row["constant_status"] == OPAQUE_SUBDEPENDENCY)
    missing_count = sum(1 for row in rows if row["constant_status"] == SOURCE_LABEL_MISSING)

    silent_generic_count = 0
    for row in rows:
        for dep in row["dependency_rows"]:
            if dep["source_found"] and dep["dependency_status"] == OPAQUE_SUBDEPENDENCY and not dep.get("named_blocker"):
                silent_generic_count += 1

    if tao_source_sha_status != "PASS_PINNED_SOURCE_SHA_MATCH" or missing_count > 0 or silent_generic_count > 0:
        status = "BASE_LEMMA_CONSTANT_AUDIT_BLOCKED_SOURCE"
    elif numeric_count == len(rows):
        status = "BASE_LEMMA_CONSTANT_AUDIT_READY"
    elif numeric_count > 0:
        status = "BASE_LEMMA_CONSTANT_AUDIT_PARTIAL_NUMERIC"
    else:
        status = "BASE_LEMMA_CONSTANT_AUDIT_REMAINS_SYMBOLIC"

    if missing_count:
        first_failed = "source label missing"
    elif silent_generic_count:
        first_failed = "silent generic dependency row"
    elif opaque_count:
        first_failed = "asymptotic or comparison subconstants remain unquantified"
    else:
        first_failed = "none"

    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "source_experiment_ids": [PARENT_MAP_ID, SYMBOLIC_ID],
        "source_artifacts": [
            {"artifact_id": PARENT_MAP_ID, "path": str(parent_path), "status": parent["status"], "sha_check": parent_sha},
            {"artifact_id": SYMBOLIC_ID, "path": str(symbolic_path), "status": symbolic["status"], "sha_check": symbolic_sha},
        ],
        "tao_source_sha_status": tao_source_sha_status,
        "tao_source_sha256": tao_sha,
        "expected_tao_source_sha256": EXPECTED_TAO_SOURCE_SHA,
        "processed_base_lemma_count": len(rows),
        "constant_seed_rows": rows,
        "numeric_seed_ready_count": numeric_count,
        "symbolic_constant_ready_count": symbolic_count,
        "opaque_subdependency_count": opaque_count,
        "source_label_missing_count": missing_count,
        "silent_generic_count": silent_generic_count,
        "first_failed_condition": first_failed,
        "partial_max_N_final_section": None,
        "global_candidate_N0": None,
        "n15_or_higher_computation_authorized": False,
        "full_synthesis_allowed": False,
        "claim_ceiling": "Base-lemma constant audit only. This packet does not provide a numeric Tao threshold, does not authorize n=15 or higher computation, and does not prove the full all-degree statement.",
    }

    outdir = proof_root / EXPERIMENT_ID
    digest = write_artifact(outdir, result, report_for(result))
    print(
        json.dumps(
            {
                "status": "EHP114_BASE_LEMMA_CONSTANT_AUDIT_EMITTED",
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
