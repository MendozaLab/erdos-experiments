"""Emit an EHP114 Tao pvar exponential-constant extraction packet.

L45 attacks the current K_origin blocker from L44. It source-maps the pvar
estimates behind poz and pocl and separates direct numeric elementary bounds
from still-unquantified exponential constants.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-TAO-PVAR-EXPONENTIAL-CONSTANT-EXTRACTION-20260507-01"
L44_ID = "EXP-MATH-EHP114-TAO-BASE-LEMMA-SEED-PROPAGATION-20260507-02"
L43_ID = "EXP-MATH-EHP114-TAO-PRIMITIVE-CONSTANT-SEED-AUDIT-20260507-01"
EXPECTED_TAO_SOURCE_SHA = "ab792bd8fe806a985e44f692ccb8305cf386677509846979bb8ff35a276814a6"

NUMERIC_PVAR_READY = "NUMERIC_PVAR_READY"
SYMBOLIC_PVAR_READY = "SYMBOLIC_PVAR_READY"
OPAQUE_PVAR = "OPAQUE_PVAR"
SOURCE_LABEL_MISSING = "SOURCE_LABEL_MISSING"

TARGET_LABELS = ["pvar", "p'-bound", "p'-bound-lower", "poz", "ppcl", "pocl", "origin-repulsion"]


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


def compact(lines: list[str], start: int, end: int, max_chars: int = 1400) -> str:
    text = " ".join(line.strip() for line in lines[start - 1 : end])
    text = re.sub(r"\s+", " ", text)
    if len(text) > max_chars:
        return text[: max_chars - 3] + "..."
    return text


def source_row(lines: list[str], module: Any, label: str) -> dict[str, Any]:
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


def pvar_constant_rows() -> list[dict[str, Any]]:
    return [
        {
            "row_id": "product_factor_upper_bound",
            "pvar_status": NUMERIC_PVAR_READY,
            "source_labels": ["pp-factor", "p'-bound"],
            "constant_expression": "|1-w| <= exp(|w|)",
            "numeric_value": {"exponent_coefficient": 1.0},
            "supports": ["p'-bound first upper estimate"],
            "named_blocker": None,
        },
        {
            "row_id": "centered_product_upper_bound",
            "pvar_status": SYMBOLIC_PVAR_READY,
            "source_labels": ["p'-bound"],
            "constant_expression": "|(1-w)exp(w)| <= exp(C_centered |w|^2)",
            "numeric_value": None,
            "supports": ["p'-bound centered upper estimate"],
            "named_blocker": "elementary centered-product constant C_centered not numerically selected",
        },
        {
            "row_id": "ppcl_to_pocl_integration",
            "pvar_status": OPAQUE_PVAR,
            "source_labels": ["ppcl", "pocl"],
            "constant_expression": "integrate |p'(z)| <= exp(C_ppcl n)n^(-n)||p||_1^(n-1) over radius O(||p||_1/n)",
            "numeric_value": None,
            "supports": ["pocl"],
            "named_blocker": "ppcl exponential constant and small-radius integration radius are not quantified",
        },
        {
            "row_id": "poz_radial_derivative_first_gate",
            "pvar_status": OPAQUE_PVAR,
            "source_labels": ["p'-bound", "poz"],
            "constant_expression": "choose C so C n - C^2 ||p||_1/|z| dominates p'-bound for |z| >= 2C||p||_1/n",
            "numeric_value": None,
            "supports": ["poz small-to-intermediate radial regime"],
            "named_blocker": "radial derivative domination constant C is large-enough but not explicit",
        },
        {
            "row_id": "poz_radial_derivative_second_gate",
            "pvar_status": OPAQUE_PVAR,
            "source_labels": ["p'-bound", "poz"],
            "constant_expression": "choose C' so C' n - 2(C')^2 ||p||_1/|z| dominates p'-bound for |z| >= 3C'||p||_1",
            "numeric_value": None,
            "supports": ["poz large-radius radial regime"],
            "named_blocker": "radial derivative domination constant C_prime is large-enough but not explicit",
        },
        {
            "row_id": "k_origin_after_pvar",
            "pvar_status": OPAQUE_PVAR,
            "source_labels": ["poz", "pocl", "origin-repulsion"],
            "constant_expression": "K_origin awaits explicit C_ppcl, C_pocl, C_poz, and C_poz_prime constants",
            "numeric_value": None,
            "supports": ["K_origin"],
            "named_blocker": "K_origin remains blocked by pvar exponential constants",
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
        "# EHP114 Tao pvar Exponential Constant Extraction\n\n"
        f"Experiment: `{EXPERIMENT_ID}`\n\n"
        "## Verdict\n\n"
        f"Status: `{result['status']}`.\n\n"
        "This packet audits the pvar constants behind poz, pocl, and K_origin. "
        "It names the elementary blockers without emitting a global threshold.\n\n"
        "## Summary\n\n"
        f"- Processed pvar rows: `{result['processed_pvar_row_count']}`\n"
        f"- Numeric-ready pvar rows: `{result['numeric_pvar_ready_count']}`\n"
        f"- Symbolic-ready pvar rows: `{result['symbolic_pvar_ready_count']}`\n"
        f"- Opaque pvar rows: `{result['opaque_pvar_count']}`\n"
        f"- K_origin after pvar: `{result['k_origin_status_after_pvar']}`\n"
        f"- Tao source SHA status: `{result['tao_source_sha_status']}`\n\n"
        "## First Failed Condition\n\n"
        f"`{result['first_failed_condition']}`\n\n"
        "## Claim Ceiling\n\n"
        f"{result['claim_ceiling']}\n"
    )


def main() -> None:
    proof_root = proof_path_root()
    l44_path = proof_root / L44_ID / f"{L44_ID}_RESULTS.json"
    l43_path = proof_root / L43_ID / f"{L43_ID}_RESULTS.json"
    l44 = read_json(l44_path)
    l43 = read_json(l43_path)
    l44_sha = verify_result_sha("base_lemma_seed_propagation", l44_path)
    l43_sha = verify_result_sha("primitive_constant_seed_audit", l43_path)

    lines, source_sha, module = source_lines()
    tao_source_sha_status = "PASS_PINNED_SOURCE_SHA_MATCH" if source_sha == EXPECTED_TAO_SOURCE_SHA else "FAIL_PINNED_SOURCE_SHA_MISMATCH"
    target_sources = [source_row(lines, module, label) for label in TARGET_LABELS]
    rows = pvar_constant_rows()

    missing_count = sum(1 for row in target_sources if not row["source_found"])
    numeric_count = sum(1 for row in rows if row["pvar_status"] == NUMERIC_PVAR_READY)
    symbolic_count = sum(1 for row in rows if row["pvar_status"] == SYMBOLIC_PVAR_READY)
    opaque_count = sum(1 for row in rows if row["pvar_status"] == OPAQUE_PVAR)
    silent_generic_count = sum(1 for row in rows if row["pvar_status"] == OPAQUE_PVAR and not row.get("named_blocker"))

    if tao_source_sha_status != "PASS_PINNED_SOURCE_SHA_MATCH" or missing_count > 0 or silent_generic_count > 0:
        status = "PVAR_CONSTANT_EXTRACTION_BLOCKED_SOURCE"
    elif opaque_count == 0:
        status = "PVAR_CONSTANT_EXTRACTION_READY"
    elif numeric_count + symbolic_count > 0:
        status = "PVAR_CONSTANT_EXTRACTION_PARTIAL"
    else:
        status = "PVAR_CONSTANT_EXTRACTION_REMAINS_OPAQUE"

    if missing_count:
        first_failed = "source label missing"
    elif silent_generic_count:
        first_failed = "silent generic pvar row"
    elif opaque_count:
        first_failed = "pvar exponential and radial-derivative constants remain unquantified"
    else:
        first_failed = "none"

    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "source_experiment_ids": [L44_ID, L43_ID],
        "source_artifacts": [
            {"artifact_id": L44_ID, "path": str(l44_path), "status": l44["status"], "sha_check": l44_sha},
            {"artifact_id": L43_ID, "path": str(l43_path), "status": l43["status"], "sha_check": l43_sha},
        ],
        "tao_source_sha_status": tao_source_sha_status,
        "tao_source_sha256": source_sha,
        "expected_tao_source_sha256": EXPECTED_TAO_SOURCE_SHA,
        "target_source_rows": target_sources,
        "target_source_missing_count": missing_count,
        "processed_pvar_row_count": len(rows),
        "pvar_constant_rows": rows,
        "numeric_pvar_ready_count": numeric_count,
        "symbolic_pvar_ready_count": symbolic_count,
        "opaque_pvar_count": opaque_count,
        "silent_generic_count": silent_generic_count,
        "k_origin_status_after_pvar": "OPAQUE_K_ORIGIN_PVAR_CONSTANTS_PENDING" if opaque_count else "K_ORIGIN_SYMBOLIC_READY",
        "first_failed_condition": first_failed,
        "partial_max_N_final_section": None,
        "global_candidate_N0": None,
        "n15_or_higher_computation_authorized": False,
        "full_synthesis_allowed": False,
        "claim_ceiling": "pvar exponential constant extraction only. This packet does not provide a numeric Tao threshold, does not authorize n=15 or higher computation, and does not prove the full all-degree statement.",
    }

    outdir = proof_root / EXPERIMENT_ID
    digest = write_artifact(outdir, result, report_for(result))
    print(
        json.dumps(
            {
                "status": "EHP114_PVAR_CONSTANT_EXTRACTION_EMITTED",
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
