"""Emit an EHP114 Tao pvar elementary inequality packet.

L46 converts the pvar blockers from L45 into elementary numeric or symbolic
inequality gates. It intentionally leaves symbolic parameters such as A_pocl,
B1, and B2 rather than inventing numeric values not present in the source.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-TAO-PVAR-ELEMENTARY-INEQUALITY-PACKET-20260507-01"
L45_ID = "EXP-MATH-EHP114-TAO-PVAR-EXPONENTIAL-CONSTANT-EXTRACTION-20260507-01"
EXPECTED_TAO_SOURCE_SHA = "ab792bd8fe806a985e44f692ccb8305cf386677509846979bb8ff35a276814a6"

NUMERIC_ELEMENTARY_READY = "NUMERIC_ELEMENTARY_READY"
SYMBOLIC_ELEMENTARY_READY = "SYMBOLIC_ELEMENTARY_READY"
OPAQUE_ELEMENTARY = "OPAQUE_ELEMENTARY"
SOURCE_LABEL_MISSING = "SOURCE_LABEL_MISSING"

TARGET_LABELS = ["pvar", "p'-bound", "poz", "ppcl", "pocl", "origin-repulsion"]


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


def elementary_rows() -> list[dict[str, Any]]:
    return [
        {
            "row_id": "centered_product_bound",
            "elementary_status": NUMERIC_ELEMENTARY_READY,
            "source_labels": ["p'-bound"],
            "inequality": "|(1-w)exp(w)| <= exp(|w|^2) for |w| <= 1/2; use exp(2|w|) otherwise",
            "named_parameters": {"C_centered_local": 1.0, "local_radius": 0.5},
            "supports": ["p'-bound centered product estimate"],
            "named_blocker": None,
        },
        {
            "row_id": "pprime_bound_coefficients",
            "elementary_status": SYMBOLIC_ELEMENTARY_READY,
            "source_labels": ["p'-bound"],
            "inequality": "|p'(z)| <= n|z|^(n-1) exp(B1 min(||p||_1/|z|, ||p||_1^2/|z|^2))",
            "named_parameters": {"B1": "max(product_factor_coefficient, centered_product_global_coefficient)"},
            "supports": ["ppcl", "poz radial gates"],
            "named_blocker": "B1 is named but not numerically minimized",
        },
        {
            "row_id": "ppcl_to_pocl_symbolic_gate",
            "elementary_status": SYMBOLIC_ELEMENTARY_READY,
            "source_labels": ["ppcl", "pocl"],
            "inequality": "if |z| <= A_pocl ||p||_1/n and |p'| <= exp(B1 n)n^(-n)||p||_1^(n-1), then |p(z)-p(0)| <= A_pocl exp(B1 n)n^(-n)||p||_1^n",
            "named_parameters": {"A_pocl": "small-radius maximum-principle/integration radius", "B1": "pprime exponent coefficient"},
            "supports": ["pocl"],
            "named_blocker": "A_pocl remains symbolic",
        },
        {
            "row_id": "poz_first_radial_gate",
            "elementary_status": SYMBOLIC_ELEMENTARY_READY,
            "source_labels": ["p'-bound", "poz"],
            "inequality": "for |z| >= 2C||p||_1/n, derivative is negative if C >= max(2, 2B1)",
            "named_parameters": {"C": "max(2,2*B1)", "B1": "pprime exponent coefficient"},
            "supports": ["poz small-to-intermediate radial regime"],
            "named_blocker": "B1 remains symbolic",
        },
        {
            "row_id": "poz_second_radial_gate",
            "elementary_status": SYMBOLIC_ELEMENTARY_READY,
            "source_labels": ["p'-bound", "poz"],
            "inequality": "for |z| >= 3C'||p||_1, derivative is negative if C' >= max(2, B2) and first radial gate is available",
            "named_parameters": {"C_prime": "max(2,B2)", "B2": "large-radius pprime exponent coefficient"},
            "supports": ["poz large-radius radial regime"],
            "named_blocker": "B2 remains symbolic and depends on first radial gate",
        },
        {
            "row_id": "k_origin_after_elementary_packet",
            "elementary_status": SYMBOLIC_ELEMENTARY_READY,
            "source_labels": ["poz", "pocl", "origin-repulsion"],
            "inequality": "K_origin = F(A_pocl, B1, B2, C, C_prime) once pvar symbolic gates are accepted",
            "named_parameters": {"A_pocl": "symbolic", "B1": "symbolic", "B2": "symbolic"},
            "supports": ["K_origin"],
            "named_blocker": "K_origin is symbolic-ready but not numeric",
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
        "# EHP114 Tao pvar Elementary Inequality Packet\n\n"
        f"Experiment: `{EXPERIMENT_ID}`\n\n"
        "## Verdict\n\n"
        f"Status: `{result['status']}`.\n\n"
        "This packet replaces pvar exponential opacity with named elementary "
        "symbolic gates. It emits no global threshold.\n\n"
        "## Summary\n\n"
        f"- Processed elementary rows: `{result['processed_elementary_row_count']}`\n"
        f"- Numeric-ready rows: `{result['numeric_elementary_ready_count']}`\n"
        f"- Symbolic-ready rows: `{result['symbolic_elementary_ready_count']}`\n"
        f"- Opaque rows: `{result['opaque_elementary_count']}`\n"
        f"- K_origin after packet: `{result['k_origin_status_after_elementary_packet']}`\n"
        f"- Tao source SHA status: `{result['tao_source_sha_status']}`\n\n"
        "## First Failed Condition\n\n"
        f"`{result['first_failed_condition']}`\n\n"
        "## Claim Ceiling\n\n"
        f"{result['claim_ceiling']}\n"
    )


def main() -> None:
    proof_root = proof_path_root()
    l45_path = proof_root / L45_ID / f"{L45_ID}_RESULTS.json"
    l45 = read_json(l45_path)
    l45_sha = verify_result_sha("pvar_exponential_constant_extraction", l45_path)

    lines, source_sha, module = source_lines()
    tao_source_sha_status = "PASS_PINNED_SOURCE_SHA_MATCH" if source_sha == EXPECTED_TAO_SOURCE_SHA else "FAIL_PINNED_SOURCE_SHA_MISMATCH"
    target_sources = [source_row(lines, module, label) for label in TARGET_LABELS]
    missing_count = sum(1 for row in target_sources if not row["source_found"])
    rows = elementary_rows()

    numeric_count = sum(1 for row in rows if row["elementary_status"] == NUMERIC_ELEMENTARY_READY)
    symbolic_count = sum(1 for row in rows if row["elementary_status"] == SYMBOLIC_ELEMENTARY_READY)
    opaque_count = sum(1 for row in rows if row["elementary_status"] == OPAQUE_ELEMENTARY)
    silent_generic_count = sum(1 for row in rows if row["elementary_status"] == OPAQUE_ELEMENTARY and not row.get("named_blocker"))

    if tao_source_sha_status != "PASS_PINNED_SOURCE_SHA_MATCH" or missing_count > 0 or silent_generic_count > 0:
        status = "PVAR_ELEMENTARY_INEQUALITY_BLOCKED_SOURCE"
    elif opaque_count == 0:
        status = "PVAR_ELEMENTARY_INEQUALITY_READY"
    elif numeric_count + symbolic_count > 0:
        status = "PVAR_ELEMENTARY_INEQUALITY_PARTIAL"
    else:
        status = "PVAR_ELEMENTARY_INEQUALITY_REMAINS_OPAQUE"

    if missing_count:
        first_failed = "source label missing"
    elif silent_generic_count:
        first_failed = "silent generic elementary row"
    elif opaque_count:
        first_failed = "opaque elementary inequality remains"
    elif symbolic_count:
        first_failed = "symbolic parameters remain to be numerically selected"
    else:
        first_failed = "none"

    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "source_experiment_ids": [L45_ID],
        "source_artifacts": [
            {"artifact_id": L45_ID, "path": str(l45_path), "status": l45["status"], "sha_check": l45_sha},
        ],
        "tao_source_sha_status": tao_source_sha_status,
        "tao_source_sha256": source_sha,
        "expected_tao_source_sha256": EXPECTED_TAO_SOURCE_SHA,
        "target_source_rows": target_sources,
        "target_source_missing_count": missing_count,
        "processed_elementary_row_count": len(rows),
        "elementary_inequality_rows": rows,
        "numeric_elementary_ready_count": numeric_count,
        "symbolic_elementary_ready_count": symbolic_count,
        "opaque_elementary_count": opaque_count,
        "silent_generic_count": silent_generic_count,
        "k_origin_status_after_elementary_packet": "K_ORIGIN_SYMBOLIC_READY_NUMERIC_PARAMETERS_PENDING" if opaque_count == 0 else "K_ORIGIN_OPAQUE",
        "first_failed_condition": first_failed,
        "partial_max_N_final_section": None,
        "global_candidate_N0": None,
        "n15_or_higher_computation_authorized": False,
        "full_synthesis_allowed": False,
        "claim_ceiling": "pvar elementary inequality packet only. This packet does not provide a numeric Tao threshold, does not authorize n=15 or higher computation, and does not prove the full all-degree statement.",
    }

    outdir = proof_root / EXPERIMENT_ID
    digest = write_artifact(outdir, result, report_for(result))
    print(
        json.dumps(
            {
                "status": "EHP114_PVAR_ELEMENTARY_INEQUALITY_PACKET_EMITTED",
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
