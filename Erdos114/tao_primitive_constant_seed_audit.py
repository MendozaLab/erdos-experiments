"""Emit an EHP114 Tao primitive constant seed audit packet.

L43 takes the opaque subdependencies identified by L42 and separates them into
usable primitive seeds, symbolic seeds, and genuinely opaque asymptotic rows.
It does not compute a global Tao threshold and does not authorize computation
for n=15 or higher.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-TAO-PRIMITIVE-CONSTANT-SEED-AUDIT-20260507-01"
L42_ID = "EXP-MATH-EHP114-TAO-BASE-LEMMA-CONSTANT-AUDIT-20260507-01"
EXPECTED_TAO_SOURCE_SHA = "ab792bd8fe806a985e44f692ccb8305cf386677509846979bb8ff35a276814a6"

NUMERIC_SEED_READY = "NUMERIC_SEED_READY"
SYMBOLIC_SEED_READY = "SYMBOLIC_SEED_READY"
OPAQUE_SEED = "OPAQUE_SEED"
SOURCE_LABEL_MISSING = "SOURCE_LABEL_MISSING"

PRIMITIVE_ROWS = [
    {
        "label": "markov-many",
        "classification": NUMERIC_SEED_READY,
        "seed_expression": "# {zeta in D(0,2||p||_1/n)} >= n/4 for n >= 4",
        "numeric_value": {"population_fraction_lower_bound": 0.25, "valid_for_n_min": 4},
        "named_blocker": None,
    },
    {
        "label": "disp-upper",
        "classification": SYMBOLIC_SEED_READY,
        "seed_expression": "triangle-inequality dispersion comparison with unit coefficients",
        "numeric_value": None,
        "named_blocker": None,
    },
    {
        "label": "semi-norm",
        "classification": SYMBOLIC_SEED_READY,
        "seed_expression": "definition ||p||_1 = sum_zeta |zeta|",
        "numeric_value": None,
        "named_blocker": None,
    },
    {
        "label": "size-def",
        "classification": SYMBOLIC_SEED_READY,
        "seed_expression": "definition ||p|| = ||p||_1 + ||p||_0",
        "numeric_value": None,
        "named_blocker": None,
    },
    {
        "label": "npz",
        "classification": NUMERIC_SEED_READY,
        "seed_expression": "|p(z)-p(0)| >= n^(-n)||p||_0^n on the unit-boundary comparison point",
        "numeric_value": {"coefficient": 1.0, "scale": "n^(-n)||p||_0^n"},
        "named_blocker": None,
    },
    {
        "label": "disp-sum",
        "classification": NUMERIC_SEED_READY,
        "seed_expression": "Disp <= (2/#A) sum_{xi in A} sum_i |z_i-xi|",
        "numeric_value": {"coefficient": 2.0},
        "named_blocker": None,
    },
    {
        "label": "multip",
        "classification": NUMERIC_SEED_READY,
        "seed_expression": "integral_E sum_j 1/|z-zeta_j| <= 2*pi*m*r[E]",
        "numeric_value": {"coefficient": "2*pi"},
        "named_blocker": None,
    },
    {
        "label": "riesz-bound",
        "classification": NUMERIC_SEED_READY,
        "seed_expression": "I_1 1_E(z0) <= 2*pi*r[E]",
        "numeric_value": {"coefficient": "2*pi"},
        "named_blocker": None,
    },
    {
        "label": "Psi-def",
        "classification": NUMERIC_SEED_READY,
        "seed_expression": "Psi(E) = (1/pi) integral_E |psi| dA",
        "numeric_value": {"coefficient": "1/pi"},
        "named_blocker": None,
    },
    {
        "label": "tridef",
        "classification": SYMBOLIC_SEED_READY,
        "seed_expression": "c_tridef(C) is the positive integral of the normalized defect on D(0,1/(2C+1))",
        "numeric_value": None,
        "named_blocker": "C-dependent integral lower bound not numerically evaluated",
    },
    {
        "label": "poz",
        "classification": OPAQUE_SEED,
        "seed_expression": "p(z)=p(0)+O(|z|^n exp(O(min(||p||_1/|z|,||p||_1^2/|z|^2))))",
        "numeric_value": None,
        "named_blocker": "exponential O-constant in pvar is not quantified",
    },
    {
        "label": "pocl",
        "classification": OPAQUE_SEED,
        "seed_expression": "p(z)=p(0)+O(exp(O(n)) n^(-n)||p||_1^n)",
        "numeric_value": None,
        "named_blocker": "exponential O-constant in pvar is not quantified",
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


def find_label_source(lines: list[str], module: Any, label: str) -> dict[str, Any]:
    try:
        start, end, env, title = module.find_env(lines, label)
        return {
            "source_found": True,
            "source_span": [start, end],
            "environment": env,
            "title": title,
            "source_snippet": compact(lines, start, end),
            "source_label_blocker": None,
        }
    except RuntimeError:
        return {
            "source_found": False,
            "source_span": None,
            "environment": None,
            "title": None,
            "source_snippet": None,
            "source_label_blocker": "source_label_missing",
        }


def primitive_seed_rows(lines: list[str], module: Any) -> list[dict[str, Any]]:
    rows = []
    for spec in PRIMITIVE_ROWS:
        source = find_label_source(lines, module, spec["label"])
        classification = spec["classification"] if source["source_found"] else SOURCE_LABEL_MISSING
        named_blocker = spec["named_blocker"] if source["source_found"] else "source_label_missing"
        rows.append(
            {
                "label": spec["label"],
                "seed_status": classification,
                "seed_expression": spec["seed_expression"],
                "numeric_value": spec["numeric_value"],
                "named_blocker": named_blocker,
                "source": source,
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
        "# EHP114 Tao Primitive Constant Seed Audit\n\n"
        f"Experiment: `{EXPERIMENT_ID}`\n\n"
        "## Verdict\n\n"
        f"Status: `{result['status']}`.\n\n"
        "This packet separates primitive Tao-side constants into numeric seeds, "
        "symbolic seeds, and opaque asymptotic rows. It does not emit a global "
        "threshold.\n\n"
        "## Summary\n\n"
        f"- Processed primitives: `{result['processed_primitive_count']}`\n"
        f"- Numeric seed rows: `{result['numeric_seed_ready_count']}`\n"
        f"- Symbolic seed rows: `{result['symbolic_seed_ready_count']}`\n"
        f"- Opaque seed rows: `{result['opaque_seed_count']}`\n"
        f"- Source-label missing rows: `{result['source_label_missing_count']}`\n"
        f"- Tao source SHA status: `{result['tao_source_sha_status']}`\n\n"
        "## First Failed Condition\n\n"
        f"`{result['first_failed_condition']}`\n\n"
        "## Claim Ceiling\n\n"
        f"{result['claim_ceiling']}\n"
    )


def main() -> None:
    proof_root = proof_path_root()
    l42_path = proof_root / L42_ID / f"{L42_ID}_RESULTS.json"
    l42 = read_json(l42_path)
    l42_sha = verify_result_sha("base_lemma_constant_audit", l42_path)

    lines, tao_sha, module = source_lines()
    tao_source_sha_status = "PASS_PINNED_SOURCE_SHA_MATCH" if tao_sha == EXPECTED_TAO_SOURCE_SHA else "FAIL_PINNED_SOURCE_SHA_MISMATCH"
    rows = primitive_seed_rows(lines, module)

    numeric_count = sum(1 for row in rows if row["seed_status"] == NUMERIC_SEED_READY)
    symbolic_count = sum(1 for row in rows if row["seed_status"] == SYMBOLIC_SEED_READY)
    opaque_count = sum(1 for row in rows if row["seed_status"] == OPAQUE_SEED)
    missing_count = sum(1 for row in rows if row["seed_status"] == SOURCE_LABEL_MISSING)
    silent_generic_count = sum(1 for row in rows if row["seed_status"] == OPAQUE_SEED and not row.get("named_blocker"))

    if tao_source_sha_status != "PASS_PINNED_SOURCE_SHA_MATCH" or missing_count > 0 or silent_generic_count > 0:
        status = "PRIMITIVE_CONSTANT_SEED_AUDIT_BLOCKED_SOURCE"
    elif opaque_count == 0:
        status = "PRIMITIVE_CONSTANT_SEED_AUDIT_READY"
    elif numeric_count + symbolic_count > 0:
        status = "PRIMITIVE_CONSTANT_SEED_AUDIT_PARTIAL"
    else:
        status = "PRIMITIVE_CONSTANT_SEED_AUDIT_REMAINS_OPAQUE"

    if missing_count:
        first_failed = "source label missing"
    elif silent_generic_count:
        first_failed = "silent generic primitive row"
    elif opaque_count:
        first_failed = "poz and pocl exponential constants remain opaque"
    else:
        first_failed = "none"

    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "source_experiment_ids": [L42_ID],
        "source_artifacts": [
            {"artifact_id": L42_ID, "path": str(l42_path), "status": l42["status"], "sha_check": l42_sha},
        ],
        "tao_source_sha_status": tao_source_sha_status,
        "tao_source_sha256": tao_sha,
        "expected_tao_source_sha256": EXPECTED_TAO_SOURCE_SHA,
        "processed_primitive_count": len(rows),
        "primitive_seed_rows": rows,
        "numeric_seed_ready_count": numeric_count,
        "symbolic_seed_ready_count": symbolic_count,
        "opaque_seed_count": opaque_count,
        "source_label_missing_count": missing_count,
        "silent_generic_count": silent_generic_count,
        "first_failed_condition": first_failed,
        "partial_max_N_final_section": None,
        "global_candidate_N0": None,
        "n15_or_higher_computation_authorized": False,
        "full_synthesis_allowed": False,
        "claim_ceiling": "Primitive constant seed audit only. This packet does not provide a numeric Tao threshold, does not authorize n=15 or higher computation, and does not prove the full all-degree statement.",
    }

    outdir = proof_root / EXPERIMENT_ID
    digest = write_artifact(outdir, result, report_for(result))
    print(
        json.dumps(
            {
                "status": "EHP114_PRIMITIVE_CONSTANT_SEED_AUDIT_EMITTED",
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
