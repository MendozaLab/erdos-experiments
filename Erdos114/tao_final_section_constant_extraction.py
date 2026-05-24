"""Emit EHP114 final-section constant-extraction bridge packets.

This is a narrow proof-engineering pass over the final section of Tao's
high-degree EHP proof. It fetches the pinned arXiv source, extracts exact TeX
snippets for the selected final-section dependencies, and emits three local
experiment-contract artifacts:

  L40A final-section constant extraction
  L40B final-section threshold subcheck
  L40C bridge decision refresh

It does not extract a global threshold, authorize n=15 computation, or emit a
full proof synthesis.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import tarfile
import urllib.request
from pathlib import Path
from typing import Any


L40A_ID = "EXP-MATH-EHP114-TAO-FINAL-SECTION-CONSTANT-EXTRACTION-20260507-01"
L40B_ID = "EXP-MATH-EHP114-TAO-FINAL-SECTION-THRESHOLD-SUBCHECK-20260507-01"
L40C_ID = "EXP-MATH-EHP114-BRIDGE-DECISION-REFRESH-FINAL-SECTION-20260507-01"

L37_ID = "EXP-MATH-EHP114-TAO-THRESHOLD-EXTRACTION-SKELETON-20260507-02"
L38_ID = "EXP-MATH-EHP114-TAO-EFFECTIVE-THRESHOLD-CHECKER-20260507-02"
FINITE_ID = "EXP-MATH-EHP114-FINITE-N-LESS-15-PROOF-PACKET-20260507-01"

ARXIV_ID = "2512.12455"
ARXIV_EPRINT_URL = f"https://arxiv.org/e-print/{ARXIV_ID}"
ARXIV_ABS_URL = f"https://arxiv.org/abs/{ARXIV_ID}"
EXPECTED_TAO_SOURCE_SHA = "ab792bd8fe806a985e44f692ccb8305cf386677509846979bb8ff35a276814a6"

TARGET_LABELS = {"ets", "pots", "inside_2", "annulus_2"}
TARGET_FINAL_LABELS = {"inside-2", "annulus-2", "outside-again"}
EXPLICIT_READY = "EXPLICIT_INEQUALITY_READY"
NAMED_BLOCKER = "NAMED_BLOCKER"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def proof_path_root() -> Path:
    return Path(__file__).resolve().parent / "proof_path"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha_path_for(path: Path) -> Path:
    name = path.name
    if not name.endswith("_RESULTS.json"):
        raise ValueError(f"result JSON must end in _RESULTS.json: {path}")
    return path.with_name(name.replace("_RESULTS.json", "_RESULTS.sha256"))


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


def fetch_arxiv_source() -> bytes:
    req = urllib.request.Request(
        ARXIV_EPRINT_URL,
        headers={"User-Agent": "MendozaLab-EHP114-final-section-audit/1.0"},
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        return response.read()


def extract_lemniscate_tex(payload: bytes) -> str:
    data = gzip.decompress(payload)
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:") as archive:
        member = archive.extractfile("lemniscate.tex")
        if member is None:
            raise RuntimeError("lemniscate.tex not found in arXiv source archive")
        return member.read().decode("utf-8", errors="replace")


def find_env(lines: list[str], label: str) -> tuple[int, int, str, str]:
    label_pat = re.compile(r"\\label\{" + re.escape(label) + r"\}")
    label_line = None
    for idx, line in enumerate(lines):
        if label_pat.search(line):
            label_line = idx
            break
    if label_line is None:
        raise RuntimeError(f"label not found: {label}")

    start = label_line
    while start >= 0 and "\\begin{" not in lines[start]:
        start -= 1
    if start < 0:
        start = label_line

    begin_match = re.search(r"\\begin\{([^}]+)\}(?:\[([^\]]*)\])?", lines[start])
    env = begin_match.group(1) if begin_match else "unknown"
    title = begin_match.group(2) if begin_match and begin_match.group(2) else ""

    end = label_line
    end_token = f"\\end{{{env}}}"
    while end < len(lines) and end_token not in lines[end]:
        end += 1
    if end >= len(lines):
        end = label_line
    return start + 1, end + 1, env, title


def compact_snippet(lines: list[str], start: int, end: int, max_chars: int = 900) -> str:
    text = " ".join(line.strip() for line in lines[start - 1 : end])
    text = re.sub(r"\s+", " ", text)
    if len(text) > max_chars:
        return text[: max_chars - 3] + "..."
    return text


def window_for_line(line: int, total: int, radius: int = 2) -> tuple[int, int]:
    return max(1, line - radius), min(total, line + radius)


def infer_named_constants(snippet: str) -> list[str]:
    candidates = [
        ("absolute constant", r"absolute constant"),
        ("large-n lower bound", r"large enough|sufficiently large"),
        ("small parameter choice", r"small enough"),
        ("C", r"\bC\b|\\C\b"),
        ("c", r"\bc\b"),
        ("K", r"\bK\b"),
        ("mu", r"\\mu|\bmu\b"),
        ("epsilon", r"\\varepsilon|\\epsilon|\beps\b"),
        ("n", r"\bn\b"),
        ("dispersion", r"dispersion|\\disp"),
    ]
    found = []
    for name, pattern in candidates:
        if re.search(pattern, snippet):
            found.append(name)
    return found


def blocker_for(dep: dict[str, Any], snippet_text: str) -> str:
    phrase = dep.get("phrase")
    if phrase == "absolute constant":
        return "absolute_constant_not_numerically_named"
    if phrase in {"sufficiently large", "large enough"}:
        return "large_n_gate_not_converted_to_integer_bound"
    if phrase == "small enough":
        return "small_parameter_gate_not_quantified"
    if dep.get("source_label") in TARGET_FINAL_LABELS:
        return "final_deficit_constant_not_explicitly_budgeted"
    if "o(" in snippet_text or "\\to" in snippet_text:
        return "asymptotic_error_not_replaced_by_explicit_bound"
    return "constant_dependency_not_explicit"


def inequality_for(dep: dict[str, Any]) -> str:
    label = dep.get("source_label") or dep.get("classified_label")
    if label == "inside-2" or label == "inside_2":
        return "Make the final inner deficit quantitatively dominate its inherited error terms."
    if label == "annulus-2" or label == "annulus_2":
        return "Make the final intermediate-region deficit quantitatively dominate bulk and shape errors."
    if label == "outside-again":
        return "Make the final outer-region comparison quantitatively absorb endpoint and remainder errors."
    if label == "pots":
        return "Convert total-size collapse to an explicit n-threshold and error budget."
    if label == "ets":
        return "Convert critical-point or dispersion collapse to an explicit n-threshold and error budget."
    return "Convert the local asymptotic dependency into an explicit inequality."


def select_target_dependencies(l37: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for dep in l37["dependency_rows"]:
        if dep.get("classified_label") in TARGET_LABELS or dep.get("source_label") in TARGET_FINAL_LABELS:
            rows.append(dep)
    rows.sort(key=lambda row: (row.get("source_line") or 10_000, row.get("dependency_id", "")))
    return rows


def enrich_dependency(dep: dict[str, Any], lines: list[str]) -> dict[str, Any]:
    if dep.get("source_line") is not None:
        source_line = int(dep["source_line"])
        start, end = window_for_line(source_line, len(lines))
        exact_line = lines[source_line - 1].strip()
        snippet_text = compact_snippet(lines, start, end)
        source_span = [start, end]
    else:
        label = dep["source_label"]
        start, end, env, title = find_env(lines, label)
        exact_line = f"\\label{{{label}}}"
        snippet_text = compact_snippet(lines, start, end)
        source_span = [start, end]
        dep = dict(dep)
        dep["source_env"] = env
        dep["source_title"] = title

    constants = infer_named_constants(snippet_text)
    status = NAMED_BLOCKER
    blocker = blocker_for(dep, snippet_text)
    parents = []
    label = dep.get("source_label") or dep.get("classified_label")
    if label in {"inside-2", "annulus-2", "outside-again", "inside_2", "annulus_2"}:
        parents = ["pots", "ets"]
    elif label == "pots":
        parents = ["inside", "annulus", "outside", "geomcontrol"]
    elif label == "ets":
        parents = ["geomcontrol", "x2-lem", "x3-lem"]

    return {
        "dependency_id": dep["dependency_id"],
        "source_line": dep.get("source_line"),
        "source_label": dep.get("source_label"),
        "classified_label": dep.get("classified_label"),
        "priority": dep.get("priority"),
        "source_span": source_span,
        "exact_proof_phrase": exact_line,
        "source_snippet": snippet_text,
        "inequality_needed": inequality_for(dep),
        "named_constants": constants,
        "dependency_parents": parents,
        "extraction_status": status,
        "named_blocker": blocker,
        "extracted_N_i": None,
    }


def write_artifact(outdir: Path, experiment_id: str, result: dict[str, Any], report: str) -> str:
    outdir.mkdir(parents=True, exist_ok=False)
    result_path = outdir / f"{experiment_id}_RESULTS.json"
    report_path = outdir / f"{experiment_id}_REPORT.md"
    sha_path = outdir / f"{experiment_id}_RESULTS.sha256"
    result_json = json.dumps(result, indent=2, sort_keys=True)
    result_path.write_text(result_json + "\n")
    report_path.write_text(report)
    digest = hashlib.sha256((result_json + "\n").encode("utf-8")).hexdigest()
    sha_path.write_text(f"{digest}  {result_path.name}\n")
    return digest


def l40a_report(result: dict[str, Any]) -> str:
    return (
        f"# EHP114 L40A Tao Final-Section Constant Extraction\n\n"
        f"Experiment: `{L40A_ID}`\n\n"
        "## Verdict\n\n"
        f"Status: `{result['status']}`.\n\n"
        "This packet extracts exact source snippets for the first final-section "
        "threshold dependencies. Every targeted row now has either an explicit "
        "inequality-ready status or a named blocker. No numerical threshold is "
        "claimed here.\n\n"
        "## Summary\n\n"
        f"- Target dependency count: `{result['target_dependency_count']}`\n"
        f"- Explicit rows: `{result['explicit_inequality_ready_count']}`\n"
        f"- Named-blocker rows: `{result['named_blocker_count']}`\n"
        f"- Source SHA status: `{result['arxiv_source']['source_sha_status']}`\n\n"
        "## Claim Ceiling\n\n"
        f"{result['claim_ceiling']}\n"
    )


def l40b_report(result: dict[str, Any]) -> str:
    return (
        f"# EHP114 L40B Final-Section Threshold Subcheck\n\n"
        f"Experiment: `{L40B_ID}`\n\n"
        "## Verdict\n\n"
        f"Status: `{result['status']}`.\n\n"
        "The subcheck only computes row-level thresholds for dependencies whose "
        "constants are already explicit. Since the selected final-section rows "
        "still have named blockers, no partial maximum threshold is emitted.\n\n"
        "## Claim Ceiling\n\n"
        f"{result['claim_ceiling']}\n"
    )


def l40c_report(result: dict[str, Any]) -> str:
    return (
        f"# EHP114 L40C Bridge Decision Refresh\n\n"
        f"Experiment: `{L40C_ID}`\n\n"
        "## Verdict\n\n"
        f"Status: `{result['status']}`.\n\n"
        "The finite side still passes, but the final-section subcheck remains "
        "opaque. Therefore this packet does not authorize higher-degree "
        "computation or full synthesis.\n\n"
        "## Next Action\n\n"
        f"{result['next_action']}\n\n"
        "## Claim Ceiling\n\n"
        f"{result['claim_ceiling']}\n"
    )


def main() -> None:
    root = repo_root()
    proof_root = proof_path_root()
    l37_path = proof_root / L37_ID / f"{L37_ID}_RESULTS.json"
    l38_path = proof_root / L38_ID / f"{L38_ID}_RESULTS.json"
    finite_path = proof_root / FINITE_ID / f"{FINITE_ID}_RESULTS.json"

    l37 = read_json(l37_path)
    l38 = read_json(l38_path)
    finite = read_json(finite_path)
    l37_sha = verify_result_sha("threshold_extraction_skeleton", l37_path)
    l38_sha = verify_result_sha("effective_threshold_checker", l38_path)
    finite_sha = verify_result_sha("finite_n_less15_packet", finite_path)

    tex = extract_lemniscate_tex(fetch_arxiv_source())
    tex_lines = tex.splitlines()
    source_sha = hashlib.sha256(tex.encode("utf-8")).hexdigest()
    source_sha_status = "PASS_PINNED_SOURCE_SHA_MATCH" if source_sha == EXPECTED_TAO_SOURCE_SHA else "FAIL_PINNED_SOURCE_SHA_MISMATCH"

    target_deps = select_target_dependencies(l37)
    extraction_rows = [enrich_dependency(dep, tex_lines) for dep in target_deps]
    explicit_count = sum(row["extraction_status"] == EXPLICIT_READY for row in extraction_rows)
    blocker_count = sum(row["extraction_status"] == NAMED_BLOCKER for row in extraction_rows)
    silent_generic_count = sum(1 for row in extraction_rows if not row.get("named_blocker") and row["extraction_status"] != EXPLICIT_READY)

    l40a_status = (
        "FINAL_SECTION_CONSTANT_EXTRACTION_READY"
        if source_sha_status == "PASS_PINNED_SOURCE_SHA_MATCH" and silent_generic_count == 0
        else "FINAL_SECTION_CONSTANT_EXTRACTION_BLOCKED"
    )
    l40a = {
        "experiment_id": L40A_ID,
        "status": l40a_status,
        "source_artifacts": [
            {"artifact_id": L37_ID, "path": str(l37_path), "status": l37["status"], "sha_check": l37_sha},
            {"artifact_id": L38_ID, "path": str(l38_path), "status": l38["status"], "sha_check": l38_sha},
            {"artifact_id": FINITE_ID, "path": str(finite_path), "status": finite["status"], "sha_check": finite_sha},
        ],
        "arxiv_source": {
            "url": ARXIV_ABS_URL,
            "source_sha256": source_sha,
            "expected_source_sha256": EXPECTED_TAO_SOURCE_SHA,
            "source_sha_status": source_sha_status,
            "source_line_count": len(tex_lines),
        },
        "target_dependency_count": len(extraction_rows),
        "explicit_inequality_ready_count": explicit_count,
        "named_blocker_count": blocker_count,
        "silent_generic_count": silent_generic_count,
        "targeted_labels": ["ets", "pots", "inside-2", "annulus-2", "outside-again"],
        "extraction_rows": extraction_rows,
        "first_failed_condition": "none" if l40a_status == "FINAL_SECTION_CONSTANT_EXTRACTION_READY" else "source SHA mismatch or silent generic row",
        "claim_ceiling": "Final-section constant extraction only. This packet names blockers and snippets but does not extract Tao's global threshold or prove the full EHP114 statement.",
    }

    explicit_rows = [row for row in extraction_rows if row["extraction_status"] == EXPLICIT_READY and row["extracted_N_i"] is not None]
    partial_max_n = max((row["extracted_N_i"] for row in explicit_rows), default=None)
    if blocker_count == len(extraction_rows):
        l40b_status = "FINAL_SECTION_REMAINS_OPAQUE"
    elif blocker_count > 0:
        l40b_status = "FINAL_SECTION_PARTIAL_EXPLICIT"
    else:
        l40b_status = "FINAL_SECTION_EXPLICIT"
    l40b = {
        "experiment_id": L40B_ID,
        "status": l40b_status,
        "source_experiment_id": L40A_ID,
        "target_dependency_count": len(extraction_rows),
        "explicit_dependency_count": len(explicit_rows),
        "named_blocker_count": blocker_count,
        "partial_max_N_final_section": partial_max_n,
        "global_candidate_N0": None,
        "subcheck_rows": [
            {
                "dependency_id": row["dependency_id"],
                "extraction_status": row["extraction_status"],
                "extracted_N_i": row["extracted_N_i"],
                "named_blocker": row["named_blocker"],
            }
            for row in extraction_rows
        ],
        "first_failed_condition": "final-section constants remain opaque" if blocker_count else "none",
        "claim_ceiling": "Final-section threshold subcheck only. It emits no global candidate_N0 and authorizes no higher-degree computation.",
    }

    l40c_status = (
        "BRIDGE_DECISION_REFRESH_FINAL_SECTION_OPAQUE"
        if l40b_status != "FINAL_SECTION_EXPLICIT"
        else "BRIDGE_DECISION_REFRESH_BACKWARD_PROPAGATION_READY"
    )
    l40c = {
        "experiment_id": L40C_ID,
        "status": l40c_status,
        "finite_n_less15_status": finite["status"],
        "finite_n_less15_pass": finite["all_n_less_15_pass"],
        "tao_threshold_status": l38["status"],
        "final_section_subcheck_status": l40b_status,
        "partial_max_N_final_section": partial_max_n,
        "n15_or_higher_computation_authorized": False,
        "full_proof_synthesis_allowed": False,
        "next_action": (
            "Resolve named blockers in final-section constants before back-propagating to inside, annulus, outside, and geomcontrol."
            if l40b_status != "FINAL_SECTION_EXPLICIT"
            else "Continue backward extraction to inside, annulus, outside, and geomcontrol."
        ),
        "first_failed_condition": "final-section constants remain opaque" if l40b_status != "FINAL_SECTION_EXPLICIT" else "none",
        "claim_ceiling": "Bridge decision refresh only. It does not prove the full all-degree statement and does not authorize n=15 or higher computation.",
    }

    artifacts = [
        (L40A_ID, l40a, l40a_report(l40a)),
        (L40B_ID, l40b, l40b_report(l40b)),
        (L40C_ID, l40c, l40c_report(l40c)),
    ]
    emitted = []
    for experiment_id, result, report in artifacts:
        outdir = proof_root / experiment_id
        digest = write_artifact(outdir, experiment_id, result, report)
        emitted.append({"experiment_id": experiment_id, "status": result["status"], "sha256": digest, "outdir": str(outdir)})

    print(json.dumps({"status": "EHP114_FINAL_SECTION_PACKETS_EMITTED", "artifact_count": len(emitted), "artifacts": emitted}, indent=2))


if __name__ == "__main__":
    main()
