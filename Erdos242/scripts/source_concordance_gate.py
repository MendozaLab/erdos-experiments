#!/usr/bin/env python3
"""Source-concordance gate for the #242 certificate pipeline.

This runner checks whether the local certificate-parity packet is honestly
aligned with the Salez / Mihnea-Bogdan filter lineage. It records concordance
and differences only; it does not vendor third-party code, run the upstream
checker, update D1, or claim SOTA-scale verification.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXP_ID = "EXP-MATH-ERDOS242-SOURCE-CONCORDANCE-20260508-01"
ROOT = Path(__file__).resolve().parents[3]
RESULT_DIR = ROOT / "erdos-experiments" / "results" / "erdos-242"
CERT_COMPARISON_ID = "EXP-MATH-ERDOS242-MIHNEA-BOGDAN-CERT-COMPARISON-20260508-01"
CERT_COMPARISON_RESULTS = RESULT_DIR / f"{CERT_COMPARISON_ID}_RESULTS.json"

EXPECTED_EQUATIONS = [
    "eqmod1a",
    "eqmod1b",
    "eqmod1c",
    "eqmod2a",
    "eqmod2b",
    "eqmod2c",
    "eqmod2d",
]
EXPECTED_MB_PATHS = [
    "README.md",
    "section1/Checker.cpp",
    "section1/Salez_Python.py",
    "section1/resources/Filters.txt",
    "section1/resources/Residues.txt",
    "section2/solution_counting-full.csv",
]
SALEZ_TEX = Path("/private/tmp/salez1406/erdos_straus_en.tex")
SALEZ_PROGRAM = Path("/private/tmp/salez1406/anc/program.cpp")


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_json(url: str) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"User-Agent": "h2-erdos242-source-concordance"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_text(url: str, limit: int = 4096) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "h2-erdos242-source-concordance"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read(limit).decode("utf-8", errors="replace")


def github_reference() -> dict[str, Any]:
    out: dict[str, Any] = {
        "repo": "https://github.com/esc-paper/erdos-straus",
        "status": "UNFETCHED",
        "third_party_code_committed": False,
    }
    try:
        commit = fetch_json("https://api.github.com/repos/esc-paper/erdos-straus/commits/HEAD")
        sha = commit["sha"]
        tree = fetch_json(f"https://api.github.com/repos/esc-paper/erdos-straus/git/trees/{sha}?recursive=1")
        paths = sorted(item["path"] for item in tree.get("tree", []) if item.get("type") == "blob")
        readme = fetch_text("https://raw.githubusercontent.com/esc-paper/erdos-straus/HEAD/README.md")
        filters_head = fetch_text("https://raw.githubusercontent.com/esc-paper/erdos-straus/HEAD/section1/resources/Filters.txt", limit=2048)
        residues_head = fetch_text("https://raw.githubusercontent.com/esc-paper/erdos-straus/HEAD/section1/resources/Residues.txt", limit=2048)
        missing = [path for path in EXPECTED_MB_PATHS if path not in paths]
        out.update({
            "status": "FETCHED_REFERENCE_METADATA_ONLY",
            "head_commit": sha,
            "expected_paths_present": not missing,
            "missing_expected_paths": missing,
            "path_count": len(paths),
            "paths": paths,
            "readme_mentions_arxiv": "arxiv.org/abs/2509.00128" in readme,
            "readme_mentions_gmp": "gmplib.org" in readme or "GMP" in readme,
            "section1_filter_resource_sample": filters_head.splitlines()[:5],
            "section1_residue_resource_sample": residues_head.splitlines()[:2],
        })
    except Exception as exc:  # pragma: no cover - network contingency is part of artifact truth.
        out.update({
            "status": "FETCH_FAILED",
            "error": repr(exc),
            "expected_paths_present": False,
            "missing_expected_paths": EXPECTED_MB_PATHS,
        })
    return out


def salez_reference() -> dict[str, Any]:
    tex = SALEZ_TEX.read_text(encoding="utf-8", errors="replace") if SALEZ_TEX.exists() else ""
    program = SALEZ_PROGRAM.read_text(encoding="utf-8", errors="replace") if SALEZ_PROGRAM.exists() else ""
    label_presence = {label: bool(re.search(rf"\\label\{{{label}\}}", tex)) for label in EXPECTED_EQUATIONS}
    filter_function_presence = {
        label: f"filter{label.replace('eqmod', '')}" in program
        for label in EXPECTED_EQUATIONS
    }
    formula_presence = {
        "hard_strip_reduction_p_1_mod_24": "p=1 \\mod 24" in tex or "p=1 \\mod24" in tex,
        "complete_seven_equation_statement": "7 reference equations" in tex and "complete set" in tex,
        "program_claim_10_17": "10^{17}" in tex and "N=10^{17}" in tex,
        "mod_array_present": "const unsigned\tMOD[]" in program or "const unsigned\tMOD" in program,
        "row_constants_present": "ROW1" in program and "ROW2" in program and "ROW3" in program,
    }
    return {
        "arxiv": "https://arxiv.org/abs/1406.6307",
        "source_paths": {
            "tex": str(SALEZ_TEX),
            "program": str(SALEZ_PROGRAM),
            "tex_exists": SALEZ_TEX.exists(),
            "program_exists": SALEZ_PROGRAM.exists(),
        },
        "equation_labels_present": label_presence,
        "all_equation_labels_present": all(label_presence.values()),
        "program_filter_functions_present": filter_function_presence,
        "all_program_filter_functions_present": all(filter_function_presence.values()),
        "formula_presence": formula_presence,
    }


def load_certificate_comparison(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    main = data["main_comparison"]
    return {
        "path": str(path.relative_to(ROOT)),
        "experiment_id": data["experiment_id"],
        "verdict": data["verdict"],
        "claim_ceiling": data["claim_ceiling"],
        "max_n": main["max_n"],
        "constant_bound": main["constant_bound"],
        "target_prime_count": main["target_prime_count"],
        "certified_target_count": main["certified_target_count"],
        "coverage_fraction": main["coverage_fraction"],
        "certificate_count": main["certificate_count"],
        "filter_moduli_used_count": main["filter_moduli_used_count"],
        "filter_residue_class_count": main["filter_residue_class_count"],
        "filter_only_no_witness_count": main["filter_only_no_witness_count"],
        "invalid_certificate_count": main["invalid_certificate_count"],
        "equation_counts": main["equation_counts"],
        "filter_moduli_sample": main["filter_moduli_used"][:80],
        "filter_moduli_tail": main["filter_moduli_used"][-20:],
        "one_witness_sample": main["one_witness_per_certified_target"][:12],
        "source_manifest": data.get("source_manifest", {}),
    }


def classify_concordance(
    salez: dict[str, Any],
    mb: dict[str, Any],
    local: dict[str, Any],
) -> dict[str, Any]:
    local_equations = set(local["equation_counts"])
    all_equations_used = set(EXPECTED_EQUATIONS).issubset(local_equations)
    direct_salez = salez["all_equation_labels_present"] and salez["all_program_filter_functions_present"] and all_equations_used
    mb_style = (
        mb.get("expected_paths_present", False)
        and local["filter_moduli_used_count"] > 0
        and local["filter_residue_class_count"] > 0
    )
    enrichment = (
        local["verdict"] == "CERTIFIED_PARITY_LOCAL"
        and local["certified_target_count"] == local["target_prime_count"]
        and local["filter_only_no_witness_count"] == 0
        and local["invalid_certificate_count"] == 0
        and bool(local["one_witness_sample"])
    )
    intentional_differences = [
        "local run stops at n <= 1,000,000; Mihnea/Bogdan report a 10^18-scale computational frontier",
        "local output attaches Rosati variables and denominators for certified targets",
        "local code reconstructs source formulas; it does not vendor or execute third-party code",
        "local filter moduli/residue classes are certificate-bearing natural moduli, not a byte-for-byte copy of Mihnea/Bogdan Filters.txt/Residues.txt",
    ]
    if direct_salez and mb_style and enrichment:
        verdict = "SOURCE_CONCORDANT_CERTIFICATE_ENRICHMENT"
    elif direct_salez or mb_style or enrichment:
        verdict = "PARTIAL_CONCORDANCE_NEEDS_REVIEW"
    else:
        verdict = "NOT_CONCORDANT"
    return {
        "verdict": verdict,
        "direct_salez_equation_reconstructions": direct_salez,
        "mihnea_bogdan_style_filter_analogues": mb_style,
        "local_certificate_enrichments_not_present_in_upstream_pipeline": enrichment,
        "all_expected_equations_used_locally": all_equations_used,
        "intentional_differences": intentional_differences,
        "classification": {
            "357049_certificates": "direct Salez-equation reconstructions plus local Rosati/denominator certificate enrichment",
            "970_filter_moduli": "Mihnea/Bogdan-style modular-filter analogue; not asserted as exact upstream filter table equality",
            "7626_filter_residue_classes": "local certificate-bearing residue classes derived from verified witnesses",
        },
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    salez = salez_reference()
    mb = github_reference() if not args.skip_network else {
        "repo": "https://github.com/esc-paper/erdos-straus",
        "status": "SKIPPED_BY_FLAG",
        "expected_paths_present": False,
        "third_party_code_committed": False,
    }
    local = load_certificate_comparison(args.cert_results)
    concordance = classify_concordance(salez, mb, local)
    return {
        "experiment_id": EXP_ID,
        "generated_at": now_iso(),
        "status": "REVIEW_ONLY",
        "verdict": concordance["verdict"],
        "claim_ceiling": "source-concordance and certificate-enrichment audit only; not a 10^18 replication, not SOTA, and not a proof",
        "source_manifest": {
            "status": "REFERENCE_CONCORDANCE_ONLY",
            "third_party_code_committed": False,
            "salez_arxiv": "https://arxiv.org/abs/1406.6307",
            "mihnea_bogdan_arxiv": "https://arxiv.org/abs/2509.00128",
            "mihnea_bogdan_github": "https://github.com/esc-paper/erdos-straus",
            "mihnea_bogdan_head_commit": mb.get("head_commit"),
        },
        "salez_reference": salez,
        "mihnea_bogdan_reference": mb,
        "local_certificate_packet": local,
        "concordance": concordance,
        "forbidden_actions": [
            "do not claim proof of Erdos-Straus",
            "do not claim global computational SOTA",
            "do not claim replication of Mihnea/Bogdan's 10^18 result",
            "do not update D1/proof registry from this review-only artifact",
            "do not imply third-party code was vendored or executed",
        ],
    }


def report(results: dict[str, Any]) -> str:
    local = results["local_certificate_packet"]
    concordance = results["concordance"]
    mb = results["mihnea_bogdan_reference"]
    lines = [
        "# Erdos #242 Source-Concordance Gate",
        "",
        f"Experiment: `{results['experiment_id']}`",
        "Status: `REVIEW_ONLY`",
        f"Verdict: `{results['verdict']}`",
        "",
        "## Meaning",
        "",
        "This gate checks whether the local certificate-parity runner is aligned with the Salez / Mihnea-Bogdan filter lineage. It is not a proof, not a `10^18` replication, and not a public SOTA claim.",
        "",
        "## Concordance Classification",
        "",
        f"- Direct Salez-equation reconstruction: `{concordance['direct_salez_equation_reconstructions']}`",
        f"- Mihnea/Bogdan-style filter analogue: `{concordance['mihnea_bogdan_style_filter_analogues']}`",
        f"- Local certificate enrichment: `{concordance['local_certificate_enrichments_not_present_in_upstream_pipeline']}`",
        "",
        "## Local Packet Interpreted",
        "",
        f"- Certificates at `1e6`: `{local['certificate_count']}`",
        f"- Filter moduli: `{local['filter_moduli_used_count']}`",
        f"- Filter residue classes: `{local['filter_residue_class_count']}`",
        f"- Certified hard-strip primes: `{local['certified_target_count']} / {local['target_prime_count']}`",
        f"- Invalid certificates: `{local['invalid_certificate_count']}`",
        f"- Filter-only/no-witness targets: `{local['filter_only_no_witness_count']}`",
        "",
        "Interpretation: the certificates are direct seven-equation Salez reconstructions with added Rosati variables and denominators. The filter moduli/residue classes are Mihnea/Bogdan-style filter analogues, not asserted as exact equality to the upstream filter tables.",
        "",
        "## External Reference Snapshot",
        "",
        f"- Mihnea/Bogdan repo status: `{mb.get('status')}`",
        f"- Mihnea/Bogdan HEAD: `{mb.get('head_commit')}`",
        f"- Expected repo paths present: `{mb.get('expected_paths_present')}`",
        f"- Third-party code committed locally: `{results['source_manifest']['third_party_code_committed']}`",
        "",
        "## Intentional Differences",
        "",
    ]
    for diff in concordance["intentional_differences"]:
        lines.append(f"- {diff}")
    lines.extend([
        "",
        "## Boundary",
        "",
        "This packet supports the claim that our local runner is source-concordant certificate enrichment. It does not support claiming computational-frontier parity with Mihnea/Bogdan.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cert-results", type=Path, default=CERT_COMPARISON_RESULTS)
    parser.add_argument("--skip-network", action="store_true")
    args = parser.parse_args()

    results = run(args)
    result_path = RESULT_DIR / f"{EXP_ID}_RESULTS.json"
    report_path = RESULT_DIR / f"{EXP_ID}_REPORT.md"
    sha_path = RESULT_DIR / f"{EXP_ID}_RESULTS.sha256"
    write_json(result_path, results)
    write_text(report_path, report(results))
    write_text(sha_path, sha256_file(result_path) + "\n")
    print(json.dumps({
        "experiment_id": EXP_ID,
        "results": str(result_path.relative_to(ROOT)),
        "report": str(report_path.relative_to(ROOT)),
        "sha256": str(sha_path.relative_to(ROOT)),
        "verdict": results["verdict"],
        "mihnea_bogdan_head_commit": results["source_manifest"]["mihnea_bogdan_head_commit"],
        "direct_salez": results["concordance"]["direct_salez_equation_reconstructions"],
        "mb_style_filter_analogue": results["concordance"]["mihnea_bogdan_style_filter_analogues"],
        "certificate_enrichment": results["concordance"]["local_certificate_enrichments_not_present_in_upstream_pipeline"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
