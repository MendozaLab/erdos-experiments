#!/usr/bin/env python3
"""Local Mihnea/Bogdan-style #242 comparison with explicit certificates.

This runner does not reproduce the 10^18 computational frontier. It implements
the first local comparison target: Salez-style modular-filter coverage for
hard-strip primes up to 1,000,000, with an explicit Rosati witness attached to
every certified target whenever available.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import salez_filter_certificates as certs
import salez_reference_equation_audit as reference_audit


EXP_ID = "EXP-MATH-ERDOS242-MIHNEA-BOGDAN-CERT-COMPARISON-20260508-01"
ROOT = Path(__file__).resolve().parents[3]
RESULT_DIR = ROOT / "erdos-experiments" / "results" / "erdos-242"


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


def timed_summary(
    *,
    max_n: int,
    constant_bound: int,
    allowed_filter_moduli: set[int] | None = None,
) -> dict[str, Any]:
    started = time.perf_counter()
    certificates = certs.generate_filter_certificates(
        max_n=max_n,
        constant_bound=constant_bound,
        hard_strip_only=True,
        allowed_filter_moduli=allowed_filter_moduli,
    )
    summary = certs.summarize_certificates(certificates, max_n=max_n, hard_strip_only=True)
    summary["verdict"] = certs.verdict_for_summary(summary)
    summary["max_n"] = max_n
    summary["constant_bound"] = constant_bound
    summary["allowed_filter_moduli"] = sorted(allowed_filter_moduli) if allowed_filter_moduli is not None else None
    summary["runtime_seconds"] = round(time.perf_counter() - started, 6)
    return summary


def unit_checks(reference_t_max: int) -> dict[str, Any]:
    reference = reference_audit.audit(reference_t_max)
    expected_equations = {"eqmod1a", "eqmod1b", "eqmod1c", "eqmod2a", "eqmod2b", "eqmod2c", "eqmod2d"}
    observed_equations = set(reference["verified_counts_by_equation"])
    return {
        "salez_example_1_t_max": reference_t_max,
        "salez_example_1_verdict": reference["verdict"],
        "salez_example_1_failure_count": reference["failure_count"],
        "all_seven_equations_present": observed_equations == expected_equations,
        "verified_counts_by_equation": reference["verified_counts_by_equation"],
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    units = unit_checks(args.reference_t_max)
    regression_100k = timed_summary(max_n=100000, constant_bound=args.regression_bound)
    main = timed_summary(max_n=args.max_n, constant_bound=args.constant_bound)

    controls: dict[str, Any] = {}
    if not args.skip_controls:
        controls["near_miss_lower_bound"] = timed_summary(max_n=args.max_n, constant_bound=args.near_miss_bound)
        controls["restricted_modulus"] = timed_summary(
            max_n=min(args.max_n, args.control_max_n),
            constant_bound=args.constant_bound,
            allowed_filter_moduli=set(args.restricted_moduli),
        )

    verdict = main["verdict"]
    if units["salez_example_1_failure_count"] or not units["all_seven_equations_present"]:
        verdict = "FILTER_PARITY_WITNESS_GAP"
    if regression_100k["verdict"] != "CERTIFIED_PARITY_LOCAL":
        verdict = "NEEDS_FILTER_COVERAGE"

    return {
        "experiment_id": EXP_ID,
        "generated_at": now_iso(),
        "status": "REVIEW_ONLY",
        "verdict": verdict,
        "claim_ceiling": "local certificate-parity comparison only; not a 10^18 verification and not a proof of Erdos-Straus",
        "comparison_target": {
            "problem": "Erdos #242 / Erdos-Straus",
            "target_strip": "prime p = 1 mod 24",
            "max_n": args.max_n,
            "constant_bound": args.constant_bound,
            "priority_axis": "clearer certificates, not speed",
        },
        "source_manifest": {
            "local_implementation_status": "RECONSTRUCTED_FROM_SOURCES",
            "third_party_code_committed": False,
            "salez_arxiv": "https://arxiv.org/abs/1406.6307",
            "salez_source_bundle_note": "arXiv source bundle includes erdos_straus_en.tex and anc/program.cpp; no third-party code was vendored.",
            "mihnea_bogdan_arxiv": "https://arxiv.org/abs/2509.00128",
            "mihnea_bogdan_github": "https://github.com/esc-paper/erdos-straus",
        },
        "unit_checks": units,
        "regression_100k": regression_100k,
        "main_comparison": main,
        "controls": controls,
        "interpretation": {
            "mihnea_bogdan_axis": "computational frontier and large-scale filter verification",
            "local_axis": "certificate-bearing modular-filter reconstruction",
            "what_this_shows": "At the local acceptance scale, every hard-strip prime certified by the filter layer also has an explicit verified Rosati witness.",
            "what_this_does_not_show": "No claim to 10^18 scale, no speed superiority, no public proof status.",
        },
        "forbidden_actions": [
            "do not claim proof of Erdos-Straus",
            "do not claim global computational SOTA",
            "do not update D1/proof registry from this review-only artifact",
            "do not treat reconstructed filters as a vendored Mihnea/Bogdan implementation",
        ],
    }


def report(results: dict[str, Any]) -> str:
    main = results["main_comparison"]
    regression = results["regression_100k"]
    controls = results["controls"]
    lines = [
        "# Erdos #242 Mihnea/Bogdan Certificate Comparison",
        "",
        f"Experiment: `{results['experiment_id']}`",
        "Status: `REVIEW_ONLY`",
        f"Verdict: `{results['verdict']}`",
        "",
        "## Meaning",
        "",
        "This is a local comparison against the Salez/Mihnea-Bogdan filter lineage on the certificate axis. It does not compete with the `10^18` computational frontier. It asks whether the filter layer can be made witness-bearing at a modest scale.",
        "",
        "## Main Result",
        "",
        f"- Target: hard-strip primes `p = 1 mod 24` up to `{main['max_n']}`",
        f"- Constant bound: `{main['constant_bound']}`",
        f"- Certified targets: `{main['certified_target_count']} / {main['target_prime_count']} = {main['coverage_fraction']:.6f}`",
        f"- Filter-only/no-witness targets: `{main['filter_only_no_witness_count']}`",
        f"- Invalid certificates: `{main['invalid_certificate_count']}`",
        f"- Certificate count: `{main['certificate_count']}`",
        f"- Filter moduli used: `{main['filter_moduli_used_count']}`",
        f"- Runtime seconds: `{main['runtime_seconds']}`",
        "",
        "## Regression And Controls",
        "",
        f"- 100k regression: `{regression['certified_target_count']} / {regression['target_prime_count']}` with verdict `{regression['verdict']}`",
    ]
    if controls:
        lower = controls["near_miss_lower_bound"]
        restricted = controls["restricted_modulus"]
        lines.extend([
            f"- Near-miss lower bound control: bound `{lower['constant_bound']}` covers `{lower['certified_target_count']} / {lower['target_prime_count']}` with verdict `{lower['verdict']}`",
            f"- Restricted-modulus control: moduli `{restricted['allowed_filter_moduli']}` covers `{restricted['certified_target_count']} / {restricted['target_prime_count']}` with verdict `{restricted['verdict']}`",
        ])
    lines.extend([
        "",
        "## Boundary",
        "",
        "Mihnea/Bogdan remain the computational-frontier comparison point. This packet's contribution is different: each local certified target carries a concrete Rosati witness and verified denominators.",
        "",
    ])
    return "\n".join(lines)


def parse_moduli(raw: str) -> list[int]:
    values = [int(part.strip()) for part in raw.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("at least one modulus is required")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=1_000_000)
    parser.add_argument("--constant-bound", type=int, default=21)
    parser.add_argument("--regression-bound", type=int, default=10)
    parser.add_argument("--near-miss-bound", type=int, default=20)
    parser.add_argument("--reference-t-max", type=int, default=1000)
    parser.add_argument("--restricted-moduli", type=parse_moduli, default=parse_moduli("3"))
    parser.add_argument("--control-max-n", type=int, default=100000)
    parser.add_argument("--skip-controls", action="store_true")
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
        "main": {
            "max_n": results["main_comparison"]["max_n"],
            "constant_bound": results["main_comparison"]["constant_bound"],
            "certified": results["main_comparison"]["certified_target_count"],
            "target": results["main_comparison"]["target_prime_count"],
            "coverage": results["main_comparison"]["coverage_fraction"],
            "filter_only_no_witness": results["main_comparison"]["filter_only_no_witness_count"],
            "invalid_certificates": results["main_comparison"]["invalid_certificate_count"],
        },
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
