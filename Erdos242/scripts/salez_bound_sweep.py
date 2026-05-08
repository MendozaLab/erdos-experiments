#!/usr/bin/env python3
"""Sweep bounded Salez seven-equation coverage by constant window."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXP_ID = "EXP-MATH-ERDOS242-SALEZ-BOUND-SWEEP-20260508-01"
ROOT = Path(__file__).resolve().parents[3]
RESULT_DIR = ROOT / "erdos-experiments" / "results" / "erdos-242"
SEARCH_SCRIPT = ROOT / "erdos-experiments" / "Erdos242" / "scripts" / "salez_general_equation_search.py"


def load_search_module() -> Any:
    spec = importlib.util.spec_from_file_location("salez_general_equation_search", SEARCH_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {SEARCH_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


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


def run(max_p: int, bounds: list[int]) -> dict[str, Any]:
    search = load_search_module()
    rows = []
    for bound in bounds:
        result = search.run(max_p=max_p, constant_bound=bound, target_hard_strip=True)
        coverage = result["bounded_salez_search"]
        rows.append({
            "constant_bound": bound,
            "covered_target_prime_count": coverage["covered_target_prime_count"],
            "target_prime_count": result["prime_population"]["target_prime_count"],
            "coverage_fraction": coverage["coverage_fraction"],
            "equation_counts": coverage["equation_counts"],
            "first_uncovered_target_primes": coverage["first_uncovered_target_primes"][:20],
        })

    first_full = next((row["constant_bound"] for row in rows if row["coverage_fraction"] == 1.0), None)
    verdict = "FULL_HARD_STRIP_COVERAGE_WITHIN_SWEEP" if first_full is not None else "PARTIAL_WITHIN_SWEEP"
    return {
        "experiment_id": EXP_ID,
        "generated_at": now_iso(),
        "status": "REVIEW_ONLY",
        "verdict": verdict,
        "claim_ceiling": "bounded coverage sweep for Salez seven-equation search; not an optimized sieve and not a proof",
        "max_p": max_p,
        "target": "prime hard strip p = 1 mod 24",
        "bounds": bounds,
        "first_full_coverage_bound": first_full,
        "rows": rows,
        "forbidden_actions": [
            "do not claim proof of Erdos-Straus",
            "do not claim global computational SOTA",
            "do not treat bounded sweep as asymptotic evidence",
        ],
    }


def report(results: dict[str, Any]) -> str:
    lines = [
        "# Erdos #242 Salez Bound Sweep",
        "",
        f"Experiment: `{results['experiment_id']}`",
        "Status: `REVIEW_ONLY`",
        f"Verdict: `{results['verdict']}`",
        "",
        "## Meaning",
        "",
        "This sweep asks whether the bounded seven-equation search only works at a large constant window, or whether hard-strip coverage appears steadily as the window grows.",
        "",
        "## Coverage By Bound",
        "",
        "| Constant bound | Covered hard-strip primes | Coverage | First uncovered targets |",
        "|---:|---:|---:|---|",
    ]
    for row in results["rows"]:
        lines.append(
            f"| {row['constant_bound']} | {row['covered_target_prime_count']} / {row['target_prime_count']} | "
            f"{row['coverage_fraction']:.6f} | `{row['first_uncovered_target_primes'][:8]}` |"
        )
    lines.extend([
        "",
        "## Boundary",
        "",
        "The sweep is a tractability diagnostic, not a theorem. It is useful because it shows how quickly the Salez-family reconstruction covers the hard strip under finite constants.",
        "",
    ])
    return "\n".join(lines)


def parse_bounds(raw: str) -> list[int]:
    bounds = [int(part.strip()) for part in raw.split(",") if part.strip()]
    if not bounds or any(bound <= 0 for bound in bounds):
        raise argparse.ArgumentTypeError("bounds must be a comma-separated list of positive integers")
    return sorted(set(bounds))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-p", type=int, default=100000)
    parser.add_argument("--bounds", type=parse_bounds, default=parse_bounds("5,10,20,40,80"))
    args = parser.parse_args()

    results = run(args.max_p, args.bounds)
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
        "first_full_coverage_bound": results["first_full_coverage_bound"],
        "rows": results["rows"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
