#!/usr/bin/env python3
"""Verify Salez's seven reference-equation example family.

Salez gives seven reference modular equations and an example family
`p = 24 * 5 * t - 23` where all seven are present. This script encodes those
seven labels, verifies the explicit decompositions for a finite range, and
records the boundary: this is a direct reference-equation alignment, not a new
global sieve or proof.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Any


EXP_ID = "EXP-MATH-ERDOS242-SALEZ-SEVEN-EQUATION-EXAMPLE-20260508-01"
ROOT = Path(__file__).resolve().parents[3]
RESULT_DIR = ROOT / "erdos-experiments" / "results" / "erdos-242"


REFERENCE_EQUATIONS = [
    {
        "id": "eqmod1a",
        "paper_label": r"\ref{eqmod1a}",
        "condition": "B + p*C = 0 mod 4*B*C*D - 1",
    },
    {
        "id": "eqmod1b",
        "paper_label": r"\ref{eqmod1b}",
        "condition": "p + E = 0 mod 4*A*B and A + B = 0 mod E",
    },
    {
        "id": "eqmod1c",
        "paper_label": r"\ref{eqmod1c}",
        "condition": "p + E + 4*B^2*D = 0 mod 4*B*D*E",
    },
    {
        "id": "eqmod2a",
        "paper_label": r"\ref{eqmod2a}",
        "condition": "p*E + 1 = 0 mod 4*A*B and A + B = 0 mod E",
    },
    {
        "id": "eqmod2b",
        "paper_label": r"\ref{eqmod2b}",
        "condition": "p + F = 0 mod 4*B*C and p*B + C = 0 mod F",
    },
    {
        "id": "eqmod2c",
        "paper_label": r"\ref{eqmod2c}",
        "condition": "p + F = 0 mod 4*B*D and 4*B^2*D + 1 = 0 mod F",
    },
    {
        "id": "eqmod2d",
        "paper_label": r"\ref{eqmod2d}",
        "condition": "p + F = 0 mod 4*C*D and p^2 + 4*C^2*D = 0 mod F",
    },
]


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


def p_of_t(t: int) -> int:
    return 120 * t - 23


def denominators(label: str, t: int) -> list[int]:
    p = p_of_t(t)
    if label == "eqmod1a":
        return [4 * p, 4 * (16 * t - 3) * p, 2 * (16 * t - 3)]
    if label == "eqmod1b":
        return [10 * (6 * t - 1) * p, 2 * (6 * t - 1) * p, 5 * (6 * t - 1)]
    if label == "eqmod1c":
        return [10 * t * p, 10 * t * (6 * t - 1) * p, 5 * (6 * t - 1)]
    if label == "eqmod2a":
        return [5 * (21 * t - 4), 2 * (21 * t - 4), 10 * (21 * t - 4) * p]
    if label == "eqmod2b":
        return [
            5 * (6 * t - 1),
            2 * (6 * t - 1) * (100 * t - 19),
            10 * (6 * t - 1) * (100 * t - 19) * p,
        ]
    if label == "eqmod2c":
        return [
            5 * (6 * t - 1),
            10 * (6 * t - 1) * (21 * t - 4),
            10 * (21 * t - 4) * p,
        ]
    if label == "eqmod2d":
        q = 120 * t * t - 43 * t + 4
        return [5 * (6 * t - 1), 10 * q, 10 * (6 * t - 1) * q * p]
    raise KeyError(label)


def verify(p: int, ds: list[int]) -> bool:
    return (
        2 < p
        and len(set(ds)) == 3
        and all(d > 0 for d in ds)
        and Fraction(4, p) == sum(Fraction(1, d) for d in ds)
    )


def audit(t_max: int) -> dict[str, Any]:
    failures = []
    examples: dict[str, Any] = {}
    counts = {row["id"]: 0 for row in REFERENCE_EQUATIONS}

    for t in range(1, t_max + 1):
        p = p_of_t(t)
        for equation in REFERENCE_EQUATIONS:
            label = equation["id"]
            ds = denominators(label, t)
            ok = verify(p, ds)
            if ok:
                counts[label] += 1
                examples.setdefault(label, {
                    "t": t,
                    "p": p,
                    "denominators": sorted(ds),
                    "p_mod_24": p % 24,
                    "p_mod_5": p % 5,
                })
            else:
                failures.append({"t": t, "p": p, "equation": label, "denominators": ds})

    return {
        "experiment_id": EXP_ID,
        "generated_at": now_iso(),
        "status": "REVIEW_ONLY",
        "verdict": "SALEZ_SEVEN_EXAMPLE_VERIFIED" if not failures else "FAIL",
        "claim_ceiling": "verifies Salez Example 1 seven-equation family; not a proof and not a full sieve implementation",
        "t_range": [1, t_max],
        "p_family": "p = 120*t - 23; hence p = 1 mod 24 and p = 2 mod 5",
        "reference_equations": REFERENCE_EQUATIONS,
        "verified_counts_by_equation": counts,
        "example_witness_by_equation": examples,
        "failure_count": len(failures),
        "first_failures": failures[:20],
        "sota_context": {
            "salez_source": "https://arxiv.org/abs/1406.6307",
            "salez_source_bundle": "arXiv source bundle contained erdos_straus_en.tex and anc/program.cpp",
            "boundary": "The seven labels are encoded from Proposition eqmod; this script verifies the paper's Example 1 family, not the whole Salez sieve.",
        },
        "forbidden_actions": [
            "do not claim proof of Erdos-Straus",
            "do not claim global computational SOTA",
            "do not treat this as D1/proof-registry evidence without review",
        ],
    }


def report(results: dict[str, Any]) -> str:
    lines = [
        "# Erdos #242 Salez Seven-Equation Example Audit",
        "",
        f"Experiment: `{results['experiment_id']}`",
        "Status: `REVIEW_ONLY`",
        f"Verdict: `{results['verdict']}`",
        "",
        "## Meaning",
        "",
        "This is the first direct alignment with Salez's seven reference-equation layer. It does not implement the full sieve, but it encodes the seven equation labels and verifies Salez's Example 1 family where all seven are present.",
        "",
        "## Result",
        "",
        f"- Family checked: `{results['p_family']}`",
        f"- t range: `{results['t_range']}`",
        f"- Failure count: {results['failure_count']}",
        f"- Verified counts by equation: `{results['verified_counts_by_equation']}`",
        "",
        "## Boundary",
        "",
        "This advances the #242 packet from basic-formula alignment to seven-equation example alignment. The next step is still larger: implement each reference equation as a general constant-coefficient search family and compare its coverage with the local residue operator.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t-max", type=int, default=1000)
    args = parser.parse_args()

    results = audit(args.t_max)
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
        "failure_count": results["failure_count"],
        "verified_counts_by_equation": results["verified_counts_by_equation"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
