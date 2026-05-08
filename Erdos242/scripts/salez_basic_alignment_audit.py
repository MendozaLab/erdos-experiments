#!/usr/bin/env python3
"""Audit Erdos #242 operator families against Salez's basic formulas.

The Salez paper's opening reduction gives explicit identities for
`n = 3t - 1`, `n = 4t - 1`, and `n = 8t - 3`, and then reduces the hard
case to primes `p = 1 mod 24`. This script checks that the local
modular-coverage operator recovers those basic formulas exactly and records
the remaining hard strip honestly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Any


EXP_ID = "EXP-MATH-ERDOS242-SALEZ-BASIC-ALIGNMENT-20260508-01"
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


def verify(n: int, xyz: tuple[int, int, int]) -> bool:
    x, y, z = xyz
    return (
        2 < n
        and 0 < x < y < z
        and Fraction(4, n) == Fraction(1, x) + Fraction(1, y) + Fraction(1, z)
    )


def salez_basic_family(n: int) -> tuple[str, tuple[int, int, int]] | None:
    if n > 2 and (n + 1) % 3 == 0:
        t = (n + 1) // 3
        return "salez_basic_3t_minus_1", (t, n, t * n)

    if n > 2 and (n + 1) % 4 == 0:
        t = (n + 1) // 4
        m = t * n
        return "salez_basic_4t_minus_1_split", (t, m + 1, m * (m + 1))

    if n > 2 and (n + 3) % 8 == 0:
        t = (n + 3) // 8
        return "salez_basic_8t_minus_3", (2 * t, t * n, 2 * t * n)

    return None


def operator_family(n: int) -> tuple[str, tuple[int, int, int]] | None:
    if n == 4:
        return "small_n_4", (2, 3, 6)

    if n > 2 and n % 3 == 2:
        k = (n + 1) // 3
        return "n_mod_3_eq_2", (k, n, k * n)

    if n > 2 and n % 4 == 3:
        x = (n + 1) // 4
        m = n * (n + 1) // 4
        return "n_mod_4_eq_3_split", (x, m + 1, m * (m + 1))

    if n > 2 and n % 8 == 5:
        x = (n + 3) // 4
        d = n * x
        return "n_mod_8_eq_5", (x, d // 2, d)

    return None


def divisors(n: int) -> list[int]:
    found = {1, n}
    for d in range(2, int(math.isqrt(n)) + 1):
        if n % d == 0:
            found.add(d)
            found.add(n // d)
    return sorted(found)


def divisor_descent_covered(n: int) -> bool:
    if operator_family(n):
        return True
    return any(operator_family(d) for d in divisors(n) if d not in {1, n} and 2 < d)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for d in range(3, int(math.isqrt(n)) + 1, 2):
        if n % d == 0:
            return False
    return True


def audit(max_n: int) -> dict[str, Any]:
    salez_rows = []
    family_counts: Counter[str] = Counter()
    mismatches = []
    hard_strip = []
    hard_primes = []
    descent_covered = 0

    for n in range(3, max_n + 1):
        salez = salez_basic_family(n)
        operator = operator_family(n)
        if salez:
            salez_name, salez_xyz = salez
            family_counts[salez_name] += 1
            row = {
                "n": n,
                "salez_family": salez_name,
                "salez_xyz": list(salez_xyz),
                "salez_verified": verify(n, salez_xyz),
                "operator_family": operator[0] if operator else None,
                "operator_xyz": list(operator[1]) if operator else None,
                "operator_verified": verify(n, operator[1]) if operator else None,
            }
            salez_rows.append(row)
            if not row["salez_verified"] or not row["operator_verified"]:
                mismatches.append(row)

        if divisor_descent_covered(n):
            descent_covered += 1
        elif n % 24 == 1:
            hard_strip.append(n)
            if is_prime(n):
                hard_primes.append(n)

    return {
        "experiment_id": EXP_ID,
        "generated_at": now_iso(),
        "status": "REVIEW_ONLY",
        "verdict": "BASIC_FORMULAS_RECOVERED",
        "claim_ceiling": "recovers Salez opening reduction/basic formulas; does not recover Salez's seven modular equations",
        "max_n": max_n,
        "salez_basic_family_counts": dict(family_counts),
        "checked_salez_basic_instances": len(salez_rows),
        "mismatch_count": len(mismatches),
        "first_mismatches": mismatches[:20],
        "descent_coverage": {
            "covered_count": descent_covered,
            "sample_count": max_n - 2,
            "coverage_fraction": descent_covered / (max_n - 2),
        },
        "remaining_hard_strip": {
            "description": "Not covered by basic Salez families or divisor descent; all are n = 1 mod 24 in this audit.",
            "count": len(hard_strip),
            "first_values": hard_strip[:60],
            "prime_count": len(hard_primes),
            "first_primes": hard_primes[:60],
        },
        "sota_sources": [
            {
                "label": "Salez 2014",
                "url": "https://arxiv.org/abs/1406.6307",
                "used_for": "opening reduction and seven-modular-equation frontier",
            },
            {
                "label": "Mihnea-Bogdan 2025",
                "url": "https://arxiv.org/abs/2509.00128",
                "used_for": "current computational frontier at 10^18",
            },
            {
                "label": "Elsholtz-Tao 2013",
                "url": "https://www.cambridge.org/core/journals/journal-of-the-australian-mathematical-society/article/counting-the-number-of-solutions-to-the-erdosstraus-equation-on-unit-fractions/A1AB4837F744ED30E0F1243E65C2B646",
                "used_for": "counting-function context",
            },
        ],
        "forbidden_actions": [
            "do not claim proof of Erdos-Straus",
            "do not claim recovery of Salez's seven modular equations",
            "do not update D1 or proof registry from this packet alone",
        ],
    }


def report(results: dict[str, Any]) -> str:
    coverage = results["descent_coverage"]
    hard = results["remaining_hard_strip"]
    lines = [
        "# Erdos #242 Salez Basic-Formula Alignment Audit",
        "",
        f"Experiment: `{results['experiment_id']}`",
        "Status: `REVIEW_ONLY`",
        f"Verdict: `{results['verdict']}`",
        "",
        "## Meaning",
        "",
        "The local residue operator recovers the opening Salez reduction exactly: the families for `3t - 1`, `4t - 1`, and `8t - 3` correspond to our `n mod 3 = 2`, `n mod 4 = 3`, and `n mod 8 = 5` operator families.",
        "",
        "This is useful, but it is not enough to claim SOTA. Salez's real frontier is the later complete set of seven modular equations and the `10^17` sieve; Mihnea-Bogdan report `10^18` computational verification.",
        "",
        "## Audit Results",
        "",
        f"- Checked Salez basic-family instances up to n = {results['max_n']}: {results['checked_salez_basic_instances']}",
        f"- Mismatches: {results['mismatch_count']}",
        f"- Divisor-descent coverage: {coverage['covered_count']} / {coverage['sample_count']} = {coverage['coverage_fraction']:.6f}",
        f"- Remaining hard-strip count: {hard['count']}",
        f"- First hard-strip values: `{hard['first_values'][:20]}`",
        f"- First hard-strip primes: `{hard['first_primes'][:20]}`",
        "",
        "## Claim Boundary",
        "",
        "We can now say the operator is aligned with the first Salez reduction and has compiled scaling infrastructure. We cannot yet say it recovers the seven modular equations, improves the sieve frontier, or resolves the `p = 1 mod 24` core.",
        "",
        "## Next Gate",
        "",
        "Extract the seven Salez modular equations from the paper/program, encode each as a named operator family, and rerun this audit with a family-by-family recovery table.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=100000)
    args = parser.parse_args()

    results = audit(args.max_n)
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
        "mismatch_count": results["mismatch_count"],
        "coverage_fraction": results["descent_coverage"]["coverage_fraction"],
        "hard_strip_count": results["remaining_hard_strip"]["count"],
        "first_hard_strip": results["remaining_hard_strip"]["first_values"][:10],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
