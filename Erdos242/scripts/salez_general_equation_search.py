#!/usr/bin/env python3
"""Bounded general search for Salez's seven modular equations.

This audit moves one step beyond the Example 1 check: it treats each of
Salez's seven reference equations as a constant-coefficient search family,
recovers Rosati variables A,B,C,D when the congruences fire, and verifies the
resulting 3-unit-fraction certificate directly.

It is still bounded and review-only. It is not Salez's optimized sieve, not a
proof of Erdős-Straus, and not a computational-frontier claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Iterable


EXP_ID = "EXP-MATH-ERDOS242-SALEZ-GENERAL-EQUATION-SEARCH-20260508-01"
ROOT = Path(__file__).resolve().parents[3]
RESULT_DIR = ROOT / "erdos-experiments" / "results" / "erdos-242"


@dataclass(frozen=True)
class Witness:
    equation: str
    constants: dict[str, int]
    rosati_case: str
    A: int
    B: int
    C: int
    D: int
    denominators: tuple[int, int, int]


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


def primes_up_to(n: int) -> list[int]:
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, int(math.isqrt(n)) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:n + 1:p] = b"\x00" * (((n - start) // p) + 1)
    return [i for i in range(2, n + 1) if sieve[i]]


def inv_mod(a: int, m: int) -> int | None:
    if math.gcd(a, m) != 1:
        return None
    return pow(a, -1, m)


def crt_pair(r1: int, m1: int, r2: int, m2: int) -> tuple[int, int] | None:
    g = math.gcd(m1, m2)
    if (r2 - r1) % g:
        return None
    lcm = m1 // g * m2
    k = ((r2 - r1) // g * pow(m1 // g, -1, m2 // g)) % (m2 // g)
    return (r1 + m1 * k) % lcm, lcm


def p_values_for_residue(max_p: int, residue: int, modulus: int, prime_set: set[int]) -> Iterable[int]:
    residue %= modulus
    start = residue if residue >= 3 else residue + ((3 - residue + modulus - 1) // modulus) * modulus
    for p in range(start, max_p + 1, modulus):
        if p in prime_set:
            yield p


def verify_solution(p: int, denominators: tuple[int, int, int]) -> bool:
    ds = sorted(denominators)
    return (
        p > 2
        and len(set(ds)) == 3
        and all(d > 0 for d in ds)
        and Fraction(4, p) == sum(Fraction(1, d) for d in ds)
    )


def rosati1_witness(p: int, equation: str, constants: dict[str, int], A: int, B: int, C: int, D: int) -> Witness | None:
    if min(A, B, C, D) <= 0:
        return None
    if 4 * A * B * C * D != A + B + p * C:
        return None
    if math.gcd(A * B * D, p) != 1:
        return None
    ds = (p * B * C * D, p * A * C * D, A * B * D)
    if not verify_solution(p, ds):
        return None
    return Witness(equation, constants, "Rosati_1", A, B, C, D, tuple(sorted(ds)))


def rosati2_witness(p: int, equation: str, constants: dict[str, int], A: int, B: int, C: int, D: int) -> Witness | None:
    if min(A, B, C, D) <= 0:
        return None
    if 4 * A * B * C * D != p * (A + B) + C:
        return None
    if math.gcd(A * B * C * D, p) != 1:
        return None
    ds = (B * C * D, A * C * D, p * A * B * D)
    if not verify_solution(p, ds):
        return None
    return Witness(equation, constants, "Rosati_2", A, B, C, D, tuple(sorted(ds)))


def search_eqmod1a(max_p: int, bound: int, prime_set: set[int]) -> Iterable[tuple[int, Witness]]:
    for B in range(1, bound + 1):
        for C in range(1, bound + 1):
            for D in range(1, bound + 1):
                modulus = 4 * B * C * D - 1
                inv_c = inv_mod(C, modulus)
                if inv_c is None:
                    continue
                residue = (-B * inv_c) % modulus
                for p in p_values_for_residue(max_p, residue, modulus, prime_set):
                    numerator = B + p * C
                    if numerator % modulus:
                        continue
                    A = numerator // modulus
                    witness = rosati1_witness(p, "eqmod1a", {"B": B, "C": C, "D": D}, A, B, C, D)
                    if witness:
                        yield p, witness


def search_eqmod1b(max_p: int, bound: int, prime_set: set[int]) -> Iterable[tuple[int, Witness]]:
    for A in range(1, bound + 1):
        for B in range(1, bound + 1):
            for E in range(1, bound + 1):
                if (A + B) % E:
                    continue
                C = (A + B) // E
                if C <= 0:
                    continue
                modulus = 4 * A * B
                residue = (-E) % modulus
                for p in p_values_for_residue(max_p, residue, modulus, prime_set):
                    D = (p + E) // modulus
                    witness = rosati1_witness(p, "eqmod1b", {"A": A, "B": B, "E": E}, A, B, C, D)
                    if witness:
                        yield p, witness


def search_eqmod1c(max_p: int, bound: int, prime_set: set[int]) -> Iterable[tuple[int, Witness]]:
    for B in range(1, bound + 1):
        for D in range(1, bound + 1):
            for E in range(1, bound + 1):
                modulus = 4 * B * D * E
                residue = (-E - 4 * B * B * D) % modulus
                for p in p_values_for_residue(max_p, residue, modulus, prime_set):
                    if (p + E) % (4 * B * D):
                        continue
                    A = (p + E) // (4 * B * D)
                    C = (p + E + 4 * B * B * D) // modulus
                    witness = rosati1_witness(p, "eqmod1c", {"B": B, "D": D, "E": E}, A, B, C, D)
                    if witness:
                        yield p, witness


def search_eqmod2a(max_p: int, bound: int, prime_set: set[int]) -> Iterable[tuple[int, Witness]]:
    for A in range(1, bound + 1):
        for B in range(1, bound + 1):
            for E in range(1, bound + 1):
                if (A + B) % E:
                    continue
                modulus = 4 * A * B
                inv_e = inv_mod(E, modulus)
                if inv_e is None:
                    continue
                residue = (-inv_e) % modulus
                C = (A + B) // E
                for p in p_values_for_residue(max_p, residue, modulus, prime_set):
                    D = (p * E + 1) // modulus
                    witness = rosati2_witness(p, "eqmod2a", {"A": A, "B": B, "E": E}, A, B, C, D)
                    if witness:
                        yield p, witness


def search_eqmod2b(max_p: int, bound: int, prime_set: set[int]) -> Iterable[tuple[int, Witness]]:
    for B in range(1, bound + 1):
        for C in range(1, bound + 1):
            for F in range(1, bound + 1):
                inv_b = inv_mod(B, F)
                if inv_b is None:
                    continue
                merged = crt_pair((-F) % (4 * B * C), 4 * B * C, (-C * inv_b) % F, F)
                if merged is None:
                    continue
                residue, modulus = merged
                for p in p_values_for_residue(max_p, residue, modulus, prime_set):
                    D = (p + F) // (4 * B * C)
                    A = (p * B + C) // F
                    witness = rosati2_witness(p, "eqmod2b", {"B": B, "C": C, "F": F}, A, B, C, D)
                    if witness:
                        yield p, witness


def search_eqmod2c(max_p: int, bound: int, prime_set: set[int]) -> Iterable[tuple[int, Witness]]:
    for B in range(1, bound + 1):
        for D in range(1, bound + 1):
            for F in range(1, bound + 1):
                if (4 * B * B * D + 1) % F:
                    continue
                modulus = 4 * B * D
                residue = (-F) % modulus
                E = (4 * B * B * D + 1) // F
                for p in p_values_for_residue(max_p, residue, modulus, prime_set):
                    C = (p + F) // modulus
                    A = C * E - B
                    witness = rosati2_witness(p, "eqmod2c", {"B": B, "D": D, "F": F}, A, B, C, D)
                    if witness:
                        yield p, witness


def search_eqmod2d(max_p: int, bound: int, prime_set: set[int]) -> Iterable[tuple[int, Witness]]:
    for C in range(1, bound + 1):
        for D in range(1, bound + 1):
            first_modulus = 4 * C * D
            for F in range(1, bound + 1):
                residue1 = (-F) % first_modulus
                for p in p_values_for_residue(max_p, residue1, first_modulus, prime_set):
                    if (p * p + 4 * C * C * D) % F:
                        continue
                    B = (p + F) // first_modulus
                    if (p * B + C) % F:
                        continue
                    A = (p * B + C) // F
                    witness = rosati2_witness(p, "eqmod2d", {"C": C, "D": D, "F": F}, A, B, C, D)
                    if witness:
                        yield p, witness


SEARCHERS: list[Callable[[int, int, set[int]], Iterable[tuple[int, Witness]]]] = [
    search_eqmod1a,
    search_eqmod1b,
    search_eqmod1c,
    search_eqmod2a,
    search_eqmod2b,
    search_eqmod2c,
    search_eqmod2d,
]


def witness_json(witness: Witness) -> dict[str, Any]:
    return {
        "equation": witness.equation,
        "constants": witness.constants,
        "rosati_case": witness.rosati_case,
        "A": witness.A,
        "B": witness.B,
        "C": witness.C,
        "D": witness.D,
        "denominators": list(witness.denominators),
    }


def local_basic_prime_covered(p: int) -> bool:
    return p % 3 == 2 or p % 4 == 3 or p % 8 == 5


def run(max_p: int, constant_bound: int, target_hard_strip: bool) -> dict[str, Any]:
    prime_list = [p for p in primes_up_to(max_p) if p > 2]
    prime_set = set(prime_list)
    target_primes = [p for p in prime_list if p % 24 == 1] if target_hard_strip else prime_list

    witnesses_by_prime: dict[int, list[Witness]] = defaultdict(list)
    first_witness_by_equation: dict[str, dict[str, Any]] = {}
    candidate_hits = 0

    for searcher in SEARCHERS:
        for p, witness in searcher(max_p, constant_bound, prime_set):
            if target_hard_strip and p % 24 != 1:
                continue
            candidate_hits += 1
            existing_keys = {
                (w.equation, tuple(sorted(w.constants.items())), w.denominators)
                for w in witnesses_by_prime[p]
            }
            key = (witness.equation, tuple(sorted(witness.constants.items())), witness.denominators)
            if key in existing_keys:
                continue
            witnesses_by_prime[p].append(witness)
            first_witness_by_equation.setdefault(witness.equation, {"p": p, **witness_json(witness)})

    covered_primes = sorted(p for p in target_primes if witnesses_by_prime.get(p))
    equation_counts = Counter()
    witness_count_by_prime: dict[str, int] = {}
    examples_by_prime: dict[str, Any] = {}

    for p in covered_primes:
        witness_count_by_prime[str(p)] = len(witnesses_by_prime[p])
        for witness in witnesses_by_prime[p]:
            equation_counts[witness.equation] += 1
        if len(examples_by_prime) < 20:
            examples_by_prime[str(p)] = [witness_json(w) for w in witnesses_by_prime[p][:5]]

    local_basic = [p for p in prime_list if local_basic_prime_covered(p)]
    hard_strip_uncovered = [p for p in target_primes if not witnesses_by_prime.get(p)]
    verdict = "SALEZ_GENERAL_BOUNDED_SEARCH_READY" if covered_primes else "NEEDS_LARGER_BOUND"

    return {
        "experiment_id": EXP_ID,
        "generated_at": now_iso(),
        "status": "REVIEW_ONLY",
        "verdict": verdict,
        "claim_ceiling": "bounded reconstruction/search over Salez's seven reference-equation families; not the optimized Salez sieve and not a proof",
        "max_p": max_p,
        "constant_bound": constant_bound,
        "target": "prime hard strip p = 1 mod 24" if target_hard_strip else "all odd primes",
        "prime_population": {
            "odd_prime_count": len(prime_list),
            "local_basic_prime_count": len(local_basic),
            "hard_strip_prime_count": len([p for p in prime_list if p % 24 == 1]),
            "target_prime_count": len(target_primes),
        },
        "bounded_salez_search": {
            "covered_target_prime_count": len(covered_primes),
            "coverage_fraction": len(covered_primes) / len(target_primes) if target_primes else 0.0,
            "candidate_hits_before_dedup": candidate_hits,
            "equation_counts": dict(equation_counts),
            "first_covered_primes": covered_primes[:60],
            "first_uncovered_target_primes": hard_strip_uncovered[:60],
            "witness_count_by_prime_sample": dict(list(witness_count_by_prime.items())[:40]),
            "examples_by_prime": examples_by_prime,
            "first_witness_by_equation": first_witness_by_equation,
        },
        "sota_boundary": {
            "salez_source": "https://arxiv.org/abs/1406.6307",
            "salez_boundary": "Salez reports a complete seven-equation framework and optimized checking to 10^17; this bounded search is a reconstruction audit only.",
            "current_computational_frontier_source": "https://arxiv.org/abs/2509.00128",
            "current_computational_frontier_boundary": "Mihnea-Bogdan report verification to 10^18; this audit does not compete with that bound.",
        },
        "forbidden_actions": [
            "do not claim proof of Erdos-Straus",
            "do not claim global computational SOTA",
            "do not treat bounded constant search as Salez's optimized sieve",
            "do not update D1/proof registry without review",
        ],
    }


def report(results: dict[str, Any]) -> str:
    pop = results["prime_population"]
    search = results["bounded_salez_search"]
    lines = [
        "# Erdos #242 Salez General-Equation Bounded Search",
        "",
        f"Experiment: `{results['experiment_id']}`",
        "Status: `REVIEW_ONLY`",
        f"Verdict: `{results['verdict']}`",
        "",
        "## Meaning",
        "",
        "This audit reconstructs Salez's seven reference equations as bounded constant-parameter search families. For each hit it rebuilds Rosati variables `A,B,C,D`, reconstructs the denominators, and verifies the identity `4/p = 1/x + 1/y + 1/z` directly.",
        "",
        "It is a stronger alignment than the Example 1 check, but it is not Salez's optimized sieve and it does not move the computational frontier.",
        "",
        "## Run Shape",
        "",
        f"- Max prime checked: `{results['max_p']}`",
        f"- Constant bound: `{results['constant_bound']}`",
        f"- Target: `{results['target']}`",
        f"- Odd primes: `{pop['odd_prime_count']}`",
        f"- Basic local-family primes: `{pop['local_basic_prime_count']}`",
        f"- Hard-strip primes `p = 1 mod 24`: `{pop['hard_strip_prime_count']}`",
        "",
        "## Bounded Search Result",
        "",
        f"- Covered target primes: `{search['covered_target_prime_count']} / {pop['target_prime_count']} = {search['coverage_fraction']:.6f}`",
        f"- Equation hit counts: `{search['equation_counts']}`",
        f"- First covered target primes: `{search['first_covered_primes'][:20]}`",
        f"- First uncovered target primes: `{search['first_uncovered_target_primes'][:20]}`",
        "",
        "## Boundary",
        "",
        "This makes the #242 lane genuinely Salez-facing: the seven equations are now executable search families, not just labels in a report. The remaining gap is to match Salez's optimized sieve behavior and parameter strategy, then decide whether our operator view simplifies or only repackages that known machinery.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-p", type=int, default=100000)
    parser.add_argument("--constant-bound", type=int, default=80)
    parser.add_argument("--all-odd-primes", action="store_true", help="Search all odd primes instead of only p = 1 mod 24 target strip.")
    args = parser.parse_args()

    results = run(args.max_p, args.constant_bound, target_hard_strip=not args.all_odd_primes)
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
        "coverage_fraction": results["bounded_salez_search"]["coverage_fraction"],
        "covered_target_prime_count": results["bounded_salez_search"]["covered_target_prime_count"],
        "target_prime_count": results["prime_population"]["target_prime_count"],
        "equation_counts": results["bounded_salez_search"]["equation_counts"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
