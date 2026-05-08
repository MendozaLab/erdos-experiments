#!/usr/bin/env python3
"""Reusable Salez-style filter certificates for Erdos #242.

This module wraps the seven-equation reconstruction in
``salez_general_equation_search.py`` and gives each verified witness a
Salez-style filter identity: a natural modulus, residue, equation label, Rosati
variables, and denominators. The filters here are certificate-bearing local
filters, not Salez's optimized high-scale sieve.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Iterable

import salez_general_equation_search as salez


@dataclass(frozen=True)
class FilterCertificate:
    p: int
    filter_modulus: int
    filter_residue: int
    equation_id: str
    constants: dict[str, int]
    rosati_case: str
    A: int
    B: int
    C: int
    D: int
    x: int
    y: int
    z: int
    verified: bool

    def key(self) -> tuple[Any, ...]:
        return (
            self.p,
            self.filter_modulus,
            self.filter_residue,
            self.equation_id,
            tuple(sorted(self.constants.items())),
            (self.x, self.y, self.z),
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "p": self.p,
            "filter_modulus": self.filter_modulus,
            "filter_residue": self.filter_residue,
            "equation_id": self.equation_id,
            "constants": self.constants,
            "rosati_case": self.rosati_case,
            "A": self.A,
            "B": self.B,
            "C": self.C,
            "D": self.D,
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "verified": self.verified,
        }


def target_primes(max_n: int, hard_strip_only: bool = True) -> list[int]:
    primes = [p for p in salez.primes_up_to(max_n) if p > 2]
    if hard_strip_only:
        return [p for p in primes if p % 24 == 1]
    return primes


def witness_filter_modulus_residue(witness: salez.Witness) -> tuple[int, int]:
    c = witness.constants
    if witness.equation == "eqmod1a":
        modulus = 4 * c["B"] * c["C"] * c["D"] - 1
        residue = (-c["B"] * pow(c["C"], -1, modulus)) % modulus
        return modulus, residue
    if witness.equation == "eqmod1b":
        modulus = 4 * c["A"] * c["B"]
        return modulus, (-c["E"]) % modulus
    if witness.equation == "eqmod1c":
        modulus = 4 * c["B"] * c["D"] * c["E"]
        return modulus, (-c["E"] - 4 * c["B"] * c["B"] * c["D"]) % modulus
    if witness.equation == "eqmod2a":
        modulus = 4 * c["A"] * c["B"]
        return modulus, (-pow(c["E"], -1, modulus)) % modulus
    if witness.equation == "eqmod2b":
        merged = salez.crt_pair(
            (-c["F"]) % (4 * c["B"] * c["C"]),
            4 * c["B"] * c["C"],
            (-c["C"] * pow(c["B"], -1, c["F"])) % c["F"],
            c["F"],
        )
        if merged is None:
            raise ValueError(f"invalid eqmod2b constants: {c}")
        residue, modulus = merged
        return modulus, residue
    if witness.equation == "eqmod2c":
        modulus = 4 * c["B"] * c["D"]
        return modulus, (-c["F"]) % modulus
    if witness.equation == "eqmod2d":
        # The first congruence is p + F = 0 mod 4CD; the quadratic condition
        # p^2 + 4C^2D = 0 mod F is already checked by the witness builder.
        modulus = 4 * c["C"] * c["D"]
        return modulus, (-c["F"]) % modulus
    raise KeyError(witness.equation)


def to_filter_certificate(p: int, witness: salez.Witness) -> FilterCertificate:
    modulus, residue = witness_filter_modulus_residue(witness)
    x, y, z = witness.denominators
    verified = (
        p % modulus == residue
        and salez.verify_solution(p, witness.denominators)
    )
    return FilterCertificate(
        p=p,
        filter_modulus=modulus,
        filter_residue=residue,
        equation_id=witness.equation,
        constants=dict(witness.constants),
        rosati_case=witness.rosati_case,
        A=witness.A,
        B=witness.B,
        C=witness.C,
        D=witness.D,
        x=x,
        y=y,
        z=z,
        verified=verified,
    )


def generate_filter_certificates(
    max_n: int,
    constant_bound: int,
    hard_strip_only: bool = True,
    allowed_filter_moduli: set[int] | None = None,
) -> list[FilterCertificate]:
    prime_set = set(target_primes(max_n, hard_strip_only=False))
    certificates: list[FilterCertificate] = []
    seen: set[tuple[Any, ...]] = set()

    for searcher in salez.SEARCHERS:
        for p, witness in searcher(max_n, constant_bound, prime_set):
            if hard_strip_only and p % 24 != 1:
                continue
            cert = to_filter_certificate(p, witness)
            if allowed_filter_moduli is not None and cert.filter_modulus not in allowed_filter_moduli:
                continue
            key = cert.key()
            if key in seen:
                continue
            seen.add(key)
            certificates.append(cert)

    certificates.sort(key=lambda c: (c.p, c.filter_modulus, c.equation_id, c.x, c.y, c.z))
    return certificates


def choose_one_certificate_per_prime(certificates: Iterable[FilterCertificate]) -> dict[int, FilterCertificate]:
    chosen: dict[int, FilterCertificate] = {}
    for cert in certificates:
        current = chosen.get(cert.p)
        if current is None:
            chosen[cert.p] = cert
            continue
        current_rank = (current.filter_modulus, current.equation_id, current.x + current.y + current.z)
        cert_rank = (cert.filter_modulus, cert.equation_id, cert.x + cert.y + cert.z)
        if cert_rank < current_rank:
            chosen[cert.p] = cert
    return chosen


def summarize_certificates(
    certificates: list[FilterCertificate],
    max_n: int,
    hard_strip_only: bool = True,
) -> dict[str, Any]:
    targets = target_primes(max_n, hard_strip_only)
    chosen = choose_one_certificate_per_prime(certificates)
    certified_targets = [p for p in targets if p in chosen]
    missing_targets = [p for p in targets if p not in chosen]
    invalid = [cert for cert in certificates if not cert.verified]

    equation_counts: defaultdict[str, int] = defaultdict(int)
    filter_residue_counts: defaultdict[int, set[int]] = defaultdict(set)
    for cert in certificates:
        equation_counts[cert.equation_id] += 1
        filter_residue_counts[cert.filter_modulus].add(cert.filter_residue)

    filter_moduli = sorted(filter_residue_counts)
    return {
        "target_prime_count": len(targets),
        "certified_target_count": len(certified_targets),
        "coverage_fraction": len(certified_targets) / len(targets) if targets else 0.0,
        "certificate_count": len(certificates),
        "invalid_certificate_count": len(invalid),
        "filter_only_no_witness_count": 0,
        "first_certified_targets": certified_targets[:80],
        "first_missing_targets": missing_targets[:80],
        "equation_counts": dict(sorted(equation_counts.items())),
        "filter_moduli_used_count": len(filter_moduli),
        "filter_moduli_used": filter_moduli,
        "filter_residue_class_count": sum(len(v) for v in filter_residue_counts.values()),
        "filter_residue_class_sample": {
            str(m): sorted(filter_residue_counts[m])[:40]
            for m in filter_moduli[:80]
        },
        "one_witness_per_certified_target": [
            chosen[p].to_json() for p in certified_targets
        ],
        "first_invalid_certificates": [cert.to_json() for cert in invalid[:20]],
    }


def verdict_for_summary(summary: dict[str, Any]) -> str:
    if summary["certified_target_count"] < summary["target_prime_count"]:
        return "NEEDS_FILTER_COVERAGE"
    if summary["invalid_certificate_count"] or summary["filter_only_no_witness_count"]:
        return "FILTER_PARITY_WITNESS_GAP"
    return "CERTIFIED_PARITY_LOCAL"
