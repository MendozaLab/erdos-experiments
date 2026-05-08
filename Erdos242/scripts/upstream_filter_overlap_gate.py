#!/usr/bin/env python3
"""Exact upstream filter-table overlap gate for Erdos #242.

This runner compares the local certificate-bearing filter surface against the
pinned Mihnea/Bogdan `Filters.txt` table. It measures overlap and semantic
differences. It does not vendor upstream code/data, execute the upstream
checker, claim exact equality, or claim 10^18 replication.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import salez_filter_certificates as certs


EXP_ID = "EXP-MATH-ERDOS242-UPSTREAM-FILTER-OVERLAP-20260508-01"
ROOT = Path(__file__).resolve().parents[3]
RESULT_DIR = ROOT / "erdos-experiments" / "results" / "erdos-242"
PINNED_COMMIT = "e36eef1815d339701b9f168fe7fa504ccfa401e8"
UPSTREAM_BASE = f"https://raw.githubusercontent.com/esc-paper/erdos-straus/{PINNED_COMMIT}/section1/resources"
FILTERS_URL = f"{UPSTREAM_BASE}/Filters.txt"
RESIDUES_URL = f"{UPSTREAM_BASE}/Residues.txt"


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


def open_url(url: str):
    request = urllib.request.Request(url, headers={"User-Agent": "h2-erdos242-overlap-gate"})
    return urllib.request.urlopen(request, timeout=120)


def parse_filter_line(line: str, line_number: int) -> tuple[int, set[int]]:
    raw_tokens = line.split()
    if not raw_tokens:
        raise ValueError(f"empty filter row at line {line_number}")
    try:
        tokens = [int(token) for token in raw_tokens]
    except ValueError as exc:
        raise ValueError(f"non-integer token at line {line_number}") from exc
    if tokens[-1] != -1:
        raise ValueError(f"filter row missing -1 terminator at line {line_number}")
    if len(tokens) < 3:
        raise ValueError(f"filter row has no residues at line {line_number}")
    modulus = tokens[0]
    if modulus <= 0:
        raise ValueError(f"non-positive modulus at line {line_number}")
    residues = set()
    for residue in tokens[1:-1]:
        if residue < 0:
            raise ValueError(f"negative residue before terminator at line {line_number}")
        residues.add(residue % modulus)
    return modulus, residues


def fetch_upstream_filters(url: str) -> dict[str, Any]:
    started = time.perf_counter()
    digest = hashlib.sha256()
    residues_by_modulus: defaultdict[int, set[int]] = defaultdict(set)
    row_count = 0
    parse_errors: list[str] = []
    first_rows: list[dict[str, Any]] = []
    byte_count = 0

    with open_url(url) as response:
        for row_count, raw in enumerate(response, start=1):
            digest.update(raw)
            byte_count += len(raw)
            line = raw.decode("utf-8", errors="replace").strip()
            if not line:
                continue
            try:
                modulus, residues = parse_filter_line(line, row_count)
            except ValueError as exc:
                parse_errors.append(str(exc))
                if len(parse_errors) >= 20:
                    break
                continue
            residues_by_modulus[modulus].update(residues)
            if len(first_rows) < 8:
                first_rows.append({
                    "line": row_count,
                    "modulus": modulus,
                    "residue_count": len(residues),
                    "residue_sample": sorted(residues)[:20],
                })

    filter_moduli = sorted(residues_by_modulus)
    pair_count = sum(len(v) for v in residues_by_modulus.values())
    return {
        "url": url,
        "fetch_status": "FETCHED",
        "bytes": byte_count,
        "sha256": digest.hexdigest(),
        "row_count": row_count,
        "parse_error_count": len(parse_errors),
        "parse_errors": parse_errors,
        "filter_modulus_count": len(filter_moduli),
        "filter_pair_count": pair_count,
        "first_rows": first_rows,
        "moduli_sample": filter_moduli[:80],
        "moduli_tail": filter_moduli[-20:],
        "residues_by_modulus": residues_by_modulus,
        "runtime_seconds": round(time.perf_counter() - started, 6),
    }


def fetch_residue_reference(url: str, sample_size: int = 40) -> dict[str, Any]:
    started = time.perf_counter()
    with open_url(url) as response:
        data = response.read()
    digest = hashlib.sha256(data).hexdigest()
    text = data.decode("utf-8", errors="replace")
    values = [int(token) for token in text.split()]
    return {
        "url": url,
        "fetch_status": "FETCHED",
        "bytes": len(data),
        "sha256": digest,
        "residue_count": len(values),
        "first_residues": values[:sample_size],
        "last_residues": values[-sample_size:],
        "runtime_seconds": round(time.perf_counter() - started, 6),
        "semantic_note": "Residues.txt is treated as a large-gap residual search surface, not as a certificate-witness table.",
    }


def local_certificate_surface(max_n: int, constant_bound: int) -> dict[str, Any]:
    started = time.perf_counter()
    certificates = certs.generate_filter_certificates(
        max_n=max_n,
        constant_bound=constant_bound,
        hard_strip_only=True,
    )
    invalid = [c for c in certificates if not c.verified or c.p % c.filter_modulus != c.filter_residue]
    pairs_by_modulus: defaultdict[int, set[int]] = defaultdict(set)
    witness_by_pair: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for cert in certificates:
        pair = (cert.filter_modulus, cert.filter_residue)
        pairs_by_modulus[cert.filter_modulus].add(cert.filter_residue)
        if len(witness_by_pair[pair]) < 5:
            witness_by_pair[pair].append(cert.to_json())

    filter_moduli = sorted(pairs_by_modulus)
    pair_count = sum(len(v) for v in pairs_by_modulus.values())
    target_primes = certs.target_primes(max_n, hard_strip_only=True)
    certified_primes = sorted({c.p for c in certificates})
    return {
        "max_n": max_n,
        "constant_bound": constant_bound,
        "certificate_count": len(certificates),
        "target_prime_count": len(target_primes),
        "certified_prime_count": len([p for p in target_primes if p in set(certified_primes)]),
        "invalid_certificate_count": len(invalid),
        "filter_modulus_count": len(filter_moduli),
        "filter_pair_count": pair_count,
        "filter_moduli_sample": filter_moduli[:80],
        "filter_moduli_tail": filter_moduli[-20:],
        "first_invalid_certificates": [c.to_json() for c in invalid[:20]],
        "pairs_by_modulus": pairs_by_modulus,
        "witness_by_pair": witness_by_pair,
        "runtime_seconds": round(time.perf_counter() - started, 6),
    }


def pair_iter(pairs_by_modulus: dict[int, set[int]]) -> Iterable[tuple[int, int]]:
    for modulus, residues in pairs_by_modulus.items():
        for residue in residues:
            yield modulus, residue


def compare_surfaces(upstream: dict[str, Any], local: dict[str, Any], sample_size: int) -> dict[str, Any]:
    upstream_pairs = set(pair_iter(upstream["residues_by_modulus"]))
    local_pairs = set(pair_iter(local["pairs_by_modulus"]))
    intersection = upstream_pairs & local_pairs
    upstream_only = sorted(upstream_pairs - local_pairs)[:sample_size]
    local_only = sorted(local_pairs - upstream_pairs)[:sample_size]

    witness_samples = []
    for pair in sorted(intersection):
        for witness in local["witness_by_pair"].get(pair, []):
            witness_samples.append(witness)
            if len(witness_samples) >= sample_size:
                break
        if len(witness_samples) >= sample_size:
            break

    return {
        "upstream_pair_count": len(upstream_pairs),
        "local_pair_count": len(local_pairs),
        "intersection_pair_count": len(intersection),
        "intersection_fraction_of_upstream": len(intersection) / len(upstream_pairs) if upstream_pairs else 0.0,
        "intersection_fraction_of_local": len(intersection) / len(local_pairs) if local_pairs else 0.0,
        "upstream_only_pair_sample": [{"modulus": m, "residue": r} for m, r in upstream_only],
        "local_only_pair_sample": [{"modulus": m, "residue": r} for m, r in local_only],
        "local_witnesses_with_upstream_pair_sample": witness_samples,
        "local_only_classification": "certificate-enriched natural moduli unless later proven equivalent to upstream filters",
        "semantic_note": "Exact equality is not expected because upstream Filters.txt is an optimized sieve table and local pairs are natural certificate-bearing filters derived from verified witnesses.",
    }


def compact_upstream_filters(upstream: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in upstream.items() if k != "residues_by_modulus"}


def compact_local_surface(local: dict[str, Any]) -> dict[str, Any]:
    return {
        k: v
        for k, v in local.items()
        if k not in {"pairs_by_modulus", "witness_by_pair"}
    }


def verdict_for(upstream: dict[str, Any], local: dict[str, Any], comparison: dict[str, Any]) -> str:
    if upstream.get("fetch_status") != "FETCHED" or upstream.get("parse_error_count"):
        return "UPSTREAM_TABLE_UNAVAILABLE"
    if local["invalid_certificate_count"] or not comparison:
        return "PARTIAL_OVERLAP_NEEDS_REVIEW"
    return "UPSTREAM_FILTER_OVERLAP_QUANTIFIED"


def run(args: argparse.Namespace) -> dict[str, Any]:
    try:
        upstream = fetch_upstream_filters(FILTERS_URL)
        residues = fetch_residue_reference(RESIDUES_URL, sample_size=args.sample_size)
    except Exception as exc:
        return {
            "experiment_id": EXP_ID,
            "generated_at": now_iso(),
            "status": "REVIEW_ONLY",
            "verdict": "UPSTREAM_TABLE_UNAVAILABLE",
            "claim_ceiling": "upstream-table availability failure; no SOTA or proof claim",
            "source_manifest": source_manifest(),
            "error": repr(exc),
            "forbidden_actions": forbidden_actions(),
        }

    local = local_certificate_surface(args.max_n, args.constant_bound)
    comparison = compare_surfaces(upstream, local, args.sample_size)
    verdict = verdict_for(upstream, local, comparison)
    return {
        "experiment_id": EXP_ID,
        "generated_at": now_iso(),
        "status": "REVIEW_ONLY",
        "verdict": verdict,
        "claim_ceiling": "exact upstream filter-table overlap audit only; not exact pipeline equality, not 10^18 replication, not SOTA, and not a proof",
        "source_manifest": source_manifest(),
        "upstream_filters": compact_upstream_filters(upstream),
        "upstream_residues_reference": residues,
        "local_certificate_surface": compact_local_surface(local),
        "overlap": comparison,
        "forbidden_actions": forbidden_actions(),
    }


def source_manifest() -> dict[str, Any]:
    return {
        "status": "REFERENCE_INPUTS_ONLY",
        "third_party_code_committed": False,
        "mihnea_bogdan_github": "https://github.com/esc-paper/erdos-straus",
        "mihnea_bogdan_pinned_commit": PINNED_COMMIT,
        "filters_url": FILTERS_URL,
        "residues_url": RESIDUES_URL,
        "local_implementation_status": "RECONSTRUCTED_FROM_SOURCES",
    }


def forbidden_actions() -> list[str]:
    return [
        "do not claim proof of Erdos-Straus",
        "do not claim exact equality to Mihnea/Bogdan's optimized filter pipeline",
        "do not claim replication of the 10^18 result",
        "do not claim global computational SOTA",
        "do not update D1/proof registry from this review-only artifact",
    ]


def report(results: dict[str, Any]) -> str:
    if results["verdict"] == "UPSTREAM_TABLE_UNAVAILABLE":
        return "\n".join([
            "# Erdos #242 Upstream Filter-Table Overlap Gate",
            "",
            f"Experiment: `{results['experiment_id']}`",
            "Status: `REVIEW_ONLY`",
            f"Verdict: `{results['verdict']}`",
            "",
            f"Error: `{results.get('error')}`",
            "",
        ])

    upstream = results["upstream_filters"]
    residues = results["upstream_residues_reference"]
    local = results["local_certificate_surface"]
    overlap = results["overlap"]
    return "\n".join([
        "# Erdos #242 Upstream Filter-Table Overlap Gate",
        "",
        f"Experiment: `{results['experiment_id']}`",
        "Status: `REVIEW_ONLY`",
        f"Verdict: `{results['verdict']}`",
        "",
        "## Meaning",
        "",
        "This gate compares our certificate-bearing filter surface with Mihnea/Bogdan's exact pinned `Filters.txt` table. It quantifies overlap and non-overlap; it does not assert exact equality or computational-frontier replication.",
        "",
        "## Upstream Table",
        "",
        f"- Pinned commit: `{PINNED_COMMIT}`",
        f"- Filter rows: `{upstream['row_count']}`",
        f"- Upstream moduli: `{upstream['filter_modulus_count']}`",
        f"- Upstream `(modulus,residue)` pairs: `{upstream['filter_pair_count']}`",
        f"- Filters SHA-256: `{upstream['sha256']}`",
        f"- Residues SHA-256: `{residues['sha256']}`",
        f"- Residues count: `{residues['residue_count']}`",
        "",
        "## Local Surface",
        "",
        f"- Max n: `{local['max_n']}`",
        f"- Constant bound: `{local['constant_bound']}`",
        f"- Certificates: `{local['certificate_count']}`",
        f"- Local moduli: `{local['filter_modulus_count']}`",
        f"- Local `(modulus,residue)` pairs: `{local['filter_pair_count']}`",
        f"- Invalid certificates: `{local['invalid_certificate_count']}`",
        "",
        "## Overlap",
        "",
        f"- Pair intersection: `{overlap['intersection_pair_count']}`",
        f"- Fraction of upstream pairs: `{overlap['intersection_fraction_of_upstream']:.8f}`",
        f"- Fraction of local pairs: `{overlap['intersection_fraction_of_local']:.8f}`",
        f"- Upstream-only pair sample: `{overlap['upstream_only_pair_sample'][:12]}`",
        f"- Local-only pair sample: `{overlap['local_only_pair_sample'][:12]}`",
        "",
        "## Boundary",
        "",
        "Local-only pairs are classified as certificate-enriched natural moduli unless later proven equivalent to upstream filters. Upstream `Residues.txt` is treated as a residual progression/search surface, not a witness table.",
        "",
    ])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=1_000_000)
    parser.add_argument("--constant-bound", type=int, default=21)
    parser.add_argument("--sample-size", type=int, default=40)
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
        "overlap": results.get("overlap", {}),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
