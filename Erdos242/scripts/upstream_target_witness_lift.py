#!/usr/bin/env python3
"""Target-level witness lift over pinned Mihnea/Bogdan filters for Erdos #242.

This runner asks a target-level question that pair overlap cannot answer:
when the pinned upstream `Filters.txt` table certifies a hard-strip prime
`p <= 1_000_000`, can the local certificate layer attach a concrete verified
Rosati witness and denominators to that same prime?

It is REVIEW_ONLY. It does not vendor upstream code/data, execute the upstream
pipeline, claim exact equality, claim 10^18 replication, or prove
Erdos-Straus.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import salez_filter_certificates as certs
import upstream_filter_overlap_gate as overlap_gate


EXP_ID = "EXP-MATH-ERDOS242-UPSTREAM-TARGET-WITNESS-LIFT-20260508-01"
ROOT = Path(__file__).resolve().parents[3]
RESULT_DIR = ROOT / "erdos-experiments" / "results" / "erdos-242"
EXPECTED_FILTERS_SHA256 = "4aafb1e2cd909439be5a32c84bff5d6cd27fdce95a91edaee219de237f2c0ffd"


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


def local_witness_surface(max_n: int, constant_bound: int) -> dict[str, Any]:
    started = time.perf_counter()
    certificates = certs.generate_filter_certificates(
        max_n=max_n,
        constant_bound=constant_bound,
        hard_strip_only=True,
    )
    chosen = certs.choose_one_certificate_per_prime(certificates)
    targets = certs.target_primes(max_n, hard_strip_only=True)
    invalid = [
        c
        for c in certificates
        if not c.verified or c.p % c.filter_modulus != c.filter_residue
    ]
    pairs: set[tuple[int, int]] = {
        (c.filter_modulus, c.filter_residue)
        for c in certificates
    }
    return {
        "max_n": max_n,
        "constant_bound": constant_bound,
        "target_primes": targets,
        "target_prime_count": len(targets),
        "certificates": certificates,
        "chosen_by_prime": chosen,
        "witnessed_targets": sorted(chosen),
        "locally_witnessed_target_count": len(chosen),
        "certificate_count": len(certificates),
        "invalid_certificate_count": len(invalid),
        "invalid_certificate_sample": [c.to_json() for c in invalid[:20]],
        "filter_modulus_count": len({c.filter_modulus for c in certificates}),
        "filter_pair_count": len(pairs),
        "runtime_seconds": round(time.perf_counter() - started, 6),
    }


def first_candidate_at_or_above(residue: int, modulus: int, lower: int) -> int:
    if residue >= lower:
        return residue
    return residue + ((lower - residue + modulus - 1) // modulus) * modulus


def scan_upstream_target_hits(
    *,
    max_n: int,
    target_set: set[int],
    sample_per_target: int,
    global_sample_size: int,
) -> dict[str, Any]:
    started = time.perf_counter()
    digest = hashlib.sha256()
    row_count = 0
    byte_count = 0
    parse_errors: list[str] = []
    moduli_seen: set[int] = set()
    pair_count_by_modulus: defaultdict[int, set[int]] = defaultdict(set)
    hit_counts_by_target: defaultdict[int, int] = defaultdict(int)
    hit_samples_by_target: dict[int, list[dict[str, int]]] = defaultdict(list)
    global_hit_samples: list[dict[str, int]] = []
    hit_validation_error_count = 0

    with overlap_gate.open_url(overlap_gate.FILTERS_URL) as response:
        for row_count, raw in enumerate(response, start=1):
            digest.update(raw)
            byte_count += len(raw)
            line = raw.decode("utf-8", errors="replace").strip()
            if not line:
                continue
            try:
                modulus, residues = overlap_gate.parse_filter_line(line, row_count)
            except ValueError as exc:
                parse_errors.append(str(exc))
                if len(parse_errors) >= 20:
                    break
                continue

            moduli_seen.add(modulus)
            pair_count_by_modulus[modulus].update(residues)
            for residue in residues:
                p = first_candidate_at_or_above(residue, modulus, 3)
                while p <= max_n:
                    if p in target_set:
                        if p % modulus != residue:
                            hit_validation_error_count += 1
                        hit_counts_by_target[p] += 1
                        sample = {
                            "p": p,
                            "line": row_count,
                            "filter_modulus": modulus,
                            "filter_residue": residue,
                        }
                        if len(hit_samples_by_target[p]) < sample_per_target:
                            hit_samples_by_target[p].append(sample)
                        if len(global_hit_samples) < global_sample_size:
                            global_hit_samples.append(sample)
                    p += modulus

    upstream_filtered_targets = sorted(hit_counts_by_target)
    filter_pair_count = sum(len(residues) for residues in pair_count_by_modulus.values())
    return {
        "url": overlap_gate.FILTERS_URL,
        "fetch_status": "FETCHED",
        "bytes": byte_count,
        "sha256": digest.hexdigest(),
        "sha256_matches_prior_overlap": digest.hexdigest() == EXPECTED_FILTERS_SHA256,
        "expected_prior_overlap_sha256": EXPECTED_FILTERS_SHA256,
        "row_count": row_count,
        "parse_error_count": len(parse_errors),
        "parse_errors": parse_errors,
        "filter_modulus_count": len(moduli_seen),
        "filter_pair_count": filter_pair_count,
        "upstream_filtered_targets": upstream_filtered_targets,
        "upstream_filtered_target_count": len(upstream_filtered_targets),
        "total_target_hit_records": sum(hit_counts_by_target.values()),
        "hit_counts_by_target": dict(sorted(hit_counts_by_target.items())),
        "hit_samples_by_target": {
            str(p): hit_samples_by_target[p]
            for p in upstream_filtered_targets
        },
        "global_hit_samples": global_hit_samples,
        "hit_validation_error_count": hit_validation_error_count,
        "parser_contract": {
            "row_shape": "modulus residues... -1",
            "malformed_rows_rejected": True,
            "repeated_moduli_union_residues": True,
        },
        "runtime_seconds": round(time.perf_counter() - started, 6),
    }


def build_target_records(
    *,
    upstream: dict[str, Any],
    local: dict[str, Any],
    sample_missing: int,
) -> dict[str, Any]:
    upstream_targets = set(upstream["upstream_filtered_targets"])
    local_targets = set(local["witnessed_targets"])
    lifted_targets = sorted(upstream_targets & local_targets)
    upstream_missing_local = sorted(upstream_targets - local_targets)
    local_without_upstream = sorted(local_targets - upstream_targets)

    chosen = local["chosen_by_prime"]
    hit_counts = upstream["hit_counts_by_target"]
    hit_samples = upstream["hit_samples_by_target"]

    target_records: list[dict[str, Any]] = []
    emitted_witness_invalid_count = 0
    for p in sorted(upstream_targets):
        witness = chosen.get(p)
        witness_json = witness.to_json() if witness is not None else None
        if witness is not None and not witness.verified:
            emitted_witness_invalid_count += 1
        target_records.append({
            "p": p,
            "upstream_hit_count": hit_counts[p],
            "upstream_hit_samples": hit_samples.get(str(p), []),
            "locally_witnessed": witness is not None,
            "local_witness": witness_json,
        })

    sample_lifted = []
    for p in lifted_targets[:sample_missing]:
        sample_lifted.append({
            "p": p,
            "upstream_hit_count": hit_counts[p],
            "upstream_hit_samples": hit_samples.get(str(p), []),
            "local_witness": chosen[p].to_json(),
        })

    return {
        "upstream_filter_target_count": len(upstream_targets),
        "locally_witnessed_target_count": len(local_targets),
        "upstream_filtered_and_locally_witnessed_count": len(lifted_targets),
        "upstream_filtered_local_witness_missing_count": len(upstream_missing_local),
        "local_witness_without_observed_upstream_hit_count": len(local_without_upstream),
        "emitted_witness_invalid_count": emitted_witness_invalid_count,
        "upstream_filtered_local_witness_missing_sample": upstream_missing_local[:sample_missing],
        "local_witness_without_observed_upstream_hit_sample": local_without_upstream[:sample_missing],
        "lifted_target_sample": sample_lifted,
        "target_records": target_records,
        "semantic_note": "The primary score is target lift, not exact (modulus,residue) equality.",
    }


def pair_overlap_support(upstream_scan: dict[str, Any], local: dict[str, Any]) -> dict[str, Any]:
    overlap_path = (
        RESULT_DIR
        / "EXP-MATH-ERDOS242-UPSTREAM-FILTER-OVERLAP-20260508-01_RESULTS.json"
    )
    support: dict[str, Any] = {
        "role": "supporting_context_only",
        "semantic_note": "Pair overlap is not the main score for this gate.",
    }
    if overlap_path.exists():
        previous = json.loads(overlap_path.read_text(encoding="utf-8"))
        overlap = previous.get("overlap", {})
        support.update({
            "source": str(overlap_path.relative_to(ROOT)),
            "previous_pair_intersection_count": overlap.get("intersection_pair_count"),
            "previous_intersection_fraction_of_upstream": overlap.get("intersection_fraction_of_upstream"),
            "previous_intersection_fraction_of_local": overlap.get("intersection_fraction_of_local"),
            "previous_upstream_pair_count": overlap.get("upstream_pair_count"),
            "previous_local_pair_count": overlap.get("local_pair_count"),
        })
    else:
        support.update({
            "source": "not_found",
            "current_upstream_pair_count": upstream_scan["filter_pair_count"],
            "current_local_pair_count": local["filter_pair_count"],
        })
    return support


def verdict_for(upstream: dict[str, Any], local: dict[str, Any], lift: dict[str, Any]) -> str:
    if upstream.get("fetch_status") != "FETCHED":
        return "UPSTREAM_FILTER_INPUT_UNAVAILABLE"
    if (
        upstream["parse_error_count"]
        or upstream["hit_validation_error_count"]
        or local["invalid_certificate_count"]
        or lift["emitted_witness_invalid_count"]
        or lift["upstream_filtered_local_witness_missing_count"]
    ):
        return "PARTIAL_TARGET_LIFT_NEEDS_REVIEW"
    return "UPSTREAM_TARGETS_WITNESS_LIFTED"


def source_manifest() -> dict[str, Any]:
    return {
        "status": "REFERENCE_INPUTS_ONLY",
        "third_party_code_committed": False,
        "mihnea_bogdan_github": "https://github.com/esc-paper/erdos-straus",
        "mihnea_bogdan_pinned_commit": overlap_gate.PINNED_COMMIT,
        "filters_url": overlap_gate.FILTERS_URL,
        "local_implementation_status": "RECONSTRUCTED_FROM_SOURCES",
        "witness_layer_status": "LOCAL_CERTIFICATE_ENRICHMENT",
    }


def forbidden_actions() -> list[str]:
    return [
        "do not claim proof of Erdos-Straus",
        "do not claim exact equality to Mihnea/Bogdan's optimized filter pipeline",
        "do not claim replication of the 10^18 result",
        "do not claim global computational SOTA",
        "do not update D1/proof registry from this review-only artifact",
    ]


def run(args: argparse.Namespace) -> dict[str, Any]:
    try:
        local = local_witness_surface(args.max_n, args.constant_bound)
        target_set = set(local["target_primes"])
        upstream = scan_upstream_target_hits(
            max_n=args.max_n,
            target_set=target_set,
            sample_per_target=args.sample_per_target,
            global_sample_size=args.sample_size,
        )
    except Exception as exc:
        return {
            "experiment_id": EXP_ID,
            "generated_at": now_iso(),
            "status": "REVIEW_ONLY",
            "verdict": "UPSTREAM_FILTER_INPUT_UNAVAILABLE",
            "claim_ceiling": "upstream target-lift input failure; no SOTA, exact-pipeline, 10^18, or proof claim",
            "source_manifest": source_manifest(),
            "error": repr(exc),
            "forbidden_actions": forbidden_actions(),
        }

    lift = build_target_records(
        upstream=upstream,
        local=local,
        sample_missing=args.sample_size,
    )
    verdict = verdict_for(upstream, local, lift)
    local_summary = {
        k: v
        for k, v in local.items()
        if k not in {"certificates", "chosen_by_prime", "target_primes"}
    }
    return {
        "experiment_id": EXP_ID,
        "generated_at": now_iso(),
        "status": "REVIEW_ONLY",
        "verdict": verdict,
        "claim_ceiling": "target-level upstream-filter witness lift only; not exact pipeline equality, not 10^18 replication, not SOTA, and not a proof",
        "comparison_target": {
            "problem": "Erdos #242 / Erdos-Straus",
            "target_strip": "prime p = 1 mod 24",
            "max_n": args.max_n,
            "constant_bound": args.constant_bound,
        },
        "source_manifest": source_manifest(),
        "upstream_target_scan": {
            k: v
            for k, v in upstream.items()
            if k not in {"upstream_filtered_targets", "hit_counts_by_target", "hit_samples_by_target"}
        },
        "local_witness_surface": local_summary,
        "target_lift": lift,
        "pair_overlap_support": pair_overlap_support(upstream, local),
        "interpretation": {
            "main_question": "For primes hit by the upstream filter table at the local scale, can the local layer attach verified Rosati witnesses?",
            "target_level_answer_field": "target_lift.upstream_filtered_and_locally_witnessed_count",
            "boundary": "Pair overlap remains supporting context because different filter moduli can certify the same prime.",
        },
        "forbidden_actions": forbidden_actions(),
    }


def report(results: dict[str, Any]) -> str:
    if results["verdict"] == "UPSTREAM_FILTER_INPUT_UNAVAILABLE":
        return "\n".join([
            "# Erdos #242 Upstream Target Witness Lift",
            "",
            f"Experiment: `{results['experiment_id']}`",
            "Status: `REVIEW_ONLY`",
            f"Verdict: `{results['verdict']}`",
            "",
            f"Error: `{results.get('error')}`",
            "",
        ])

    upstream = results["upstream_target_scan"]
    local = results["local_witness_surface"]
    lift = results["target_lift"]
    pair = results["pair_overlap_support"]
    return "\n".join([
        "# Erdos #242 Upstream Target Witness Lift",
        "",
        f"Experiment: `{results['experiment_id']}`",
        "Status: `REVIEW_ONLY`",
        f"Verdict: `{results['verdict']}`",
        "",
        "## Meaning",
        "",
        "This gate moves from residue-pair overlap to target-level meaning. It asks whether a hard-strip prime hit by Mihnea/Bogdan's pinned upstream filter table can be lifted locally to an explicit Rosati witness with verified denominators.",
        "",
        "## Upstream Target Scan",
        "",
        f"- Pinned commit: `{overlap_gate.PINNED_COMMIT}`",
        f"- Filter rows: `{upstream['row_count']}`",
        f"- Upstream moduli: `{upstream['filter_modulus_count']}`",
        f"- Upstream `(modulus,residue)` pairs: `{upstream['filter_pair_count']}`",
        f"- Filters SHA-256: `{upstream['sha256']}`",
        f"- SHA matches prior overlap gate: `{upstream['sha256_matches_prior_overlap']}`",
        f"- Upstream-filtered hard-strip targets: `{lift['upstream_filter_target_count']}`",
        f"- Total upstream target-hit records: `{upstream['total_target_hit_records']}`",
        f"- Upstream hit validation errors: `{upstream['hit_validation_error_count']}`",
        "",
        "## Local Witness Lift",
        "",
        f"- Max n: `{local['max_n']}`",
        f"- Constant bound: `{local['constant_bound']}`",
        f"- Hard-strip targets: `{local['target_prime_count']}`",
        f"- Locally witnessed targets: `{lift['locally_witnessed_target_count']}`",
        f"- Upstream-filtered and locally witnessed: `{lift['upstream_filtered_and_locally_witnessed_count']}`",
        f"- Upstream-filtered but local witness missing: `{lift['upstream_filtered_local_witness_missing_count']}`",
        f"- Local witness but no observed upstream hit: `{lift['local_witness_without_observed_upstream_hit_count']}`",
        f"- Invalid emitted witnesses: `{lift['emitted_witness_invalid_count']}`",
        f"- Local certificates: `{local['certificate_count']}`",
        f"- Local filter pairs: `{local['filter_pair_count']}`",
        "",
        "## Supporting Pair Context",
        "",
        f"- Prior pair intersection: `{pair.get('previous_pair_intersection_count')}`",
        f"- Prior local-pair overlap fraction: `{pair.get('previous_intersection_fraction_of_local')}`",
        f"- Prior upstream-pair overlap fraction: `{pair.get('previous_intersection_fraction_of_upstream')}`",
        "",
        "## Boundary",
        "",
        "This is not `10^18` replication. It is not exact pipeline equality. It is not a proof of Erdős-Straus. It is a local certificate-enrichment layer over upstream filter hits.",
        "",
    ])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=1_000_000)
    parser.add_argument("--constant-bound", type=int, default=21)
    parser.add_argument("--sample-size", type=int, default=40)
    parser.add_argument("--sample-per-target", type=int, default=3)
    args = parser.parse_args()

    results = run(args)
    result_path = RESULT_DIR / f"{EXP_ID}_RESULTS.json"
    report_path = RESULT_DIR / f"{EXP_ID}_REPORT.md"
    sha_path = RESULT_DIR / f"{EXP_ID}_RESULTS.sha256"
    write_json(result_path, results)
    write_text(report_path, report(results))
    write_text(sha_path, sha256_file(result_path) + "\n")
    lift = results.get("target_lift", {})
    upstream = results.get("upstream_target_scan", {})
    print(json.dumps({
        "experiment_id": EXP_ID,
        "results": str(result_path.relative_to(ROOT)),
        "report": str(report_path.relative_to(ROOT)),
        "sha256": str(sha_path.relative_to(ROOT)),
        "verdict": results["verdict"],
        "upstream_filtered_targets": lift.get("upstream_filter_target_count"),
        "lifted_targets": lift.get("upstream_filtered_and_locally_witnessed_count"),
        "missing_local_witnesses": lift.get("upstream_filtered_local_witness_missing_count"),
        "filters_sha256": upstream.get("sha256"),
        "filters_sha256_matches_prior_overlap": upstream.get("sha256_matches_prior_overlap"),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
