#!/usr/bin/env python3
"""Review-only theorem packet for the Erdos #242 Salez filter lemmas.

This runner connects the local certificate generator to the Lean theorem names
in ``Erdos242SalezFilters.lean``. It reports whether every locally emitted
certificate is covered by one of the seven formal Salez-style sufficiency
lemmas, and whether the current 1e6 hard-strip certificate packet remains
complete.

It is not a computational-frontier claim, not 10^18 replication, and not a
proof of Erdos-Straus.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import salez_filter_certificates as certs


EXP_ID = "EXP-MATH-ERDOS242-SALEZ-FILTER-THEOREM-PACKET-20260508-01"
ROOT = Path(__file__).resolve().parents[3]
RESULT_DIR = ROOT / "erdos-experiments" / "results" / "erdos-242"
LEAN_PATH = ROOT / "erdos-experiments" / "Erdos242" / "lean" / "Erdos242SalezFilters.lean"

THEOREM_MANIFEST: dict[str, dict[str, str]] = {
    "eqmod1a": {
        "lean_theorem": "salez_eqmod1a_clearedWitness",
        "rosati_case": "Rosati_1",
        "condition_shape": "B + p*C + A = 4*A*B*C*D",
    },
    "eqmod1b": {
        "lean_theorem": "salez_eqmod1b_clearedWitness",
        "rosati_case": "Rosati_1",
        "condition_shape": "A+B = C*E and p+E = 4*A*B*D",
    },
    "eqmod1c": {
        "lean_theorem": "salez_eqmod1c_clearedWitness",
        "rosati_case": "Rosati_1",
        "condition_shape": "p+E = 4*A*B*D and 4*B*D*E*C = p+E+4*B*B*D",
    },
    "eqmod2a": {
        "lean_theorem": "salez_eqmod2a_clearedWitness",
        "rosati_case": "Rosati_2",
        "condition_shape": "A+B = C*E and p*E+1 = 4*A*B*D",
    },
    "eqmod2b": {
        "lean_theorem": "salez_eqmod2b_clearedWitness",
        "rosati_case": "Rosati_2",
        "condition_shape": "p+F = 4*B*C*D and p*B+C = A*F",
    },
    "eqmod2c": {
        "lean_theorem": "salez_eqmod2c_clearedWitness",
        "rosati_case": "Rosati_2",
        "condition_shape": "A+B = C*E, p+F = 4*B*D*C, and 4*B*B*D+1 = E*F",
    },
    "eqmod2d": {
        "lean_theorem": "salez_eqmod2d_clearedWitness",
        "rosati_case": "Rosati_2",
        "condition_shape": "p+F = 4*B*C*D and p*B+C = A*F",
    },
}


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


def lean_manifest_status() -> dict[str, Any]:
    text = LEAN_PATH.read_text(encoding="utf-8")
    expected_theorems = [entry["lean_theorem"] for entry in THEOREM_MANIFEST.values()]
    present = [name for name in expected_theorems if f"theorem {name}" in text]
    missing = [name for name in expected_theorems if name not in present]
    forbidden_hits = []
    for token in ("sorry", "admit", "axiom"):
        for line_number, line in enumerate(text.splitlines(), start=1):
            if token in line:
                forbidden_hits.append({"token": token, "line": line_number, "text": line.strip()})
    return {
        "lean_file": str(LEAN_PATH.relative_to(ROOT)),
        "lean_file_sha256": sha256_file(LEAN_PATH),
        "expected_theorem_count": len(expected_theorems),
        "present_theorem_count": len(present),
        "present_theorems": present,
        "missing_theorems": missing,
        "forbidden_token_hits": forbidden_hits,
    }


def certificate_theorem_coverage(max_n: int, constant_bound: int) -> dict[str, Any]:
    started = time.perf_counter()
    certificates = certs.generate_filter_certificates(
        max_n=max_n,
        constant_bound=constant_bound,
        hard_strip_only=True,
    )
    targets = certs.target_primes(max_n, hard_strip_only=True)
    chosen = certs.choose_one_certificate_per_prime(certificates)
    invalid = [c for c in certificates if not c.verified]

    equation_counts = Counter(c.equation_id for c in certificates)
    rosati_counts = Counter(c.rosati_case for c in certificates)
    known_labels = set(THEOREM_MANIFEST)
    unknown_label_counts = {
        label: count
        for label, count in sorted(equation_counts.items())
        if label not in known_labels
    }

    first_certificate_by_equation: dict[str, Any] = {}
    for cert in certificates:
        first_certificate_by_equation.setdefault(cert.equation_id, cert.to_json())

    return {
        "max_n": max_n,
        "constant_bound": constant_bound,
        "target_strip": "prime p = 1 mod 24",
        "target_prime_count": len(targets),
        "certified_target_count": len([p for p in targets if p in chosen]),
        "missing_target_count": len([p for p in targets if p not in chosen]),
        "certificate_count": len(certificates),
        "invalid_certificate_count": len(invalid),
        "equation_counts": dict(sorted(equation_counts.items())),
        "rosati_case_counts": dict(sorted(rosati_counts.items())),
        "unknown_label_counts": unknown_label_counts,
        "all_certificates_have_formal_label": not unknown_label_counts,
        "all_targets_have_witness": all(p in chosen for p in targets),
        "first_invalid_certificate_sample": [c.to_json() for c in invalid[:20]],
        "first_certificate_by_equation": first_certificate_by_equation,
        "runtime_seconds": round(time.perf_counter() - started, 6),
    }


def verdict_for(lean_status: dict[str, Any], coverage: dict[str, Any]) -> str:
    if (
        lean_status["missing_theorems"]
        or lean_status["forbidden_token_hits"]
        or coverage["invalid_certificate_count"]
        or not coverage["all_certificates_have_formal_label"]
        or not coverage["all_targets_have_witness"]
    ):
        return "FORMAL_PACKET_NEEDS_REVIEW"
    return "SALEZ_FILTER_LEMMAS_FORMALIZED_LOCAL"


def source_manifest() -> dict[str, Any]:
    return {
        "status": "REVIEW_ONLY_LOCAL_THEOREM_PACKET",
        "third_party_code_committed": False,
        "salez_arxiv": "https://arxiv.org/abs/1406.6307",
        "mihnea_bogdan_arxiv": "https://arxiv.org/abs/2509.00128",
        "mihnea_bogdan_github": "https://github.com/esc-paper/erdos-straus",
        "local_implementation_status": "FORMAL_CERTIFICATE_LAYER_RECONSTRUCTED_FROM_SOURCES",
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    lean_status = lean_manifest_status()
    coverage = certificate_theorem_coverage(args.max_n, args.constant_bound)
    verdict = verdict_for(lean_status, coverage)
    return {
        "experiment_id": EXP_ID,
        "generated_at": now_iso(),
        "status": "REVIEW_ONLY",
        "verdict": verdict,
        "claim_ceiling": "formal certificate infrastructure for seven Salez-style filter lemmas; not a new computational bound, not 10^18 replication, and not a proof of Erdos-Straus",
        "source_manifest": source_manifest(),
        "theorem_manifest": THEOREM_MANIFEST,
        "lean_manifest_status": lean_status,
        "certificate_theorem_coverage": coverage,
        "interpretation": {
            "what_this_adds": "The local certificate generator now has a Lean-checked algebraic sufficiency layer for every emitted Salez-style equation label.",
            "what_this_does_not_add": "No new global verification bound and no exact formalization of the upstream Mihnea/Bogdan sieve table.",
            "public_worthy_axis": "useful verification infrastructure and a clean formal theorem subset",
        },
        "forbidden_actions": [
            "do not claim proof of Erdos-Straus",
            "do not claim beyond-10^18 verification",
            "do not claim exact equality to Mihnea/Bogdan's optimized pipeline",
            "do not update D1/proof registry from this review-only packet",
        ],
    }


def report(results: dict[str, Any]) -> str:
    lean = results["lean_manifest_status"]
    coverage = results["certificate_theorem_coverage"]
    lines = [
        "# Erdos #242 Salez Filter Theorem Packet",
        "",
        f"Experiment: `{results['experiment_id']}`",
        "Status: `REVIEW_ONLY`",
        f"Verdict: `{results['verdict']}`",
        "",
        "## Meaning",
        "",
        "This packet turns the local #242 certificate layer into formal certificate infrastructure. The Lean file proves seven Salez-style sufficient-condition lemmas that produce cleared Erdős-Straus witnesses.",
        "",
        "## Lean Manifest",
        "",
        f"- Lean file: `{lean['lean_file']}`",
        f"- Lean SHA-256: `{lean['lean_file_sha256']}`",
        f"- Expected theorem labels present: `{lean['present_theorem_count']} / {lean['expected_theorem_count']}`",
        f"- Forbidden Lean tokens: `{len(lean['forbidden_token_hits'])}`",
        "",
        "## Certificate Coverage",
        "",
        f"- Target: `{coverage['target_strip']}`, `p <= {coverage['max_n']}`",
        f"- Constant bound: `{coverage['constant_bound']}`",
        f"- Certified targets: `{coverage['certified_target_count']} / {coverage['target_prime_count']}`",
        f"- Missing targets: `{coverage['missing_target_count']}`",
        f"- Certificates: `{coverage['certificate_count']}`",
        f"- Invalid certificates: `{coverage['invalid_certificate_count']}`",
        f"- Unknown equation labels: `{coverage['unknown_label_counts']}`",
        f"- Equation counts: `{coverage['equation_counts']}`",
        "",
        "## Boundary",
        "",
        "This is formal certificate infrastructure, not a new computational bound. It is not `10^18` replication, not exact upstream pipeline equality, and not a proof of Erdős-Straus.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=1_000_000)
    parser.add_argument("--constant-bound", type=int, default=21)
    args = parser.parse_args()

    results = run(args)
    result_path = RESULT_DIR / f"{EXP_ID}_RESULTS.json"
    report_path = RESULT_DIR / f"{EXP_ID}_REPORT.md"
    sha_path = RESULT_DIR / f"{EXP_ID}_RESULTS.sha256"
    write_json(result_path, results)
    write_text(report_path, report(results))
    write_text(sha_path, sha256_file(result_path) + "\n")
    coverage = results["certificate_theorem_coverage"]
    print(json.dumps({
        "experiment_id": EXP_ID,
        "results": str(result_path.relative_to(ROOT)),
        "report": str(report_path.relative_to(ROOT)),
        "sha256": str(sha_path.relative_to(ROOT)),
        "verdict": results["verdict"],
        "certificates": coverage["certificate_count"],
        "certified_targets": coverage["certified_target_count"],
        "target_primes": coverage["target_prime_count"],
        "invalid_certificates": coverage["invalid_certificate_count"],
        "missing_targets": coverage["missing_target_count"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
