#!/usr/bin/env python3
"""Local n=15 slice pilot plus Tao bridge ledger for Erdős #114.

This run is intentionally narrow. It does not launch the full n=15 branch-and-
bound engine. It binds three selected n=14 seed slices to their n=15 analog
work requests, keeps Tao's large-degree constants as a ledger, and emits one
review-only artifact triplet.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-N15-LOCAL-SLICE-TAO-BRIDGE-20260508-01"
CLAIM_CEILING = (
    "We have a finite certified packet through degree 14 and an independent "
    "large-degree theorem by Tao, but the effective finite bridge is still open."
)

ROOT = Path(__file__).resolve().parents[1]
ERDOS114 = ROOT / "Erdos114"
PROOF_PATH = ERDOS114 / "proof_path"
VALIDATED = ERDOS114 / "validated_length"
RESULTS_114 = ROOT / "results" / "erdos-114"
OUTDIR = PROOF_PATH / EXPERIMENT_ID
DOWNLOADS_STORY_DIR = Path.home() / "Downloads" / "Eratosthenes_Ahmes_Stories"

SOURCES = {
    "finite_less15_packet": PROOF_PATH
    / "EXP-MATH-EHP114-FINITE-N-LESS-15-PROOF-PACKET-20260507-01"
    / "EXP-MATH-EHP114-FINITE-N-LESS-15-PROOF-PACKET-20260507-01_RESULTS.json",
    "n15_interval_contract": PROOF_PATH
    / "EXP-MATH-EHP114-N15-INTERVAL-PILOT-20260508-01"
    / "EXP-MATH-EHP114-N15-INTERVAL-PILOT-20260508-01_RESULTS.json",
    "n15_envelope": PROOF_PATH
    / "EXP-MATH-EHP114-N15-PILOT-ENVELOPE-20260508-01"
    / "EXP-MATH-EHP114-N15-PILOT-ENVELOPE-20260508-01_RESULTS.json",
    "n15_fourier_diagnostic": RESULTS_114
    / "EXP-MATH-EHP114-N15-FOURIER-HESSIAN-20260502-02_RESULTS.json",
    "n14_low_margin_local_cell": VALIDATED
    / "EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-07-03-20260506-01"
    / "EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-07-03-20260506-01_RESULTS.json",
    "n14_regular_residual_cell": VALIDATED
    / "EXP-MATH-EHP114-N14-VARIABLE-REGULAR-RESIDUAL-CLOSURE-CELL-06-03-20260506-01"
    / "EXP-MATH-EHP114-N14-VARIABLE-REGULAR-RESIDUAL-CLOSURE-CELL-06-03-20260506-01_RESULTS.json",
    "n14_hard_collar_cell": VALIDATED
    / "EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-02-03-20260506-01"
    / "EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-02-03-20260506-01_RESULTS.json",
    "tao_blocker_ledger": PROOF_PATH
    / "EXP-MATH-EHP114-TAO-FINAL-SECTION-BLOCKER-LEDGER-20260508-01"
    / "EXP-MATH-EHP114-TAO-FINAL-SECTION-BLOCKER-LEDGER-20260508-01_RESULTS.json",
    "tao_final_constant_extraction": PROOF_PATH
    / "EXP-MATH-EHP114-TAO-FINAL-SECTION-CONSTANT-EXTRACTION-20260507-01"
    / "EXP-MATH-EHP114-TAO-FINAL-SECTION-CONSTANT-EXTRACTION-20260507-01_RESULTS.json",
    "tao_parent_constant_map": PROOF_PATH
    / "EXP-MATH-EHP114-TAO-PARENT-CONSTANT-SOURCE-MAP-20260507-01"
    / "EXP-MATH-EHP114-TAO-PARENT-CONSTANT-SOURCE-MAP-20260507-01_RESULTS.json",
    "tao_effective_checker": PROOF_PATH
    / "EXP-MATH-EHP114-TAO-EFFECTIVE-THRESHOLD-CHECKER-20260507-01"
    / "EXP-MATH-EHP114-TAO-EFFECTIVE-THRESHOLD-CHECKER-20260507-01_RESULTS.json",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_sha(path: Path) -> str:
    digest = sha256_file(path)
    sha_path = path.with_name(path.name.replace("_RESULTS.json", "_RESULTS.sha256"))
    sha_path.write_text(f"{digest}  {path.name}\n", encoding="utf-8")
    return digest


def source_fingerprints() -> list[dict[str, Any]]:
    rows = []
    for role, path in sorted(SOURCES.items()):
        rows.append(
            {
                "role": role,
                "path": str(path),
                "exists": path.exists(),
                "byte_count": path.stat().st_size if path.exists() else None,
                "sha256": sha256_file(path) if path.exists() else None,
            }
        )
    return rows


def require_sources() -> None:
    missing = [str(path) for path in SOURCES.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("missing required source artifacts: " + ", ".join(missing))


def selected_slice(role: str, source_key: str, n15_request: str) -> dict[str, Any]:
    source = read_json(SOURCES[source_key])
    return {
        "role": role,
        "source_path": str(SOURCES[source_key]),
        "source_sha256": sha256_file(SOURCES[source_key]),
        "n14_seed_cell": source.get("cell_tag"),
        "n14_seed_status": source.get("status"),
        "n14_margin_to_cap": source.get("margin_to_cap"),
        "n14_first_failed_condition": source.get("first_failed_condition"),
        "n15_slice_status": "SELECTED_FOR_LOCAL_REDUCTION_NOT_EXECUTED",
        "n15_request": n15_request,
        "why_this_slice": (
            "This is one of the three smallest local probes that can test whether "
            "n14 cell/collar machinery transfers to degree 15 without launching "
            "the full search."
        ),
        "acceptance_gate_for_future_run": [
            "write a fresh n15 slice RESULTS/REPORT/SHA triplet",
            "record interval widths and margin-to-cap if a length bound is emitted",
            "if branch-and-bound is invoked, record nonzero evaluation count",
            "do not emit a global n15 certificate row from a slice result",
        ],
    }


def build_slice_section() -> dict[str, Any]:
    fourier = read_json(SOURCES["n15_fourier_diagnostic"])
    envelope = read_json(SOURCES["n15_envelope"])
    return {
        "section": "n15_local_slice_pilot",
        "status": "NEEDS_N15_SLICE_RUNNERS",
        "full_n15_branch_and_bound_launched": False,
        "source_envelope_decision": envelope.get("decision"),
        "source_envelope_stop_go": envelope.get("stop_go_rule_result"),
        "n15_existing_diagnostic": {
            "path": str(SOURCES["n15_fourier_diagnostic"]),
            "sha256": sha256_file(SOURCES["n15_fourier_diagnostic"]),
            "status": fourier.get("status"),
            "degree": fourier.get("degree"),
            "rigorous": bool(fourier.get("method", {}).get("rigorous"))
            if isinstance(fourier.get("method"), dict)
            else False,
        },
        "selected_slices": [
            selected_slice(
                "low_margin_local_cell",
                "n14_low_margin_local_cell",
                "instantiate n15 local-cell source generation and residual closure accounting from CELL-07-03",
            ),
            selected_slice(
                "regular_residual_analog",
                "n14_regular_residual_cell",
                "instantiate n15 regular residual closure from CELL-06-03 and demand a positive margin",
            ),
            selected_slice(
                "hard_collar_root_isolation_analog",
                "n14_hard_collar_cell",
                "instantiate n15 collar/root-isolation pilot from CELL-02-03 and classify the wall-separation failure",
            ),
        ],
        "local_run_decision": "NEEDS_REDUCTION",
        "reason": (
            "There is no n15-specific slice runner in the current local code. "
            "The general n15 engine would be a full branch-and-bound launch, so it was not run."
        ),
    }


def final_section_rows() -> list[dict[str, Any]]:
    extraction = read_json(SOURCES["tao_final_constant_extraction"])
    wanted = {"inside-2", "annulus-2", "outside-again"}
    rows = []
    for row in extraction.get("extraction_rows", []):
        if row.get("source_label") not in wanted:
            continue
        rows.append(
            {
                "source_label": row.get("source_label"),
                "dependency_id": row.get("dependency_id"),
                "priority": row.get("priority"),
                "extraction_status": row.get("extraction_status"),
                "extracted_N_i": row.get("extracted_N_i"),
                "named_blocker": row.get("named_blocker"),
                "named_constants": row.get("named_constants", []),
                "dependency_parents": row.get("dependency_parents", []),
                "source_span": row.get("source_span"),
                "next_required_work": row.get("inequality_needed"),
            }
        )
    rows.sort(key=lambda r: ["inside-2", "annulus-2", "outside-again"].index(r["source_label"]))
    return rows


def intermediate_rows_by_label() -> dict[str, dict[str, Any]]:
    extraction = read_json(SOURCES["tao_final_constant_extraction"])
    wanted = {"pots", "ets"}
    out = {}
    for row in extraction.get("extraction_rows", []):
        label = row.get("classified_label")
        if label not in wanted:
            continue
        out.setdefault(
            label,
            {
                "label": label,
                "dependency_ids": [],
                "source_spans": [],
                "extraction_statuses": set(),
                "named_blockers": set(),
                "named_constants": set(),
                "next_required_work": [],
            },
        )
        target = out[label]
        target["dependency_ids"].append(row.get("dependency_id"))
        target["source_spans"].append(row.get("source_span"))
        target["extraction_statuses"].add(row.get("extraction_status"))
        target["named_blockers"].add(row.get("named_blocker"))
        for constant in row.get("named_constants", []):
            target["named_constants"].add(constant)
        if row.get("inequality_needed"):
            target["next_required_work"].append(row.get("inequality_needed"))
    cleaned = {}
    for label, row in out.items():
        cleaned[label] = {
            "label": row["label"],
            "dependency_ids": row["dependency_ids"],
            "source_spans": row["source_spans"],
            "extraction_statuses": sorted(x for x in row["extraction_statuses"] if x),
            "named_blockers": sorted(x for x in row["named_blockers"] if x),
            "named_constants": sorted(x for x in row["named_constants"] if x),
            "next_required_work": sorted(set(row["next_required_work"])),
        }
    return cleaned


def parent_rows_by_label() -> dict[str, dict[str, Any]]:
    parent = read_json(SOURCES["tao_parent_constant_map"])
    out = {}
    for row in parent.get("parent_source_rows", []):
        out[row.get("label")] = {
            "label": row.get("label"),
            "title": row.get("title"),
            "environment": row.get("environment"),
            "source_span": row.get("source_span"),
            "quantification_status": row.get("quantification_status"),
        }
    return out


def build_tao_section() -> dict[str, Any]:
    checker = read_json(SOURCES["tao_effective_checker"])
    parent = read_json(SOURCES["tao_parent_constant_map"])
    parents = parent_rows_by_label()
    rows = final_section_rows()
    parent_labels = sorted({p for row in rows for p in row.get("dependency_parents", [])})
    intermediate = intermediate_rows_by_label()
    parent_map = [parents[label] for label in parent_labels if label in parents]
    parent_map.extend(intermediate[label] for label in parent_labels if label in intermediate)
    return {
        "section": "tao_bridge_ledger",
        "status": "NO_FINITE_BRIDGE_YET",
        "candidate_N0": checker.get("candidate_N0"),
        "dependency_count": checker.get("dependency_count"),
        "explicit_dependency_count": checker.get("explicit_dependency_count"),
        "opaque_dependency_count": checker.get("opaque_dependency_count"),
        "first_failed_condition": checker.get("first_failed_condition"),
        "target_blockers": rows,
        "immediate_parent_rows_to_quantify": parent_map,
        "constant_rows_still_opaque": [
            {
                "constant_symbol": row.get("constant_symbol"),
                "needed_for": row.get("needed_for"),
                "role": row.get("role"),
                "source_parents": row.get("source_parents"),
            }
            for row in parent.get("constant_rows", [])
            if row.get("quantification_status") == "PARENT_CONSTANTS_OPAQUE"
        ],
        "ledger_decision": "OPAQUE_BLOCKERS_REMAIN",
        "reason": (
            "The final labels inside-2, annulus-2, and outside-again are mapped, "
            "but their absolute constants and large-n gates remain nonnumeric."
        ),
    }


def build_result() -> dict[str, Any]:
    require_sources()
    slice_section = build_slice_section()
    tao_section = build_tao_section()
    return {
        "experiment_id": EXPERIMENT_ID,
        "created_unix": int(time.time()),
        "generated_by": "erdos_atlas_autoresearch_librarian",
        "origin": "auto-research",
        "persona": "Eratosthenes of Cyrene",
        "short_name": "Eratosthenes",
        "scribe": "Ahmes",
        "story_writer": "Ahmes",
        "problem_id": 114,
        "promotion_state": "review_only",
        "claim_ceiling": CLAIM_CEILING,
        "final_decision": "NEEDS_REDUCTION_AND_TAO_CONSTANTS",
        "full_n15_authorized": False,
        "tao_finite_bridge_authorized": False,
        "n15_slice_section": slice_section,
        "tao_bridge_section": tao_section,
        "source_fingerprints": source_fingerprints(),
        "forbidden_writes_observed": False,
        "forbidden_claims": [
            "global n15 certificate",
            "finite/high-degree bridge closure",
            "all-degree conclusion",
            "formal proof-status promotion",
        ],
        "next_executable_experiments": [
            "EXP-MATH-EHP114-N15-SLICE-CELL-07-03-LOCAL-<new-id>",
            "EXP-MATH-EHP114-N15-REGULAR-RESIDUAL-CELL-06-03-LOCAL-<new-id>",
            "EXP-MATH-EHP114-N15-COLLAR-ROOT-ISOLATION-CELL-02-03-LOCAL-<new-id>",
            "EXP-MATH-EHP114-TAO-PARENT-CONSTANT-NUMERICIZATION-<new-id>",
        ],
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    slice_rows = result["n15_slice_section"]["selected_slices"]
    tao_rows = result["tao_bridge_section"]["target_blockers"]
    body = [
        f"# {EXPERIMENT_ID} Report",
        "",
        "## Claim Ceiling",
        "",
        CLAIM_CEILING,
        "",
        "## Decision",
        "",
        f"- Final decision: `{result['final_decision']}`",
        f"- Full n=15 authorized: `{str(result['full_n15_authorized']).lower()}`",
        f"- Tao finite bridge authorized: `{str(result['tao_finite_bridge_authorized']).lower()}`",
        "",
        "## n=15 Local Slice Pilot",
        "",
        f"- Status: `{result['n15_slice_section']['status']}`",
        f"- Local run decision: `{result['n15_slice_section']['local_run_decision']}`",
        "- Full branch-and-bound launched: `false`",
        "",
    ]
    for row in slice_rows:
        body.append(
            f"- `{row['role']}` from `{row['n14_seed_cell']}`: "
            f"`{row['n15_slice_status']}`; next: {row['n15_request']}"
        )
    body.extend(
        [
            "",
            "## Tao Bridge Ledger",
            "",
            f"- Status: `{result['tao_bridge_section']['status']}`",
            f"- Candidate N0: `{result['tao_bridge_section']['candidate_N0']}`",
            f"- Opaque dependencies: `{result['tao_bridge_section']['opaque_dependency_count']}`",
            "",
        ]
    )
    for row in tao_rows:
        body.append(
            f"- `{row['source_label']}`: `{row['extraction_status']}`; "
            f"needed: {row['next_required_work']}"
        )
    body.extend(
        [
            "",
            "## Meaning",
            "",
            "The useful next work is no longer broad exploration. It is to write "
            "three n=15 slice runners and numericize the Tao parent constants. "
            "Until one of those artifacts changes, the finite bridge remains open.",
            "",
        ]
    )
    path.write_text("\n".join(body), encoding="utf-8")


def write_story(result: dict[str, Any]) -> Path:
    DOWNLOADS_STORY_DIR.mkdir(parents=True, exist_ok=True)
    path = DOWNLOADS_STORY_DIR / f"AHMES_STORY_{EXPERIMENT_ID}.md"
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing Ahmes story: {path}")
    story = f"""# Ahmes Story: {EXPERIMENT_ID}

Eratosthenes did not try to force degree 15 through the full machine. He took
the three smallest doors from degree 14 and marked what each would have to prove
at degree 15: the low-margin local cell, the regular residual cell, and the
hard collar/root-isolation cell.

He also kept Tao's large-degree theorem in its proper place. The labels
`inside-2`, `annulus-2`, and `outside-again` are now the bridge ledger, but no
integer threshold has emerged from them.

## Receipt

- Run ID: `{EXPERIMENT_ID}`
- Decision: `{result["final_decision"]}`
- Full n=15 authorized: `{str(result["full_n15_authorized"]).lower()}`
- Tao finite bridge authorized: `{str(result["tao_finite_bridge_authorized"]).lower()}`
- Claim ceiling: {CLAIM_CEILING}

## Next Doors

- n15 local-cell runner for `CELL-07-03`
- n15 regular-residual runner for `CELL-06-03`
- n15 collar/root-isolation runner for `CELL-02-03`
- Tao parent-constant numericization for `inside-2`, `annulus-2`, `outside-again`

This is a reading companion. The JSON/report/hash triplet is the source of truth.
"""
    path.write_text(story, encoding="utf-8")
    return path


def materialize(result: dict[str, Any], *, write_downloads_story: bool) -> dict[str, str]:
    if OUTDIR.exists():
        raise FileExistsError(f"refusing to overwrite existing run folder: {OUTDIR}")
    OUTDIR.mkdir(parents=True)
    result_path = OUTDIR / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = OUTDIR / f"{EXPERIMENT_ID}_REPORT.md"
    write_json(result_path, result)
    write_report(result, report_path)
    digest = write_sha(result_path)
    story_path = write_story(result) if write_downloads_story else None
    return {
        "results": str(result_path),
        "report": str(report_path),
        "sha256": str(result_path.with_name(f"{EXPERIMENT_ID}_RESULTS.sha256")),
        "digest": digest,
        "story": str(story_path) if story_path else "",
    }


def main() -> int:
    global EXPERIMENT_ID, OUTDIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-id", default=EXPERIMENT_ID)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--no-downloads-story", action="store_true")
    args = parser.parse_args()
    EXPERIMENT_ID = args.experiment_id
    OUTDIR = PROOF_PATH / EXPERIMENT_ID
    result = build_result()
    if args.dry_run:
        print(json.dumps({"experiment_id": EXPERIMENT_ID, "would_write": str(OUTDIR), "result": result}, indent=2, sort_keys=True))
        return 0
    result["materialized_paths"] = materialize(
        result,
        write_downloads_story=not args.no_downloads_story,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
