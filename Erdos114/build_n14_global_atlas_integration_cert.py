#!/usr/bin/env python3
"""EHP114 n=14 global atlas integration certificate.

Composes 64 already-emitted local-cell certificate packets (8x8 root-affine
parameter grid, with hard-cell (6,4) using the LOCAL-HARD-CELL-CERTIFICATE
naming convention) into one global integration artifact.

This is a bookkeeping orchestration layer only. It does not recompute interval
geometry. Per the IN_SEARCH gate #2, the integration must prove:

  - 64 cells covered: every (sub_i, sub_j) in 0..7 x 0..7 is present
  - per-cell pass status: each local-cell packet has status
    LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF and pass=true
  - per-cell SHA integrity: each cell's own RESULTS.json SHA-256 matches
    the sidecar .sha256 file
  - per-cell child-source SHA integrity: each cell's claimed source SHAs
    are PASS (already verified inside the cell, but re-asserted here)
  - cross-cell source-id uniqueness: no two cells claim the same L24/L27/
    L28/L29/L31 source experiment ID
  - cap consistency: every cell uses the same exact_length_cap
  - global length bound: max over cells of total_validated_length_upper
    is below the exact_length_cap
  - bookkeeping zeros: every cell has zero ownership_duplicate, zero
    source_filter_mismatch, zero candidate_count_drift, zero source_sha_fail

The output is a single integration certificate at
EXP-MATH-EHP114-N14-GLOBAL-ATLAS-INTEGRATION-CERT-20260508-01.

Claim ceiling: this is a global atlas integration certificate for n=14 only.
It is not a proof of Erdos #114, not a result for n != 14, and not a Tao
bridge. The all-degree conjecture remains open.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
VALIDATED = SCRIPT_DIR / "validated_length"

EXACT_LENGTH_CAP = 20.672796062619668

OUT_ID = "EXP-MATH-EHP114-N14-GLOBAL-ATLAS-INTEGRATION-CERT-20260508-01"

PASS_LOCAL_STATUS = "LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF"
PASS_GLOBAL_STATUS = "GLOBAL_ATLAS_INTEGRATION_PASS_NOT_FULL_PROOF"
FAIL_GLOBAL_STATUS = "GLOBAL_ATLAS_INTEGRATION_BLOCKED"

CLAIM_CEILING = (
    "Global n=14 atlas integration certificate only. Not a proof of "
    "Erdos #114, not a result for n != 14, and not a Tao-bridge or "
    "all-degree statement. The all-degree EHP conjecture remains open."
)

EXPECTED_SOURCE_ROLES = (
    "l24_regular_monotone",
    "l27_pprime_p8",
    "l28_pprime_root_location",
    "l29_sharp_pprime_root_location",
    "l31_residual_chain_integration",
)


def cell_tag(sub_i: int, sub_j: int) -> str:
    return f"CELL-{sub_i:02d}-{sub_j:02d}"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def find_canonical_packet(sub_i: int, sub_j: int) -> Path:
    """Locate the canonical local-cell certificate packet directory for (i,j).

    Hard cell (6,4) uses LOCAL-HARD-CELL-CERTIFICATE-PACKET-* (no CELL suffix).
    Other cells use LOCAL-CELL-CERTIFICATE-PACKET-CELL-ii-jj-*; if multiple
    versions exist (-01, -02, ...), we pick the highest-numbered version on the
    most recent date.
    """
    if (sub_i, sub_j) == (6, 4):
        candidates = sorted(VALIDATED.glob("EXP-MATH-EHP114-N14-LOCAL-HARD-CELL-CERTIFICATE-PACKET-*"))
        if not candidates:
            raise FileNotFoundError("hard-cell packet missing for (6,4)")
        return candidates[-1]

    tag = cell_tag(sub_i, sub_j)
    pattern = f"EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-{tag}-*"
    candidates = sorted(VALIDATED.glob(pattern))
    if not candidates:
        raise FileNotFoundError(f"no local-cell packet found for {tag}")
    return candidates[-1]


def load_packet(pkt_dir: Path) -> dict[str, Any]:
    rj = pkt_dir / f"{pkt_dir.name}_RESULTS.json"
    if not rj.exists():
        raise FileNotFoundError(f"RESULTS.json missing in {pkt_dir.name}")
    with rj.open("r", encoding="utf-8") as f:
        return json.load(f)


def verify_own_sha(pkt_dir: Path) -> tuple[bool, str, str]:
    """Recompute the cell's own RESULTS.json SHA-256, compare to .sha256 sidecar.

    Returns (pass, actual, expected_or_message).
    """
    rj = pkt_dir / f"{pkt_dir.name}_RESULTS.json"
    sha = pkt_dir / f"{pkt_dir.name}_RESULTS.sha256"
    if not sha.exists():
        return (False, "", f"sidecar missing: {sha.name}")
    actual = sha256_file(rj)
    expected_line = sha.read_text(encoding="utf-8").strip()
    expected = expected_line.split()[0] if expected_line else ""
    return (actual == expected, actual, expected)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Emit n=14 global atlas integration certificate."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run all checks and print the would-be report; do not write artifacts.",
    )
    args = parser.parse_args(argv)

    out_dir = VALIDATED / OUT_ID
    if out_dir.exists() and not args.dry_run:
        # Strict: never overwrite an existing experiment ID.
        print(
            f"REFUSE: output directory already exists: {out_dir}\n"
            "Per H2 platform invariants Rule 3, results are immutable. "
            "Bump the experiment ID suffix to retry.",
            file=sys.stderr,
        )
        return 2

    cells: list[dict[str, Any]] = []
    blockers: list[dict[str, Any]] = []
    seen_source_ids: dict[str, str] = {}
    cap_values_seen: set[float] = set()

    for i in range(8):
        for j in range(8):
            tag = cell_tag(i, j)
            entry: dict[str, Any] = {
                "cell_tag": tag,
                "sub_i": i,
                "sub_j": j,
                "checks": {},
            }
            try:
                pkt_dir = find_canonical_packet(i, j)
            except FileNotFoundError as exc:
                entry["packet_dir"] = None
                entry["packet_status"] = "MISSING"
                entry["error"] = str(exc)
                blockers.append({"cell_tag": tag, "kind": "missing_packet", "detail": str(exc)})
                cells.append(entry)
                continue

            entry["packet_dir"] = pkt_dir.name
            try:
                pkt = load_packet(pkt_dir)
            except (FileNotFoundError, json.JSONDecodeError) as exc:
                entry["packet_status"] = "UNREADABLE"
                entry["error"] = str(exc)
                blockers.append({"cell_tag": tag, "kind": "unreadable_packet", "detail": str(exc)})
                cells.append(entry)
                continue

            entry["packet_id"] = pkt.get("experiment_id")
            entry["packet_status"] = pkt.get("status")
            entry["total_validated_length_upper"] = pkt.get("total_validated_length_upper")
            entry["margin_to_cap"] = pkt.get("margin_to_cap")
            entry["exact_length_cap"] = pkt.get("exact_length_cap")
            entry["pass"] = pkt.get("local_hard_cell_certificate_pass", False)

            cap = pkt.get("exact_length_cap")
            if cap is not None:
                cap_values_seen.add(float(cap))

            # Check 1: status
            ok_status = pkt.get("status") == PASS_LOCAL_STATUS and pkt.get(
                "local_hard_cell_certificate_pass"
            ) is True
            entry["checks"]["status_pass"] = ok_status
            if not ok_status:
                blockers.append(
                    {
                        "cell_tag": tag,
                        "kind": "cell_status_not_pass",
                        "detail": f"status={pkt.get('status')!r}, pass={pkt.get('local_hard_cell_certificate_pass')!r}",
                    }
                )

            # Check 2: bookkeeping zeros
            zeros_ok = (
                pkt.get("source_sha_fail_count", 1) == 0
                and pkt.get("ownership_duplicate_count", 1) == 0
                and pkt.get("source_filter_mismatch_count", 1) == 0
                and pkt.get("candidate_count_drift_count", 1) == 0
                and pkt.get("l31_failure_count", 1) == 0
            )
            entry["checks"]["bookkeeping_zeros"] = zeros_ok
            if not zeros_ok:
                blockers.append(
                    {
                        "cell_tag": tag,
                        "kind": "cell_bookkeeping_nonzero",
                        "detail": (
                            f"sha_fail={pkt.get('source_sha_fail_count')} "
                            f"own_dup={pkt.get('ownership_duplicate_count')} "
                            f"src_mis={pkt.get('source_filter_mismatch_count')} "
                            f"cand_drift={pkt.get('candidate_count_drift_count')} "
                            f"l31_fail={pkt.get('l31_failure_count')}"
                        ),
                    }
                )

            # Check 3: length below cap
            length_upper = pkt.get("total_validated_length_upper")
            length_ok = (
                isinstance(length_upper, (int, float))
                and isinstance(cap, (int, float))
                and length_upper < cap
            )
            entry["checks"]["length_below_cap"] = bool(length_ok)
            if not length_ok:
                blockers.append(
                    {
                        "cell_tag": tag,
                        "kind": "cell_length_at_or_above_cap",
                        "detail": f"length_upper={length_upper}, cap={cap}",
                    }
                )

            # Check 4: own RESULTS.json SHA matches its sidecar
            own_sha_pass, own_actual, own_expected = verify_own_sha(pkt_dir)
            entry["checks"]["own_results_sha_pass"] = own_sha_pass
            entry["own_results_sha_actual"] = own_actual
            entry["own_results_sha_expected"] = own_expected
            if not own_sha_pass:
                blockers.append(
                    {
                        "cell_tag": tag,
                        "kind": "cell_own_sha_mismatch",
                        "detail": f"actual={own_actual} expected={own_expected}",
                    }
                )

            # Check 5: per-cell child-source SHA statuses all PASS
            sha_rows = pkt.get("source_sha_status", []) or []
            child_pass = all(row.get("sha_status") == "PASS" for row in sha_rows) and bool(sha_rows)
            entry["checks"]["child_source_sha_all_pass"] = child_pass
            entry["child_source_sha_row_count"] = len(sha_rows)
            if not child_pass:
                blockers.append(
                    {
                        "cell_tag": tag,
                        "kind": "cell_child_source_sha_not_all_pass",
                        "detail": f"row_count={len(sha_rows)}",
                    }
                )

            # Check 6: collect source experiment IDs for cross-cell uniqueness.
            # Skip None / empty values: some cells declare optional slots
            # (e.g. residual_universe_accounting) with a null placeholder when
            # the variable-universe accounting is not invoked. Those are not
            # ownership claims and must not register as collisions.
            source_ids = pkt.get("source_experiment_ids", {}) or {}
            entry["source_experiment_ids"] = source_ids
            for role, sid in source_ids.items():
                if not sid:
                    continue
                if sid in seen_source_ids:
                    blockers.append(
                        {
                            "cell_tag": tag,
                            "kind": "cross_cell_source_id_collision",
                            "detail": (
                                f"source id {sid!r} (role={role}) already claimed by "
                                f"{seen_source_ids[sid]}"
                            ),
                        }
                    )
                else:
                    seen_source_ids[sid] = f"{tag}:{role}"

            # Check 7: declared roles cover the expected L24/L27/L28/L29/L31 set
            missing_roles = [r for r in EXPECTED_SOURCE_ROLES if r not in source_ids]
            entry["checks"]["all_expected_roles_present"] = not missing_roles
            if missing_roles:
                blockers.append(
                    {
                        "cell_tag": tag,
                        "kind": "cell_missing_source_roles",
                        "detail": f"missing={missing_roles}",
                    }
                )

            cells.append(entry)

    # Aggregate
    processed = sum(1 for c in cells if c.get("packet_status") == PASS_LOCAL_STATUS)
    pass_count = sum(
        1 for c in cells if all(c.get("checks", {}).get(k) for k in c.get("checks", {}))
    )
    lengths = [
        c["total_validated_length_upper"]
        for c in cells
        if isinstance(c.get("total_validated_length_upper"), (int, float))
    ]
    margins = [
        c["margin_to_cap"]
        for c in cells
        if isinstance(c.get("margin_to_cap"), (int, float))
    ]

    coverage_complete = len(cells) == 64 and all(c.get("packet_dir") for c in cells)
    cap_consistent = len(cap_values_seen) == 1 and EXACT_LENGTH_CAP in cap_values_seen
    if not cap_consistent:
        blockers.append(
            {
                "cell_tag": None,
                "kind": "cap_inconsistent_across_cells",
                "detail": f"cap_values_seen={sorted(cap_values_seen)}",
            }
        )

    global_max_length = max(lengths) if lengths else None
    global_min_margin = min(margins) if margins else None
    global_pass = (
        coverage_complete
        and cap_consistent
        and pass_count == 64
        and not blockers
        and isinstance(global_max_length, (int, float))
        and global_max_length < EXACT_LENGTH_CAP
    )

    status = PASS_GLOBAL_STATUS if global_pass else FAIL_GLOBAL_STATUS
    first_failed_condition: str
    if global_pass:
        first_failed_condition = "none"
    elif not coverage_complete:
        first_failed_condition = "atlas_coverage_incomplete"
    elif not cap_consistent:
        first_failed_condition = "exact_length_cap_inconsistent"
    elif blockers:
        first_failed_condition = blockers[0]["kind"]
    else:
        first_failed_condition = "unknown"

    results: dict[str, Any] = {
        "experiment_id": OUT_ID,
        "status": status,
        "claim_ceiling": CLAIM_CEILING,
        "atlas_grid_dim_i": 8,
        "atlas_grid_dim_j": 8,
        "expected_cell_count": 64,
        "discovered_cell_count": sum(1 for c in cells if c.get("packet_dir")),
        "passing_cell_count": pass_count,
        "blocker_count": len(blockers),
        "first_failed_condition": first_failed_condition,
        "exact_length_cap": EXACT_LENGTH_CAP,
        "cap_consistent_across_cells": cap_consistent,
        "cap_values_seen": sorted(cap_values_seen),
        "global_max_length_upper": global_max_length,
        "global_min_margin_to_cap": global_min_margin,
        "global_length_below_cap": (
            isinstance(global_max_length, (int, float))
            and global_max_length < EXACT_LENGTH_CAP
        ),
        "tightest_cells": sorted(
            [
                {
                    "cell_tag": c["cell_tag"],
                    "total_validated_length_upper": c.get("total_validated_length_upper"),
                    "margin_to_cap": c.get("margin_to_cap"),
                }
                for c in cells
                if isinstance(c.get("margin_to_cap"), (int, float))
            ],
            key=lambda x: x["margin_to_cap"],
        )[:5],
        "loosest_cells": sorted(
            [
                {
                    "cell_tag": c["cell_tag"],
                    "total_validated_length_upper": c.get("total_validated_length_upper"),
                    "margin_to_cap": c.get("margin_to_cap"),
                }
                for c in cells
                if isinstance(c.get("margin_to_cap"), (int, float))
            ],
            key=lambda x: -x["margin_to_cap"],
        )[:5],
        "cross_cell_unique_source_ids": len(seen_source_ids),
        "cross_cell_collision_count": sum(
            1 for b in blockers if b["kind"] == "cross_cell_source_id_collision"
        ),
        "blockers": blockers,
        "cells": cells,
        "proof_obligations": [
            "atlas covers 8x8 root-affine grid (64 cells; cell (6,4) uses hard-cell packet)",
            "every cell packet status equals LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF",
            "every cell packet bookkeeping counters are zero",
            "every cell total_validated_length_upper < exact_length_cap",
            "every cell RESULTS.json SHA-256 matches its sidecar",
            "every cell child-source SHA status is PASS",
            "no source experiment ID is claimed by two different cells",
            "exact_length_cap is identical across all 64 cells",
        ],
        "next_dependency": (
            "n=15 boundary slice atlas (Track B in plan yes-wondrous-blum.md), "
            "Tao threshold effectivization (Track A), and finite-degree certificate "
            "index for n=3..14. This integration cert is necessary but not "
            "sufficient for an all-degree EHP proof."
        ),
        "timestamp_unix": str(int(time.time())),
        "schema_version": "1.0",
    }

    if args.dry_run:
        summary = {
            k: results[k]
            for k in (
                "status",
                "discovered_cell_count",
                "passing_cell_count",
                "blocker_count",
                "first_failed_condition",
                "global_max_length_upper",
                "global_min_margin_to_cap",
                "global_length_below_cap",
                "cross_cell_unique_source_ids",
                "cross_cell_collision_count",
                "cap_consistent_across_cells",
            )
        }
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0 if global_pass else 1

    out_dir.mkdir(parents=True, exist_ok=False)
    rj_path = out_dir / f"{OUT_ID}_RESULTS.json"
    rj_path.write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    sha = sha256_file(rj_path)
    sha_path = out_dir / f"{OUT_ID}_RESULTS.sha256"
    sha_path.write_text(f"{sha}  {rj_path.name}\n", encoding="utf-8")

    report_lines = [
        f"# EHP114 n=14 Global Atlas Integration Certificate",
        "",
        f"Experiment: `{OUT_ID}`",
        "",
        "## Verdict",
        "",
        f"- Status: `{status}`",
        f"- Discovered cell count: `{results['discovered_cell_count']}` of `64`",
        f"- Passing cell count: `{results['passing_cell_count']}` of `64`",
        f"- Blocker count: `{results['blocker_count']}`",
        f"- First failed condition: `{first_failed_condition}`",
        f"- Cap consistent across cells: `{cap_consistent}`",
        f"- Global max length upper: `{global_max_length}`",
        f"- Global min margin to cap: `{global_min_margin}`",
        f"- Exact length cap: `{EXACT_LENGTH_CAP}`",
        f"- Cross-cell unique source IDs: `{results['cross_cell_unique_source_ids']}`",
        f"- Cross-cell source ID collisions: `{results['cross_cell_collision_count']}`",
        "",
        "## Tightest Cells (smallest margin to cap)",
        "",
        "| Cell | Length upper | Margin to cap |",
        "|---|---:|---:|",
    ]
    for c in results["tightest_cells"]:
        report_lines.append(
            f"| `{c['cell_tag']}` | `{c['total_validated_length_upper']}` | `{c['margin_to_cap']}` |"
        )
    report_lines.extend(
        [
            "",
            "## Loosest Cells (largest margin to cap)",
            "",
            "| Cell | Length upper | Margin to cap |",
            "|---|---:|---:|",
        ]
    )
    for c in results["loosest_cells"]:
        report_lines.append(
            f"| `{c['cell_tag']}` | `{c['total_validated_length_upper']}` | `{c['margin_to_cap']}` |"
        )
    report_lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This certificate composes 64 already-emitted local-cell certificate "
            "packets into one global integration artifact. It does not recompute "
            "interval geometry. The pass condition is bookkeeping integrity: "
            "every per-cell packet must carry a passing local certificate, every "
            "cell artifact's own RESULTS.json SHA-256 must match its sidecar, every "
            "cell must use the same exact length cap, no source experiment ID may "
            "be claimed by two cells, and the maximum local length upper across "
            "the atlas must remain below the exact length cap.",
            "",
            "## Claim Ceiling",
            "",
            CLAIM_CEILING,
            "",
            "## Next Dependencies",
            "",
            f"- {results['next_dependency']}",
            "",
        ]
    )
    report_path = out_dir / f"{OUT_ID}_REPORT.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")

    print(f"WROTE {rj_path}")
    print(f"WROTE {sha_path}")
    print(f"WROTE {report_path}")
    print(f"STATUS={status}")
    print(f"GLOBAL_MAX_LENGTH_UPPER={global_max_length}")
    print(f"GLOBAL_MIN_MARGIN_TO_CAP={global_min_margin}")
    return 0 if global_pass else 1


if __name__ == "__main__":
    sys.exit(main())
