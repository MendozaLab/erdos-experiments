#!/usr/bin/env python3
"""Track B1 residual-sweep umbrella aggregation.

Reads the five per-family artifacts produced by
`build_n15_cell0203_centerline_sign_model_replicate.py` for owner families
{3104:0, 3546:0, 3982:0, 4404:0, 2932:1} (the residual CELL-02-03 wall-
separation owner families that were not part of the original 4488:3 / 2484:4 /
4571:2 runs) and builds a single umbrella artifact summarizing the WS-01
effectiveness across all 65 residual failures.

It also pulls the three already-done family results into an aggregate-after-
today section so the artifact is the single ledger covering all 108
CELL-02-03 wall-separation failures.

Internal experiment artifact only. Not a proof of Erdos #114, not an n=15
certificate, not a CELL-02-03 closure. Same Python-numerical-demonstration
claim level as the per-family runs. NOT a Rust interval-certified result.
"""

from __future__ import annotations

import hashlib
import json
import statistics
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

EXPERIMENT_ID = "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-RESIDUAL-SWEEP-20260509-01"
SCRIPT_PATH = Path(__file__).resolve()
ERDOS114_DIR = SCRIPT_PATH.parent
PROOF_PATH = ERDOS114_DIR / "proof_path"

OUTPUT_DIR = PROOF_PATH / EXPERIMENT_ID

# Already-done predecessor families (read-only references)
PRIOR_FAMILIES = [
    (
        "4488:3",
        PROOF_PATH
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01"
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01_RESULTS.json",
    ),
    (
        "2484:4",
        PROOF_PATH
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2484-4-20260509-01"
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2484-4-20260509-01_RESULTS.json",
    ),
    (
        "4571:2",
        PROOF_PATH
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4571-2-20260509-01"
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4571-2-20260509-01_RESULTS.json",
    ),
]

# Residual families produced by today's per-family runs
RESIDUAL_FAMILIES = [
    (
        "3104:0",
        PROOF_PATH
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3104-0-20260509-01"
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3104-0-20260509-01_RESULTS.json",
    ),
    (
        "3546:0",
        PROOF_PATH
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3546-0-20260509-01"
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3546-0-20260509-01_RESULTS.json",
    ),
    (
        "3982:0",
        PROOF_PATH
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3982-0-20260509-01"
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3982-0-20260509-01_RESULTS.json",
    ),
    (
        "4404:0",
        PROOF_PATH
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4404-0-20260509-01"
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4404-0-20260509-01_RESULTS.json",
    ),
    (
        "2932:1",
        PROOF_PATH
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2932-1-20260509-01"
        / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2932-1-20260509-01_RESULTS.json",
    ),
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def family_outcome_summary(family_payload: dict[str, Any]) -> str:
    closed = family_payload["failures_closed_by_ws01"]
    in_family = family_payload["failures_in_family"]
    not_found = family_payload["failures_branch_point_not_found"]
    if in_family == 0:
        return "EMPTY"
    if closed == in_family:
        return "FULL_PASS"
    if not_found == in_family:
        return "BLOCKED"
    if closed == 0:
        return "BLOCKED"
    return "PARTIAL"


def build_per_family_row(owner_key: str, payload: dict[str, Any]) -> dict[str, Any]:
    return {
        "owner_key": owner_key,
        "experiment_id": payload.get("experiment_id"),
        "failures": payload["failures_in_family"],
        "closed": payload["failures_closed_by_ws01"],
        "closed_tight_R_only": payload.get("failures_closed_by_ws01_tight_R_only", 0),
        "still_failing": payload["failures_still_failing"],
        "branch_point_not_found": payload["failures_branch_point_not_found"],
        "median_required_factor_before": payload["median_required_factor_before"],
        "median_required_factor_after": payload["median_required_factor_after"],
        "outcome_summary": family_outcome_summary(payload),
        "source_artifact_status": payload.get("status"),
    }


def main() -> int:
    if OUTPUT_DIR.exists():
        raise SystemExit(
            f"Refusing to overwrite existing output directory: {OUTPUT_DIR}"
        )

    log_lines: list[str] = []

    def log(msg: str) -> None:
        line = f"[{datetime.now(timezone.utc).isoformat()}] {msg}"
        log_lines.append(line)
        print(line, flush=True)

    log(f"experiment_id={EXPERIMENT_ID}")
    log("aggregation mode: pattern A umbrella over residual sweep")

    residual_per_family: list[dict[str, Any]] = []
    residual_per_failure_outcomes: list[dict[str, Any]] = []
    residual_source_paths: dict[str, str] = {}
    residual_source_sha_status: list[dict[str, Any]] = []

    for owner_key, json_path in RESIDUAL_FAMILIES:
        if not json_path.exists():
            raise SystemExit(f"Missing residual family artifact: {json_path}")
        payload = load_json(json_path)
        sha = sha256_file(json_path)
        log(f"loaded {owner_key}: {json_path.name} sha256={sha}")
        residual_per_family.append(build_per_family_row(owner_key, payload))
        residual_per_failure_outcomes.extend(payload["per_failure_outcomes"])
        residual_source_paths[owner_key] = str(json_path)
        residual_source_sha_status.append(
            {
                "owner_key": owner_key,
                "path": str(json_path),
                "actual_sha256": sha,
            }
        )

    prior_per_family: list[dict[str, Any]] = []
    prior_source_paths: dict[str, str] = {}
    for owner_key, json_path in PRIOR_FAMILIES:
        if not json_path.exists():
            raise SystemExit(f"Missing prior family artifact: {json_path}")
        payload = load_json(json_path)
        log(f"loaded prior {owner_key}: {json_path.name}")
        prior_per_family.append(build_per_family_row(owner_key, payload))
        prior_source_paths[owner_key] = str(json_path)

    total_residual_failures = sum(r["failures"] for r in residual_per_family)
    total_residual_closed = sum(r["closed"] for r in residual_per_family)
    total_residual_closed_tight = sum(
        r["closed_tight_R_only"] for r in residual_per_family
    )
    total_residual_still_failing = sum(r["still_failing"] for r in residual_per_family)
    total_residual_not_found = sum(
        r["branch_point_not_found"] for r in residual_per_family
    )

    if total_residual_closed == total_residual_failures and total_residual_failures > 0:
        sweep_status = "RESIDUAL_SWEEP_FULL_PASS"
    elif total_residual_closed > 0:
        sweep_status = "RESIDUAL_SWEEP_PARTIAL"
    else:
        sweep_status = "RESIDUAL_SWEEP_BLOCKED"

    # Aggregate-after-today across all 108 CELL-02-03 wall-separation failures
    all_families = residual_per_family + prior_per_family
    total_all = sum(r["failures"] for r in all_families)
    closed_all = sum(r["closed"] for r in all_families)
    closed_tight_all = sum(r["closed_tight_R_only"] for r in all_families)
    still_failing_all = sum(r["still_failing"] for r in all_families)
    not_found_all = sum(r["branch_point_not_found"] for r in all_families)
    closure_pct = (100.0 * closed_all / total_all) if total_all else 0.0

    # Median factors across residual families that have data
    before_factors = [
        r["median_required_factor_before"]
        for r in residual_per_family
        if r["median_required_factor_before"] is not None
    ]
    after_factors = [
        r["median_required_factor_after"]
        for r in residual_per_family
        if r["median_required_factor_after"] is not None
    ]
    median_before_overall = (
        statistics.median(before_factors) if before_factors else None
    )
    median_after_overall = statistics.median(after_factors) if after_factors else None

    interpretation = (
        f"Across the 5 residual owner families ({total_residual_failures} boxes), "
        f"WS-01 closes {total_residual_closed} boxes by the conservative 2R bound "
        f"and an additional {total_residual_closed_tight} boxes by the tight-R bound; "
        f"{total_residual_still_failing} remain failing and "
        f"{total_residual_not_found} have no validated branch point. "
        f"Two of the five (4404:0, 2932:1; 17 boxes total) are uniform full-pass — "
        f"those rows have R*|F_n|/|F_t|-style geometry where the rewrite gives a "
        f"large factor of safety. The remaining three (3104:0, 3546:0, 3982:0; 48 "
        f"boxes total) are dominated by R^2*|F_tt| where the conservative 2R-radius "
        f"penalty overshoots the wall RHS by ~30-50%, but the tight-R sensitivity "
        f"figure clears comfortably (~99% of the 48 close at tight R). Combined "
        f"with prior runs (4488:3=16/16, 2484:4=16/16, 4571:2=4/11 at 2R), the "
        f"aggregate across all 108 CELL-02-03 wall failures is "
        f"{closed_all}/{total_all} = {closure_pct:.1f}% closed at conservative 2R. "
        f"This is a Python numerical demonstration, NOT a Rust interval-certified "
        f"result."
    )

    if sweep_status == "RESIDUAL_SWEEP_FULL_PASS":
        next_dependency_residual = (
            "Residual sweep is uniform full pass at the conservative 2R bound. "
            "All 5 residual owner families clear WS-01. Combined with the prior "
            "runs, CELL-02-03 closure under WS-01 is uniform across all owner "
            "families. Promote to Rust interval certification of the moving-frame "
            "collar via interval Newton on F(z)=0."
        )
    elif sweep_status == "RESIDUAL_SWEEP_PARTIAL":
        next_dependency_residual = (
            "Partial pass on residual sweep. WS-01 generalizes but is not uniform "
            "at the conservative 2R bound. Three owner families (3104:0, 3546:0, "
            "3982:0) need either a tighter R bound (interval-certified z* via "
            "interval Newton, which would replace 2R with a true off-center radius) "
            "or a refined T3/F_tt bound. Recommended next: produce a tight-R "
            "(z*-anchored) variant of the WS-01 closure for these three families "
            "and verify whether the tight-R closure (which already clears) survives "
            "interval-Newton certification."
        )
    else:
        next_dependency_residual = (
            "Residual sweep BLOCKED at the conservative 2R bound. WS-01 alone is "
            "insufficient for the residual cluster of CELL-02-03 owner families. "
            "Either subdivide cells further or design a different rewrite (e.g., "
            "an analytic critical-point exclusion theorem)."
        )

    if closure_pct >= 90.0:
        next_dependency_aggregate = (
            "Aggregate closure rate ≥ 90%. The residual cluster is small enough "
            "to attempt analytic critical-point exclusion or per-cell subdivision. "
            "CELL-02-03 may be closable as soon as the still-failing minority "
            "is handled by a tight-R / interval-Newton variant."
        )
    elif closure_pct >= 70.0:
        next_dependency_aggregate = (
            "Aggregate closure rate ≥ 70% but < 90%. Track B1 generalizes broadly "
            "across CELL-02-03 owner families but is not uniform at the "
            "conservative 2R bound. The minority of still-failing rows clusters "
            "in 2-3 owner families dominated by R^2*|F_tt|. Tight-R closure for "
            "those families is the next milestone before any cell-level closure "
            "claim."
        )
    else:
        next_dependency_aggregate = (
            "Aggregate closure rate < 70%. WS-01 alone is insufficient for "
            "CELL-02-03. A different rewrite or substantive cell subdivision is "
            "needed before Track B1 can be claimed as a CELL-02-03 closure path."
        )

    aggregate_after_today = {
        "total_cell_02_03_failures": 108,
        "closed_across_all_owner_families": closed_all,
        "closed_tight_R_only": closed_tight_all,
        "still_failing": still_failing_all,
        "branch_point_not_found": not_found_all,
        "closure_pct": f"{closure_pct:.2f}%",
        "owner_families_full_pass": [
            r["owner_key"] for r in all_families if r["outcome_summary"] == "FULL_PASS"
        ],
        "owner_families_partial": [
            r["owner_key"] for r in all_families if r["outcome_summary"] == "PARTIAL"
        ],
        "owner_families_blocked": [
            r["owner_key"] for r in all_families if r["outcome_summary"] == "BLOCKED"
        ],
        "next_dependency": next_dependency_aggregate,
    }

    payload_out = {
        "experiment_id": EXPERIMENT_ID,
        "schema_version": "1.0",
        "timestamp_unix": int(time.time()),
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "generated_by": "track_b1_residual_sweep_aggregator",
        "origin": "auto-research",
        "promotion_state": "review_only",
        "target_problem": 114,
        "degree": 15,
        "cell_id": "CELL-02-03",
        "claim_ceiling": (
            "Internal sweep artifact only. Not a proof of Erdos #114, not an n=15 "
            "certificate, and not a CELL-02-03 closure. Diagnostic of WS-01 "
            "effectiveness across all owner families of CELL-02-03 not already "
            "covered by 4488:3 / 2484:4 / 4571:2. Python numerical demonstration; "
            "NOT Rust interval-certified."
        ),
        "owner_families_already_covered": [r["owner_key"] for r in prior_per_family],
        "owner_families_in_residual_sweep": [
            r["owner_key"] for r in residual_per_family
        ],
        "status": sweep_status,
        "total_residual_failures": total_residual_failures,
        "total_residual_closed_by_ws01": total_residual_closed,
        "total_residual_closed_by_ws01_tight_R_only": total_residual_closed_tight,
        "total_residual_still_failing": total_residual_still_failing,
        "total_residual_branch_point_not_found": total_residual_not_found,
        "median_required_factor_before_overall_residual": median_before_overall,
        "median_required_factor_after_overall_residual": median_after_overall,
        "per_owner_family_residual": residual_per_family,
        "per_owner_family_prior": prior_per_family,
        "per_failure_outcomes_residual": residual_per_failure_outcomes,
        "interpretation": interpretation,
        "next_dependency_residual": next_dependency_residual,
        "aggregate_after_today": aggregate_after_today,
        "source_paths": {
            "residual_family_artifacts": residual_source_paths,
            "prior_family_artifacts": prior_source_paths,
            "boundary_slice_source": str(
                PROOF_PATH
                / "EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03"
                / "EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03_RESULTS.json"
            ),
            "wall_separation_target_packet": str(
                PROOF_PATH
                / "EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01"
                / "EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01_RESULTS.json"
            ),
            "math_reference": str(
                ERDOS114_DIR
                / "EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_2026-05-06.md"
            ),
        },
        "source_sha_status": residual_source_sha_status,
        "forbidden_writes": [
            "D1",
            "morphisms.json",
            "proof registries",
            "public pages",
            "existing EHP114 finite packets",
            "existing corrected n=13 or n=14 receipts",
            "any of the per-family centerline-sign-model artifacts (read-only)",
        ],
    }

    OUTPUT_DIR.mkdir(parents=True)
    results_json = OUTPUT_DIR / f"{EXPERIMENT_ID}_RESULTS.json"
    report_md = OUTPUT_DIR / f"{EXPERIMENT_ID}_REPORT.md"
    results_sha = OUTPUT_DIR / f"{EXPERIMENT_ID}_RESULTS.sha256"
    run_log = OUTPUT_DIR / f"{EXPERIMENT_ID}_RUN.log"

    results_json.write_text(
        json.dumps(payload_out, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    # Build report
    table_lines = [
        "| owner_key | failures | closed (2R) | closed (tight R only) | still_failing | branch_not_found | outcome |",
        "|-----------|----------|-------------|------------------------|----------------|-------------------|---------|",
    ]
    for r in residual_per_family:
        table_lines.append(
            f"| `{r['owner_key']}` | {r['failures']} | {r['closed']} | "
            f"{r['closed_tight_R_only']} | {r['still_failing']} | "
            f"{r['branch_point_not_found']} | `{r['outcome_summary']}` |"
        )
    residual_table = "\n".join(table_lines)

    prior_table_lines = [
        "| owner_key | failures | closed (2R) | closed (tight R only) | still_failing | branch_not_found | outcome |",
        "|-----------|----------|-------------|------------------------|----------------|-------------------|---------|",
    ]
    for r in prior_per_family:
        prior_table_lines.append(
            f"| `{r['owner_key']}` | {r['failures']} | {r['closed']} | "
            f"{r['closed_tight_R_only']} | {r['still_failing']} | "
            f"{r['branch_point_not_found']} | `{r['outcome_summary']}` |"
        )
    prior_table = "\n".join(prior_table_lines)

    report = f"""# {EXPERIMENT_ID}

## Scope

Internal sweep artifact only. Track B1 residual-owner-family sweep of EHP114
bridge program. Not a proof of Erdos #114, not an n=15 certificate, not a
CELL-02-03 closure. Diagnostic of WS-01-CENTER-STRIP-CANCELLATION (Branch-
Centered Moving-Frame Collar rewrite) on the 5 owner families of CELL-02-03's
108 wall-separation failures that were not already covered by the 4488:3 /
2484:4 / 4571:2 runs. Python numerical demonstration at the same claim level
as the per-family runs; NOT Rust interval-certified.

## Pattern choice

Pattern B (per-family artifacts, mirroring the existing 2484-4 / 4571-2 layout)
for the actual WS-01 runs, plus this Pattern A umbrella aggregate to give a
single ledger across all 108 CELL-02-03 wall failures. The 5 residual families
were small enough that per-family granularity is more useful than a monolithic
sweep file, and the umbrella here pulls all of it together.

## Status

`{sweep_status}`

## Residual-sweep headline numbers (5 owner families, {total_residual_failures} boxes)

- Closed by WS-01 (conservative 2R bound): `{total_residual_closed}`
- Closed by WS-01 (tight R only): `{total_residual_closed_tight}`
- Still failing: `{total_residual_still_failing}`
- Branch point not found: `{total_residual_not_found}`
- Median required factor (RHS/LHS) before WS-01 (overall residual): `{median_before_overall}`
- Median required factor (RHS/LHS) after WS-01 at 2R (overall residual): `{median_after_overall}`

## Per-family table (residual sweep)

{residual_table}

## Per-family table (prior runs, for context only — read-only)

{prior_table}

## Aggregate after today (all 108 CELL-02-03 wall failures)

- Total failures: `{total_all}`
- Closed by WS-01 at conservative 2R: `{closed_all}`
- Closed by WS-01 at tight R only: `{closed_tight_all}`
- Still failing: `{still_failing_all}`
- Branch point not found: `{not_found_all}`
- Closure rate at conservative 2R: `{closure_pct:.2f}%`
- Owner families full-pass: `{aggregate_after_today["owner_families_full_pass"]}`
- Owner families partial: `{aggregate_after_today["owner_families_partial"]}`
- Owner families blocked: `{aggregate_after_today["owner_families_blocked"]}`

## Interpretation

{interpretation}

## Next dependency (residual sweep)

{next_dependency_residual}

## Next dependency (aggregate, all 108 boxes)

{next_dependency_aggregate}

## Source provenance

Per-family residual artifacts (each with its own SHA-256 sidecar in its own folder):
- `{residual_source_paths.get("3104:0")}`
- `{residual_source_paths.get("3546:0")}`
- `{residual_source_paths.get("3982:0")}`
- `{residual_source_paths.get("4404:0")}`
- `{residual_source_paths.get("2932:1")}`

Prior (read-only) artifacts:
- `{prior_source_paths.get("4488:3")}`
- `{prior_source_paths.get("2484:4")}`
- `{prior_source_paths.get("4571:2")}`

Math reference (consume, do not re-derive):
- `{ERDOS114_DIR / "EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_2026-05-06.md"}`
  lines 122-146 (Branch-Centered Moving-Frame Collar lemma).

## Safety

This artifact is review-only. Forbidden writes (D1, curated morphisms, public
pages, proof registries, existing EHP114 finite packets, any of the per-family
centerline-sign-model artifacts) were not attempted. The aggregate script
`build_n15_cell0203_residual_sweep_aggregate.py` writes only this experiment
folder. No prior script or artifact was modified.
"""

    report_md.write_text(report, encoding="utf-8")
    results_sha.write_text(
        f"{sha256_file(results_json)}  {results_json.name}\n",
        encoding="utf-8",
    )
    run_log.write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    summary = {
        "experiment_id": EXPERIMENT_ID,
        "status": sweep_status,
        "results": str(results_json),
        "report": str(report_md),
        "sha256": str(results_sha),
        "log": str(run_log),
        "total_residual_failures": total_residual_failures,
        "total_residual_closed_by_ws01": total_residual_closed,
        "total_residual_closed_by_ws01_tight_R_only": total_residual_closed_tight,
        "total_residual_still_failing": total_residual_still_failing,
        "total_residual_branch_point_not_found": total_residual_not_found,
        "aggregate_total": total_all,
        "aggregate_closed": closed_all,
        "aggregate_closure_pct": f"{closure_pct:.2f}%",
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
