#!/usr/bin/env python3
"""Build a review-only wall-separation target packet for EHP114 n=15 CELL-02-03.

This is an analysis artifact over the latest adaptive boundary-slice run. It
does not certify degree 15, promote proof status, or mutate prior packets.
"""

from __future__ import annotations

import hashlib
import json
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01"
SOURCE_EXPERIMENT_ID = "EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03"
TARGET_PROBLEM = 114
DEGREE = 15
CELL_ID = "CELL-02-03"

SCRIPT_PATH = Path(__file__).resolve()
ERDOS114_DIR = SCRIPT_PATH.parent
PROOF_PATH = ERDOS114_DIR / "proof_path"
SOURCE_DIR = PROOF_PATH / SOURCE_EXPERIMENT_ID
SOURCE_JSON = SOURCE_DIR / f"{SOURCE_EXPERIMENT_ID}_RESULTS.json"
SOURCE_SHA = SOURCE_DIR / f"{SOURCE_EXPERIMENT_ID}_RESULTS.sha256"
OUTPUT_DIR = PROOF_PATH / EXPERIMENT_ID
RESULTS_JSON = OUTPUT_DIR / f"{EXPERIMENT_ID}_RESULTS.json"
REPORT_MD = OUTPUT_DIR / f"{EXPERIMENT_ID}_REPORT.md"
RESULTS_SHA = OUTPUT_DIR / f"{EXPERIMENT_ID}_RESULTS.sha256"
STORY_DIR = Path("/Users/kenbengoetxea/Downloads/Eratosthenes_Ahmes_Stories")
STORY_MD = STORY_DIR / f"AHMES_STORY_{EXPERIMENT_ID}.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def quantiles(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"min": None, "p25": None, "median": None, "p75": None, "max": None}
    ordered = sorted(values)

    def nearest(q: float) -> float:
        if len(ordered) == 1:
            return ordered[0]
        idx = round(q * (len(ordered) - 1))
        return ordered[idx]

    return {
        "min": ordered[0],
        "p25": nearest(0.25),
        "median": statistics.median(ordered),
        "p75": nearest(0.75),
        "max": ordered[-1],
    }


def safe_div(num: float, den: float) -> float | None:
    if den == 0:
        return None
    return num / den


def source_hash_from_sidecar() -> str | None:
    if not SOURCE_SHA.exists():
        return None
    text = SOURCE_SHA.read_text(encoding="utf-8").strip()
    return text.split()[0] if text else None


def flatten_wall_rows(source: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in source.get("remaining_unresolved", []):
        if item.get("reason") != "third_order_wall_separation_failed":
            continue
        ineq = item.get("inequality", {})
        lhs = float(ineq.get("wall_lhs", 0.0))
        rhs = float(ineq.get("wall_rhs", 0.0))
        center = float(ineq.get("center_strip_f_abs_upper", 0.0))
        remainder = float(ineq.get("normal_remainder_bound", 0.0))
        fn_lower = float(ineq.get("fn_abs_lower", 0.0))
        gap = lhs - rhs
        remaining_center_budget = rhs - remainder
        rows.append(
            {
                "source_index": item.get("source_index"),
                "source_ownership_key": item.get("source_ownership_key"),
                "split_path": item.get("split_path"),
                "wall_lhs": lhs,
                "wall_rhs": rhs,
                "gap": gap,
                "margin": float(ineq.get("margin", rhs - lhs)),
                "center_strip_f_abs_upper": center,
                "normal_remainder_bound": remainder,
                "fn_abs_lower": fn_lower,
                "normal_radius": float(ineq.get("normal_radius", 0.0)),
                "tangent_radius": float(ineq.get("tangent_radius", 0.0)),
                "center_strip_share_of_lhs": safe_div(center, lhs),
                "remainder_share_of_lhs": safe_div(remainder, lhs),
                "required_lhs_factor": safe_div(rhs, lhs),
                "required_center_strip_factor_after_remainder": (
                    safe_div(max(0.0, remaining_center_budget), center)
                ),
                "center_budget_negative": remaining_center_budget < 0.0,
            }
        )
    return rows


def flatten_critical_rows(source: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in source.get("remaining_unresolved", []):
        if item.get("reason") != "third_order_critical_exclusion_failed":
            continue
        ineq = item.get("inequality", {})
        rows.append(
            {
                "source_index": item.get("source_index"),
                "source_ownership_key": item.get("source_ownership_key"),
                "split_path": item.get("split_path"),
                "critical_margin": float(ineq.get("critical_exclusion_margin", 0.0)),
                "critical_point_upper_bound": float(ineq.get("critical_point_upper_bound", 0.0)),
                "fn_abs_lower": float(ineq.get("fn_abs_lower", 0.0)),
            }
        )
    return rows


def top_wall_groups(rows: list[dict[str, Any]], limit: int = 8) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[str(row["source_ownership_key"])].append(row)
    summaries = []
    for key, items in groups.items():
        worst = max(items, key=lambda r: r["gap"])
        summaries.append(
            {
                "source_ownership_key": key,
                "wall_failure_count": len(items),
                "max_gap": worst["gap"],
                "worst_source_index": worst["source_index"],
                "worst_split_path": worst["split_path"],
                "worst_lhs": worst["wall_lhs"],
                "worst_rhs": worst["wall_rhs"],
                "median_required_lhs_factor": statistics.median(
                    [r["required_lhs_factor"] for r in items if r["required_lhs_factor"] is not None]
                ),
                "median_center_strip_share": statistics.median(
                    [
                        r["center_strip_share_of_lhs"]
                        for r in items
                        if r["center_strip_share_of_lhs"] is not None
                    ]
                ),
            }
        )
    summaries.sort(key=lambda row: (row["max_gap"], row["wall_failure_count"]), reverse=True)
    return summaries[:limit]


def build_results(source: dict[str, Any], source_hash: str) -> dict[str, Any]:
    wall_rows = flatten_wall_rows(source)
    critical_rows = flatten_critical_rows(source)
    reasons = Counter(item.get("reason", "unknown") for item in source.get("remaining_unresolved", []))

    required_lhs_factors = [r["required_lhs_factor"] for r in wall_rows if r["required_lhs_factor"] is not None]
    center_shares = [
        r["center_strip_share_of_lhs"] for r in wall_rows if r["center_strip_share_of_lhs"] is not None
    ]
    remainder_shares = [
        r["remainder_share_of_lhs"] for r in wall_rows if r["remainder_share_of_lhs"] is not None
    ]
    gaps = [r["gap"] for r in wall_rows]
    center_budget_negative_count = sum(1 for r in wall_rows if r["center_budget_negative"])
    top_groups = top_wall_groups(wall_rows)

    median_center_share = statistics.median(center_shares) if center_shares else None
    dominant_mode = (
        "CENTER_STRIP_DOMINATED_WALL_SEPARATION"
        if median_center_share is not None and median_center_share >= 0.95
        else "MIXED_WALL_SEPARATION"
    )

    generated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    results = {
        "experiment_id": EXPERIMENT_ID,
        "generated_at": generated_at,
        "generated_by": "erdos_atlas_autoresearch_librarian",
        "origin": "auto-research",
        "promotion_state": "review_only",
        "persona": "Eratosthenes of Cyrene",
        "short_name": "Eratosthenes",
        "scribe": "Ahmes",
        "story_writer": "Ahmes",
        "mathematician_persona": "David Hilbert",
        "target_problem": TARGET_PROBLEM,
        "degree": DEGREE,
        "cell_id": CELL_ID,
        "source_experiment_id": SOURCE_EXPERIMENT_ID,
        "source_results_path": str(SOURCE_JSON),
        "source_results_sha256": source_hash,
        "claim_ceiling": (
            "This packet diagnoses a boundary-slice obstruction and proposes "
            "review-only next experiments. It is not a degree-15 certification."
        ),
        "final_verdict": "WALL_SEPARATION_TARGET_READY",
        "source_summary": {
            "source_final_verdict": source.get("final_verdict"),
            "source_processed_region_count": source.get("processed_region_count"),
            "source_remaining_unresolved_region_count": source.get("remaining_unresolved_region_count"),
            "unresolved_reason_counts": dict(sorted(reasons.items())),
        },
        "wall_separation_analysis": {
            "wall_failure_count": len(wall_rows),
            "critical_exclusion_failure_count": len(critical_rows),
            "dominant_failure_mode": dominant_mode,
            "gap_lhs_minus_rhs": quantiles(gaps),
            "required_lhs_factor_rhs_over_lhs": quantiles(required_lhs_factors),
            "center_strip_share_of_lhs": quantiles(center_shares),
            "normal_remainder_share_of_lhs": quantiles(remainder_shares),
            "center_budget_negative_count": center_budget_negative_count,
            "interpretation": (
                "The wall inequality is failing mostly because the center-strip absolute bound is far "
                "larger than the admissible wall threshold; the normal remainder is not the main term."
            ),
        },
        "top_wall_owner_groups": top_groups,
        "representative_wall_failures": sorted(wall_rows, key=lambda r: r["gap"], reverse=True)[:12],
        "target_inequalities": [
            {
                "target_id": "WS-01-CENTER-STRIP-CANCELLATION",
                "status": "FORMAL_TARGET_READY",
                "question": (
                    "Can the centerline term be bounded by signed or cancellation-aware structure "
                    "instead of the raw absolute upper bound?"
                ),
                "acceptance_rule": (
                    "On the 108 current wall failures, replace the raw center-strip bound with a "
                    "certified bound that brings wall_lhs below wall_rhs for at least one dominant "
                    "owner group without weakening interval containment."
                ),
                "evidence_from_packet": (
                    "Median center-strip share of wall_lhs is recorded above; if it is near one, "
                    "the blocker is the centerline estimate rather than collar remainder."
                ),
            },
            {
                "target_id": "WS-02-SIGNED-WALL-ENDPOINTS",
                "status": "NEEDS_EVIDENCE",
                "question": (
                    "Can the normal-wall test use signed endpoint intervals or opposite-side "
                    "separation instead of sup_r |F(0,r)|?"
                ),
                "acceptance_rule": (
                    "Run a new immutable slice experiment that records signed wall endpoints for "
                    "the current wall-failure boxes and converts a nonzero subset to excluded or "
                    "certified status."
                ),
                "failure_boundary": (
                    "If signed endpoint intervals still overlap zero with the same scale as the raw "
                    "center strip, the route needs a different local model."
                ),
            },
            {
                "target_id": "WS-03-DOMINANT-OWNER-LOCAL-MODEL",
                "status": "SCOUT",
                "question": (
                    "Can the worst owner groups be explained by a local chart artifact, a root "
                    "near-wall geometry, or a genuine collar obstruction?"
                ),
                "acceptance_rule": (
                    "Build a focused packet for the largest-gap source_ownership_key and report "
                    "whether a chart change, deeper split, or analytic bound reduces the max gap."
                ),
                "top_owner_groups": [row["source_ownership_key"] for row in top_groups[:3]],
            },
        ],
        "recommended_next_experiments": [
            {
                "experiment_id": "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260508-01",
                "purpose": "Replace raw center-strip absolute bound with a signed or cancellation-aware centerline model.",
                "go_rule": "Proceed if it reduces wall_lhs by the factor recorded in required_lhs_factor for any top owner group.",
                "stop_rule": "Stop if signed intervals remain centered at zero with no owner-specific separation.",
            },
            {
                "experiment_id": "EXP-MATH-EHP114-N15-CELL-02-03-WALL-ENDPOINT-SPLIT-20260508-01",
                "purpose": "Record direct signed endpoint separation on normal walls for the 108 current wall failures.",
                "go_rule": "Proceed to another adaptive slice if endpoint signs separate on a nonzero subset.",
                "stop_rule": "Do not launch full n=15 if endpoint signs remain unresolved across dominant groups.",
            },
            {
                "experiment_id": "EXP-MATH-EHP114-N15-CELL-02-03-OWNER-LOCAL-MODEL-20260508-01",
                "purpose": "Analyze the worst source_ownership_key group before broadening to CELL-07-03 or CELL-06-03.",
                "go_rule": "Proceed if the local model identifies a reusable reduction rule.",
                "stop_rule": "If the local model is owner-specific and non-transferable, mark this route as needing analytic reduction.",
            },
        ],
        "forbidden_writes": [
            "D1",
            "morphisms.json",
            "proof registries",
            "public pages",
            "existing EHP114 finite packets",
            "existing corrected n=13 or n=14 receipts",
        ],
    }
    return results


def build_report(results: dict[str, Any]) -> str:
    analysis = results["wall_separation_analysis"]
    required = analysis["required_lhs_factor_rhs_over_lhs"]
    center = analysis["center_strip_share_of_lhs"]
    remainder = analysis["normal_remainder_share_of_lhs"]
    top = results["top_wall_owner_groups"][:3]

    top_lines = "\n".join(
        f"- `{row['source_ownership_key']}`: {row['wall_failure_count']} failures, "
        f"max gap {row['max_gap']:.6g}, median required factor {row['median_required_lhs_factor']:.6g}"
        for row in top
    )
    next_lines = "\n".join(
        f"- `{row['experiment_id']}`: {row['purpose']}" for row in results["recommended_next_experiments"]
    )

    return f"""# {EXPERIMENT_ID}

## Meaning

This packet turns the latest `CELL-02-03` adaptive slice into a narrow target:
the blocker is now wall separation, and the wall inequality is dominated by the
center-strip absolute bound rather than by the third-order collar remainder.

This is review-only. It is not a degree-15 certification and it does not change
any existing Erdős #114 packet.

## Source

- Source experiment: `{results['source_experiment_id']}`
- Source SHA-256: `{results['source_results_sha256']}`
- Source verdict: `{results['source_summary']['source_final_verdict']}`
- Remaining unresolved reasons: `{results['source_summary']['unresolved_reason_counts']}`

## Wall-Separation Diagnosis

- Wall failures: `{analysis['wall_failure_count']}`
- Critical-exclusion failures still present: `{analysis['critical_exclusion_failure_count']}`
- Dominant mode: `{analysis['dominant_failure_mode']}`
- Required lhs factor, median: `{required['median']}`
- Center-strip share of lhs, median: `{center['median']}`
- Normal-remainder share of lhs, median: `{remainder['median']}`

The practical read is simple: deeper interval subdivision helped move the first
blocker from critical exclusion to wall separation, but the present wall test is
too conservative because it pays for `sup_r |F(0,r)|` on the center strip. The
normal remainder is small. The next useful experiment is therefore not broad
compute; it is a signed/cancellation-aware wall test.

## Top Owner Groups

{top_lines}

## Target Inequalities

1. `WS-01-CENTER-STRIP-CANCELLATION`: replace the raw absolute center-strip
   bound with a certified signed or cancellation-aware bound.
2. `WS-02-SIGNED-WALL-ENDPOINTS`: record signed normal-wall endpoint intervals
   for the current wall-failure boxes.
3. `WS-03-DOMINANT-OWNER-LOCAL-MODEL`: explain the largest-gap owner group
   before broadening this method to other cells.

## Recommended Next Experiments

{next_lines}

## Safety

Forbidden writes remain D1, curated morphisms, public pages, proof registries,
and existing #114 finite packets. This artifact is an auto-research target
packet only.
"""


def build_story(results: dict[str, Any]) -> str:
    analysis = results["wall_separation_analysis"]
    return f"""# Ahmes Story: {EXPERIMENT_ID}

Ahmes records that Eratosthenes returned to the hard `CELL-02-03` boundary
slice not to widen the search, but to read the obstruction. The refined slice
had already shown that blind subdivision was not the main answer. The first
wall now stood in a precise place: the center strip was too large under the
absolute-value estimate, while the collar remainder was small.

So the next research move is no longer "run more boxes." It is to ask a sharper
mathematical question. Does the centerline have a sign, cancellation, or local
model that the current interval test cannot see? If yes, the same boxes may
start to separate. If no, this route needs a new reduction before degree 15 is
worth a full launch.

## Receipts

- Run ID: `{EXPERIMENT_ID}`
- Target: Erdős #114, degree 15, `CELL-02-03`
- Source run: `{results['source_experiment_id']}`
- Source SHA-256: `{results['source_results_sha256']}`
- Wall failures: `{analysis['wall_failure_count']}`
- Critical-exclusion failures: `{analysis['critical_exclusion_failure_count']}`
- Dominant failure mode: `{analysis['dominant_failure_mode']}`
- Median required lhs factor: `{analysis['required_lhs_factor_rhs_over_lhs']['median']}`
- Median center-strip share of lhs: `{analysis['center_strip_share_of_lhs']['median']}`
- Final verdict: `{results['final_verdict']}`

Review-only warning: this story is a companion reading artifact. The source of
truth is the run folder and its SHA-256 sidecar. No proof status, curated edge,
D1 row, public page, or prior finite packet was changed.
"""


def main() -> int:
    if OUTPUT_DIR.exists():
        raise SystemExit(f"Refusing to overwrite existing output directory: {OUTPUT_DIR}")
    if STORY_MD.exists():
        raise SystemExit(f"Refusing to overwrite existing story file: {STORY_MD}")
    if not SOURCE_JSON.exists():
        raise SystemExit(f"Missing source results: {SOURCE_JSON}")

    computed_source_hash = sha256_file(SOURCE_JSON)
    recorded_source_hash = source_hash_from_sidecar()
    if recorded_source_hash and recorded_source_hash != computed_source_hash:
        raise SystemExit(
            "Source SHA sidecar mismatch: "
            f"recorded={recorded_source_hash} computed={computed_source_hash}"
        )

    source = load_json(SOURCE_JSON)
    results = build_results(source, computed_source_hash)

    OUTPUT_DIR.mkdir(parents=True)
    STORY_DIR.mkdir(parents=True, exist_ok=True)
    write_json(RESULTS_JSON, results)
    REPORT_MD.write_text(build_report(results), encoding="utf-8")
    RESULTS_SHA.write_text(f"{sha256_file(RESULTS_JSON)}  {RESULTS_JSON.name}\n", encoding="utf-8")
    STORY_MD.write_text(build_story(results), encoding="utf-8")

    print(json.dumps({
        "experiment_id": EXPERIMENT_ID,
        "results": str(RESULTS_JSON),
        "report": str(REPORT_MD),
        "sha256": str(RESULTS_SHA),
        "story": str(STORY_MD),
        "final_verdict": results["final_verdict"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
