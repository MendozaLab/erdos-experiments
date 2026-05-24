#!/usr/bin/env python3
"""Baseline-selection diagnostic for Erdos #20 core closure.

Reads only saved #20 artifacts and writes only inside baseline_selection.
It does not enumerate new sunflower families, update D1, touch scorecards,
or write public-facing documents.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Any


EXPERIMENT_ID = "EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01"
OUT_DIR = Path(__file__).resolve().parent
ROOT = OUT_DIR.parents[0]

DEEPENING_RESULTS = ROOT / "experimental_deepening" / "EXP-MATH-ERDOS20-CORE-CLOSURE-DEEPENING-20260505-01_RESULTS.json"
DEEPENING_PACKET = ROOT / "experimental_deepening" / "ERDOS20_CORE_CLOSURE_DEEPENING_PACKET_2026-05-05.md"
FORMAL_TARGETS = ROOT / "formal_core_closure" / "ERDOS20_FORMAL_CORE_CLOSURE_TARGETS_2026-05-05.md"
LITERATURE_GATE = ROOT / "Q1_LITERATURE_GATE_2026-04-17.md"
APRIL_001 = ROOT / "EXP-MATH-ERDOS20-SUNFLOWER-001_REPORT_2026-04-17.md"
APRIL_002 = ROOT / "EXP-MATH-ERDOS20-SUNFLOWER-002_REPORT.md"

INPUTS = [DEEPENING_RESULTS, DEEPENING_PACKET, FORMAL_TARGETS, LITERATURE_GATE, APRIL_001, APRIL_002]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def score_sum(scores: dict[str, int]) -> int:
    return sum(scores.values())


def round_float(value: float | None, places: int = 6) -> float | None:
    if value is None or not math.isfinite(value):
        return None
    return round(value, places)


def extract_stability(deepening: dict[str, Any]) -> dict[str, Any]:
    diagnostics = deepening.get("diagnostics", [])
    fixed_w3_s2 = next(
        (row for row in diagnostics if row.get("series") == "w=3, fixed core size s=2"),
        None,
    )
    strongest_w3 = next(
        (row for row in diagnostics if row.get("series") == "w=3 strongest per-core"),
        None,
    )
    geometry_drift_count = sum(1 for row in diagnostics if row.get("classification") == "GEOMETRY_DRIFT")
    stable_count = sum(1 for row in diagnostics if row.get("classification") == "STABLE_BY_SLOPE_ONLY")
    return {
        "best_current_series": fixed_w3_s2,
        "strongest_w3_series": strongest_w3,
        "geometry_drift_count": geometry_drift_count,
        "stable_by_slope_count": stable_count,
        "diagnostic_count": len(diagnostics),
        "observations_count": len(deepening.get("observations", [])),
    }


def build_candidates(stability: dict[str, Any]) -> list[dict[str, Any]]:
    """Score candidates on fixed, auditable criteria.

    Scale: 0-2 per criterion. The point is not to certify truth; it is to pick
    the next experiment's numerator/baseline without confusing computability
    with Leg-4 meaning.
    """
    fixed_w3 = stability.get("best_current_series") or {}
    observed_scores = {
        "data_available_now": 2,
        "core_channel_specificity": 2,
        "predeclared_denominator": 0,
        "mendoza_floor_binding": 0,
        "ahs_literature_control": 0,
        "geometry_drift_resistance": 1 if fixed_w3.get("classification") == "STABLE_BY_SLOPE_ONLY" else 0,
        "next_run_operationality": 2,
    }
    floor_scores = {
        "data_available_now": 1,
        "core_channel_specificity": 2,
        "predeclared_denominator": 1,
        "mendoza_floor_binding": 2,
        "ahs_literature_control": 1,
        "geometry_drift_resistance": 1,
        "next_run_operationality": 1,
    }
    ahs_scores = {
        "data_available_now": 0,
        "core_channel_specificity": 1,
        "predeclared_denominator": 1,
        "mendoza_floor_binding": 0,
        "ahs_literature_control": 2,
        "geometry_drift_resistance": 1,
        "next_run_operationality": 0,
    }
    return [
        {
            "candidate": "observed_I_core_local",
            "quantity": "I_core_local(s,m)",
            "scores": observed_scores,
            "total_score": score_sum(observed_scores),
            "rank_note": "Best immediate numerator, but not a Leg-4 baseline because it lacks a denominator.",
        },
        {
            "candidate": "floor_normalized_I_core_local",
            "quantity": "I_core_local(s,m) / I_floor(s,m)",
            "scores": floor_scores,
            "total_score": score_sum(floor_scores),
            "rank_note": "Best Leg-4 target if I_floor is predeclared before the run; preserves the measured core channel.",
        },
        {
            "candidate": "ahs_normalized_closure_cost",
            "quantity": "I_core_local(s,m) / I_AHS(s,m) or closure excess over an AHS-style construction baseline",
            "scores": ahs_scores,
            "total_score": score_sum(ahs_scores),
            "rank_note": "Necessary external control, but not ready as the primary numerator because no AHS construction artifact exists here.",
        },
    ]


def report_markdown(results: dict[str, Any]) -> str:
    rows = []
    for candidate in results["candidate_ranking"]:
        rows.append(
            "| {candidate} | `{quantity}` | {score} | {note} |".format(
                candidate=candidate["candidate"],
                quantity=candidate["quantity"],
                score=candidate["total_score"],
                note=candidate["rank_note"],
            )
        )
    best = results["recommendation"]
    fixed = results["stability_summary"].get("best_current_series") or {}
    return "\n".join(
        [
            f"# {EXPERIMENT_ID} Report",
            "",
            "**Scope:** saved-artifact diagnostic only. No new enumeration, no D1, no scorecard, no git, no public docs.",
            "",
            f"**Recommendation:** `{best['candidate']}` as the next-run baseline, with observed `I_core_local(s,m)` as the numerator and AHS as a secondary external control.",
            "",
            "## Ranking",
            "",
            "| Candidate | Quantity | Score | Meaning |",
            "|---|---:|---:|---|",
            *rows,
            "",
            "## Current Signal Used",
            "",
            f"The best saved precursor remains `{fixed.get('series')}` with classification `{fixed.get('classification')}`, slope `{fixed.get('loglog_slope_vs_N')}`, and late-window CV `{fixed.get('late3_cv')}`.",
            "",
            "## Claim Ceiling",
            "",
            results["claim_ceiling"],
            "",
        ]
    )


def main() -> None:
    deepening = load_json(DEEPENING_RESULTS)
    stability = extract_stability(deepening)
    candidates = build_candidates(stability)
    leg4_sorted = sorted(
        candidates,
        key=lambda row: (
            row["scores"]["mendoza_floor_binding"],
            row["scores"]["core_channel_specificity"],
            row["scores"]["next_run_operationality"],
            row["total_score"],
        ),
        reverse=True,
    )
    results = {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #20 sunflower core closure",
        "scope": "baseline selection over saved May 5 core-closure artifacts and April literature/experiment reports",
        "guardrails": {
            "read_existing_artifacts_only": True,
            "no_new_enumeration": True,
            "no_d1": True,
            "no_scorecard": True,
            "no_git": True,
            "write_scope": "erdos-experiments/Erdos20/baseline_selection/",
        },
        "input_metadata": [
            {
                "path": str(path.relative_to(ROOT.parents[1])),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in INPUTS
        ],
        "stability_summary": stability,
        "candidate_ranking": sorted(candidates, key=lambda row: row["total_score"], reverse=True),
        "leg4_suitability_order": leg4_sorted,
        "recommendation": leg4_sorted[0],
        "recommended_next_run_quantity": {
            "primary": "I_core_local(s,m) / I_floor(s,m)",
            "numerator": "fixed-core-size local closure information I_core_local(s,m), starting with w=3, s=2, m near m* and across m* +/- delta",
            "baseline": "predeclared Mendoza-floor information denominator I_floor(s,m)",
            "secondary_control": "AHS-style construction-normalized closure cost, reported separately as a literature control rather than merged into the floor denominator",
        },
        "claim_ceiling": "A-axis A0; shadow signature, not universal law; no theorem, no lower-bound progress, no Leg-4 pass.",
    }

    result_path = OUT_DIR / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = OUT_DIR / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = OUT_DIR / f"{EXPERIMENT_ID}_RESULTS.sha256"

    result_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_path.write_text(report_markdown(results), encoding="utf-8")
    sha_path.write_text(f"{sha256_file(result_path)}  {result_path.name}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
