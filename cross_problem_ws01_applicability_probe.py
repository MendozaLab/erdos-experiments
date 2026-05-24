#!/usr/bin/env python3
"""WS-01 applicability probe -- problem-agnostic version of the
Branch-Centered Moving-Frame Collar diagnostic from Erdos #114.

This is Track B2 of the Tao-gap-closing plan. It generalizes the
problem-hardcoded
  Erdos114/build_n15_cell0203_centerline_sign_model.py
into a probe that accepts an Erdos problem ID and a boundary-slice JSON
in the same schema, runs the WS-01 algorithm (Branch-Centered Moving-Frame
Collar wall-separation rewrite) problem-agnostically, and outputs a verdict.

Reference implementation (read-only, not modified by this script):
  Erdos114/build_n15_cell0203_centerline_sign_model.py
Wall-separation primitives reference (read-only):
  Erdos114/build_n15_cell0203_wall_separation_target.py
Math reference for the rewrite:
  Erdos114/EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_2026-05-06.md
  (Branch-Centered Moving-Frame Collar lemma, lines 122-146)

INPUT SCHEMA (per failure-box, taken from the actual #114 boundary-slice
RESULTS.json; the canonical mapping the probe uses):

  item['reason']                       == 'third_order_wall_separation_failed'
  item['source_ownership_key']         (string; e.g. '4488:3') -- optional filter
  item['source_index']                 (int; for reporting)
  item['split_path']                   (string; for reporting)
  item['x_interval'], ['y_interval']   (list[float,float]; for reporting)
  item['tangent'], ['normal']          (dict with 'x','y'; for reporting)
  item['f_interval']                   {'lo': float, 'hi': float}
                                       interval enclosure of F over the box
  item['fn_full_first_order_interval'] {'lo': float, 'hi': float} (recorded)
  item['ft_full_first_order_interval'] {'lo': float, 'hi': float} (recorded)
  item['inequality']                   dict with:
    'center_fn_abs_lower'              float  -- lower bound on |F_n| at box center
    'tangent_radius_upper'             float  -- the R parameter
    'f_tt_abs_upper'                   float  -- upper bound on |F_tt|
    'third_directional_upper'          float  -- upper bound on T_3
    'normal_remainder_bound'           float  -- R_n
    'wall_lhs'                         float  -- old wall LHS
    'wall_rhs'                         float  -- S * lower(|F_n|)
    'center_strip_f_abs_upper'         float  -- old C0 hull (for reporting)

NAME MAPPING (this probe's vocabulary -> actual schema):
  F_iv_lower            <- f_interval['lo']
  F_iv_upper            <- f_interval['hi']
  F_tt_abs_upper        <- inequality['f_tt_abs_upper']
  t3_upper              <- inequality['third_directional_upper']
  r_n_upper             <- inequality['normal_remainder_bound']
  s_factor * f_n_abs_lower  <- inequality['wall_rhs']
  f_n_abs_lower         <- inequality['center_fn_abs_lower']
  tangent_radius_upper  <- inequality['tangent_radius_upper']  (the R)
  wall_lhs_before       <- inequality['wall_lhs']

VERDICT TAXONOMY (per-box outcome):
  BRANCH_EXISTS                     IVT validated (F_iv straddles 0)
  BRANCH_POINT_NOT_FOUND            F_iv does NOT straddle 0 (or grad vanishes)
  BRANCH_POINT_NOT_WELL_DEFINED     f_n_abs_lower <= 0 (gradient near-zero)
  CLOSED_BY_WS01                    new_LHS at conservative 2R < wall_rhs
  CLOSED_BY_WS01_TIGHT_R_ONLY       new_LHS at tight R < wall_rhs only
  STILL_FAILING_AFTER_WS01          even tight-R bound fails

VERDICT TAXONOMY (aggregate):
  WS01_APPLICABLE                   branch validation rate >= 0.95 AND
                                    closure rate at tight R >= 0.50
  WS01_PARTIAL                      branch validation rate >= 0.95 AND
                                    closure rate at tight R >= 0.10
  WS01_INAPPLICABLE_STRUCTURAL      branch validation rate <  0.95
                                    (architecture does not apply at all)
  WS01_INAPPLICABLE_QUANTITATIVE    branch validation rate >= 0.95 but
                                    closure rate at tight R <  0.10
                                    (architecture applies but constants
                                    do not favor closure)

CLAIM CEILING: this probe is a DIAGNOSTIC, NOT a proof. It only reports
whether the WS-01 architecture appears applicable; it does not certify
the wall separation, the boundary slice, or the underlying Erdos problem.

REGRESSION TEST: invoke with --regression. Loads the canonical #114
CELL-02-03 owner-family-4488:3 16 boxes and verifies 16/16 close at the
conservative 2R bound, matching the 2026-05-09 reference run.

USAGE:
  # Regression test (no other args needed; exits non-zero on failure):
  python3 cross_problem_ws01_applicability_probe.py --regression

  # Probe a new problem:
  python3 cross_problem_ws01_applicability_probe.py \
      --problem-id 1038 \
      --input-json path/to/boundary_slice.json \
      --output-dir path/to/output_dir \
      [--owner-family-filter 4488:3]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


# -----------------------------------------------------------------------------
# Reference paths for regression test (read-only).
# -----------------------------------------------------------------------------

SCRIPT_PATH = Path(__file__).resolve()
ERDOS_EXP_DIR = SCRIPT_PATH.parent
ERDOS114_DIR = ERDOS_EXP_DIR / "Erdos114"
ERDOS114_PROOF_PATH = ERDOS114_DIR / "proof_path"

REGRESSION_BOUNDARY_SLICE_JSON = (
    ERDOS114_PROOF_PATH
    / "EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03"
    / "EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03_RESULTS.json"
)
REGRESSION_OWNER_FAMILY = "4488:3"
REGRESSION_EXPECTED_TOTAL = 16
REGRESSION_EXPECTED_CLOSED_2R = 16
REGRESSION_EXPECTED_CLOSED_TIGHT_ONLY = 0
REGRESSION_EXPECTED_BRANCH_NOT_FOUND = 0


# -----------------------------------------------------------------------------
# I/O helpers (matching the reference implementation's conventions).
# -----------------------------------------------------------------------------

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
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


# -----------------------------------------------------------------------------
# Schema validation.
# -----------------------------------------------------------------------------

REQUIRED_TOP_FIELDS = (
    "f_interval",
    "inequality",
    "reason",
)

REQUIRED_INEQ_FIELDS = (
    "center_fn_abs_lower",
    "tangent_radius_upper",
    "f_tt_abs_upper",
    "third_directional_upper",
    "normal_remainder_bound",
    "wall_lhs",
    "wall_rhs",
)


def validate_box_schema(item: dict[str, Any], idx: int) -> list[str]:
    """Return list of schema-error strings; empty list means valid."""
    errors: list[str] = []
    for k in REQUIRED_TOP_FIELDS:
        if k not in item:
            errors.append(f"box[{idx}] missing top-level field '{k}'")
    if "f_interval" in item:
        fiv = item["f_interval"]
        if not isinstance(fiv, dict) or "lo" not in fiv or "hi" not in fiv:
            errors.append(f"box[{idx}] f_interval must have 'lo' and 'hi'")
    if "inequality" in item:
        ineq = item["inequality"]
        if not isinstance(ineq, dict):
            errors.append(f"box[{idx}] inequality must be a dict")
        else:
            for k in REQUIRED_INEQ_FIELDS:
                if k not in ineq:
                    errors.append(f"box[{idx}] inequality missing field '{k}'")
    return errors


# -----------------------------------------------------------------------------
# Core WS-01 classification (problem-agnostic; algorithm copied from
# build_n15_cell0203_centerline_sign_model.classify_box and de-hardcoded).
# -----------------------------------------------------------------------------

def classify_box(item: dict[str, Any]) -> dict[str, Any]:
    """Apply WS-01 (Branch-Centered Moving-Frame Collar) rewrite to one box.

    Algorithm (verbatim equivalent to the reference implementation):
      Branch-point validation:
        (a) f_interval straddles zero  (IVT => F has a zero in the box)
        (b) center_fn_abs_lower > 0    (grad F nonvanishing => smooth curve)
        BRANCH_POINT_NOT_WELL_DEFINED if (b) fails outright
        BRANCH_POINT_NOT_FOUND otherwise if (a) fails or (a) AND (b) both fail
      Conservative wall LHS at 2R (worst case z*-off-center):
        new_LHS_2R = 0.5 * f_tt_up * (2R)^2 + (T3/6) * (2R)^3 + R_n
      Tight wall LHS at R (best case, z* near box center):
        new_LHS_R  = 0.5 * f_tt_up *  R^2  + (T3/6) *  R^3  + R_n
      RHS unchanged:
        wall_rhs = inequality['wall_rhs']
      Outcome:
        CLOSED_BY_WS01                if new_LHS_2R < wall_rhs
        CLOSED_BY_WS01_TIGHT_R_ONLY   else if new_LHS_R < wall_rhs
        STILL_FAILING_AFTER_WS01      otherwise
    """
    ineq = item["inequality"]
    f_iv = item["f_interval"]

    f_lo = float(f_iv["lo"])
    f_hi = float(f_iv["hi"])
    fn_center_lo = float(ineq["center_fn_abs_lower"])
    R = float(ineq["tangent_radius_upper"])
    f_tt_up = float(ineq["f_tt_abs_upper"])
    T3 = float(ineq["third_directional_upper"])
    Rn = float(ineq["normal_remainder_bound"])
    wall_lhs_before = float(ineq["wall_lhs"])
    wall_rhs = float(ineq["wall_rhs"])
    C0_old = float(ineq.get("center_strip_f_abs_upper", 0.0))

    # Branch-point existence test (interval-arithmetic only).
    f_iv_contains_zero = (f_lo <= 0.0 <= f_hi)
    grad_well_defined = fn_center_lo > 0.0

    # Conservative-2R vs tight-R Taylor bounds (LINEAR-IN-R term gone).
    new_lhs_R = 0.5 * f_tt_up * R * R + (T3 / 6.0) * R * R * R + Rn
    new_lhs_2R = (
        0.5 * f_tt_up * (2.0 * R) ** 2 + (T3 / 6.0) * (2.0 * R) ** 3 + Rn
    )

    if not grad_well_defined:
        outcome = "BRANCH_POINT_NOT_WELL_DEFINED"
        branch_validated = False
    elif not f_iv_contains_zero:
        outcome = "BRANCH_POINT_NOT_FOUND"
        branch_validated = False
    else:
        branch_validated = True
        if new_lhs_2R < wall_rhs:
            outcome = "CLOSED_BY_WS01"
        elif new_lhs_R < wall_rhs:
            outcome = "CLOSED_BY_WS01_TIGHT_R_ONLY"
        else:
            outcome = "STILL_FAILING_AFTER_WS01"

    required_factor_before = (
        wall_rhs / wall_lhs_before if wall_lhs_before > 0 else None
    )
    required_factor_after_2R = (
        wall_rhs / new_lhs_2R if new_lhs_2R > 0 else None
    )
    required_factor_after_R = (
        wall_rhs / new_lhs_R if new_lhs_R > 0 else None
    )

    return {
        "source_index": item.get("source_index"),
        "source_ownership_key": item.get("source_ownership_key"),
        "split_path": item.get("split_path"),
        "x_interval": item.get("x_interval"),
        "y_interval": item.get("y_interval"),
        "branch_validated": branch_validated,
        "f_interval_contains_zero": f_iv_contains_zero,
        "grad_F_well_defined": grad_well_defined,
        "tangent_radius_R": R,
        "f_tt_abs_upper": f_tt_up,
        "third_directional_upper_T3": T3,
        "normal_remainder_bound_Rn": Rn,
        "wall_threshold_rhs": wall_rhs,
        "wall_lhs_before": wall_lhs_before,
        "C0_raw_absolute_hull_old": C0_old,
        "wall_lhs_after_R": new_lhs_R,
        "wall_lhs_after_2R": new_lhs_2R,
        "required_factor_before_rhs_over_lhs": required_factor_before,
        "required_factor_after_rhs_over_lhs_R": required_factor_after_R,
        "required_factor_after_rhs_over_lhs_2R": required_factor_after_2R,
        "outcome": outcome,
    }


# -----------------------------------------------------------------------------
# Aggregation and verdict.
# -----------------------------------------------------------------------------

def collect_wall_failures(
    boundary_slice: dict[str, Any],
    owner_family_filter: str | None = None,
) -> list[dict[str, Any]]:
    """Pull all 'third_order_wall_separation_failed' rows, optionally
    filtered by source_ownership_key."""
    rows = []
    for item in boundary_slice.get("remaining_unresolved", []):
        if item.get("reason") != "third_order_wall_separation_failed":
            continue
        if owner_family_filter is not None and item.get(
            "source_ownership_key"
        ) != owner_family_filter:
            continue
        rows.append(item)
    return rows


def aggregate_verdict(classified: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(classified)
    if n == 0:
        return {
            "n_boxes": 0,
            "verdict": "WS01_INPUT_EMPTY",
            "branch_point_validation_rate": None,
            "closure_rate_at_2R": None,
            "closure_rate_at_R_tight": None,
            "structural_difference_rate": None,
            "median_required_factor_before": None,
            "median_required_factor_after_2R": None,
            "n_closed_2R": 0,
            "n_closed_tight_R_only": 0,
            "n_still_failing": 0,
            "n_branch_not_found": 0,
            "n_branch_not_well_defined": 0,
        }

    n_closed_2R = sum(1 for c in classified if c["outcome"] == "CLOSED_BY_WS01")
    n_closed_tight = sum(
        1 for c in classified if c["outcome"] == "CLOSED_BY_WS01_TIGHT_R_ONLY"
    )
    n_still = sum(
        1 for c in classified if c["outcome"] == "STILL_FAILING_AFTER_WS01"
    )
    n_branch_not_found = sum(
        1 for c in classified if c["outcome"] == "BRANCH_POINT_NOT_FOUND"
    )
    n_branch_ill = sum(
        1
        for c in classified
        if c["outcome"] == "BRANCH_POINT_NOT_WELL_DEFINED"
    )

    branch_validated = sum(1 for c in classified if c["branch_validated"])
    branch_rate = branch_validated / n
    closure_2R = n_closed_2R / n
    closure_tight = (n_closed_2R + n_closed_tight) / n
    structural_diff = (n_branch_not_found + n_branch_ill) / n

    if branch_rate < 0.95:
        verdict = "WS01_INAPPLICABLE_STRUCTURAL"
    elif closure_tight >= 0.50:
        verdict = "WS01_APPLICABLE"
    elif closure_tight >= 0.10:
        verdict = "WS01_PARTIAL"
    else:
        verdict = "WS01_INAPPLICABLE_QUANTITATIVE"

    before_factors = [
        c["required_factor_before_rhs_over_lhs"]
        for c in classified
        if c["required_factor_before_rhs_over_lhs"] is not None
    ]
    after_factors_2R = [
        c["required_factor_after_rhs_over_lhs_2R"]
        for c in classified
        if c["required_factor_after_rhs_over_lhs_2R"] is not None
    ]
    median_before = (
        statistics.median(before_factors) if before_factors else None
    )
    median_after_2R = (
        statistics.median(after_factors_2R) if after_factors_2R else None
    )

    return {
        "n_boxes": n,
        "verdict": verdict,
        "branch_point_validation_rate": branch_rate,
        "closure_rate_at_2R": closure_2R,
        "closure_rate_at_R_tight": closure_tight,
        "structural_difference_rate": structural_diff,
        "median_required_factor_before": median_before,
        "median_required_factor_after_2R": median_after_2R,
        "n_closed_2R": n_closed_2R,
        "n_closed_tight_R_only": n_closed_tight,
        "n_still_failing": n_still,
        "n_branch_not_found": n_branch_not_found,
        "n_branch_not_well_defined": n_branch_ill,
    }


# -----------------------------------------------------------------------------
# Artifact assembly.
# -----------------------------------------------------------------------------

def build_experiment_id(problem_id: int, run_date: str) -> str:
    return f"EXP-MATH-ERDOS-{problem_id}-WS01-APPLICABILITY-PROBE-{run_date}-01"


def build_results(
    problem_id: int,
    experiment_id: str,
    input_json_path: Path,
    input_json_sha: str,
    owner_family_filter: str | None,
    classified: list[dict[str, Any]],
    schema_errors: list[str],
) -> dict[str, Any]:
    agg = aggregate_verdict(classified)
    return {
        "experiment_id": experiment_id,
        "schema_version": "1.0",
        "timestamp_unix": int(time.time()),
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "generated_by": "cross_problem_ws01_applicability_probe",
        "origin": "auto-research",
        "promotion_state": "review_only",
        "target_problem": problem_id,
        "owner_family_filter": owner_family_filter,
        "claim_ceiling": (
            "Internal experiment artifact only. This is a DIAGNOSTIC, NOT "
            "a proof. The probe reports whether the WS-01 architecture "
            "(Branch-Centered Moving-Frame Collar wall-separation rewrite) "
            "is applicable to the supplied boundary-slice failure-box data. "
            "It does not certify the wall separation, the boundary slice, "
            "any cell decomposition, or the underlying Erdos problem."
        ),
        "wall_test_rewrite": {
            "name": "WS-01-CENTER-STRIP-CANCELLATION",
            "old_LHS": "sup_r |F(0,r)| + R_n",
            "new_LHS": "0.5 * |F_tt|_upper * R^2 + (T3/6) * R^3 + R_n",
            "RHS": "S * lower(|F_n|)",
            "math_reference": (
                "Erdos114/EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_"
                "2026-05-06.md lines 122-146"
            ),
            "z_star_validation": (
                "interval-arithmetic only: f_interval contains 0 (IVT) AND "
                "center_fn_abs_lower > 0 (grad F nonvanishing)"
            ),
            "F_t_at_z_star": (
                "zero by construction: t is grad F rotated 90 deg, so F "
                "restricted to zero curve is constant 0"
            ),
            "R_choice": (
                "R = tangent_radius_upper from source data (tight bound); "
                "2R conservative bound for off-center z*"
            ),
        },
        "input_source": {
            "path": str(input_json_path),
            "sha256": input_json_sha,
        },
        "schema_validation": {
            "schema_errors": schema_errors,
            "validation_status": "PASS" if not schema_errors else "FAIL",
        },
        "aggregate": agg,
        "per_failure_outcomes": classified,
        "interpretation": _interpret(agg),
        "next_action": _next_action(agg),
        "forbidden_writes": [
            "D1",
            "morphisms.json",
            "proof registries",
            "public pages",
            "existing EHP114 finite packets",
            "existing corrected n=13 or n=14 receipts",
            "reference implementation script",
        ],
    }


def _interpret(agg: dict[str, Any]) -> str:
    if agg["n_boxes"] == 0:
        return (
            "Probe ran on 0 failure boxes. Either the supplied JSON had no "
            "'third_order_wall_separation_failed' rows or the owner-family "
            "filter excluded all rows. No verdict possible."
        )
    verdict = agg["verdict"]
    n = agg["n_boxes"]
    if verdict == "WS01_APPLICABLE":
        return (
            f"On {n} wall-separation failure boxes, WS-01 closes "
            f"{agg['n_closed_2R']} at the conservative 2R bound and an "
            f"additional {agg['n_closed_tight_R_only']} at the tight-R "
            f"bound. Branch-point validation rate "
            f"{agg['branch_point_validation_rate']:.2%}, closure rate at "
            f"tight R {agg['closure_rate_at_R_tight']:.2%}. The architecture "
            f"appears applicable to this problem's boundary-slice failures."
        )
    if verdict == "WS01_PARTIAL":
        return (
            f"On {n} wall-separation failure boxes, WS-01 closes "
            f"{agg['n_closed_2R']} at 2R and {agg['n_closed_tight_R_only']} "
            f"more at tight R. Branch-point validation rate "
            f"{agg['branch_point_validation_rate']:.2%}, closure rate at "
            f"tight R {agg['closure_rate_at_R_tight']:.2%}. The architecture "
            f"applies but does not clear a majority of boxes."
        )
    if verdict == "WS01_INAPPLICABLE_STRUCTURAL":
        return (
            f"On {n} wall-separation failure boxes, only "
            f"{agg['branch_point_validation_rate']:.2%} have a validated "
            f"interval branch point (IVT + grad-F nonvanishing). The WS-01 "
            f"architecture itself does not apply -- most boxes are not "
            f"sitting on the F=0 curve in a way that admits the moving-frame "
            f"collar reduction."
        )
    if verdict == "WS01_INAPPLICABLE_QUANTITATIVE":
        return (
            f"On {n} wall-separation failure boxes, "
            f"{agg['branch_point_validation_rate']:.2%} have a validated "
            f"interval branch point, but only "
            f"{agg['closure_rate_at_R_tight']:.2%} of them clear the wall "
            f"threshold even at the tight-R bound. The architecture applies "
            f"but the local constants (|F_tt|, T3, R) do not favor closure."
        )
    return f"Unknown verdict: {verdict}"


def _next_action(agg: dict[str, Any]) -> str:
    if agg["n_boxes"] == 0:
        return (
            "No boxes processed. Re-supply a boundary-slice JSON with at "
            "least one third_order_wall_separation_failed row, or relax "
            "the owner-family filter."
        )
    verdict = agg["verdict"]
    if verdict == "WS01_APPLICABLE":
        return (
            "Proceed with WS-01 for this problem. Run owner-family-by-"
            "owner-family probes to confirm uniform closure across the "
            "full wall-failure set, then attempt Rust-level interval "
            "certification of the moving-frame collar."
        )
    if verdict == "WS01_PARTIAL":
        return (
            "WS-01 partially applicable. Identify which owner families "
            "close vs. which do not; consider tighter R partitions or a "
            "hybrid decomposition that uses WS-01 only on the favorable "
            "subset."
        )
    if verdict == "WS01_INAPPLICABLE_STRUCTURAL":
        return (
            "Do NOT pursue WS-01 for this problem. The boundary-slice "
            "geometry does not present validated branch points; the "
            "underlying analytic reduction differs from #114's CELL-02-03 "
            "case. Look for a different architecture."
        )
    if verdict == "WS01_INAPPLICABLE_QUANTITATIVE":
        return (
            "WS-01 architecture applies but the constants do not favor "
            "closure. Try smaller boxes (reduce R) or a sharper |F_tt| / "
            "T3 estimate before re-running the probe."
        )
    return "No recommendation; unknown verdict."


def build_report(results: dict[str, Any]) -> str:
    agg = results["aggregate"]
    n = agg["n_boxes"]
    rows = results["per_failure_outcomes"]
    table_lines = [
        "| # | src_idx | owner | R | wall_LHS_before | wall_LHS_after_2R | "
        "wall_rhs | outcome |",
        "|---|---------|-------|---|-----------------|-------------------|"
        "----------|---------|",
    ]
    for i, c in enumerate(rows):
        owner = c.get("source_ownership_key") or ""
        table_lines.append(
            f"| {i} | {c['source_index']} | `{owner}` | "
            f"{c['tangent_radius_R']:.4g} | {c['wall_lhs_before']:.4g} | "
            f"{c['wall_lhs_after_2R']:.4g} | "
            f"{c['wall_threshold_rhs']:.4g} | `{c['outcome']}` |"
        )
    table = "\n".join(table_lines)

    schema = results["schema_validation"]
    schema_block = (
        "All rows passed schema validation."
        if schema["validation_status"] == "PASS"
        else "Schema errors detected:\n" + "\n".join(
            f"- {e}" for e in schema["schema_errors"]
        )
    )

    return f"""# {results["experiment_id"]}

## Scope

Diagnostic, not a proof. Cross-problem applicability probe for the WS-01
(Branch-Centered Moving-Frame Collar) wall-separation rewrite, generalized
from Erdos #114's CELL-02-03 owner-family-4488:3 reference run.

## Verdict

`{agg["verdict"]}`

## Headline numbers

- Failure boxes processed: `{n}`
- Branch-point validation rate: `{agg["branch_point_validation_rate"]}`
- Closure rate at conservative 2R: `{agg["closure_rate_at_2R"]}`
- Closure rate at tight R (cumulative): `{agg["closure_rate_at_R_tight"]}`
- Structural-difference rate: `{agg["structural_difference_rate"]}`
- Closed at 2R: `{agg["n_closed_2R"]}`
- Closed at tight R only: `{agg["n_closed_tight_R_only"]}`
- Still failing: `{agg["n_still_failing"]}`
- Branch point not found: `{agg["n_branch_not_found"]}`
- Branch point not well-defined: `{agg["n_branch_not_well_defined"]}`
- Median required factor (RHS/LHS) before WS-01: \
`{agg["median_required_factor_before"]}`
- Median required factor (RHS/LHS) after WS-01 at 2R: \
`{agg["median_required_factor_after_2R"]}`

## Schema validation

{schema_block}

## What WS-01 does

The original wall test bounds `sup_r |F(0,r)|` by the raw absolute interval
hull of `|F|` on the center strip. That hull carries the linear-in-R
coefficient `|F_t(box-midpoint)|`, which can be O(1).

WS-01 picks a validated zero point `z*` on the F=0 curve in each box and
uses the tangent direction `t` to that curve as the new coordinate. By
construction `F(z*) = 0` and `F_t(z*) = 0`, so the linear-in-R term
vanishes:

```
new_LHS = 0.5 |F_tt| R^2 + (T3/6) R^3 + R_n
```

Branch-point validation per box uses two interval-arithmetic facts:

1. `f_interval` straddles zero (IVT => F has a zero in the box).
2. `center_fn_abs_lower > 0` (grad F nonvanishing => zero curve is a smooth
   1-manifold => `t` is well-defined).

`F_t(z*) = 0` then holds by construction, not by any numerical estimate.

`R` is taken as `tangent_radius_upper` from the source data. The probe
reports the conservative `2R` figure as the canonical pass condition
(worst case z* near a box corner) and the tight `R` figure as a sensitivity
diagnostic.

## Per-failure outcomes

{table}

## Interpretation

{results["interpretation"]}

## Next action

{results["next_action"]}

## Source provenance

- Input JSON: `{results["input_source"]["path"]}`
- SHA-256: `{results["input_source"]["sha256"]}`

## Safety

This artifact is review-only. Forbidden writes (D1, curated morphisms,
public pages, proof registries, existing EHP114 finite packets,
reference implementation script) were not attempted. The probe writes
only this experiment folder.
"""


# -----------------------------------------------------------------------------
# Probe entry point (problem-agnostic).
# -----------------------------------------------------------------------------

def run_probe(
    problem_id: int,
    input_json_path: Path,
    output_dir: Path,
    owner_family_filter: str | None,
    run_date: str | None = None,
) -> dict[str, Any]:
    if not input_json_path.exists():
        raise SystemExit(f"Missing input JSON: {input_json_path}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(
            f"Refusing to overwrite non-empty output directory: {output_dir}"
        )

    if run_date is None:
        run_date = datetime.now(timezone.utc).strftime("%Y%m%d")
    experiment_id = build_experiment_id(problem_id, run_date)

    boundary_slice = load_json(input_json_path)
    input_json_sha = sha256_file(input_json_path)

    failures = collect_wall_failures(boundary_slice, owner_family_filter)

    schema_errors: list[str] = []
    for i, item in enumerate(failures):
        schema_errors.extend(validate_box_schema(item, i))

    if schema_errors:
        # Still produce an artifact, but mark the verdict accordingly.
        classified = []
    else:
        classified = [classify_box(item) for item in failures]

    results = build_results(
        problem_id=problem_id,
        experiment_id=experiment_id,
        input_json_path=input_json_path,
        input_json_sha=input_json_sha,
        owner_family_filter=owner_family_filter,
        classified=classified,
        schema_errors=schema_errors,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    results_json = output_dir / f"{experiment_id}_RESULTS.json"
    report_md = output_dir / f"{experiment_id}_REPORT.md"
    results_sha = output_dir / f"{experiment_id}_RESULTS.sha256"

    write_json(results_json, results)
    report_md.write_text(build_report(results), encoding="utf-8")
    results_sha.write_text(
        f"{sha256_file(results_json)}  {results_json.name}\n",
        encoding="utf-8",
    )

    return {
        "experiment_id": experiment_id,
        "results_json": str(results_json),
        "report_md": str(report_md),
        "results_sha": str(results_sha),
        "verdict": results["aggregate"]["verdict"],
        "aggregate": results["aggregate"],
    }


# -----------------------------------------------------------------------------
# Regression test against the #114 CELL-02-03 4488:3 reference run.
# -----------------------------------------------------------------------------

def run_regression() -> int:
    """Verify the probe matches the 2026-05-09 reference run on the 16 boxes
    of CELL-02-03 owner family 4488:3. Returns 0 on PASS, 1 on FAIL."""
    print("[regression] loading reference boundary-slice...")
    if not REGRESSION_BOUNDARY_SLICE_JSON.exists():
        print(
            f"[regression] FAIL: reference data missing: "
            f"{REGRESSION_BOUNDARY_SLICE_JSON}",
            file=sys.stderr,
        )
        return 1

    boundary_slice = load_json(REGRESSION_BOUNDARY_SLICE_JSON)
    failures = collect_wall_failures(
        boundary_slice, owner_family_filter=REGRESSION_OWNER_FAMILY
    )
    print(f"[regression] collected {len(failures)} failure boxes")
    if len(failures) != REGRESSION_EXPECTED_TOTAL:
        print(
            f"[regression] FAIL: expected {REGRESSION_EXPECTED_TOTAL} boxes, "
            f"got {len(failures)}",
            file=sys.stderr,
        )
        return 1

    schema_errors: list[str] = []
    for i, item in enumerate(failures):
        schema_errors.extend(validate_box_schema(item, i))
    if schema_errors:
        print(
            f"[regression] FAIL: schema errors: {schema_errors}",
            file=sys.stderr,
        )
        return 1

    classified = [classify_box(item) for item in failures]
    agg = aggregate_verdict(classified)

    print(f"[regression] verdict: {agg['verdict']}")
    print(f"[regression] n_closed_2R: {agg['n_closed_2R']}")
    print(f"[regression] n_closed_tight_R_only: {agg['n_closed_tight_R_only']}")
    print(f"[regression] n_still_failing: {agg['n_still_failing']}")
    print(f"[regression] n_branch_not_found: {agg['n_branch_not_found']}")
    print(
        f"[regression] n_branch_not_well_defined: "
        f"{agg['n_branch_not_well_defined']}"
    )

    expected = {
        "n_boxes": REGRESSION_EXPECTED_TOTAL,
        "n_closed_2R": REGRESSION_EXPECTED_CLOSED_2R,
        "n_closed_tight_R_only": REGRESSION_EXPECTED_CLOSED_TIGHT_ONLY,
        "n_branch_not_found": REGRESSION_EXPECTED_BRANCH_NOT_FOUND,
    }
    failures_list: list[str] = []
    for k, v in expected.items():
        if agg[k] != v:
            failures_list.append(
                f"  {k}: expected {v}, got {agg[k]}"
            )

    if failures_list:
        print(
            "[regression] FAIL:\n" + "\n".join(failures_list),
            file=sys.stderr,
        )
        return 1

    if agg["verdict"] != "WS01_APPLICABLE":
        print(
            f"[regression] FAIL: verdict expected WS01_APPLICABLE, "
            f"got {agg['verdict']}",
            file=sys.stderr,
        )
        return 1

    print(
        "[regression] PASS: 16/16 boxes close at conservative 2R; "
        "verdict WS01_APPLICABLE matches 2026-05-09 reference run."
    )
    return 0


# -----------------------------------------------------------------------------
# CLI.
# -----------------------------------------------------------------------------

def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="WS-01 applicability probe (problem-agnostic).",
    )
    parser.add_argument(
        "--problem-id",
        type=int,
        help="Erdos problem ID (e.g., 1038, 20, 233).",
    )
    parser.add_argument(
        "--input-json",
        type=Path,
        help="Boundary-slice failure-box JSON (same schema as #114).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory (created if missing; must not be non-empty).",
    )
    parser.add_argument(
        "--owner-family-filter",
        type=str,
        default=None,
        help="Optional source_ownership_key filter (e.g., '4488:3').",
    )
    parser.add_argument(
        "--regression",
        action="store_true",
        help=(
            "Run the #114 CELL-02-03 owner-family-4488:3 regression test "
            "and exit. No other args needed."
        ),
    )
    args = parser.parse_args(argv)

    if args.regression:
        return run_regression()

    missing = [
        flag
        for flag, val in (
            ("--problem-id", args.problem_id),
            ("--input-json", args.input_json),
            ("--output-dir", args.output_dir),
        )
        if val is None
    ]
    if missing:
        parser.error(f"missing required arguments: {', '.join(missing)}")

    summary = run_probe(
        problem_id=args.problem_id,
        input_json_path=args.input_json,
        output_dir=args.output_dir,
        owner_family_filter=args.owner_family_filter,
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
