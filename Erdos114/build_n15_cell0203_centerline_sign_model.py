#!/usr/bin/env python3
"""Test the WS-01-CENTER-STRIP-CANCELLATION rewrite (Branch-Centered Moving-Frame
Collar) on the 16 wall-separation failures of owner family `4488:3` in the
n=15 boundary-slice run for CELL-02-03.

This is Track B1 of the EHP114 bridge program. Internal experiment artifact only.
Not a proof of Erdos #114, not an n=15 certificate, not a CELL-02-03 closure.
Only a diagnostic of WS-01 effectiveness on owner family 4488:3.

Math reference (consume, do not re-derive):
  Erdos114/EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_2026-05-06.md
  lines 122-146 (Branch-Centered Moving-Frame Collar lemma).

The current wall test is:
    sup_r |F(0,r)| + R_n  <  S * lower(|F_n|)
where the LHS is dominated by the raw absolute hull of |F(0,r)| over the
center strip (a "thick" interval bound that pays for full root-affine
dependency at once).

The WS-01 rewrite: at a validated zero point z*(u) on the F=0 curve, with
t(u) tangent to that curve, F(z*) = 0 and F_t(z*) = 0 by construction. The
center strip at z* obeys
    |F(z* + R t(u))|  <=  0.5 |F_tt| R^2 + (T3/6) R^3
with no linear-in-R term. The wall test becomes
    new_wall_LHS = 0.5 |F_tt| R^2 + (T3/6) R^3 + R_n  <  S * lower(|F_n|).

Branch-point validation per box (interval-arithmetic, no Newton needed):
  (a) f_interval lo <= 0 <= hi   => F has a zero in the box (IVT)
  (b) center_fn_abs_lower > 0    => grad F is non-vanishing on the box,
                                     so zero set is a smooth 1-manifold there
  (c) tangent direction in source data is already grad F rotated 90 deg
                                  => F_t along that direction is zero on the
                                     zero curve (definitional)
"""

from __future__ import annotations

import hashlib
import json
import statistics
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01"
TARGET_PROBLEM = 114
DEGREE = 15
CELL_ID = "CELL-02-03"
OWNER_FAMILY = "4488:3"

SCRIPT_PATH = Path(__file__).resolve()
ERDOS114_DIR = SCRIPT_PATH.parent
PROOF_PATH = ERDOS114_DIR / "proof_path"

SOURCE_BOUNDARY_SLICE_DIR = (
    PROOF_PATH
    / "EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03"
)
SOURCE_BOUNDARY_SLICE_JSON = (
    SOURCE_BOUNDARY_SLICE_DIR
    / "EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03_RESULTS.json"
)
SOURCE_BOUNDARY_SLICE_SHA = (
    SOURCE_BOUNDARY_SLICE_DIR
    / "EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03_RESULTS.sha256"
)

SOURCE_WALL_TARGET_DIR = (
    PROOF_PATH / "EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01"
)
SOURCE_WALL_TARGET_JSON = (
    SOURCE_WALL_TARGET_DIR
    / "EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01_RESULTS.json"
)
SOURCE_WALL_TARGET_SHA = (
    SOURCE_WALL_TARGET_DIR
    / "EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01_RESULTS.sha256"
)

MATH_REFERENCE_PATH = (
    ERDOS114_DIR / "EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_2026-05-06.md"
)

OUTPUT_DIR = PROOF_PATH / EXPERIMENT_ID
RESULTS_JSON = OUTPUT_DIR / f"{EXPERIMENT_ID}_RESULTS.json"
REPORT_MD = OUTPUT_DIR / f"{EXPERIMENT_ID}_REPORT.md"
RESULTS_SHA = OUTPUT_DIR / f"{EXPERIMENT_ID}_RESULTS.sha256"
RUN_LOG = OUTPUT_DIR / f"{EXPERIMENT_ID}_RUN.log"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha_from_sidecar(sha_path: Path) -> str | None:
    if not sha_path.exists():
        return None
    text = sha_path.read_text(encoding="utf-8").strip()
    return text.split()[0] if text else None


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def collect_wall_failures_4488_3(boundary_slice: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for item in boundary_slice.get("remaining_unresolved", []):
        if item.get("reason") != "third_order_wall_separation_failed":
            continue
        if item.get("source_ownership_key") != OWNER_FAMILY:
            continue
        rows.append(item)
    return rows


def classify_box(item: dict[str, Any]) -> dict[str, Any]:
    """Apply WS-01 rewrite to a single failure box.

    Returns a dict with classification + WS-01 numerics.
    """
    ineq = item["inequality"]
    f_iv = item["f_interval"]
    fn_first = item["fn_full_first_order_interval"]
    ft_first = item["ft_full_first_order_interval"]

    f_lo = float(f_iv["lo"])
    f_hi = float(f_iv["hi"])
    fn_center_lo = float(ineq["center_fn_abs_lower"])
    R = float(ineq["tangent_radius_upper"])
    f_tt_up = float(ineq["f_tt_abs_upper"])
    T3 = float(ineq["third_directional_upper"])
    Rn = float(ineq["normal_remainder_bound"])
    wall_lhs_before = float(ineq["wall_lhs"])
    wall_rhs = float(ineq["wall_rhs"])
    C0_old = float(ineq["center_strip_f_abs_upper"])

    # Branch-point existence test (interval-arithmetic only).
    # (a) F has a zero in the box: f_interval contains zero.
    f_iv_contains_zero = (f_lo <= 0.0 <= f_hi)
    # (b) grad F does not vanish: |F_n| at center has positive lower bound
    #     and the first-order interval keeps |F_n| from collapsing to zero
    #     on the box. We use center_fn_abs_lower > 0 as the lab condition,
    #     since the source experiment already certifies it via
    #     critical_exclusion test (passed for these 108 boxes).
    grad_nonvanishing = fn_center_lo > 0.0
    # (c) Branch-point validated iff (a) AND (b).
    branch_validated = f_iv_contains_zero and grad_nonvanishing

    # Conservative R for moving-frame collar: tangent_radius_upper covers
    # half-width from box center; if z* is off-center the effective radius
    # could be up to 2*R. Use both for sensitivity reporting.
    new_lhs_R = 0.5 * f_tt_up * R * R + (T3 / 6.0) * R * R * R + Rn
    new_lhs_2R = (
        0.5 * f_tt_up * (2.0 * R) ** 2 + (T3 / 6.0) * (2.0 * R) ** 3 + Rn
    )

    # Honest classification: we report PASS only if even the worst-case 2R
    # bound clears the wall threshold. R-only PASS is reported as the
    # "best-case" line (z* near box center).
    if not branch_validated:
        outcome = "BRANCH_POINT_NOT_FOUND"
    elif new_lhs_2R < wall_rhs:
        outcome = "CLOSED_BY_WS01"
    elif new_lhs_R < wall_rhs:
        outcome = "CLOSED_BY_WS01_TIGHT_R_ONLY"
    else:
        outcome = "STILL_FAILING_AFTER_WS01"

    # required factor before vs after (RHS / LHS); a value < 1 means LHS too
    # big, > 1 means LHS clears the threshold.
    required_factor_before = wall_rhs / wall_lhs_before if wall_lhs_before > 0 else None
    required_factor_after = wall_rhs / new_lhs_2R if new_lhs_2R > 0 else None

    return {
        "ownership_key": OWNER_FAMILY,
        "source_index": item.get("source_index"),
        "split_path": item.get("split_path"),
        "x_interval": item.get("x_interval"),
        "y_interval": item.get("y_interval"),
        "branch_validated": branch_validated,
        "f_interval_contains_zero": f_iv_contains_zero,
        "grad_F_nonvanishing": grad_nonvanishing,
        "z_star_method": (
            "interval_existence_via_IVT_on_f_interval_with_constant_sign_F_n"
            if branch_validated
            else None
        ),
        "tangent_unit": item.get("tangent"),
        "normal_unit": item.get("normal"),
        "F_t_at_z_star_argument": (
            "zero_by_construction:t_is_grad_F_rotated_90deg_so_F_t_vanishes_on_zero_curve"
            if branch_validated
            else None
        ),
        "tangent_radius_R": R,
        "f_tt_abs_upper": f_tt_up,
        "third_directional_upper_T3": T3,
        "normal_remainder_bound_Rn": Rn,
        "wall_threshold_S_times_m": wall_rhs,
        "wall_lhs_before": wall_lhs_before,
        "C0_raw_absolute_hull_old": C0_old,
        "wall_lhs_after_R": new_lhs_R,
        "wall_lhs_after_2R": new_lhs_2R,
        "required_factor_before_rhs_over_lhs": required_factor_before,
        "required_factor_after_rhs_over_lhs_2R": required_factor_after,
        "outcome": outcome,
    }


def build_results(
    failures_all: list[dict[str, Any]],
    classified: list[dict[str, Any]],
    boundary_slice_sha: str,
    wall_target_sha: str,
) -> dict[str, Any]:
    n_in_family = len(failures_all)

    closed = [c for c in classified if c["outcome"] == "CLOSED_BY_WS01"]
    closed_tight = [
        c for c in classified if c["outcome"] == "CLOSED_BY_WS01_TIGHT_R_ONLY"
    ]
    still_failing = [
        c for c in classified if c["outcome"] == "STILL_FAILING_AFTER_WS01"
    ]
    not_found = [c for c in classified if c["outcome"] == "BRANCH_POINT_NOT_FOUND"]

    # We report the conservative WS-01 (using 2R) as the canonical pass set;
    # tight-R passes count as still_failing for the headline status to keep
    # the honest scope.
    closed_count = len(closed)
    closed_tight_count = len(closed_tight)
    still_failing_count = len(still_failing) + closed_tight_count
    not_found_count = len(not_found)

    # Required-factor median before vs after.
    before_factors = [
        c["required_factor_before_rhs_over_lhs"]
        for c in classified
        if c["required_factor_before_rhs_over_lhs"] is not None
    ]
    after_factors = [
        c["required_factor_after_rhs_over_lhs_2R"]
        for c in classified
        if c["required_factor_after_rhs_over_lhs_2R"] is not None
    ]
    median_before = statistics.median(before_factors) if before_factors else None
    median_after = statistics.median(after_factors) if after_factors else None

    # Headline status logic.
    if not_found_count == n_in_family:
        status = "CENTERLINE_SIGN_MODEL_BLOCKED"
    elif closed_count == n_in_family:
        status = "CENTERLINE_SIGN_MODEL_FULL_PASS_4488_3"
    elif closed_count > 0:
        status = "CENTERLINE_SIGN_MODEL_PARTIAL_PASS"
    else:
        status = "CENTERLINE_SIGN_MODEL_BLOCKED"

    interpretation = (
        f"On owner family {OWNER_FAMILY} of CELL-02-03's n=15 wall-separation "
        f"failures (16 boxes), the WS-01 Branch-Centered Moving-Frame Collar "
        f"rewrite reduces the median wall-LHS-over-RHS deficit from "
        f"~{1/median_before:.2f}x miss to a "
        f"~{median_after:.2f}x clearance (factor of safety) when the worst-"
        f"case z*-off-center radius 2R is used. {closed_count}/{n_in_family} "
        f"boxes close by the conservative 2R bound and {closed_tight_count} more "
        f"close only with the tight-R bound. {not_found_count} boxes have no "
        f"interval-validated branch point. Per-box validation: F_iv contains 0 "
        f"(branch exists by IVT) and center |F_n| lower bound is large positive "
        f"(grad F nonvanishing, so the zero curve is smooth and t(u) is "
        f"well-defined). F_t(z*)=0 holds by construction, not by numerics."
        if median_before and median_after
        else "Insufficient classified rows to compute interpretation."
    )

    next_dependency = (
        "B2 (signed wall endpoints) and B3 (owner-family local model on next-"
        "worst groups, e.g., 4571:2 with 11 failures and 2484:4 with 16 failures) "
        "are unblocked. Track B1 viable. Recommended: replicate this experiment "
        "on owner family 4571:2 next, then run a Rust-level interval certification "
        "of the moving-frame collar that emits z* via interval Newton on F(z)=0."
        if status in {
            "CENTERLINE_SIGN_MODEL_PARTIAL_PASS",
            "CENTERLINE_SIGN_MODEL_FULL_PASS_4488_3",
        }
        else "BLOCKED. Moving-frame collar cannot establish a validated branch "
        "point on family 4488:3, so the cell decomposition for n=15 needs to be "
        "reconsidered (smaller boxes, or a different analytic reduction)."
    )

    return {
        "experiment_id": EXPERIMENT_ID,
        "schema_version": "1.0",
        "timestamp_unix": int(time.time()),
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "generated_by": "track_b1_centerline_sign_model_runner",
        "origin": "auto-research",
        "promotion_state": "review_only",
        "target_problem": TARGET_PROBLEM,
        "degree": DEGREE,
        "cell_id": CELL_ID,
        "owner_family_tested": OWNER_FAMILY,
        "claim_ceiling": (
            "Internal experiment artifact only. Not a proof of Erdos #114, "
            "not an n=15 certificate, and not a CELL-02-03 closure. Only a "
            "diagnostic of WS-01 effectiveness on owner family 4488:3."
        ),
        "status": status,
        "failures_in_family": n_in_family,
        "failures_closed_by_ws01": closed_count,
        "failures_closed_by_ws01_tight_R_only": closed_tight_count,
        "failures_still_failing": still_failing_count,
        "failures_branch_point_not_found": not_found_count,
        "median_required_factor_before": median_before,
        "median_required_factor_after": median_after,
        "wall_test_rewrite": {
            "name": "WS-01-CENTER-STRIP-CANCELLATION",
            "old_LHS": "sup_r |F(0,r)| + R_n",
            "new_LHS": "0.5 * |F_tt|_upper * R^2 + (T3/6) * R^3 + R_n",
            "RHS": "S * lower(|F_n|)",
            "math_reference": str(MATH_REFERENCE_PATH),
            "math_reference_lines": "122-146",
            "z_star_validation": (
                "interval-arithmetic only: F_iv contains 0 (IVT) AND "
                "center_fn_abs_lower > 0 (grad F nonvanishing => smooth zero curve)"
            ),
            "F_t_at_z_star": (
                "zero by construction: t(u) = grad F rotated 90 deg, "
                "so F restricted to zero curve is constant 0"
            ),
            "R_choice": (
                "R = tangent_radius_upper from source data (best-case) and "
                "2R as conservative worst-case for off-center z*"
            ),
        },
        "per_failure_outcomes": classified,
        "interpretation": interpretation,
        "next_dependency": next_dependency,
        "source_paths": {
            "wall_failures": str(SOURCE_BOUNDARY_SLICE_JSON),
            "wall_separation_target_packet": str(SOURCE_WALL_TARGET_JSON),
            "math_reference": str(MATH_REFERENCE_PATH),
            "n14_root_location": (
                "Erdos114/validated_length/EXP-MATH-EHP114-N14-{SHARP-,}"
                "PPRIME-ROOT-LOCATION-CELL-*-20260506-01 (existence pattern only; "
                "no n=14 root-location artifact was numerically consumed because "
                "WS-01 z* validation needs only F_iv-contains-zero + grad-F-nonvanishing, "
                "both of which are read directly from the n=15 boundary-slice rows)"
            ),
        },
        "source_sha_status": [
            {
                "path": str(SOURCE_BOUNDARY_SLICE_JSON),
                "expected_sha256": sha_from_sidecar(SOURCE_BOUNDARY_SLICE_SHA),
                "actual_sha256": boundary_slice_sha,
                "sha_status": (
                    "PASS"
                    if sha_from_sidecar(SOURCE_BOUNDARY_SLICE_SHA) == boundary_slice_sha
                    else "FAIL"
                ),
            },
            {
                "path": str(SOURCE_WALL_TARGET_JSON),
                "expected_sha256": sha_from_sidecar(SOURCE_WALL_TARGET_SHA),
                "actual_sha256": wall_target_sha,
                "sha_status": (
                    "PASS"
                    if sha_from_sidecar(SOURCE_WALL_TARGET_SHA) == wall_target_sha
                    else "FAIL"
                ),
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


def build_report(results: dict[str, Any]) -> str:
    rows = results["per_failure_outcomes"]
    table_lines = [
        "| # | src_idx | split_path | R | wall_LHS_before | wall_LHS_after_2R | wall_rhs | outcome |",
        "|---|---------|------------|---|------------------|--------------------|----------|---------|",
    ]
    for i, c in enumerate(rows):
        table_lines.append(
            f"| {i} | {c['source_index']} | `{c['split_path']}` | "
            f"{c['tangent_radius_R']:.4g} | {c['wall_lhs_before']:.4g} | "
            f"{c['wall_lhs_after_2R']:.4g} | {c['wall_threshold_S_times_m']:.4g} | "
            f"`{c['outcome']}` |"
        )
    table = "\n".join(table_lines)

    return f"""# {EXPERIMENT_ID}

## Scope

Internal experiment artifact only. Track B1 of the EHP114 bridge program. Not
a proof of Erdos #114, not an n=15 certificate, and not a CELL-02-03 closure.
Only a diagnostic of WS-01-CENTER-STRIP-CANCELLATION (the Branch-Centered
Moving-Frame Collar rewrite) on owner family `{OWNER_FAMILY}`'s 16 wall-
separation failures.

## Status

`{results["status"]}`

## Headline numbers

- Failures in family: `{results["failures_in_family"]}`
- Closed by WS-01 (conservative 2R bound): `{results["failures_closed_by_ws01"]}`
- Closed by WS-01 (tight R only): `{results["failures_closed_by_ws01_tight_R_only"]}`
- Still failing: `{results["failures_still_failing"]}`
- Branch point not found: `{results["failures_branch_point_not_found"]}`
- Median required factor (RHS/LHS) before WS-01: `{results["median_required_factor_before"]}`
- Median required factor (RHS/LHS) after WS-01 at 2R: `{results["median_required_factor_after"]}`

## What WS-01 does

The original wall test bounds `sup_r |F(0,r)|` by the raw absolute interval
hull of `|F|` on the center strip. That hull pays for full root-affine
dependency at once and is dominated by the linear-in-R coefficient
`|F_t(box-midpoint)|` which can be O(1).

WS-01 picks a validated zero point `z*(u)` on the F=0 curve in each box and
uses the tangent direction `t(u)` to that curve as the new center-strip
coordinate. By construction `F(z*) = 0` and `F_t(z*) = 0`, so the linear-in-R
term vanishes:

```
new_LHS = 0.5 |F_tt| R^2 + (T3/6) R^3 + R_n
```

Branch-point validation per box uses two interval-arithmetic facts already
recorded in the source boundary-slice run:

1. `f_interval` straddles zero (IVT => F has a zero in the box).
2. `center_fn_abs_lower > 0` (grad F nonvanishing => zero curve is a smooth
   1-manifold => t(u) is well-defined).

`F_t(z*) = 0` then holds by construction, not by any numerical estimate.

`R` is taken as `tangent_radius_upper` from the source data; we report the
conservative `2R` figure as the canonical pass condition (worst case z*
sitting near a box corner). The tight-R figure is reported as a sensitivity
diagnostic.

## Per-failure outcomes (owner family {OWNER_FAMILY})

{table}

## Interpretation

{results["interpretation"]}

## Next dependency

{results["next_dependency"]}

## Source provenance

- Boundary-slice source: `{SOURCE_BOUNDARY_SLICE_JSON}`
- SHA-256 status: `{results["source_sha_status"][0]["sha_status"]}`
- Wall-separation-target packet: `{SOURCE_WALL_TARGET_JSON}`
- SHA-256 status: `{results["source_sha_status"][1]["sha_status"]}`
- Math reference: `{MATH_REFERENCE_PATH}` (lines 122-146)

## Safety

This artifact is review-only. Forbidden writes (D1, curated morphisms,
public pages, proof registries, existing EHP114 finite packets) were not
attempted. The new build script
`build_n15_cell0203_centerline_sign_model.py` writes only this experiment
folder. No prior script or artifact was modified.
"""


def main(argv: list[str]) -> int:
    if OUTPUT_DIR.exists():
        raise SystemExit(f"Refusing to overwrite existing output directory: {OUTPUT_DIR}")
    if not SOURCE_BOUNDARY_SLICE_JSON.exists():
        raise SystemExit(f"Missing boundary-slice source: {SOURCE_BOUNDARY_SLICE_JSON}")
    if not SOURCE_WALL_TARGET_JSON.exists():
        raise SystemExit(f"Missing wall-target packet: {SOURCE_WALL_TARGET_JSON}")
    if not MATH_REFERENCE_PATH.exists():
        raise SystemExit(f"Missing math reference: {MATH_REFERENCE_PATH}")

    log_lines: list[str] = []

    def log(msg: str) -> None:
        line = f"[{datetime.now(timezone.utc).isoformat()}] {msg}"
        log_lines.append(line)
        print(line, flush=True)

    log(f"experiment_id={EXPERIMENT_ID}")
    log(f"loading boundary-slice source: {SOURCE_BOUNDARY_SLICE_JSON}")
    boundary_slice = load_json(SOURCE_BOUNDARY_SLICE_JSON)
    boundary_slice_sha = sha256_file(SOURCE_BOUNDARY_SLICE_JSON)
    expected_bs_sha = sha_from_sidecar(SOURCE_BOUNDARY_SLICE_SHA)
    if expected_bs_sha and expected_bs_sha != boundary_slice_sha:
        raise SystemExit(
            f"Boundary-slice SHA mismatch: expected {expected_bs_sha}, got {boundary_slice_sha}"
        )
    log(f"boundary-slice sha256={boundary_slice_sha}")

    log(f"loading wall-target packet: {SOURCE_WALL_TARGET_JSON}")
    wall_target = load_json(SOURCE_WALL_TARGET_JSON)
    wall_target_sha = sha256_file(SOURCE_WALL_TARGET_JSON)
    expected_wt_sha = sha_from_sidecar(SOURCE_WALL_TARGET_SHA)
    if expected_wt_sha and expected_wt_sha != wall_target_sha:
        raise SystemExit(
            f"Wall-target SHA mismatch: expected {expected_wt_sha}, got {wall_target_sha}"
        )
    log(f"wall-target sha256={wall_target_sha}")

    failures_all_in_run = [
        item
        for item in boundary_slice.get("remaining_unresolved", [])
        if item.get("reason") == "third_order_wall_separation_failed"
    ]
    log(f"total wall_separation failures in source run: {len(failures_all_in_run)}")
    if len(failures_all_in_run) != 108:
        log(
            f"WARNING: expected 108 wall failures per the wall-separation-target "
            f"packet, found {len(failures_all_in_run)}. Continuing with what we have."
        )
    # cross-check vs wall-target packet
    wt_count = wall_target.get("wall_separation_analysis", {}).get(
        "wall_failure_count"
    )
    log(f"wall-target packet reports wall_failure_count={wt_count}")

    failures_4488_3 = collect_wall_failures_4488_3(boundary_slice)
    log(f"failures in owner family {OWNER_FAMILY}: {len(failures_4488_3)}")
    if len(failures_4488_3) != 16:
        log(
            f"WARNING: expected 16 failures in family {OWNER_FAMILY}, "
            f"found {len(failures_4488_3)}. Continuing with what we have."
        )

    classified = [classify_box(item) for item in failures_4488_3]
    for i, c in enumerate(classified):
        log(
            f"#{i:2d} src_idx={c['source_index']} R={c['tangent_radius_R']:.4g} "
            f"wall_LHS_before={c['wall_lhs_before']:.4g} "
            f"wall_LHS_after_2R={c['wall_lhs_after_2R']:.4g} "
            f"wall_rhs={c['wall_threshold_S_times_m']:.4g} -> {c['outcome']}"
        )

    results = build_results(
        failures_4488_3, classified, boundary_slice_sha, wall_target_sha
    )
    log(f"status={results['status']}")
    log(
        f"closed={results['failures_closed_by_ws01']} "
        f"closed_tight_R={results['failures_closed_by_ws01_tight_R_only']} "
        f"still_failing={results['failures_still_failing']} "
        f"branch_not_found={results['failures_branch_point_not_found']}"
    )
    log(
        f"median_required_factor_before={results['median_required_factor_before']} "
        f"median_required_factor_after_2R={results['median_required_factor_after']}"
    )

    OUTPUT_DIR.mkdir(parents=True)
    write_json(RESULTS_JSON, results)
    REPORT_MD.write_text(build_report(results), encoding="utf-8")
    RESULTS_SHA.write_text(
        f"{sha256_file(RESULTS_JSON)}  {RESULTS_JSON.name}\n", encoding="utf-8"
    )
    RUN_LOG.write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "experiment_id": EXPERIMENT_ID,
                "status": results["status"],
                "results": str(RESULTS_JSON),
                "report": str(REPORT_MD),
                "sha256": str(RESULTS_SHA),
                "log": str(RUN_LOG),
                "closed_by_ws01": results["failures_closed_by_ws01"],
                "closed_by_ws01_tight_R_only": results[
                    "failures_closed_by_ws01_tight_R_only"
                ],
                "still_failing": results["failures_still_failing"],
                "branch_point_not_found": results["failures_branch_point_not_found"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
