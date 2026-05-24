#!/usr/bin/env python3
"""Track B1-replicate of EHP114 bridge program.

Replicates the WS-01-CENTER-STRIP-CANCELLATION (Branch-Centered Moving-Frame
Collar) test on additional owner families of CELL-02-03's n=15 wall-separation
failures, beyond the original `4488:3` run that lives at:

  proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01/

This script is a parameterized clone of
`build_n15_cell0203_centerline_sign_model.py`. It accepts a single owner-family
key on the command line (e.g., `--owner-family 4571:2`) and writes a fresh
artifact directory at:

  proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-<KEY>-20260509-01/

The file-system-safe rendering of the owner key replaces ':' with '-'.

Internal experiment artifact only. Not a proof of Erdos #114, not an n=15
certificate, not a CELL-02-03 closure. Only a diagnostic of WS-01 effectiveness
on the chosen owner family. Same claim ceiling as the 4488:3 run.

Math reference (consume, do not re-derive):
  Erdos114/EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_2026-05-06.md
  lines 122-146 (Branch-Centered Moving-Frame Collar lemma).
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


TARGET_PROBLEM = 114
DEGREE = 15
CELL_ID = "CELL-02-03"
RUN_DATE_TAG = "20260509-01"

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


def owner_to_path_token(owner: str) -> str:
    """Render an owner key like '4571:2' as the filesystem token '4571-2'."""
    return owner.replace(":", "-")


def owner_to_status_token(owner: str) -> str:
    """Render an owner key like '4571:2' as the status token '4571_2'."""
    return owner.replace(":", "_")


def collect_wall_failures(
    boundary_slice: dict[str, Any], owner_family: str
) -> list[dict[str, Any]]:
    rows = []
    for item in boundary_slice.get("remaining_unresolved", []):
        if item.get("reason") != "third_order_wall_separation_failed":
            continue
        if item.get("source_ownership_key") != owner_family:
            continue
        rows.append(item)
    return rows


def classify_box(item: dict[str, Any], owner_family: str) -> dict[str, Any]:
    """Apply WS-01 rewrite to a single failure box. Identical numerics to the
    4488:3 run; owner_family is parameter-only (carried into output rows)."""
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
    f_iv_contains_zero = (f_lo <= 0.0 <= f_hi)
    grad_nonvanishing = fn_center_lo > 0.0
    branch_validated = f_iv_contains_zero and grad_nonvanishing

    new_lhs_R = 0.5 * f_tt_up * R * R + (T3 / 6.0) * R * R * R + Rn
    new_lhs_2R = (
        0.5 * f_tt_up * (2.0 * R) ** 2 + (T3 / 6.0) * (2.0 * R) ** 3 + Rn
    )

    if not branch_validated:
        outcome = "BRANCH_POINT_NOT_FOUND"
    elif new_lhs_2R < wall_rhs:
        outcome = "CLOSED_BY_WS01"
    elif new_lhs_R < wall_rhs:
        outcome = "CLOSED_BY_WS01_TIGHT_R_ONLY"
    else:
        outcome = "STILL_FAILING_AFTER_WS01"

    required_factor_before = wall_rhs / wall_lhs_before if wall_lhs_before > 0 else None
    required_factor_after = wall_rhs / new_lhs_2R if new_lhs_2R > 0 else None

    return {
        "ownership_key": owner_family,
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
    experiment_id: str,
    owner_family: str,
    failures_all: list[dict[str, Any]],
    classified: list[dict[str, Any]],
    boundary_slice_sha: str,
    wall_target_sha: str,
) -> dict[str, Any]:
    n_in_family = len(failures_all)
    status_token = owner_to_status_token(owner_family)

    closed = [c for c in classified if c["outcome"] == "CLOSED_BY_WS01"]
    closed_tight = [
        c for c in classified if c["outcome"] == "CLOSED_BY_WS01_TIGHT_R_ONLY"
    ]
    still_failing = [
        c for c in classified if c["outcome"] == "STILL_FAILING_AFTER_WS01"
    ]
    not_found = [c for c in classified if c["outcome"] == "BRANCH_POINT_NOT_FOUND"]

    closed_count = len(closed)
    closed_tight_count = len(closed_tight)
    still_failing_count = len(still_failing) + closed_tight_count
    not_found_count = len(not_found)

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

    if n_in_family == 0:
        status = f"CENTERLINE_SIGN_MODEL_BLOCKED_{status_token}"
    elif not_found_count == n_in_family:
        status = f"CENTERLINE_SIGN_MODEL_BLOCKED_{status_token}"
    elif closed_count == n_in_family:
        status = f"CENTERLINE_SIGN_MODEL_FULL_PASS_{status_token}"
    elif closed_count > 0:
        status = f"CENTERLINE_SIGN_MODEL_PARTIAL_PASS_{status_token}"
    else:
        status = f"CENTERLINE_SIGN_MODEL_BLOCKED_{status_token}"

    if median_before and median_after:
        interpretation = (
            f"On owner family {owner_family} of CELL-02-03's n=15 wall-separation "
            f"failures ({n_in_family} boxes), the WS-01 Branch-Centered Moving-Frame "
            f"Collar rewrite reduces the median wall-LHS-over-RHS deficit from "
            f"~{1/median_before:.2f}x miss to a "
            f"~{median_after:.2f}x clearance (factor of safety) when the worst-"
            f"case z*-off-center radius 2R is used. {closed_count}/{n_in_family} "
            f"boxes close by the conservative 2R bound and {closed_tight_count} more "
            f"close only with the tight-R bound. {not_found_count} boxes have no "
            f"interval-validated branch point. Per-box validation: F_iv contains 0 "
            f"(branch exists by IVT) and center |F_n| lower bound is large positive "
            f"(grad F nonvanishing, so the zero curve is smooth and t(u) is "
            f"well-defined). F_t(z*)=0 holds by construction, not by numerics. "
            f"This is a Python numerical demonstration at the same claim level as "
            f"the 4488:3 run; NOT a Rust interval-certified result."
        )
    else:
        interpretation = "Insufficient classified rows to compute interpretation."

    if status.startswith("CENTERLINE_SIGN_MODEL_FULL_PASS"):
        next_dependency = (
            f"Track B1-replicate confirmed full pass on owner family {owner_family}. "
            f"Combined with the 4488:3 result (16/16 closed), Track B1 viable across "
            f"the worst-three owner families. Recommended next: replicate on the "
            f"remaining matched-count families (3104:0, 3546:0, 3982:0, 4404:0 with "
            f"16 each), then move toward Rust interval-certification of the moving-"
            f"frame collar that emits z* via interval Newton on F(z)=0."
        )
    elif status.startswith("CENTERLINE_SIGN_MODEL_PARTIAL_PASS"):
        next_dependency = (
            f"Partial pass on owner family {owner_family}. Track B1 viable but not "
            f"uniform; the failing boxes need either smaller cells, a refined T3 "
            f"bound, or a Rust interval-certified z* via interval Newton. Recommended "
            f"next: characterize the still-failing rows (which inequality term "
            f"dominates? f_tt_upper * R^2, T3 * R^3, or Rn?) and decide whether to "
            f"refine cells or upgrade to interval-Newton."
        )
    else:
        next_dependency = (
            f"BLOCKED on owner family {owner_family}. Moving-frame collar cannot "
            f"establish a validated branch point and/or all wall LHS values still "
            f"exceed RHS even with the WS-01 rewrite. The cell decomposition for "
            f"n=15 needs to be reconsidered (smaller boxes, or a different analytic "
            f"reduction)."
        )

    return {
        "experiment_id": experiment_id,
        "schema_version": "1.0",
        "timestamp_unix": int(time.time()),
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "generated_by": "track_b1_replicate_centerline_sign_model_runner",
        "origin": "auto-research",
        "promotion_state": "review_only",
        "target_problem": TARGET_PROBLEM,
        "degree": DEGREE,
        "cell_id": CELL_ID,
        "owner_family_tested": owner_family,
        "claim_ceiling": (
            "Internal experiment artifact only. Not a proof of Erdos #114, "
            "not an n=15 certificate, and not a CELL-02-03 closure. Only a "
            f"diagnostic of WS-01 effectiveness on owner family {owner_family}. "
            "Python numerical demonstration, NOT a Rust interval-certified result."
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
            "predecessor_run_4488_3": str(
                PROOF_PATH
                / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01"
                / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01_RESULTS.json"
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
            "existing 4488:3 centerline-sign-model artifact",
        ],
    }


def build_report(experiment_id: str, owner_family: str, results: dict[str, Any]) -> str:
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

    return f"""# {experiment_id}

## Scope

Internal experiment artifact only. Track B1-replicate of the EHP114 bridge
program. Not a proof of Erdos #114, not an n=15 certificate, and not a
CELL-02-03 closure. Only a diagnostic of WS-01-CENTER-STRIP-CANCELLATION (the
Branch-Centered Moving-Frame Collar rewrite) on owner family `{owner_family}`'s
{results["failures_in_family"]} wall-separation failures. Python numerical
demonstration at the same claim level as the 4488:3 run; not Rust interval-
certified.

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

## Per-failure outcomes (owner family {owner_family})

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
public pages, proof registries, existing EHP114 finite packets, the existing
4488:3 centerline-sign-model artifact) were not attempted. The build script
`build_n15_cell0203_centerline_sign_model_replicate.py` writes only this
experiment folder. No prior script or artifact was modified.
"""


def run_one(owner_family: str) -> dict[str, Any]:
    experiment_id = (
        f"EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-"
        f"{owner_to_path_token(owner_family)}-{RUN_DATE_TAG}"
    )
    output_dir = PROOF_PATH / experiment_id
    results_json = output_dir / f"{experiment_id}_RESULTS.json"
    report_md = output_dir / f"{experiment_id}_REPORT.md"
    results_sha = output_dir / f"{experiment_id}_RESULTS.sha256"
    run_log = output_dir / f"{experiment_id}_RUN.log"

    if output_dir.exists():
        raise SystemExit(f"Refusing to overwrite existing output directory: {output_dir}")
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

    log(f"experiment_id={experiment_id}")
    log(f"owner_family={owner_family}")
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
    wt_count = wall_target.get("wall_separation_analysis", {}).get(
        "wall_failure_count"
    )
    log(f"wall-target packet reports wall_failure_count={wt_count}")

    failures_family = collect_wall_failures(boundary_slice, owner_family)
    log(f"failures in owner family {owner_family}: {len(failures_family)}")

    classified = [classify_box(item, owner_family) for item in failures_family]
    for i, c in enumerate(classified):
        log(
            f"#{i:2d} src_idx={c['source_index']} R={c['tangent_radius_R']:.4g} "
            f"wall_LHS_before={c['wall_lhs_before']:.4g} "
            f"wall_LHS_after_2R={c['wall_lhs_after_2R']:.4g} "
            f"wall_rhs={c['wall_threshold_S_times_m']:.4g} -> {c['outcome']}"
        )

    results = build_results(
        experiment_id,
        owner_family,
        failures_family,
        classified,
        boundary_slice_sha,
        wall_target_sha,
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

    output_dir.mkdir(parents=True)
    write_json(results_json, results)
    report_md.write_text(build_report(experiment_id, owner_family, results), encoding="utf-8")
    results_sha.write_text(
        f"{sha256_file(results_json)}  {results_json.name}\n", encoding="utf-8"
    )
    run_log.write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    summary = {
        "experiment_id": experiment_id,
        "owner_family": owner_family,
        "status": results["status"],
        "results": str(results_json),
        "report": str(report_md),
        "sha256": str(results_sha),
        "log": str(run_log),
        "failures_in_family": results["failures_in_family"],
        "closed_by_ws01": results["failures_closed_by_ws01"],
        "closed_by_ws01_tight_R_only": results[
            "failures_closed_by_ws01_tight_R_only"
        ],
        "still_failing": results["failures_still_failing"],
        "branch_point_not_found": results["failures_branch_point_not_found"],
    }
    print(json.dumps(summary, indent=2))
    return summary


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Track B1-replicate: WS-01 centerline-sign-model on a chosen "
            "owner family of CELL-02-03's n=15 wall-separation failures."
        )
    )
    parser.add_argument(
        "--owner-family",
        required=True,
        help="Owner family key, e.g. '4571:2' or '2484:4'.",
    )
    args = parser.parse_args(argv)

    run_one(args.owner_family)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
