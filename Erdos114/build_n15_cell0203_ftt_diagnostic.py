#!/usr/bin/env python3
"""F_tt distribution diagnostic on residual WS-01 wall-separation failures.

Internal diagnostic only. Not a proof of Erdos #114, not a CELL-02-03 closure.
Produces the F_tt distribution comparison the three reviewers (Perplexity,
Gemini, Claude) requested before the Modal Tier-2 Rust-port spend is committed.

Reviewer hypothesis under test
--------------------------------
The 8 residual wall-separation failures in CELL-02-03 (the 7 still failing at
any R bound from owner family 4571:2, plus 1 still failing at tight R from
owner family 3982:0) might NOT be recoverable by interval-Newton because the
boxes contain or are adjacent to a critical point of F where |F_tt| itself
nearly vanishes. If true, the moving-frame collar architecture is inherently
insufficient for those boxes — the linear-in-R term goes away by frame
construction (good), but the quadratic 0.5 |F_tt| R^2 term ALSO becomes small
or sign-flipped (bad), and the wall test fails for a deeper structural reason
than any z*-placement penalty.

What the existing Python pipeline records
------------------------------------------
The upstream boundary-slice runner (Rust binary
`ehp114_n15_boundary_slice_cell_02_03`, invoked by
`scripts/erdos-114/modal_n15_boundary_slice_cell_02_03.py`) emits per-box
`f_tt_abs_upper` (an interval upper bound on |F_tt|). It does NOT emit a
lower bound `f_tt_abs_lower`. The downstream Python build scripts
(`build_n15_cell0203_centerline_sign_model.py` and `..._replicate.py`) consume
that upper bound directly and never recompute |F_tt|.

Consequence for the diagnostic
-------------------------------
A `|F_tt|` UPPER bound by itself cannot directly confirm or refute the
critical-point hypothesis in the strong form the reviewers stated (interval
straddling zero requires both a lower and upper bound). It can, however,
disclose the WEAK form: if failing boxes had systematically SMALLER upper
bounds than closing boxes, the critical-point hypothesis would gain credence
because |F_tt| <= small upper bound is consistent with |F_tt| approaching
zero. If failing boxes have LARGER upper bounds than closing boxes, the
critical-point hypothesis loses credence because the second derivative is
clearly NOT shrinking on the failing boxes — the wall test is failing for
some other reason (most likely: small wall RHS / small `S * lower(|F_n|)`,
or large T3 cubic remainder, or large normal_remainder_bound R_n).

This script does NOT recompute |F_tt| lower bound. Recovering that would
require modifying the upstream Rust binary and re-running it, which is
outside scope for an honest day-of diagnostic. The verdict rules below
reflect that.

Verdict rules
--------------
- CONFIRMED: failing-box `f_tt_abs_upper` median is < 0.5 * closing-box
  median, AND failing-box minimum upper bound is < 0.25 * closing-box
  minimum upper bound. (Failing boxes have markedly smaller upper bounds,
  consistent with |F_tt| collapsing toward zero.)
- REFUTED: failing-box `f_tt_abs_upper` median is > 1.25 * closing-box
  median, AND every failing-box upper bound > closing-box minimum.
  (Failing boxes have UPPER bounds clearly bounded away from zero —
  even in the worst case |F_tt| does not vanish on these boxes.)
- INCONCLUSIVE: anything else, including the case where lower-bound
  data is unavailable and the reviewer hypothesis cannot be tested in
  its strong form.

Honest limitation: REFUTED here means "the upper bound does not collapse
toward zero on the failing boxes." It does NOT prove |F_tt| is uniformly
bounded BELOW on those boxes. To prove that, the Rust binary needs an
`f_tt_abs_lower` emitter. The Track B / Tier-2 Rust port should add that
emitter as a small free-rider task.

Same Python-numerical claim level as the existing WS-01 work.
NOT Rust interval-certified.
"""

from __future__ import annotations

import hashlib
import json
import statistics
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-N15-CELL-02-03-FTT-DIAGNOSTIC-20260509-01"
TARGET_PROBLEM = 114
DEGREE = 15
CELL_ID = "CELL-02-03"

SCRIPT_PATH = Path(__file__).resolve()
ERDOS114_DIR = SCRIPT_PATH.parent
PROOF_PATH = ERDOS114_DIR / "proof_path"

OUTPUT_DIR = PROOF_PATH / EXPERIMENT_ID
RESULTS_JSON = OUTPUT_DIR / f"{EXPERIMENT_ID}_RESULTS.json"
REPORT_MD = OUTPUT_DIR / f"{EXPERIMENT_ID}_REPORT.md"
RESULTS_SHA = OUTPUT_DIR / f"{EXPERIMENT_ID}_RESULTS.sha256"


SOURCE_4571_2 = (
    PROOF_PATH
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4571-2-20260509-01"
)
SOURCE_4571_2_JSON = (
    SOURCE_4571_2
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4571-2-20260509-01_RESULTS.json"
)
SOURCE_4571_2_SHA = (
    SOURCE_4571_2
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4571-2-20260509-01_RESULTS.sha256"
)

SOURCE_3982_0 = (
    PROOF_PATH
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3982-0-20260509-01"
)
SOURCE_3982_0_JSON = (
    SOURCE_3982_0
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3982-0-20260509-01_RESULTS.json"
)
SOURCE_3982_0_SHA = (
    SOURCE_3982_0
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3982-0-20260509-01_RESULTS.sha256"
)

SOURCE_4488_3 = (
    PROOF_PATH
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01"
)
SOURCE_4488_3_JSON = (
    SOURCE_4488_3
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01_RESULTS.json"
)
SOURCE_4488_3_SHA = (
    SOURCE_4488_3
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01_RESULTS.sha256"
)

SOURCE_2484_4 = (
    PROOF_PATH
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2484-4-20260509-01"
)
SOURCE_2484_4_JSON = (
    SOURCE_2484_4
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2484-4-20260509-01_RESULTS.json"
)
SOURCE_2484_4_SHA = (
    SOURCE_2484_4
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2484-4-20260509-01_RESULTS.sha256"
)

SOURCE_2932_1 = (
    PROOF_PATH
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2932-1-20260509-01"
)
SOURCE_2932_1_JSON = (
    SOURCE_2932_1
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2932-1-20260509-01_RESULTS.json"
)
SOURCE_2932_1_SHA = (
    SOURCE_2932_1
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2932-1-20260509-01_RESULTS.sha256"
)

SOURCE_4404_0 = (
    PROOF_PATH
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4404-0-20260509-01"
)
SOURCE_4404_0_JSON = (
    SOURCE_4404_0
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4404-0-20260509-01_RESULTS.json"
)
SOURCE_4404_0_SHA = (
    SOURCE_4404_0
    / "EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4404-0-20260509-01_RESULTS.sha256"
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


def collect_failing_boxes() -> list[dict[str, Any]]:
    """The 8 residual failures (= boxes that don't pass the conservative
    FULL_PASS-2R criterion): 7 from 4571:2 (6 STILL_FAILING + 1 TIGHT_R_ONLY)
    plus 1 from 3982:0 (STILL_FAILING_AFTER_WS01).

    Note: per the build script, `failures_still_failing` counts both
    STILL_FAILING_AFTER_WS01 and CLOSED_BY_WS01_TIGHT_R_ONLY. They are
    "still failing at any R bound" in the sense that the 2R conservative
    bound does not clear them. The 1 box from 3982:0 grouped here is the
    truly STILL_FAILING_AFTER_WS01 entry — the 11 TIGHT_R_ONLY boxes from
    3982:0 are excluded because the user's prompt explicitly groups them
    separately ("1 still failing at tight R" = the harder 3982:0 case)."""
    failing: list[dict[str, Any]] = []

    d4571 = load_json(SOURCE_4571_2_JSON)
    for pfo in d4571["per_failure_outcomes"]:
        if pfo["outcome"] in ("STILL_FAILING_AFTER_WS01", "CLOSED_BY_WS01_TIGHT_R_ONLY"):
            failing.append(
                {
                    "owner_key": "4571:2",
                    "src_idx": pfo["source_index"],
                    "split_path": pfo.get("split_path"),
                    "outcome": pfo["outcome"],
                    "f_tt_abs_upper": pfo["f_tt_abs_upper"],
                    "tangent_radius_R": pfo["tangent_radius_R"],
                    "third_directional_upper_T3": pfo["third_directional_upper_T3"],
                    "wall_lhs_after_2R": pfo["wall_lhs_after_2R"],
                    "wall_threshold_S_times_m": pfo["wall_threshold_S_times_m"],
                    "required_factor_after_rhs_over_lhs_2R": pfo[
                        "required_factor_after_rhs_over_lhs_2R"
                    ],
                    "x_interval": pfo.get("x_interval"),
                    "y_interval": pfo.get("y_interval"),
                }
            )

    d3982 = load_json(SOURCE_3982_0_JSON)
    for pfo in d3982["per_failure_outcomes"]:
        if pfo["outcome"] == "STILL_FAILING_AFTER_WS01":
            failing.append(
                {
                    "owner_key": "3982:0",
                    "src_idx": pfo["source_index"],
                    "split_path": pfo.get("split_path"),
                    "outcome": pfo["outcome"],
                    "f_tt_abs_upper": pfo["f_tt_abs_upper"],
                    "tangent_radius_R": pfo["tangent_radius_R"],
                    "third_directional_upper_T3": pfo["third_directional_upper_T3"],
                    "wall_lhs_after_2R": pfo["wall_lhs_after_2R"],
                    "wall_threshold_S_times_m": pfo["wall_threshold_S_times_m"],
                    "required_factor_after_rhs_over_lhs_2R": pfo[
                        "required_factor_after_rhs_over_lhs_2R"
                    ],
                    "x_interval": pfo.get("x_interval"),
                    "y_interval": pfo.get("y_interval"),
                }
            )

    return failing


def collect_closing_boxes_baseline() -> list[dict[str, Any]]:
    """Representative subset of CLOSED_BY_WS01 (FULL_PASS-2R) boxes from the
    four closing families: 4488:3 (16), 2484:4, 2932:1, 4404:0."""
    closing: list[dict[str, Any]] = []
    for owner_key, src_path in [
        ("4488:3", SOURCE_4488_3_JSON),
        ("2484:4", SOURCE_2484_4_JSON),
        ("2932:1", SOURCE_2932_1_JSON),
        ("4404:0", SOURCE_4404_0_JSON),
    ]:
        d = load_json(src_path)
        owner_labeled = d.get("owner_family_tested", owner_key)
        for pfo in d["per_failure_outcomes"]:
            if pfo["outcome"] == "CLOSED_BY_WS01":
                closing.append(
                    {
                        "owner_key": owner_labeled,
                        "src_idx": pfo["source_index"],
                        "split_path": pfo.get("split_path"),
                        "outcome": pfo["outcome"],
                        "f_tt_abs_upper": pfo["f_tt_abs_upper"],
                        "tangent_radius_R": pfo["tangent_radius_R"],
                        "third_directional_upper_T3": pfo[
                            "third_directional_upper_T3"
                        ],
                        "wall_lhs_after_2R": pfo["wall_lhs_after_2R"],
                        "wall_threshold_S_times_m": pfo[
                            "wall_threshold_S_times_m"
                        ],
                        "required_factor_after_rhs_over_lhs_2R": pfo[
                            "required_factor_after_rhs_over_lhs_2R"
                        ],
                    }
                )
    return closing


def select_representative_subset(boxes: list[dict[str, Any]], n: int = 16) -> list[dict[str, Any]]:
    """Stratified pick across owner families to keep the closing-baseline at
    ~16 entries while preserving family coverage."""
    by_owner: dict[str, list[dict[str, Any]]] = {}
    for b in boxes:
        by_owner.setdefault(b["owner_key"], []).append(b)
    # take min(4, len(family)) per family, then pad to n
    rep: list[dict[str, Any]] = []
    for owner, items in by_owner.items():
        rep.extend(items[: min(4, len(items))])
    if len(rep) < n:
        # pad with extras
        seen = {(b["owner_key"], b["src_idx"], b.get("split_path")) for b in rep}
        for owner, items in by_owner.items():
            for it in items[min(4, len(items)) :]:
                key = (it["owner_key"], it["src_idx"], it.get("split_path"))
                if key not in seen and len(rep) < n:
                    rep.append(it)
                    seen.add(key)
    return rep[:n]


def compute_stats(values: list[float]) -> dict[str, Any]:
    if not values:
        return {"n": 0}
    sorted_v = sorted(values)
    return {
        "n": len(values),
        "min": sorted_v[0],
        "max": sorted_v[-1],
        "median": statistics.median(values),
        "mean": statistics.mean(values),
        "p25": sorted_v[len(values) // 4],
        "p75": sorted_v[(3 * len(values)) // 4],
    }


def render_verdict(comparison: dict[str, Any]) -> str:
    fail_med = comparison["failing_F_tt_median"]
    close_med = comparison["closing_F_tt_median"]
    fail_min = comparison["failing_F_tt_min"]
    close_min = comparison["closing_F_tt_min"]
    ratio = comparison["ratio_of_medians_failing_over_closing"]

    if fail_med < 0.5 * close_med and fail_min < 0.25 * close_min:
        return "FTT_DIAGNOSTIC_CRITICAL_POINT_HYPOTHESIS_CONFIRMED"
    if fail_med > 1.25 * close_med and fail_min > close_min:
        return "FTT_DIAGNOSTIC_CRITICAL_POINT_HYPOTHESIS_REFUTED"
    return "FTT_DIAGNOSTIC_INCONCLUSIVE"


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Source SHA verification
    sha_status = []
    for label, jpath, spath in [
        ("4571_2", SOURCE_4571_2_JSON, SOURCE_4571_2_SHA),
        ("3982_0", SOURCE_3982_0_JSON, SOURCE_3982_0_SHA),
        ("4488_3", SOURCE_4488_3_JSON, SOURCE_4488_3_SHA),
        ("2484_4", SOURCE_2484_4_JSON, SOURCE_2484_4_SHA),
        ("2932_1", SOURCE_2932_1_JSON, SOURCE_2932_1_SHA),
        ("4404_0", SOURCE_4404_0_JSON, SOURCE_4404_0_SHA),
    ]:
        rec_sha = sha_from_sidecar(spath)
        live_sha = sha256_file(jpath) if jpath.exists() else None
        sha_status.append(
            {
                "owner_key_label": label,
                "json_path": str(jpath),
                "sidecar_sha": rec_sha,
                "live_sha": live_sha,
                "match": (rec_sha == live_sha) if rec_sha and live_sha else None,
            }
        )

    failing_boxes = collect_failing_boxes()
    closing_all = collect_closing_boxes_baseline()
    closing_rep = select_representative_subset(closing_all, n=16)

    fail_ftt = [b["f_tt_abs_upper"] for b in failing_boxes]
    close_ftt = [b["f_tt_abs_upper"] for b in closing_all]

    fail_stats = compute_stats(fail_ftt)
    close_stats = compute_stats(close_ftt)

    comparison_stats = {
        "failing_F_tt_median": fail_stats.get("median"),
        "closing_F_tt_median": close_stats.get("median"),
        "failing_F_tt_min": fail_stats.get("min"),
        "closing_F_tt_min": close_stats.get("min"),
        "failing_F_tt_max": fail_stats.get("max"),
        "closing_F_tt_max": close_stats.get("max"),
        "failing_F_tt_mean": fail_stats.get("mean"),
        "closing_F_tt_mean": close_stats.get("mean"),
        "ratio_of_medians_failing_over_closing": (
            fail_stats["median"] / close_stats["median"]
            if close_stats.get("median")
            else None
        ),
        "any_failing_box_F_tt_straddles_zero": (
            None
        ),  # cannot determine: only upper bound recorded by upstream Rust binary
        "lower_bound_data_available": False,
        "lower_bound_data_unavailable_reason": (
            "Upstream Rust binary (ehp114_n15_boundary_slice_cell_02_03) emits "
            "only f_tt_abs_upper, not f_tt_abs_lower. Recovering the lower bound "
            "requires modifying that binary and re-running the boundary-slice. "
            "Out of scope for an honest day-of Python diagnostic."
        ),
    }

    status = render_verdict(comparison_stats)

    # Build the failing-box records with explicit straddle-zero field (= None
    # when only upper bound is known)
    failing_records = []
    for b in failing_boxes:
        ftt_up = b["f_tt_abs_upper"]
        failing_records.append(
            {
                "owner_key": b["owner_key"],
                "src_idx": b["src_idx"],
                "split_path": b["split_path"],
                "F_tt_interval": [None, ftt_up],
                "F_tt_upper_bound": ftt_up,
                "F_tt_lower_bound_known": False,
                "F_tt_straddles_zero": None,
                "F_tt_magnitude_upper": ftt_up,
                "tangent_radius_R": b["tangent_radius_R"],
                "third_directional_upper_T3": b["third_directional_upper_T3"],
                "wall_lhs_after_2R": b["wall_lhs_after_2R"],
                "wall_threshold_S_times_m": b["wall_threshold_S_times_m"],
                "required_factor_after_rhs_over_lhs_2R": b[
                    "required_factor_after_rhs_over_lhs_2R"
                ],
                "x_interval": b["x_interval"],
                "y_interval": b["y_interval"],
            }
        )

    closing_records = []
    for b in closing_rep:
        ftt_up = b["f_tt_abs_upper"]
        closing_records.append(
            {
                "owner_key": b["owner_key"],
                "src_idx": b["src_idx"],
                "split_path": b["split_path"],
                "F_tt_interval": [None, ftt_up],
                "F_tt_upper_bound": ftt_up,
                "F_tt_magnitude_upper": ftt_up,
            }
        )

    # Verdict prose
    fail_med = comparison_stats["failing_F_tt_median"]
    close_med = comparison_stats["closing_F_tt_median"]
    ratio = comparison_stats["ratio_of_medians_failing_over_closing"]
    if status == "FTT_DIAGNOSTIC_CRITICAL_POINT_HYPOTHESIS_REFUTED":
        verdict_prose = (
            f"REFUTED (weak form). The 8 residual failing boxes have a median "
            f"|F_tt| upper bound of {fail_med:.1f}, vs {close_med:.1f} for "
            f"closing-family boxes (ratio {ratio:.2f}x). Failing boxes do NOT "
            f"have collapsing |F_tt| upper bounds — if anything they are "
            f"larger. The wall test is failing for some other reason "
            f"(small wall RHS / S * lower(|F_n|), large T3 cubic remainder, "
            f"or large normal_remainder R_n). Modal Tier-2 Rust port + "
            f"interval-Newton recovery story remains plausible. STRONG-FORM "
            f"REFUTATION (proving |F_tt| bounded uniformly BELOW on the "
            f"failing boxes) requires adding f_tt_abs_lower to the upstream "
            f"Rust binary — recommend folding that into Tier-2 as a small "
            f"free-rider task."
        )
    elif status == "FTT_DIAGNOSTIC_CRITICAL_POINT_HYPOTHESIS_CONFIRMED":
        verdict_prose = (
            f"CONFIRMED. Failing boxes have systematically smaller |F_tt| "
            f"upper bounds (median {fail_med:.1f} vs {close_med:.1f}, "
            f"ratio {ratio:.2f}x). Critical-point hypothesis gains support. "
            f"Modal Tier-2 Rust port should NOT be the next dependency; "
            f"need symmetry-reduction or analytic critical-point exclusion."
        )
    else:
        verdict_prose = (
            f"INCONCLUSIVE on strong form (lower bound unavailable). Weak-form "
            f"signal: failing-box |F_tt| upper-bound median {fail_med:.1f} vs "
            f"closing {close_med:.1f}, ratio {ratio:.2f}x. Recommend Tier-2 "
            f"Rust port include an f_tt_abs_lower emitter to settle the "
            f"strong form before committing to symmetry-reduction work."
        )

    # Implications block
    implications = {}
    if status == "FTT_DIAGNOSTIC_CRITICAL_POINT_HYPOTHESIS_REFUTED":
        implications = {
            "modal_tier_2_spend_justified": "conditional",
            "next_dependency": (
                "Tier-2 Rust port with f_tt_abs_lower emitter folded in. "
                "Lower bound proves strong-form refutation; without it, "
                "the weak-form refutation here only shows the upper bound "
                "is not collapsing, which is necessary but not sufficient."
            ),
        }
    elif status == "FTT_DIAGNOSTIC_CRITICAL_POINT_HYPOTHESIS_CONFIRMED":
        implications = {
            "modal_tier_2_spend_justified": False,
            "next_dependency": (
                "Symmetry-reduction (Gemini's suggestion) or analytic "
                "critical-point exclusion. Rust port + interval-Newton "
                "alone will not recover these 8 boxes."
            ),
        }
    else:
        implications = {
            "modal_tier_2_spend_justified": "conditional",
            "next_dependency": (
                "Modify upstream Rust binary to emit f_tt_abs_lower. Re-run "
                "boundary-slice on failing boxes only. Then re-run this "
                "diagnostic to settle the strong form."
            ),
        }

    payload = {
        "experiment_id": EXPERIMENT_ID,
        "schema_version": "1.0",
        "status": status,
        "claim_ceiling": (
            "Internal diagnostic only. Not a proof of Erdos #114, not a "
            "CELL-02-03 closure. Just the F_tt upper-bound distribution on "
            "8 residual failing boxes vs closing-family baseline. Same "
            "Python-numerical claim level as the existing WS-01 work; not "
            "Rust interval-certified. Lower bound unavailable upstream."
        ),
        "hypothesis": (
            "8 residual WS-01 failures contain or are adjacent to critical "
            "points where |F_tt| vanishes, making moving-frame collar "
            "inherently insufficient."
        ),
        "data_availability": {
            "f_tt_abs_upper": "AVAILABLE (recorded by upstream Rust binary, "
            "consumed verbatim by Python build scripts)",
            "f_tt_abs_lower": "NOT_AVAILABLE (upstream Rust binary does not "
            "emit this field; recomputing requires modifying that binary)",
        },
        "failing_boxes": failing_records,
        "closing_boxes_baseline": closing_records,
        "comparison_stats": comparison_stats,
        "verdict": verdict_prose,
        "implications": implications,
        "source_paths": {
            "owner_4571_2": str(SOURCE_4571_2_JSON),
            "owner_3982_0": str(SOURCE_3982_0_JSON),
            "owner_4488_3_aka_default_run": str(SOURCE_4488_3_JSON),
            "owner_2484_4": str(SOURCE_2484_4_JSON),
            "owner_2932_1": str(SOURCE_2932_1_JSON),
            "owner_4404_0": str(SOURCE_4404_0_JSON),
        },
        "source_sha_status": sha_status,
        "target_problem": TARGET_PROBLEM,
        "degree": DEGREE,
        "cell_id": CELL_ID,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generated_by": str(SCRIPT_PATH),
        "timestamp_unix": int(time.time()),
    }

    write_json(RESULTS_JSON, payload)

    # SHA256 sidecar
    rj_sha = sha256_file(RESULTS_JSON)
    RESULTS_SHA.write_text(f"{rj_sha}  {RESULTS_JSON.name}\n", encoding="utf-8")

    # Report
    lines: list[str] = []
    lines.append(f"# {EXPERIMENT_ID}")
    lines.append("")
    lines.append("## Status")
    lines.append("")
    lines.append(f"- {status}")
    lines.append("")
    lines.append("## Claim Ceiling")
    lines.append("")
    lines.append(payload["claim_ceiling"])
    lines.append("")
    lines.append("## Hypothesis Under Test")
    lines.append("")
    lines.append(payload["hypothesis"])
    lines.append("")
    lines.append("## Data Availability")
    lines.append("")
    lines.append(
        "- `f_tt_abs_upper`: AVAILABLE — recorded per box by the upstream "
        "Rust binary `ehp114_n15_boundary_slice_cell_02_03` and consumed "
        "verbatim by the Python build scripts."
    )
    lines.append(
        "- `f_tt_abs_lower`: NOT AVAILABLE — the upstream Rust binary does "
        "not emit this field. Recovering it requires modifying the binary "
        "and re-running the boundary-slice. Out of scope for this Python "
        "diagnostic."
    )
    lines.append("")
    lines.append(
        "Consequence: this diagnostic can only test the WEAK form of the "
        "critical-point hypothesis (does the upper bound collapse on failing "
        "boxes?). The STRONG form (does |F_tt| straddle zero?) needs the "
        "lower bound and is reported as `null` in the output."
    )
    lines.append("")
    lines.append("## Comparison: Failing vs Closing |F_tt| Upper Bounds")
    lines.append("")
    lines.append("| Group | n | min | median | mean | max |")
    lines.append("|---|---|---|---|---|---|")
    lines.append(
        f"| Failing (residual at FULL_PASS-2R) | {fail_stats['n']} | "
        f"{fail_stats['min']:.1f} | {fail_stats['median']:.1f} | "
        f"{fail_stats['mean']:.1f} | {fail_stats['max']:.1f} |"
    )
    lines.append(
        f"| Closing (FULL_PASS-2R) | {close_stats['n']} | "
        f"{close_stats['min']:.1f} | {close_stats['median']:.1f} | "
        f"{close_stats['mean']:.1f} | {close_stats['max']:.1f} |"
    )
    lines.append(
        f"| Ratio (failing / closing) | — | "
        f"{fail_stats['min']/close_stats['min']:.2f} | "
        f"{fail_stats['median']/close_stats['median']:.2f} | "
        f"{fail_stats['mean']/close_stats['mean']:.2f} | "
        f"{fail_stats['max']/close_stats['max']:.2f} |"
    )
    lines.append("")
    lines.append("## Per-Failing-Box |F_tt| Upper Bound")
    lines.append("")
    lines.append(
        "| owner_key | src_idx | split_path | F_tt upper | wall_lhs_after_2R | "
        "wall_RHS | req_factor (RHS/LHS) |"
    )
    lines.append("|---|---|---|---|---|---|---|")
    for b in failing_records:
        lines.append(
            f"| {b['owner_key']} | {b['src_idx']} | "
            f"{b['split_path']} | {b['F_tt_upper_bound']:.1f} | "
            f"{b['wall_lhs_after_2R']:.6e} | "
            f"{b['wall_threshold_S_times_m']:.6e} | "
            f"{b['required_factor_after_rhs_over_lhs_2R']:.4f} |"
        )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(verdict_prose)
    lines.append("")
    lines.append("## Implications for Modal Tier-2 Spend")
    lines.append("")
    lines.append(
        f"- **Spend justified:** {implications['modal_tier_2_spend_justified']}"
    )
    lines.append(f"- **Next dependency:** {implications['next_dependency']}")
    lines.append("")
    lines.append("## Source Paths")
    lines.append("")
    for k, v in payload["source_paths"].items():
        lines.append(f"- `{k}`: `{v}`")
    lines.append("")
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[ok] wrote {RESULTS_JSON}")
    print(f"[ok] wrote {REPORT_MD}")
    print(f"[ok] wrote {RESULTS_SHA}")
    print(f"[ok] status: {status}")
    print(f"[ok] failing F_tt upper median: {fail_stats['median']:.1f}")
    print(f"[ok] closing F_tt upper median: {close_stats['median']:.1f}")
    print(f"[ok] ratio (failing/closing) of medians: {ratio:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
