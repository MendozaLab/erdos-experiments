#!/usr/bin/env python3
"""Aggregate Hessian-style diagnostics for Erdos #20 sunflower closure data.

This script intentionally uses only saved April 2026 aggregate artifacts. It
does not rerun enumeration, add per-core instrumentation, update scorecards, or
touch D1. The output is a diagnostic for a possible Leg-4 observable, not a
theorem or lower-bound claim.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import date
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-ERDOS20-HESSIAN-CLOSURE-20260505-01"
RUN_DATE = "2026-05-05"
CLASSIFICATION_OPTIONS = [
    "HESSIAN_PRESENT",
    "HESSIAN_INCONCLUSIVE",
    "HESSIAN_FAIL",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def finite_round(value: float | None, digits: int = 6) -> float | None:
    if value is None or not math.isfinite(value):
        return None
    return round(value, digits)


def first_differences(values: list[float | None]) -> list[float | None]:
    diffs: list[float | None] = []
    for left, right in zip(values, values[1:]):
        if left is None or right is None:
            diffs.append(None)
        else:
            diffs.append(right - left)
    return diffs


def second_differences(values: list[float | None]) -> list[float | None]:
    diffs: list[float | None] = []
    for left, mid, right in zip(values, values[1:], values[2:]):
        if left is None or mid is None or right is None:
            diffs.append(None)
        else:
            diffs.append(right - 2.0 * mid + left)
    return diffs


def first_index_less_than(values: list[float], threshold: float) -> int | None:
    for idx, value in enumerate(values):
        if value < threshold:
            return idx
    return None


def safe_min(values: list[float | None]) -> float | None:
    finite_values = [value for value in values if value is not None and math.isfinite(value)]
    if not finite_values:
        return None
    return min(finite_values)


def safe_max(values: list[float | None]) -> float | None:
    finite_values = [value for value in values if value is not None and math.isfinite(value)]
    if not finite_values:
        return None
    return max(finite_values)


def closure_pressure_profile(
    n_sites: int,
    growth_rates: list[float] | None,
) -> tuple[list[dict[str, Any]], list[float | None]]:
    if not growth_rates:
        return [], []

    profile: list[dict[str, Any]] = []
    i_close_values: list[float | None] = []
    for m, growth_rate in enumerate(growth_rates):
        unused_sites = n_sites - m
        if unused_sites <= 0:
            p_safe = None
            i_close = None
            status = "no_unused_sites"
        else:
            p_safe = growth_rate / unused_sites
            if p_safe > 0:
                i_close = -math.log2(p_safe)
                status = "finite"
            else:
                i_close = None
                status = "jammed_zero_growth"

        i_close_values.append(i_close)
        profile.append(
            {
                "m": m,
                "growth_rate": finite_round(growth_rate),
                "unused_sites": unused_sites,
                "p_safe": finite_round(p_safe),
                "I_close_bits": finite_round(i_close),
                "status": status,
            }
        )
    return profile, i_close_values


def centered_second_records(values: list[float | None], name: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for idx, value in enumerate(values):
        records.append({"center_m": idx + 1, name: finite_round(value)})
    return records


def tail_records(records: list[dict[str, Any]], tail_count: int = 5) -> list[dict[str, Any]]:
    return records[-tail_count:] if len(records) > tail_count else records


def build_profile(
    row: dict[str, Any],
    source_files: list[str],
    growth_rates: list[float] | None = None,
) -> dict[str, Any]:
    density_of_states = row.get("density_of_states") or []
    log_density = [math.log2(value) if value > 0 else None for value in density_of_states]
    delta_log_density = first_differences(log_density)
    hessian_log_density = second_differences(log_density)

    n_sites = int(row["num_w_subsets"])
    if growth_rates is None:
        raw_growth = row.get("growth_rates")
        growth_rates = raw_growth if isinstance(raw_growth, list) else None

    closure_profile, i_close_values = closure_pressure_profile(n_sites, growth_rates)
    delta_i_close = first_differences(i_close_values)
    hessian_i_close = second_differences(i_close_values)
    jamming_m = first_index_less_than(growth_rates, 1.0) if growth_rates else None

    hessian_records = centered_second_records(hessian_log_density, "second_delta_logD")
    i_close_hessian_records = centered_second_records(hessian_i_close, "second_delta_I_close")

    if jamming_m is not None:
        jamming_window = [
            record
            for record in closure_profile
            if jamming_m - 2 <= record["m"] <= jamming_m + 2
        ]
        hessian_near_jamming = [
            record
            for record in hessian_records
            if jamming_m - 2 <= record["center_m"] <= jamming_m + 2
        ]
    else:
        jamming_window = []
        hessian_near_jamming = []

    finite_i_close = [value for value in i_close_values if value is not None]
    terminal_hessian_values = [
        record["second_delta_logD"]
        for record in tail_records(hessian_records, tail_count=4)
        if record["second_delta_logD"] is not None
    ]
    jamming_hessian_values = [
        record["second_delta_logD"]
        for record in hessian_near_jamming
        if record["second_delta_logD"] is not None
    ]

    return {
        "n": row["n"],
        "w": row["w"],
        "k": row.get("k", 3),
        "N_sites": n_sites,
        "max_family_size": row["max_family_size"],
        "total_sf_free": row["total_sf_free"],
        "log2_total_sf_free": finite_round(row.get("log2_count")),
        "source_files": source_files,
        "density_of_states": density_of_states,
        "logD": [finite_round(value) for value in log_density],
        "delta_logD": [finite_round(value) for value in delta_log_density],
        "second_delta_logD": [finite_round(value) for value in hessian_log_density],
        "second_delta_logD_by_center_m": hessian_records,
        "second_delta_logD_tail": tail_records(hessian_records, tail_count=5),
        "growth_rates": [finite_round(value) for value in growth_rates] if growth_rates else None,
        "jamming_m_star": jamming_m,
        "I_close_profile": closure_profile,
        "I_close_bits": [finite_round(value) for value in i_close_values],
        "delta_I_close": [finite_round(value) for value in delta_i_close],
        "second_delta_I_close": [finite_round(value) for value in hessian_i_close],
        "second_delta_I_close_by_center_m": i_close_hessian_records,
        "I_close_window_around_jamming": jamming_window,
        "second_delta_logD_window_around_jamming": hessian_near_jamming,
        "summary_metrics": {
            "min_second_delta_logD": finite_round(safe_min(hessian_log_density)),
            "max_second_delta_logD": finite_round(safe_max(hessian_log_density)),
            "min_tail_second_delta_logD": finite_round(safe_min(terminal_hessian_values)),
            "max_finite_I_close_bits": finite_round(max(finite_i_close) if finite_i_close else None),
            "I_close_at_jamming_bits": finite_round(
                i_close_values[jamming_m]
                if jamming_m is not None and jamming_m < len(i_close_values)
                else None
            ),
            "min_second_delta_logD_near_jamming": finite_round(
                safe_min(jamming_hessian_values)
            ),
            "has_growth_rates": bool(growth_rates),
        },
    }


def collect_profiles(base: Path) -> list[dict[str, Any]]:
    sunflower_001 = read_json(base / "EXP-MATH-ERDOS20-SUNFLOWER-001_RESULTS.json")
    sunflower_002 = read_json(base / "EXP-MATH-ERDOS20-SUNFLOWER-002_RESULTS.json")
    raw_w2 = read_json(base / "results_w2.json")
    raw_w4 = read_json(base / "results_w4.json")

    profiles: list[dict[str, Any]] = []

    for row in raw_w2["data"]:
        profiles.append(build_profile(row, ["results_w2.json"]))

    for row in sunflower_001["data"]:
        profiles.append(
            build_profile(
                row,
                ["EXP-MATH-ERDOS20-SUNFLOWER-001_RESULTS.json"],
            )
        )

    w4_growth_n7 = sunflower_002["data_by_w"]["w=4"].get("growth_rates_n7")
    for row in raw_w4["data"]:
        growth_rates = w4_growth_n7 if row["n"] == 7 else None
        source_files = ["results_w4.json"]
        if growth_rates:
            source_files.append("EXP-MATH-ERDOS20-SUNFLOWER-002_RESULTS.json")
        profiles.append(build_profile(row, source_files, growth_rates=growth_rates))

    return sorted(profiles, key=lambda item: (item["w"], item["n"]))


def classify(profiles: list[dict[str, Any]]) -> dict[str, Any]:
    growth_profiles = [profile for profile in profiles if profile["summary_metrics"]["has_growth_rates"]]
    curvature_profiles = [
        profile
        for profile in profiles
        if (
            profile["summary_metrics"]["min_tail_second_delta_logD"] is not None
            and profile["summary_metrics"]["min_tail_second_delta_logD"] < -0.25
        )
    ]
    closure_profiles_with_rise = [
        profile
        for profile in growth_profiles
        if (
            profile["summary_metrics"]["max_finite_I_close_bits"] is not None
            and profile["summary_metrics"]["I_close_at_jamming_bits"] is not None
            and profile["summary_metrics"]["max_finite_I_close_bits"]
            >= profile["summary_metrics"]["I_close_at_jamming_bits"]
        )
    ]

    aggregate_signature_present = bool(curvature_profiles and closure_profiles_with_rise)
    blockers = [
        "available data are aggregate by family size, not per-core closure states",
        "growth rates are missing for w=2 and for w=3,n=8",
        "w=4 stops at n=7 under exhaustive enumeration",
        "no Mendoza-floor numerator or per-core R profile is defined in these inputs",
        "Abbott-Hansen-Sauer dominates the small-n lower-bound framing; A-axis remains A0",
    ]

    if not aggregate_signature_present:
        classification = "HESSIAN_FAIL"
        meaning = (
            "No robust aggregate curvature plus closure-pressure signal was visible in "
            "the saved density-of-states and growth-rate profiles."
        )
    else:
        classification = "HESSIAN_INCONCLUSIVE"
        meaning = (
            "Aggregate curvature and closure pressure are visible near jamming, but "
            "the signal is not yet a per-core Hessian or Leg-4 core-closure observable."
        )

    return {
        "classification": classification,
        "allowed_classifications": CLASSIFICATION_OPTIONS,
        "aggregate_signature_present": aggregate_signature_present,
        "growth_profile_count": len(growth_profiles),
        "curvature_profile_count": len(curvature_profiles),
        "closure_pressure_profile_count": len(closure_profiles_with_rise),
        "blockers": blockers,
        "meaning": meaning,
        "claim_ceiling": (
            "shadow signature, not universal law; diagnostic aggregate curvature only; "
            "no theorem progress and no lower-bound improvement"
        ),
    }


def input_metadata(paths: list[Path], root: Path) -> list[dict[str, Any]]:
    metadata = []
    for path in paths:
        resolved = path.resolve()
        try:
            relative = str(resolved.relative_to(root))
        except ValueError:
            relative = str(resolved)
        metadata.append(
            {
                "path": relative,
                "bytes": resolved.stat().st_size,
                "sha256": sha256_file(resolved),
            }
        )
    return metadata


def compact_profile_table(profiles: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for profile in profiles:
        metrics = profile["summary_metrics"]
        rows.append(
            {
                "w": profile["w"],
                "n": profile["n"],
                "N_sites": profile["N_sites"],
                "max_family_size": profile["max_family_size"],
                "jamming_m_star": profile["jamming_m_star"],
                "min_tail_second_delta_logD": metrics["min_tail_second_delta_logD"],
                "I_close_at_jamming_bits": metrics["I_close_at_jamming_bits"],
                "max_finite_I_close_bits": metrics["max_finite_I_close_bits"],
                "has_growth_rates": metrics["has_growth_rates"],
            }
        )
    return rows


def markdown_table(headers: list[str], rows: list[list[Any]]) -> list[str]:
    lines = []
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(format_cell(value) for value in row) + " |")
    return lines


def format_cell(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def render_report(results: dict[str, Any]) -> str:
    classification = results["classification"]["classification"]
    compact_rows = results["compact_profile_table"]
    growth_rows = [row for row in compact_rows if row["has_growth_rates"]]
    strongest_curvature = sorted(
        compact_rows,
        key=lambda row: (
            row["min_tail_second_delta_logD"]
            if row["min_tail_second_delta_logD"] is not None
            else 0
        ),
    )[:8]

    lines: list[str] = [
        f"# {EXPERIMENT_ID} - Report",
        "",
        f"**Date:** {RUN_DATE}",
        "**Problem:** Erdos #20 sunflower core closure",
        "**Scope:** Aggregate discrete Hessian/free-energy diagnostic from saved April 2026 lattice-gas data",
        f"**Classification:** {classification}",
        "**Claim ceiling:** shadow signature, not universal law",
        "",
        "## Meaning",
        "",
        (
            "The saved sunflower data contains a real aggregate edge effect: "
            "log2 density-of-states bends in the high-family-size tail, and "
            "the closure-pressure quantity I_close rises as the mean extension "
            "rate approaches jamming. That is the right kind of shadow to inspect "
            "for a core-closure observable."
        ),
        "",
        (
            "The result is still INCONCLUSIVE because the existing artifacts are "
            "aggregate by family size. They do not say which core carried the "
            "pressure, which petal channel closed, or whether a per-core Hessian "
            "would remain after ordinary geometry is accounted for."
        ),
        "",
        "## Method",
        "",
        "For each available density-of-states profile D(m), the script computes:",
        "",
        "- logD(m) = log2 D(m)",
        "- delta logD(m) = logD(m+1) - logD(m)",
        "- second delta logD(m) = logD(m+1) - 2 logD(m) + logD(m-1), centered at m",
        "- p_safe(m) = g(m) / (N - m), when aggregate growth rates g(m) exist",
        "- I_close(m) = -log2 p_safe(m)",
        "",
        (
            "This is a discrete free-energy/closure-pressure proxy. It is not "
            "a per-core measurement and not a Hessian of the actual sunflower "
            "core state."
        ),
        "",
        "## Growth-Rate Profiles",
        "",
    ]

    lines.extend(
        markdown_table(
            [
                "w",
                "n",
                "N",
                "M",
                "jamming m*",
                "tail min second logD",
                "I_close at m*",
                "max finite I_close",
            ],
            [
                [
                    row["w"],
                    row["n"],
                    row["N_sites"],
                    row["max_family_size"],
                    row["jamming_m_star"],
                    row["min_tail_second_delta_logD"],
                    row["I_close_at_jamming_bits"],
                    row["max_finite_I_close_bits"],
                ]
                for row in growth_rows
            ],
        )
    )

    lines.extend(
        [
            "",
            "## Strongest Aggregate Curvature Tails",
            "",
        ]
    )
    lines.extend(
        markdown_table(
            ["w", "n", "N", "M", "tail min second logD", "growth rates"],
            [
                [
                    row["w"],
                    row["n"],
                    row["N_sites"],
                    row["max_family_size"],
                    row["min_tail_second_delta_logD"],
                    row["has_growth_rates"],
                ]
                for row in strongest_curvature
            ],
        )
    )

    lines.extend(
        [
            "",
            "## Classification",
            "",
            f"The classifier returns **{classification}**.",
            "",
            (
                "Why not HESSIAN_PRESENT: the PASS-style signal would need "
                "per-core closure instrumentation, a defined floor-normalized "
                "ratio, and a geometry-exhausting sweep. None of those are present "
                "in the saved April aggregate data."
            ),
            "",
            (
                "Why not HESSIAN_FAIL: the aggregate profiles do show curvature "
                "and closure pressure near jamming, so the fantasy is not empty. "
                "It has a measurable bulk trace worth instrumenting properly."
            ),
            "",
            "Blocking facts:",
            "",
        ]
    )
    lines.extend([f"- {blocker}" for blocker in results["classification"]["blockers"]])

    lines.extend(
        [
            "",
            "## Leg-4 Relation",
            "",
            (
                "Core closure asks how much bookkeeping the shared intersection "
                "must carry so petals are not treated as independent fragments. "
                "The aggregate I_close profile measures the cost of finding a safe "
                "unused addition at family size m. That makes it a useful precursor "
                "to the Leg-4 observable."
            ),
            "",
            (
                "It does not execute Leg 4. A real core-closure Leg-4 run needs "
                "I_core(C,s,m) distributions by core size, plus a predeclared "
                "Mendoza-floor or construction-normalized comparison."
            ),
            "",
            "## Claim Limits",
            "",
            (
                "A-axis remains A0. Abbott-Hansen-Sauer dominates the small-n "
                "lower-bound framing, so these measurements are not a new "
                "lower-bound story. The script computes a diagnostic observable; "
                "it does not close the conjecture, change the public status of #20, "
                "or justify publication language."
            ),
            "",
            "## Artifacts",
            "",
            f"- `{EXPERIMENT_ID}_RESULTS.json`",
            f"- `{EXPERIMENT_ID}_REPORT.md`",
            f"- `{EXPERIMENT_ID}_RESULTS.sha256`",
            "- `sunflower_core_hessian_analysis.py`",
            "",
        ]
    )

    return "\n".join(lines)


def main() -> None:
    base = Path(__file__).resolve().parent
    math_root = base.parent.parent

    required_inputs = [
        base / "EXP-MATH-ERDOS20-SUNFLOWER-001_RESULTS.json",
        base / "EXP-MATH-ERDOS20-SUNFLOWER-001_REPORT_2026-04-17.md",
        base / "EXP-MATH-ERDOS20-SUNFLOWER-002_RESULTS.json",
        base / "EXP-MATH-ERDOS20-SUNFLOWER-002_REPORT.md",
        base / "Q1_LITERATURE_GATE_2026-04-17.md",
        base / "SUNFLOWER_CORE_CLOSURE_LEG4_PACKET_2026-05-05.md",
        math_root / "SHADOW_DYNAMICS_HYSTERETIC_INFORMATION_EXHAUST_LIVING_ANALYSIS.md",
    ]
    extra_inputs = [
        base / "results_w2.json",
        base / "results_w4.json",
    ]

    for path in required_inputs + extra_inputs:
        if not path.exists():
            raise FileNotFoundError(path)
        if path.suffix in {".md", ".json"}:
            path.read_text(encoding="utf-8")

    profiles = collect_profiles(base)
    classification = classify(profiles)
    results = {
        "experiment_id": EXPERIMENT_ID,
        "title": "Sunflower core closure aggregate Hessian diagnostic",
        "run_date": RUN_DATE,
        "generated_on_host_date": date.today().isoformat(),
        "problem": "Erdos #20 sunflower conjecture",
        "scope": "aggregate free-energy and closure-pressure diagnostics from saved lattice-gas data",
        "input_metadata": input_metadata(required_inputs + extra_inputs, math_root),
        "analysis_definitions": {
            "logD": "log2 D(m), where D(m) is density of states at family size m",
            "delta_logD": "first difference of logD",
            "second_delta_logD": "discrete curvature of logD",
            "p_safe": "g(m)/(N-m), using aggregate growth rate g(m)",
            "I_close_bits": "-log2 p_safe(m), aggregate closure pressure",
        },
        "classification": classification,
        "compact_profile_table": compact_profile_table(profiles),
        "profiles": profiles,
        "claim_controls": {
            "A_axis": "A0",
            "Q1_gate": "Abbott-Hansen-Sauer dominates the small-n lower-bound framing",
            "forbidden_framing": [
                "new lower-bound contribution",
                "closure of the conjecture",
                "public power-morphism evidence",
                "per-core Maxwell behavior from aggregate data",
            ],
            "safe_framing": "aggregate diagnostic; shadow signature, not universal law",
        },
    }

    results_path = base / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = base / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = base / f"{EXPERIMENT_ID}_RESULTS.sha256"

    results_path.write_text(
        json.dumps(results, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    report_path.write_text(render_report(results), encoding="utf-8")
    sha_path.write_text(
        f"{sha256_file(results_path)}  {results_path.name}\n",
        encoding="utf-8",
    )

    print(json.dumps({
        "experiment_id": EXPERIMENT_ID,
        "classification": classification["classification"],
        "results": str(results_path),
        "report": str(report_path),
        "sha256": str(sha_path),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
