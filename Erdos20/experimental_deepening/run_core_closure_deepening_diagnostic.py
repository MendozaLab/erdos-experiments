#!/usr/bin/env python3
"""Secondary diagnostic for Erdos #20 core-closure artifacts.

This script reads only existing #20 result artifacts. It does not enumerate new
families, update D1, touch scorecards, or write outside experimental_deepening.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean, pstdev
from typing import Any


EXPERIMENT_ID = "EXP-MATH-ERDOS20-CORE-CLOSURE-DEEPENING-20260505-01"
ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = Path(__file__).resolve().parent

INPUTS = [
    ROOT / "SUNFLOWER_CORE_CLOSURE_LEG4_PACKET_2026-05-05.md",
    ROOT / "EXP-MATH-ERDOS20-HESSIAN-CLOSURE-20260505-01_RESULTS.json",
    ROOT / "EXP-MATH-ERDOS20-PER-CORE-CLOSURE-20260505-01_RESULTS.json",
    ROOT / "EXP-MATH-ERDOS20-RUST-PER-CORE-CLOSURE-20260505-01_RESULTS.json",
    ROOT / "EXP-MATH-ERDOS20-RUST-PER-CORE-CLOSURE-W3N8-20260505-02_RESULTS.json",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def round_float(value: float | None, places: int = 6) -> float | None:
    if value is None or not math.isfinite(value):
        return None
    return round(value, places)


def linreg_slope(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    x_bar = mean(xs)
    y_bar = mean(ys)
    den = sum((x - x_bar) ** 2 for x in xs)
    if den == 0:
        return None
    return sum((x - x_bar) * (y - y_bar) for x, y in zip(xs, ys)) / den


def cv(values: list[float]) -> float | None:
    positive = [value for value in values if value is not None and math.isfinite(value)]
    if len(positive) < 2:
        return None
    avg = mean(positive)
    if avg == 0:
        return None
    return pstdev(positive) / avg


def trend_label(slope: float | None) -> str:
    if slope is None:
        return "INSUFFICIENT_POINTS"
    abs_slope = abs(slope)
    if abs_slope <= 0.15:
        return "STABLE_BY_SLOPE_ONLY"
    if abs_slope > 0.25:
        return "GEOMETRY_DRIFT"
    return "BORDERLINE_DRIFT"


def build_python_observations(per_core: dict[str, Any]) -> list[dict[str, Any]]:
    observations: list[dict[str, Any]] = []
    for run in per_core.get("runs", []):
        summary = run.get("summary", {})
        target_m = summary.get("target_m")
        rows = []
        for row in run.get("per_core_by_size", []):
            if row.get("m") == target_m:
                rows.append(
                    {
                        "core_size_s": row.get("core_size_s"),
                        "petal_size": row.get("petal_size"),
                        "I_core_local_bits": row.get("I_core_local_weighted_bits"),
                        "candidate_count": row.get("candidate_count"),
                        "local_valid_count": row.get("local_valid_count"),
                    }
                )
        observations.append(
            {
                "source": "python_exact",
                "w": run.get("w"),
                "n": run.get("n"),
                "N_sites": run.get("N_sites"),
                "total_sf_free": run.get("total_sf_free"),
                "target_m": target_m,
                "jamming_m_star": run.get("jamming_m_star"),
                "aggregate_I_close_bits": summary.get("aggregate_I_close_bits_at_target"),
                "all_unused_I_close_bits": summary.get("all_unused_I_close_bits_at_target"),
                "strongest_I_core_local_bits": summary.get("strongest_I_core_local_weighted_bits_at_target"),
                "strongest_core_size_s": summary.get("strongest_core_size_at_target"),
                "per_core_rows": rows,
            }
        )
    return observations


def build_rust_observations(*rust_results: dict[str, Any]) -> list[dict[str, Any]]:
    observations: list[dict[str, Any]] = []
    for result in rust_results:
        for target in result.get("targets", []):
            if target.get("status") != "PASS":
                continue
            data = target.get("data", {})
            rows = []
            for row in data.get("core_rows", []):
                rows.append(
                    {
                        "core_size_s": row.get("core_size_s"),
                        "petal_size": row.get("petal_size"),
                        "I_core_local_bits": row.get("I_core_local_bits"),
                        "candidate_count": row.get("candidate_count"),
                        "local_valid_count": row.get("local_valid_count"),
                    }
                )
            strongest = max(
                rows,
                key=lambda row: row.get("I_core_local_bits")
                if row.get("I_core_local_bits") is not None
                else -1.0,
                default={},
            )
            observations.append(
                {
                    "source": "rust_exact",
                    "w": data.get("w"),
                    "n": data.get("n"),
                    "N_sites": data.get("site_count"),
                    "total_sf_free": data.get("total_sf_free"),
                    "target_m": data.get("target_m"),
                    "jamming_m_star": data.get("auto_jamming_m"),
                    "aggregate_I_close_bits": data.get("aggregate_I_close_bits_at_target"),
                    "all_unused_I_close_bits": None,
                    "strongest_I_core_local_bits": strongest.get("I_core_local_bits"),
                    "strongest_core_size_s": strongest.get("core_size_s"),
                    "per_core_rows": rows,
                }
            )
    return observations


def dedupe_observations(observations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key: dict[tuple[int, int], dict[str, Any]] = {}
    source_rank = {"python_exact": 1, "rust_exact": 2}
    for obs in observations:
        key = (obs["w"], obs["n"])
        existing = by_key.get(key)
        if existing is None or source_rank[obs["source"]] >= source_rank[existing["source"]]:
            by_key[key] = obs
    return sorted(by_key.values(), key=lambda item: (item["w"], item["n"]))


def series_diagnostics(observations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    diagnostics: list[dict[str, Any]] = []
    by_w: dict[int, list[dict[str, Any]]] = defaultdict(list)
    by_ws: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for obs in observations:
        if obs.get("strongest_I_core_local_bits") is not None:
            by_w[obs["w"]].append(obs)
        for row in obs.get("per_core_rows", []):
            value = row.get("I_core_local_bits")
            if value is not None and value > 0:
                by_ws[(obs["w"], row["core_size_s"])].append(
                    {
                        "n": obs["n"],
                        "N_sites": obs["N_sites"],
                        "I_core_local_bits": value,
                        "source": obs["source"],
                    }
                )

    for w, rows in sorted(by_w.items()):
        usable = [row for row in rows if row["strongest_I_core_local_bits"] and row["strongest_I_core_local_bits"] > 0]
        if len(usable) < 3:
            diagnostics.append(
                {
                    "series": f"w={w} strongest per-core",
                    "point_count": len(usable),
                    "classification": "INSUFFICIENT_POINTS",
                    "meaning": "Fewer than three nonzero jamming observations after dedupe.",
                }
            )
            continue
        xs = [math.log(row["N_sites"]) for row in usable]
        ys = [math.log(row["strongest_I_core_local_bits"]) for row in usable]
        slope = linreg_slope(xs, ys)
        late_values = [row["strongest_I_core_local_bits"] for row in usable[-3:]]
        diagnostics.append(
            {
                "series": f"w={w} strongest per-core",
                "point_count": len(usable),
                "n_window": [row["n"] for row in usable],
                "N_window": [row["N_sites"] for row in usable],
                "values_bits": [round_float(row["strongest_I_core_local_bits"]) for row in usable],
                "loglog_slope_vs_N": round_float(slope),
                "late3_cv": round_float(cv(late_values)),
                "classification": trend_label(slope),
                "meaning": "Slope is a geometry-drift screen, not a floor-normalized Leg-4 verdict.",
            }
        )

    for (w, s), rows in sorted(by_ws.items()):
        rows = sorted(rows, key=lambda row: row["n"])
        if len(rows) < 3:
            continue
        xs = [math.log(row["N_sites"]) for row in rows]
        ys = [math.log(row["I_core_local_bits"]) for row in rows]
        slope = linreg_slope(xs, ys)
        diagnostics.append(
            {
                "series": f"w={w}, fixed core size s={s}",
                "point_count": len(rows),
                "n_window": [row["n"] for row in rows],
                "N_window": [row["N_sites"] for row in rows],
                "values_bits": [round_float(row["I_core_local_bits"]) for row in rows],
                "loglog_slope_vs_N": round_float(slope),
                "late3_cv": round_float(cv([row["I_core_local_bits"] for row in rows[-3:]])),
                "classification": trend_label(slope),
                "meaning": "Fixed-s series is closer to the needed Leg-4 observable than strongest-core switching.",
            }
        )

    return diagnostics


def write_report(result: dict[str, Any], report_path: Path) -> None:
    lines = [
        f"# {EXPERIMENT_ID} Report",
        "",
        "## Status",
        "",
        "Secondary deterministic diagnostic from saved Erdos #20 artifacts only. No new family enumeration, no D1, no scorecard, no public docs.",
        "",
        "## Verdict",
        "",
        f"- Classification: `{result['classification']}`",
        f"- Claim ceiling: `{result['claim_ceiling']}`",
        "- Meaning: existing data support a measurable per-core channel, but not a Leg-4 pass or theorem progress.",
        "",
        "## Drift Screen",
        "",
        "| Series | points | slope vs N | late CV | class |",
        "|---|---:|---:|---:|---|",
    ]
    for diag in result["diagnostics"]:
        lines.append(
            "| {series} | {point_count} | {slope} | {cv} | {classification} |".format(
                series=diag["series"],
                point_count=diag["point_count"],
                slope=diag.get("loglog_slope_vs_N"),
                cv=diag.get("late3_cv"),
                classification=diag["classification"],
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The w=3 fixed-core-size series is the useful next lane because it now spans n=5..8 after the Rust W3N8 artifact. The screen is still only geometry-drift analysis: it does not define a Mendoza-floor numerator, compare Abbott-Hansen-Sauer constructions, or show a theorem beyond encoding.",
            "",
            "## Next Measurable",
            "",
            "Pre-register fixed core-size local closure cost `I_core_local(s,m*)` with a defect/hysteresis window around `m*`, then run a symmetry-reduced or Rust exact sweep where feasible. The immediate target is not a sunflower bound; it is whether the per-core channel stays stable once ordinary ambient geometry is accounted for.",
            "",
            "## Source Boundary",
            "",
            "Inputs are listed with SHA-256 hashes in the companion results JSON.",
        ]
    )
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    result_path = OUT_DIR / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = OUT_DIR / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = OUT_DIR / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {path}")

    hessian = load_json(INPUTS[1])
    per_core = load_json(INPUTS[2])
    rust_01 = load_json(INPUTS[3])
    rust_w3n8 = load_json(INPUTS[4])

    observations = dedupe_observations(
        build_python_observations(per_core) + build_rust_observations(rust_01, rust_w3n8)
    )
    diagnostics = series_diagnostics(observations)
    drift_count = sum(1 for diag in diagnostics if diag["classification"] == "GEOMETRY_DRIFT")
    stable_count = sum(1 for diag in diagnostics if diag["classification"] == "STABLE_BY_SLOPE_ONLY")
    classification = (
        "DEEPENING_GEOMETRY_DRIFT_PRESENT"
        if drift_count
        else "DEEPENING_STABLE_PRECURSOR_PRESENT"
        if stable_count
        else "DEEPENING_INCONCLUSIVE"
    )

    result = {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #20 sunflower core closure",
        "scope": "secondary drift/CV screen over saved aggregate, per-core, and Rust exact artifacts",
        "input_metadata": [
            {
                "path": str(path.relative_to(ROOT.parent.parent)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in INPUTS
        ],
        "classification": classification,
        "claim_ceiling": "A-axis A0; shadow signature, not universal law; no theorem, no lower-bound progress, no Leg-4 pass",
        "source_classifications": {
            "hessian": hessian.get("classification", {}).get("classification"),
            "per_core": per_core.get("classification", {}).get("classification"),
            "rust_01": rust_01.get("classification"),
            "rust_w3n8": rust_w3n8.get("classification"),
        },
        "observations": observations,
        "diagnostics": diagnostics,
        "next_leg4_measurable": {
            "quantity": "fixed-core-size local closure cost with defect/hysteresis window",
            "symbol": "I_core_local(s,m) near m* and across m* +/- delta",
            "why": "It keeps the channel tied to a specified core size instead of letting the strongest core switch with n.",
            "not_yet": [
                "no Mendoza-floor numerator",
                "no Abbott-Hansen-Sauer construction baseline",
                "no three-n geometry-exhausting sweep beyond w=3",
                "no theorem beyond encoding",
            ],
        },
        "guardrails": {
            "read_existing_artifacts_only": True,
            "no_d1": True,
            "no_scorecard": True,
            "no_public_docs": True,
            "no_git": True,
            "no_downloads": True,
        },
    }

    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(result, report_path)
    sha_path.write_text(f"{sha256_file(result_path)}  {result_path.name}\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "experiment_id": EXPERIMENT_ID,
                "classification": classification,
                "result": str(result_path),
                "report": str(report_path),
                "sha256": sha_path.read_text(encoding="utf-8").split()[0],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
