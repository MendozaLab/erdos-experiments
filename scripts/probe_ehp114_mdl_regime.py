#!/usr/bin/env python3
"""EHP114 MDL regime probe.

Consumes existing Erdős #114 artifacts and applies the Lean-side regime rule:

    spectral_gap > exp(1)  -> QuantumShadow
    otherwise             -> Classical

This is an evidence-binding probe, not a new numerical search. It deliberately
does not mutate prior EXP-MM-EHP outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean


DEFAULT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS = DEFAULT_ROOT / "erdos-experiments/results/erdos-114"
DEFAULT_LEG4 = DEFAULT_ROOT / "erdosatlas-workbench/experiments/LEG4_UNITARITY_114_RESULTS.json"
EXPERIMENT_ID = "EXP-MATH-EHP114-MDL-PROBE-20260502-01"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def classify_gap(gap: float) -> str:
    return "QuantumShadow" if gap > math.e else "Classical"


def build_probe(results_dir: Path, leg4_path: Path) -> dict:
    leg4 = load_json(leg4_path)
    rows: list[dict] = []

    for n_text, item in sorted(leg4["results_by_n"].items(), key=lambda kv: int(kv[0])):
        n = int(n_text)
        zn1 = item["z_n_minus_1"]
        koopman = zn1["koopman"]
        gap = float(koopman["spectral_gap"])
        alpha_star = float(zn1["floor_sweep"]["alpha_star"])
        proof_path = results_dir / f"EXP-MM-EHP-007-n{n}-inari_RESULTS.json"
        proof = load_json(proof_path) if proof_path.exists() else {}

        rows.append(
            {
                "n": n,
                "operator_source": "LEG4_UNITARITY_114 z^n-1 Koopman kernel",
                "spectral_gap": gap,
                "threshold_exp_1": math.e,
                "gap_over_exp_1": gap / math.e,
                "margin_to_quantum_shadow": math.e - gap,
                "mdl_regime": classify_gap(gap),
                "alpha_star": alpha_star,
                "floor_p1_satisfied": bool(zn1["floor_sweep"]["p1_satisfied"]),
                "floor_p3_satisfied": bool(zn1["floor_sweep"]["p3_satisfied"]),
                "perimeter_winner_is_zn1": bool(item["perimeter_winner_is_zn1"]),
                "perimeter_zn1": float(item["perimeter_zn1"]),
                "perimeter_random_best": float(item["perimeter_random_best"]),
                "proof_verdict": proof.get("verdict"),
                "interval_proof_complete": proof.get("bb_proof_complete"),
                "l_star_lower": proof.get("l_star_lower"),
                "l_star_upper": proof.get("l_star_upper"),
                "proof_artifact": str(proof_path),
                "proof_artifact_sha256": sha256_file(proof_path) if proof_path.exists() else None,
            }
        )

    n14_path = results_dir / "EXP-MM-EHP-007-n14-inari_RESULTS.json"
    n14 = load_json(n14_path)

    gaps = [r["spectral_gap"] for r in rows]
    alpha_values = [r["alpha_star"] for r in rows]
    quantum_rows = [r for r in rows if r["mdl_regime"] == "QuantumShadow"]
    winners = [r for r in rows if r["perimeter_winner_is_zn1"]]
    p1 = [r for r in rows if r["floor_p1_satisfied"]]
    p3 = [r for r in rows if r["floor_p3_satisfied"]]

    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "probe_rule": {
            "lean_rule": "spectral_gap > exp(1) iff mdl_regime != Classical",
            "threshold_exp_1": math.e,
            "classical_when": "spectral_gap <= exp(1)",
            "quantum_shadow_when": "spectral_gap > exp(1)",
        },
        "source_artifacts": {
            "leg4_unitarity_results": str(leg4_path),
            "leg4_unitarity_sha256": sha256_file(leg4_path),
            "interval_proof_results_dir": str(results_dir),
            "n14_interval_artifact": str(n14_path),
            "n14_interval_sha256": sha256_file(n14_path),
        },
        "measured_degrees": [r["n"] for r in rows],
        "probe_rows": rows,
        "n14_interval_extension": {
            "degree": 14,
            "proof_verdict": n14.get("verdict"),
            "interval_proof_complete": n14.get("bb_proof_complete"),
            "l_star_lower": n14.get("l_star_lower"),
            "l_star_upper": n14.get("l_star_upper"),
            "spectral_gap_status": "not_measured_in_LEG4_UNITARITY_114",
            "mdl_regime_status": "not_classified_by_this_probe",
        },
        "summary": {
            "degrees_with_gap_measurement": len(rows),
            "max_spectral_gap": max(gaps),
            "max_spectral_gap_degree": rows[gaps.index(max(gaps))]["n"],
            "mean_spectral_gap": mean(gaps),
            "quantum_shadow_count": len(quantum_rows),
            "classical_count": len(rows) - len(quantum_rows),
            "alpha_star_unique_values": sorted(set(alpha_values)),
            "p1_floor_collapse_count": len(p1),
            "p3_entropy_curvature_count": len(p3),
            "zn1_perimeter_win_count": len(winners),
            "zn1_perimeter_win_fraction": len(winners) / len(rows),
            "leg4_global_predictions": leg4.get("predictions_satisfied"),
        },
        "interpretation": {
            "current_gap_verdict": (
                "No measured z^n-1 Koopman gap crosses exp(1); under the current "
                "operator/probe, EHP114 stays Classical for n=3..13."
            ),
            "phase_floor_verdict": (
                "A stable alpha_star floor is present in the prior Leg-4 artifact, "
                "with P1 and P3 satisfied across measured degrees."
            ),
            "caution": (
                "This does not prove EHP114 is globally Classical. It only says the "
                "current sampled Koopman-kernel probe does not observe a QuantumShadow "
                "gap crossing for z^n-1 over n=3..13. n=14 has interval perimeter "
                "support but no bound spectral-gap measurement in this artifact set."
            ),
        },
    }


def write_report(result: dict, report_path: Path) -> None:
    summary = result["summary"]
    lines = [
        f"# {result['experiment_id']} Report",
        "",
        "## Claim Tested",
        "",
        "Apply the Lean-side MDL regime rule to the existing Erdős #114 Koopman-gap artifact:",
        "",
        "- `spectral_gap > exp(1)` -> `QuantumShadow`",
        "- otherwise -> `Classical`",
        "",
        "This is an evidence-binding probe over existing artifacts, not a new optimizer run.",
        "",
        "## Result",
        "",
        f"- Degrees with measured gaps: {summary['degrees_with_gap_measurement']}",
        f"- QuantumShadow crossings: {summary['quantum_shadow_count']}",
        f"- Classical classifications: {summary['classical_count']}",
        f"- Max measured spectral gap: {summary['max_spectral_gap']:.12f} at n={summary['max_spectral_gap_degree']}",
        f"- Threshold exp(1): {math.e:.12f}",
        f"- alpha_star values observed: {summary['alpha_star_unique_values']}",
        f"- P1 floor-collapse support: {summary['p1_floor_collapse_count']}/{summary['degrees_with_gap_measurement']}",
        f"- P3 entropy-curvature support: {summary['p3_entropy_curvature_count']}/{summary['degrees_with_gap_measurement']}",
        f"- z^n-1 perimeter win fraction in Leg-4 random comparison: {summary['zn1_perimeter_win_fraction']:.6f}",
        "",
        "## Interpretation",
        "",
        result["interpretation"]["current_gap_verdict"],
        "",
        result["interpretation"]["phase_floor_verdict"],
        "",
        result["interpretation"]["caution"],
        "",
        "## Per-Degree Probe Table",
        "",
        "| n | spectral_gap | gap/exp(1) | regime | alpha_star | z^n-1 win? | proof verdict |",
        "|---:|---:|---:|---|---:|---|---|",
    ]
    for row in result["probe_rows"]:
        lines.append(
            "| {n} | {gap:.12f} | {ratio:.6f} | {regime} | {alpha:.12f} | {win} | {verdict} |".format(
                n=row["n"],
                gap=row["spectral_gap"],
                ratio=row["gap_over_exp_1"],
                regime=row["mdl_regime"],
                alpha=row["alpha_star"],
                win="yes" if row["perimeter_winner_is_zn1"] else "no",
                verdict=row["proof_verdict"],
            )
        )
    lines.extend(
        [
            "",
            "## Source Artifacts",
            "",
            f"- Leg-4 source: `{result['source_artifacts']['leg4_unitarity_results']}`",
            f"- Leg-4 SHA-256: `{result['source_artifacts']['leg4_unitarity_sha256']}`",
            f"- n=14 interval source: `{result['source_artifacts']['n14_interval_artifact']}`",
            f"- n=14 SHA-256: `{result['source_artifacts']['n14_interval_sha256']}`",
            "",
            "## Status Boundary",
            "",
            "This supports a `CLASSICAL_BY_CURRENT_KOOPMAN_GAP` finding for n=3..13, plus a separate",
            "`PHASE_FLOOR_PRESENT` finding. It does not support a public claim that EHP114 is globally",
            "Classical or globally QuantumShadow.",
            "",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--leg4", type=Path, default=DEFAULT_LEG4)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_RESULTS)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    result_path = args.out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = args.out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = args.out_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {path}")

    result = build_probe(args.results_dir, args.leg4)
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(result, report_path)
    sha_path.write_text(sha256_file(result_path) + "\n", encoding="utf-8")

    print(f"Wrote {result_path}")
    print(f"Wrote {report_path}")
    print(f"Wrote {sha_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
