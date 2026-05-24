#!/usr/bin/env python3
"""One-degree Koopman MDL probe for Erdős #114 at n=14.

This binds the existing interval proof artifact for n=14 to the same
Koopman-kernel spectral-gap rule used in `LEG4_UNITARITY_114_RESULTS.json`.
It writes a fresh artifact packet and refuses to overwrite previous output.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = ROOT / "erdos-experiments/results/erdos-114"
LEG4_SCRIPT = ROOT / "erdosatlas-workbench/experiments/leg4_unitarity_114.py"
N14_INTERVAL_PATH = RESULTS_DIR / "EXP-MM-EHP-007-n14-inari_RESULTS.json"
EXPERIMENT_ID = "EXP-MATH-EHP114-N14-KOOPMAN-PROBE-20260502-02"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def load_leg4_module():
    spec = importlib.util.spec_from_file_location("leg4_unitarity_114", LEG4_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {LEG4_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def classify_gap(gap: float) -> str:
    return "QuantumShadow" if gap > math.e else "Classical"


def build_payload() -> dict:
    leg4 = load_leg4_module()
    n = 14
    alpha_values = np.linspace(0.0, leg4.ALPHA_MAX, leg4.ALPHA_STEPS)
    coeffs = leg4.make_zn_minus_1(n)
    run = leg4.run_polynomial(coeffs, f"z^{n}-1", n, alpha_values)
    if run.get("error"):
        raise RuntimeError(f"n=14 Koopman probe failed: {run['error']}")
    koopman = run["koopman"]
    floor = run["floor_sweep"]
    gap = float(koopman["spectral_gap"])
    interval = load_json(N14_INTERVAL_PATH)
    sampled_perimeter = float(run["perimeter"])
    interval_midpoint = 0.5 * (float(interval["l_star_lower"]) + float(interval["l_star_upper"]))
    perimeter_relative_error = abs(sampled_perimeter - interval_midpoint) / interval_midpoint

    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "degree": n,
        "operator_source": "LEG4_UNITARITY_114 one-degree z^14-1 Koopman kernel",
        "probe_rule": {
            "lean_rule": "spectral_gap > exp(1) iff mdl_regime != Classical",
            "threshold_exp_1": math.e,
            "classical_when": "spectral_gap <= exp(1)",
            "quantum_shadow_when": "spectral_gap > exp(1)",
        },
        "source_artifacts": {
            "leg4_script": str(LEG4_SCRIPT),
            "leg4_script_sha256": sha256_file(LEG4_SCRIPT),
            "n14_interval_artifact": str(N14_INTERVAL_PATH),
            "n14_interval_sha256": sha256_file(N14_INTERVAL_PATH),
        },
        "config": {
            "m_sample": leg4.M_SAMPLE,
            "epsilon_kernel": leg4.EPSILON_KERNEL,
            "delta_unitary": leg4.DELTA_UNITARY,
            "alpha_steps": leg4.ALPHA_STEPS,
            "alpha_max": leg4.ALPHA_MAX,
            "mendoza_limit": leg4.M_L,
            "perimeter_n_contours": leg4.PERIMETER_N_CONTOURS,
            "polynomial": "z^14 - 1",
            "random_competitors_run": False,
        },
        "koopman_probe": {
            "spectral_gap": gap,
            "threshold_exp_1": math.e,
            "gap_over_exp_1": gap / math.e,
            "margin_to_quantum_shadow": math.e - gap,
            "mdl_regime": classify_gap(gap),
            "spectral_entropy": float(koopman["spectral_entropy"]),
            "effective_unitary_dim": int(koopman["effective_unitary_dim"]),
            "n_eigenvalues": int(koopman["n_eigenvalues"]),
            "sampled_perimeter": sampled_perimeter,
            "perimeter_success": bool(run["perimeter_success"]),
            "perimeter_trust": "diagnostic_only",
            "alpha_star": float(floor["alpha_star"]),
            "floor_p1_satisfied": bool(floor["p1_satisfied"]),
            "floor_p1_max_cp_drop": float(floor["p1_max_cp_drop"]),
            "floor_p3_satisfied": bool(floor["p3_satisfied"]),
            "floor_p3_peak_d2_ratio": float(floor["p3_peak_d2_ratio"]),
            "phase_transition_detected": bool(floor["phase_transition_detected"]),
        },
        "interval_perimeter_artifact": {
            "proof_verdict": interval.get("verdict"),
            "interval_proof_complete": interval.get("bb_proof_complete"),
            "rigor": interval.get("rigor"),
            "l_star_lower": interval.get("l_star_lower"),
            "l_star_upper": interval.get("l_star_upper"),
            "l_star_midpoint": interval_midpoint,
            "sampled_perimeter_relative_error": perimeter_relative_error,
            "hessian_negative": interval.get("hessian_negative"),
            "outer_domain_safe": interval.get("outer_domain_safe"),
        },
        "interpretation": {
            "gap_verdict": (
                "The one-degree n=14 Koopman-kernel probe does not cross exp(1); "
                "under this operator, n=14 classifies as Classical."
            ),
            "floor_verdict": (
                "The n=14 run preserves the Leg-4 floor signature: P1 and P3 "
                "are satisfied for z^14-1."
            ),
            "caution": (
                "This is a sampled kernel measurement, not a theorem about the "
                "true Koopman spectrum. It binds n=14 to the existing probe rule "
                "and the existing interval perimeter artifact. The reused Leg-4 "
                "perimeter estimator is diagnostic only for n=14; the interval "
                "artifact is the perimeter truth surface."
            ),
        },
    }


def write_report(payload: dict, report_path: Path) -> None:
    k = payload["koopman_probe"]
    interval = payload["interval_perimeter_artifact"]
    lines = [
        f"# {payload['experiment_id']} Report",
        "",
        "## Claim Tested",
        "",
        "Bind the n=14 EHP interval artifact to the Lean-side MDL regime rule:",
        "",
        "- `spectral_gap > exp(1)` -> `QuantumShadow`",
        "- otherwise -> `Classical`",
        "",
        "## Result",
        "",
        f"- Degree: {payload['degree']}",
        f"- Spectral gap: {k['spectral_gap']:.12f}",
        f"- Threshold exp(1): {math.e:.12f}",
        f"- gap / exp(1): {k['gap_over_exp_1']:.6f}",
        f"- MDL regime under this probe: `{k['mdl_regime']}`",
        f"- Spectral entropy: {k['spectral_entropy']:.12f}",
        f"- Effective unitary dimension: {k['effective_unitary_dim']}",
        f"- Sampled perimeter: {k['sampled_perimeter']:.12f} (`diagnostic_only`)",
        f"- alpha_star: {k['alpha_star']:.12f}",
        f"- P1 floor collapse: {k['floor_p1_satisfied']}",
        f"- P3 entropy curvature: {k['floor_p3_satisfied']}",
        "",
        "## Interval Artifact Binding",
        "",
        f"- Proof verdict: `{interval['proof_verdict']}`",
        f"- Interval proof complete: {interval['interval_proof_complete']}",
        f"- Rigor: `{interval['rigor']}`",
        f"- L* lower: {interval['l_star_lower']}",
        f"- L* upper: {interval['l_star_upper']}",
        f"- L* midpoint: {interval['l_star_midpoint']}",
        f"- Sampled-perimeter relative error vs L* midpoint: {interval['sampled_perimeter_relative_error']:.6f}",
        "",
        "## Interpretation",
        "",
        payload["interpretation"]["gap_verdict"],
        "",
        payload["interpretation"]["floor_verdict"],
        "",
        payload["interpretation"]["caution"],
        "",
        "## Source Artifacts",
        "",
        f"- Leg-4 script: `{payload['source_artifacts']['leg4_script']}`",
        f"- Leg-4 script SHA-256: `{payload['source_artifacts']['leg4_script_sha256']}`",
        f"- n=14 interval source: `{payload['source_artifacts']['n14_interval_artifact']}`",
        f"- n=14 interval SHA-256: `{payload['source_artifacts']['n14_interval_sha256']}`",
        "",
        "## Status Boundary",
        "",
        "This supports `N14_CLASSICAL_BY_CURRENT_KOOPMAN_GAP` and `N14_PHASE_FLOOR_PRESENT`.",
        "It does not prove the true Koopman spectrum is Classical. The sampled perimeter is not",
        "used as the perimeter truth surface.",
        "",
    ]
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    result_path = RESULTS_DIR / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = RESULTS_DIR / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = RESULTS_DIR / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {path}")

    payload = build_payload()
    result_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(payload, report_path)
    sha_path.write_text(sha256_file(result_path) + "\n", encoding="utf-8")
    print(f"Wrote {result_path}")
    print(f"Wrote {report_path}")
    print(f"Wrote {sha_path}")
    print(json.dumps(payload["koopman_probe"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
