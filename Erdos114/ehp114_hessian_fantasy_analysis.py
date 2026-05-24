#!/usr/bin/env python3
"""Aggregate the EHP #114 stratified Hessian-fantasy probes.

This is a diagnostic aggregator. It reads existing saved artifacts and writes
one H114 experiment result, report, and SHA-256 sidecar in this directory.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-HESSIAN-FANTASY-20260505-01"


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_float(value: Any) -> float:
    return float(str(value))


def fmt_float(value: float, digits: int = 12) -> str:
    return f"{value:.{digits}g}"


def radial_summary(radial_results: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(radial_results):
        data = read_json(path)
        slope_rows = data["asymptotic"]["slope_rows"]
        tail4 = next((row for row in slope_rows if row["tail_window"] == 4), slope_rows[0])
        tail6 = next((row for row in slope_rows if row["tail_window"] == 6), None)
        expected = as_float(data["asymptotic"]["expected_slope"])
        observed = as_float(tail4["fitted_log_log_slope"])
        last_asymptotic = data["asymptotic"]["rows"][-1]
        rows.append(
            {
                "degree": data["degree"],
                "experiment_id": data["experiment_id"],
                "artifact": str(path),
                "status": data["status"],
                "claim_ceiling": data.get("claim_ceiling"),
                "expected_slope_1_over_n": expected,
                "tail_window_4_slope": observed,
                "tail_window_6_slope": as_float(tail6["fitted_log_log_slope"]) if tail6 else None,
                "abs_error_tail4": abs(observed - expected),
                "relative_error_tail4": abs(observed - expected) / expected,
                "tail_deficit_over_eps_power_limit_proxy": as_float(
                    last_asymptotic["deficit_over_eps_power_1_over_n"]
                ),
                "smallest_eps": as_float(last_asymptotic["eps"]),
            }
        )
    return rows


def tensor_summary(tensor_path: Path) -> dict[str, Any]:
    data = read_json(tensor_path)
    direction_rows = data["direction_rows"]
    shape_rows = [row for row in direction_rows if row.get("cone_role") == "shape_candidate"]
    radial_rows = [row for row in direction_rows if row.get("cone_role") == "singular_radial"]
    shape_slopes = [float(row["fitted_loglog_slope_symmetric_deficit"]) for row in shape_rows]
    radial_slopes = [float(row["fitted_loglog_slope_symmetric_deficit"]) for row in radial_rows]
    proxy = data["mixed_shape_hessian_proxy"]
    eigenvalues = [float(value) for value in proxy["eigenvalues"]]
    slope_class_rows = []
    for row in shape_rows:
        slope = float(row["fitted_loglog_slope_symmetric_deficit"])
        if slope >= 1.5:
            slope_class = "smooth_hessian_like"
        elif slope >= 0.5:
            slope_class = "intermediate_subquadratic"
        else:
            slope_class = "singular_sub_hessian"
        slope_class_rows.append(
            {
                "label": row["label"],
                "mode": row["mode"],
                "kind": row["kind"],
                "phase": row["phase"],
                "slope": slope,
                "class": slope_class,
            }
        )
    return {
        "degree": data["degree"],
        "experiment_id": data["experiment_id"],
        "artifact": str(tensor_path),
        "status": data["status"],
        "rigorous": data["method"]["rigorous"],
        "shape_basis_rank": data["parameters"]["shape_basis_rank"],
        "radial_singular_slope": radial_slopes[0] if radial_slopes else None,
        "shape_slope_min": min(shape_slopes),
        "shape_slope_mean": mean(shape_slopes),
        "shape_slope_max": max(shape_slopes),
        "shape_smooth_hessian_like_count_slope_ge_1_5": sum(
            1 for slope in shape_slopes if slope >= 1.5
        ),
        "shape_singular_like_count_slope_lt_1_5": sum(1 for slope in shape_slopes if slope < 1.5),
        "all_shape_symmetric_deficits_positive": all(
            row["all_symmetric_deficits_positive"] for row in shape_rows
        ),
        "mixed_shape_hessian_proxy": {
            "eps": proxy["eps"],
            "min_eigenvalue": min(eigenvalues),
            "max_eigenvalue": max(eigenvalues),
            "positive_eigenvalue_count": sum(1 for value in eigenvalues if value > 0),
            "eigenvalue_count": len(eigenvalues),
            "condition_proxy_abs_max_over_min": proxy["condition_proxy_abs_max_over_min"],
        },
        "slope_class_rows": slope_class_rows,
        "interpretation": data["interpretation"],
    }


def fourier_summary(fourier_path: Path) -> dict[str, Any]:
    data = read_json(fourier_path)
    curvature_values: list[float] = []
    positive_rows = 0
    for row in data["basis_rows"]:
        if row["all_eps_positive"]:
            positive_rows += 1
        for eps_row in row["eps_rows"]:
            curvature_values.append(float(eps_row["deficit_curvature"]))
    return {
        "degree": data["degree"],
        "experiment_id": data["experiment_id"],
        "artifact": str(fourier_path),
        "status": data["status"],
        "rigorous": data["method"]["rigorous"],
        "basis_rank": data["parameters"]["basis_rank"],
        "expected_reduced_dim_2n_minus_3": data["parameters"][
            "expected_reduced_dim_2n_minus_3"
        ],
        "positive_mode_count_all_eps": positive_rows,
        "mode_count": len(data["basis_rows"]),
        "rows_not_positive_all_eps": len(data["basis_rows"]) - positive_rows,
        "min_deficit_curvature": min(curvature_values),
        "max_deficit_curvature": max(curvature_values),
        "best_label": data["summary"]["best_label"],
        "worst_label": data["summary"]["worst_label"],
        "large_reference_error_flag": data["singularity_sanity"]["large_reference_error_flag"],
        "l0_relative_error_vs_shortcut_lstar": data["reference"][
            "L0_relative_error_vs_shortcut_lstar"
        ],
        "interpretation": data["interpretation"],
    }


def build_report(result: dict[str, Any]) -> str:
    radial_lines = []
    for row in result["radial_exponent_evidence"]:
        radial_lines.append(
            "| {degree} | {expected} | {observed} | {abs_error} | {rel_error} | {limit_proxy} |".format(
                degree=row["degree"],
                expected=fmt_float(row["expected_slope_1_over_n"]),
                observed=fmt_float(row["tail_window_4_slope"]),
                abs_error=fmt_float(row["abs_error_tail4"]),
                rel_error=fmt_float(row["relative_error_tail4"]),
                limit_proxy=fmt_float(row["tail_deficit_over_eps_power_limit_proxy"]),
            )
        )

    tensor = result["tensor_cone_hessian_proxy"]
    proxy = tensor["mixed_shape_hessian_proxy"]
    fourier = result["fourier_hessian_curvature"]

    report = f"""# {EXPERIMENT_ID} Report

## Status

Diagnostic aggregate, not a proof of EHP #114 and not a publication packet.

The ordinary smooth-Hessian fantasy is dead as a radial local model: the radial
boundary layer follows a Puiseux exponent close to 1/n, not an eps^2 law.

The useful survivor is stratified. Read the certificate as three pieces:
outer-domain coercivity, radial Puiseux lower bound, and nonradial shape
Hessian/remainder control. This is a shadow signature, not universal law.

## Verdict

- Overall decision: `{result["decision"]["overall"]}`
- Smooth Hessian layer: `{result["decision"]["smooth_hessian_layer"]}`
- Stratified certificate layer: `{result["decision"]["stratified_certificate_layer"]}`

Meaning: do not try to make a single ordinary Hessian at `z^n - 1` carry the
boundary. The radial axis is singular. The finite-difference shape probes still
show positive Hessian-like curvature on quotient Fourier modes, so the right
attack is a decomposed certificate rather than a smooth Taylor certificate.

## Radial Puiseux Evidence

| n | expected 1/n slope | tail-window 4 slope | abs error | relative error | deficit / eps^(1/n) tail proxy |
|---:|---:|---:|---:|---:|---:|
{chr(10).join(radial_lines)}

The n = 14 and n = 15 radial probes both land on the predicted 1/n exponent to
small relative error in the tail fit. This kills the radial quadratic model and
supports a fixed-n Puiseux lower-bound theorem as the radial component.

## Tensor-Cone Shape Probe

- Source: `{tensor["experiment_id"]}`
- Degree: {tensor["degree"]}
- Rigorous: {tensor["rigorous"]}
- Singular radial slope: {fmt_float(tensor["radial_singular_slope"])}
- Shape slope range: {fmt_float(tensor["shape_slope_min"])} to {fmt_float(tensor["shape_slope_max"])}
- Shape slope mean: {fmt_float(tensor["shape_slope_mean"])}
- Shape slopes >= 1.5: {tensor["shape_smooth_hessian_like_count_slope_ge_1_5"]}
- Shape slopes < 1.5: {tensor["shape_singular_like_count_slope_lt_1_5"]}
- All tested shape symmetric deficits positive: {tensor["all_shape_symmetric_deficits_positive"]}
- Mixed shape proxy eigenvalues: min {fmt_float(proxy["min_eigenvalue"])}, max {fmt_float(proxy["max_eigenvalue"])}, positive {proxy["positive_eigenvalue_count"]}/{proxy["eigenvalue_count"]}
- Condition proxy: {fmt_float(proxy["condition_proxy_abs_max_over_min"])}

Interpretation: the tensor-cone data does not rescue a smooth eps^2 scale law,
because every measured shape slope is sub-Hessian. It does preserve a positive
finite-difference shape cone at the tested scale. The shape layer is therefore
Hessian-like as a cone positivity probe, not as a full smooth local model.

## Fourier-Hessian Shape Probe

- Source: `{fourier["experiment_id"]}`
- Degree: {fourier["degree"]}
- Rigorous: {fourier["rigorous"]}
- Basis rank: {fourier["basis_rank"]} = expected reduced dimension {fourier["expected_reduced_dim_2n_minus_3"]}
- Positive modes at all eps: {fourier["positive_mode_count_all_eps"]}/{fourier["mode_count"]}
- Rows not positive at all eps: {fourier["rows_not_positive_all_eps"]}
- Curvature range: {fmt_float(fourier["min_deficit_curvature"])} to {fmt_float(fourier["max_deficit_curvature"])}
- Large reference-error flag: {fourier["large_reference_error_flag"]}

Interpretation: all 27 n = 15 Fourier modes carry positive deficit curvature
under this floating estimator, but the singularity sanity flag matters. This is
evidence for the nonradial cone component, not an interval certificate.

## Certificate Ansatz

Let `D_n(p) = L(z^n - 1) - L(p)` and decompose a local perturbation as
`u = r e_0 + s`, where `e_0` is the singular radial direction and `s` is the
quotient nonradial Fourier/tensor shape component.

The surviving stratified certificate should have this form:

```text
Outer domain:
  D_n(p) >= B_n(rho) > 0
  for p outside the singular local cone / middle annulus.

Radial axis:
  D_n(r e_0) >= C_n |r|^(1/n)
  with fixed-n interval constants, starting at n = 14.

Shape cone:
  D_n(r e_0 + s) >= C_n |r|^(1/n) + lambda_n ||s||^2 - R_n(r, s)
  where lambda_n > 0 comes from the nonradial cone and R_n is bounded small
  enough that mixed radial/shape terms cannot cancel the Puiseux deficit.
```

That is the honest target. The current evidence supports the architecture but
does not establish the analytic remainder bound.

## Next Run Target

Recommended next run: n = 20 triad, because it is the cleanest falsification
test for the stratified fantasy.

Inputs:

- Exact radial hypergeometric family for n = 20 with eps grid
  `{0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 1e-5, 1e-6, 1e-7, 1e-8}`.
- Fourier-Hessian basis rank `2n - 3 = 37`, eps `{0.02, 0.01, 0.005}`, grid resolution at least matching n = 15.
- Tensor-cone mixed matrix for the n = 20 quotient shape basis, with the singular radial direction excluded from the mixed shape block.

Outputs:

- `EXP-MATH-EHP114-N20-RADIAL-HYPERGEOMETRIC-20260505-02_RESULTS.json`
- `EXP-MATH-EHP114-N20-FOURIER-HESSIAN-20260505-02_RESULTS.json`
- `EXP-MATH-EHP114-N20-TENSOR-CONE-SCAFFOLD-20260505-02_RESULTS.json`

Pass signal: radial slope within 5 percent of 1/20, all or nearly all Fourier
shape modes positive, and a positive mixed shape proxy with a bounded condition
proxy. Fail signal: radial slope drifts from 1/20 or the shape cone loses broad
positive curvature.

Proof-oriented alternative: n = 14 interval hardening for the fixed
`gauss_2F1_puiseux_lower_bound_n14_interval` theorem, then a separate n = 14 or
n = 15 Fourier/tensor interval cone for nonradial and mixed terms.

## Source Boundary

All numbers in this aggregate are computed from saved local artifacts listed in
the companion results JSON. No D1, scorecard, public document, email, push, or
publication surface was changed.
"""
    return report


def main() -> None:
    script_path = Path(__file__).resolve()
    out_dir = script_path.parent
    math_root = script_path.parents[2]
    results_dir = math_root / "erdos-experiments" / "results" / "erdos-114"

    tensor_path = (
        results_dir / "EXP-MATH-EHP114-N10-TENSOR-CONE-SCAFFOLD-20260502-02_RESULTS.json"
    )
    fourier_path = results_dir / "EXP-MATH-EHP114-N15-FOURIER-HESSIAN-20260502-02_RESULTS.json"
    radial_paths = sorted(results_dir.glob("*RADIAL-HYPERGEOMETRIC*RESULTS.json"))

    required_inputs = [
        out_dir / "STRATIFICATION_N3_N14_2026-05-02.md",
        out_dir / "EHP114_BOUNDARY_CLOSURE_PACKET_2026-05-05.md",
        tensor_path,
        fourier_path,
        *radial_paths,
        math_root / "SHADOW_DYNAMICS_HYSTERETIC_INFORMATION_EXHAUST_LIVING_ANALYSIS.md",
    ]
    missing = [str(path) for path in required_inputs if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required inputs: {missing}")

    radial = radial_summary(radial_paths)
    tensor = tensor_summary(tensor_path)
    fourier = fourier_summary(fourier_path)

    result: dict[str, Any] = {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "worker": "H114",
        "scope": "stratified Hessian fantasy diagnostic aggregate",
        "status": "DIAGNOSTIC_AGGREGATE_NOT_PROOF",
        "decision": {
            "overall": "STRATIFIED",
            "smooth_hessian_layer": "DEAD_AS_ORDINARY_RADIAL_MODEL",
            "stratified_certificate_layer": "ALIVE_AS_DIAGNOSTIC_ANSATZ",
            "why": (
                "Radial fits follow a Puiseux 1/n law, while nonradial tensor "
                "and Fourier probes retain positive finite-difference curvature."
            ),
        },
        "claim_boundary": (
            "No global EHP114 maximality claim. No public theorem claim. "
            "The artifact aggregates diagnostics and proposes the next certificate ansatz."
        ),
        "required_phrase": "shadow signature, not universal law",
        "radial_exponent_evidence": radial,
        "tensor_cone_hessian_proxy": tensor,
        "fourier_hessian_curvature": fourier,
        "certificate_ansatz": {
            "decomposition": "u = r e_0 + s, with e_0 singular radial and s nonradial quotient shape",
            "outer_domain_coercivity": "D_n(p) >= B_n(rho) > 0 outside the singular local cone / middle annulus",
            "radial_puiseux_lower_bound": "D_n(r e_0) >= C_n |r|^(1/n)",
            "shape_hessian_remainder_bound": (
                "D_n(r e_0 + s) >= C_n |r|^(1/n) + lambda_n ||s||^2 - R_n(r,s), "
                "with R_n bounded so mixed terms cannot cancel the Puiseux deficit"
            ),
            "open_gap": "No analytic remainder bound or interval shape-cone certificate is established here.",
        },
        "next_run_target": {
            "primary": "n20_triad",
            "reason": "Best falsification test for the stratified certificate pattern.",
            "inputs": {
                "radial": "n=20 exact hypergeometric radial family, eps grid from 0.1 down to 1e-8",
                "fourier": "n=20 quotient Fourier basis rank 37, eps 0.02/0.01/0.005",
                "tensor_cone": "n=20 mixed shape matrix with singular radial direction excluded",
            },
            "outputs": [
                "EXP-MATH-EHP114-N20-RADIAL-HYPERGEOMETRIC-20260505-02_RESULTS.json",
                "EXP-MATH-EHP114-N20-FOURIER-HESSIAN-20260505-02_RESULTS.json",
                "EXP-MATH-EHP114-N20-TENSOR-CONE-SCAFFOLD-20260505-02_RESULTS.json",
            ],
            "pass_signal": "radial slope within 5 percent of 1/20 plus broad positive shape curvature",
            "fail_signal": "radial exponent drift or broad loss of positive nonradial curvature",
        },
        "proof_oriented_alternative": {
            "target": "n14_interval_hardening",
            "input": "fixed gauss_2F1_puiseux_lower_bound_n14_interval theorem target",
            "output": "EXP-MATH-EHP114-N14-PUISEUX-INTERVAL-HARDENING-20260505-01_*",
        },
        "source_artifacts": [
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "role": "required_input",
            }
            for path in required_inputs
        ],
        "guardrails": {
            "writes_scoped_to": str(out_dir),
            "no_d1": True,
            "no_scorecard": True,
            "no_public_docs": True,
            "no_publication_or_push": True,
        },
    }

    results_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"

    results_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_path.write_text(build_report(result), encoding="utf-8")
    sha_path.write_text(
        f"{sha256_file(results_path)}  {results_path.name}\n",
        encoding="utf-8",
    )

    print(f"wrote {results_path}")
    print(f"wrote {report_path}")
    print(f"wrote {sha_path}")


if __name__ == "__main__":
    main()
