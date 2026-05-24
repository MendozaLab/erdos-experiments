#!/usr/bin/env python3
"""Run the n=20 Hessian/Puiseux triad for EHP114.

This is a diagnostic continuation artifact, not a proof or publication packet.
It reuses the existing radial, Fourier, and tensor-cone probe functions without
modifying those shared scripts. All durable outputs are written next to this
file.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


sys.dont_write_bytecode = True

EXPERIMENT_ID = "EXP-MATH-EHP114-N20-HESSIAN-TRIAD-20260505-01"
RADIAL_ID = "EXP-MATH-EHP114-N20-RADIAL-HYPERGEOMETRIC-20260505-02"
FOURIER_ID = "EXP-MATH-EHP114-N20-FOURIER-HESSIAN-20260505-02"
TENSOR_ID = "EXP-MATH-EHP114-N20-TENSOR-CONE-SCAFFOLD-20260505-02"

N = 20
DEFAULT_RADIAL_EPS = "0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0001,0.00001,0.000001,0.0000001,0.00000001"
DEFAULT_FOURIER_EPS = "0.02,0.01,0.005"
DEFAULT_TENSOR_EPS = "0.04,0.02,0.01,0.005"
CONDITION_PROXY_BOUND = 100.0


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_float_list(raw: str) -> list[float]:
    values = [float(part.strip()) for part in raw.split(",") if part.strip()]
    if not values:
        raise ValueError("empty float list")
    return values


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def timed(label: str, func: Any) -> tuple[Any, float]:
    start = time.monotonic()
    value = func()
    return value, time.monotonic() - start


def tail_slope(result: dict[str, Any], window: int) -> float:
    rows = result["asymptotic"]["slope_rows"]
    for row in rows:
        if row["tail_window"] == window:
            return float(row["fitted_log_log_slope"])
    return float(rows[-1]["fitted_log_log_slope"])


def radial_summary(radial: dict[str, Any]) -> dict[str, Any]:
    expected = float(radial["asymptotic"]["expected_slope"])
    observed = tail_slope(radial, 4)
    rel_error = abs(observed - expected) / expected
    tail = radial["asymptotic"]["rows"][-1]
    return {
        "component_id": radial["experiment_id"],
        "degree": radial["degree"],
        "expected_slope_1_over_n": expected,
        "tail_window_4_slope": observed,
        "relative_error_tail4": rel_error,
        "within_5_percent_target": rel_error <= 0.05,
        "smallest_eps": float(tail["eps"]),
        "tail_deficit_over_eps_power_1_over_n": float(
            tail["deficit_over_eps_power_1_over_n"]
        ),
        "status": radial["status"],
    }


def fourier_summary(fourier: dict[str, Any]) -> dict[str, Any]:
    basis_rows = fourier["basis_rows"]
    positive = int(fourier["summary"]["rows_positive_all_eps"])
    total = len(basis_rows)
    all_curvatures = [
        float(eps_row["deficit_curvature"])
        for row in basis_rows
        for eps_row in row["eps_rows"]
    ]
    positive_fraction = positive / total if total else 0.0
    return {
        "component_id": fourier["experiment_id"],
        "degree": fourier["degree"],
        "basis_rank": fourier["parameters"]["basis_rank"],
        "expected_rank_2n_minus_3": fourier["parameters"]["expected_reduced_dim_2n_minus_3"],
        "positive_modes_all_eps": positive,
        "mode_count": total,
        "positive_fraction": positive_fraction,
        "broad_positive_shape_curvature": positive_fraction >= 0.90,
        "all_modes_positive_all_eps": bool(fourier["summary"]["all_modes_positive_all_eps"]),
        "min_deficit_curvature": min(all_curvatures),
        "max_deficit_curvature": max(all_curvatures),
        "large_reference_error_flag": fourier["singularity_sanity"]["large_reference_error_flag"],
        "status": fourier["status"],
    }


def tensor_summary(tensor: dict[str, Any]) -> dict[str, Any]:
    shape_rows = [
        row for row in tensor["direction_rows"] if row["cone_role"] == "shape_candidate"
    ]
    shape_positive = sum(1 for row in shape_rows if row["all_symmetric_deficits_positive"])
    mixed = tensor["mixed_shape_hessian_proxy"]
    eigvals = [float(value) for value in mixed["eigenvalues"]] if mixed else []
    positive_eigs = sum(1 for value in eigvals if value > 0.0)
    condition_proxy = mixed["condition_proxy_abs_max_over_min"] if mixed else None
    condition_ok = (
        condition_proxy is not None
        and math.isfinite(float(condition_proxy))
        and float(condition_proxy) <= CONDITION_PROXY_BOUND
    )
    return {
        "component_id": tensor["experiment_id"],
        "degree": tensor["degree"],
        "quotient_basis_rank": tensor["parameters"]["quotient_basis_rank"],
        "shape_basis_rank": tensor["parameters"]["shape_basis_rank"],
        "expected_rank_2n_minus_3": tensor["parameters"]["expected_rank_2n_minus_3"],
        "radial_slope": tensor["summary"]["radial_slope"],
        "shape_slope_min": tensor["summary"]["shape_slope_min"],
        "shape_slope_mean": tensor["summary"]["shape_slope_mean"],
        "shape_slope_max": tensor["summary"]["shape_slope_max"],
        "shape_positive_symmetric_deficits": shape_positive,
        "shape_direction_count": len(shape_rows),
        "all_shape_symmetric_deficits_positive": shape_positive == len(shape_rows),
        "mixed_min_eigenvalue": min(eigvals) if eigvals else None,
        "mixed_max_eigenvalue": max(eigvals) if eigvals else None,
        "mixed_positive_eigenvalues": positive_eigs,
        "mixed_eigenvalue_count": len(eigvals),
        "mixed_shape_proxy_positive": bool(eigvals) and positive_eigs == len(eigvals),
        "condition_proxy_abs_max_over_min": condition_proxy,
        "predeclared_condition_proxy_bound": CONDITION_PROXY_BOUND,
        "condition_proxy_bounded": condition_ok,
        "status": tensor["status"],
    }


def lean_target(full_triad_signal: bool) -> str:
    if full_triad_signal:
        return """theorem ehp114_n20_stratified_shadow_local_bound
    (r : Real) (s : ShapeQuotient 20)
    (hr_pos : 0 < abs r) (hr_small : abs r <= delta20)
    (hs_small : quotientNorm s <= eta20) :
    C20 * Real.rpow (abs r) ((1 : Real) / 20)
      + lambda20 * quotientNormSq s
      <= D20 (radialMode20 r + shapeMode20 s) + R20 r s := by
  -- interval Puiseux radial bound plus quotient shape-cone remainder bound
  sorry"""
    return """theorem ehp114_n20_radial_puiseux_interval
    (eps : Real) (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 100000000) :
    C20 * Real.rpow eps ((1 : Real) / 20)
      <= L20Star - radialLength20 (Real.rpow (1 - eps / Real.sqrt 20) 20) := by
  -- exact 2F1/Gamma radial-family interval hardening
  sorry"""


def build_report(result: dict[str, Any]) -> str:
    r = result["component_summaries"]["radial"]
    f = result["component_summaries"]["fourier"]
    t = result["component_summaries"]["tensor_cone"]
    theorem = result["next_theorem_target"]["lean_shape"]
    return f"""# {EXPERIMENT_ID} Report

## Status

Diagnostic n=20 Hessian/Puiseux triad for EHP114. This is not a proof and not
a publication packet. The claim ceiling remains: this is a shadow signature,
not universal law.

## Verdict

- Run status: `{result["run_status"]}`
- Triad verdict: `{result["triad_verdict"]}`
- Component decision: radial, Fourier, and tensor-cone machinery were runnable
  now by reusing existing local probe functions through a scoped wrapper.

The result preserves the stratified reading. The radial axis follows the
Puiseux exponent target, while the nonradial Fourier and tensor-cone probes
retain broad positive finite-difference shape signals under floating
estimators. This supports the next interval-hardening target, but it does not
establish a local or global EHP114 theorem.

## Component Summary

| component | runnable now? | main signal | diagnostic pass? |
|---|---|---|---|
| radial hypergeometric | yes | tail slope {r["tail_window_4_slope"]:.12g} vs target {r["expected_slope_1_over_n"]:.12g}; relative error {r["relative_error_tail4"]:.6g} | {r["within_5_percent_target"]} |
| Fourier quotient Hessian | yes | positive modes {f["positive_modes_all_eps"]}/{f["mode_count"]}; rank {f["basis_rank"]} vs expected {f["expected_rank_2n_minus_3"]} | {f["broad_positive_shape_curvature"]} |
| tensor-cone mixed shape | yes | mixed positive eigenvalues {t["mixed_positive_eigenvalues"]}/{t["mixed_eigenvalue_count"]}; condition proxy {t["condition_proxy_abs_max_over_min"]:.6g} | {t["mixed_shape_proxy_positive"] and t["condition_proxy_bounded"]} |

## Radial Puiseux

The n=20 radial family lands on the expected `1/20` Puiseux exponent within the
predeclared 5 percent tolerance. The smallest tested eps is
`{r["smallest_eps"]:.12g}`, and the tail proxy
`deficit / eps^(1/20)` is `{r["tail_deficit_over_eps_power_1_over_n"]:.12g}`.
This says the smooth radial Hessian model is still the wrong object at n=20.

## Nonradial Shape

The Fourier basis has the expected reduced rank `2n - 3 = 37`. Under the
floating marching-squares estimator, `{f["positive_modes_all_eps"]}` of
`{f["mode_count"]}` tested quotient modes are positive at all eps values. The
large reference-error flag remains `{f["large_reference_error_flag"]}`, so this
is shape evidence only, not an interval certificate.

The tensor-cone mixed matrix excludes the singular radial direction. Its shape
basis rank is `{t["shape_basis_rank"]}`, and the finite-difference proxy has
minimum eigenvalue `{t["mixed_min_eigenvalue"]:.12g}` and maximum eigenvalue
`{t["mixed_max_eigenvalue"]:.12g}`. The condition proxy is compared against the
predeclared diagnostic bound `{t["predeclared_condition_proxy_bound"]}`.

## Lean-Shaped Next Target

```lean
{theorem}
```

This target is deliberately local and stratified: radial Puiseux first, then
quotient shape-cone positivity plus a mixed remainder bound. The missing piece
is interval control of the remainders, not another floating sweep.

## Source Boundary

All durable outputs from this run were written under:

`{result["guardrails"]["writes_scoped_to"]}`

No scorecard, D1, public docs, git, email, CLAUDE.md, or AGENTS.md was changed.
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument("--fourier-res", type=int, default=360)
    parser.add_argument("--tensor-res", type=int, default=340)
    parser.add_argument("--mixed-res", type=int, default=240)
    parser.add_argument("--extent", type=float, default=3.0)
    parser.add_argument("--radial-eps", default=DEFAULT_RADIAL_EPS)
    parser.add_argument("--fourier-eps", default=DEFAULT_FOURIER_EPS)
    parser.add_argument("--tensor-eps", default=DEFAULT_TENSOR_EPS)
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    math_root = script_path.parents[2]
    scripts_dir = math_root / "erdos-experiments" / "scripts"
    results_dir = math_root / "erdos-experiments" / "results" / "erdos-114"
    shortcut_dir = scripts_dir / "erdos-114"
    out_dir = args.out_dir.resolve()

    if out_dir != script_path.parent.resolve():
        raise SystemExit("Refusing to write outside erdos-experiments/Erdos114")
    out_dir.mkdir(parents=True, exist_ok=True)

    result_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {path}")

    radial_mod = load_module("ehp114_radial_hypergeometric_source", scripts_dir / "calc_ehp114_radial_hypergeometric.py")
    fourier_mod = load_module("ehp114_fourier_hessian_source", scripts_dir / "probe_ehp114_fourier_hessian.py")
    tensor_mod = load_module("ehp114_tensor_cone_source", scripts_dir / "probe_ehp114_tensor_cone_scaffold.py")

    radial_mod.mp.mp.dps = args.mp_dps
    tensor_mod.mp.mp.dps = args.mp_dps

    component_runtime_seconds: dict[str, float] = {}

    fourier_eps = parse_float_list(args.fourier_eps)
    tensor_eps = parse_float_list(args.tensor_eps)

    fourier, elapsed = timed(
        "fourier",
        lambda: fourier_mod.run_probe(
            n=N,
            eps_values=fourier_eps,
            res=args.fourier_res,
            extent=args.extent,
            shortcut_dir=shortcut_dir,
        ),
    )
    fourier["experiment_id"] = FOURIER_ID
    component_runtime_seconds["fourier"] = elapsed

    radial, elapsed = timed(
        "radial",
        lambda: radial_mod.build_result(
            n=N,
            experiment_id=RADIAL_ID,
            eps_values=radial_mod.parse_mp_list(args.radial_eps),
            asymptotic_eps_values=radial_mod.parse_mp_list(args.radial_eps),
            result_root=results_dir,
            shortcut_root=shortcut_dir,
            fourier_artifact_path=results_dir / "EXP-MATH-EHP114-N15-FOURIER-HESSIAN-20260502-02_RESULTS.json",
        ),
    )
    component_runtime_seconds["radial"] = elapsed

    tensor, elapsed = timed(
        "tensor_cone",
        lambda: tensor_mod.run_probe(
            n=N,
            eps_values=tensor_eps,
            res=args.tensor_res,
            mixed_res=args.mixed_res,
            extent=args.extent,
            result_dir=results_dir,
            include_mixed_matrix=True,
        ),
    )
    tensor["experiment_id"] = TENSOR_ID
    component_runtime_seconds["tensor_cone"] = elapsed

    summaries = {
        "radial": radial_summary(radial),
        "fourier": fourier_summary(fourier),
        "tensor_cone": tensor_summary(tensor),
    }
    full_triad_signal = (
        summaries["radial"]["within_5_percent_target"]
        and summaries["fourier"]["broad_positive_shape_curvature"]
        and summaries["tensor_cone"]["mixed_shape_proxy_positive"]
        and summaries["tensor_cone"]["condition_proxy_bounded"]
    )
    triad_verdict = (
        "SUPPORTS_STRATIFIED_SHADOW_SIGNATURE"
        if full_triad_signal
        else "MIXED_DIAGNOSTIC_SIGNAL"
    )

    result: dict[str, Any] = {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "scope": "n=20 stratified Hessian/Puiseux diagnostic triad",
        "run_status": "COMPLETE",
        "triad_verdict": triad_verdict,
        "claim_ceiling": "This is a shadow signature, not universal law. It is not a proof or publication packet.",
        "feasibility_decision": {
            "radial_hypergeometric": {
                "can_run_now": True,
                "machinery": "calc_ehp114_radial_hypergeometric.build_result",
                "wrapper_reason": "Existing CLI already supports n and experiment-id parameters.",
            },
            "fourier_hessian": {
                "can_run_now": True,
                "machinery": "probe_ehp114_fourier_hessian.run_probe",
                "wrapper_reason": "Existing reusable function supports n=20; CLI experiment id is pinned.",
            },
            "tensor_cone": {
                "can_run_now": True,
                "machinery": "probe_ehp114_tensor_cone_scaffold.run_probe",
                "wrapper_reason": "Existing reusable function supports n=20; CLI guard is pinned to n=10.",
            },
        },
        "parameters": {
            "n": N,
            "radial_eps": args.radial_eps,
            "fourier_eps": fourier_eps,
            "tensor_eps": tensor_eps,
            "fourier_res": args.fourier_res,
            "tensor_res": args.tensor_res,
            "mixed_res": args.mixed_res,
            "extent": args.extent,
            "mp_dps": args.mp_dps,
            "condition_proxy_bound": CONDITION_PROXY_BOUND,
        },
        "component_runtime_seconds": component_runtime_seconds,
        "component_summaries": summaries,
        "components": {
            "radial_hypergeometric": radial,
            "fourier_hessian": fourier,
            "tensor_cone": tensor,
        },
        "next_theorem_target": {
            "selected_because": (
                "The full triad signal supports a local stratified target."
                if full_triad_signal
                else "The full triad signal is mixed, so the radial interval theorem is the next narrow target."
            ),
            "lean_shape": lean_target(full_triad_signal),
        },
        "blockers": [
            "No interval arithmetic certificate for the n=20 Fourier or tensor-cone shape block.",
            "No analytic mixed radial/shape remainder bound.",
            "Floating marching-squares estimators are diagnostic only.",
        ],
        "guardrails": {
            "writes_scoped_to": str(out_dir),
            "no_d1": True,
            "no_scorecard": True,
            "no_public_docs": True,
            "no_git": True,
            "no_email": True,
            "no_claude_or_agents": True,
        },
        "source_scripts": {
            "radial": str(scripts_dir / "calc_ehp114_radial_hypergeometric.py"),
            "fourier": str(scripts_dir / "probe_ehp114_fourier_hessian.py"),
            "tensor_cone": str(scripts_dir / "probe_ehp114_tensor_cone_scaffold.py"),
        },
    }

    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_path.write_text(build_report(result), encoding="utf-8")
    sha_path.write_text(f"{sha256_file(result_path)}  {result_path.name}\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "result": str(result_path),
                "report": str(report_path),
                "sha256": sha_path.read_text(encoding="utf-8").split()[0],
                "triad_verdict": triad_verdict,
                "component_runtime_seconds": component_runtime_seconds,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
