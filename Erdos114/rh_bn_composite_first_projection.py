#!/usr/bin/env python3
"""
EXP-MATH-RH-BN-COMPOSITE-FIRST-PROJECTION-20260506-01

Finite composite-first projection diagnostic for the Beurling-Nyman MDL lane.

Claim ceiling:
    INTERNAL / METHOD-SHAPING ONLY.
    This is a finite projection/residual experiment. It is not RH support, not
    a zeta formalization, and not a spectral-zeta proof mechanism.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import lstsq, svd

from rh_beurling_nyman_mdl_probe import (
    ACTIVE_REL_TOL,
    DEFAULT_QUANT_BITS,
    DICTIONARIES,
    DictionarySpec,
    fractional_part,
    theta_values,
)


EXPERIMENT_ID = "EXP-MATH-RH-BN-COMPOSITE-FIRST-PROJECTION-20260506-01"
CLAIM_CEILING = (
    "INTERNAL / METHOD-SHAPING ONLY: finite Beurling-Nyman projection "
    "diagnostic; no RH claim, no zeta formalization, and no asymptotic claim."
)
DEFAULT_N_VALUES = (8, 16, 24, 32)
GRID_SPECS = (
    {"name": "legendre_2048", "rule": "Gauss-Legendre on [0,1]", "nodes": 2048},
    {"name": "midpoint_2048", "rule": "equal-weight midpoint grid on [0,1]", "nodes": 2048},
)
PRIMARY_BITS = 16
ORDER_TOL = 1e-8


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def finite_float(x: float) -> float | str:
    if math.isfinite(x):
        return float(x)
    if math.isnan(x):
        return "NaN"
    return "Infinity" if x > 0 else "-Infinity"


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    limit = int(math.isqrt(n))
    for d in range(3, limit + 1, 2):
        if n % d == 0:
            return False
    return True


def combination_bits(n: int, k: int) -> int:
    if k <= 0 or k >= n:
        return 1
    return math.ceil(math.log2(math.comb(n, k)))


def build_grid(name: str, nodes: int) -> tuple[np.ndarray, np.ndarray]:
    if name.startswith("legendre"):
        raw_nodes, raw_weights = leggauss(nodes)
        return 0.5 * (raw_nodes + 1.0), 0.5 * raw_weights
    if name.startswith("midpoint"):
        x = (np.arange(nodes, dtype=np.float64) + 0.5) / nodes
        w = np.full(nodes, 1.0 / nodes, dtype=np.float64)
        return x, w
    raise ValueError(f"unknown grid {name}")


def build_design(theta: np.ndarray, x: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    basis = fractional_part(theta[None, :] / x[:, None])
    sqrt_w = np.sqrt(w)
    return basis * sqrt_w[:, None], sqrt_w


def index_groups(n: int) -> dict[str, list[int]]:
    unit = [0]
    primes = [j - 1 for j in range(2, n + 1) if is_prime(j)]
    composites = [j - 1 for j in range(2, n + 1) if not is_prime(j)]
    return {
        "unit": unit,
        "prime": primes,
        "composite": composites,
        "unit_prime": unit + primes,
        "unit_composite": unit + composites,
        "all": list(range(n)),
    }


def active_count(coeffs: np.ndarray) -> int:
    if coeffs.size == 0:
        return 0
    max_abs = float(np.max(np.abs(coeffs)))
    return int(np.sum(np.abs(coeffs) > ACTIVE_REL_TOL * max(1.0, max_abs)))


def bit_summary(spec: DictionarySpec, subset_size: int, active: int, coeff_bits: int, stage_count: int) -> dict[str, int]:
    support_bits = combination_bits(subset_size, active)
    coefficient_bits = active * coeff_bits
    stage_header_bits = 8 * stage_count
    total = spec.header_bits + stage_header_bits + support_bits + coefficient_bits
    return {
        "dictionary_header_bits": spec.header_bits,
        "stage_header_bits": stage_header_bits,
        "support_pattern_bits": support_bits,
        "coefficient_bits": coefficient_bits,
        "total_description_bits": total,
    }


def quantization_summary(
    spec: DictionarySpec,
    coeffs: np.ndarray,
    singular_values: np.ndarray,
    subset_size: int,
    target_l2: float,
    coeff_bits: int,
    stage_count: int = 1,
) -> dict[str, Any]:
    active = active_count(coeffs)
    scale = max(1.0, float(np.max(np.abs(coeffs))) if coeffs.size else 1.0)
    step = scale * (2.0 ** (-coeff_bits))
    sigma_max = float(singular_values[0]) if singular_values.size else 0.0
    penalty = sigma_max * math.sqrt(active) * step / 2.0
    return {
        "coeff_bits": coeff_bits,
        "active_count": active,
        "quantization_step": step,
        "quantization_penalty_l2_bound": penalty,
        "quantization_penalty_relative_bound": penalty / target_l2 if target_l2 else float("nan"),
        **bit_summary(spec, subset_size, active, coeff_bits, stage_count),
    }


def fit_subset(
    A: np.ndarray,
    y: np.ndarray,
    cols: list[int],
    label: str,
    spec: DictionarySpec,
    coeff_bits: Iterable[int],
) -> dict[str, Any]:
    if not cols:
        residual = y.copy()
        coeffs = np.zeros(0, dtype=np.float64)
        singular_values = np.zeros(0, dtype=np.float64)
        rank = 0
        fitted = np.zeros_like(y)
    else:
        sub = A[:, cols]
        coeffs, _, rank, singular_values = lstsq(sub, y, lapack_driver="gelsd")
        fitted = sub @ coeffs
        residual = y - fitted
    residual_l2 = float(np.linalg.norm(residual))
    target_l2 = float(np.linalg.norm(y))
    condition_number = (
        float(singular_values[0] / singular_values[-1])
        if singular_values.size and singular_values[-1] > 0
        else float("inf")
    )
    q = [
        quantization_summary(spec, coeffs, singular_values, len(cols), target_l2, bits)
        for bits in coeff_bits
    ]
    return {
        "label": label,
        "columns": [c + 1 for c in cols],
        "column_count": len(cols),
        "rank": int(rank),
        "active_count": active_count(coeffs),
        "coefficients_l1": float(np.linalg.norm(coeffs, 1)),
        "coefficients_l2": float(np.linalg.norm(coeffs)),
        "coefficients_linf": float(np.max(np.abs(coeffs))) if coeffs.size else 0.0,
        "singular_values": [float(x) for x in singular_values],
        "sigma_max": float(singular_values[0]) if singular_values.size else 0.0,
        "sigma_min": float(singular_values[-1]) if singular_values.size else 0.0,
        "condition_number": finite_float(condition_number),
        "gram_eigen_max": float(singular_values[0] ** 2) if singular_values.size else 0.0,
        "gram_eigen_min": float(singular_values[-1] ** 2) if singular_values.size else 0.0,
        "residual_l2": residual_l2,
        "relative_residual": residual_l2 / target_l2 if target_l2 else float("nan"),
        "negative_log2_relative_residual": (
            -math.log2(residual_l2 / target_l2) if residual_l2 > 0 and target_l2 else None
        ),
        "primary_quantization": next(item for item in q if item["coeff_bits"] == PRIMARY_BITS),
        "quantized": q,
        "_coefficients": coeffs,
        "_fitted": fitted,
        "_residual": residual,
        "_singular_values_np": singular_values,
    }


def spectral_capture(A: np.ndarray, residual: np.ndarray, cols: list[int]) -> dict[str, Any]:
    if not cols:
        return {
            "added_column_count": 0,
            "rank": 0,
            "capturable_energy": 0.0,
            "capturable_fraction_of_residual": 0.0,
            "top1_fraction_of_capturable_energy": None,
            "top3_fraction_of_capturable_energy": None,
        }
    sub = A[:, cols]
    u, singular_values, _vh = svd(sub, full_matrices=False)
    rank = int(np.sum(singular_values > np.finfo(float).eps * max(sub.shape) * singular_values[0])) if singular_values.size else 0
    coords = u[:, :rank].T @ residual if rank else np.zeros(0, dtype=np.float64)
    energies = np.sort(coords ** 2)[::-1]
    capturable = float(np.sum(energies))
    residual_energy = float(np.dot(residual, residual))
    return {
        "added_column_count": len(cols),
        "rank": rank,
        "singular_values": [float(x) for x in singular_values],
        "capturable_energy": capturable,
        "capturable_fraction_of_residual": capturable / residual_energy if residual_energy else 0.0,
        "top1_fraction_of_capturable_energy": float(energies[0] / capturable) if capturable and len(energies) else None,
        "top3_fraction_of_capturable_energy": float(np.sum(energies[:3]) / capturable) if capturable and len(energies) else None,
    }


def sequential_fit(
    A: np.ndarray,
    y: np.ndarray,
    first_cols: list[int],
    second_cols: list[int],
    label: str,
    spec: DictionarySpec,
) -> dict[str, Any]:
    first = fit_subset(A, y, first_cols, f"{label}_stage1", spec, (PRIMARY_BITS,))
    second = fit_subset(A, first["_residual"], second_cols, f"{label}_stage2", spec, (PRIMARY_BITS,))
    final_residual = first["_residual"] - second["_fitted"]
    target_l2 = float(np.linalg.norm(y))
    final_l2 = float(np.linalg.norm(final_residual))
    q1 = first["primary_quantization"]
    q2 = second["primary_quantization"]
    total_penalty = q1["quantization_penalty_l2_bound"] + q2["quantization_penalty_l2_bound"]
    total_bits = (
        q1["dictionary_header_bits"]
        + q1["stage_header_bits"]
        + q1["support_pattern_bits"]
        + q1["coefficient_bits"]
        + q2["stage_header_bits"]
        + q2["support_pattern_bits"]
        + q2["coefficient_bits"]
    )
    return {
        "label": label,
        "stage1_label": first["label"],
        "stage2_label": second["label"],
        "stage1_relative_residual": first["relative_residual"],
        "stage2_relative_residual_on_stage1_residual": second["relative_residual"],
        "final_residual_l2": final_l2,
        "final_relative_residual": final_l2 / target_l2 if target_l2 else float("nan"),
        "negative_log2_final_relative_residual": (
            -math.log2(final_l2 / target_l2) if final_l2 > 0 and target_l2 else None
        ),
        "residual_drop_from_stage1": first["residual_l2"] - final_l2,
        "residual_drop_from_stage1_relative_to_target": (first["residual_l2"] - final_l2) / target_l2 if target_l2 else float("nan"),
        "residual_drop_per_added_column": (first["residual_l2"] - final_l2) / max(1, len(second_cols)),
        "residual_drop_per_16bit_coefficient_bit": (first["residual_l2"] - final_l2) / max(1, q2["coefficient_bits"]),
        "total_active_count": first["active_count"] + second["active_count"],
        "total_description_bits_16": total_bits,
        "quantization_penalty_l2_bound_16": total_penalty,
        "certified_upper_relative_residual_16": (final_l2 + total_penalty) / target_l2 if target_l2 else float("nan"),
        "stage2_projection_spectrum": spectral_capture(A, first["_residual"], second_cols),
    }


def public_fit(fit: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in fit.items() if not k.startswith("_")}


def run_experiment(n_values: Iterable[int]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for grid in GRID_SPECS:
        x, w = build_grid(grid["name"], grid["nodes"])
        for n in n_values:
            groups = index_groups(n)
            for spec in DICTIONARIES:
                theta = theta_values(spec, n)
                A, y = build_design(theta, x, w)
                unit_composite = fit_subset(A, y, groups["unit_composite"], "unit_plus_composite", spec, DEFAULT_QUANT_BITS)
                unit_prime = fit_subset(A, y, groups["unit_prime"], "unit_plus_prime", spec, DEFAULT_QUANT_BITS)
                all_fit = fit_subset(A, y, groups["all"], "all_columns", spec, DEFAULT_QUANT_BITS)
                composite_first = sequential_fit(A, y, groups["unit_composite"], groups["prime"], "composite_then_prime", spec)
                prime_first = sequential_fit(A, y, groups["unit_prime"], groups["composite"], "prime_then_composite", spec)

                order_gap = composite_first["final_relative_residual"] - prime_first["final_relative_residual"]
                if abs(order_gap) <= ORDER_TOL:
                    order_winner = "tie"
                elif order_gap < 0:
                    order_winner = "composite_first"
                else:
                    order_winner = "prime_first"

                rows.append(
                    {
                        "grid": grid,
                        "dictionary": spec.name,
                        "theta_rule": spec.theta_rule,
                        "N": n,
                        "unit_count": len(groups["unit"]),
                        "prime_column_count": len(groups["prime"]),
                        "composite_column_count": len(groups["composite"]),
                        "fits": {
                            "unit_plus_composite": public_fit(unit_composite),
                            "unit_plus_prime": public_fit(unit_prime),
                            "all_columns": public_fit(all_fit),
                        },
                        "sequential": {
                            "composite_then_prime": composite_first,
                            "prime_then_composite": prime_first,
                        },
                        "order_gap_composite_minus_prime": order_gap,
                        "order_winner": order_winner,
                        "standalone_bulk_winner": (
                            "composite"
                            if unit_composite["relative_residual"] < unit_prime["relative_residual"] - ORDER_TOL
                            else "prime"
                            if unit_prime["relative_residual"] < unit_composite["relative_residual"] - ORDER_TOL
                            else "tie"
                        ),
                        "prime_after_composite_drop_per_column": composite_first["residual_drop_per_added_column"],
                        "composite_after_prime_drop_per_column": prime_first["residual_drop_per_added_column"],
                        "prime_after_composite_drop_per_16bit": composite_first["residual_drop_per_16bit_coefficient_bit"],
                        "composite_after_prime_drop_per_16bit": prime_first["residual_drop_per_16bit_coefficient_bit"],
                    }
                )
    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_unix": int(time.time()),
        "claim_ceiling": CLAIM_CEILING,
        "model_note": "COMPOSITE_FIRST_PROJECTION_MODEL_2026-05-06.md",
        "unit_policy": "j=1 is included in first-stage standalone prime/composite fits; second-stage sequential fits add only the complementary prime or composite columns.",
        "n_values": list(n_values),
        "grid_specs": GRID_SPECS,
        "dictionaries": [asdict(spec) for spec in DICTIONARIES],
        "primary_quantization_bits": PRIMARY_BITS,
        "rows": rows,
    }


def summarize(results: dict[str, Any]) -> dict[str, Any]:
    rows = results["rows"]
    counts = {
        "composite_first": sum(1 for row in rows if row["order_winner"] == "composite_first"),
        "prime_first": sum(1 for row in rows if row["order_winner"] == "prime_first"),
        "tie": sum(1 for row in rows if row["order_winner"] == "tie"),
    }
    bulk_counts = {
        "composite": sum(1 for row in rows if row["standalone_bulk_winner"] == "composite"),
        "prime": sum(1 for row in rows if row["standalone_bulk_winner"] == "prime"),
        "tie": sum(1 for row in rows if row["standalone_bulk_winner"] == "tie"),
    }
    comparable = len(rows) - counts["tie"]
    dominant = "tie"
    if counts["composite_first"] > counts["prime_first"]:
        dominant = "composite_first"
    elif counts["prime_first"] > counts["composite_first"]:
        dominant = "prime_first"
    dominant_count = counts.get(dominant, 0) if dominant != "tie" else counts["tie"]

    by_grid: dict[str, dict[str, int]] = {}
    by_n: dict[str, dict[str, int]] = {}
    for row in rows:
        grid = row["grid"]["name"]
        by_grid.setdefault(grid, {"composite_first": 0, "prime_first": 0, "tie": 0})
        by_grid[grid][row["order_winner"]] += 1
        n_key = str(row["N"])
        by_n.setdefault(n_key, {"composite_first": 0, "prime_first": 0, "tie": 0})
        by_n[n_key][row["order_winner"]] += 1

    grid_agrees = dominant != "tie" and all(v[dominant] > max(v["tie"], v["prime_first" if dominant == "composite_first" else "composite_first"]) for v in by_grid.values())
    n_values_with_dominant = (
        sum(1 for v in by_n.values() if dominant != "tie" and v[dominant] > max(v["tie"], v["prime_first" if dominant == "composite_first" else "composite_first"]))
    )
    stable = bool(
        dominant != "tie"
        and comparable > 0
        and dominant_count / comparable >= 0.75
        and grid_agrees
        and n_values_with_dominant >= 3
    )

    median_abs_gap = float(np.median([abs(row["order_gap_composite_minus_prime"]) for row in rows]))
    median_prime_marginal = float(np.median([row["prime_after_composite_drop_per_column"] for row in rows]))
    median_composite_marginal = float(np.median([row["composite_after_prime_drop_per_column"] for row in rows]))

    return {
        "row_count": len(rows),
        "order_winner_counts": counts,
        "standalone_bulk_winner_counts": bulk_counts,
        "dominant_order": dominant,
        "dominant_order_fraction_of_comparable": dominant_count / comparable if comparable else 0.0,
        "median_abs_order_gap_relative_residual": median_abs_gap,
        "order_counts_by_grid": by_grid,
        "order_counts_by_N": by_n,
        "n_values_with_dominant_order": n_values_with_dominant,
        "stable_asymmetry": stable,
        "stable_asymmetry_status": (
            "STABLE_FINITE_PROJECTION_ASYMMETRY"
            if stable
            else "NO_STABLE_FINITE_PROJECTION_ASYMMETRY"
        ),
        "median_prime_after_composite_drop_per_column": median_prime_marginal,
        "median_composite_after_prime_drop_per_column": median_composite_marginal,
    }


def markdown_report(results: dict[str, Any]) -> str:
    summary = results["summary"]
    selected = []
    for row in results["rows"]:
        if row["N"] in (8, 16, 32) and row["dictionary"] == "geometric":
            selected.append(
                "| {grid} | {N} | {winner} | {comp:.6g} | {prime:.6g} | {allfit:.6g} | {gap:+.3e} |".format(
                    grid=row["grid"]["name"],
                    N=row["N"],
                    winner=row["order_winner"],
                    comp=row["sequential"]["composite_then_prime"]["final_relative_residual"],
                    prime=row["sequential"]["prime_then_composite"]["final_relative_residual"],
                    allfit=row["fits"]["all_columns"]["relative_residual"],
                    gap=row["order_gap_composite_minus_prime"],
                )
            )

    return f"""# {EXPERIMENT_ID} Report

## Verdict

- Status: `{summary["stable_asymmetry_status"]}`
- Rows: `{summary["row_count"]}`
- Dominant order: `{summary["dominant_order"]}`
- Order counts: `{summary["order_winner_counts"]}`
- Standalone bulk counts: `{summary["standalone_bulk_winner_counts"]}`
- Median absolute order gap in relative residual: `{summary["median_abs_order_gap_relative_residual"]:.6g}`

## Meaning

This experiment tests a finite projection model for the RH-MDL lane. The
question is whether composite-indexed Beurling-Nyman columns absorb bulk
residual first, while prime-indexed columns supply harder marginal directions.

The result is a finite linear-algebra diagnostic. It is not a zeta
formalization, not a theorem about primes, and not a public RH claim.

## Claim Ceiling

{CLAIM_CEILING}

## Model

- Target: `chi_(0,1] = 1` on `(0,1]`.
- Basis: `rho_theta(x) = fractional_part(theta / x)`.
- Column index split: `1` is unit, prime indices are prime columns, composite
  indices greater than `1` are composite columns.
- Grids: `{", ".join(grid["name"] for grid in GRID_SPECS)}`.
- Dictionary sizes: `{results["n_values"]}`.
- Dictionaries: `{", ".join(spec["name"] for spec in results["dictionaries"])}`.

## Selected Geometric Rows

| grid | N | order winner | composite-then-prime residual | prime-then-composite residual | all-column residual | composite-minus-prime gap |
|---|---:|---|---:|---:|---:|---:|
{chr(10).join(selected)}

## Stability Checks

The protocol required the order asymmetry to be stable across at least three
dictionary sizes, two quadrature grids, and both orderings. The observed status
is:

```text
{summary["stable_asymmetry_status"]}
```

Grid counts:

```json
{json.dumps(summary["order_counts_by_grid"], indent=2)}
```

N counts:

```json
{json.dumps(summary["order_counts_by_N"], indent=2)}
```

## Spectral Diagnostic Boundary

Each row stores finite projection-spectrum data for the second-stage residual:
capturable energy, the fraction of residual energy captured, and top-mode
concentration. This is only a finite matrix diagnostic. It should not be
identified with the Riemann zeta function or with a Hilbert-Polya operator.

## Artifact Boundary

Generated artifacts:

- `{EXPERIMENT_ID}_RESULTS.json`
- `{EXPERIMENT_ID}_REPORT.md`
- `{EXPERIMENT_ID}_RESULTS.sha256`

No D1, scorecard, public page, git staging, commit, Zenodo, arXiv, or publisher
surface was updated.
"""


def write_outputs(results: dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    result_path = output_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = output_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = output_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite existing immutable artifact: {path}")

    summary = summarize(results)
    payload = {**results, "summary": summary}
    results_bytes = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8")
    result_path.write_bytes(results_bytes)
    digest = hashlib.sha256(results_bytes).hexdigest()
    sha_path.write_text(f"{digest}  {result_path.name}\n", encoding="utf-8")
    report_path.write_text(markdown_report(payload), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=here, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results = run_experiment(DEFAULT_N_VALUES)
    write_outputs(results, args.output_dir)
    print(
        json.dumps(
            {
                "experiment_id": EXPERIMENT_ID,
                "status": summarize(results)["stable_asymmetry_status"],
                "output_dir": str(args.output_dir),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
