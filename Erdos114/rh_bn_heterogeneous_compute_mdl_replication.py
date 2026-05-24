#!/usr/bin/env python3
"""
EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-REPLICATION-20260507-01

Replication gate for the finite heterogeneous-compute MDL diagnostic.

Claim ceiling:
    INTERNAL / METHOD-SHAPING ONLY.
    This is a finite operation-weighted Beurling-Nyman MDL robustness check.
    It is not an RH claim, not a zeta formalization, and not a
    quantum-mechanical claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import date
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import lstsq

from rh_beurling_nyman_mdl_probe import (
    ACTIVE_REL_TOL,
    DICTIONARIES,
    DictionarySpec,
    fractional_part,
    theta_values,
)
from rh_bn_prime_harness_mdl import (
    RESIDUAL_TOLERANCES,
    composite_address_bits,
    composite_prefix,
    factor_cost_ordered,
    factor_harness_bits,
    factorization,
    factorized_address_bits,
    flat_address_bits,
    int_bits_at_most,
    integer_prefix,
    primes_up_to,
)


EXPERIMENT_ID = "EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-REPLICATION-20260507-01"
CLAIM_CEILING = (
    "INTERNAL / METHOD-SHAPING ONLY: finite heterogeneous-compute "
    "Beurling-Nyman MDL replication gate; no RH claim, no zeta "
    "formalization, and no quantum claim."
)
MODEL_NOTES = (
    "HETEROGENEOUS_COMPUTE_MDL_MODEL_2026-05-06.md",
    "HETEROGENEOUS_COMPUTE_MDL_EXPERIMENTAL_NOTE_2026-05-06.md",
)
SOURCE_EXPERIMENT = "EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01"
DEFAULT_N_VALUES = (16, 24, 32, 48, 64, 96)
GRID_SPECS = (
    {"name": "legendre_2048", "rule": "Gauss-Legendre on [0,1]", "nodes": 2048},
    {"name": "midpoint_2048", "rule": "equal-weight midpoint grid on [0,1]", "nodes": 2048},
    {
        "name": "shifted_midpoint_2048",
        "rule": "equal-weight midpoint grid shifted to quarter-cell nodes",
        "nodes": 2048,
    },
)
SUPPORT_BUDGETS = (4, 8, 12, 16, 24, 32, 48, 64, 96)
COEFFICIENT_BITS = (8, 16, 24)
BASELINE_ENCODINGS = {"flat_index", "composite_index", "factorized_reusable_harness"}
ABLATION_CHANNELS = ("prime_lookup", "multiply", "exponent", "factor_tree_depth")
STAGE_HEADER_BITS = 8
FACTOR_TREE_DEPTH_WEIGHT = 1
MULTIPLY_WEIGHT = 2
EXPONENT_WEIGHT = 2
ORDER_TOL = 1e-9


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


def build_grid(name: str, nodes: int) -> tuple[np.ndarray, np.ndarray]:
    if name.startswith("legendre"):
        raw_nodes, raw_weights = leggauss(nodes)
        return 0.5 * (raw_nodes + 1.0), 0.5 * raw_weights
    if name == "midpoint_2048" or name.startswith("midpoint"):
        x = (np.arange(nodes, dtype=np.float64) + 0.5) / nodes
        w = np.full(nodes, 1.0 / nodes, dtype=np.float64)
        return x, w
    if name.startswith("shifted_midpoint"):
        x = (np.arange(nodes, dtype=np.float64) + 0.25) / nodes
        w = np.full(nodes, 1.0 / nodes, dtype=np.float64)
        return x, w
    raise ValueError(f"unknown grid {name}")


def build_design(theta: np.ndarray, x: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    basis = fractional_part(theta[None, :] / x[:, None])
    sqrt_w = np.sqrt(w)
    return basis * sqrt_w[:, None], sqrt_w


def active_count(coeffs: np.ndarray) -> int:
    if coeffs.size == 0:
        return 0
    max_abs = float(np.max(np.abs(coeffs)))
    return int(np.sum(np.abs(coeffs) > ACTIVE_REL_TOL * max(1.0, max_abs)))


def fit_columns_base(A: np.ndarray, y: np.ndarray, cols_zero: list[int]) -> dict[str, Any]:
    if not cols_zero:
        coeffs = np.zeros(0, dtype=np.float64)
        residual = y.copy()
        rank = 0
        singular_values = np.zeros(0, dtype=np.float64)
    else:
        sub = A[:, cols_zero]
        coeffs, _, rank, singular_values = lstsq(sub, y, lapack_driver="gelsd")
        residual = y - sub @ coeffs
    residual_l2 = float(np.linalg.norm(residual))
    target_l2 = float(np.linalg.norm(y))
    sigma_max = float(singular_values[0]) if singular_values.size else 0.0
    sigma_min = float(singular_values[-1]) if singular_values.size else 0.0
    cond = sigma_max / sigma_min if sigma_min > 0 else float("inf")
    return {
        "coeffs": coeffs,
        "rank": int(rank),
        "residual_l2": residual_l2,
        "relative_residual": residual_l2 / target_l2 if target_l2 else float("nan"),
        "target_l2": target_l2,
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "condition_number": finite_float(cond),
    }


def fit_metrics(base: dict[str, Any], coeff_bits: int) -> dict[str, Any]:
    coeffs = base["coeffs"]
    active = active_count(coeffs)
    scale = max(1.0, float(np.max(np.abs(coeffs))) if coeffs.size else 1.0)
    step = scale * (2.0 ** (-coeff_bits))
    quant_penalty = base["sigma_max"] * math.sqrt(active) * step / 2.0
    target_l2 = base["target_l2"]
    certified = (base["residual_l2"] + quant_penalty) / target_l2 if target_l2 else float("nan")
    info = -math.log2(certified) if certified > 0 else None
    return {
        "rank": base["rank"],
        "active_count": active,
        "residual_l2": base["residual_l2"],
        "relative_residual": base["relative_residual"],
        "sigma_max": base["sigma_max"],
        "sigma_min": base["sigma_min"],
        "condition_number": base["condition_number"],
        "coefficient_l1": float(np.linalg.norm(coeffs, 1)),
        "coefficient_l2": float(np.linalg.norm(coeffs)),
        "coefficient_linf": float(np.max(np.abs(coeffs))) if coeffs.size else 0.0,
        "coefficient_precision_bits": coeff_bits,
        "quantization_step": step,
        "quantization_penalty_l2_bound": quant_penalty,
        "certified_upper_relative_residual": certified,
        f"certified_upper_relative_residual_{coeff_bits}": certified,
        "information_bits_certified": info,
    }


def total_bits(spec: DictionarySpec, address_bits: int, active: int, coeff_bits: int) -> int:
    return spec.header_bits + STAGE_HEADER_BITS + address_bits + active * coeff_bits


def information_bits(row: dict[str, Any]) -> float | None:
    residual = row.get("certified_upper_relative_residual")
    if isinstance(residual, (int, float)) and residual > 0:
        return -math.log2(float(residual))
    return None


def attach_objective(row: dict[str, Any], compute_cost_bits: int, channels: dict[str, int]) -> dict[str, Any]:
    info = information_bits(row)
    row["compute_cost_bits"] = compute_cost_bits
    row["heterogeneous_cost_channels"] = channels
    row["certified_residual_information_bits"] = info
    row["total_objective_bits"] = compute_cost_bits - info if info is not None else None
    return row


def base_row(
    spec: DictionarySpec,
    n: int,
    grid: str,
    support_strategy: str,
    scenario: str,
    encoding: str,
    cols: list[int],
    fit: dict[str, Any],
    address_bits: int,
    include_harness: bool,
) -> dict[str, Any]:
    coeff_bits = int(fit["coefficient_precision_bits"])
    bits = total_bits(spec, address_bits, fit["active_count"], coeff_bits)
    return {
        "dictionary": spec.name,
        "theta_rule": spec.theta_rule,
        "grid": grid,
        "N": n,
        "support_budget": len(cols),
        "support_strategy": support_strategy,
        "scenario": scenario,
        "encoding": encoding,
        "channel_profile": scenario,
        "include_prime_harness_setup_cost": include_harness,
        "columns": cols,
        "address_bits": int(address_bits),
        "dictionary_header_bits": spec.header_bits,
        "stage_header_bits": STAGE_HEADER_BITS,
        "coefficient_bits": fit["active_count"] * coeff_bits,
        "total_description_bits": bits,
        "description_bits_per_information_bit": (
            bits / fit["information_bits_certified"]
            if fit["information_bits_certified"] and fit["information_bits_certified"] > 0
            else None
        ),
        **fit,
    }


def baseline_channels(encoding: str, address_bits: int, row: dict[str, Any]) -> dict[str, int]:
    coefficient_bits = int(row["coefficient_bits"])
    if encoding in {"flat_index", "composite_index"}:
        return {
            "direct_index": int(address_bits),
            "prime_lookup": 0,
            "multiply": 0,
            "exponent": 0,
            "factor_tree_depth": 0,
            "coefficient_move": coefficient_bits,
            "shared_harness": 0,
        }
    return {
        "direct_index": 0,
        "prime_lookup": int(address_bits),
        "multiply": 0,
        "exponent": 0,
        "factor_tree_depth": 0,
        "coefficient_move": coefficient_bits,
        "shared_harness": 0,
    }


def operation_channels(
    indices: list[int],
    n: int,
    active: int,
    coeff_bits: int,
    include_harness: bool,
    drop_channels: set[str] | None = None,
) -> dict[str, int]:
    drop_channels = drop_channels or set()
    distinct_primes: set[int] = set()
    multiply_ops = 0
    exponent_ops = 0
    max_depth = 0
    for j in indices:
        factors = factorization(j)
        if not factors:
            continue
        distinct_primes.update(factors)
        multiplicity = sum(factors.values())
        multiply_ops += max(0, multiplicity - 1)
        exponent_ops += sum(max(0, exp - 1) for exp in factors.values())
        max_depth = max(max_depth, multiplicity)
    prime_lookup_unit = int_bits_at_most(len(primes_up_to(n)))
    channels = {
        "direct_index": 0,
        "prime_lookup": len(distinct_primes) * prime_lookup_unit,
        "multiply": multiply_ops * MULTIPLY_WEIGHT,
        "exponent": exponent_ops * EXPONENT_WEIGHT,
        "factor_tree_depth": max_depth * FACTOR_TREE_DEPTH_WEIGHT,
        "coefficient_move": active * coeff_bits,
        "shared_harness": factor_harness_bits(n) if include_harness else 0,
    }
    for channel in drop_channels:
        channels[channel] = 0
    return channels


def cost_from_channels(spec: DictionarySpec, channels: dict[str, int]) -> int:
    return spec.header_bits + STAGE_HEADER_BITS + sum(channels.values())


def add_baseline_row(
    rows: list[dict[str, Any]],
    spec: DictionarySpec,
    n: int,
    grid: str,
    support_strategy: str,
    scenario: str,
    encoding: str,
    cols: list[int],
    fit: dict[str, Any],
    address_bits: int,
) -> None:
    row = base_row(
        spec,
        n,
        grid,
        support_strategy,
        scenario,
        encoding,
        cols,
        fit,
        address_bits,
        False,
    )
    channels = baseline_channels(encoding, address_bits, row)
    rows.append(attach_objective(row, int(row["total_description_bits"]), channels))


def add_heterogeneous_row(
    rows: list[dict[str, Any]],
    spec: DictionarySpec,
    n: int,
    grid: str,
    support_strategy: str,
    scenario: str,
    encoding: str,
    cols: list[int],
    fit: dict[str, Any],
    include_harness: bool,
    drop_channels: set[str] | None = None,
) -> None:
    channels = operation_channels(
        cols,
        n,
        fit["active_count"],
        int(fit["coefficient_precision_bits"]),
        include_harness,
        drop_channels,
    )
    compute_cost = cost_from_channels(spec, channels)
    address_bits = sum(channels.values()) - channels["coefficient_move"]
    row = base_row(
        spec,
        n,
        grid,
        support_strategy,
        scenario,
        encoding,
        cols,
        fit,
        address_bits,
        include_harness,
    )
    row["total_description_bits"] = compute_cost
    rows.append(attach_objective(row, compute_cost, channels))


def add_scenario_baselines(
    rows: list[dict[str, Any]],
    spec: DictionarySpec,
    n: int,
    grid: str,
    support_strategy: str,
    scenario: str,
    cols: list[int],
    fit: dict[str, Any],
) -> None:
    if support_strategy == "integer_prefix":
        add_baseline_row(
            rows,
            spec,
            n,
            grid,
            support_strategy,
            scenario,
            "flat_index",
            cols,
            fit,
            flat_address_bits(cols, n),
        )
        add_baseline_row(
            rows,
            spec,
            n,
            grid,
            support_strategy,
            scenario,
            "factorized_reusable_harness",
            cols,
            fit,
            factorized_address_bits(cols, n, include_harness=False),
        )
    elif support_strategy == "composite_prefix":
        bits = composite_address_bits(cols, n)
        if bits is not None:
            add_baseline_row(
                rows,
                spec,
                n,
                grid,
                support_strategy,
                scenario,
                "composite_index",
                cols,
                fit,
                bits,
            )
    elif support_strategy == "factor_cost_ordered":
        add_baseline_row(
            rows,
            spec,
            n,
            grid,
            support_strategy,
            scenario,
            "factorized_reusable_harness",
            cols,
            fit,
            factorized_address_bits(cols, n, include_harness=False),
        )
    else:
        raise ValueError(f"unknown support strategy {support_strategy}")


def add_rows_for_support(
    rows: list[dict[str, Any]],
    spec: DictionarySpec,
    n: int,
    grid: str,
    A: np.ndarray,
    y: np.ndarray,
    support_strategy: str,
    cols: list[int],
) -> None:
    cols = sorted(dict.fromkeys(cols))
    fit_base = fit_columns_base(A, y, [j - 1 for j in cols])
    for coeff_bits in COEFFICIENT_BITS:
        fit = fit_metrics(fit_base, coeff_bits)

        add_scenario_baselines(rows, spec, n, grid, support_strategy, "base", cols, fit)
        add_heterogeneous_row(
            rows,
            spec,
            n,
            grid,
            support_strategy,
            "base",
            "heterogeneous_compute_reusable",
            cols,
            fit,
            include_harness=False,
        )
        add_heterogeneous_row(
            rows,
            spec,
            n,
            grid,
            support_strategy,
            "base",
            "heterogeneous_compute_with_harness",
            cols,
            fit,
            include_harness=True,
        )

        for channel in ABLATION_CHANNELS:
            scenario = f"drop_{channel}"
            add_scenario_baselines(rows, spec, n, grid, support_strategy, scenario, cols, fit)
            add_heterogeneous_row(
                rows,
                spec,
                n,
                grid,
                support_strategy,
                scenario,
                f"heterogeneous_compute_reusable_drop_{channel}",
                cols,
                fit,
                include_harness=False,
                drop_channels={channel},
            )


def support_budgets(n: int) -> list[int]:
    values = sorted({m for m in SUPPORT_BUDGETS if m <= n} | {n})
    return [m for m in values if m > 0]


def tolerance_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, str, int, int, float], list[dict[str, Any]]] = {}
    for row in rows:
        for tol in RESIDUAL_TOLERANCES:
            key = (
                row["scenario"],
                row["grid"],
                row["dictionary"],
                row["N"],
                row["coefficient_precision_bits"],
                tol,
            )
            groups.setdefault(key, []).append(row)

    table = []
    for (scenario, grid, dictionary, n, coeff_bits, tol), group_rows in sorted(groups.items()):
        candidates = [
            row
            for row in group_rows
            if row["certified_upper_relative_residual"] <= tol
            and row["total_objective_bits"] is not None
        ]
        if not candidates:
            continue
        winner = min(candidates, key=lambda row: row["total_objective_bits"])
        by_encoding = {}
        for row in candidates:
            encoding = row["encoding"]
            best = by_encoding.get(encoding)
            if best is None or row["total_objective_bits"] < best["total_objective_bits"]:
                by_encoding[encoding] = row
        table.append(
            {
                "scenario": scenario,
                "grid": grid,
                "dictionary": dictionary,
                "N": n,
                "coefficient_precision_bits": coeff_bits,
                "tolerance": tol,
                "winner_encoding": winner["encoding"],
                "winner_support_strategy": winner["support_strategy"],
                "winner_support_budget": winner["support_budget"],
                "winner_compute_cost_bits": winner["compute_cost_bits"],
                "winner_total_objective_bits": winner["total_objective_bits"],
                "winner_certified_upper_relative_residual": winner[
                    "certified_upper_relative_residual"
                ],
                "encoding_best": [
                    {
                        "encoding": enc,
                        "support_strategy": row["support_strategy"],
                        "support_budget": row["support_budget"],
                        "compute_cost_bits": row["compute_cost_bits"],
                        "total_objective_bits": row["total_objective_bits"],
                        "certified_upper_relative_residual": row[
                            "certified_upper_relative_residual"
                        ],
                    }
                    for enc, row in sorted(by_encoding.items())
                ],
            }
        )
    return table


def hetero_encoding_for_scenario(scenario: str) -> str:
    if scenario == "base":
        return "heterogeneous_compute_reusable"
    if scenario.startswith("drop_"):
        return f"heterogeneous_compute_reusable_{scenario}"
    raise ValueError(f"unknown scenario {scenario}")


def same_support_summary(rows: list[dict[str, Any]], scenario: str) -> dict[str, Any]:
    hetero_encoding = hetero_encoding_for_scenario(scenario)
    groups: dict[tuple[str, str, int, int, str, int], dict[str, dict[str, Any]]] = {}
    for row in rows:
        if row["scenario"] != scenario:
            continue
        key = (
            row["grid"],
            row["dictionary"],
            row["N"],
            row["support_budget"],
            row["support_strategy"],
            row["coefficient_precision_bits"],
        )
        groups.setdefault(key, {})[row["encoding"]] = row

    comparable = 0
    wins = 0
    savings = []
    for encs in groups.values():
        baselines = [row for enc, row in encs.items() if enc in BASELINE_ENCODINGS]
        hetero = encs.get(hetero_encoding)
        if not baselines or hetero is None:
            continue
        comparable += 1
        best_baseline = min(baselines, key=lambda row: row["total_objective_bits"])
        delta = best_baseline["total_objective_bits"] - hetero["total_objective_bits"]
        savings.append(delta)
        if delta > ORDER_TOL:
            wins += 1
    return {
        "same_support_comparable_count": comparable,
        "heterogeneous_same_support_win_count": wins,
        "heterogeneous_same_support_win_fraction": wins / comparable if comparable else None,
        "median_objective_savings_vs_best_baseline": (
            float(np.median(savings)) if savings else None
        ),
    }


def scenario_tolerance_summary(tolerances: list[dict[str, Any]], scenario: str) -> dict[str, Any]:
    hetero_encoding = hetero_encoding_for_scenario(scenario)
    scenario_rows = [row for row in tolerances if row["scenario"] == scenario]
    winner_counts: dict[str, int] = {}
    hetero_by_grid = {spec["name"]: 0 for spec in GRID_SPECS}
    hetero_by_n = {str(n): 0 for n in DEFAULT_N_VALUES}
    hetero_by_coeff_bits = {str(bits): 0 for bits in COEFFICIENT_BITS}
    for row in scenario_rows:
        winner = row["winner_encoding"]
        winner_counts[winner] = winner_counts.get(winner, 0) + 1
        if winner == hetero_encoding:
            hetero_by_grid[row["grid"]] = hetero_by_grid.get(row["grid"], 0) + 1
            n_key = str(row["N"])
            hetero_by_n[n_key] = hetero_by_n.get(n_key, 0) + 1
            bit_key = str(row["coefficient_precision_bits"])
            hetero_by_coeff_bits[bit_key] = hetero_by_coeff_bits.get(bit_key, 0) + 1
    hetero_wins = winner_counts.get(hetero_encoding, 0)
    return {
        "scenario": scenario,
        "heterogeneous_encoding": hetero_encoding,
        "tolerance_row_count": len(scenario_rows),
        "winner_counts": winner_counts,
        "heterogeneous_tolerance_win_count": hetero_wins,
        "heterogeneous_tolerance_win_fraction": (
            hetero_wins / len(scenario_rows) if scenario_rows else None
        ),
        "heterogeneous_wins_by_grid": hetero_by_grid,
        "heterogeneous_wins_by_N": hetero_by_n,
        "heterogeneous_wins_by_coefficient_bits": hetero_by_coeff_bits,
    }


def ablation_passes(summary: dict[str, Any]) -> bool:
    tolerance_fraction = summary.get("heterogeneous_tolerance_win_fraction") or 0.0
    same_support_median = summary.get("median_objective_savings_vs_best_baseline")
    grid_wins = summary.get("heterogeneous_wins_by_grid", {})
    n_wins = summary.get("heterogeneous_wins_by_N", {})
    return (
        tolerance_fraction >= 0.50
        and same_support_median is not None
        and same_support_median > 0
        and sum(1 for value in grid_wins.values() if value > 0) >= 2
        and sum(1 for value in n_wins.values() if value > 0) >= 3
    )


def base_replicates(summary: dict[str, Any]) -> bool:
    tolerance_fraction = summary.get("heterogeneous_tolerance_win_fraction") or 0.0
    same_support_median = summary.get("median_objective_savings_vs_best_baseline")
    grid_wins = summary.get("heterogeneous_wins_by_grid", {})
    n_wins = summary.get("heterogeneous_wins_by_N", {})
    large_n_win = n_wins.get("64", 0) > 0 or n_wins.get("96", 0) > 0
    return (
        tolerance_fraction >= 0.50
        and all(grid_wins.get(spec["name"], 0) > 0 for spec in GRID_SPECS)
        and sum(1 for value in n_wins.values() if value > 0) >= 3
        and large_n_win
        and same_support_median is not None
        and same_support_median > 0
    )


def classify(scenario_summaries: dict[str, dict[str, Any]]) -> tuple[str, list[str]]:
    reasons: list[str] = []
    base = scenario_summaries["base"]
    base_ok = base_replicates(base)
    ablation_ok = {
        scenario: ablation_passes(summary)
        for scenario, summary in scenario_summaries.items()
        if scenario.startswith("drop_")
    }
    ablation_pass_count = sum(1 for ok in ablation_ok.values() if ok)
    reasons.append(
        "base_replicates="
        f"{base_ok}; ablation_pass_count={ablation_pass_count}/{len(ablation_ok)}"
    )
    if base_ok and ablation_pass_count >= 3:
        return "REPLICATED_HETEROGENEOUS_COMPUTE_SIGNAL", reasons
    if base_ok:
        failed = sorted(scenario for scenario, ok in ablation_ok.items() if not ok)
        reasons.append(f"failed_ablation_scenarios={failed}")
        return "CHANNEL_SENSITIVE_HETEROGENEOUS_COMPUTE_SIGNAL", reasons
    if (base.get("heterogeneous_tolerance_win_count") or 0) > 0:
        return "MIXED_HETEROGENEOUS_COMPUTE_SIGNAL", reasons
    if any(
        (summary.get("heterogeneous_tolerance_win_count") or 0) > 0
        for summary in scenario_summaries.values()
    ):
        return "MIXED_HETEROGENEOUS_COMPUTE_SIGNAL", reasons
    return "FAILED_REPLICATION", reasons


def summarize(rows: list[dict[str, Any]], tolerances: list[dict[str, Any]]) -> dict[str, Any]:
    scenario_names = ["base"] + [f"drop_{channel}" for channel in ABLATION_CHANNELS]
    scenario_summaries: dict[str, dict[str, Any]] = {}
    for scenario in scenario_names:
        combined = scenario_tolerance_summary(tolerances, scenario)
        combined.update(same_support_summary(rows, scenario))
        scenario_summaries[scenario] = combined
    status, reasons = classify(scenario_summaries)
    return {
        "status": status,
        "classification_reasons": reasons,
        "row_count": len(rows),
        "tolerance_row_count": len(tolerances),
        "scenario_summaries": scenario_summaries,
    }


def run(n_values: Iterable[int]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for grid_spec in GRID_SPECS:
        x, w = build_grid(grid_spec["name"], grid_spec["nodes"])
        for n in n_values:
            for spec in DICTIONARIES:
                theta = theta_values(spec, n)
                A, y = build_design(theta, x, w)
                for m in support_budgets(n):
                    add_rows_for_support(
                        rows,
                        spec,
                        n,
                        grid_spec["name"],
                        A,
                        y,
                        "integer_prefix",
                        integer_prefix(m),
                    )
                    comp_cols = composite_prefix(n, m)
                    if len(comp_cols) >= 2:
                        add_rows_for_support(
                            rows,
                            spec,
                            n,
                            grid_spec["name"],
                            A,
                            y,
                            "composite_prefix",
                            comp_cols,
                        )
                    add_rows_for_support(
                        rows,
                        spec,
                        n,
                        grid_spec["name"],
                        A,
                        y,
                        "factor_cost_ordered",
                        factor_cost_ordered(n, m),
                    )
    tolerances = tolerance_table(rows)
    return {
        "experiment_id": EXPERIMENT_ID,
        "date": date.today().isoformat(),
        "claim_ceiling": CLAIM_CEILING,
        "model_notes": MODEL_NOTES,
        "source_experiment": SOURCE_EXPERIMENT,
        "cost_channels": {
            "direct_index": "direct finite lookup of selected dictionary indices",
            "prime_lookup": "lookup of distinct prime tokens used by selected indices",
            "multiply": f"{MULTIPLY_WEIGHT} bit-equivalent units per factor composition",
            "exponent": f"{EXPONENT_WEIGHT} bit-equivalent units per repeated-factor operation",
            "factor_tree_depth": f"{FACTOR_TREE_DEPTH_WEIGHT} bit-equivalent unit per max factor-tree level",
            "coefficient_move": "coefficient precision sweep over 8, 16, and 24 bits",
            "shared_harness": "optional setup cost for declaring the reusable prime harness",
        },
        "objective": "total_objective_bits = compute_cost_bits - certified_residual_information_bits",
        "grids": GRID_SPECS,
        "n_values": list(n_values),
        "support_budgets": SUPPORT_BUDGETS,
        "coefficient_bits": COEFFICIENT_BITS,
        "residual_tolerances": RESIDUAL_TOLERANCES,
        "summary": summarize(rows, tolerances),
        "tolerance_table": tolerances,
        "rows": rows,
    }


def markdown_report(results: dict[str, Any]) -> str:
    summary = results["summary"]
    scenarios = summary["scenario_summaries"]
    base = scenarios["base"]
    lines = [
        f"# {EXPERIMENT_ID} Report",
        "",
        "## Verdict",
        "",
        f"- Status: `{summary['status']}`",
        f"- Rows: `{summary['row_count']}`",
        f"- Tolerance rows: `{summary['tolerance_row_count']}`",
        f"- Base heterogeneous tolerance wins: `{base['heterogeneous_tolerance_win_count']}` / `{base['tolerance_row_count']}`",
        f"- Base same-support wins: `{base['heterogeneous_same_support_win_count']}` / `{base['same_support_comparable_count']}`",
        f"- Base median objective savings vs best baseline: `{base['median_objective_savings_vs_best_baseline']}`",
        "",
        "Classification reasons:",
        "",
        "```json",
        json.dumps(summary["classification_reasons"], indent=2, sort_keys=True),
        "```",
        "",
        "## Meaning",
        "",
        "This replication tries to break the positive heterogeneous-compute",
        "signal by adding larger dictionary sizes, a shifted grid, coefficient",
        "precision sensitivity, and one-channel-at-a-time ablations.",
        "",
        "The result is still a finite operation-weighted diagnostic. It is not",
        "a theorem about primes, not a zeta formalization, and not a",
        "quantum-mechanical claim.",
        "",
        "## Claim Ceiling",
        "",
        results["claim_ceiling"],
        "",
        "## Base Scenario",
        "",
        "```json",
        json.dumps(base, indent=2, sort_keys=True),
        "```",
        "",
        "## Ablation Scenarios",
        "",
        "| scenario | tolerance wins | tolerance rows | same-support wins | comparable | median savings |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for scenario in sorted(name for name in scenarios if name.startswith("drop_")):
        row = scenarios[scenario]
        lines.append(
            f"| `{scenario}` | `{row['heterogeneous_tolerance_win_count']}` | "
            f"`{row['tolerance_row_count']}` | "
            f"`{row['heterogeneous_same_support_win_count']}` | "
            f"`{row['same_support_comparable_count']}` | "
            f"`{row['median_objective_savings_vs_best_baseline']}` |"
        )
    lines.extend(
        [
            "",
            "## Base Tolerance Winners",
            "",
            "| grid | dictionary | N | coeff bits | tolerance | winner | support | objective | residual |",
            "|---|---|---:|---:|---:|---|---|---:|---:|",
        ]
    )
    base_rows = [row for row in results["tolerance_table"] if row["scenario"] == "base"]
    for row in base_rows[:80]:
        lines.append(
            "| {grid} | {dictionary} | {N} | {coefficient_precision_bits} | "
            "{tolerance:.2f} | {winner_encoding} | {winner_support_strategy} | "
            "{winner_total_objective_bits:.3f} | "
            "{winner_certified_upper_relative_residual:.6g} |".format(**row)
        )
    remaining = len(base_rows) - 80
    if remaining > 0:
        lines.append(f"| ... | ... | ... | ... | ... | {remaining} more rows | ... | ... | ... |")
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "A replicated result here means only that the finite reusable",
            "operation-weighted cost model survived this robustness gate against",
            "the listed finite baselines and ablations. It does not establish",
            "an infinite-dimensional theorem, dictionary invariance, or a",
            "result about zeta zeros.",
            "",
            "## Artifact Boundary",
            "",
            "Generated artifacts:",
            "",
            f"- `{EXPERIMENT_ID}_RESULTS.json`",
            f"- `{EXPERIMENT_ID}_REPORT.md`",
            f"- `{EXPERIMENT_ID}_RESULTS.sha256`",
            "",
            "No D1, scorecard, public page, git staging, commit, Zenodo, arXiv, or publisher surface was updated.",
            "",
        ]
    )
    return "\n".join(lines)


def write_outputs(results: dict[str, Any], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    json_path = outdir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = outdir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = outdir / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for path in (json_path, report_path, sha_path):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite existing artifact: {path}")
    json_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    report_path.write_text(markdown_report(results))
    sha_path.write_text(f"{sha256_file(json_path)}  {json_path.name}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Output directory for immutable experiment artifacts.",
    )
    parser.add_argument(
        "--n-values",
        type=int,
        nargs="*",
        default=list(DEFAULT_N_VALUES),
        help="Dictionary sizes to test.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results = run(args.n_values)
    write_outputs(results, args.outdir)
    print(
        json.dumps(
            {
                "experiment_id": EXPERIMENT_ID,
                "status": results["summary"]["status"],
                "output_dir": str(args.outdir),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
