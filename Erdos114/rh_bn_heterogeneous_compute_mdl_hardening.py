#!/usr/bin/env python3
"""
EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-HARDENING-20260507-01

Hardening gate for the finite heterogeneous-compute MDL diagnostic.

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
from dataclasses import dataclass, asdict
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
    factor_harness_bits,
    factorization,
    factorized_address_bits,
    flat_address_bits,
    int_bits_at_most,
    primes_up_to,
)


EXPERIMENT_ID = "EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-HARDENING-20260507-01"
CLAIM_CEILING = (
    "INTERNAL / METHOD-SHAPING ONLY: finite heterogeneous-compute "
    "Beurling-Nyman MDL hardening gate; no RH claim, no zeta formalization, "
    "and no quantum claim."
)
SOURCE_EXPERIMENTS = (
    "EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01",
    "EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-REPLICATION-20260507-01",
)
MODEL_NOTES = (
    "HETEROGENEOUS_COMPUTE_MDL_MODEL_2026-05-06.md",
    "HETEROGENEOUS_COMPUTE_MDL_EXPERIMENTAL_NOTE_2026-05-06.md",
)
DEFAULT_N_VALUES = (16, 32, 48, 64, 96)
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
RIDGE_ALPHAS = (1e-6, 1e-3)
OMP_ALPHA = 1e-6
STAGE_HEADER_BITS = 8
ORDER_TOL = 1e-9


@dataclass(frozen=True)
class WeightProfile:
    name: str
    prime_lookup_scale: int
    multiply_weight: int
    exponent_weight: int
    factor_tree_depth_weight: int


WEIGHT_PROFILES = (
    WeightProfile("base", 1, 2, 2, 1),
    WeightProfile("all_unit", 1, 1, 1, 1),
    WeightProfile("low_ops", 1, 1, 1, 0),
    WeightProfile("high_ops", 1, 4, 4, 2),
    WeightProfile("prime_heavy", 3, 2, 2, 1),
    WeightProfile("operator_heavy", 1, 6, 6, 3),
    WeightProfile("free_arithmetic", 0, 0, 0, 0),
)
BASELINE_ENCODINGS = {"flat_index", "factorized_reusable_harness"}


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
        "information_bits_certified": info,
    }


def ridge_coefficients(A: np.ndarray, y: np.ndarray, alpha: float) -> np.ndarray:
    gram = A.T @ A
    rhs = A.T @ y
    return np.linalg.solve(gram + alpha * np.eye(gram.shape[0]), rhs)


def ridge_top_support(A: np.ndarray, y: np.ndarray, alpha: float, m: int) -> list[int]:
    coeffs = ridge_coefficients(A, y, alpha)
    order = np.argsort(-np.abs(coeffs), kind="mergesort")[:m]
    return sorted(int(idx) + 1 for idx in order)


def ridge_refit_residual(A: np.ndarray, y: np.ndarray, cols_zero: list[int], alpha: float) -> np.ndarray:
    if not cols_zero:
        return y.copy()
    sub = A[:, cols_zero]
    gram = sub.T @ sub
    rhs = sub.T @ y
    coeffs = np.linalg.solve(gram + alpha * np.eye(gram.shape[0]), rhs)
    return y - sub @ coeffs


def omp_order(A: np.ndarray, y: np.ndarray, alpha: float) -> list[int]:
    selected: list[int] = []
    remaining = set(range(A.shape[1]))
    residual = y.copy()
    while remaining:
        correlations = A.T @ residual
        best = max(remaining, key=lambda idx: (abs(float(correlations[idx])), -idx))
        selected.append(best)
        remaining.remove(best)
        residual = ridge_refit_residual(A, y, selected, alpha)
    return [idx + 1 for idx in selected]


def support_budgets(n: int) -> list[int]:
    values = sorted({m for m in SUPPORT_BUDGETS if m <= n} | {n})
    return [m for m in values if m > 0]


def regularized_supports(A: np.ndarray, y: np.ndarray, n: int) -> dict[str, dict[int, list[int]]]:
    supports: dict[str, dict[int, list[int]]] = {}
    for alpha in RIDGE_ALPHAS:
        strategy = f"ridge_top_alpha_{alpha:g}"
        supports[strategy] = {m: ridge_top_support(A, y, alpha, m) for m in support_budgets(n)}
    order = omp_order(A, y, OMP_ALPHA)
    strategy = f"omp_ridge_alpha_{OMP_ALPHA:g}"
    supports[strategy] = {m: sorted(order[:m]) for m in support_budgets(n)}
    return supports


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
    weight_profile: str,
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
        "weight_profile": weight_profile,
        "encoding": encoding,
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
    if encoding == "flat_index":
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
    profile: WeightProfile,
    include_harness: bool,
) -> dict[str, int]:
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
    return {
        "direct_index": 0,
        "prime_lookup": len(distinct_primes) * prime_lookup_unit * profile.prime_lookup_scale,
        "multiply": multiply_ops * profile.multiply_weight,
        "exponent": exponent_ops * profile.exponent_weight,
        "factor_tree_depth": max_depth * profile.factor_tree_depth_weight,
        "coefficient_move": active * coeff_bits,
        "shared_harness": factor_harness_bits(n) if include_harness else 0,
    }


def cost_from_channels(spec: DictionarySpec, channels: dict[str, int]) -> int:
    return spec.header_bits + STAGE_HEADER_BITS + sum(channels.values())


def add_baseline_rows(
    rows: list[dict[str, Any]],
    spec: DictionarySpec,
    n: int,
    grid: str,
    support_strategy: str,
    profile: WeightProfile,
    cols: list[int],
    fit: dict[str, Any],
) -> None:
    for encoding, address_bits in (
        ("flat_index", flat_address_bits(cols, n)),
        ("factorized_reusable_harness", factorized_address_bits(cols, n, include_harness=False)),
    ):
        row = base_row(
            spec,
            n,
            grid,
            support_strategy,
            profile.name,
            encoding,
            cols,
            fit,
            address_bits,
            False,
        )
        rows.append(attach_objective(row, int(row["total_description_bits"]), baseline_channels(encoding, address_bits, row)))


def add_heterogeneous_rows(
    rows: list[dict[str, Any]],
    spec: DictionarySpec,
    n: int,
    grid: str,
    support_strategy: str,
    profile: WeightProfile,
    cols: list[int],
    fit: dict[str, Any],
) -> None:
    for encoding, include_harness in (
        ("heterogeneous_compute_reusable", False),
        ("heterogeneous_compute_with_harness", True),
    ):
        channels = operation_channels(
            cols,
            n,
            fit["active_count"],
            int(fit["coefficient_precision_bits"]),
            profile,
            include_harness,
        )
        compute_cost = cost_from_channels(spec, channels)
        address_bits = sum(channels.values()) - channels["coefficient_move"]
        row = base_row(
            spec,
            n,
            grid,
            support_strategy,
            profile.name,
            encoding,
            cols,
            fit,
            address_bits,
            include_harness,
        )
        row["total_description_bits"] = compute_cost
        rows.append(attach_objective(row, compute_cost, channels))


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
        for profile in WEIGHT_PROFILES:
            add_baseline_rows(rows, spec, n, grid, support_strategy, profile, cols, fit)
            add_heterogeneous_rows(rows, spec, n, grid, support_strategy, profile, cols, fit)


def tolerance_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, str, str, int, int, float], list[dict[str, Any]]] = {}
    for row in rows:
        for tol in RESIDUAL_TOLERANCES:
            key = (
                row["weight_profile"],
                row["support_strategy"],
                row["grid"],
                row["dictionary"],
                row["N"],
                row["coefficient_precision_bits"],
                tol,
            )
            groups.setdefault(key, []).append(row)

    table = []
    for (profile, support_strategy, grid, dictionary, n, coeff_bits, tol), group_rows in sorted(groups.items()):
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
                "weight_profile": profile,
                "support_strategy": support_strategy,
                "grid": grid,
                "dictionary": dictionary,
                "N": n,
                "coefficient_precision_bits": coeff_bits,
                "tolerance": tol,
                "winner_encoding": winner["encoding"],
                "winner_support_budget": winner["support_budget"],
                "winner_compute_cost_bits": winner["compute_cost_bits"],
                "winner_total_objective_bits": winner["total_objective_bits"],
                "winner_certified_upper_relative_residual": winner[
                    "certified_upper_relative_residual"
                ],
                "encoding_best": [
                    {
                        "encoding": enc,
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


def same_support_summary(rows: list[dict[str, Any]], profile: str) -> dict[str, Any]:
    groups: dict[tuple[str, str, int, int, str, int], dict[str, dict[str, Any]]] = {}
    for row in rows:
        if row["weight_profile"] != profile:
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
    reusable_wins = 0
    with_harness_wins = 0
    reusable_savings = []
    harness_savings = []
    for encs in groups.values():
        baselines = [row for enc, row in encs.items() if enc in BASELINE_ENCODINGS]
        reusable = encs.get("heterogeneous_compute_reusable")
        with_harness = encs.get("heterogeneous_compute_with_harness")
        if not baselines or reusable is None or with_harness is None:
            continue
        comparable += 1
        best_baseline = min(baselines, key=lambda row: row["total_objective_bits"])
        reusable_delta = best_baseline["total_objective_bits"] - reusable["total_objective_bits"]
        harness_delta = best_baseline["total_objective_bits"] - with_harness["total_objective_bits"]
        reusable_savings.append(reusable_delta)
        harness_savings.append(harness_delta)
        if reusable_delta > ORDER_TOL:
            reusable_wins += 1
        if harness_delta > ORDER_TOL:
            with_harness_wins += 1
    return {
        "same_support_comparable_count": comparable,
        "heterogeneous_reusable_same_support_win_count": reusable_wins,
        "heterogeneous_reusable_same_support_win_fraction": reusable_wins / comparable if comparable else None,
        "heterogeneous_with_harness_same_support_win_count": with_harness_wins,
        "heterogeneous_with_harness_same_support_win_fraction": with_harness_wins / comparable if comparable else None,
        "median_reusable_objective_savings_vs_best_baseline": (
            float(np.median(reusable_savings)) if reusable_savings else None
        ),
        "median_with_harness_objective_savings_vs_best_baseline": (
            float(np.median(harness_savings)) if harness_savings else None
        ),
    }


def profile_tolerance_summary(tolerances: list[dict[str, Any]], profile: str) -> dict[str, Any]:
    profile_rows = [row for row in tolerances if row["weight_profile"] == profile]
    winner_counts: dict[str, int] = {}
    wins_by_grid = {spec["name"]: 0 for spec in GRID_SPECS}
    wins_by_n = {str(n): 0 for n in DEFAULT_N_VALUES}
    wins_by_support = {}
    wins_by_coeff_bits = {str(bits): 0 for bits in COEFFICIENT_BITS}
    for row in profile_rows:
        winner = row["winner_encoding"]
        winner_counts[winner] = winner_counts.get(winner, 0) + 1
        if winner == "heterogeneous_compute_reusable":
            wins_by_grid[row["grid"]] = wins_by_grid.get(row["grid"], 0) + 1
            n_key = str(row["N"])
            wins_by_n[n_key] = wins_by_n.get(n_key, 0) + 1
            wins_by_support[row["support_strategy"]] = wins_by_support.get(row["support_strategy"], 0) + 1
            bit_key = str(row["coefficient_precision_bits"])
            wins_by_coeff_bits[bit_key] = wins_by_coeff_bits.get(bit_key, 0) + 1
    wins = winner_counts.get("heterogeneous_compute_reusable", 0)
    return {
        "weight_profile": profile,
        "tolerance_row_count": len(profile_rows),
        "winner_counts": winner_counts,
        "heterogeneous_reusable_tolerance_win_count": wins,
        "heterogeneous_reusable_tolerance_win_fraction": (
            wins / len(profile_rows) if profile_rows else None
        ),
        "heterogeneous_reusable_wins_by_grid": wins_by_grid,
        "heterogeneous_reusable_wins_by_N": wins_by_n,
        "heterogeneous_reusable_wins_by_support_strategy": wins_by_support,
        "heterogeneous_reusable_wins_by_coefficient_bits": wins_by_coeff_bits,
    }


def profile_passes(summary: dict[str, Any]) -> bool:
    tolerance_fraction = summary.get("heterogeneous_reusable_tolerance_win_fraction") or 0.0
    same_support_median = summary.get("median_reusable_objective_savings_vs_best_baseline")
    grid_wins = summary.get("heterogeneous_reusable_wins_by_grid", {})
    n_wins = summary.get("heterogeneous_reusable_wins_by_N", {})
    support_wins = summary.get("heterogeneous_reusable_wins_by_support_strategy", {})
    large_n_win = n_wins.get("64", 0) > 0 or n_wins.get("96", 0) > 0
    return (
        tolerance_fraction >= 0.50
        and all(grid_wins.get(spec["name"], 0) > 0 for spec in GRID_SPECS)
        and sum(1 for value in n_wins.values() if value > 0) >= 3
        and large_n_win
        and len(support_wins) >= 2
        and same_support_median is not None
        and same_support_median > 0
    )


def independent_factorization(n: int) -> dict[int, int]:
    factors: dict[int, int] = {}
    value = n
    candidate = 2
    while candidate <= value:
        while value % candidate == 0:
            factors[candidate] = factors.get(candidate, 0) + 1
            value //= candidate
        candidate += 1
    return factors


def independent_operation_channels(
    indices: list[int],
    n: int,
    active: int,
    coeff_bits: int,
    profile: WeightProfile,
    include_harness: bool,
) -> dict[str, int]:
    distinct_primes: set[int] = set()
    multiply_ops = 0
    exponent_ops = 0
    max_depth = 0
    for j in indices:
        factors = independent_factorization(j)
        if not factors:
            continue
        distinct_primes.update(factors)
        multiplicity = sum(factors.values())
        multiply_ops += max(0, multiplicity - 1)
        exponent_ops += sum(max(0, exp - 1) for exp in factors.values())
        max_depth = max(max_depth, multiplicity)
    prime_lookup_unit = int_bits_at_most(len(primes_up_to(n)))
    return {
        "direct_index": 0,
        "prime_lookup": len(distinct_primes) * prime_lookup_unit * profile.prime_lookup_scale,
        "multiply": multiply_ops * profile.multiply_weight,
        "exponent": exponent_ops * profile.exponent_weight,
        "factor_tree_depth": max_depth * profile.factor_tree_depth_weight,
        "coefficient_move": active * coeff_bits,
        "shared_harness": factor_harness_bits(n) if include_harness else 0,
    }


def audit_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    errors: list[str] = []
    profile_map = {profile.name: profile for profile in WEIGHT_PROFILES}
    for idx, row in enumerate(rows):
        cols = row["columns"]
        n = row["N"]
        if cols != sorted(set(cols)):
            errors.append(f"row {idx}: columns are not sorted unique")
        if any(col < 1 or col > n for col in cols):
            errors.append(f"row {idx}: column outside 1..N")
        residual = row["certified_upper_relative_residual"]
        if residual > 0 and row["certified_residual_information_bits"] is not None:
            recomputed_info = -math.log2(residual)
            if abs(recomputed_info - row["certified_residual_information_bits"]) > 1e-9:
                errors.append(f"row {idx}: information bits mismatch")
        if row["encoding"].startswith("heterogeneous_compute"):
            profile = profile_map[row["weight_profile"]]
            expected = independent_operation_channels(
                cols,
                n,
                row["active_count"],
                row["coefficient_precision_bits"],
                profile,
                row["include_prime_harness_setup_cost"],
            )
            if expected != row["heterogeneous_cost_channels"]:
                errors.append(f"row {idx}: channel mismatch")
            expected_cost = row["dictionary_header_bits"] + row["stage_header_bits"] + sum(expected.values())
            if expected_cost != row["compute_cost_bits"]:
                errors.append(f"row {idx}: compute cost mismatch")
    return {
        "audit_status": "PASS" if not errors else "FAIL",
        "checked_rows": len(rows),
        "error_count": len(errors),
        "errors": errors[:20],
    }


def classify(profile_summaries: dict[str, dict[str, Any]], audit: dict[str, Any]) -> tuple[str, list[str]]:
    reasons: list[str] = []
    profile_pass = {name: profile_passes(summary) for name, summary in profile_summaries.items()}
    pass_count = sum(1 for ok in profile_pass.values() if ok)
    base_ok = profile_pass.get("base", False)
    reasons.append(f"base_passes={base_ok}; profile_pass_count={pass_count}/{len(profile_pass)}")
    reasons.append(f"audit_status={audit['audit_status']}; audit_error_count={audit['error_count']}")
    if audit["audit_status"] != "PASS":
        return "BLOCKED_BY_IMPLEMENTATION_AUDIT", reasons
    if base_ok and pass_count >= 4:
        return "HARDENED_HETEROGENEOUS_COMPUTE_SIGNAL", reasons
    if base_ok:
        failed = sorted(name for name, ok in profile_pass.items() if not ok)
        reasons.append(f"failed_weight_profiles={failed}")
        return "WEIGHT_SENSITIVE_HETEROGENEOUS_COMPUTE_SIGNAL", reasons
    if any((summary.get("heterogeneous_reusable_tolerance_win_count") or 0) > 0 for summary in profile_summaries.values()):
        return "MIXED_HETEROGENEOUS_COMPUTE_SIGNAL", reasons
    return "FAILED_HARDENING", reasons


def summarize(rows: list[dict[str, Any]], tolerances: list[dict[str, Any]], audit: dict[str, Any]) -> dict[str, Any]:
    profile_summaries: dict[str, dict[str, Any]] = {}
    for profile in WEIGHT_PROFILES:
        combined = profile_tolerance_summary(tolerances, profile.name)
        combined.update(same_support_summary(rows, profile.name))
        profile_summaries[profile.name] = combined
    status, reasons = classify(profile_summaries, audit)
    return {
        "status": status,
        "classification_reasons": reasons,
        "row_count": len(rows),
        "tolerance_row_count": len(tolerances),
        "profile_summaries": profile_summaries,
        "implementation_audit": audit,
    }


def run(n_values: Iterable[int]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for grid_spec in GRID_SPECS:
        x, w = build_grid(grid_spec["name"], grid_spec["nodes"])
        for n in n_values:
            for spec in DICTIONARIES:
                theta = theta_values(spec, n)
                A, y = build_design(theta, x, w)
                supports_by_strategy = regularized_supports(A, y, n)
                for support_strategy, supports_by_budget in supports_by_strategy.items():
                    for cols in supports_by_budget.values():
                        add_rows_for_support(rows, spec, n, grid_spec["name"], A, y, support_strategy, cols)
    tolerances = tolerance_table(rows)
    audit = audit_rows(rows)
    return {
        "experiment_id": EXPERIMENT_ID,
        "date": date.today().isoformat(),
        "claim_ceiling": CLAIM_CEILING,
        "source_experiments": SOURCE_EXPERIMENTS,
        "model_notes": MODEL_NOTES,
        "objective": "total_objective_bits = compute_cost_bits - certified_residual_information_bits",
        "grids": GRID_SPECS,
        "n_values": list(n_values),
        "support_budgets": SUPPORT_BUDGETS,
        "coefficient_bits": COEFFICIENT_BITS,
        "ridge_alphas": RIDGE_ALPHAS,
        "omp_alpha": OMP_ALPHA,
        "weight_profiles": [asdict(profile) for profile in WEIGHT_PROFILES],
        "residual_tolerances": RESIDUAL_TOLERANCES,
        "summary": summarize(rows, tolerances, audit),
        "tolerance_table": tolerances,
        "rows": rows,
    }


def markdown_report(results: dict[str, Any]) -> str:
    summary = results["summary"]
    profiles = summary["profile_summaries"]
    base = profiles["base"]
    lines = [
        f"# {EXPERIMENT_ID} Report",
        "",
        "## Verdict",
        "",
        f"- Status: `{summary['status']}`",
        f"- Rows: `{summary['row_count']}`",
        f"- Tolerance rows: `{summary['tolerance_row_count']}`",
        f"- Base reusable tolerance wins: `{base['heterogeneous_reusable_tolerance_win_count']}` / `{base['tolerance_row_count']}`",
        f"- Base same-support wins: `{base['heterogeneous_reusable_same_support_win_count']}` / `{base['same_support_comparable_count']}`",
        f"- Base median reusable savings vs best baseline: `{base['median_reusable_objective_savings_vs_best_baseline']}`",
        f"- Implementation audit: `{summary['implementation_audit']['audit_status']}`",
        "",
        "Classification reasons:",
        "",
        "```json",
        json.dumps(summary["classification_reasons"], indent=2, sort_keys=True),
        "```",
        "",
        "## Meaning",
        "",
        "This hardening gate tests whether the heterogeneous-compute signal",
        "survives regularized active-set supports, operation-weight sensitivity,",
        "and an independent implementation audit.",
        "",
        "The result is still a finite operation-weighted diagnostic. It is not",
        "a theorem about primes, not a zeta formalization, and not a quantum",
        "claim.",
        "",
        "## Claim Ceiling",
        "",
        results["claim_ceiling"],
        "",
        "## Weight Profiles",
        "",
        "| profile | tolerance wins | tolerance rows | same-support wins | comparable | median savings |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for profile in sorted(profiles):
        row = profiles[profile]
        lines.append(
            f"| `{profile}` | `{row['heterogeneous_reusable_tolerance_win_count']}` | "
            f"`{row['tolerance_row_count']}` | "
            f"`{row['heterogeneous_reusable_same_support_win_count']}` | "
            f"`{row['same_support_comparable_count']}` | "
            f"`{row['median_reusable_objective_savings_vs_best_baseline']}` |"
        )
    lines.extend(
        [
            "",
            "## Implementation Audit",
            "",
            "```json",
            json.dumps(summary["implementation_audit"], indent=2, sort_keys=True),
            "```",
            "",
            "## Base Tolerance Winners",
            "",
            "| support | grid | dictionary | N | coeff bits | tolerance | winner | objective | residual |",
            "|---|---|---|---:|---:|---:|---|---:|---:|",
        ]
    )
    base_rows = [row for row in results["tolerance_table"] if row["weight_profile"] == "base"]
    for row in base_rows[:80]:
        lines.append(
            "| {support_strategy} | {grid} | {dictionary} | {N} | "
            "{coefficient_precision_bits} | {tolerance:.2f} | {winner_encoding} | "
            "{winner_total_objective_bits:.3f} | "
            "{winner_certified_upper_relative_residual:.6g} |".format(**row)
        )
    remaining = len(base_rows) - 80
    if remaining > 0:
        lines.append(f"| ... | ... | ... | ... | ... | ... | {remaining} more rows | ... | ... |")
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "A hardened result here means only that the finite reusable",
            "operation-weighted cost model survived this stricter finite gate.",
            "It does not establish an infinite-dimensional theorem, dictionary",
            "invariance, or a result about zeta zeros.",
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
