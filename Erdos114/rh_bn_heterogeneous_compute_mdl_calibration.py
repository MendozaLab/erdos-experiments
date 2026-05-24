#!/usr/bin/env python3
"""
EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-CALIBRATION-20260507-01

Calibration and negative-control gate for the finite heterogeneous-compute
MDL diagnostic.

Claim ceiling:
    INTERNAL / METHOD-SHAPING ONLY.
    This is a finite operation-weighted Beurling-Nyman MDL calibration check.
    It is not an RH claim, not a zeta formalization, and not a
    quantum-mechanical claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from rh_beurling_nyman_mdl_probe import DICTIONARIES, DictionarySpec, theta_values
from rh_bn_prime_harness_mdl import (
    RESIDUAL_TOLERANCES,
    factor_harness_bits,
    factorization,
    factorized_address_bits,
    flat_address_bits,
    int_bits_at_most,
    primes_up_to,
)
from rh_bn_heterogeneous_compute_mdl_hardening import (
    COEFFICIENT_BITS,
    GRID_SPECS,
    OMP_ALPHA,
    RIDGE_ALPHAS,
    STAGE_HEADER_BITS,
    SUPPORT_BUDGETS,
    build_design,
    build_grid,
    fit_columns_base,
    fit_metrics,
    regularized_supports,
    support_budgets,
)


EXPERIMENT_ID = "EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-CALIBRATION-20260507-01"
CLAIM_CEILING = (
    "INTERNAL / METHOD-SHAPING ONLY: finite heterogeneous-compute "
    "Beurling-Nyman MDL calibration gate; no RH claim, no zeta "
    "formalization, and no quantum claim."
)
SOURCE_EXPERIMENTS = (
    "EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01",
    "EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-REPLICATION-20260507-01",
    "EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-HARDENING-20260507-01",
)
MODEL_NOTES = (
    "HETEROGENEOUS_COMPUTE_MDL_MODEL_2026-05-06.md",
    "HETEROGENEOUS_COMPUTE_MDL_EXPERIMENTAL_NOTE_2026-05-06.md",
)
DEFAULT_N_VALUES = (16, 32, 48, 64, 96)
BASELINE_ENCODINGS = {"flat_index", "factorized_reusable_harness"}
CALIBRATED_PROFILES = {"circuit_unit", "prefix_decode", "symbolic_runtime_proxy"}
STRESS_PROFILES = {"stress_prime_heavy", "stress_operator_heavy"}
NEGATIVE_CONTROLS = {"shuffled_indices", "random_support", "synthetic_nonfactorized"}
CONTROL_SIMILARITY_RATIO = 0.80
ORDER_TOL = 1e-9
RANDOM_SEED = 11420260507


@dataclass(frozen=True)
class WeightProfile:
    name: str
    prime_lookup_scale: int
    multiply_weight: int
    exponent_weight: int
    factor_tree_depth_weight: int
    calibration_source: str


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_seed(*parts: object) -> int:
    payload = "|".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little") % (2**32)


def median_runtime_ns(fn, repeats: int = 7) -> float:
    values = []
    for _ in range(repeats):
        start = time.perf_counter_ns()
        fn()
        values.append(time.perf_counter_ns() - start)
    return float(np.median(values))


def runtime_proxy_profile(max_n: int) -> tuple[WeightProfile, dict[str, Any]]:
    values = list(range(2, max_n + 1))
    factorizations = {j: factorization(j) for j in values}
    prime_ops = sum(max(1, len(factors)) for factors in factorizations.values())
    multiply_ops = sum(max(1, sum(factors.values()) - 1) for factors in factorizations.values())
    exponent_ops = sum(
        max(1, sum(max(0, exp - 1) for exp in factors.values()))
        for factors in factorizations.values()
    )
    depth_ops = len(values)

    def run_prime_lookup() -> None:
        total = 0
        for j in values:
            total += len(factorization(j))
        if total < 0:
            raise RuntimeError("unreachable")

    def run_multiply() -> None:
        total = 0
        for factors in factorizations.values():
            product = 1
            for prime, exp in factors.items():
                for _ in range(exp):
                    product *= prime
            total += product
        if total < 0:
            raise RuntimeError("unreachable")

    def run_exponent() -> None:
        total = 0
        for factors in factorizations.values():
            for prime, exp in factors.items():
                total += prime ** exp
        if total < 0:
            raise RuntimeError("unreachable")

    def run_depth() -> None:
        total = 0
        for factors in factorizations.values():
            total += sum(factors.values())
        if total < 0:
            raise RuntimeError("unreachable")

    timings = {
        "prime_lookup_ns_per_op": median_runtime_ns(run_prime_lookup) / prime_ops,
        "multiply_ns_per_op": median_runtime_ns(run_multiply) / multiply_ops,
        "exponent_ns_per_op": median_runtime_ns(run_exponent) / exponent_ops,
        "factor_tree_depth_ns_per_op": median_runtime_ns(run_depth) / depth_ops,
    }
    unit = min(value for value in timings.values() if value > 0)
    weights = {
        key.replace("_ns_per_op", ""): max(1, min(12, int(round(value / unit))))
        for key, value in timings.items()
    }
    profile = WeightProfile(
        "symbolic_runtime_proxy",
        weights["prime_lookup"],
        weights["multiply"],
        weights["exponent"],
        weights["factor_tree_depth"],
        "local median Python symbolic factorization/reconstruction timings",
    )
    return profile, {"max_n": max_n, "timings": timings, "weights": asdict(profile)}


def weight_profiles(max_n: int) -> tuple[tuple[WeightProfile, ...], dict[str, Any]]:
    runtime_profile, runtime_meta = runtime_proxy_profile(max_n)
    profiles = (
        WeightProfile("circuit_unit", 1, 1, 1, 1, "simple arithmetic circuit step count"),
        WeightProfile("prefix_decode", 2, 1, 2, 3, "finite prefix decode proxy"),
        runtime_profile,
        WeightProfile("stress_prime_heavy", 4, 2, 2, 1, "failure probe with expensive prime lookup"),
        WeightProfile("stress_operator_heavy", 1, 6, 6, 3, "failure probe with expensive arithmetic operations"),
    )
    return profiles, {"symbolic_runtime_proxy": runtime_meta}


def synthetic_design_like(A: np.ndarray, key: str) -> np.ndarray:
    rng = np.random.default_rng(stable_seed(RANDOM_SEED, "synthetic", key))
    raw = rng.normal(size=A.shape)
    raw_norms = np.linalg.norm(raw, axis=0)
    target_norms = np.linalg.norm(A, axis=0)
    raw_norms = np.where(raw_norms > 0, raw_norms, 1.0)
    return raw * (target_norms / raw_norms)[None, :]


def shuffled_cost_labels(cols: list[int], n: int, key: str) -> list[int]:
    rng = np.random.default_rng(stable_seed(RANDOM_SEED, "shuffle", key, n))
    perm = rng.permutation(np.arange(1, n + 1))
    return sorted(int(perm[col - 1]) for col in cols)


def random_supports(n: int, key: str) -> dict[str, dict[int, list[int]]]:
    rng = np.random.default_rng(stable_seed(RANDOM_SEED, "random_support", key, n))
    supports = {}
    by_budget = {}
    for m in support_budgets(n):
        by_budget[m] = sorted(int(value) for value in rng.choice(np.arange(1, n + 1), size=m, replace=False))
    supports["random_matched_support_budget"] = by_budget
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
    control_type: str,
    support_strategy: str,
    profile: WeightProfile,
    encoding: str,
    support_cols: list[int],
    cost_labels: list[int],
    fit: dict[str, Any],
    address_bits: int,
    include_harness: bool,
) -> dict[str, Any]:
    coeff_bits = int(fit["coefficient_precision_bits"])
    bits = total_bits(spec, address_bits, fit["active_count"], coeff_bits)
    return {
        "dictionary": spec.name,
        "theta_rule": spec.theta_rule,
        "control_type": control_type,
        "grid": grid,
        "N": n,
        "support_budget": len(support_cols),
        "support_strategy": support_strategy,
        "weight_profile": profile.name,
        "calibration_source": profile.calibration_source,
        "encoding": encoding,
        "include_prime_harness_setup_cost": include_harness,
        "columns": support_cols,
        "cost_labels": cost_labels,
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
    labels: list[int],
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
    for j in labels:
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


def add_rows_for_support(
    rows: list[dict[str, Any]],
    spec: DictionarySpec,
    n: int,
    grid: str,
    control_type: str,
    support_strategy: str,
    A: np.ndarray,
    y: np.ndarray,
    support_cols: list[int],
    cost_labels: list[int],
    profiles: tuple[WeightProfile, ...],
) -> None:
    support_cols = sorted(dict.fromkeys(support_cols))
    cost_labels = sorted(dict.fromkeys(cost_labels))
    if len(support_cols) != len(cost_labels):
        raise ValueError("support columns and cost labels must have the same unique size")
    fit_base = fit_columns_base(A, y, [j - 1 for j in support_cols])
    for coeff_bits in COEFFICIENT_BITS:
        fit = fit_metrics(fit_base, coeff_bits)
        for profile in profiles:
            for encoding, address_bits in (
                ("flat_index", flat_address_bits(cost_labels, n)),
                (
                    "factorized_reusable_harness",
                    factorized_address_bits(cost_labels, n, include_harness=False),
                ),
            ):
                row = base_row(
                    spec,
                    n,
                    grid,
                    control_type,
                    support_strategy,
                    profile,
                    encoding,
                    support_cols,
                    cost_labels,
                    fit,
                    address_bits,
                    False,
                )
                rows.append(attach_objective(row, int(row["total_description_bits"]), baseline_channels(encoding, address_bits, row)))
            for encoding, include_harness in (
                ("heterogeneous_compute_reusable", False),
                ("heterogeneous_compute_with_harness", True),
            ):
                channels = operation_channels(
                    cost_labels,
                    n,
                    fit["active_count"],
                    coeff_bits,
                    profile,
                    include_harness,
                )
                compute_cost = cost_from_channels(spec, channels)
                address_bits = sum(channels.values()) - channels["coefficient_move"]
                row = base_row(
                    spec,
                    n,
                    grid,
                    control_type,
                    support_strategy,
                    profile,
                    encoding,
                    support_cols,
                    cost_labels,
                    fit,
                    address_bits,
                    include_harness,
                )
                row["total_description_bits"] = compute_cost
                rows.append(attach_objective(row, compute_cost, channels))


def control_cases(A: np.ndarray, y: np.ndarray, n: int, key: str) -> list[tuple[str, np.ndarray, dict[str, dict[int, list[int]]], bool]]:
    real_supports = regularized_supports(A, y, n)
    synthetic = synthetic_design_like(A, key)
    synthetic_supports = regularized_supports(synthetic, y, n)
    return [
        ("real", A, real_supports, False),
        ("shuffled_indices", A, real_supports, True),
        ("random_support", A, random_supports(n, key), False),
        ("synthetic_nonfactorized", synthetic, synthetic_supports, False),
    ]


def tolerance_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, str, str, str, int, int, float], list[dict[str, Any]]] = {}
    for row in rows:
        for tol in RESIDUAL_TOLERANCES:
            key = (
                row["control_type"],
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
    for (control, profile, support_strategy, grid, dictionary, n, coeff_bits, tol), group_rows in sorted(groups.items()):
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
                "control_type": control,
                "weight_profile": profile,
                "support_strategy": support_strategy,
                "grid": grid,
                "dictionary": dictionary,
                "N": n,
                "coefficient_precision_bits": coeff_bits,
                "tolerance": tol,
                "winner_encoding": winner["encoding"],
                "winner_support_budget": winner["support_budget"],
                "winner_total_objective_bits": winner["total_objective_bits"],
                "winner_certified_upper_relative_residual": winner[
                    "certified_upper_relative_residual"
                ],
                "encoding_best": [
                    {
                        "encoding": enc,
                        "support_budget": row["support_budget"],
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


def same_support_summary(rows: list[dict[str, Any]], control: str, profile: str) -> dict[str, Any]:
    groups: dict[tuple[str, str, int, int, str, int], dict[str, dict[str, Any]]] = {}
    for row in rows:
        if row["control_type"] != control or row["weight_profile"] != profile:
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
    reusable_savings = []
    for encs in groups.values():
        baselines = [row for enc, row in encs.items() if enc in BASELINE_ENCODINGS]
        reusable = encs.get("heterogeneous_compute_reusable")
        if not baselines or reusable is None:
            continue
        comparable += 1
        best_baseline = min(baselines, key=lambda row: row["total_objective_bits"])
        delta = best_baseline["total_objective_bits"] - reusable["total_objective_bits"]
        reusable_savings.append(delta)
        if delta > ORDER_TOL:
            reusable_wins += 1
    return {
        "same_support_comparable_count": comparable,
        "heterogeneous_reusable_same_support_win_count": reusable_wins,
        "heterogeneous_reusable_same_support_win_fraction": reusable_wins / comparable if comparable else None,
        "median_reusable_objective_savings_vs_best_baseline": (
            float(np.median(reusable_savings)) if reusable_savings else None
        ),
    }


def control_profile_tolerance_summary(tolerances: list[dict[str, Any]], control: str, profile: str) -> dict[str, Any]:
    selected = [
        row
        for row in tolerances
        if row["control_type"] == control and row["weight_profile"] == profile
    ]
    winner_counts: dict[str, int] = {}
    wins_by_grid = {spec["name"]: 0 for spec in GRID_SPECS}
    wins_by_n = {str(n): 0 for n in DEFAULT_N_VALUES}
    wins_by_support = {}
    for row in selected:
        winner = row["winner_encoding"]
        winner_counts[winner] = winner_counts.get(winner, 0) + 1
        if winner == "heterogeneous_compute_reusable":
            wins_by_grid[row["grid"]] = wins_by_grid.get(row["grid"], 0) + 1
            n_key = str(row["N"])
            wins_by_n[n_key] = wins_by_n.get(n_key, 0) + 1
            wins_by_support[row["support_strategy"]] = wins_by_support.get(row["support_strategy"], 0) + 1
    wins = winner_counts.get("heterogeneous_compute_reusable", 0)
    return {
        "control_type": control,
        "weight_profile": profile,
        "tolerance_row_count": len(selected),
        "winner_counts": winner_counts,
        "heterogeneous_reusable_tolerance_win_count": wins,
        "heterogeneous_reusable_tolerance_win_fraction": wins / len(selected) if selected else None,
        "heterogeneous_reusable_wins_by_grid": wins_by_grid,
        "heterogeneous_reusable_wins_by_N": wins_by_n,
        "heterogeneous_reusable_wins_by_support_strategy": wins_by_support,
    }


def profile_passes(summary: dict[str, Any]) -> bool:
    tolerance_fraction = summary.get("heterogeneous_reusable_tolerance_win_fraction") or 0.0
    same_support_median = summary.get("median_reusable_objective_savings_vs_best_baseline")
    grid_wins = summary.get("heterogeneous_reusable_wins_by_grid", {})
    n_wins = summary.get("heterogeneous_reusable_wins_by_N", {})
    large_n_win = n_wins.get("64", 0) > 0 or n_wins.get("96", 0) > 0
    return (
        tolerance_fraction >= 0.50
        and all(grid_wins.get(spec["name"], 0) > 0 for spec in GRID_SPECS)
        and sum(1 for value in n_wins.values() if value > 0) >= 3
        and large_n_win
        and same_support_median is not None
        and same_support_median > 0
    )


def control_is_similar(real: dict[str, Any], control: dict[str, Any]) -> bool:
    real_fraction = real.get("heterogeneous_reusable_tolerance_win_fraction") or 0.0
    control_fraction = control.get("heterogeneous_reusable_tolerance_win_fraction") or 0.0
    real_median = real.get("median_reusable_objective_savings_vs_best_baseline")
    control_median = control.get("median_reusable_objective_savings_vs_best_baseline")
    if real_median is None or real_median <= 0 or control_median is None:
        return False
    return (
        control_fraction >= CONTROL_SIMILARITY_RATIO * real_fraction
        and control_median >= CONTROL_SIMILARITY_RATIO * real_median
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


def independent_channels(
    labels: list[int],
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
    for j in labels:
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


def audit_rows(rows: list[dict[str, Any]], profiles: tuple[WeightProfile, ...]) -> dict[str, Any]:
    errors: list[str] = []
    profile_map = {profile.name: profile for profile in profiles}
    for idx, row in enumerate(rows):
        cols = row["columns"]
        labels = row["cost_labels"]
        n = row["N"]
        if len(cols) != len(labels):
            errors.append(f"row {idx}: column/cost-label size mismatch")
        if cols != sorted(set(cols)):
            errors.append(f"row {idx}: columns are not sorted unique")
        if labels != sorted(set(labels)):
            errors.append(f"row {idx}: cost labels are not sorted unique")
        if any(col < 1 or col > n for col in cols + labels):
            errors.append(f"row {idx}: column or cost label outside 1..N")
        residual = row["certified_upper_relative_residual"]
        if residual > 0 and row["certified_residual_information_bits"] is not None:
            recomputed_info = -math.log2(residual)
            if abs(recomputed_info - row["certified_residual_information_bits"]) > 1e-9:
                errors.append(f"row {idx}: information bits mismatch")
        if row["encoding"].startswith("heterogeneous_compute"):
            profile = profile_map[row["weight_profile"]]
            expected = independent_channels(
                labels,
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


def classify(profile_summaries: dict[str, dict[str, dict[str, Any]]], audit: dict[str, Any]) -> tuple[str, list[str]]:
    reasons: list[str] = []
    if audit["audit_status"] != "PASS":
        return "BLOCKED_BY_IMPLEMENTATION_AUDIT", [f"audit_error_count={audit['error_count']}"]

    calibrated_real_passes = []
    contaminated_controls = {}
    for profile in CALIBRATED_PROFILES:
        real = profile_summaries[profile]["real"]
        if not profile_passes(real):
            continue
        calibrated_real_passes.append(profile)
        for control in NEGATIVE_CONTROLS:
            control_summary = profile_summaries[profile][control]
            if profile_passes(control_summary) and control_is_similar(real, control_summary):
                contaminated_controls.setdefault(profile, []).append(control)

    stress_real_passes = [
        profile for profile in STRESS_PROFILES if profile_passes(profile_summaries[profile]["real"])
    ]
    reasons.append(f"calibrated_real_passes={calibrated_real_passes}")
    reasons.append(f"stress_real_passes={stress_real_passes}")
    reasons.append(f"contaminated_controls={contaminated_controls}")
    if calibrated_real_passes and not contaminated_controls:
        return "CALIBRATED_HETEROGENEOUS_COMPUTE_SIGNAL", reasons
    if calibrated_real_passes and contaminated_controls:
        return "CONTROL_CONTAMINATED_SIGNAL", reasons
    if stress_real_passes:
        return "WEIGHT_DEPENDENT_SIGNAL", reasons
    return "FAILED_CALIBRATION", reasons


def summarize(
    rows: list[dict[str, Any]],
    tolerances: list[dict[str, Any]],
    profiles: tuple[WeightProfile, ...],
    audit: dict[str, Any],
) -> dict[str, Any]:
    control_names = ("real", "shuffled_indices", "random_support", "synthetic_nonfactorized")
    profile_summaries: dict[str, dict[str, dict[str, Any]]] = {}
    for profile in profiles:
        by_control = {}
        for control in control_names:
            combined = control_profile_tolerance_summary(tolerances, control, profile.name)
            combined.update(same_support_summary(rows, control, profile.name))
            by_control[control] = combined
        profile_summaries[profile.name] = by_control
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
    n_values = tuple(n_values)
    profiles, calibration_metadata = weight_profiles(max(n_values))
    rows: list[dict[str, Any]] = []
    for grid_spec in GRID_SPECS:
        x, w = build_grid(grid_spec["name"], grid_spec["nodes"])
        for n in n_values:
            for spec in DICTIONARIES:
                key = f"{grid_spec['name']}:{spec.name}:{n}"
                theta = theta_values(spec, n)
                A, y = build_design(theta, x, w)
                for control_type, control_A, supports_by_strategy, shuffled in control_cases(A, y, n, key):
                    for support_strategy, supports_by_budget in supports_by_strategy.items():
                        for cols in supports_by_budget.values():
                            labels = shuffled_cost_labels(cols, n, key) if shuffled else cols
                            add_rows_for_support(
                                rows,
                                spec,
                                n,
                                grid_spec["name"],
                                control_type,
                                support_strategy,
                                control_A,
                                y,
                                cols,
                                labels,
                                profiles,
                            )
    tolerances = tolerance_table(rows)
    audit = audit_rows(rows, profiles)
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
        "weight_profiles": [asdict(profile) for profile in profiles],
        "calibration_metadata": calibration_metadata,
        "negative_controls": sorted(NEGATIVE_CONTROLS),
        "control_similarity_ratio": CONTROL_SIMILARITY_RATIO,
        "residual_tolerances": RESIDUAL_TOLERANCES,
        "summary": summarize(rows, tolerances, profiles, audit),
        "tolerance_table": tolerances,
        "rows": rows,
    }


def markdown_report(results: dict[str, Any]) -> str:
    summary = results["summary"]
    profiles = summary["profile_summaries"]
    lines = [
        f"# {EXPERIMENT_ID} Report",
        "",
        "## Verdict",
        "",
        f"- Status: `{summary['status']}`",
        f"- Rows: `{summary['row_count']}`",
        f"- Tolerance rows: `{summary['tolerance_row_count']}`",
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
        "This calibration gate tests whether the heterogeneous-compute signal",
        "survives externally motivated operation weights and negative controls.",
        "",
        "The result is still a finite operation-weighted diagnostic. It is not",
        "a theorem about primes, not a zeta formalization, and not a quantum",
        "claim.",
        "",
        "## Claim Ceiling",
        "",
        results["claim_ceiling"],
        "",
        "## Real Dictionary Profiles",
        "",
        "| profile | real tolerance wins | tolerance rows | real same-support wins | comparable | median savings |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for profile in sorted(profiles):
        real = profiles[profile]["real"]
        lines.append(
            f"| `{profile}` | `{real['heterogeneous_reusable_tolerance_win_count']}` | "
            f"`{real['tolerance_row_count']}` | "
            f"`{real['heterogeneous_reusable_same_support_win_count']}` | "
            f"`{real['same_support_comparable_count']}` | "
            f"`{real['median_reusable_objective_savings_vs_best_baseline']}` |"
        )
    lines.extend(
        [
            "",
            "## Negative Controls",
            "",
            "| profile | control | tolerance wins | tolerance rows | median savings |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for profile in sorted(profiles):
        for control in sorted(NEGATIVE_CONTROLS):
            row = profiles[profile][control]
            lines.append(
                f"| `{profile}` | `{control}` | "
                f"`{row['heterogeneous_reusable_tolerance_win_count']}` | "
                f"`{row['tolerance_row_count']}` | "
                f"`{row['median_reusable_objective_savings_vs_best_baseline']}` |"
            )
    lines.extend(
        [
            "",
            "## Calibration Metadata",
            "",
            "```json",
            json.dumps(results["calibration_metadata"], indent=2, sort_keys=True),
            "```",
            "",
            "## Implementation Audit",
            "",
            "```json",
            json.dumps(summary["implementation_audit"], indent=2, sort_keys=True),
            "```",
            "",
            "## Boundary",
            "",
            "A positive result here means only that this finite cost model survived",
            "the stated calibration and control gate. A contaminated result means",
            "the negative controls look too similar and the interpretation must be",
            "downgraded to cost-model artifact risk.",
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
