#!/usr/bin/env python3
"""
EXP-MATH-RH-BN-PRIME-HARNESS-AMORTIZED-20260506-01

Finite amortization diagnostic for the prime-harness MDL model.

Claim ceiling:
    INTERNAL / METHOD-SHAPING ONLY.
    This tests finite description cost across a small bundle of related
    fractional-part approximation tasks. It is not an RH claim, not a zeta
    formalization, and not a quantum-mechanical claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import date
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
from numpy.polynomial.legendre import leggauss

from rh_beurling_nyman_mdl_probe import DICTIONARIES, fractional_part, theta_values
from rh_bn_prime_harness_mdl import (
    PRIMARY_BITS,
    CLAIM_CEILING as PRIOR_CLAIM_CEILING,
    composite_address_bits,
    composite_prefix,
    factor_cost_ordered,
    factor_harness_bits,
    factorized_address_bits,
    fit_columns,
    flat_address_bits,
    integer_prefix,
    row_for_encoding,
)


EXPERIMENT_ID = "EXP-MATH-RH-BN-PRIME-HARNESS-AMORTIZED-20260506-01"
SOURCE_EXPERIMENT_ID = "EXP-MATH-RH-BN-PRIME-HARNESS-MDL-20260506-01"
CLAIM_CEILING = (
    "INTERNAL / METHOD-SHAPING ONLY: finite amortized encoding-cost diagnostic; "
    "no RH claim, no zeta formalization, and no quantum claim."
)
DEFAULT_N_VALUES = (16, 32, 48)
GRID_SPECS = (
    {"name": "legendre_2048", "rule": "Gauss-Legendre on [0,1]", "nodes": 2048},
    {"name": "midpoint_2048", "rule": "equal-weight midpoint grid on [0,1]", "nodes": 2048},
)
SUPPORT_BUDGETS = (4, 8, 12, 16, 24, 32, 48)
RESIDUAL_TOLERANCES = (0.25, 0.20, 0.15, 0.12, 0.10, 0.08, 0.06)


TargetFn = Callable[[np.ndarray], np.ndarray]


TARGETS: tuple[tuple[str, str, TargetFn], ...] = (
    ("bn_constant", "f(x) = 1", lambda x: np.ones_like(x)),
    ("sqrt_x", "f(x) = sqrt(x)", lambda x: np.sqrt(x)),
    ("linear_x", "f(x) = x", lambda x: x),
    ("left_half_step", "f(x) = 1[x <= 1/2]", lambda x: (x <= 0.5).astype(float)),
    ("right_half_step", "f(x) = 1[x > 1/2]", lambda x: (x > 0.5).astype(float)),
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_grid(name: str, nodes: int) -> tuple[np.ndarray, np.ndarray]:
    if name.startswith("legendre"):
        raw_nodes, raw_weights = leggauss(nodes)
        return 0.5 * (raw_nodes + 1.0), 0.5 * raw_weights
    if name.startswith("midpoint"):
        x = (np.arange(nodes, dtype=np.float64) + 0.5) / nodes
        w = np.full(nodes, 1.0 / nodes, dtype=np.float64)
        return x, w
    raise ValueError(f"unknown grid {name}")


def build_weighted_design(theta: np.ndarray, x: np.ndarray, w: np.ndarray) -> np.ndarray:
    sqrt_w = np.sqrt(w)
    basis = fractional_part(theta[None, :] / x[:, None])
    return basis * sqrt_w[:, None]


def weighted_target(target_fn: TargetFn, x: np.ndarray, w: np.ndarray) -> np.ndarray:
    return target_fn(x) * np.sqrt(w)


def support_budgets(n: int) -> list[int]:
    values = sorted({m for m in SUPPORT_BUDGETS if m <= n} | {n})
    return [m for m in values if m > 0]


def add_integer_prefix_rows(
    rows: list[dict[str, Any]],
    spec,
    n: int,
    grid: str,
    target_name: str,
    target_description: str,
    A: np.ndarray,
    y: np.ndarray,
    m: int,
) -> None:
    cols = integer_prefix(m)
    fit = fit_columns(A, y, [j - 1 for j in cols])
    flat = row_for_encoding(
        spec,
        n,
        grid,
        "integer_prefix",
        "flat_index",
        cols,
        fit,
        flat_address_bits(cols, n),
        False,
    )
    factor = row_for_encoding(
        spec,
        n,
        grid,
        "integer_prefix",
        "factorized_reusable_harness",
        cols,
        fit,
        factorized_address_bits(cols, n, include_harness=False),
        False,
    )
    for row in (flat, factor):
        row["target"] = target_name
        row["target_description"] = target_description
        rows.append(row)


def add_composite_prefix_row(
    rows: list[dict[str, Any]],
    spec,
    n: int,
    grid: str,
    target_name: str,
    target_description: str,
    A: np.ndarray,
    y: np.ndarray,
    m: int,
) -> None:
    cols = composite_prefix(n, m)
    if len(cols) < 2:
        return
    address_bits = composite_address_bits(cols, n)
    if address_bits is None:
        return
    fit = fit_columns(A, y, [j - 1 for j in cols])
    row = row_for_encoding(
        spec,
        n,
        grid,
        "composite_prefix",
        "composite_index",
        cols,
        fit,
        address_bits,
        False,
    )
    row["target"] = target_name
    row["target_description"] = target_description
    rows.append(row)


def add_factor_cost_row(
    rows: list[dict[str, Any]],
    spec,
    n: int,
    grid: str,
    target_name: str,
    target_description: str,
    A: np.ndarray,
    y: np.ndarray,
    m: int,
) -> None:
    cols = factor_cost_ordered(n, m)
    fit = fit_columns(A, y, [j - 1 for j in cols])
    row = row_for_encoding(
        spec,
        n,
        grid,
        "factor_cost_ordered",
        "factorized_reusable_harness",
        cols,
        fit,
        factorized_address_bits(cols, n, include_harness=False),
        False,
    )
    row["target"] = target_name
    row["target_description"] = target_description
    rows.append(row)


def best_for_task(
    rows: list[dict[str, Any]],
    target: str,
    tolerance: float,
    encoding: str,
) -> dict[str, Any] | None:
    candidates = [
        row
        for row in rows
        if row["target"] == target
        and row["encoding"] == encoding
        and row["certified_upper_relative_residual_16"] <= tolerance
    ]
    if not candidates:
        return None
    return min(candidates, key=lambda row: row["total_description_bits"])


def bundle_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, int], list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault((row["grid"], row["dictionary"], row["N"]), []).append(row)

    bundles = []
    target_names = [name for name, _desc, _fn in TARGETS]
    for (grid, dictionary, n), group_rows in sorted(groups.items()):
        common_tasks = []
        flat_total = 0
        factor_total = 0
        composite_total = 0
        composite_common = True
        task_winners: dict[str, int] = {}
        for target in target_names:
            for tolerance in RESIDUAL_TOLERANCES:
                flat = best_for_task(group_rows, target, tolerance, "flat_index")
                factor = best_for_task(
                    group_rows, target, tolerance, "factorized_reusable_harness"
                )
                if flat is None or factor is None:
                    continue
                composite = best_for_task(group_rows, target, tolerance, "composite_index")
                flat_bits = flat["total_description_bits"]
                factor_bits = factor["total_description_bits"]
                flat_total += flat_bits
                factor_total += factor_bits
                if composite is None:
                    composite_common = False
                else:
                    composite_total += composite["total_description_bits"]
                winner = "flat_index" if flat_bits <= factor_bits else "factorized_reusable_harness"
                task_winners[winner] = task_winners.get(winner, 0) + 1
                common_tasks.append(
                    {
                        "target": target,
                        "tolerance": tolerance,
                        "flat_bits": flat_bits,
                        "factor_reusable_bits": factor_bits,
                        "factor_minus_flat_bits": factor_bits - flat_bits,
                        "best_flat_support": flat["support_strategy"],
                        "best_factor_support": factor["support_strategy"],
                    }
                )
        if not common_tasks:
            continue
        harness_bits = factor_harness_bits(n)
        factor_shared_total = factor_total + harness_bits
        composite_total_value = composite_total if composite_common else None
        bundles.append(
            {
                "grid": grid,
                "dictionary": dictionary,
                "N": n,
                "common_task_count": len(common_tasks),
                "harness_setup_bits": harness_bits,
                "flat_total_bits": flat_total,
                "factor_reusable_total_bits_without_setup": factor_total,
                "factor_shared_total_bits": factor_shared_total,
                "composite_total_bits_if_all_tasks_representable": composite_total_value,
                "factor_reusable_savings_vs_flat_without_setup": flat_total - factor_total,
                "factor_shared_savings_vs_flat": flat_total - factor_shared_total,
                "factor_shared_wins": factor_shared_total < flat_total,
                "task_winner_counts_flat_vs_factor_without_setup": task_winners,
                "tasks": common_tasks,
            }
        )
    return bundles


def classify(bundles: list[dict[str, Any]]) -> str:
    if not bundles:
        return "NO_COMPARABLE_BUNDLES"
    shared_wins = sum(1 for bundle in bundles if bundle["factor_shared_wins"])
    win_fraction = shared_wins / len(bundles)
    median_savings = float(
        np.median([bundle["factor_shared_savings_vs_flat"] for bundle in bundles])
    )
    reusable_positive = sum(
        1 for bundle in bundles if bundle["factor_reusable_savings_vs_flat_without_setup"] > 0
    )
    if win_fraction >= 0.50 and median_savings > 0:
        return "PRIME_HARNESS_AMORTIZATION_SIGNAL"
    if reusable_positive / len(bundles) >= 0.50:
        return "REUSABLE_ONLY_NO_SETUP_RECOVERY"
    return "NO_AMORTIZED_PRIME_HARNESS_SIGNAL"


def summarize(rows: list[dict[str, Any]], bundles: list[dict[str, Any]]) -> dict[str, Any]:
    shared_wins = sum(1 for bundle in bundles if bundle["factor_shared_wins"])
    reusable_positive = sum(
        1 for bundle in bundles if bundle["factor_reusable_savings_vs_flat_without_setup"] > 0
    )
    savings = [bundle["factor_shared_savings_vs_flat"] for bundle in bundles]
    reusable_savings = [
        bundle["factor_reusable_savings_vs_flat_without_setup"] for bundle in bundles
    ]
    best = max(bundles, key=lambda bundle: bundle["factor_shared_savings_vs_flat"]) if bundles else None
    return {
        "status": classify(bundles),
        "row_count": len(rows),
        "bundle_count": len(bundles),
        "shared_harness_win_count": shared_wins,
        "shared_harness_win_fraction": shared_wins / len(bundles) if bundles else None,
        "reusable_without_setup_positive_count": reusable_positive,
        "reusable_without_setup_positive_fraction": (
            reusable_positive / len(bundles) if bundles else None
        ),
        "median_shared_harness_savings_vs_flat": float(np.median(savings)) if savings else None,
        "median_reusable_savings_vs_flat_without_setup": (
            float(np.median(reusable_savings)) if reusable_savings else None
        ),
        "best_shared_harness_bundle": (
            {
                "grid": best["grid"],
                "dictionary": best["dictionary"],
                "N": best["N"],
                "common_task_count": best["common_task_count"],
                "factor_shared_savings_vs_flat": best["factor_shared_savings_vs_flat"],
                "harness_setup_bits": best["harness_setup_bits"],
            }
            if best
            else None
        ),
    }


def run(n_values: Iterable[int]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for grid_spec in GRID_SPECS:
        x, w = build_grid(grid_spec["name"], grid_spec["nodes"])
        for n in n_values:
            for spec in DICTIONARIES:
                theta = theta_values(spec, n)
                A = build_weighted_design(theta, x, w)
                for target_name, target_description, target_fn in TARGETS:
                    y = weighted_target(target_fn, x, w)
                    for m in support_budgets(n):
                        add_integer_prefix_rows(
                            rows,
                            spec,
                            n,
                            grid_spec["name"],
                            target_name,
                            target_description,
                            A,
                            y,
                            m,
                        )
                        add_composite_prefix_row(
                            rows,
                            spec,
                            n,
                            grid_spec["name"],
                            target_name,
                            target_description,
                            A,
                            y,
                            m,
                        )
                        add_factor_cost_row(
                            rows,
                            spec,
                            n,
                            grid_spec["name"],
                            target_name,
                            target_description,
                            A,
                            y,
                            m,
                        )
    bundles = bundle_rows(rows)
    return {
        "experiment_id": EXPERIMENT_ID,
        "date": date.today().isoformat(),
        "claim_ceiling": CLAIM_CEILING,
        "source_experiment_id": SOURCE_EXPERIMENT_ID,
        "source_claim_ceiling": PRIOR_CLAIM_CEILING,
        "model_note": "PRIME_HARNESS_MDL_MODEL_2026-05-06.md",
        "targets": [
            {"name": name, "description": description}
            for name, description, _target_fn in TARGETS
        ],
        "grids": GRID_SPECS,
        "n_values": list(n_values),
        "support_budgets": SUPPORT_BUDGETS,
        "residual_tolerances": RESIDUAL_TOLERANCES,
        "primary_coeff_bits": PRIMARY_BITS,
        "summary": summarize(rows, bundles),
        "bundles": bundles,
        "rows": rows,
    }


def markdown_report(results: dict[str, Any]) -> str:
    summary = results["summary"]
    lines = [
        f"# {EXPERIMENT_ID} Report",
        "",
        "## Verdict",
        "",
        f"- Status: `{summary['status']}`",
        f"- Rows: `{summary['row_count']}`",
        f"- Bundles: `{summary['bundle_count']}`",
        f"- Shared-harness wins: `{summary['shared_harness_win_count']}`",
        (
            "- Reusable-without-setup positive bundles: "
            f"`{summary['reusable_without_setup_positive_count']}`"
        ),
        (
            "- Median shared-harness savings vs flat: "
            f"`{summary['median_shared_harness_savings_vs_flat']}` bits"
        ),
        (
            "- Median reusable savings before setup: "
            f"`{summary['median_reusable_savings_vs_flat_without_setup']}` bits"
        ),
        "",
        "Best shared-harness bundle:",
        "",
        "```json",
        json.dumps(summary["best_shared_harness_bundle"], indent=2, sort_keys=True),
        "```",
        "",
        "## Meaning",
        "",
        "This experiment tests whether a prime-address harness becomes useful when",
        "its setup cost is shared across a bundle of related finite approximation",
        "tasks. The targets are simple weighted functions over the same",
        "fractional-part dictionary, including the standard constant target.",
        "",
        "The result is still finite encoding accounting. It is not a theorem about",
        "primes, not a zeta formalization, and not a quantum-mechanical claim.",
        "",
        "## Claim Ceiling",
        "",
        results["claim_ceiling"],
        "",
        "## Bundle Results",
        "",
        "| grid | dictionary | N | tasks | harness bits | flat bits | shared factor bits | savings | win |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for bundle in results["bundles"]:
        lines.append(
            "| {grid} | {dictionary} | {N} | {common_task_count} | {harness_setup_bits} | "
            "{flat_total_bits} | {factor_shared_total_bits} | "
            "{factor_shared_savings_vs_flat} | {factor_shared_wins} |".format(**bundle)
        )
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "A win here would mean only that this finite cost model rewards a shared",
            "factorization-address layer across several approximation tasks. A loss",
            "means the naive prime-harness code is still too expensive under this",
            "finite accounting model.",
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
