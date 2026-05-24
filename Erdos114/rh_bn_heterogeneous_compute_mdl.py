#!/usr/bin/env python3
"""
EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01

Finite heterogeneous-compute MDL diagnostic for the Beurling-Nyman lane.

Claim ceiling:
    INTERNAL / METHOD-SHAPING ONLY.
    This is a finite operation-weighted cost experiment. It is not an RH claim,
    not a zeta formalization, and not a quantum-mechanical claim.
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

from rh_beurling_nyman_mdl_probe import DICTIONARIES, DictionarySpec, theta_values
from rh_bn_prime_harness_mdl import (
    GRID_SPECS,
    PRIMARY_BITS,
    RESIDUAL_TOLERANCES,
    SUPPORT_BUDGETS,
    build_design,
    build_grid,
    composite_address_bits,
    composite_prefix,
    factor_cost_ordered,
    factor_harness_bits,
    factorization,
    factorized_address_bits,
    fit_columns,
    flat_address_bits,
    int_bits_at_most,
    integer_prefix,
    primes_up_to,
    row_for_encoding,
)


EXPERIMENT_ID = "EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01"
CLAIM_CEILING = (
    "INTERNAL / METHOD-SHAPING ONLY: finite operation-weighted "
    "Beurling-Nyman MDL diagnostic; no RH claim, no zeta formalization, and "
    "no quantum claim."
)
MODEL_NOTE = "HETEROGENEOUS_COMPUTE_MDL_MODEL_2026-05-06.md"
SOURCE_EXPERIMENTS = (
    "EXP-MATH-RH-BN-COMPOSITE-FIRST-PROJECTION-20260506-01",
    "EXP-MATH-RH-BN-PRIME-HARNESS-MDL-20260506-01",
    "EXP-MATH-RH-BN-PRIME-HARNESS-AMORTIZED-20260506-01",
)
DEFAULT_N_VALUES = (8, 16, 24, 32, 48)
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


def information_bits(row: dict[str, Any]) -> float | None:
    residual = row.get("certified_upper_relative_residual_16")
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
    if encoding == "composite_index":
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


def operation_channels(indices: list[int], n: int, active: int, include_harness: bool) -> dict[str, int]:
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
        "prime_lookup": len(distinct_primes) * prime_lookup_unit,
        "multiply": multiply_ops * MULTIPLY_WEIGHT,
        "exponent": exponent_ops * EXPONENT_WEIGHT,
        "factor_tree_depth": max_depth * FACTOR_TREE_DEPTH_WEIGHT,
        "coefficient_move": active * PRIMARY_BITS,
        "shared_harness": factor_harness_bits(n) if include_harness else 0,
    }


def cost_from_channels(spec: DictionarySpec, channels: dict[str, int]) -> int:
    return spec.header_bits + STAGE_HEADER_BITS + sum(channels.values())


def add_baseline_row(
    rows: list[dict[str, Any]],
    spec: DictionarySpec,
    n: int,
    grid: str,
    support_strategy: str,
    encoding: str,
    cols: list[int],
    fit: dict[str, Any],
    address_bits: int,
) -> None:
    row = row_for_encoding(
        spec,
        n,
        grid,
        support_strategy,
        encoding,
        cols,
        fit,
        address_bits,
        False,
    )
    channels = baseline_channels(encoding, address_bits, row)
    rows.append(attach_objective(row, int(row["total_description_bits"]), channels))


def add_heterogeneous_rows(
    rows: list[dict[str, Any]],
    spec: DictionarySpec,
    n: int,
    grid: str,
    support_strategy: str,
    cols: list[int],
    fit: dict[str, Any],
) -> None:
    for encoding, include_harness in (
        ("heterogeneous_compute_reusable", False),
        ("heterogeneous_compute_with_harness", True),
    ):
        channels = operation_channels(cols, n, fit["active_count"], include_harness)
        compute_cost = cost_from_channels(spec, channels)
        row = row_for_encoding(
            spec,
            n,
            grid,
            support_strategy,
            encoding,
            cols,
            fit,
            sum(channels.values()) - channels["coefficient_move"],
            include_harness,
        )
        row["total_description_bits"] = compute_cost
        row["address_bits"] = sum(channels.values()) - channels["coefficient_move"]
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
    fit = fit_columns(A, y, [j - 1 for j in cols])
    if support_strategy == "integer_prefix":
        add_baseline_row(
            rows,
            spec,
            n,
            grid,
            support_strategy,
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
            "factorized_reusable_harness",
            cols,
            fit,
            factorized_address_bits(cols, n, include_harness=False),
        )
    else:
        raise ValueError(f"unknown support strategy {support_strategy}")
    add_heterogeneous_rows(rows, spec, n, grid, support_strategy, cols, fit)


def support_budgets(n: int) -> list[int]:
    values = sorted({m for m in SUPPORT_BUDGETS if m <= n} | {n})
    return [m for m in values if m > 0]


def tolerance_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, int, float], list[dict[str, Any]]] = {}
    for row in rows:
        for tol in RESIDUAL_TOLERANCES:
            key = (row["grid"], row["dictionary"], row["N"], tol)
            groups.setdefault(key, []).append(row)

    table = []
    for (grid, dictionary, n, tol), group_rows in sorted(groups.items()):
        candidates = [
            row
            for row in group_rows
            if row["certified_upper_relative_residual_16"] <= tol
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
                "grid": grid,
                "dictionary": dictionary,
                "N": n,
                "tolerance": tol,
                "winner_encoding": winner["encoding"],
                "winner_support_strategy": winner["support_strategy"],
                "winner_support_budget": winner["support_budget"],
                "winner_compute_cost_bits": winner["compute_cost_bits"],
                "winner_total_objective_bits": winner["total_objective_bits"],
                "winner_certified_upper_relative_residual_16": winner[
                    "certified_upper_relative_residual_16"
                ],
                "encoding_best": [
                    {
                        "encoding": enc,
                        "support_strategy": row["support_strategy"],
                        "support_budget": row["support_budget"],
                        "compute_cost_bits": row["compute_cost_bits"],
                        "total_objective_bits": row["total_objective_bits"],
                        "certified_upper_relative_residual_16": row[
                            "certified_upper_relative_residual_16"
                        ],
                    }
                    for enc, row in sorted(by_encoding.items())
                ],
            }
        )
    return table


def same_support_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    groups: dict[tuple[str, str, int, int, str], dict[str, dict[str, Any]]] = {}
    for row in rows:
        key = (
            row["grid"],
            row["dictionary"],
            row["N"],
            row["support_budget"],
            row["support_strategy"],
        )
        groups.setdefault(key, {})[row["encoding"]] = row

    comparable = 0
    reusable_wins = 0
    with_harness_wins = 0
    reusable_savings = []
    with_harness_savings = []
    for encs in groups.values():
        baselines = [
            encs.get("flat_index"),
            encs.get("composite_index"),
            encs.get("factorized_reusable_harness"),
        ]
        baselines = [row for row in baselines if row is not None]
        reusable = encs.get("heterogeneous_compute_reusable")
        with_harness = encs.get("heterogeneous_compute_with_harness")
        if not baselines or reusable is None or with_harness is None:
            continue
        comparable += 1
        best_baseline = min(baselines, key=lambda row: row["total_objective_bits"])
        reusable_delta = best_baseline["total_objective_bits"] - reusable["total_objective_bits"]
        harness_delta = best_baseline["total_objective_bits"] - with_harness["total_objective_bits"]
        reusable_savings.append(reusable_delta)
        with_harness_savings.append(harness_delta)
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
            float(np.median(with_harness_savings)) if with_harness_savings else None
        ),
    }


def classify(tolerances: list[dict[str, Any]], same_support: dict[str, Any]) -> str:
    if not tolerances:
        return "NO_COMPARABLE_TOLERANCE_ROWS"
    winner_counts: dict[str, int] = {}
    grid_winners: dict[str, set[str]] = {}
    n_winners: dict[int, set[str]] = {}
    for row in tolerances:
        winner = row["winner_encoding"]
        winner_counts[winner] = winner_counts.get(winner, 0) + 1
        grid_winners.setdefault(row["grid"], set()).add(winner)
        n_winners.setdefault(row["N"], set()).add(winner)
    total = len(tolerances)
    hetero_total = (
        winner_counts.get("heterogeneous_compute_reusable", 0)
        + winner_counts.get("heterogeneous_compute_with_harness", 0)
    )
    hetero_grids = {
        grid
        for grid, winners in grid_winners.items()
        if "heterogeneous_compute_reusable" in winners
        or "heterogeneous_compute_with_harness" in winners
    }
    hetero_n = {
        n
        for n, winners in n_winners.items()
        if "heterogeneous_compute_reusable" in winners
        or "heterogeneous_compute_with_harness" in winners
    }
    same_support_fraction = same_support.get("heterogeneous_reusable_same_support_win_fraction") or 0.0
    if hetero_total / total >= 0.50 and len(hetero_grids) >= 2 and len(hetero_n) >= 3:
        if same_support_fraction >= 0.50:
            return "STABLE_HETEROGENEOUS_COMPUTE_SIGNAL"
        return "TOLERANCE_ONLY_HETEROGENEOUS_COMPUTE_SIGNAL"
    if hetero_total > 0:
        return "MIXED_HETEROGENEOUS_COMPUTE_SIGNAL"
    return "NO_HETEROGENEOUS_COMPUTE_SIGNAL"


def summarize(rows: list[dict[str, Any]], tolerances: list[dict[str, Any]]) -> dict[str, Any]:
    winner_counts: dict[str, int] = {}
    grid_counts: dict[str, dict[str, int]] = {}
    n_counts: dict[str, dict[str, int]] = {}
    for row in tolerances:
        winner = row["winner_encoding"]
        winner_counts[winner] = winner_counts.get(winner, 0) + 1
        grid_counts.setdefault(row["grid"], {})[winner] = grid_counts.setdefault(row["grid"], {}).get(winner, 0) + 1
        n_key = str(row["N"])
        n_counts.setdefault(n_key, {})[winner] = n_counts.setdefault(n_key, {}).get(winner, 0) + 1
    same_support = same_support_summary(rows)
    status = classify(tolerances, same_support)
    return {
        "status": status,
        "row_count": len(rows),
        "tolerance_row_count": len(tolerances),
        "tolerance_winner_counts": winner_counts,
        "tolerance_winner_counts_by_grid": grid_counts,
        "tolerance_winner_counts_by_N": n_counts,
        **same_support,
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
        "model_note": MODEL_NOTE,
        "source_experiments": SOURCE_EXPERIMENTS,
        "cost_channels": {
            "direct_index": "direct finite lookup of selected dictionary indices",
            "prime_lookup": "lookup of distinct prime tokens used by selected indices",
            "multiply": f"{MULTIPLY_WEIGHT} bit-equivalent units per factor composition",
            "exponent": f"{EXPONENT_WEIGHT} bit-equivalent units per repeated-factor operation",
            "factor_tree_depth": f"{FACTOR_TREE_DEPTH_WEIGHT} bit-equivalent unit per max factor-tree level",
            "coefficient_move": f"{PRIMARY_BITS} bits per active coefficient",
            "shared_harness": "optional setup cost for declaring the reusable prime harness",
        },
        "objective": "total_objective_bits = compute_cost_bits - certified_residual_information_bits",
        "grids": GRID_SPECS,
        "n_values": list(n_values),
        "support_budgets": SUPPORT_BUDGETS,
        "residual_tolerances": RESIDUAL_TOLERANCES,
        "primary_coeff_bits": PRIMARY_BITS,
        "summary": summarize(rows, tolerances),
        "tolerance_table": tolerances,
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
        f"- Tolerance rows: `{summary['tolerance_row_count']}`",
        (
            "- Heterogeneous reusable same-support wins: "
            f"`{summary['heterogeneous_reusable_same_support_win_count']}` / "
            f"`{summary['same_support_comparable_count']}`"
        ),
        (
            "- Heterogeneous with-harness same-support wins: "
            f"`{summary['heterogeneous_with_harness_same_support_win_count']}` / "
            f"`{summary['same_support_comparable_count']}`"
        ),
        (
            "- Median reusable objective savings vs best baseline: "
            f"`{summary['median_reusable_objective_savings_vs_best_baseline']}`"
        ),
        (
            "- Median with-harness objective savings vs best baseline: "
            f"`{summary['median_with_harness_objective_savings_vs_best_baseline']}`"
        ),
        "",
        "Tolerance winner counts:",
        "",
        "```json",
        json.dumps(summary["tolerance_winner_counts"], indent=2, sort_keys=True),
        "```",
        "",
        "## Meaning",
        "",
        "This experiment tests the TDP-inspired correction to the finite RH-MDL",
        "cost model: when flat bits and factorization strings fail, charge",
        "arithmetic operations as heterogeneous compute channels.",
        "",
        "The result is a finite operation-weighted diagnostic. It is not a theorem",
        "about primes, not a zeta formalization, and not a quantum-mechanical",
        "claim.",
        "",
        "## Claim Ceiling",
        "",
        results["claim_ceiling"],
        "",
        "## Objective",
        "",
        "```text",
        results["objective"],
        "```",
        "",
        "## Tolerance Winners",
        "",
        "| grid | dictionary | N | tolerance | winner | support | objective | residual |",
        "|---|---|---:|---:|---|---|---:|---:|",
    ]
    for row in results["tolerance_table"][:60]:
        lines.append(
            "| {grid} | {dictionary} | {N} | {tolerance:.2f} | {winner_encoding} | "
            "{winner_support_strategy} | {winner_total_objective_bits:.3f} | "
            "{winner_certified_upper_relative_residual_16:.6g} |".format(**row)
        )
    remaining = len(results["tolerance_table"]) - 60
    if remaining > 0:
        lines.append(f"| ... | ... | ... | ... | {remaining} more rows | ... | ... | ... |")
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "A positive result here means only that this finite operation-weighted",
            "cost model beats the listed finite baselines under the stated objective.",
            "It does not establish an asymptotic law, dictionary invariance, or a",
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
