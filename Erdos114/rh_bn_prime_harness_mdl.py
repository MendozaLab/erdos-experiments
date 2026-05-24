#!/usr/bin/env python3
"""
EXP-MATH-RH-BN-PRIME-HARNESS-MDL-20260506-01

Finite prime-harness MDL diagnostic for the Beurling-Nyman lane.

Claim ceiling:
    INTERNAL / METHOD-SHAPING ONLY.
    This is a finite encoding-cost experiment. It is not an RH claim, not a
    zeta formalization, and not a quantum-mechanical claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import asdict
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


EXPERIMENT_ID = "EXP-MATH-RH-BN-PRIME-HARNESS-MDL-20260506-01"
CLAIM_CEILING = (
    "INTERNAL / METHOD-SHAPING ONLY: finite Beurling-Nyman encoding-cost "
    "diagnostic; no RH claim, no zeta formalization, and no quantum claim."
)
DEFAULT_N_VALUES = (8, 16, 24, 32, 48)
GRID_SPECS = (
    {"name": "legendre_2048", "rule": "Gauss-Legendre on [0,1]", "nodes": 2048},
    {"name": "midpoint_2048", "rule": "equal-weight midpoint grid on [0,1]", "nodes": 2048},
)
SUPPORT_BUDGETS = (4, 8, 12, 16, 24, 32, 48)
RESIDUAL_TOLERANCES = (0.25, 0.20, 0.15, 0.12, 0.10, 0.08, 0.06)
PRIMARY_BITS = 16
STAGE_HEADER_BITS = 8


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


def primes_up_to(n: int) -> list[int]:
    return [j for j in range(2, n + 1) if is_prime(j)]


def factorization(n: int) -> dict[int, int]:
    if n < 2:
        return {}
    out: dict[int, int] = {}
    d = 2
    value = n
    while d * d <= value:
        while value % d == 0:
            out[d] = out.get(d, 0) + 1
            value //= d
        d += 1 if d == 2 else 2
    if value > 1:
        out[value] = out.get(value, 0) + 1
    return out


def int_bits_at_most(n: int) -> int:
    return max(1, math.ceil(math.log2(max(2, n + 1))))


def exponent_bits(exp: int) -> int:
    return 1 + int(math.floor(math.log2(max(1, exp))))


def factor_address_bits(j: int, n: int) -> int:
    if j == 1:
        return 1
    prime_count = max(1, len(primes_up_to(n)))
    prime_token_bits = int_bits_at_most(prime_count)
    factors = factorization(j)
    # One leading type bit plus one separator bit per factor keeps the code
    # prefix-readable enough for this finite diagnostic.
    return 1 + sum(prime_token_bits + exponent_bits(exp) + 1 for exp in factors.values())


def factor_harness_bits(n: int) -> int:
    # Cost to declare the prime-address harness itself if it is not reusable.
    return len(primes_up_to(n)) * int_bits_at_most(n)


def flat_address_bits(indices: list[int], n: int) -> int:
    return len(indices) * int_bits_at_most(n)


def composite_address_bits(indices: list[int], n: int) -> int | None:
    composites = [j for j in range(2, n + 1) if not is_prime(j)]
    bucket_bits = int_bits_at_most(len(composites))
    total = 0
    for j in indices:
        if j == 1:
            total += 1
        elif is_prime(j):
            return None
        else:
            total += bucket_bits
    return total


def factorized_address_bits(indices: list[int], n: int, include_harness: bool) -> int:
    total = sum(factor_address_bits(j, n) for j in indices)
    if include_harness:
        total += factor_harness_bits(n)
    return total


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


def active_count(coeffs: np.ndarray) -> int:
    if coeffs.size == 0:
        return 0
    max_abs = float(np.max(np.abs(coeffs)))
    return int(np.sum(np.abs(coeffs) > ACTIVE_REL_TOL * max(1.0, max_abs)))


def fit_columns(A: np.ndarray, y: np.ndarray, cols_zero: list[int]) -> dict[str, Any]:
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
    active = active_count(coeffs)
    scale = max(1.0, float(np.max(np.abs(coeffs))) if coeffs.size else 1.0)
    step = scale * (2.0 ** (-PRIMARY_BITS))
    quant_penalty = sigma_max * math.sqrt(active) * step / 2.0
    certified = (residual_l2 + quant_penalty) / target_l2 if target_l2 else float("nan")
    return {
        "rank": int(rank),
        "active_count": active,
        "residual_l2": residual_l2,
        "relative_residual": residual_l2 / target_l2 if target_l2 else float("nan"),
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "condition_number": finite_float(cond),
        "coefficient_l1": float(np.linalg.norm(coeffs, 1)),
        "coefficient_l2": float(np.linalg.norm(coeffs)),
        "coefficient_linf": float(np.max(np.abs(coeffs))) if coeffs.size else 0.0,
        "quantization_step_16": step,
        "quantization_penalty_l2_bound_16": quant_penalty,
        "certified_upper_relative_residual_16": certified,
        "information_bits_certified": -math.log2(certified) if certified > 0 else None,
    }


def total_bits(spec: DictionarySpec, address_bits: int, active: int) -> int:
    return spec.header_bits + STAGE_HEADER_BITS + address_bits + active * PRIMARY_BITS


def row_for_encoding(
    spec: DictionarySpec,
    n: int,
    grid: str,
    support_strategy: str,
    encoding: str,
    cols: list[int],
    fit: dict[str, Any],
    address_bits: int,
    include_harness: bool,
) -> dict[str, Any]:
    bits = total_bits(spec, address_bits, fit["active_count"])
    info_bits = fit["information_bits_certified"]
    return {
        "dictionary": spec.name,
        "theta_rule": spec.theta_rule,
        "grid": grid,
        "N": n,
        "support_budget": len(cols),
        "support_strategy": support_strategy,
        "encoding": encoding,
        "include_prime_harness_setup_cost": include_harness,
        "columns": cols,
        "address_bits": address_bits,
        "dictionary_header_bits": spec.header_bits,
        "stage_header_bits": STAGE_HEADER_BITS,
        "coefficient_bits": fit["active_count"] * PRIMARY_BITS,
        "total_description_bits": bits,
        "description_bits_per_information_bit": (
            bits / info_bits if info_bits and info_bits > 0 else None
        ),
        **fit,
    }


def support_budgets(n: int) -> list[int]:
    values = sorted({m for m in SUPPORT_BUDGETS if m <= n} | {n})
    return [m for m in values if m > 0]


def integer_prefix(m: int) -> list[int]:
    return list(range(1, m + 1))


def composite_prefix(n: int, m: int) -> list[int]:
    cols = [1]
    for j in range(2, n + 1):
        if not is_prime(j):
            cols.append(j)
        if len(cols) >= m:
            break
    return cols


def factor_cost_ordered(n: int, m: int) -> list[int]:
    costs = [(factor_address_bits(j, n), j) for j in range(1, n + 1)]
    return sorted(j for _cost, j in sorted(costs)[:m])


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
    cols_zero = [j - 1 for j in cols]
    fit = fit_columns(A, y, cols_zero)

    if support_strategy == "integer_prefix":
        rows.append(
            row_for_encoding(
                spec,
                n,
                grid,
                support_strategy,
                "flat_index",
                cols,
                fit,
                flat_address_bits(cols, n),
                False,
            )
        )
        rows.append(
            row_for_encoding(
                spec,
                n,
                grid,
                support_strategy,
                "factorized_reusable_harness",
                cols,
                fit,
                factorized_address_bits(cols, n, include_harness=False),
                False,
            )
        )
        rows.append(
            row_for_encoding(
                spec,
                n,
                grid,
                support_strategy,
                "factorized_with_harness",
                cols,
                fit,
                factorized_address_bits(cols, n, include_harness=True),
                True,
            )
        )
    elif support_strategy == "composite_prefix":
        bits = composite_address_bits(cols, n)
        if bits is not None:
            rows.append(
                row_for_encoding(
                    spec,
                    n,
                    grid,
                    support_strategy,
                    "composite_index",
                    cols,
                    fit,
                    bits,
                    False,
                )
            )
    elif support_strategy == "factor_cost_ordered":
        rows.append(
            row_for_encoding(
                spec,
                n,
                grid,
                support_strategy,
                "factorized_reusable_harness",
                cols,
                fit,
                factorized_address_bits(cols, n, include_harness=False),
                False,
            )
        )
        rows.append(
            row_for_encoding(
                spec,
                n,
                grid,
                support_strategy,
                "factorized_with_harness",
                cols,
                fit,
                factorized_address_bits(cols, n, include_harness=True),
                True,
            )
        )
    else:
        raise ValueError(f"unknown support strategy {support_strategy}")


def tolerance_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, int, float], list[dict[str, Any]]] = {}
    for row in rows:
        for tol in RESIDUAL_TOLERANCES:
            key = (row["grid"], row["dictionary"], row["N"], tol)
            groups.setdefault(key, []).append(row)

    table = []
    for (grid, dictionary, n, tol), group_rows in sorted(groups.items()):
        winners = []
        by_encoding = sorted({r["encoding"] for r in group_rows})
        for encoding in by_encoding:
            candidates = [
                r
                for r in group_rows
                if r["encoding"] == encoding
                and r["certified_upper_relative_residual_16"] <= tol
            ]
            if candidates:
                best = min(candidates, key=lambda r: r["total_description_bits"])
                winners.append(
                    {
                        "encoding": encoding,
                        "support_strategy": best["support_strategy"],
                        "support_budget": best["support_budget"],
                        "total_description_bits": best["total_description_bits"],
                        "certified_upper_relative_residual_16": best[
                            "certified_upper_relative_residual_16"
                        ],
                        "description_bits_per_information_bit": best[
                            "description_bits_per_information_bit"
                        ],
                    }
                )
        if winners:
            overall = min(winners, key=lambda r: r["total_description_bits"])
            table.append(
                {
                    "grid": grid,
                    "dictionary": dictionary,
                    "N": n,
                    "tolerance": tol,
                    "winner_encoding": overall["encoding"],
                    "winner_support_strategy": overall["support_strategy"],
                    "winner_total_description_bits": overall["total_description_bits"],
                    "encoding_best": winners,
                }
            )
    return table


def summarize(rows: list[dict[str, Any]], tolerances: list[dict[str, Any]]) -> dict[str, Any]:
    same_support_pairs: dict[tuple[str, str, int, int], dict[str, dict[str, Any]]] = {}
    for row in rows:
        if row["support_strategy"] != "integer_prefix":
            continue
        key = (row["grid"], row["dictionary"], row["N"], row["support_budget"])
        same_support_pairs.setdefault(key, {})[row["encoding"]] = row

    reusable_wins = 0
    charged_wins = 0
    comparable = 0
    reusable_savings = []
    charged_savings = []
    for encs in same_support_pairs.values():
        flat = encs.get("flat_index")
        reusable = encs.get("factorized_reusable_harness")
        charged = encs.get("factorized_with_harness")
        if flat and reusable and charged:
            comparable += 1
            reusable_delta = flat["total_description_bits"] - reusable["total_description_bits"]
            charged_delta = flat["total_description_bits"] - charged["total_description_bits"]
            reusable_savings.append(reusable_delta)
            charged_savings.append(charged_delta)
            if reusable_delta > 0:
                reusable_wins += 1
            if charged_delta > 0:
                charged_wins += 1

    winner_counts: dict[str, int] = {}
    for row in tolerances:
        winner_counts[row["winner_encoding"]] = winner_counts.get(row["winner_encoding"], 0) + 1

    return {
        "row_count": len(rows),
        "tolerance_row_count": len(tolerances),
        "same_support_comparable_count": comparable,
        "same_support_factorized_reusable_win_count": reusable_wins,
        "same_support_factorized_reusable_win_fraction": reusable_wins / comparable if comparable else None,
        "same_support_factorized_with_harness_win_count": charged_wins,
        "same_support_factorized_with_harness_win_fraction": charged_wins / comparable if comparable else None,
        "median_reusable_address_bit_savings_vs_flat": (
            float(np.median(reusable_savings)) if reusable_savings else None
        ),
        "median_charged_address_bit_savings_vs_flat": (
            float(np.median(charged_savings)) if charged_savings else None
        ),
        "tolerance_winner_counts": winner_counts,
        "status": classify_status(
            comparable,
            reusable_wins,
            charged_wins,
            winner_counts,
        ),
    }


def classify_status(
    comparable: int,
    reusable_wins: int,
    charged_wins: int,
    winner_counts: dict[str, int],
) -> str:
    if comparable == 0:
        return "NO_COMPARABLE_ROWS"
    reusable_fraction = reusable_wins / comparable
    charged_fraction = charged_wins / comparable
    tolerance_total = sum(winner_counts.values())
    tolerance_factor_wins = (
        winner_counts.get("factorized_reusable_harness", 0)
        + winner_counts.get("factorized_with_harness", 0)
    )
    tolerance_fraction = tolerance_factor_wins / tolerance_total if tolerance_total else 0.0
    if reusable_fraction >= 0.60 and tolerance_fraction >= 0.50:
        if charged_fraction >= 0.50:
            return "PRIME_HARNESS_SIGNAL_WITH_SETUP_COST"
        return "PRIME_HARNESS_SIGNAL_REUSABLE_ONLY"
    if reusable_fraction >= 0.60:
        return "PRIME_HARNESS_SAME_SUPPORT_ONLY"
    return "NO_PRIME_HARNESS_MDL_SIGNAL"


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
    summary = summarize(rows, tolerances)
    return {
        "experiment_id": EXPERIMENT_ID,
        "date": date.today().isoformat(),
        "claim_ceiling": CLAIM_CEILING,
        "model_note": "PRIME_HARNESS_MDL_MODEL_2026-05-06.md",
        "source_reversal": "EXP-MATH-RH-BN-COMPOSITE-FIRST-PROJECTION-20260506-01",
        "encoding_models": {
            "flat_index": "Direct index code, ceil(log2(N + 1)) bits per selected dictionary index.",
            "composite_index": "Unit plus composite-rank code; prime indices are not representable.",
            "factorized_reusable_harness": "Prime-factorization address code with the prime harness treated as reusable setup.",
            "factorized_with_harness": "Prime-factorization address code charging the prime harness once per row.",
        },
        "grids": GRID_SPECS,
        "n_values": list(n_values),
        "support_budgets": SUPPORT_BUDGETS,
        "residual_tolerances": RESIDUAL_TOLERANCES,
        "primary_coeff_bits": PRIMARY_BITS,
        "summary": summary,
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
        f"- Same-support comparable rows: `{summary['same_support_comparable_count']}`",
        (
            "- Factorized reusable wins on same support: "
            f"`{summary['same_support_factorized_reusable_win_count']}`"
        ),
        (
            "- Factorized with-harness wins on same support: "
            f"`{summary['same_support_factorized_with_harness_win_count']}`"
        ),
        (
            "- Median reusable savings vs flat direct index: "
            f"`{summary['median_reusable_address_bit_savings_vs_flat']}` bits"
        ),
        (
            "- Median charged savings vs flat direct index: "
            f"`{summary['median_charged_address_bit_savings_vs_flat']}` bits"
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
        "This experiment tests the revised prime-harness interpretation after the",
        "composite-first projection diagnostic reversed. The question is whether",
        "prime factorization works as a cheaper finite addressing layer for",
        "Beurling-Nyman dictionary indices than flat integer addressing.",
        "",
        "The result is an encoding-cost diagnostic. It is not a theorem about",
        "primes, not a zeta formalization, and not a quantum-mechanical claim.",
        "",
        "## Claim Ceiling",
        "",
        results["claim_ceiling"],
        "",
        "## Cost Models",
        "",
        "- `flat_index`: direct index address for selected dictionary columns.",
        "- `composite_index`: unit plus composite-rank addressing only.",
        "- `factorized_reusable_harness`: factorization addressing after the prime harness is already shared.",
        "- `factorized_with_harness`: factorization addressing plus setup cost for the prime harness.",
        "",
        "## Tolerance Winners",
        "",
        "| grid | dictionary | N | tolerance | winner | support | bits |",
        "|---|---|---:|---:|---|---|---:|",
    ]
    for row in results["tolerance_table"][:40]:
        lines.append(
            "| {grid} | {dictionary} | {N} | {tolerance:.2f} | {winner_encoding} | "
            "{winner_support_strategy} | {winner_total_description_bits} |".format(**row)
        )
    if len(results["tolerance_table"]) > 40:
        lines.append(f"| ... | ... | ... | ... | {len(results['tolerance_table']) - 40} more rows | ... | ... |")
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "If the reusable factorized code wins but the charged code loses, the",
            "finite message is amortization: the prime-address layer helps only",
            "when the harness is already part of the shared dictionary protocol.",
            "That is useful MDL structure, but not a free-description claim.",
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
    digest = sha256_file(json_path)
    sha_path.write_text(f"{digest}  {json_path.name}\n")


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
