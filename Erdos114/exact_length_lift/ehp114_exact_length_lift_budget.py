#!/usr/bin/env python3
"""Budget prototype for lifting EHP114 n=14 SUBDIV8 to exact length.

This script reads the current root-affine SUBDIV8 certificate and computes the
remaining additive/relative slack available to a future exact or validated
lemniscate-length enclosure. It does not compute exact lemniscate length.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from statistics import mean


EXPERIMENT_ID = "EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01"
MATH_ROOT = Path(__file__).resolve().parents[3]
EHP114_ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    EHP114_ROOT
    / "low_dim_cone"
    / "EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01_RESULTS.json"
)
OUT_DIR = Path(__file__).resolve().parent
RESULTS_PATH = OUT_DIR / f"{EXPERIMENT_ID}_RESULTS.json"
REPORT_PATH = OUT_DIR / f"{EXPERIMENT_ID}_REPORT.md"
SHA_PATH = OUT_DIR / f"{EXPERIMENT_ID}_RESULTS.sha256"


def load_source() -> dict:
    with SOURCE.open("r", encoding="utf-8") as f:
        return json.load(f)


def build_result(source: dict) -> dict:
    rows = source["rows"]
    eps = float(source["parameters"]["eps"])
    target = 12.0 * math.pow(eps, 1.0 / 14.0)
    lstar_lower = float(source["lstar_lower"])
    exact_length_cap = lstar_lower - target

    enriched = []
    for row in rows:
        length_upper = float(row["length_upper"])
        margin_lower = float(row["margin_lower"])
        inferred_cap = length_upper + margin_lower
        enriched.append(
            {
                "sub_i": row["sub_i"],
                "sub_j": row["sub_j"],
                "u0_interval": row["u0_interval"],
                "u1_interval": row["u1_interval"],
                "marching_length_upper": length_upper,
                "margin_lower": margin_lower,
                "exact_length_upper_cap": inferred_cap,
                "allowed_additive_exact_length_error": margin_lower,
                "allowed_relative_exact_length_error_vs_marching": margin_lower
                / length_upper,
                "active_cells": row["active_cells"],
                "uncertain_corner_cells": row["uncertain_corner_cells"],
            }
        )

    worst_additive = min(enriched, key=lambda r: r["allowed_additive_exact_length_error"])
    worst_relative = min(
        enriched, key=lambda r: r["allowed_relative_exact_length_error_vs_marching"]
    )

    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": "BUDGET_ONLY_NOT_EXACT_LENGTH_CERTIFICATE",
        "claim_ceiling": (
            "Prototype budget extraction from the current marching-squares "
            "SUBDIV8 certificate. It does not compute exact lemniscate length, "
            "does not validate the marching-squares functional against exact "
            "length, and does not prove Erdos #114."
        ),
        "source_experiment_id": source["experiment_id"],
        "source_path": str(SOURCE.relative_to(MATH_ROOT)),
        "source_claim_ceiling": source.get("claim_ceiling"),
        "eps": eps,
        "lstar_lower": lstar_lower,
        "target_reserve": target,
        "exact_length_upper_cap_uniform": exact_length_cap,
        "subcell_count": len(rows),
        "source_failure_count": source["failure_count"],
        "source_max_marching_length_upper": source["max_length_upper"],
        "source_min_margin_lower": source["min_margin_lower"],
        "additive_error_budget": {
            "min": min(r["allowed_additive_exact_length_error"] for r in enriched),
            "max": max(r["allowed_additive_exact_length_error"] for r in enriched),
            "mean": mean(r["allowed_additive_exact_length_error"] for r in enriched),
            "worst_cell": worst_additive,
        },
        "relative_error_budget_vs_marching": {
            "min": min(
                r["allowed_relative_exact_length_error_vs_marching"] for r in enriched
            ),
            "max": max(
                r["allowed_relative_exact_length_error_vs_marching"] for r in enriched
            ),
            "mean": mean(
                r["allowed_relative_exact_length_error_vs_marching"] for r in enriched
            ),
            "worst_cell": worst_relative,
        },
        "sufficient_bridge_condition": (
            "For every SUBDIV8 subcell C, prove a validated exact lemniscate "
            "length upper enclosure L_exact(C) <= L_ms(C) + E(C), with "
            "E(C) <= margin_lower(C). A uniform additive bridge with "
            f"E <= {worst_additive['allowed_additive_exact_length_error']:.15g} "
            "would preserve the current scalar reserve on the whole selected cell."
        ),
        "coarea_bridge_inputs_still_missing": [
            "a lower bound on |grad(|p|)| or |p'| on the relevant level-set collar",
            "a validated collar width tau around |p| = 1",
            "an interval enclosure for the area of the collar or a direct implicit-curve length enclosure",
            "a proof that the numerical contour extraction overcounts exact length by at most E(C)",
        ],
        "rows": enriched,
    }
    return result


def write_report(result: dict) -> None:
    add = result["additive_error_budget"]
    rel = result["relative_error_budget_vs_marching"]
    worst = add["worst_cell"]
    lines = [
        f"# {EXPERIMENT_ID} Report",
        "",
        "## Meaning",
        "",
        "This is a budget prototype for the exact-length lift after the Python SUBDIV8 pass.",
        "It does not compute exact lemniscate length and does not certify the marching-squares oracle as proof-grade length.",
        "",
        "## Source",
        "",
        f"- Source experiment: `{result['source_experiment_id']}`",
        f"- Source path: `{result['source_path']}`",
        f"- Source failures: `{result['source_failure_count']}`",
        "",
        "## Budget",
        "",
        f"- `eps`: `{result['eps']}`",
        f"- `lstar_lower`: `{result['lstar_lower']}`",
        f"- `target_reserve`: `{result['target_reserve']}`",
        f"- uniform exact length cap: `{result['exact_length_upper_cap_uniform']}`",
        f"- source max marching length upper: `{result['source_max_marching_length_upper']}`",
        f"- minimum additive exact-length error budget: `{add['min']}`",
        f"- mean additive exact-length error budget: `{add['mean']}`",
        f"- minimum relative error budget vs marching upper: `{rel['min']}`",
        "",
        "Worst additive-budget cell:",
        "",
        f"- subcell: `({worst['sub_i']}, {worst['sub_j']})`",
        f"- `u0_interval`: `{worst['u0_interval']}`",
        f"- `u1_interval`: `{worst['u1_interval']}`",
        f"- marching length upper: `{worst['marching_length_upper']}`",
        f"- allowed additive exact-length error: `{worst['allowed_additive_exact_length_error']}`",
        "",
        "## Sufficient Bridge Condition",
        "",
        result["sufficient_bridge_condition"],
        "",
        "## Missing For Proof Grade",
        "",
    ]
    for item in result["coarea_bridge_inputs_still_missing"]:
        lines.append(f"- {item}")
    lines.extend(
        [
            "",
            "Claim ceiling: budget-only prototype; not exact lemniscate certification, not a Lean theorem, and not a proof of Erdos #114.",
            "",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    source = load_source()
    result = build_result(source)
    RESULTS_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(result)
    digest = hashlib.sha256(RESULTS_PATH.read_bytes()).hexdigest()
    SHA_PATH.write_text(f"{digest}  {RESULTS_PATH.name}\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ["experiment_id", "status"]}, indent=2))
    print(f"wrote {RESULTS_PATH}")
    print(f"wrote {REPORT_PATH}")
    print(f"wrote {SHA_PATH}")


if __name__ == "__main__":
    main()
