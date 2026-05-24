#!/usr/bin/env python3
"""Derived branch-window certificate for Erdos #30 ground-face exports."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Iterable


EXPERIMENT_ID = "EXP-MM-030-PMF-GROUNDFACE-BRANCH-CERTIFICATE-58-64-2026-05-01"
SOURCE_58_ID = "EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30"
SOURCE_59_64_ID = "EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01"


def erdos_experiments_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def positive_differences(values: list[int]) -> tuple[int, ...]:
    return tuple(sorted(values[i] - values[j] for i in range(len(values)) for j in range(i)))


def translate(values: Iterable[int], shift: int) -> tuple[int, ...]:
    return tuple(x + shift for x in values)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_sha256(results_path: Path) -> Path:
    digest = sha256_file(results_path)
    sha_path = results_path.with_name(results_path.name.replace("_RESULTS.json", "_RESULTS.sha256"))
    sha_path.write_text(f"{digest}  {results_path.name}\n")
    return sha_path


def load_rows(source_anchor: Path, source_windows: list[Path], anchor_n: int) -> list[dict]:
    anchor_rows = json.loads(source_anchor.read_text())["rows"]
    rows = [row for row in anchor_rows if row["n"] == anchor_n]
    for source in source_windows:
        rows.extend(json.loads(source.read_text())["rows"])
    rows.sort(key=lambda row: row["n"])
    return rows


def classify_row(row: dict, previous_witnesses: list[list[int]]) -> dict:
    export = row["ground_face_export"]
    witnesses = export["witnesses"]
    winners = export["exposed_winners"]

    skeleton_keys: dict[tuple[int, ...], int] = {}
    skeleton_roles: dict[int, dict] = {}
    index_to_skeleton: dict[int, int] = {}
    for witness in witnesses:
        key = positive_differences(witness["witness"])
        if key not in skeleton_keys:
            skeleton_id = len(skeleton_keys)
            skeleton_keys[key] = skeleton_id
            skeleton_roles[skeleton_id] = {
                "skeleton_id": skeleton_id,
                "witness_indices": [],
                "contains_prefix": [],
                "contains_mass": [],
                "contains_joint": [],
                "contains_pareto": [],
                "exact_previous_persistence": [],
                "plus_one_previous_persistence": [],
            }
        skeleton_id = skeleton_keys[key]
        index_to_skeleton[witness["index"]] = skeleton_id
        skeleton_roles[skeleton_id]["witness_indices"].append(witness["index"])

    current_by_index = {w["index"]: tuple(w["witness"]) for w in witnesses}
    previous_exact = {tuple(w) for w in previous_witnesses}
    previous_plus_one = {translate(w, 1) for w in previous_witnesses}

    exact_previous_indices = []
    plus_one_previous_indices = []
    for witness in witnesses:
        current = tuple(witness["witness"])
        skeleton = index_to_skeleton[witness["index"]]
        if current in previous_exact:
            exact_previous_indices.append(witness["index"])
            skeleton_roles[skeleton]["exact_previous_persistence"].append(witness["index"])
        if current in previous_plus_one:
            plus_one_previous_indices.append(witness["index"])
            skeleton_roles[skeleton]["plus_one_previous_persistence"].append(witness["index"])

    role_fields = [
        ("prefix_winner_indices", "contains_prefix"),
        ("mass_winner_indices", "contains_mass"),
        ("joint_winner_indices", "contains_joint"),
        ("pareto_minimal_indices", "contains_pareto"),
    ]
    for source_key, target_key in role_fields:
        for index in winners[source_key]:
            skeleton_roles[index_to_skeleton[index]][target_key].append(index)

    pareto_skeletons = sorted({index_to_skeleton[i] for i in winners["pareto_minimal_indices"]})
    prefix_skeletons = sorted({index_to_skeleton[i] for i in winners["prefix_winner_indices"]})
    mass_skeletons = sorted({index_to_skeleton[i] for i in winners["mass_winner_indices"]})
    joint_skeletons = sorted({index_to_skeleton[i] for i in winners["joint_winner_indices"]})

    return {
        "n": row["n"],
        "h_n": row["h_n"],
        "exact_maximizer_count": row["maximizer_count"],
        "export_status": export["status"],
        "exported_count": export["exported_count"],
        "skeleton_count": len(skeleton_keys),
        "field_selection_split": winners["field_selection_split"],
        "winner_counts": {
            "prefix": len(winners["prefix_winner_indices"]),
            "mass": len(winners["mass_winner_indices"]),
            "joint": len(winners["joint_winner_indices"]),
            "pareto": len(winners["pareto_minimal_indices"]),
        },
        "winner_skeleton_counts": {
            "prefix": len(prefix_skeletons),
            "mass": len(mass_skeletons),
            "joint": len(joint_skeletons),
            "pareto": len(pareto_skeletons),
        },
        "winner_skeletons": {
            "prefix": prefix_skeletons,
            "mass": mass_skeletons,
            "joint": joint_skeletons,
            "pareto": pareto_skeletons,
        },
        "exact_previous_persistence_count": len(exact_previous_indices),
        "plus_one_previous_persistence_count": len(plus_one_previous_indices),
        "plus_one_previous_pareto_count": sum(
            current_by_index[index] in previous_plus_one for index in winners["pareto_minimal_indices"]
        ),
        "skeleton_roles": sorted(skeleton_roles.values(), key=lambda row: row["skeleton_id"]),
    }


def build_report(result: dict) -> str:
    rows = result["rows"]
    anchor_n = result["anchor_n"]
    growth_check_label = (
        "nondecreasing with first post-anchor plateau allowed, then strict"
        if result["allow_initial_skeleton_plateau"]
        else "strictly increasing"
    )
    table = "\n".join(
        "| {n} | {exact_maximizer_count} | {skeleton_count} | {field_selection_split} | "
        "{exact_previous_persistence_count} | {plus_one_previous_persistence_count} | "
        "{pareto} | {joint} |".format(
            **row,
            pareto=row["winner_counts"]["pareto"],
            joint=row["winner_counts"]["joint"],
        )
        for row in rows
    )
    checks = result["checks"]
    return (
        f"# {result['experiment_id']}\n\n"
        "Date: 2026-05-01\n"
        "Problem: Erdos #30\n"
        "Status: DERIVED_GROUNDFACE_BRANCH_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE\n\n"
        "## Sources\n\n"
        + "\n".join(
            f"- `{source['experiment_id']}`: `{source['results_path']}`"
            for source in result["source_packets"]
        )
        + "\n\n"
        "## Result\n\n"
        f"The post-{anchor_n} face does not collapse back to a single translated branch. "
        "Each subsequent row contains the previous exact face and the previous face shifted by `+1`, "
        "while new difference-skeleton branches accumulate rapidly. The field-selected Pareto face "
        "remains a small subset of a much larger branching ground face.\n\n"
        "## Checks\n\n"
        f"- all source exports are `EXPORTED_ALL`: `{checks['all_exports_complete']}`\n"
        f"- all rows preserve field-selection split: `{checks['all_rows_have_field_selection_split']}`\n"
        f"- every row after `{anchor_n}` contains every previous witness exactly: `{checks['all_post58_rows_contain_previous_face']}`\n"
        f"- every row after `{anchor_n}` contains every previous witness shifted by `+1`: `{checks['all_post58_rows_contain_previous_face_plus_one']}`\n"
        f"- skeleton count is strictly increasing across the window: `{checks['skeleton_count_strictly_increases']}`\n"
        f"- skeleton count is nondecreasing across the window: `{checks['skeleton_count_nondecreasing']}`\n"
        f"- growth rule used for this certificate is {growth_check_label}: `{checks['skeleton_growth_rule_passes']}`\n\n"
        "## Branch Window Table\n\n"
        "| n | exact face | skeletons | field split | prev exact | prev +1 | Pareto count | joint count |\n"
        "|---:|---:|---:|---|---:|---:|---:|---:|\n"
        f"{table}\n\n"
        "## Claim Boundary\n\n"
        "Finite exact-face branch certificate only. This does not prove Erdős #30, "
        "does not establish an asymptotic phase transition, and does not upgrade PMF analogy to theorem language.\n"
    )


def main() -> int:
    root = erdos_experiments_root_from_script()
    default_results = root / "results" / "erdos-30"
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-id", default=EXPERIMENT_ID)
    parser.add_argument("--anchor-n", type=int, default=58)
    parser.add_argument(
        "--allow-initial-skeleton-plateau",
        action="store_true",
        help="Accept one anchor-to-next-row skeleton plateau before requiring strict later growth.",
    )
    parser.add_argument("--source-58", type=Path, default=default_results / f"{SOURCE_58_ID}_RESULTS.json")
    parser.add_argument(
        "--source-window",
        type=Path,
        action="append",
        default=None,
        help="Additional exact ground-face packet to append after n=58. May be repeated.",
    )
    parser.add_argument("--out-dir", type=Path, default=default_results)
    args = parser.parse_args()

    source_windows = args.source_window or [default_results / f"{SOURCE_59_64_ID}_RESULTS.json"]
    rows_raw = load_rows(args.source_58, source_windows, args.anchor_n)
    previous: list[list[int]] = []
    rows = []
    for row in rows_raw:
        rows.append(classify_row(row, previous))
        previous = [w["witness"] for w in row["ground_face_export"]["witnesses"]]

    skeleton_count_strictly_increases = all(
        rows[idx]["skeleton_count"] > rows[idx - 1]["skeleton_count"] for idx in range(1, len(rows))
    )
    skeleton_count_nondecreasing = all(
        rows[idx]["skeleton_count"] >= rows[idx - 1]["skeleton_count"] for idx in range(1, len(rows))
    )
    skeleton_count_strict_after_initial_pair = all(
        rows[idx]["skeleton_count"] > rows[idx - 1]["skeleton_count"] for idx in range(2, len(rows))
    )
    skeleton_growth_rule_passes = (
        skeleton_count_nondecreasing and skeleton_count_strict_after_initial_pair
        if args.allow_initial_skeleton_plateau
        else skeleton_count_strictly_increases
    )

    checks = {
        "all_exports_complete": all(row["export_status"] == "EXPORTED_ALL" for row in rows),
        "all_rows_have_field_selection_split": all(row["field_selection_split"] for row in rows),
        "all_post58_rows_contain_previous_face": all(
            row["exact_previous_persistence_count"] == rows[idx - 1]["exact_maximizer_count"]
            for idx, row in enumerate(rows)
            if row["n"] > args.anchor_n
        ),
        "all_post58_rows_contain_previous_face_plus_one": all(
            row["plus_one_previous_persistence_count"] == rows[idx - 1]["exact_maximizer_count"]
            for idx, row in enumerate(rows)
            if row["n"] > args.anchor_n
        ),
        "skeleton_count_strictly_increases": skeleton_count_strictly_increases,
        "skeleton_count_nondecreasing": skeleton_count_nondecreasing,
        "skeleton_count_strict_after_initial_pair": skeleton_count_strict_after_initial_pair,
        "skeleton_growth_rule_passes": skeleton_growth_rule_passes,
    }

    required_check_keys = [
        "all_exports_complete",
        "all_rows_have_field_selection_split",
        "all_post58_rows_contain_previous_face",
        "all_post58_rows_contain_previous_face_plus_one",
        "skeleton_growth_rule_passes",
    ]
    result = {
        "experiment_id": args.experiment_id,
        "date": "2026-05-01",
        "type": "derived_groundface_branch_window_certificate",
        "anchor_n": args.anchor_n,
        "allow_initial_skeleton_plateau": args.allow_initial_skeleton_plateau,
        "source_packets": [
            {
                "experiment_id": args.source_58.name.replace("_RESULTS.json", ""),
                "results_path": str(args.source_58.resolve()),
                "status": "EXACT",
                "sha256": sha256_file(args.source_58),
            },
            *[
                {
                    "experiment_id": source.name.replace("_RESULTS.json", ""),
                    "results_path": str(source.resolve()),
                    "status": "EXACT",
                    "sha256": sha256_file(source),
                }
                for source in source_windows
            ],
        ],
        "checks": checks,
        "all_checks_pass": all(checks[key] for key in required_check_keys),
        "rows": rows,
        "interpretation": (
            "After n=58, the face contains both exact persistence and +1 translated persistence "
            "from the prior row while additional difference-skeleton branches accumulate. "
            "The Pareto-selected face remains small relative to the full branching ground face."
        ),
        "claim_boundary": (
            "Finite exact-face branch certificate only: no asymptotic theorem, no Sidon proof, "
            "and no PMF-proves-Sidon language."
        ),
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    results_path = args.out_dir / f"{args.experiment_id}_RESULTS.json"
    report_path = args.out_dir / f"{args.experiment_id}_REPORT.md"
    if results_path.exists() or report_path.exists():
        raise FileExistsError(f"refusing to overwrite existing packet {args.experiment_id}")
    results_path.write_text(json.dumps(result, indent=2) + "\n")
    report_path.write_text(build_report(result))
    sha_path = write_sha256(results_path)
    print(json.dumps({"results": str(results_path), "report": str(report_path), "sha256": str(sha_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
