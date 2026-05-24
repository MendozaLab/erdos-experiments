#!/usr/bin/env python3
"""Derived n=58 branch certificate for the Erdos #30 face export.

This script derives finite checks from an existing exact packet. It does not
run a new Sidon search or widen the enumeration window.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Iterable


EXPERIMENT_ID = "EXP-MM-030-PMF-N58-BRANCH-CERTIFICATE-2026-04-30"
SOURCE_ID = "EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30"


def erdos_experiments_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def positive_differences(values: list[int]) -> list[int]:
    return sorted(values[i] - values[j] for i in range(len(values)) for j in range(i))


def is_interval_sidon(values: list[int]) -> bool:
    diffs = positive_differences(values)
    return len(diffs) == len(set(diffs))


def canonical(values: Iterable[int]) -> str:
    return json.dumps(list(values))


def set_distance(left: Iterable[int], right: Iterable[int]) -> int:
    a = set(left)
    b = set(right)
    return len(a - b) + len(b - a)


def translate(values: list[int], shift: int) -> list[int]:
    return [x + shift for x in values]


def translate_mod(values: list[int], shift: int, modulus: int) -> list[int]:
    return sorted((x + shift) % modulus for x in values)


def reflect_about(values: list[int], center: int) -> list[int]:
    return sorted(center - x for x in values)


def reflect_mod(values: list[int], center: int, modulus: int) -> list[int]:
    return sorted((center - x) % modulus for x in values)


def unit_group(modulus: int) -> list[int]:
    return [a for a in range(modulus) if math.gcd(a, modulus) == 1]


def affine_mod(values: list[int], multiplier: int, translate_by: int, modulus: int) -> list[int]:
    return sorted(((multiplier * x) + translate_by) % modulus for x in values)


def find_relations(
    left: list[int],
    right: list[int],
    *,
    moduli: list[int],
    reflection_centers: list[int],
    translation_span: int,
) -> dict:
    relations: dict[str, list[dict]] = {
        "integer_translations": [],
        "integer_reflections": [],
        "mod_translations": [],
        "mod_reflections": [],
        "mod_affine_maps": [],
    }
    right_key = canonical(right)
    for shift in range(-translation_span, translation_span + 1):
        if shift and canonical(translate(left, shift)) == right_key:
            relations["integer_translations"].append({"shift": shift})
    for center in reflection_centers:
        if canonical(reflect_about(left, center)) == right_key:
            relations["integer_reflections"].append({"center": center})
    for modulus in moduli:
        right_mod_key = canonical(sorted(x % modulus for x in right))
        for shift in range(modulus):
            if canonical(translate_mod(left, shift, modulus)) == right_mod_key:
                relations["mod_translations"].append({"modulus": modulus, "shift": shift})
        for center in range(modulus):
            if canonical(reflect_mod(left, center, modulus)) == right_mod_key:
                relations["mod_reflections"].append({"modulus": modulus, "center": center})
        for multiplier in unit_group(modulus):
            for translate_by in range(modulus):
                if canonical(affine_mod(left, multiplier, translate_by, modulus)) == right_mod_key:
                    relations["mod_affine_maps"].append(
                        {
                            "modulus": modulus,
                            "multiplier": multiplier,
                            "translate": translate_by,
                        }
                    )
    return relations


def classify_skeletons(witnesses: list[dict]) -> tuple[list[dict], list[dict]]:
    skeleton_by_key: dict[str, int] = {}
    skeletons: list[dict] = []
    classified: list[dict] = []
    for witness in witnesses:
        diffs = positive_differences(witness["witness"])
        key = canonical(diffs)
        if key not in skeleton_by_key:
            skeleton_by_key[key] = len(skeleton_by_key)
            skeletons.append(
                {
                    "skeleton_id": skeleton_by_key[key],
                    "positive_differences": diffs,
                    "witness_indices": [],
                }
            )
        skeleton_id = skeleton_by_key[key]
        skeletons[skeleton_id]["witness_indices"].append(witness["index"])
        classified.append({**witness, "skeleton_id": skeleton_id, "is_interval_sidon": is_interval_sidon(witness["witness"])})
    return classified, skeletons


def row_by_n(raw: dict, n: int) -> dict:
    return next(row for row in raw["rows"] if row["n"] == n)


def write_sha256(results_path: Path) -> Path:
    digest = hashlib.sha256(results_path.read_bytes()).hexdigest()
    sha_path = results_path.with_name(results_path.name.replace("_RESULTS.json", "_RESULTS.sha256"))
    sha_path.write_text(f"{digest}  {results_path.name}\n")
    return sha_path


def main() -> int:
    erdos_root = erdos_experiments_root_from_script()
    default_source = erdos_root / "results" / "erdos-30" / f"{SOURCE_ID}_RESULTS.json"
    default_out = erdos_root / "results" / "erdos-30"

    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=default_source)
    parser.add_argument("--out-dir", type=Path, default=default_out)
    args = parser.parse_args()

    source_path = args.source.resolve()
    out_dir = args.out_dir.resolve()
    raw = json.loads(source_path.read_text())
    row57 = row_by_n(raw, 57)
    row58 = row_by_n(raw, 58)
    w57 = row57["ground_face_export"]["witnesses"]
    w58 = row58["ground_face_export"]["witnesses"]
    c58, skeletons58 = classify_skeletons(w58)

    by57 = {w["index"]: w for w in w57}
    by58 = {w["index"]: w for w in c58}
    by57_key = {canonical(w["witness"]): w["index"] for w in w57}

    n58_pareto = row58["ground_face_export"]["exposed_winners"]["pareto_minimal_indices"]
    n58_prefix = row58["ground_face_export"]["exposed_winners"]["prefix_winner_indices"]
    n58_mass = row58["ground_face_export"]["exposed_winners"]["mass_winner_indices"]
    n58_joint = row58["ground_face_export"]["exposed_winners"]["joint_winner_indices"]

    persistence_from_57 = []
    for witness in c58:
        source_index = by57_key.get(canonical(witness["witness"]))
        if source_index is not None:
            persistence_from_57.append({"n57_index": source_index, "n58_index": witness["index"]})

    pair_7_9_relations = find_relations(
        by58[7]["witness"],
        by58[9]["witness"],
        moduli=[57, 58],
        reflection_centers=[57, 58],
        translation_span=5,
    )
    n57_5_to_n58_7 = canonical(by57[5]["witness"]) == canonical(by58[7]["witness"])
    n57_5_to_n58_9_relations = find_relations(
        by57[5]["witness"],
        by58[9]["witness"],
        moduli=[57, 58],
        reflection_centers=[57, 58],
        translation_span=5,
    )
    n57_3_to_n58_9_relations = find_relations(
        by57[3]["witness"],
        by58[9]["witness"],
        moduli=[57, 58],
        reflection_centers=[57, 58],
        translation_span=5,
    )

    pairwise_skeleton_distances = []
    for i, left in enumerate(c58):
        for right in c58[i + 1 :]:
            pairwise_skeleton_distances.append(
                {
                    "left_index": left["index"],
                    "right_index": right["index"],
                    "left_skeleton": left["skeleton_id"],
                    "right_skeleton": right["skeleton_id"],
                    "site_symmetric_difference": set_distance(left["witness"], right["witness"]),
                    "difference_skeleton_symmetric_difference": set_distance(
                        positive_differences(left["witness"]),
                        positive_differences(right["witness"]),
                    ),
                }
            )

    skeleton_roles = []
    for skeleton in skeletons58:
        indices = set(skeleton["witness_indices"])
        skeleton_roles.append(
            {
                "skeleton_id": skeleton["skeleton_id"],
                "witness_indices": skeleton["witness_indices"],
                "contains_pareto": sorted(indices & set(n58_pareto)),
                "contains_prefix": sorted(indices & set(n58_prefix)),
                "contains_mass": sorted(indices & set(n58_mass)),
                "contains_joint": sorted(indices & set(n58_joint)),
                "persists_from_n57": sorted(
                    p["n58_index"] for p in persistence_from_57 if p["n58_index"] in indices
                ),
            }
        )

    branch_skeleton_ids = sorted(
        role["skeleton_id"] for role in skeleton_roles if not role["persists_from_n57"]
    )
    pareto_skeleton_ids = sorted({by58[index]["skeleton_id"] for index in n58_pareto})
    prefix_skeleton_ids = sorted({by58[index]["skeleton_id"] for index in n58_prefix})

    checks = {
        "all_n58_exported_witnesses_are_interval_sidon": all(w["is_interval_sidon"] for w in c58),
        "n58_has_two_difference_skeletons": len(skeletons58) == 2,
        "pareto_7_and_9_share_skeleton": by58[7]["skeleton_id"] == by58[9]["skeleton_id"],
        "pareto_9_is_plus_one_translation_of_pareto_7": bool(pair_7_9_relations["integer_translations"])
        and pair_7_9_relations["integer_translations"][0]["shift"] == 1,
        "n58_index7_persists_from_n57_index5": n57_5_to_n58_7,
        "n58_index9_is_n57_index5_plus_one": bool(n57_5_to_n58_9_relations["integer_translations"])
        and n57_5_to_n58_9_relations["integer_translations"][0]["shift"] == 1,
        "new_skeleton_branch_exists_at_n58": bool(branch_skeleton_ids),
        "new_skeleton_branch_is_prefix_side_not_pareto_side": (
            bool(set(prefix_skeleton_ids) & set(branch_skeleton_ids))
            and not bool(set(pareto_skeleton_ids) & set(branch_skeleton_ids))
        ),
    }

    result = {
        "experiment_id": EXPERIMENT_ID,
        "date": "2026-04-30",
        "type": "derived_n58_branch_certificate_from_exact_packet",
        "source_packet": {
            "experiment_id": SOURCE_ID,
            "results_path": str(source_path),
            "status": "EXACT",
        },
        "n": 58,
        "h_n": row58["h_n"],
        "exact_maximizer_count": row58["maximizer_count"],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "n58_winners": {
            "prefix_winner_indices": n58_prefix,
            "mass_winner_indices": n58_mass,
            "joint_winner_indices": n58_joint,
            "pareto_minimal_indices": n58_pareto,
        },
        "n58_skeletons": skeletons58,
        "n58_skeleton_roles": skeleton_roles,
        "n58_witnesses": c58,
        "persistence_from_n57": persistence_from_57,
        "relations": {
            "n58_7_to_n58_9": pair_7_9_relations,
            "n57_5_to_n58_7_exact_persistence": n57_5_to_n58_7,
            "n57_5_to_n58_9": n57_5_to_n58_9_relations,
            "n57_3_to_n58_9": n57_3_to_n58_9_relations,
        },
        "pairwise_skeleton_distances": pairwise_skeleton_distances,
        "interpretation": (
            "n=58 is the first local branch row, but the Pareto face remains on "
            "the original translated-embedding skeleton. Candidate 7 persists "
            "from n=57 index 5; candidate 9 is candidate 7 shifted by +1. The "
            "new second skeleton appears on the prefix side through witness index 2."
        ),
        "claim_boundary": (
            "Finite branch certificate only: no asymptotic claim, no #30 proof, "
            "and no upgrade from Singer-modulus lock-on to Singer PDS containment."
        ),
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    results_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    results_path.write_text(json.dumps(result, indent=2) + "\n")
    sha_path = write_sha256(results_path)

    skeleton_table = "\n".join(
        f"| {role['skeleton_id']} | {role['witness_indices']} | {role['persists_from_n57']} | "
        f"{role['contains_prefix']} | {role['contains_mass']} | {role['contains_joint']} | {role['contains_pareto']} |"
        for role in skeleton_roles
    )
    report_path.write_text(
        f"# {EXPERIMENT_ID}\n\n"
        "Date: 2026-04-30\n"
        "Problem: Erdos #30\n"
        "Status: DERIVED_N58_BRANCH_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE\n\n"
        "## Source\n\n"
        f"Derived from `{SOURCE_ID}`. No new enumeration was run.\n\n"
        "## Result\n\n"
        "The `n=58` row is the first local branch row, but the branch is not the "
        "Pareto face. Pareto candidates `7` and `9` share the original skeleton, "
        "and `9 = 7 + 1`. The new second skeleton appears on the prefix side.\n\n"
        "## Checks\n\n"
        f"- all `n=58` exported witnesses are interval Sidon: `{checks['all_n58_exported_witnesses_are_interval_sidon']}`\n"
        f"- `n=58` has two difference skeletons: `{checks['n58_has_two_difference_skeletons']}`\n"
        f"- Pareto `7` and `9` share a skeleton: `{checks['pareto_7_and_9_share_skeleton']}`\n"
        f"- Pareto `9 = 7 + 1`: `{checks['pareto_9_is_plus_one_translation_of_pareto_7']}`\n"
        f"- `n=58` index `7` persists from `n=57` index `5`: `{checks['n58_index7_persists_from_n57_index5']}`\n"
        f"- `n=58` index `9 = n=57 index 5 + 1`: `{checks['n58_index9_is_n57_index5_plus_one']}`\n"
        f"- new skeleton branch exists at `n=58`: `{checks['new_skeleton_branch_exists_at_n58']}`\n"
        f"- new branch is prefix-side, not Pareto-side: `{checks['new_skeleton_branch_is_prefix_side_not_pareto_side']}`\n\n"
        f"All checks pass: `{result['all_checks_pass']}`\n\n"
        "## Skeleton Roles\n\n"
        "| skeleton | witnesses | persists from n=57 | prefix | mass | joint | Pareto |\n"
        "|---:|---|---|---|---|---|---|\n"
        f"{skeleton_table}\n\n"
        "## Three-Stage Mechanism\n\n"
        "```text\n"
        "56: single exposed embedding\n"
        "57: translated-embedding handoff inside one Singer-modulus skeleton\n"
        "58: first local skeleton branch appears, but Pareto remains on translated chain\n"
        "```\n\n"
        "## Claim Boundary\n\n"
        f"{result['claim_boundary']}\n"
    )

    print(json.dumps({"results": str(results_path), "report": str(report_path), "sha256": str(sha_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
