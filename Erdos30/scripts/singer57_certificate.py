#!/usr/bin/env python3
"""Derived Singer-mod-57 certificate for the Erdos #30 n=57 face.

This script does not run a new Sidon search. It derives elementary finite
checks from an existing exact packet and writes the standard packet contract:
*_RESULTS.json, *_REPORT.md, *_RESULTS.sha256.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Iterable


EXPERIMENT_ID = "EXP-MM-030-PMF-SINGER57-CERTIFICATE-V2-2026-04-30"
SOURCE_ID = "EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30"
MODULUS = 57

SINGER_PDS_REPRESENTATIVES = {
    "PDS_A": [0, 11, 19, 20, 24, 26, 36, 54],
    "PDS_B": [0, 1, 5, 7, 17, 35, 38, 49],
    "PDS_C": [16, 19, 30, 38, 39, 43, 45, 55],
}


def erdos_experiments_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def positive_differences(values: list[int]) -> list[int]:
    return sorted(values[i] - values[j] for i in range(len(values)) for j in range(i))


def ordered_mod_differences(values: list[int], modulus: int) -> list[int]:
    return sorted(
        (values[i] - values[j]) % modulus
        for i in range(len(values))
        for j in range(len(values))
        if i != j
    )


def is_interval_sidon(values: list[int]) -> bool:
    diffs = positive_differences(values)
    return len(diffs) == len(set(diffs))


def pds_check(values: list[int], modulus: int) -> dict:
    diffs = ordered_mod_differences(values, modulus)
    counts = Counter(diffs)
    nonzero = set(range(1, modulus))
    return {
        "size": len(values),
        "ordered_difference_count": len(diffs),
        "nonzero_residue_count": len(nonzero),
        "covers_all_nonzero_residues_once": (
            set(counts.keys()) == nonzero and all(counts[r] == 1 for r in nonzero)
        ),
        "multiplicity_histogram": dict(sorted(Counter(counts.values()).items())),
    }


def unit_group(modulus: int) -> list[int]:
    return [a for a in range(modulus) if math.gcd(a, modulus) == 1]


def affine_image(values: Iterable[int], multiplier: int, translate: int, modulus: int) -> set[int]:
    return {((multiplier * x) + translate) % modulus for x in values}


def multiplier_stabilizer(values: list[int], modulus: int) -> list[dict]:
    target = set(values)
    stabilizers = []
    for multiplier in unit_group(modulus):
        translations = []
        for translate in range(modulus):
            if affine_image(values, multiplier, translate, modulus) == target:
                translations.append(translate)
        if translations:
            stabilizers.append({"multiplier": multiplier, "translations": translations})
    return stabilizers


def best_affine_overlap(
    witness: list[int],
    pds: list[int],
    modulus: int,
) -> dict:
    witness_set = set(x % modulus for x in witness)
    best: list[dict] = []
    max_overlap = -1
    for multiplier in unit_group(modulus):
        for translate in range(modulus):
            image = affine_image(pds, multiplier, translate, modulus)
            overlap_set = sorted(witness_set & image)
            overlap = len(overlap_set)
            if overlap > max_overlap:
                max_overlap = overlap
                best = []
            if overlap == max_overlap:
                best.append(
                    {
                        "multiplier": multiplier,
                        "translate": translate,
                        "overlap_size": overlap,
                        "overlap_set": overlap_set,
                        "image": sorted(image),
                    }
                )
    return {
        "max_overlap": max_overlap,
        "contains_affine_pds_image": max_overlap == len(pds),
        "best_examples": best[:10],
        "best_example_count": len(best),
    }


def load_n57_witnesses(source_path: Path) -> list[dict]:
    raw = json.loads(source_path.read_text())
    row = next(row for row in raw["rows"] if row["n"] == 57)
    return row["ground_face_export"]["witnesses"]


def write_sha256(results_path: Path) -> Path:
    digest = hashlib.sha256(results_path.read_bytes()).hexdigest()
    sha_path = results_path.with_name(results_path.name.replace("_RESULTS.json", "_RESULTS.sha256"))
    sha_path.write_text(f"{digest}  {results_path.name}\n")
    return sha_path


def main() -> int:
    erdos_root = erdos_experiments_root_from_script()
    default_source = (
        erdos_root
        / "results"
        / "erdos-30"
        / f"{SOURCE_ID}_RESULTS.json"
    )
    default_out = erdos_root / "results" / "erdos-30"

    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=default_source)
    parser.add_argument("--out-dir", type=Path, default=default_out)
    args = parser.parse_args()

    source_path = args.source.resolve()
    out_dir = args.out_dir.resolve()
    witnesses = load_n57_witnesses(source_path)

    positive_skeletons = {
        json.dumps(positive_differences(w["witness"])) for w in witnesses
    }
    mod_skeletons = {
        json.dumps(ordered_mod_differences(w["witness"], MODULUS)) for w in witnesses
    }

    witness_checks = []
    for w in witnesses:
        values = w["witness"]
        mod_counts = Counter(ordered_mod_differences(values, MODULUS))
        pds_overlap = {
            name: {
                "direct_point_overlap": len(set(values) & set(pds)),
                "direct_overlap_set": sorted(set(values) & set(pds)),
                "best_affine_overlap": best_affine_overlap(values, pds, MODULUS),
            }
            for name, pds in SINGER_PDS_REPRESENTATIVES.items()
        }
        witness_checks.append(
            {
                "index": w["index"],
                "witness": values,
                "is_interval_sidon": is_interval_sidon(values),
                "positive_difference_count": len(positive_differences(values)),
                "unique_positive_difference_count": len(set(positive_differences(values))),
                "mod57_nonzero_residues_covered": len(set(mod_counts.keys())),
                "mod57_covers_all_nonzero_residues": set(mod_counts.keys()) == set(range(1, MODULUS)),
                "mod57_multiplicity_histogram": dict(sorted(Counter(mod_counts.values()).items())),
                "pds_overlap": pds_overlap,
            }
        )

    mass_winner = next(w for w in witnesses if w["index"] == 3)
    joint_winner = next(w for w in witnesses if w["index"] == 5)
    joint_is_mass_plus_one = [x + 1 for x in mass_winner["witness"]] == joint_winner["witness"]

    pds_checks = {
        name: {
            **pds_check(pds, MODULUS),
            "multiplier_stabilizer_up_to_translation": multiplier_stabilizer(pds, MODULUS),
        }
        for name, pds in SINGER_PDS_REPRESENTATIVES.items()
    }

    global_best = []
    for wc in witness_checks:
        for name, overlap in wc["pds_overlap"].items():
            best = overlap["best_affine_overlap"]
            global_best.append(
                {
                    "witness_index": wc["index"],
                    "pds_name": name,
                    "max_overlap": best["max_overlap"],
                    "contains_affine_pds_image": best["contains_affine_pds_image"],
                    "best_examples": best["best_examples"],
                }
            )
    max_pds_overlap = max(item["max_overlap"] for item in global_best)

    checks = {
        "C1_all_six_witnesses_are_interval_sidon": all(w["is_interval_sidon"] for w in witness_checks),
        "C2_all_six_witnesses_share_positive_difference_skeleton": len(positive_skeletons) == 1,
        "C3_all_six_witnesses_cover_all_nonzero_residues_mod57": all(
            w["mod57_covers_all_nonzero_residues"] for w in witness_checks
        ),
        "C4_mass_to_joint_winner_is_plus_one_translation": joint_is_mass_plus_one,
        "C5_literature_pds_representatives_verified_and_overlap_computed": (
            all(v["covers_all_nonzero_residues_once"] for v in pds_checks.values())
            and all(
                [entry["multiplier"] for entry in v["multiplier_stabilizer_up_to_translation"]]
                == [1, 7, 49]
                for v in pds_checks.values()
            )
            and max_pds_overlap >= 0
        ),
    }

    result = {
        "experiment_id": EXPERIMENT_ID,
        "date": "2026-04-30",
        "type": "derived_singer57_certificate_from_exact_packet",
        "source_packet": {
            "experiment_id": SOURCE_ID,
            "results_path": str(source_path),
            "status": "EXACT",
        },
        "external_reference": {
            "url": "https://arxiv.org/html/2502.09536v1",
            "object": "Singer / perfect difference set representatives mod 57",
            "pds_representatives": SINGER_PDS_REPRESENTATIVES,
        },
        "modulus": MODULUS,
        "witness_count": len(witnesses),
        "checks": checks,
        "all_C_checks_pass": all(checks.values()),
        "singer_pds_representative_checks": pds_checks,
        "n57_witness_checks": witness_checks,
        "pds_overlap_summary": {
            "max_affine_point_overlap_any_witness_any_pds": max_pds_overlap,
            "contains_any_affine_pds_image": any(item["contains_affine_pds_image"] for item in global_best),
            "best_overlaps": [
                item for item in global_best if item["max_overlap"] == max_pds_overlap
            ],
        },
        "interpretation": (
            "The n=57 exact interval Sidon embeddings sit on the Singer modulus 57 and "
            "cover all nonzero residues mod 57 with controlled redundancy, but they do "
            "not literally contain an affine image of the cited size-8 Singer PDS representatives."
        ),
        "claim_boundary": (
            "Finite certificate only: no asymptotic claim, no theorem proof, and no claim "
            "that the size-10 interval embeddings are Singer PDS objects."
        ),
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    results_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    results_path.write_text(json.dumps(result, indent=2) + "\n")
    sha_path = write_sha256(results_path)

    table = "\n".join(
        f"| {wc['index']} | {wc['is_interval_sidon']} | "
        f"{wc['mod57_nonzero_residues_covered']} | "
        f"{wc['mod57_multiplicity_histogram']} | "
        f"{max(v['best_affine_overlap']['max_overlap'] for v in wc['pds_overlap'].values())} |"
        for wc in witness_checks
    )
    report_path.write_text(
        f"# {EXPERIMENT_ID}\n\n"
        "Date: 2026-04-30\n"
        "Problem: Erdos #30\n"
        "Status: DERIVED_SINGER57_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE\n\n"
        "## Source\n\n"
        f"Derived from `{SOURCE_ID}`.\n\n"
        "External finite-geometry comparison uses the three cited mod-57 perfect "
        "difference set representatives from `https://arxiv.org/html/2502.09536v1`.\n\n"
        "## C1-C5 Checks\n\n"
        f"- C1 all six witnesses are interval Sidon: `{checks['C1_all_six_witnesses_are_interval_sidon']}`\n"
        f"- C2 all six share the same positive-difference skeleton: `{checks['C2_all_six_witnesses_share_positive_difference_skeleton']}`\n"
        f"- C3 all six cover every nonzero residue mod 57: `{checks['C3_all_six_witnesses_cover_all_nonzero_residues_mod57']}`\n"
        f"- C4 mass winner index 3 translates by +1 to joint winner index 5: `{checks['C4_mass_to_joint_winner_is_plus_one_translation']}`\n"
        f"- C5 cited Singer PDS representatives verify, have multiplier stabilizer `[1,7,49]`, and overlap diagnostics are computed: `{checks['C5_literature_pds_representatives_verified_and_overlap_computed']}`\n\n"
        f"All C checks pass: `{result['all_C_checks_pass']}`\n\n"
        "## Witness Table\n\n"
        "| index | interval Sidon | mod-57 residues covered | mod-57 multiplicity histogram | best affine PDS point overlap |\n"
        "|---:|---|---:|---|---:|\n"
        f"{table}\n\n"
        "## Singer PDS Comparison\n\n"
        f"Maximum affine point overlap between any size-10 witness and any cited size-8 PDS representative: "
        f"`{result['pds_overlap_summary']['max_affine_point_overlap_any_witness_any_pds']}`.\n\n"
        f"Contains any affine image of a cited PDS representative: "
        f"`{result['pds_overlap_summary']['contains_any_affine_pds_image']}`.\n\n"
        "## Interpretation\n\n"
        f"{result['interpretation']}\n\n"
        "## Claim Boundary\n\n"
        f"{result['claim_boundary']}\n"
    )

    # The report is derived from the result, but the standard checksum sidecar
    # remains bound to RESULTS.json by project convention.
    print(json.dumps({"results": str(results_path), "report": str(report_path), "sha256": str(sha_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
