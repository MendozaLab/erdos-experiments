#!/usr/bin/env python3
"""Phi/self-replication litmus for the Erdos #30 ground-face branch window."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable


EXPERIMENT_ID = "EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-58-64-2026-05-01"
SOURCE_58_ID = "EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30"
SOURCE_59_64_ID = "EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01"
PHI = (1.0 + math.sqrt(5.0)) / 2.0


def erdos_experiments_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_sha256(results_path: Path) -> Path:
    digest = sha256_file(results_path)
    sha_path = results_path.with_name(results_path.name.replace("_RESULTS.json", "_RESULTS.sha256"))
    sha_path.write_text(f"{digest}  {results_path.name}\n")
    return sha_path


def witness_key(values: Iterable[int]) -> tuple[int, ...]:
    return tuple(values)


def translate_key(values: Iterable[int], shift: int) -> tuple[int, ...]:
    return tuple(x + shift for x in values)


def positive_differences(values: Iterable[int]) -> tuple[int, ...]:
    xs = list(values)
    return tuple(sorted(xs[i] - xs[j] for i in range(len(xs)) for j in range(i)))


def load_rows(source_anchor: Path, source_windows: list[Path], anchor_n: int) -> list[dict]:
    anchor_rows = json.loads(source_anchor.read_text())["rows"]
    rows = [row for row in anchor_rows if row["n"] == anchor_n]
    for source in source_windows:
        rows.extend(json.loads(source.read_text())["rows"])
    rows.sort(key=lambda row: row["n"])
    return rows


def face_keys(row: dict) -> set[tuple[int, ...]]:
    return {witness_key(w["witness"]) for w in row["ground_face_export"]["witnesses"]}


def skeleton_keys(row: dict) -> set[tuple[int, ...]]:
    return {positive_differences(w["witness"]) for w in row["ground_face_export"]["witnesses"]}


def ratio(value: int, previous: int | None) -> float | None:
    if previous in (None, 0):
        return None
    return value / previous


def recurrence_rows(series: list[dict], field: str) -> list[dict]:
    out = []
    for idx in range(2, len(series)):
        current = series[idx][field]
        prev = series[idx - 1][field]
        prev2 = series[idx - 2][field]
        fib_pred = prev + prev2
        phi_pred = PHI * prev
        out.append(
            {
                "n": series[idx]["n"],
                "field": field,
                "value": current,
                "fibonacci_prediction": fib_pred,
                "fibonacci_residual": current - fib_pred,
                "fibonacci_relative_residual": (current - fib_pred) / current if current else None,
                "phi_prediction_from_previous": phi_pred,
                "phi_residual": current - phi_pred,
                "phi_relative_residual": (current - phi_pred) / current if current else None,
            }
        )
    return out


def build_report(result: dict) -> str:
    rows = result["rows"]
    anchor_n = result["anchor_n"]
    table = "\n".join(
        "| {n} | {face_count} | {inherited_union_count} | {inherited_overlap_count} | "
        "{new_face_count} | {skeleton_count} | {new_skeleton_count} | {face_growth_ratio} | {skeleton_growth_ratio} |".format(
            **{
                **row,
                "face_growth_ratio": "NA" if row["face_growth_ratio"] is None else f"{row['face_growth_ratio']:.6f}",
                "skeleton_growth_ratio": "NA"
                if row["skeleton_growth_ratio"] is None
                else f"{row['skeleton_growth_ratio']:.6f}",
            }
        )
        for row in rows
    )
    checks = result["checks"]
    verdict = result["verdict"]
    return (
        f"# {result['experiment_id']}\n\n"
        "Date: 2026-05-01\n"
        "Problem: Erdos #30\n"
        "Status: DERIVED_PHI_REPLICATION_LITMUS / EXACT_PACKET_BACKED / INTERPRETIVE\n\n"
        "## Question\n\n"
        f"Does the post-{anchor_n} branch fan show phi/Fibonacci-mediated self-replication, "
        "or only a more general inherited-plus-translated replication mechanism?\n\n"
        "## Result\n\n"
        f"{verdict}\n\n"
        "## Checks\n\n"
        f"- every post-{anchor_n} face contains previous face union previous face shifted by `+1`: "
        f"`{checks['post58_contains_inherited_union']}`\n"
        f"- finite ratios are within tolerance of phi: `{checks['face_growth_ratios_phi_like']}`\n"
        f"- Fibonacci residuals are small relative to exact face counts: "
        f"`{checks['face_fibonacci_residuals_small']}`\n"
        f"- skeleton ratios are within tolerance of phi: `{checks['skeleton_growth_ratios_phi_like']}`\n\n"
        "## Replication Table\n\n"
        "| n | face | inherited union | inherited overlap | new face | skeletons | new skeletons | face ratio | skeleton ratio |\n"
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n"
        f"{table}\n\n"
        "## Claim Boundary\n\n"
        "This packet can support self-replication language in the finite exact-face sense. "
        "It does not support phi/Fibonacci mediation unless the explicit ratio and recurrence checks pass. "
        "No asymptotic claim and no Erdős #30 proof claim are made.\n"
    )


def main() -> int:
    root = erdos_experiments_root_from_script()
    default_results = root / "results" / "erdos-30"
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-id", default=EXPERIMENT_ID)
    parser.add_argument("--anchor-n", type=int, default=58)
    parser.add_argument("--source-58", type=Path, default=default_results / f"{SOURCE_58_ID}_RESULTS.json")
    parser.add_argument(
        "--source-window",
        type=Path,
        action="append",
        default=None,
        help="Additional exact ground-face packet to append after n=58. May be repeated.",
    )
    parser.add_argument("--out-dir", type=Path, default=default_results)
    parser.add_argument("--phi-tolerance", type=float, default=0.15)
    parser.add_argument("--fib-relative-tolerance", type=float, default=0.15)
    args = parser.parse_args()

    source_windows = args.source_window or [default_results / f"{SOURCE_59_64_ID}_RESULTS.json"]
    raw_rows = load_rows(args.source_58, source_windows, args.anchor_n)
    rows = []
    prev_face: set[tuple[int, ...]] | None = None
    prev_skeletons: set[tuple[int, ...]] | None = None

    for raw in raw_rows:
        current_face = face_keys(raw)
        current_skeletons = skeleton_keys(raw)
        previous_exact = prev_face or set()
        previous_plus_one = {translate_key(face, 1) for face in previous_exact}
        inherited_union = previous_exact | previous_plus_one
        inherited_overlap = previous_exact & previous_plus_one
        inherited_present = current_face & inherited_union
        new_face = current_face - inherited_union

        previous_skeletons = prev_skeletons or set()
        new_skeletons = current_skeletons - previous_skeletons

        row = {
            "n": raw["n"],
            "face_count": len(current_face),
            "skeleton_count": len(current_skeletons),
            "previous_face_count": len(previous_exact),
            "previous_plus_one_count": len(previous_plus_one),
            "inherited_union_count": len(inherited_union),
            "inherited_present_count": len(inherited_present),
            "inherited_overlap_count": len(inherited_overlap),
            "new_face_count": len(new_face),
            "previous_skeleton_count": len(previous_skeletons),
            "new_skeleton_count": len(new_skeletons),
            "contains_inherited_union": inherited_union <= current_face if previous_exact else None,
        }
        if rows:
            row["face_growth_ratio"] = ratio(row["face_count"], rows[-1]["face_count"])
            row["skeleton_growth_ratio"] = ratio(row["skeleton_count"], rows[-1]["skeleton_count"])
            row["new_face_growth_ratio"] = ratio(row["new_face_count"], rows[-1]["new_face_count"])
            row["new_skeleton_growth_ratio"] = ratio(row["new_skeleton_count"], rows[-1]["new_skeleton_count"])
        else:
            row["face_growth_ratio"] = None
            row["skeleton_growth_ratio"] = None
            row["new_face_growth_ratio"] = None
            row["new_skeleton_growth_ratio"] = None
        rows.append(row)

        prev_face = current_face
        prev_skeletons = current_skeletons

    face_recurrence = recurrence_rows(rows, "face_count")
    skeleton_recurrence = recurrence_rows(rows, "skeleton_count")
    new_face_recurrence = recurrence_rows(rows, "new_face_count")
    new_skeleton_recurrence = recurrence_rows(rows, "new_skeleton_count")

    post58_rows = [row for row in rows if row["n"] > args.anchor_n]
    phi_like_face = all(
        row["face_growth_ratio"] is not None and abs(row["face_growth_ratio"] - PHI) <= args.phi_tolerance
        for row in post58_rows
    )
    phi_like_skeleton = all(
        row["skeleton_growth_ratio"] is not None
        and abs(row["skeleton_growth_ratio"] - PHI) <= args.phi_tolerance
        for row in post58_rows
    )
    fib_like_face = all(
        rec["fibonacci_relative_residual"] is not None
        and abs(rec["fibonacci_relative_residual"]) <= args.fib_relative_tolerance
        for rec in face_recurrence
    )

    checks = {
        "post58_contains_inherited_union": all(row["contains_inherited_union"] for row in post58_rows),
        "face_growth_ratios_phi_like": phi_like_face,
        "face_fibonacci_residuals_small": fib_like_face,
        "skeleton_growth_ratios_phi_like": phi_like_skeleton,
    }

    if checks["post58_contains_inherited_union"] and not (
        checks["face_growth_ratios_phi_like"]
        or checks["face_fibonacci_residuals_small"]
        or checks["skeleton_growth_ratios_phi_like"]
    ):
        verdict = (
            "The finite window supports inherited-plus-translated self-replication, "
            "but it does not support phi/Fibonacci mediation at the current tolerance."
        )
    elif checks["post58_contains_inherited_union"]:
        verdict = (
            "The finite window supports inherited-plus-translated self-replication and has "
            "some phi/Fibonacci-compatible diagnostics, but this remains finite evidence only."
        )
    else:
        verdict = (
            "The finite window does not support the inherited-plus-translated replication rule."
        )

    result = {
        "experiment_id": args.experiment_id,
        "date": "2026-05-01",
        "type": "derived_phi_replication_litmus",
        "anchor_n": args.anchor_n,
        "phi": PHI,
        "tolerances": {
            "phi_growth_abs_tolerance": args.phi_tolerance,
            "fibonacci_relative_residual_tolerance": args.fib_relative_tolerance,
        },
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
        "verdict": verdict,
        "rows": rows,
        "recurrence_tests": {
            "face_count": face_recurrence,
            "skeleton_count": skeleton_recurrence,
            "new_face_count": new_face_recurrence,
            "new_skeleton_count": new_skeleton_recurrence,
        },
        "claim_boundary": (
            "Finite exact-face replication litmus only. Supports self-replication only if the "
            "inherited union check passes. Supports phi language only if explicit ratio or recurrence "
            "checks pass. No asymptotic theorem or Sidon proof claim."
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
