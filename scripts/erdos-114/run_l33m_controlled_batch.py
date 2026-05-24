#!/usr/bin/env python3
"""EHP114 L33M controlled remaining-cell batch runner.

This is an orchestration layer only. The proof-facing certificates are still
emitted by the Rust binaries in this crate. The runner enforces the L33M rule:
process missing n=14 cells in lexicographic order and stop on the first blocker
instead of filling gaps with prose.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
VALIDATED = REPO_ROOT / "Erdos114" / "validated_length"
TARGET_EXPERIMENT_ID = "EXP-MATH-EHP114-N14-CONTROLLED-REMAINING-CELL-BATCH-20260506-01"
TARGET_OUTDIR = VALIDATED / TARGET_EXPERIMENT_ID
CAP = 20.672796062619668
CLAIM_CEILING = (
    "local n=14 cell-batch diagnostic only; not a proof of Erdos #114, "
    "not a global n=14 certificate, and not an exact lemniscate-length certificate"
)
PASS_PACKET_STATUS = "LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF"
PASS_BATCH_STATUS = "CONTROLLED_REMAINING_CELL_BATCH_PASS_NOT_GLOBAL_PROOF"
FAIL_BATCH_STATUS = "CONTROLLED_REMAINING_CELL_BATCH_BLOCKED"


@dataclass
class Step:
    name: str
    artifact_id: str
    status: str | None = None
    result_path: str | None = None


@dataclass
class CellRun:
    sub_i: int
    sub_j: int
    cell_tag: str
    status: str = "NOT_STARTED"
    steps: list[Step] = field(default_factory=list)
    slab_total: float | None = None
    source_total: float | None = None
    final_total: float | None = None
    margin_to_cap: float | None = None
    repair_targets: list[str] = field(default_factory=list)
    first_failed_condition: str = "none"
    packet_path: str | None = None


def cell_tag(sub_i: int, sub_j: int) -> str:
    return f"CELL-{sub_i:02d}-{sub_j:02d}"


def exp_id(kind: str, tag: str) -> str:
    return f"EXP-MATH-EHP114-N14-{kind}-{tag}-20260506-01"


def result_path_for(artifact_id: str) -> Path:
    return VALIDATED / artifact_id / f"{artifact_id}_RESULTS.json"


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def dump_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bin_path(name: str) -> Path:
    return SCRIPT_DIR / "target" / "release" / name


def run_cmd(cmd: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        log.flush()
        proc = subprocess.run(
            cmd,
            cwd=SCRIPT_DIR,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"command failed with exit {proc.returncode}: {' '.join(cmd)}")


def ensure_artifact(
    artifact_id: str,
    cmd: list[str],
    log_path: Path,
) -> dict[str, Any]:
    path = result_path_for(artifact_id)
    if not path.exists():
        run_cmd(cmd, log_path)
    if not path.exists():
        raise RuntimeError(f"artifact did not produce results JSON: {path}")
    return load_json(path)


def add_step(cell: CellRun, name: str, artifact_id: str, result: dict[str, Any]) -> None:
    cell.steps.append(
        Step(
            name=name,
            artifact_id=artifact_id,
            status=str(result.get("status", "UNKNOWN")),
            result_path=str(result_path_for(artifact_id)),
        )
    )


def existing_passed_cells() -> set[tuple[int, int]]:
    passed: set[tuple[int, int]] = {(6, 4)}
    for path in VALIDATED.glob("EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-*-20260506-*/"
                               "EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-*-20260506-*_RESULTS.json"):
        try:
            data = load_json(path)
        except Exception:
            continue
        if data.get("status") != PASS_PACKET_STATUS:
            continue
        subcell = data.get("subcell", {})
        try:
            passed.add((int(subcell["sub_i"]), int(subcell["sub_j"])))
        except Exception:
            continue
    return passed


def top_repair_targets(sharpening: dict[str, Any], max_targets: int) -> list[str]:
    stats = sharpening.get("source_stats", {})
    total = float(
        sharpening.get(
            "source_total_validated_length_upper",
            stats.get("total_validated_length_upper", stats.get("sum_branch_length_upper", 0.0)),
        )
    )
    excess = max(0.0, total - CAP)
    rows = sharpening.get("top_branches_by_length") or stats.get("top_branches_by_length", [])
    targets: list[str] = []
    captured = 0.0
    for row in rows:
        key = str(row.get("ownership_key", ""))
        length = float(row.get("length_upper", 0.0))
        if not key:
            continue
        targets.append(key)
        captured += length
        if captured > excess + 0.25 or len(targets) >= max_targets:
            break
    if not targets and rows:
        key = str(rows[0].get("ownership_key", ""))
        if key:
            targets.append(key)
    return targets


def run_cell(sub_i: int, sub_j: int, max_repair_targets: int, log_path: Path) -> CellRun:
    tag = cell_tag(sub_i, sub_j)
    cell = CellRun(sub_i=sub_i, sub_j=sub_j, cell_tag=tag)

    slab_id = exp_id("SLAB-VALIDATED-LENGTH", tag)
    slab = ensure_artifact(
        slab_id,
        [
            str(bin_path("ehp114_n14_validated_length_slab")),
            "--sub-i",
            str(sub_i),
            "--sub-j",
            str(sub_j),
            "--experiment-id",
            slab_id,
            "--outdir",
            str(VALIDATED / slab_id),
            "--z-subdivision",
            "32",
            "--local-bound-subdivision",
            "6",
        ],
        log_path,
    )
    add_step(cell, "slab_source", slab_id, slab)
    cell.slab_total = float(slab.get("total_validated_length_upper", 0.0))

    source_id = slab_id
    source_path = result_path_for(slab_id)
    source_total = cell.slab_total

    if source_total is None or source_total > CAP:
        sharp_id = exp_id(f"{tag}-SLAB-SOURCE-SHARPENING", "")
        sharp_id = sharp_id.replace("--20260506", "-20260506")
        sharpening = ensure_artifact(
            sharp_id,
            [
                str(bin_path("ehp114_n14_cell00_slab_source_sharpening")),
                "--sub-i",
                str(sub_i),
                "--sub-j",
                str(sub_j),
                "--experiment-id",
                sharp_id,
                "--source",
                str(source_path),
                "--outdir",
                str(VALIDATED / sharp_id),
            ],
            log_path,
        )
        add_step(cell, "source_sharpening", sharp_id, sharpening)
        targets = top_repair_targets(sharpening, max_repair_targets)
        cell.repair_targets = targets
        if not targets:
            cell.status = FAIL_BATCH_STATUS
            cell.first_failed_condition = "source_over_cap_but_no_repair_targets"
            return cell

        repair_id = exp_id(f"{tag}-HIGH-SLOPE-SLAB-REPAIR", "")
        repair_id = repair_id.replace("--20260506", "-20260506")
        repair = ensure_artifact(
            repair_id,
            [
                str(bin_path("ehp114_n14_cell00_high_slope_slab_repair")),
                "--sub-i",
                str(sub_i),
                "--sub-j",
                str(sub_j),
                "--experiment-id",
                repair_id,
                "--source",
                str(source_path),
                "--targets",
                ",".join(targets),
                "--outdir",
                str(VALIDATED / repair_id),
            ],
            log_path,
        )
        add_step(cell, "high_slope_repair", repair_id, repair)
        if int(repair.get("target_unresolved_group_count", 1)) != 0:
            cell.status = FAIL_BATCH_STATUS
            cell.first_failed_condition = "high_slope_repair_left_unresolved_groups"
            return cell
        if float(repair.get("adjusted_total_if_replaced", CAP + 1.0)) > CAP:
            cell.status = FAIL_BATCH_STATUS
            cell.first_failed_condition = "high_slope_repair_still_over_cap"
            return cell

        overlay_id = exp_id(f"{tag}-REPAIRED-SLAB-SOURCE-OVERLAY", "")
        overlay_id = overlay_id.replace("--20260506", "-20260506")
        overlay = ensure_artifact(
            overlay_id,
            [
                str(bin_path("ehp114_n14_cell00_repaired_slab_source_overlay")),
                "--sub-i",
                str(sub_i),
                "--sub-j",
                str(sub_j),
                "--experiment-id",
                overlay_id,
                "--source",
                str(source_path),
                "--repair",
                str(result_path_for(repair_id)),
                "--outdir",
                str(VALIDATED / overlay_id),
            ],
            log_path,
        )
        add_step(cell, "repaired_source_overlay", overlay_id, overlay)
        source_id = overlay_id
        source_path = result_path_for(overlay_id)
        source_total = float(overlay.get("total_validated_length_upper", 0.0))

    cell.source_total = source_total

    smoke_id = exp_id("PER-CELL-SOURCE-GENERATION-SMOKE", tag)
    smoke = ensure_artifact(
        smoke_id,
        [
            str(bin_path("ehp114_n14_per_cell_source_generation_smoke")),
            "--sub-i",
            str(sub_i),
            "--sub-j",
            str(sub_j),
            "--experiment-id",
            smoke_id,
            "--source",
            str(source_path),
            "--outdir",
            str(VALIDATED / smoke_id),
        ],
        log_path,
    )
    add_step(cell, "source_smoke", smoke_id, smoke)
    if smoke.get("status") != "PER_CELL_SOURCE_CHAIN_READY":
        cell.status = FAIL_BATCH_STATUS
        cell.first_failed_condition = str(smoke.get("first_failed_condition", "source_smoke_failed"))
        return cell

    branch_id = exp_id("BRANCH-ISOLATION-COLLAR-ATLAS", tag)
    branch = ensure_artifact(
        branch_id,
        [
            str(bin_path("ehp114_n14_branch_isolation_collar_atlas")),
            "--sub-i",
            str(sub_i),
            "--sub-j",
            str(sub_j),
            "--experiment-id",
            branch_id,
            "--source",
            str(source_path),
            "--outdir",
            str(VALIDATED / branch_id),
            "--branch-limit",
            "16",
            "--max-depth",
            "2",
        ],
        log_path,
    )
    add_step(cell, "branch_atlas", branch_id, branch)

    normal_id = exp_id("NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT", tag)
    normal = ensure_artifact(
        normal_id,
        [
            str(bin_path("ehp114_n14_normal_collar_critical_exclusion_pilot")),
            "--sub-i",
            str(sub_i),
            "--sub-j",
            str(sub_j),
            "--experiment-id",
            normal_id,
            "--source",
            str(result_path_for(branch_id)),
            "--outdir",
            str(VALIDATED / normal_id),
        ],
        log_path,
    )
    add_step(cell, "normal_collar", normal_id, normal)

    third_id = exp_id("THIRD-ORDER-COLLAR-REMAINDER-PILOT", tag)
    third = ensure_artifact(
        third_id,
        [
            str(bin_path("ehp114_n14_third_order_collar_remainder_pilot")),
            "--sub-i",
            str(sub_i),
            "--sub-j",
            str(sub_j),
            "--experiment-id",
            third_id,
            "--source",
            str(result_path_for(normal_id)),
            "--outdir",
            str(VALIDATED / third_id),
        ],
        log_path,
    )
    add_step(cell, "third_order", third_id, third)

    global_id = exp_id("GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET", tag)
    global_result = ensure_artifact(
        global_id,
        [
            str(bin_path("ehp114_n14_global_critical_point_exclusion_target")),
            "--sub-i",
            str(sub_i),
            "--sub-j",
            str(sub_j),
            "--experiment-id",
            global_id,
            "--source",
            str(result_path_for(third_id)),
            "--outdir",
            str(VALIDATED / global_id),
        ],
        log_path,
    )
    add_step(cell, "global_critical_target", global_id, global_result)

    l24_id = exp_id("VARIABLE-REGULAR-RESIDUAL-CLOSURE", tag)
    l24 = ensure_artifact(
        l24_id,
        [
            str(bin_path("ehp114_n14_regular_residual_decomposition")),
            "--sub-i",
            str(sub_i),
            "--sub-j",
            str(sub_j),
            "--experiment-id",
            l24_id,
            "--source",
            str(result_path_for(global_id)),
            "--outdir",
            str(VALIDATED / l24_id),
            "--wall-mode",
            "taylor-y-monotone",
            "--max-depth",
            "5",
            "--x-adaptive-depth",
            "8",
        ],
        log_path,
    )
    add_step(cell, "regular_residual_closure", l24_id, l24)
    if l24.get("status") != "REGULAR_RESIDUAL_PASS_NOT_GLOBAL_PROOF":
        cell.status = FAIL_BATCH_STATUS
        cell.first_failed_condition = str(l24.get("first_failed_condition", "regular_residual_failed"))
        return cell

    l27_id = exp_id("CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8", tag)
    l27 = ensure_artifact(
        l27_id,
        [
            str(bin_path("ehp114_n14_critical_candidate_affine_gradient")),
            "--sub-i",
            str(sub_i),
            "--sub-j",
            str(sub_j),
            "--experiment-id",
            l27_id,
            "--source",
            str(result_path_for(global_id)),
            "--outdir",
            str(VALIDATED / l27_id),
            "--test-mode",
            "pprime",
            "--param-subdivision",
            "8",
        ],
        log_path,
    )
    add_step(cell, "direct_pprime_p8", l27_id, l27)

    l28_id = exp_id("PPRIME-ROOT-LOCATION", tag)
    l28 = ensure_artifact(
        l28_id,
        [
            str(bin_path("ehp114_n14_pprime_root_location")),
            "--sub-i",
            str(sub_i),
            "--sub-j",
            str(sub_j),
            "--experiment-id",
            l28_id,
            "--source",
            str(result_path_for(l27_id)),
            "--outdir",
            str(VALIDATED / l28_id),
            "--param-subdivision",
            "8",
        ],
        log_path,
    )
    add_step(cell, "pprime_root_location", l28_id, l28)

    l29_id = exp_id("SHARP-PPRIME-ROOT-LOCATION", tag)
    l29 = ensure_artifact(
        l29_id,
        [
            str(bin_path("ehp114_n14_sharp_pprime_root_location")),
            "--sub-i",
            str(sub_i),
            "--sub-j",
            str(sub_j),
            "--experiment-id",
            l29_id,
            "--source",
            str(result_path_for(l28_id)),
            "--outdir",
            str(VALIDATED / l29_id),
            "--param-subdivision",
            "8",
            "--max-depth",
            "4",
        ],
        log_path,
    )
    add_step(cell, "sharp_pprime_root_location", l29_id, l29)
    if l29.get("status") != "SHARP_PPRIME_ROOT_LOCATION_PASS_NOT_GLOBAL_PROOF":
        cell.status = FAIL_BATCH_STATUS
        cell.first_failed_condition = str(l29.get("first_failed_condition", "sharp_pprime_failed"))
        return cell

    l31_id = exp_id("RESIDUAL-CHAIN-INTEGRATION-CERT", tag)
    l31 = ensure_artifact(
        l31_id,
        [
            str(bin_path("ehp114_n14_residual_chain_integration_cert")),
            "--sub-i",
            str(sub_i),
            "--sub-j",
            str(sub_j),
            "--experiment-id",
            l31_id,
            "--l24",
            str(result_path_for(l24_id)),
            "--l27",
            str(result_path_for(l27_id)),
            "--l28",
            str(result_path_for(l28_id)),
            "--l29",
            str(result_path_for(l29_id)),
            "--outdir",
            str(VALIDATED / l31_id),
        ],
        log_path,
    )
    add_step(cell, "residual_chain_integration", l31_id, l31)
    if l31.get("status") != "RESIDUAL_CHAIN_INTEGRATION_PASS_NOT_GLOBAL_PROOF":
        cell.status = FAIL_BATCH_STATUS
        cell.first_failed_condition = str(l31.get("first_failed_condition", "residual_chain_failed"))
        return cell

    residual_universe_path: Path | None = None
    if int(l31.get("expected_residual_chain_total_count", 64)) != 64:
        residual_universe_id = exp_id("RESIDUAL-UNIVERSE-ACCOUNTING-CERT", tag)
        residual_universe = ensure_artifact(
            residual_universe_id,
            [
                str(bin_path("ehp114_n14_residual_universe_accounting_cert")),
                "--sub-i",
                str(sub_i),
                "--sub-j",
                str(sub_j),
                "--experiment-id",
                residual_universe_id,
                "--branch",
                str(result_path_for(branch_id)),
                "--normal",
                str(result_path_for(normal_id)),
                "--third",
                str(result_path_for(third_id)),
                "--global",
                str(result_path_for(global_id)),
                "--l24",
                str(result_path_for(l24_id)),
                "--l27",
                str(result_path_for(l27_id)),
                "--l28",
                str(result_path_for(l28_id)),
                "--l29",
                str(result_path_for(l29_id)),
                "--l31",
                str(result_path_for(l31_id)),
                "--outdir",
                str(VALIDATED / residual_universe_id),
            ],
            log_path,
        )
        add_step(cell, "residual_universe_accounting", residual_universe_id, residual_universe)
        if residual_universe.get("status") != "RESIDUAL_UNIVERSE_ACCOUNTING_PASS_NOT_GLOBAL_PROOF":
            cell.status = FAIL_BATCH_STATUS
            cell.first_failed_condition = str(
                residual_universe.get("first_failed_condition", "residual_universe_accounting_failed")
            )
            return cell
        residual_universe_path = result_path_for(residual_universe_id)

    packet_id = exp_id("LOCAL-CELL-CERTIFICATE-PACKET", tag)
    packet_cmd = [
        str(bin_path("ehp114_n14_local_hard_cell_certificate_packet")),
        "--sub-i",
        str(sub_i),
        "--sub-j",
        str(sub_j),
        "--experiment-id",
        packet_id,
        "--l31",
        str(result_path_for(l31_id)),
        "--outdir",
        str(VALIDATED / packet_id),
    ]
    if residual_universe_path is not None:
        packet_cmd.extend(["--residual-universe", str(residual_universe_path)])
    packet = ensure_artifact(
        packet_id,
        packet_cmd,
        log_path,
    )
    add_step(cell, "local_cell_packet", packet_id, packet)
    cell.packet_path = str(result_path_for(packet_id))
    cell.final_total = float(packet.get("total_validated_length_upper", 0.0))
    cell.margin_to_cap = float(packet.get("margin_to_cap", 0.0))
    if packet.get("status") != PASS_PACKET_STATUS:
        cell.status = FAIL_BATCH_STATUS
        cell.first_failed_condition = str(packet.get("first_failed_condition", "local_packet_failed"))
        return cell

    cell.status = PASS_PACKET_STATUS
    return cell


def cell_to_json(cell: CellRun) -> dict[str, Any]:
    return {
        "subcell": {"sub_i": cell.sub_i, "sub_j": cell.sub_j},
        "cell_tag": cell.cell_tag,
        "status": cell.status,
        "slab_total": cell.slab_total,
        "source_total": cell.source_total,
        "final_total": cell.final_total,
        "margin_to_cap": cell.margin_to_cap,
        "repair_targets": cell.repair_targets,
        "first_failed_condition": cell.first_failed_condition,
        "packet_path": cell.packet_path,
        "steps": [
            {
                "name": step.name,
                "artifact_id": step.artifact_id,
                "status": step.status,
                "result_path": step.result_path,
            }
            for step in cell.steps
        ],
    }


def write_report(result: dict[str, Any], report_path: Path) -> None:
    rows = []
    for cell in result["processed_cells"]:
        rows.append(
            f"| `{cell['cell_tag']}` | `{cell['status']}` | "
            f"`{cell.get('final_total')}` | `{cell.get('margin_to_cap')}` | "
            f"`{','.join(cell.get('repair_targets') or []) or 'none'}` |"
        )
    table = "\n".join(rows) if rows else "| none | none | none | none | none |"
    report = f"""# EHP114 n=14 Controlled Remaining-Cell Batch

Experiment: `{result['experiment_id']}`

## Verdict

- Status: `{result['status']}`
- Processed cell count: `{result['processed_cell_count']}`
- Passing cell count: `{result['passing_cell_count']}`
- Stopped on blocker: `{result['stopped_on_blocker']}`
- First failed condition: `{result['first_failed_condition']}`
- Remaining cell count after this run: `{result['remaining_cell_count_after_run']}`
- Exact length cap: `{result['exact_length_cap']}`

| Cell | Status | Final length upper | Margin to cap | Repair targets |
|---|---|---:|---:|---|
{table}

## Interpretation

This batch is a proof-track control gate. It runs the same per-cell source,
repair, residual-chain, and local-packet machinery used by the representative
cells, and it stops on the first blocker. A pass here means only that the
processed local cells received theorem-packet-shaped artifacts under the local
n=14 contract.

## Claim Ceiling

{result['claim_ceiling']}
"""
    report_path.write_text(report, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=TARGET_OUTDIR)
    parser.add_argument("--experiment-id", default=TARGET_EXPERIMENT_ID)
    parser.add_argument("--max-cells", type=int, default=0, help="0 means process all missing cells")
    parser.add_argument("--max-repair-targets", type=int, default=10)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    result_path = args.outdir / f"{args.experiment_id}_RESULTS.json"
    report_path = args.outdir / f"{args.experiment_id}_REPORT.md"
    sha_path = args.outdir / f"{args.experiment_id}_RESULTS.sha256"
    log_path = args.outdir / f"{args.experiment_id}_RUN.log"

    if result_path.exists() or report_path.exists() or sha_path.exists():
        raise SystemExit(f"refusing to overwrite existing L33M artifact in {args.outdir}")

    passed_before = existing_passed_cells()
    all_cells = [(i, j) for i in range(8) for j in range(8)]
    missing = [cell for cell in all_cells if cell not in passed_before]
    if args.max_cells > 0:
        missing = missing[: args.max_cells]

    processed: list[CellRun] = []
    started = time.time()
    blocker = "none"
    try:
        for sub_i, sub_j in missing:
            cell = run_cell(sub_i, sub_j, args.max_repair_targets, log_path)
            processed.append(cell)
            if cell.status != PASS_PACKET_STATUS:
                blocker = f"{cell.cell_tag}: {cell.first_failed_condition}"
                break
    except Exception as exc:
        if processed:
            processed[-1].status = FAIL_BATCH_STATUS
            processed[-1].first_failed_condition = str(exc)
            blocker = f"{processed[-1].cell_tag}: {exc}"
        else:
            blocker = str(exc)

    passing = [cell for cell in processed if cell.status == PASS_PACKET_STATUS]
    status = PASS_BATCH_STATUS if len(passing) == len(processed) and blocker == "none" else FAIL_BATCH_STATUS
    passed_after = passed_before | {(cell.sub_i, cell.sub_j) for cell in passing}
    remaining_after = [cell for cell in all_cells if cell not in passed_after]
    result = {
        "experiment_id": args.experiment_id,
        "status": status,
        "started_unix": int(started),
        "completed_unix": int(time.time()),
        "duration_seconds": time.time() - started,
        "existing_passed_cell_count_before_run": len(passed_before),
        "targeted_missing_cell_count": len(missing),
        "processed_cell_count": len(processed),
        "passing_cell_count": len(passing),
        "stopped_on_blocker": blocker != "none",
        "first_failed_condition": blocker,
        "remaining_cell_count_after_run": len(remaining_after),
        "remaining_cells_after_run": [{"sub_i": i, "sub_j": j, "cell_tag": cell_tag(i, j)} for i, j in remaining_after],
        "processed_cells": [cell_to_json(cell) for cell in processed],
        "exact_length_cap": CAP,
        "claim_ceiling": CLAIM_CEILING,
    }
    dump_json(result_path, result)
    write_report(result, report_path)
    digest = sha256_file(result_path)
    sha_path.write_text(f"{digest}  {result_path.name}\n", encoding="utf-8")
    print(json.dumps({
        "status": status,
        "processed_cell_count": len(processed),
        "passing_cell_count": len(passing),
        "first_failed_condition": blocker,
        "remaining_cell_count_after_run": len(remaining_after),
        "result_path": str(result_path),
    }, indent=2))
    return 0 if status == PASS_BATCH_STATUS else 2


if __name__ == "__main__":
    raise SystemExit(main())
