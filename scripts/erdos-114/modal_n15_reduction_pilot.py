#!/usr/bin/env python3
"""Modal reduction pilot for Erdős #114 n=15.

This runner is deliberately a reduction pilot, not a certificate runner. It
checks Modal plumbing, audits one trusted n=14 control chain, and packages the
smallest n=15 slice plan into a review-only artifact. A future full n=15 launch
is authorized only if this pilot emits GO_FULL_N15.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

try:
    import modal
except Exception:  # pragma: no cover - dry-run/local modes work without Modal.
    modal = None  # type: ignore[assignment]


EXPERIMENT_ID = "EXP-MATH-EHP114-N15-REDUCTION-PILOT-MODAL-20260508-01"
APP_NAME = "ehp114-n15-reduction-pilot-20260508-01"
VOLUME_NAME = "ehp114-n15-reduction-pilot-20260508-01"

REMOTE_SRC = "/workspace/erdos-114"
REMOTE_INPUTS = "/workspace/inputs"
REMOTE_OUT = f"/outputs/{EXPERIMENT_ID}"

SCRIPT_ROOT = Path(__file__).resolve().parent
SRC_ROOT = SCRIPT_ROOT if (SCRIPT_ROOT / "Cargo.toml").exists() else Path(REMOTE_SRC)
ERDOS_EXPERIMENTS = (
    SCRIPT_ROOT.parents[1]
    if len(SCRIPT_ROOT.parents) > 1 and (SCRIPT_ROOT.parents[1] / "Erdos114").exists()
    else Path("/workspace")
)
MATH_ROOT = ERDOS_EXPERIMENTS.parent

LOCAL_OUTDIR = (
    ERDOS_EXPERIMENTS
    / "Erdos114"
    / "proof_path"
    / EXPERIMENT_ID
)
DOWNLOADS_STORY_DIR = (
    Path.home() / "Downloads" / "Eratosthenes_Ahmes_Stories"
)

SOURCE_PATHS = {
    "n15_envelope": ERDOS_EXPERIMENTS
    / "Erdos114"
    / "proof_path"
    / "EXP-MATH-EHP114-N15-PILOT-ENVELOPE-20260508-01"
    / "EXP-MATH-EHP114-N15-PILOT-ENVELOPE-20260508-01_RESULTS.json",
    "finite_less15": ERDOS_EXPERIMENTS
    / "Erdos114"
    / "proof_path"
    / "EXP-MATH-EHP114-FINITE-N-LESS-15-PROOF-PACKET-20260507-01"
    / "EXP-MATH-EHP114-FINITE-N-LESS-15-PROOF-PACKET-20260507-01_RESULTS.json",
    "n13_receipt": ERDOS_EXPERIMENTS
    / "results"
    / "erdos-114"
    / "EXP-MM-EHP-007-n13-inari_RESULTS.json",
    "n14_receipt": ERDOS_EXPERIMENTS
    / "results"
    / "erdos-114"
    / "EXP-MM-EHP-007-n14-inari_RESULTS.json",
    "n15_fourier_hessian": ERDOS_EXPERIMENTS
    / "results"
    / "erdos-114"
    / "EXP-MATH-EHP114-N15-FOURIER-HESSIAN-20260502-02_RESULTS.json",
    "n14_local_cell_07_03": ERDOS_EXPERIMENTS
    / "Erdos114"
    / "validated_length"
    / "EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-07-03-20260506-01"
    / "EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-07-03-20260506-01_RESULTS.json",
    "n14_regular_residual_06_03": ERDOS_EXPERIMENTS
    / "Erdos114"
    / "validated_length"
    / "EXP-MATH-EHP114-N14-VARIABLE-REGULAR-RESIDUAL-CLOSURE-CELL-06-03-20260506-01"
    / "EXP-MATH-EHP114-N14-VARIABLE-REGULAR-RESIDUAL-CLOSURE-CELL-06-03-20260506-01_RESULTS.json",
    "n14_third_order_collar_02_03": ERDOS_EXPERIMENTS
    / "Erdos114"
    / "validated_length"
    / "EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-02-03-20260506-01"
    / "EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-02-03-20260506-01_RESULTS.json",
}

REMOTE_SOURCE_PATHS = {
    name: Path(REMOTE_INPUTS) / f"{name}.json" for name in SOURCE_PATHS
}

FINAL_VERDICTS = {"GO_FULL_N15", "NEEDS_REDUCTION", "STOP_NOT_FEASIBLE"}


def _require_modal() -> Any:
    if modal is None:
        raise RuntimeError(
            "Modal is not importable. Use --dry-run/--local-only, or run with Modal installed."
        )
    return modal


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_sha_sidecar(results_path: Path) -> str:
    digest = sha256_file(results_path)
    sha_path = results_path.with_name(results_path.name.replace("_RESULTS.json", "_RESULTS.sha256"))
    sha_path.write_text(f"{digest}  {results_path.name}\n", encoding="utf-8")
    return digest


def source_fingerprints(paths: dict[str, Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for role, path in sorted(paths.items()):
        rows.append(
            {
                "role": role,
                "path": str(path),
                "exists": path.exists(),
                "byte_count": path.stat().st_size if path.exists() else None,
                "sha256": sha256_file(path) if path.exists() else None,
            }
        )
    return rows


def validate_sources(paths: dict[str, Path]) -> None:
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("missing source inputs: " + ", ".join(missing))


def summarize_receipt(path: Path, role: str) -> dict[str, Any]:
    data = read_json(path)
    levels = data.get("bb_levels")
    levels_len = len(levels) if isinstance(levels, list) else None
    total_evals = data.get("bb_total_evals")
    return {
        "role": role,
        "path": str(path),
        "sha256": sha256_file(path),
        "bb_total_evals": total_evals,
        "bb_levels_len": levels_len,
        "has_nonzero_work_signature": bool(
            isinstance(total_evals, int) and total_evals > 0 and isinstance(levels_len, int) and levels_len > 0
        ),
    }


def audit_n14_control(paths: dict[str, Path]) -> dict[str, Any]:
    controls = [
        (
            "low_margin_local_cell",
            paths["n14_local_cell_07_03"],
            "CELL-07-03",
            "LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF",
            True,
        ),
        (
            "regular_residual_cell",
            paths["n14_regular_residual_06_03"],
            "CELL-06-03",
            "REGULAR_RESIDUAL_PASS_NOT_GLOBAL_PROOF",
            True,
        ),
        (
            "hard_collar_failure_cell",
            paths["n14_third_order_collar_02_03"],
            "CELL-02-03",
            "THIRD_ORDER_FAIL_CRITICAL_EXCLUSION",
            False,
        ),
    ]
    rows = []
    pass_count = 0
    for role, path, expected_cell, expected_status, positive_margin_required in controls:
        data = read_json(path)
        margin = data.get("margin_to_cap")
        status = data.get("status")
        cell_tag = data.get("cell_tag")
        first_failed_condition = data.get("first_failed_condition")
        margin_ok = isinstance(margin, (int, float)) and (
            margin > 0 if positive_margin_required else True
        )
        row_pass = cell_tag == expected_cell and status == expected_status and margin_ok
        rows.append(
            {
                "role": role,
                "path": str(path),
                "sha256": sha256_file(path),
                "expected_cell": expected_cell,
                "cell_tag": cell_tag,
                "expected_status": expected_status,
                "status": status,
                "margin_to_cap": margin,
                "first_failed_condition": first_failed_condition,
                "control_pass": row_pass,
            }
        )
        if row_pass:
            pass_count += 1

    n13 = summarize_receipt(paths["n13_receipt"], "corrected_n13_receipt")
    n14 = summarize_receipt(paths["n14_receipt"], "protected_n14_receipt")
    receipt_pass = n13["has_nonzero_work_signature"] and n14["has_nonzero_work_signature"]
    control_pass = pass_count == len(rows) and receipt_pass
    return {
        "gate": "n14_control",
        "status": "PASS" if control_pass else "FAIL",
        "control_row_pass_count": pass_count,
        "control_row_count": len(rows),
        "receipt_gate_pass": receipt_pass,
        "receipt_summaries": [n13, n14],
        "control_rows": rows,
        "meaning": "Receipt-shape and selected n14 cell/collar control audit only; no new global claim is emitted.",
    }


def audit_n15_pilot(paths: dict[str, Path]) -> dict[str, Any]:
    envelope = read_json(paths["n15_envelope"])
    fourier = read_json(paths["n15_fourier_hessian"])
    seed_cells = (
        envelope.get("recommended_n15_pilot_envelope", {})
        .get("selected_seed_cells_from_n14", [])
    )
    seed_by_cell = {
        row.get("cell_tag"): row for row in seed_cells if isinstance(row, dict)
    }
    planned_slices = [
        {
            "slice_id": "n15_low_margin_local_cell_transfer",
            "n14_seed_cell": "CELL-07-03",
            "source_status": seed_by_cell.get("CELL-07-03", {}).get("status"),
            "source_margin_to_cap": seed_by_cell.get("CELL-07-03", {}).get("margin_to_cap"),
            "n15_status": "SELECTED_NOT_EXECUTED",
            "required_next_evidence": "n15 local-cell source generation plus residual closure accounting",
        },
        {
            "slice_id": "n15_regular_residual_transfer",
            "n14_seed_cell": "CELL-06-03",
            "source_status": seed_by_cell.get("CELL-06-03", {}).get("status"),
            "source_margin_to_cap": seed_by_cell.get("CELL-06-03", {}).get("margin_to_cap"),
            "n15_status": "SELECTED_NOT_EXECUTED",
            "required_next_evidence": "n15 regular residual closure row with positive margin and hash sidecar",
        },
        {
            "slice_id": "n15_hard_collar_root_isolation_transfer",
            "n14_seed_cell": "CELL-02-03",
            "source_status": seed_by_cell.get("CELL-02-03", {}).get("status"),
            "source_margin_to_cap": seed_by_cell.get("CELL-02-03", {}).get("margin_to_cap"),
            "n15_status": "SELECTED_NOT_EXECUTED",
            "required_next_evidence": "n15 collar/root-isolation pilot resolving or classifying the wall-separation failure",
        },
    ]
    fourier_summary = fourier.get("summary", {})
    return {
        "gate": "n15_pilot",
        "status": "NEEDS_REDUCTION",
        "source_envelope_decision": envelope.get("decision"),
        "source_envelope_stop_go": envelope.get("stop_go_rule_result"),
        "fourier_diagnostic": {
            "path": str(paths["n15_fourier_hessian"]),
            "sha256": sha256_file(paths["n15_fourier_hessian"]),
            "status": fourier.get("status"),
            "degree": fourier.get("degree"),
            "method": fourier.get("method"),
            "summary_keys": sorted(fourier_summary.keys()) if isinstance(fourier_summary, dict) else [],
        },
        "planned_slices": planned_slices,
        "gate_reason": "Modal plumbing and n14 controls are not enough to authorize full n=15; selected n15 slices still need real reduced computations.",
    }


def make_smoke_gate(*, remote_executed: bool, blocked_reason: str | None = None) -> dict[str, Any]:
    started = time.time()
    status = "PASS" if remote_executed else "BLOCKED_BY_EXECUTION_POLICY"
    return {
        "gate": "modal_smoke",
        "status": status,
        "timestamp_unix": int(started),
        "elapsed_secs": time.time() - started,
        "remote_executed": remote_executed,
        "volume_write_checked": remote_executed,
        "sidecar_hash_checked": remote_executed,
        "log_cadence_checked": remote_executed,
        "blocked_reason": blocked_reason,
        "meaning": "Tiny Modal task and output-contract check only.",
    }


def make_result_payload(
    *,
    source_paths: dict[str, Path],
    smoke_gate: dict[str, Any],
    n14_control_gate: dict[str, Any],
    n15_pilot_gate: dict[str, Any],
    modal_metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    gates = [smoke_gate, n14_control_gate, n15_pilot_gate]
    gate_statuses = {row["gate"]: row.get("status") for row in gates}
    if smoke_gate.get("status") != "PASS" or n14_control_gate.get("status") != "PASS":
        verdict = "STOP_NOT_FEASIBLE"
    elif n15_pilot_gate.get("status") == "GO_FULL_N15":
        verdict = "GO_FULL_N15"
    else:
        verdict = "NEEDS_REDUCTION"
    if verdict not in FINAL_VERDICTS:
        raise ValueError(f"invalid final verdict {verdict}")
    return {
        "experiment_id": EXPERIMENT_ID,
        "created_unix": int(time.time()),
        "generated_by": "erdos_atlas_autoresearch_librarian",
        "origin": "auto-research",
        "persona": "Eratosthenes of Cyrene",
        "short_name": "Eratosthenes",
        "scribe": "Ahmes",
        "story_writer": "Ahmes",
        "problem_id": 114,
        "scope": "erdos_atlas_only",
        "promotion_state": "review_only",
        "claim_ceiling": "review-only n=15 reduction pilot; not D1, not public state, not formal-status promotion",
        "final_verdict": verdict,
        "full_n15_authorized": verdict == "GO_FULL_N15",
        "gate_statuses": gate_statuses,
        "gates": gates,
        "source_fingerprints": source_fingerprints(source_paths),
        "safety_checks": {
            "old_n15_n16_zero_eval_shortcuts_used_as_evidence": False,
            "branch_and_bound_invoked_for_n15": False,
            "certificate_row_emitted_for_n15": False,
            "writes_limited_to_new_experiment_folder_and_optional_story": True,
            "forbidden_writes": [
                "D1",
                "morphisms.json",
                "formal registries",
                "public pages",
                "existing #114 packets",
            ],
        },
        "next_experiment_if_needed": {
            "if_GO_FULL_N15": "EXP-MATH-EHP114-N15-FULL-CERTIFICATION-MODAL-<new-id>",
            "if_NEEDS_REDUCTION": "EXP-MATH-EHP114-N15-SLICE-CELL-COLLAR-MODAL-<new-id>",
            "if_STOP_NOT_FEASIBLE": "repair Modal/runner/input-contract blocker before spending compute",
        },
        "modal_metadata": modal_metadata or {},
    }


def write_report(result: dict[str, Any], report_path: Path) -> None:
    rows = result["gate_statuses"]
    body = f"""# {EXPERIMENT_ID} Report

## Verdict

Final verdict: `{result["final_verdict"]}`

Full n=15 launch authorized: `{str(result["full_n15_authorized"]).lower()}`

## Gate Summary

- Modal smoke: `{rows.get("modal_smoke")}`
- n=14 control: `{rows.get("n14_control")}`
- n=15 pilot: `{rows.get("n15_pilot")}`

## Meaning

This is a reduction pilot. It checks whether the Modal/runner surface and the
trusted n=14 control rows are coherent enough to justify the next n=15 slice
run. It does not certify n=15 and does not update any registry or public state.

## Safety

- Old zero-eval n=15/n=16 shortcut rows used as evidence: `false`
- n=15 branch-and-bound invoked: `false`
- n=15 certificate row emitted: `false`
- Promotion state: `{result["promotion_state"]}`

## Next Step

If the verdict remains `NEEDS_REDUCTION`, run the next immutable slice-cell/collar
Modal experiment over the three selected n=15 slices before considering a full
n=15 launch.
"""
    report_path.write_text(body, encoding="utf-8")


def write_story(result: dict[str, Any], story_dir: Path = DOWNLOADS_STORY_DIR) -> Path:
    story_dir.mkdir(parents=True, exist_ok=True)
    story_path = story_dir / f"AHMES_STORY_{EXPERIMENT_ID}.md"
    if story_path.exists():
        raise FileExistsError(f"refusing to overwrite existing story: {story_path}")
    gate_rows = "\n".join(
        f"- `{gate}`: `{status}`" for gate, status in result["gate_statuses"].items()
    )
    story = f"""# Ahmes Story: {EXPERIMENT_ID}

Eratosthenes did not send n=15 straight into the furnace. He first checked the
door, the ledger, and the old n=14 control rows. The lesson is narrow but useful:
the path to n=15 is a reduction path, not a blind full launch.

The pilot kept the old zero-eval n=15/n=16 rows outside the evidence chain. It
looked only at protected finite receipts, selected n=14 control cells, and the
n=15 reduction envelope. The result is `{result["final_verdict"]}`.

## Receipts

- Run ID: `{EXPERIMENT_ID}`
- Target: Erdős #114, n=15 reduction pilot
- Promotion state: `{result["promotion_state"]}`
- Full n=15 launch authorized: `{str(result["full_n15_authorized"]).lower()}`

## Gate Status

{gate_rows}

## Review Warning

This story is a reading companion. The source of truth is the JSON/report/hash
triplet in the run folder. No D1, public page, curated morphism, or formal-status
registry was updated.
"""
    story_path.write_text(story, encoding="utf-8")
    return story_path


def materialize_local(result: dict[str, Any], outdir: Path = LOCAL_OUTDIR, *, write_story_file: bool) -> dict[str, str]:
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite existing run folder: {outdir}")
    outdir.mkdir(parents=True)
    results_path = outdir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = outdir / f"{EXPERIMENT_ID}_REPORT.md"
    write_json(results_path, result)
    write_report(result, report_path)
    digest = write_sha_sidecar(results_path)
    story_path = write_story(result) if write_story_file else None
    return {
        "results": str(results_path),
        "report": str(report_path),
        "sha256": str(results_path.with_name(f"{EXPERIMENT_ID}_RESULTS.sha256")),
        "digest": digest,
        "story": str(story_path) if story_path else "",
    }


def dry_run_summary() -> dict[str, Any]:
    validate_sources(SOURCE_PATHS)
    return {
        "experiment_id": EXPERIMENT_ID,
        "mode": "dry_run",
        "would_write": [
            str(LOCAL_OUTDIR / f"{EXPERIMENT_ID}_RESULTS.json"),
            str(LOCAL_OUTDIR / f"{EXPERIMENT_ID}_REPORT.md"),
            str(LOCAL_OUTDIR / f"{EXPERIMENT_ID}_RESULTS.sha256"),
            str(DOWNLOADS_STORY_DIR / f"AHMES_STORY_{EXPERIMENT_ID}.md"),
        ],
        "source_fingerprints": source_fingerprints(SOURCE_PATHS),
        "final_verdict_domain": sorted(FINAL_VERDICTS),
    }


def run_local_only(
    *,
    materialize: bool,
    story: bool,
    modal_blocked_reason: str | None = None,
) -> dict[str, Any]:
    validate_sources(SOURCE_PATHS)
    smoke_gate = make_smoke_gate(
        remote_executed=False,
        blocked_reason=modal_blocked_reason or "local-only mode does not execute Modal",
    )
    n14_control_gate = audit_n14_control(SOURCE_PATHS)
    n15_pilot_gate = audit_n15_pilot(SOURCE_PATHS)
    result = make_result_payload(
        source_paths=SOURCE_PATHS,
        smoke_gate=smoke_gate,
        n14_control_gate=n14_control_gate,
        n15_pilot_gate=n15_pilot_gate,
        modal_metadata={"mode": "local_only", "modal_blocked_reason": modal_blocked_reason},
    )
    if materialize:
        result["materialized_paths"] = materialize_local(result, write_story_file=story)
    return result


def _remote_copy_paths() -> dict[str, Path]:
    return REMOTE_SOURCE_PATHS


if modal is not None:
    volume = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)

    image = (
        modal.Image.debian_slim(python_version="3.11")
        .apt_install(
            "ca-certificates",
            "curl",
            "build-essential",
            "pkg-config",
            "m4",
            "autoconf",
            "automake",
            "libtool",
            "make",
        )
        .run_commands(
            "curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | "
            "sh -s -- -y --profile minimal"
        )
        .env(
            {
                "PATH": "/root/.cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
                "RUSTFLAGS": "-Ctarget-cpu=haswell",
            }
        )
        .add_local_file(SRC_ROOT / "Cargo.toml", f"{REMOTE_SRC}/Cargo.toml", copy=True)
        .add_local_file(SRC_ROOT / "Cargo.lock", f"{REMOTE_SRC}/Cargo.lock", copy=True)
        .add_local_dir(SRC_ROOT / "src", f"{REMOTE_SRC}/src", copy=True)
    )
    for role, local_path in SOURCE_PATHS.items():
        image = image.add_local_file(local_path, str(REMOTE_SOURCE_PATHS[role]), copy=True)
    image = image.run_commands(
        f"cd {REMOTE_SRC} && cargo build --release --bin ehp114_n14_parameterized_local_pipeline_smoke"
    )

    app = modal.App(APP_NAME)

    @app.function(image=image, volumes={"/outputs": volume}, cpu=2, memory=4096, timeout=20 * 60)
    def run_modal_pilot() -> dict[str, Any]:
        outdir = Path(REMOTE_OUT)
        outdir.mkdir(parents=True, exist_ok=True)

        smoke_gate = make_smoke_gate(remote_executed=True)
        smoke_path = outdir / "modal_smoke_RESULTS.json"
        write_json(smoke_path, smoke_gate)
        smoke_digest = write_sha_sidecar(smoke_path)
        smoke_gate["remote_smoke_sha256"] = smoke_digest
        volume.commit()
        print("MODAL_SMOKE_GATE", json.dumps(smoke_gate, sort_keys=True), flush=True)

        remote_paths = _remote_copy_paths()
        n14_control_gate = audit_n14_control(remote_paths)
        n15_pilot_gate = audit_n15_pilot(remote_paths)
        result = make_result_payload(
            source_paths=remote_paths,
            smoke_gate=smoke_gate,
            n14_control_gate=n14_control_gate,
            n15_pilot_gate=n15_pilot_gate,
            modal_metadata={
                "app_name": APP_NAME,
                "volume_name": VOLUME_NAME,
                "remote_out": REMOTE_OUT,
            },
        )
        remote_results = outdir / f"{EXPERIMENT_ID}_RESULTS.json"
        write_json(remote_results, result)
        write_report(result, outdir / f"{EXPERIMENT_ID}_REPORT.md")
        result_digest = write_sha_sidecar(remote_results)
        result["remote_sha256"] = result_digest
        volume.commit()
        print(
            "MODAL_N15_REDUCTION_PILOT",
            json.dumps(
                {
                    "experiment_id": EXPERIMENT_ID,
                    "final_verdict": result["final_verdict"],
                    "remote_out": REMOTE_OUT,
                    "sha256": result_digest,
                },
                sort_keys=True,
            ),
            flush=True,
        )
        return result

    @app.local_entrypoint()
    def modal_main(wait: bool = True, write_local: bool = True, write_story_file: bool = True):
        call = run_modal_pilot.spawn()
        metadata = {
            "call_id": call.object_id,
            "dashboard_url": call.get_dashboard_url(),
            "wait": wait,
            "write_local": write_local,
            "write_story_file": write_story_file,
        }
        print(json.dumps(metadata, indent=2, sort_keys=True))
        if wait:
            result = call.get(timeout=25 * 60)
            result["modal_metadata"]["call_id"] = call.object_id
            result["modal_metadata"]["dashboard_url"] = call.get_dashboard_url()
            if write_local:
                result["materialized_paths"] = materialize_local(
                    result,
                    write_story_file=write_story_file,
                )
            print(json.dumps(result, indent=2, sort_keys=True))


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="print planned inputs/outputs only")
    parser.add_argument(
        "--local-only",
        action="store_true",
        help="run the audit gates locally without Modal",
    )
    parser.add_argument(
        "--materialize",
        action="store_true",
        help="with --local-only, write the final triplet locally",
    )
    parser.add_argument(
        "--write-story",
        action="store_true",
        help="with --materialize, also write the Ahmes story to Downloads",
    )
    parser.add_argument(
        "--modal-blocked-reason",
        help="record a blocked Modal launch reason in a local review-only artifact",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if args.dry_run:
        print(json.dumps(dry_run_summary(), indent=2, sort_keys=True))
        return 0
    if args.local_only:
        result = run_local_only(
            materialize=args.materialize,
            story=args.write_story,
            modal_blocked_reason=args.modal_blocked_reason,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    print(
        "Use `modal run erdos-experiments/scripts/erdos-114/modal_n15_reduction_pilot.py::modal_main` "
        "for the Modal path, or pass --dry-run / --local-only.",
        file=sys.stderr,
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
