#!/usr/bin/env python3
"""Modal runner for the EHP114 n=15 CELL-02-03 boundary-slice pilot.

This is a review-only reduction experiment. It runs one degree-15 collar probe
against the hard CELL-02-03 residual boxes inherited from the n=14 audit chain.
It does not launch full n=15 certification and does not update public or proof
state.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

try:
    import modal
except Exception:  # pragma: no cover - dry-run works without Modal.
    modal = None  # type: ignore[assignment]


EXPERIMENT_ID = "EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-MODAL-20260508-01"
APP_NAME = "ehp114-n15-boundary-slice-cell-02-03-20260508-01"
VOLUME_NAME = "ehp114-n15-boundary-slice-cell-02-03-20260508-01"

REMOTE_SRC = "/workspace/erdos-114"
REMOTE_INPUTS = "/workspace/inputs"
REMOTE_OUT = f"/outputs/{EXPERIMENT_ID}"
REMOTE_SOURCE = Path(REMOTE_INPUTS) / "n14_cell_02_03_third_order_source.json"

SCRIPT_ROOT = Path(__file__).resolve().parent
SRC_ROOT = SCRIPT_ROOT if (SCRIPT_ROOT / "Cargo.toml").exists() else Path(REMOTE_SRC)
ERDOS_EXPERIMENTS = (
    SCRIPT_ROOT.parents[1]
    if len(SCRIPT_ROOT.parents) > 1 and (SCRIPT_ROOT.parents[1] / "Erdos114").exists()
    else Path("/workspace")
)

LOCAL_OUTDIR = ERDOS_EXPERIMENTS / "Erdos114" / "proof_path" / EXPERIMENT_ID
DOWNLOADS_STORY_DIR = Path.home() / "Downloads" / "Eratosthenes_Ahmes_Stories"

SOURCE_PATHS = {
    "n14_cell_02_03_third_order_source": ERDOS_EXPERIMENTS
    / "Erdos114"
    / "validated_length"
    / "EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-02-03-20260506-01"
    / "EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-02-03-20260506-01_RESULTS.json",
    "n15_fourier_hessian_context": ERDOS_EXPERIMENTS
    / "results"
    / "erdos-114"
    / "EXP-MATH-EHP114-N15-FOURIER-HESSIAN-20260502-02_RESULTS.json",
    "n15_local_slice_tao_bridge_context": ERDOS_EXPERIMENTS
    / "Erdos114"
    / "proof_path"
    / "EXP-MATH-EHP114-N15-LOCAL-SLICE-TAO-BRIDGE-20260508-02"
    / "EXP-MATH-EHP114-N15-LOCAL-SLICE-TAO-BRIDGE-20260508-02_RESULTS.json",
    "n15_reduction_modal_context": ERDOS_EXPERIMENTS
    / "Erdos114"
    / "proof_path"
    / "EXP-MATH-EHP114-N15-REDUCTION-PILOT-MODAL-REMOTE-PASS-20260508-01"
    / "EXP-MATH-EHP114-N15-REDUCTION-PILOT-MODAL-20260508-01_RESULTS.json",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def source_fingerprints() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for role, path in sorted(SOURCE_PATHS.items()):
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


def validate_sources() -> None:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("missing source inputs: " + ", ".join(missing))


def dry_run_summary(region_limit: int) -> dict[str, Any]:
    validate_sources()
    return {
        "experiment_id": EXPERIMENT_ID,
        "mode": "dry_run",
        "app_name": APP_NAME,
        "volume_name": VOLUME_NAME,
        "remote_out": REMOTE_OUT,
        "region_limit": region_limit,
        "would_write": [
            str(LOCAL_OUTDIR / f"{EXPERIMENT_ID}_RESULTS.json"),
            str(LOCAL_OUTDIR / f"{EXPERIMENT_ID}_REPORT.md"),
            str(LOCAL_OUTDIR / f"{EXPERIMENT_ID}_RESULTS.sha256"),
            str(DOWNLOADS_STORY_DIR / f"AHMES_STORY_{EXPERIMENT_ID}.md"),
        ],
        "source_fingerprints": source_fingerprints(),
    }


def write_story(result: dict[str, Any], call_metadata: dict[str, Any]) -> Path:
    DOWNLOADS_STORY_DIR.mkdir(parents=True, exist_ok=True)
    story_path = DOWNLOADS_STORY_DIR / f"AHMES_STORY_{EXPERIMENT_ID}.md"
    if story_path.exists():
        raise FileExistsError(f"refusing to overwrite existing story: {story_path}")
    story = f"""# Ahmes Story: {EXPERIMENT_ID}

Eratosthenes took the hard CELL-02-03 collar slice to the accelerator, but not as
a full degree-15 launch. The task was narrower: replay the troublesome n=14
residual boxes as a degree-15 Fourier/root-slice boundary test and ask whether
the local collar machinery still has traction.

The verdict was `{result.get("final_verdict")}`. That is a reduction signal,
not public state. The next move is controlled by the blocker recorded in the
JSON packet, not by enthusiasm for a full n=15 run.

## Receipts

- Run ID: `{EXPERIMENT_ID}`
- Modal app: `{call_metadata.get("app_name")}`
- Modal call ID: `{call_metadata.get("call_id")}`
- Modal dashboard: `{call_metadata.get("dashboard_url")}`
- Modal volume: `{VOLUME_NAME}`
- Remote output: `{REMOTE_OUT}`
- Local mirror: `{LOCAL_OUTDIR}`
- Result SHA-256: `{call_metadata.get("local_sha256")}`
- Final verdict: `{result.get("final_verdict")}`
- Collar status: `{result.get("collar_status")}`
- Processed regions: `{result.get("processed_region_count")}`
- Remaining unresolved regions: `{result.get("remaining_unresolved_region_count")}`

## Review Warning

This story is a reading companion. The source of truth is the JSON/report/hash
triplet in the run folder. No D1, public page, curated morphism, existing #114
packet, or formal-status registry was updated.
"""
    story_path.write_text(story, encoding="utf-8")
    return story_path


def materialize_local(
    remote_payload: dict[str, Any],
    *,
    call_metadata: dict[str, Any],
    write_story_file: bool,
) -> dict[str, str]:
    if LOCAL_OUTDIR.exists():
        raise FileExistsError(f"refusing to overwrite existing run folder: {LOCAL_OUTDIR}")
    LOCAL_OUTDIR.mkdir(parents=True)
    result = remote_payload["result"]
    report_text = remote_payload["report_text"]

    result_path = LOCAL_OUTDIR / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = LOCAL_OUTDIR / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = LOCAL_OUTDIR / f"{EXPERIMENT_ID}_RESULTS.sha256"

    write_json(result_path, result)
    report_path.write_text(report_text, encoding="utf-8")
    digest = sha256_file(result_path)
    if digest != remote_payload["result_sha256"]:
        raise RuntimeError(
            f"local result hash {digest} does not match remote {remote_payload['result_sha256']}"
        )
    sha_path.write_text(f"{digest}  {result_path.name}\n", encoding="utf-8")
    call_metadata["local_sha256"] = digest
    story_path = write_story(result, call_metadata) if write_story_file else None
    return {
        "results": str(result_path),
        "report": str(report_path),
        "sha256": str(sha_path),
        "digest": digest,
        "story": str(story_path) if story_path else "",
    }


if modal is not None:
    volume = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)
    image = (
        modal.Image.debian_slim(python_version="3.13")
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
        .add_local_file(SOURCE_PATHS["n14_cell_02_03_third_order_source"], str(REMOTE_SOURCE), copy=True)
        .run_commands(
            f"cd {REMOTE_SRC} && cargo build --release --bin ehp114_n15_boundary_slice_cell_02_03"
        )
    )

    app = modal.App(APP_NAME)

    @app.function(image=image, volumes={"/outputs": volume}, cpu=4, memory=8192, timeout=2 * 60 * 60)
    def run_boundary_slice(region_limit: int) -> dict[str, Any]:
        outdir = Path(REMOTE_OUT)
        outdir.mkdir(parents=True, exist_ok=True)

        smoke = {
            "gate": "modal_smoke",
            "status": "PASS",
            "timestamp_unix": int(time.time()),
            "volume_write_checked": True,
            "sidecar_hash_checked": True,
            "meaning": "Tiny Modal output-contract gate before the boundary-slice runner.",
        }
        smoke_path = outdir / "modal_smoke_RESULTS.json"
        write_json(smoke_path, smoke)
        smoke_sha = sha256_file(smoke_path)
        (outdir / "modal_smoke_RESULTS.sha256").write_text(
            f"{smoke_sha}  {smoke_path.name}\n",
            encoding="utf-8",
        )
        print("MODAL_SMOKE_GATE", json.dumps(smoke, sort_keys=True), flush=True)

        cmd = [
            f"{REMOTE_SRC}/target/release/ehp114_n15_boundary_slice_cell_02_03",
            "--source",
            str(REMOTE_SOURCE),
            "--outdir",
            REMOTE_OUT,
            "--experiment-id",
            EXPERIMENT_ID,
            "--sub-i",
            "2",
            "--sub-j",
            "3",
            "--region-limit",
            str(region_limit),
        ]
        completed = subprocess.run(cmd, cwd=REMOTE_SRC, check=True, text=True, capture_output=True)
        print(completed.stdout, flush=True)
        if completed.stderr:
            print(completed.stderr, flush=True)

        result_path = outdir / f"{EXPERIMENT_ID}_RESULTS.json"
        report_path = outdir / f"{EXPERIMENT_ID}_REPORT.md"
        sha_path = outdir / f"{EXPERIMENT_ID}_RESULTS.sha256"
        result = json.loads(result_path.read_text(encoding="utf-8"))
        report_text = report_path.read_text(encoding="utf-8")
        result_sha = sha256_file(result_path)
        expected_sha = sha_path.read_text(encoding="utf-8").split()[0]
        if result_sha != expected_sha:
            raise RuntimeError(f"remote sha mismatch: {result_sha} != {expected_sha}")
        volume.commit()
        print(
            "MODAL_N15_BOUNDARY_SLICE",
            json.dumps(
                {
                    "experiment_id": EXPERIMENT_ID,
                    "final_verdict": result.get("final_verdict"),
                    "collar_status": result.get("collar_status"),
                    "remote_out": REMOTE_OUT,
                    "sha256": result_sha,
                },
                sort_keys=True,
            ),
            flush=True,
        )
        return {
            "result": result,
            "report_text": report_text,
            "result_sha256": result_sha,
            "remote_paths": {
                "results": str(result_path),
                "report": str(report_path),
                "sha256": str(sha_path),
                "smoke": str(smoke_path),
            },
        }

    @app.local_entrypoint()
    def modal_main(
        wait: bool = True,
        write_local: bool = True,
        write_story_file: bool = True,
        region_limit: int = 0,
    ):
        call = run_boundary_slice.spawn(region_limit)
        call_metadata = {
            "app_name": APP_NAME,
            "call_id": call.object_id,
            "dashboard_url": call.get_dashboard_url(),
            "wait": wait,
            "write_local": write_local,
            "write_story_file": write_story_file,
            "region_limit": region_limit,
        }
        print(json.dumps(call_metadata, indent=2, sort_keys=True))
        if wait:
            payload = call.get(timeout=2 * 60 * 60 + 10 * 60)
            if write_local:
                payload["materialized_paths"] = materialize_local(
                    payload,
                    call_metadata=call_metadata,
                    write_story_file=write_story_file,
                )
            print(json.dumps(payload, indent=2, sort_keys=True))


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="print planned inputs/outputs only")
    parser.add_argument("--region-limit", type=int, default=0, help="0 means process all source regions")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if args.dry_run:
        print(json.dumps(dry_run_summary(args.region_limit), indent=2, sort_keys=True))
        return 0
    print(
        "Use `modal run erdos-experiments/scripts/erdos-114/modal_n15_boundary_slice_cell_02_03.py::modal_main` "
        "for the Modal path, or pass --dry-run.",
        file=sys.stderr,
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
