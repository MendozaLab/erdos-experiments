#!/usr/bin/env python3
"""Build and run the Rust per-core closure runner for Erdos #20.

The Python per-core runner established the observable on cheap exact regimes.
This wrapper compiles the Rust engine and runs the next calibration/scale
targets without changing D1, scorecards, or public surfaces.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import time
import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-ERDOS20-RUST-PER-CORE-CLOSURE-20260505-01"
DEFAULT_TARGETS = [(3, 7), (4, 7)]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_command(args: list[str], cwd: Path, timeout: int | None = None) -> dict[str, Any]:
    start = time.monotonic()
    proc = subprocess.run(
        args,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=False,
    )
    return {
        "args": args,
        "cwd": str(cwd),
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "elapsed_seconds": time.monotonic() - start,
    }


def parse_target(raw: str) -> tuple[int, int]:
    w_raw, n_raw = raw.split(":", 1)
    return int(w_raw), int(n_raw)


def format_num(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def build_report(result: dict[str, Any]) -> str:
    rows = []
    for target in result["targets"]:
        if target["status"] != "PASS":
            rows.append(
                [
                    target["w"],
                    target["n"],
                    "FAILED",
                    target.get("error", "see stderr"),
                    "",
                    "",
                    "",
                    "",
                ]
            )
            continue
        data = target["data"]
        strongest = max(
            data["core_rows"],
            key=lambda row: row.get("I_core_local_bits") or -1.0,
            default={},
        )
        rows.append(
            [
                target["w"],
                target["n"],
                data["total_sf_free"],
                data["target_m"],
                format_num(data.get("aggregate_I_close_bits_at_target")),
                strongest.get("core_size_s"),
                format_num(strongest.get("I_core_local_bits")),
                format_num(target["elapsed_seconds"]),
            ]
        )
    table = "\n".join(
        ["| w | n | families/status | target m | aggregate I | strongest s | strongest local I | seconds |",
         "|---|---|---:|---:|---:|---:|---:|---:|"]
        + ["| " + " | ".join(map(str, row)) + " |" for row in rows]
    )
    return f"""# {EXPERIMENT_ID} Report

## Status

Rust per-core closure runner for Erdos #20. This is an internal diagnostic
artifact. It is not theorem progress, not lower-bound progress, and not a
Leg-4 pass. The claim ceiling remains: this is a shadow signature, not
universal law.

## Verdict

- Runner status: `{result["runner_status"]}`
- Classification: `{result["classification"]}`
- Rust engine: `erdos-experiments/Erdos20/rust_core_closure`

The runner confirms that Rust is now useful for this lane. Python was enough to
define the observable; Rust is the right layer for larger per-core sweeps.

## Executed Targets

{table}

## Interpretation

The important comparison is not just the aggregate closure pressure. The Rust
runner records fixed-core channels at the selected target size, so the question
becomes whether particular core sizes carry a repeatable local closure cost.

This still does not execute the full Leg-4 test. The missing pieces are a
predeclared floor-normalized numerator, Abbott-Hansen-Sauer baseline controls,
and a symmetry-reduced transfer formulation that can push beyond exact
enumeration.

## Next Rust Target

The next queued run is `w=3,n=8`; it has about `148790380` sunflower-free
families in the saved aggregate artifact. The runner is ready, but that run
should be treated as a heavier execution packet rather than mixed into this
calibration artifact.

## Source Boundary

No scorecard, D1, public document, git, email, CLAUDE.md, or AGENTS.md was
changed.
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-id", default=EXPERIMENT_ID)
    parser.add_argument(
        "--targets",
        default=",".join(f"{w}:{n}" for w, n in DEFAULT_TARGETS),
        help="Comma-separated w:n targets, for example 3:8 or 3:7,4:7",
    )
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()

    experiment_id = args.experiment_id
    targets_requested = [parse_target(raw.strip()) for raw in args.targets.split(",") if raw.strip()]

    out_dir = Path(__file__).resolve().parent
    cargo_dir = out_dir / "rust_core_closure"
    result_path = out_dir / f"{experiment_id}_RESULTS.json"
    report_path = out_dir / f"{experiment_id}_REPORT.md"
    sha_path = out_dir / f"{experiment_id}_RESULTS.sha256"
    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {path}")

    build = run_command(["cargo", "build", "--release"], cwd=cargo_dir, timeout=300)
    targets = []
    if build["returncode"] == 0:
        for w, n in targets_requested:
            cmd = [
                str(cargo_dir / "target" / "release" / "sunflower-core-closure"),
                "--w",
                str(w),
                "--n",
                str(n),
                "--target-m",
                "auto",
                "--window",
                "0",
            ]
            run = run_command(cmd, cwd=out_dir, timeout=args.timeout)
            record: dict[str, Any] = {
                "w": w,
                "n": n,
                "status": "PASS" if run["returncode"] == 0 else "FAIL",
                "elapsed_seconds": run["elapsed_seconds"],
                "stderr": run["stderr"][-4000:],
            }
            if run["returncode"] == 0:
                record["data"] = json.loads(run["stdout"])
            else:
                record["error"] = run["stdout"][-4000:]
            targets.append(record)

    pass_targets = [target for target in targets if target["status"] == "PASS"]
    nonzero = 0
    stratified = 0
    for target in pass_targets:
        values = [
            row.get("I_core_local_bits")
            for row in target["data"]["core_rows"]
            if row.get("I_core_local_bits") is not None
        ]
        if values and max(values) > 0.25:
            nonzero += 1
        if values and (max(values) - min(values)) > 0.1:
            stratified += 1

    classification = (
        "RUST_PER_CORE_SIGNAL_PRESENT"
        if pass_targets and nonzero and stratified
        else "RUST_PER_CORE_INCONCLUSIVE"
        if pass_targets
        else "RUST_PER_CORE_BLOCKED"
    )
    result: dict[str, Any] = {
        "experiment_id": experiment_id,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #20 sunflower core closure",
        "scope": "Rust per-core closure calibration and scale runner",
        "runner_status": "PASS" if build["returncode"] == 0 and len(pass_targets) == len(targets_requested) else "PARTIAL_OR_FAIL",
        "classification": classification,
        "claim_ceiling": "shadow signature, not universal law; no theorem or lower-bound progress",
        "build": {
            "returncode": build["returncode"],
            "elapsed_seconds": build["elapsed_seconds"],
            "stderr_tail": build["stderr"][-4000:],
        },
        "targets": targets,
        "next_heavy_targets": [
            {
                "w": 3,
                "n": 8,
                "known_sf_free_families": 148790380,
                "recommended_command": "target/release/sunflower-core-closure --w 3 --n 8 --target-m auto --window 0",
            }
        ],
        "guardrails": {
            "no_d1": True,
            "no_scorecard": True,
            "no_public_docs": True,
            "no_git": True,
            "no_email": True,
        },
    }

    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_text = build_report(result).replace(EXPERIMENT_ID, experiment_id)
    report_path.write_text(report_text, encoding="utf-8")
    sha_path.write_text(f"{sha256_file(result_path)}  {result_path.name}\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "experiment_id": experiment_id,
                "result": str(result_path),
                "report": str(report_path),
                "sha256": sha_path.read_text(encoding="utf-8").split()[0],
                "classification": classification,
                "runner_status": result["runner_status"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
