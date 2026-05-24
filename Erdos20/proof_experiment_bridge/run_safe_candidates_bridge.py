#!/usr/bin/env python3
"""Candidate-level safeCandidates bridge run for Erdos #20.

This binds the formal `safeCandidates` slot to real enumerator rows for one
small exact regime. It is a bridge artifact only: no theorem, no lower-bound
claim, no Leg-4 pass.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-ERDOS20-SAFE-CANDIDATES-BRIDGE-20260506-02"
W = 3
N = 7
K = 3
CORE_SIZE_S = 2
M_DELTA = 1
I_FLOOR_BITS = 1.0
I_AHS_BITS = math.log2(math.sqrt(10.0))


BASE = Path(__file__).resolve().parent
ERDOS20 = BASE.parent
sys.path.insert(0, str(ERDOS20))

from sunflower_per_core_closure_analysis import (  # noqa: E402
    PerCoreEnumerator,
    bit_count,
    mask_tuple,
    round_finite,
    safe_log2_ratio,
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_mask_indices(mask: int) -> list[int]:
    out: list[int] = []
    idx = 0
    while mask:
        if mask & 1:
            out.append(idx)
        mask >>= 1
        idx += 1
    return out


def first_unsafe_witnesses(
    enumerator: PerCoreEnumerator,
    family_indices: list[int],
    candidate_idx: int,
    core_idx: int,
    limit: int = 3,
) -> list[dict[str, Any]]:
    witnesses: list[dict[str, Any]] = []
    candidate_bit = 1 << candidate_idx
    for pos, left_idx in enumerate(family_indices):
        for right_idx in family_indices[pos + 1 :]:
            if enumerator.pair_core_index[left_idx][right_idx] != core_idx:
                continue
            if not (enumerator.pair_block_mask[left_idx][right_idx] & candidate_bit):
                continue
            witnesses.append(
                {
                    "left_site_index": left_idx,
                    "left_site": mask_tuple(enumerator.sites[left_idx], enumerator.n),
                    "right_site_index": right_idx,
                    "right_site": mask_tuple(enumerator.sites[right_idx], enumerator.n),
                    "candidate_site_index": candidate_idx,
                    "candidate_site": mask_tuple(enumerator.sites[candidate_idx], enumerator.n),
                    "core": mask_tuple(enumerator.cores[core_idx], enumerator.n),
                }
            )
            if len(witnesses) >= limit:
                return witnesses
    return witnesses


def build_candidate_rows(
    enumerator: PerCoreEnumerator,
    family_indices: list[int],
    family_mask: int,
    blocked_global: int,
    blocked_by_core: list[int],
    core_idx: int,
) -> dict[str, Any]:
    unused_mask = enumerator.all_site_mask & ~family_mask
    candidate_mask = enumerator.core_candidate_masks[core_idx] & unused_mask
    safe_mask = candidate_mask & ~blocked_by_core[core_idx]
    global_safe_mask = candidate_mask & ~blocked_global
    candidate_count = bit_count(candidate_mask)
    safe_count = bit_count(safe_mask)
    global_safe_count = bit_count(global_safe_mask)
    i_core = safe_log2_ratio(safe_count, candidate_count)
    family_sites = [
        {
            "site_index": idx,
            "site": mask_tuple(enumerator.sites[idx], enumerator.n),
        }
        for idx in family_indices
    ]
    rows = []
    for candidate_idx in iter_mask_indices(candidate_mask):
        candidate_bit = 1 << candidate_idx
        local_safe = bool(safe_mask & candidate_bit)
        global_safe = bool(global_safe_mask & candidate_bit)
        rows.append(
            {
                "candidate_site_index": candidate_idx,
                "candidate_extension": mask_tuple(enumerator.sites[candidate_idx], enumerator.n),
                "safe_extension": local_safe,
                "global_safe_extension": global_safe,
                "unsafe_witness_pairs": []
                if local_safe
                else first_unsafe_witnesses(enumerator, family_indices, candidate_idx, core_idx),
            }
        )
    return {
        "core_index": core_idx,
        "core": mask_tuple(enumerator.cores[core_idx], enumerator.n),
        "core_mask": enumerator.cores[core_idx],
        "core_size_s": enumerator.core_sizes[core_idx],
        "family_sites": family_sites,
        "candidate_count": candidate_count,
        "safe_count": safe_count,
        "global_safe_count": global_safe_count,
        "I_core_local_bits": round_finite(i_core),
        "I_floor_bits": I_FLOOR_BITS,
        "floor_ratio": round_finite(i_core / I_FLOOR_BITS if i_core is not None else None),
        "I_AHS_bits": round_finite(I_AHS_BITS),
        "ahs_ratio": round_finite(i_core / I_AHS_BITS if i_core is not None else None),
        "safeCandidates": [
            row["candidate_extension"] for row in rows if row["safe_extension"]
        ],
        "candidate_rows": rows,
    }


def choose_core_sample(
    enumerator: PerCoreEnumerator,
    family_indices: list[int],
    family_mask: int,
    blocked_global: int,
    blocked_by_core: list[int],
) -> dict[str, Any] | None:
    best: tuple[float, int, dict[str, Any]] | None = None
    for core_idx in enumerator.core_size_to_indices[CORE_SIZE_S]:
        record = build_candidate_rows(
            enumerator, family_indices, family_mask, blocked_global, blocked_by_core, core_idx
        )
        if record["candidate_count"] <= 0:
            continue
        if record["safe_count"] <= 0:
            continue
        if record["safe_count"] >= record["candidate_count"]:
            continue
        # Require nontrivial local closure: both safe and unsafe candidates exist.
        score = record["I_core_local_bits"] or 0.0
        if best is None or score > best[0]:
            best = (score, core_idx, record)
    return None if best is None else best[2]


def sample_families(enumerator: PerCoreEnumerator, target_ms: set[int]) -> dict[str, Any]:
    samples: dict[int, dict[str, Any]] = {}
    family_indices: list[int] = []
    family_mask = 0
    blocked_global = 0
    blocked_by_core = [0 for _ in enumerator.cores]

    def record_if_needed() -> None:
        m = len(family_indices)
        if m not in target_ms or m in samples:
            return
        sample = choose_core_sample(
            enumerator, family_indices, family_mask, blocked_global, blocked_by_core
        )
        if sample is not None:
            samples[m] = {
                "m": m,
                "family_size": m,
                **sample,
            }

    def backtrack(start: int) -> bool:
        nonlocal family_mask, blocked_global
        record_if_needed()
        if target_ms.issubset(set(samples)):
            return True
        for site_idx in range(start, enumerator.N):
            site_bit = 1 << site_idx
            if family_mask & site_bit:
                continue
            if blocked_global & site_bit:
                continue
            old_family_mask = family_mask
            old_blocked_global = blocked_global
            changed_cores: list[tuple[int, int]] = []
            new_blocked_global = blocked_global
            for existing_idx in family_indices:
                core_idx = enumerator.pair_core_index[existing_idx][site_idx]
                block_mask = enumerator.pair_block_mask[existing_idx][site_idx]
                old_core_mask = blocked_by_core[core_idx]
                new_core_mask = old_core_mask | block_mask
                if new_core_mask != old_core_mask:
                    changed_cores.append((core_idx, old_core_mask))
                    blocked_by_core[core_idx] = new_core_mask
                new_blocked_global |= block_mask
            family_indices.append(site_idx)
            family_mask = old_family_mask | site_bit
            blocked_global = new_blocked_global
            if backtrack(site_idx + 1):
                return True
            blocked_global = old_blocked_global
            family_mask = old_family_mask
            family_indices.pop()
            for core_idx, old_core_mask in reversed(changed_cores):
                blocked_by_core[core_idx] = old_core_mask
        return False

    backtrack(0)
    return {str(key): samples[key] for key in sorted(samples)}


def build_report(result: dict[str, Any]) -> str:
    sample_lines = []
    for m, sample in result["samples_by_m"].items():
        sample_lines.append(
            f"| {m} | {sample['core']} | {sample['candidate_count']} | "
            f"{sample['safe_count']} | {sample['I_core_local_bits']} | "
            f"{sample['floor_ratio']} | {sample['ahs_ratio']} |"
        )
    sample_table = "\n".join(
        [
            "| m | core | candidates | safe | I_core_local | floor ratio | AHS ratio |",
            "|---:|---|---:|---:|---:|---:|---:|",
            *sample_lines,
        ]
    )
    return f"""# {EXPERIMENT_ID} Report

## Status

- Status: `{result['status']}`
- Target: `w={W}, n={N}, k={K}, core_size_s={CORE_SIZE_S}`
- Claim ceiling: A0, shadow signature, not universal law

## What This Binds

This bridge emits candidate-level rows for the Lean-shaped relation:

```text
safeCandidates = candidateExtensions.filter (SafeExtension core family)
```

Each unsafe candidate includes concrete witness pairs from the family when
available. This is still bookkeeping evidence, not theorem progress.

## Ratios

The run reports both controls requested after third-party review:

- `floor_ratio = I_core_local / I_floor`, with `I_floor = 1 bit`
- `ahs_ratio = I_core_local / log2(sqrt(10))`, using the Q1 gate's
  Abbott-Hansen-Sauer base `c_3 >= sqrt(10)` as a construction/literature
  control, not an exact local construction model

{sample_table}

## Interpretation

This turns the prior aggregate/per-core summaries into proof-facing rows. It
does not reproduce an Abbott-Hansen-Sauer construction, does not prove a
sunflower theorem, and does not upgrade the A-axis.
"""


def main() -> int:
    result_path = BASE / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = BASE / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = BASE / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {path}")

    start = time.monotonic()
    enumerator = PerCoreEnumerator(n=N, w=W, k=K)
    full_run = enumerator.run()
    m_star = full_run["jamming_m_star"]
    if m_star is None:
        target_ms = {full_run["summary"]["target_m"]}
    else:
        target_ms = {
            m for m in (m_star - M_DELTA, m_star, m_star + M_DELTA) if m >= 0
        }
    samples = sample_families(enumerator, target_ms)
    status = (
        "SAFE_CANDIDATES_BRIDGE_PASS_A0"
        if str(m_star) in samples and samples
        else "SAFE_CANDIDATES_BRIDGE_FAIL_NO_TARGET_SAMPLE"
    )
    result: dict[str, Any] = {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "problem": "Erdos #20 sunflower conjecture",
        "scope": "candidate-level safeCandidates bridge sample; not theorem progress",
        "parameters": {
            "w": W,
            "n": N,
            "k": K,
            "core_size_s": CORE_SIZE_S,
            "m_delta": M_DELTA,
            "target_ms": sorted(target_ms),
            "jamming_m_star": m_star,
        },
        "normalizations": {
            "I_floor_bits": I_FLOOR_BITS,
            "I_floor_definition": "one predeclared bit of local closure information; preliminary Leg-4 denominator",
            "I_AHS_bits": round_finite(I_AHS_BITS),
            "I_AHS_definition": "log2(sqrt(10)); Q1 literature-gate construction/base control, not a reproduced AHS local construction",
        },
        "samples_by_m": samples,
        "full_run_summary": full_run["summary"],
        "claim_ceiling": "A0; shadow signature, not universal law; no theorem or lower-bound progress",
        "elapsed_seconds": round_finite(time.monotonic() - start, digits=3),
    }
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_path.write_text(build_report(result), encoding="utf-8")
    sha_path.write_text(f"{sha256_file(result_path)}  {result_path.name}\n", encoding="utf-8")
    print(json.dumps({
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "target_ms": sorted(target_ms),
        "sample_count": len(samples),
        "result": str(result_path),
        "sha256": sha256_file(result_path),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
