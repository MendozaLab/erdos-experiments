#!/usr/bin/env python3
"""Per-core closure instrumentation for Erdos #20 sunflower-free families.

The prior Hessian run measured aggregate density-of-states curvature. This
script re-enumerates cheap exact calibration targets and records the missing
core channel: for each core C, how many unused petal extensions through C stay
locally open, and how many stay globally safe.

This is an internal diagnostic. It is not a theorem, lower-bound result, or
public Leg-4 pass.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import date
from itertools import combinations
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-ERDOS20-PER-CORE-CLOSURE-20260505-01"
RUN_DATE = "2026-05-05"

DEFAULT_TARGETS = [
    (2, 6),
    (2, 7),
    (2, 8),
    (2, 9),
    (2, 10),
    (3, 4),
    (3, 5),
    (3, 6),
    (3, 7),
    (4, 5),
    (4, 6),
]

KNOWN_EXPENSIVE_TARGETS = [
    {
        "w": 3,
        "n": 8,
        "known_sf_free_families": 148790380,
        "reason": "per-core instrumentation over 148790380 exact families is not cheap in Python",
    },
    {
        "w": 4,
        "n": 7,
        "known_sf_free_families": 35333735,
        "reason": "per-core instrumentation over 35333735 exact families is not cheap in Python",
    },
    {
        "w": 4,
        "n": 8,
        "known_sf_free_families": None,
        "reason": "the prior exhaustive aggregate run did not finish for n=8,w=4",
    },
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def round_finite(value: float | None, digits: int = 6) -> float | None:
    if value is None or not math.isfinite(value):
        return None
    rounded = round(value, digits)
    return 0.0 if rounded == 0 else rounded


def safe_log2_ratio(numerator: int | float, denominator: int | float) -> float | None:
    if denominator <= 0 or numerator <= 0:
        return None
    return -math.log2(numerator / denominator)


def mask_from_indices(indices: tuple[int, ...]) -> int:
    mask = 0
    for idx in indices:
        mask |= 1 << idx
    return mask


def mask_tuple(mask: int, n: int) -> list[int]:
    return [idx + 1 for idx in range(n) if mask & (1 << idx)]


def bit_count(mask: int) -> int:
    return mask.bit_count()


def quantiles_from_distribution(
    distribution: dict[tuple[int, int], int],
) -> dict[str, Any]:
    finite_values: list[tuple[float, int]] = []
    finite_count = 0
    blocked_count = 0
    total_count = 0
    sum_bits = 0.0

    for (candidate_count, valid_count), count in distribution.items():
        total_count += count
        if candidate_count <= 0:
            continue
        if valid_count <= 0:
            blocked_count += count
            continue
        bits = -math.log2(valid_count / candidate_count)
        finite_values.append((bits, count))
        finite_count += count
        sum_bits += bits * count

    finite_values.sort(key=lambda item: item[0])

    def quantile(q: float) -> float | None:
        if finite_count == 0:
            return None
        threshold = max(1, math.ceil(q * finite_count))
        running = 0
        for bits, count in finite_values:
            running += count
            if running >= threshold:
                return bits
        return finite_values[-1][0]

    return {
        "finite_channel_count": finite_count,
        "blocked_channel_count": blocked_count,
        "total_channel_count": total_count,
        "blocked_channel_fraction": round_finite(
            blocked_count / total_count if total_count else None
        ),
        "mean_finite_I_core_bits": round_finite(
            sum_bits / finite_count if finite_count else None
        ),
        "median_finite_I_core_bits": round_finite(quantile(0.5)),
        "p90_finite_I_core_bits": round_finite(quantile(0.9)),
    }


@dataclass
class ChannelAccumulator:
    channel_count: int = 0
    candidate_count: int = 0
    local_valid_count: int = 0
    global_valid_count: int = 0
    distribution: dict[tuple[int, int], int] = field(default_factory=lambda: defaultdict(int))

    def add(self, candidates: int, local_valid: int, global_valid: int) -> None:
        if candidates <= 0:
            return
        self.channel_count += 1
        self.candidate_count += candidates
        self.local_valid_count += local_valid
        self.global_valid_count += global_valid
        self.distribution[(candidates, local_valid)] += 1

    def to_record(self, n: int, w: int, m: int, s: int) -> dict[str, Any]:
        local_i = safe_log2_ratio(self.local_valid_count, self.candidate_count)
        global_i = safe_log2_ratio(self.global_valid_count, self.candidate_count)
        quantiles = quantiles_from_distribution(self.distribution)
        return {
            "n": n,
            "w": w,
            "m": m,
            "core_size_s": s,
            "petal_size": w - s,
            "channel_count": self.channel_count,
            "candidate_count": self.candidate_count,
            "local_valid_count": self.local_valid_count,
            "global_valid_count": self.global_valid_count,
            "local_p_safe_weighted": round_finite(
                self.local_valid_count / self.candidate_count
                if self.candidate_count
                else None
            ),
            "global_p_safe_weighted": round_finite(
                self.global_valid_count / self.candidate_count
                if self.candidate_count
                else None
            ),
            "I_core_local_weighted_bits": round_finite(local_i),
            "I_core_global_weighted_bits": round_finite(global_i),
            **quantiles,
        }


class PerCoreEnumerator:
    def __init__(self, n: int, w: int, k: int = 3) -> None:
        if k != 3:
            raise ValueError("Only k=3 is implemented for this per-core runner")
        self.n = n
        self.w = w
        self.k = k
        self.sites = [mask_from_indices(combo) for combo in combinations(range(n), w)]
        self.N = len(self.sites)
        self.all_site_mask = (1 << self.N) - 1
        self.cores: list[int] = []
        self.core_index: dict[int, int] = {}
        self.core_sizes: list[int] = []
        self.core_candidate_masks: list[int] = []
        self.core_size_to_indices: dict[int, list[int]] = defaultdict(list)
        self.pair_core_index = [[-1] * self.N for _ in range(self.N)]
        self.pair_block_mask = [[0] * self.N for _ in range(self.N)]
        self._build_cores()
        self._build_pair_blockers()

    def _build_cores(self) -> None:
        for s in range(self.w):
            for combo in combinations(range(self.n), s):
                core = mask_from_indices(combo)
                idx = len(self.cores)
                self.cores.append(core)
                self.core_index[core] = idx
                self.core_sizes.append(s)
                self.core_size_to_indices[s].append(idx)

                candidates = 0
                for site_idx, site in enumerate(self.sites):
                    if (site & core) == core:
                        candidates |= 1 << site_idx
                self.core_candidate_masks.append(candidates)

    def _build_pair_blockers(self) -> None:
        for left in range(self.N):
            left_site = self.sites[left]
            for right in range(left + 1, self.N):
                right_site = self.sites[right]
                core = left_site & right_site
                core_idx = self.core_index[core]
                block_mask = 0
                for candidate_idx, candidate in enumerate(self.sites):
                    if candidate_idx == left or candidate_idx == right:
                        continue
                    if (
                        (candidate & core) == core
                        and (candidate & left_site) == core
                        and (candidate & right_site) == core
                    ):
                        block_mask |= 1 << candidate_idx
                self.pair_core_index[left][right] = core_idx
                self.pair_core_index[right][left] = core_idx
                self.pair_block_mask[left][right] = block_mask
                self.pair_block_mask[right][left] = block_mask

    def run(self) -> dict[str, Any]:
        max_m_possible = self.N
        density = [0 for _ in range(max_m_possible + 1)]
        all_extension_candidates = [0 for _ in range(max_m_possible + 1)]
        all_extension_valid = [0 for _ in range(max_m_possible + 1)]
        channel_stats: list[list[ChannelAccumulator]] = [
            [ChannelAccumulator() for _ in range(self.w)] for _ in range(max_m_possible + 1)
        ]
        per_core_totals: dict[tuple[int, int], ChannelAccumulator] = {
            (m, core_idx): ChannelAccumulator()
            for m in range(max_m_possible + 1)
            for core_idx in range(len(self.cores))
        }

        family_indices: list[int] = []
        family_mask = 0
        blocked_global = 0
        blocked_by_core = [0 for _ in self.cores]

        def record_family() -> None:
            m = len(family_indices)
            density[m] += 1
            unused_mask = self.all_site_mask & ~family_mask
            global_valid_mask = unused_mask & ~blocked_global
            all_extension_candidates[m] += bit_count(unused_mask)
            all_extension_valid[m] += bit_count(global_valid_mask)

            for core_idx, candidate_mask in enumerate(self.core_candidate_masks):
                candidates_mask = candidate_mask & unused_mask
                candidate_count = bit_count(candidates_mask)
                if candidate_count == 0:
                    continue
                local_valid_count = bit_count(candidates_mask & ~blocked_by_core[core_idx])
                global_valid_count = bit_count(candidates_mask & ~blocked_global)
                s = self.core_sizes[core_idx]
                channel_stats[m][s].add(
                    candidate_count,
                    local_valid_count,
                    global_valid_count,
                )
                per_core_totals[(m, core_idx)].add(
                    candidate_count,
                    local_valid_count,
                    global_valid_count,
                )

        def backtrack(start: int) -> None:
            nonlocal family_mask, blocked_global
            record_family()

            for site_idx in range(start, self.N):
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
                    core_idx = self.pair_core_index[existing_idx][site_idx]
                    block_mask = self.pair_block_mask[existing_idx][site_idx]
                    old_core_mask = blocked_by_core[core_idx]
                    new_core_mask = old_core_mask | block_mask
                    if new_core_mask != old_core_mask:
                        changed_cores.append((core_idx, old_core_mask))
                        blocked_by_core[core_idx] = new_core_mask
                    new_blocked_global |= block_mask

                family_indices.append(site_idx)
                family_mask = old_family_mask | site_bit
                blocked_global = new_blocked_global
                backtrack(site_idx + 1)

                blocked_global = old_blocked_global
                family_mask = old_family_mask
                family_indices.pop()
                for core_idx, old_core_mask in reversed(changed_cores):
                    blocked_by_core[core_idx] = old_core_mask

        t0 = time.perf_counter()
        backtrack(0)
        elapsed = time.perf_counter() - t0

        max_family_size = max(idx for idx, count in enumerate(density) if count)
        density = density[: max_family_size + 1]
        all_extension_candidates = all_extension_candidates[: max_family_size + 1]
        all_extension_valid = all_extension_valid[: max_family_size + 1]

        aggregate_profiles = []
        for m, family_count in enumerate(density):
            canonical_growth = None
            canonical_i = None
            if m + 1 < len(density) and family_count > 0:
                canonical_growth = density[m + 1] / family_count
                canonical_i = safe_log2_ratio(canonical_growth, self.N - m)

            all_i = safe_log2_ratio(all_extension_valid[m], all_extension_candidates[m])
            aggregate_profiles.append(
                {
                    "m": m,
                    "family_count": family_count,
                    "canonical_growth_rate_D_next_over_D": round_finite(canonical_growth),
                    "aggregate_I_close_bits": round_finite(canonical_i),
                    "all_unused_candidate_count": all_extension_candidates[m],
                    "all_unused_valid_count": all_extension_valid[m],
                    "all_unused_I_close_bits": round_finite(all_i),
                }
            )

        jamming_m = None
        for row in aggregate_profiles:
            growth = row["canonical_growth_rate_D_next_over_D"]
            if growth is not None and growth < 1.0:
                jamming_m = row["m"]
                break

        per_core_by_size = []
        for m in range(max_family_size + 1):
            for s in range(self.w):
                acc = channel_stats[m][s]
                if acc.channel_count:
                    per_core_by_size.append(acc.to_record(self.n, self.w, m, s))

        representative_core_profiles = self._representative_core_profiles(
            per_core_totals,
            max_family_size,
        )

        return {
            "n": self.n,
            "w": self.w,
            "k": self.k,
            "N_sites": self.N,
            "core_count_by_size": {
                str(s): len(self.core_size_to_indices[s]) for s in range(self.w)
            },
            "petal_channels": [
                {
                    "core_size_s": s,
                    "petal_size": self.w - s,
                    "meaning": "candidate set is C union P with C fixed and P inside the complement",
                }
                for s in range(self.w)
            ],
            "total_sf_free": sum(density),
            "max_family_size": max_family_size,
            "density_of_states": density,
            "jamming_m_star": jamming_m,
            "elapsed_seconds": round_finite(elapsed, digits=3),
            "aggregate_profiles": aggregate_profiles,
            "per_core_by_size": per_core_by_size,
            "representative_core_profiles": representative_core_profiles,
            "summary": summarize_run(self.n, self.w, self.N, jamming_m, aggregate_profiles, per_core_by_size),
        }

    def _representative_core_profiles(
        self,
        per_core_totals: dict[tuple[int, int], ChannelAccumulator],
        max_family_size: int,
    ) -> list[dict[str, Any]]:
        records = []
        for s in range(self.w):
            indices = self.core_size_to_indices[s]
            if not indices:
                continue
            for core_idx in indices[: min(3, len(indices))]:
                profile = []
                for m in range(max_family_size + 1):
                    acc = per_core_totals[(m, core_idx)]
                    if not acc.channel_count:
                        continue
                    profile.append(
                        {
                            "m": m,
                            "candidate_count": acc.candidate_count,
                            "local_valid_count": acc.local_valid_count,
                            "I_core_local_weighted_bits": round_finite(
                                safe_log2_ratio(acc.local_valid_count, acc.candidate_count)
                            ),
                        }
                    )
                records.append(
                    {
                        "core_id": mask_tuple(self.cores[core_idx], self.n),
                        "core_mask": self.cores[core_idx],
                        "core_size_s": s,
                        "petal_size": self.w - s,
                        "profile": profile,
                    }
                )
        return records


def summarize_run(
    n: int,
    w: int,
    N: int,
    jamming_m: int | None,
    aggregate_profiles: list[dict[str, Any]],
    per_core_by_size: list[dict[str, Any]],
) -> dict[str, Any]:
    if jamming_m is None:
        target_m = max(row["m"] for row in aggregate_profiles)
    else:
        target_m = jamming_m

    rows_at_target = [row for row in per_core_by_size if row["m"] == target_m]
    strongest = max(
        rows_at_target,
        key=lambda row: row["I_core_local_weighted_bits"] or -1.0,
        default=None,
    )
    aggregate_at_target = next(
        (row for row in aggregate_profiles if row["m"] == target_m),
        None,
    )

    local_values = [
        row["I_core_local_weighted_bits"]
        for row in rows_at_target
        if row["I_core_local_weighted_bits"] is not None
    ]
    local_range = max(local_values) - min(local_values) if local_values else None

    return {
        "target_m": target_m,
        "target_is_jamming_m": jamming_m is not None and target_m == jamming_m,
        "aggregate_I_close_bits_at_target": (
            aggregate_at_target["aggregate_I_close_bits"] if aggregate_at_target else None
        ),
        "all_unused_I_close_bits_at_target": (
            aggregate_at_target["all_unused_I_close_bits"] if aggregate_at_target else None
        ),
        "strongest_core_size_at_target": (
            strongest["core_size_s"] if strongest else None
        ),
        "strongest_petal_size_at_target": strongest["petal_size"] if strongest else None,
        "strongest_I_core_local_weighted_bits_at_target": (
            strongest["I_core_local_weighted_bits"] if strongest else None
        ),
        "I_core_local_weighted_range_by_core_size_at_target": round_finite(local_range),
        "comparison_note": (
            "per-core local I isolates closure through a fixed core; aggregate I_close "
            "uses D(m+1)/D(m) and is included only as a continuity check"
        ),
        "ambient_state_count": N,
        "n": n,
        "w": w,
    }


def input_metadata(paths: list[Path], root: Path) -> list[dict[str, Any]]:
    records = []
    for path in paths:
        abs_path = path.resolve()
        try:
            relative = str(abs_path.relative_to(root))
        except ValueError:
            relative = str(abs_path)
        records.append(
            {
                "path": relative,
                "bytes": abs_path.stat().st_size,
                "sha256": sha256_file(abs_path),
            }
        )
    return records


def classify(runs: list[dict[str, Any]]) -> dict[str, Any]:
    executed = bool(runs)
    target_runs = [run for run in runs if run["summary"]["target_is_jamming_m"]]
    runs_with_nonzero_core_i = [
        run
        for run in target_runs
        if (
            run["summary"]["strongest_I_core_local_weighted_bits_at_target"] is not None
            and run["summary"]["strongest_I_core_local_weighted_bits_at_target"] > 0.25
        )
    ]
    runs_with_stratification = [
        run
        for run in target_runs
        if (
            run["summary"]["I_core_local_weighted_range_by_core_size_at_target"] is not None
            and run["summary"]["I_core_local_weighted_range_by_core_size_at_target"] > 0.1
        )
    ]

    if not executed:
        classification = "PER_CORE_BLOCKED"
        meaning = "No exact per-core run was executed."
    elif runs_with_nonzero_core_i and runs_with_stratification:
        classification = "PER_CORE_SIGNAL_PRESENT"
        meaning = (
            "Exact cheap regimes show a nonzero per-core closure-pressure channel, "
            "and the channel depends on core size near the aggregate jamming point."
        )
    else:
        classification = "PER_CORE_INCONCLUSIVE"
        meaning = (
            "Per-core instrumentation ran, but the observed core channel was too weak "
            "or too unstratified to separate from aggregate jamming."
        )

    return {
        "classification": classification,
        "allowed_classifications": [
            "PER_CORE_SIGNAL_PRESENT",
            "PER_CORE_INCONCLUSIVE",
            "PER_CORE_BLOCKED",
        ],
        "meaning": meaning,
        "executed_target_count": len(runs),
        "jamming_target_count": len(target_runs),
        "nonzero_core_signal_target_count": len(runs_with_nonzero_core_i),
        "stratified_target_count": len(runs_with_stratification),
        "leg4_status": (
            "precursor only; no floor-normalized ratio, no geometry-exhausting sweep, "
            "and no theorem-status change"
        ),
        "claim_ceiling": "shadow signature, not universal law",
        "A_axis": "A0",
    }


def compact_run_table(runs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for run in runs:
        summary = run["summary"]
        rows.append(
            {
                "w": run["w"],
                "n": run["n"],
                "N_sites": run["N_sites"],
                "total_sf_free": run["total_sf_free"],
                "max_family_size": run["max_family_size"],
                "jamming_m_star": run["jamming_m_star"],
                "aggregate_I_close_bits_at_target": summary["aggregate_I_close_bits_at_target"],
                "all_unused_I_close_bits_at_target": summary["all_unused_I_close_bits_at_target"],
                "strongest_core_size_at_target": summary["strongest_core_size_at_target"],
                "strongest_I_core_local_weighted_bits_at_target": summary[
                    "strongest_I_core_local_weighted_bits_at_target"
                ],
                "elapsed_seconds": run["elapsed_seconds"],
            }
        )
    return rows


def markdown_table(headers: list[str], rows: list[list[Any]]) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(format_cell(value) for value in row) + " |")
    return lines


def format_cell(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def selected_core_size_rows(runs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selected = []
    wanted = {(3, 7), (4, 6), (2, 10)}
    for run in runs:
        if (run["w"], run["n"]) not in wanted:
            continue
        target_m = run["summary"]["target_m"]
        for row in run["per_core_by_size"]:
            if row["m"] == target_m:
                selected.append(row)
    return selected


def render_report(results: dict[str, Any]) -> str:
    classification = results["classification"]["classification"]
    rows = results["compact_run_table"]
    core_rows = selected_core_size_rows(results["runs"])

    lines: list[str] = [
        f"# {EXPERIMENT_ID} - Report",
        "",
        f"**Date:** {RUN_DATE}",
        "**Problem:** Erdos #20 sunflower core closure",
        "**Scope:** Exact per-core closure instrumentation on cheap calibration regimes",
        f"**Classification:** {classification}",
        "**Claim ceiling:** shadow signature, not universal law",
        "",
        "## Meaning",
        "",
        (
            "The aggregate Hessian run could see jamming, but it could not say where "
            "the closure pressure lived. This run adds that missing channel for the "
            "cheap exact regimes: a core is an actual subset C, and a petal extension "
            "is a candidate w-set of the form C union P."
        ),
        "",
        (
            "The per-core signal is present in this limited sense: near the aggregate "
            "jamming point, fixed-core channels become locally expensive, and the "
            "expense depends on core size. That makes the earlier aggregate trace more "
            "specific. It is still only a precursor, because the expensive regimes and "
            "the floor-normalized comparison are not executed here."
        ),
        "",
        "## Operational Definitions",
        "",
        "- `core_id`: the sorted 1-based elements of a subset C of [n].",
        "- `petal channel`: fixed C with candidate petals P in the complement of C, so the candidate set is C union P.",
        "- `I_core_local(C,s,m)`: `-log2(local valid petal extensions / unused petal extensions through C)` at family size m.",
        "- `I_core_global(C,s,m)`: the same denominator, but the numerator requires the candidate extension to be globally sunflower-free.",
        "- `aggregate I_close`: the prior continuity check using `D(m+1)/D(m)` divided by unused ambient sites.",
        "",
        "## Executed Targets",
        "",
    ]

    lines.extend(
        markdown_table(
            [
                "w",
                "n",
                "N",
                "families",
                "M",
                "m*",
                "aggregate I",
                "all-unused I",
                "strongest s",
                "strongest I_core",
                "seconds",
            ],
            [
                [
                    row["w"],
                    row["n"],
                    row["N_sites"],
                    row["total_sf_free"],
                    row["max_family_size"],
                    row["jamming_m_star"],
                    row["aggregate_I_close_bits_at_target"],
                    row["all_unused_I_close_bits_at_target"],
                    row["strongest_core_size_at_target"],
                    row["strongest_I_core_local_weighted_bits_at_target"],
                    row["elapsed_seconds"],
                ]
                for row in rows
            ],
        )
    )

    lines.extend(
        [
            "",
            "## Core-Size Slice At Target m",
            "",
        ]
    )
    lines.extend(
        markdown_table(
            [
                "w",
                "n",
                "m",
                "s",
                "petal",
                "channels",
                "local I",
                "global I",
                "median finite I",
                "blocked frac",
            ],
            [
                [
                    row["w"],
                    row["n"],
                    row["m"],
                    row["core_size_s"],
                    row["petal_size"],
                    row["channel_count"],
                    row["I_core_local_weighted_bits"],
                    row["I_core_global_weighted_bits"],
                    row["median_finite_I_core_bits"],
                    row["blocked_channel_fraction"],
                ]
                for row in core_rows
            ],
        )
    )

    lines.extend(
        [
            "",
            "## Classification",
            "",
            f"The classifier returns **{classification}**.",
            "",
            (
                "This means the per-core channel is no longer blocked by missing "
                "instrumentation for the cheap regimes. It does not mean a Leg-4 pass. "
                "A pass would require the floor-normalized numerator, geometry controls, "
                "and a larger symmetry-reduced sweep."
            ),
            "",
            "The comparison to aggregate `I_close` should be read narrowly. Aggregate `I_close` asks how the density of states changes with family size. `I_core` asks how a fixed core channel closes as petals accumulate. Agreement in the late regime is a useful precursor; mismatch is expected because the denominators are different.",
            "",
            "## Blockers",
            "",
        ]
    )
    lines.extend(
        [
            "- `w=3,n=8` and `w=4,n=7` are too large for this exact Python per-core pass.",
            "- `w=4,n=8` was already beyond the prior exhaustive aggregate run.",
            "- No Mendoza-floor or construction-normalized numerator is predeclared here.",
            "- No Abbott-Hansen-Sauer baseline family is instrumented as a control.",
            "- The runner is exact enumeration, not the transfer-matrix rewrite needed for scale.",
        ]
    )

    lines.extend(
        [
            "",
            "## Claim Limits",
            "",
            (
                "A-axis remains A0. The artifact measures an internal diagnostic "
                "channel. It does not change the status of Erdos #20, does not give "
                "new lower-bound progress, and does not justify public power-morphism "
                "language."
            ),
            "",
            "## Artifacts",
            "",
            f"- `{EXPERIMENT_ID}_RESULTS.json`",
            f"- `{EXPERIMENT_ID}_REPORT.md`",
            f"- `{EXPERIMENT_ID}_RESULTS.sha256`",
            "- `sunflower_per_core_closure_analysis.py`",
            "",
        ]
    )

    return "\n".join(lines)


def parse_targets(raw_targets: list[str] | None) -> list[tuple[int, int]]:
    if not raw_targets:
        return DEFAULT_TARGETS
    parsed = []
    for raw in raw_targets:
        try:
            w_raw, n_raw = raw.split(",", 1)
            parsed.append((int(w_raw), int(n_raw)))
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"target must look like w,n; got {raw!r}"
            ) from exc
    return parsed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target",
        action="append",
        help="Run one target as w,n. May be repeated. Defaults to the cheap calibration set.",
    )
    args = parser.parse_args()

    base = Path(__file__).resolve().parent
    math_root = base.parent.parent
    targets = parse_targets(args.target)

    input_paths = [
        base / "EXP-MATH-ERDOS20-HESSIAN-CLOSURE-20260505-01_RESULTS.json",
        base / "EXP-MATH-ERDOS20-HESSIAN-CLOSURE-20260505-01_REPORT.md",
        base / "EXP-MATH-ERDOS20-SUNFLOWER-001_RESULTS.json",
        base / "EXP-MATH-ERDOS20-SUNFLOWER-002_RESULTS.json",
        base / "SUNFLOWER_CORE_CLOSURE_LEG4_PACKET_2026-05-05.md",
        base / "sunflower_transfer_matrix.cpp",
    ]
    for path in input_paths:
        if not path.exists():
            raise FileNotFoundError(path)
        if path.suffix == ".json":
            read_json(path)
        elif path.suffix in {".md", ".cpp"}:
            path.read_text(encoding="utf-8")

    runs = []
    for w, n in targets:
        enumerator = PerCoreEnumerator(n=n, w=w, k=3)
        runs.append(enumerator.run())

    results = {
        "experiment_id": EXPERIMENT_ID,
        "title": "Sunflower per-core closure instrumentation",
        "run_date": RUN_DATE,
        "generated_on_host_date": date.today().isoformat(),
        "problem": "Erdos #20 sunflower conjecture",
        "scope": "exact per-core closure instrumentation on cheap calibration regimes",
        "input_metadata": input_metadata(input_paths, math_root),
        "executed_targets": [{"w": w, "n": n, "k": 3} for w, n in targets],
        "expensive_targets_not_executed": KNOWN_EXPENSIVE_TARGETS,
        "analysis_definitions": {
            "core_id": "sorted 1-based elements of C",
            "petal_channel": "candidate set equals fixed core C union a petal P in the complement",
            "I_core_local_bits": "-log2(local valid petal extensions through C / unused petal extensions through C)",
            "I_core_global_bits": "-log2(globally safe petal extensions through C / unused petal extensions through C)",
            "aggregate_I_close_bits": "-log2((D(m+1)/D(m))/(N-m)); continuity check against the prior aggregate run",
        },
        "classification": classify(runs),
        "compact_run_table": compact_run_table(runs),
        "runs": runs,
        "claim_controls": {
            "A_axis": "A0",
            "safe_framing": "per-core diagnostic; shadow signature, not universal law",
            "forbidden_framing": [
                "new lower-bound contribution",
                "status change for the conjecture",
                "public power-morphism evidence",
                "universal law",
            ],
        },
    }

    results_path = base / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = base / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = base / f"{EXPERIMENT_ID}_RESULTS.sha256"

    results_path.write_text(
        json.dumps(results, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    report_path.write_text(render_report(results), encoding="utf-8")
    sha_path.write_text(f"{sha256_file(results_path)}  {results_path.name}\n", encoding="utf-8")

    print(f"Wrote {results_path.name}")
    print(f"Wrote {report_path.name}")
    print(f"Wrote {sha_path.name}")
    print(results["classification"]["classification"])


if __name__ == "__main__":
    main()
