#!/usr/bin/env python3
"""Branch-count census for the EHP114 n=14 hard-cell slab artifact.

This is a diagnostic post-processor. It does not certify length and does not
prove Erdős #114. It answers the next cheap question: do the current slabs
look like one owned branch per x-slab, or are unresolved/multi-branch tubes
still present before interval Newton/Bernstein work?
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-N14-BRANCH-COUNT-CENSUS-HARD-CELL-20260506-01"
SOURCE_EXPERIMENT_ID = "EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03"
BASE = Path(__file__).resolve().parent
SOURCE_DIR = BASE / SOURCE_EXPERIMENT_ID
SOURCE_RESULTS = SOURCE_DIR / f"{SOURCE_EXPERIMENT_ID}_RESULTS.json"
OUT_DIR = BASE / EXPERIMENT_ID
COARSE_BUCKET_COUNT = 100


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bucket_for_ix(ix: int, x_slab_count: int) -> int:
    return min(COARSE_BUCKET_COUNT - 1, (ix * COARSE_BUCKET_COUNT) // x_slab_count)


def compact_interval(interval: list[float] | None) -> list[float] | None:
    if interval is None:
        return None
    return [round(float(interval[0]), 12), round(float(interval[1]), 12)]


def summarize_examples(records: list[dict[str, Any]], limit: int = 20) -> list[dict[str, Any]]:
    out = []
    for row in records[:limit]:
        out.append(
            {
                "ix": row["ix"],
                "accepted_count": row["accepted_count"],
                "unresolved_count": row["unresolved_count"],
                "total_count": row["total_count"],
                "x_interval": compact_interval(row.get("x_interval")),
                "accepted_y_intervals": [
                    compact_interval(item) for item in row.get("accepted_y_intervals", [])[:5]
                ],
                "unresolved_y_intervals": [
                    compact_interval(item) for item in row.get("unresolved_y_intervals", [])[:5]
                ],
                "unresolved_reasons": row.get("unresolved_reasons", {}),
            }
        )
    return out


def build_report(result: dict[str, Any]) -> str:
    dist = result["fine_slab_total_branch_count_distribution"]
    return f"""# {EXPERIMENT_ID} Report

## Verdict

- Status: `{result['status']}`
- Source: `{SOURCE_EXPERIMENT_ID}`
- Fine x-slabs: `{result['fine_x_slab_count']}`
- Active fine x-slabs: `{result['active_fine_slab_count']}`
- Active slabs with exactly one accepted branch and no unresolved tube: `{result['active_single_certified_slab_count']}`
- Fine slabs with unresolved tubes: `{result['fine_slabs_with_unresolved_count']}`
- Fine slabs with multiple branch groups: `{result['fine_slabs_with_multiple_groups_count']}`
- Coarse buckets: `{COARSE_BUCKET_COUNT}`

## Distribution

Fine-slab total branch/tube count distribution:

```json
{json.dumps(dist, indent=2, sort_keys=True)}
```

## Meaning

The z32 slab artifact already puts the certified accepted branch length below
the hard-cell cap, but this census shows the single-branch assumption is not
yet certified across the domain. The next proof-facing step is therefore not
another global length pass. It is an interval Newton or Bernstein isolation
resolver over the unresolved tubes named by this census.

## Claim Ceiling

This is a local branch-count diagnostic. It is not a proof of Erdős #114 and
not a global n=14 proof.
"""


def main() -> int:
    if not SOURCE_RESULTS.exists():
        raise SystemExit(f"Missing source artifact: {SOURCE_RESULTS}")
    OUT_DIR.mkdir(parents=True, exist_ok=False)
    result_path = OUT_DIR / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = OUT_DIR / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = OUT_DIR / f"{EXPERIMENT_ID}_RESULTS.sha256"

    data = json.loads(SOURCE_RESULTS.read_text(encoding="utf-8"))
    params = data["parameters"]
    x_slab_count = int(params["x_slab_count"])

    by_ix: dict[int, dict[str, Any]] = {}
    reason_counts: Counter[str] = Counter()

    def row_for_ix(ix: int) -> dict[str, Any]:
        if ix not in by_ix:
            by_ix[ix] = {
                "ix": ix,
                "accepted_count": 0,
                "unresolved_count": 0,
                "accepted_y_intervals": [],
                "unresolved_y_intervals": [],
                "unresolved_reasons": Counter(),
                "x_interval": None,
            }
        return by_ix[ix]

    for branch in data.get("branches", []):
        ix = int(branch["ix"])
        row = row_for_ix(ix)
        row["accepted_count"] += 1
        row["accepted_y_intervals"].append(branch["y_interval"])
        row["x_interval"] = row["x_interval"] or branch["x_interval"]

    for tube in data.get("unresolved_branches", []):
        ix = int(tube["ix"])
        reason = tube.get("reason", "unknown")
        row = row_for_ix(ix)
        row["unresolved_count"] += 1
        row["unresolved_y_intervals"].append(tube["y_interval"])
        row["unresolved_reasons"][reason] += 1
        row["x_interval"] = row["x_interval"] or tube["x_interval"]
        reason_counts[reason] += 1

    active_rows = []
    for ix, row in by_ix.items():
        row["total_count"] = row["accepted_count"] + row["unresolved_count"]
        row["unresolved_reasons"] = dict(row["unresolved_reasons"])
        active_rows.append(row)
    active_rows.sort(key=lambda item: item["ix"])

    distribution = Counter(str(row["total_count"]) for row in active_rows)
    single_certified = [
        row for row in active_rows if row["accepted_count"] == 1 and row["unresolved_count"] == 0
    ]
    unresolved_rows = [row for row in active_rows if row["unresolved_count"] > 0]
    multi_rows = [row for row in active_rows if row["total_count"] > 1]
    zero_accept_unresolved = [
        row for row in active_rows if row["accepted_count"] == 0 and row["unresolved_count"] > 0
    ]

    buckets: dict[int, dict[str, Any]] = {
        idx: {
            "bucket": idx,
            "fine_slab_start": math.floor(idx * x_slab_count / COARSE_BUCKET_COUNT),
            "fine_slab_end_exclusive": math.floor((idx + 1) * x_slab_count / COARSE_BUCKET_COUNT),
            "active_fine_slab_count": 0,
            "max_total_count": 0,
            "unresolved_fine_slab_count": 0,
            "multi_group_fine_slab_count": 0,
        }
        for idx in range(COARSE_BUCKET_COUNT)
    }
    for row in active_rows:
        bucket = buckets[bucket_for_ix(row["ix"], x_slab_count)]
        bucket["active_fine_slab_count"] += 1
        bucket["max_total_count"] = max(bucket["max_total_count"], row["total_count"])
        if row["unresolved_count"] > 0:
            bucket["unresolved_fine_slab_count"] += 1
        if row["total_count"] > 1:
            bucket["multi_group_fine_slab_count"] += 1

    buckets_with_risk = [
        bucket
        for bucket in buckets.values()
        if bucket["unresolved_fine_slab_count"] or bucket["multi_group_fine_slab_count"]
    ]

    status = (
        "BRANCH_CENSUS_PASS_SINGLE_CERTIFIED_BRANCH"
        if len(single_certified) == len(active_rows)
        else "BRANCH_CENSUS_FAIL_AMBIGUOUS_ROOT_ISOLATION"
    )
    result: dict[str, Any] = {
        "experiment_id": EXPERIMENT_ID,
        "source_experiment_id": SOURCE_EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "source_status": data["status"],
        "source_total_validated_length_upper": data["total_validated_length_upper"],
        "source_exact_length_cap": data["length_budget"]["exact_length_cap"],
        "source_margin_to_cap": data["margin_to_cap"],
        "fine_x_slab_count": x_slab_count,
        "active_fine_slab_count": len(active_rows),
        "active_single_certified_slab_count": len(single_certified),
        "fine_slabs_with_unresolved_count": len(unresolved_rows),
        "fine_slabs_with_multiple_groups_count": len(multi_rows),
        "fine_slabs_with_unresolved_only_count": len(zero_accept_unresolved),
        "fine_slab_total_branch_count_distribution": dict(sorted(distribution.items(), key=lambda kv: int(kv[0]))),
        "unresolved_reason_counts": dict(reason_counts),
        "coarse_bucket_count": COARSE_BUCKET_COUNT,
        "coarse_buckets_with_risk_count": len(buckets_with_risk),
        "coarse_buckets_with_risk": buckets_with_risk,
        "example_unresolved_slabs": summarize_examples(unresolved_rows),
        "example_multi_group_slabs": summarize_examples(multi_rows),
        "claim_ceiling": "Local branch-count diagnostic only. Not a proof of Erdős #114 and not a global n=14 proof.",
        "next_blocker": "Run interval Newton/Krawczyk or Bernstein isolation on unresolved tubes; do not rerun global length until unresolved tubes are resolved or excluded.",
    }
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_path.write_text(build_report(result), encoding="utf-8")
    sha_path.write_text(f"{sha256_file(result_path)}  {result_path.name}\n", encoding="utf-8")
    print(json.dumps({
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "active_fine_slab_count": len(active_rows),
        "single_certified": len(single_certified),
        "unresolved_slabs": len(unresolved_rows),
        "multi_group_slabs": len(multi_rows),
        "result": str(result_path),
        "sha256": sha256_file(result_path),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
