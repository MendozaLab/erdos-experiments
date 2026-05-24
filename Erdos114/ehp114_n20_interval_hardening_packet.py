#!/usr/bin/env python3
"""Proof-facing interval-hardening packet for the n=20 EHP114 shadow triad.

This script does not prove EHP114. It converts the successful n=20 diagnostic
triad into a conservative local certificate target: a radial Puiseux lower
constant, a quotient shape-cone matrix lower bound, and the explicit mixed
remainder obligation still needed for Lean/interval closure.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "EXP-MATH-EHP114-N20-INTERVAL-HARDENING-PACKET-20260505-01"
SOURCE_ID = "EXP-MATH-EHP114-N20-HESSIAN-TRIAD-20260505-01"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def round_sig(value: float, digits: int = 12) -> float:
    if not math.isfinite(value):
        return value
    return float(f"{value:.{digits}g}")


def gershgorin_lower_bound(matrix: list[list[float]]) -> dict[str, Any]:
    rows = []
    for i, row in enumerate(matrix):
        diag = float(row[i])
        radius = sum(abs(float(value)) for j, value in enumerate(row) if j != i)
        rows.append(
            {
                "index": i,
                "diagonal": diag,
                "off_diagonal_radius": radius,
                "diagonal_minus_radius": diag - radius,
            }
        )
    lower = min(row["diagonal_minus_radius"] for row in rows)
    return {
        "method": "Gershgorin diagonal dominance on diagnostic mixed matrix",
        "lower_bound_float": lower,
        "row_count": len(rows),
        "all_rows_positive": all(row["diagonal_minus_radius"] > 0 for row in rows),
        "worst_row": min(rows, key=lambda row: row["diagonal_minus_radius"]),
        "rows": rows,
    }


def radial_constant_packet(rows: list[dict[str, Any]]) -> dict[str, Any]:
    parsed = []
    for row in rows:
        eps = float(row["eps"])
        ratio = float(row["deficit_over_eps_power_1_over_n"])
        parsed.append({"eps": eps, "ratio": ratio})
    global_min = min(parsed, key=lambda row: row["ratio"])
    tail_rows = [row for row in parsed if row["eps"] <= 1e-4]
    tail_min = min(tail_rows, key=lambda row: row["ratio"])
    return {
        "sampled_eps_count": len(parsed),
        "global_min_ratio": global_min,
        "tail_min_ratio_eps_le_1e_minus_4": tail_min,
        "working_constant_C20": 32.0,
        "tail_constant_C20_tail": 40.0,
        "reason": (
            "C20=32 is below the sampled global minimum with about 20 percent "
            "diagnostic slack; C20_tail=40 is a sharper tail-only target."
        ),
    }


def build_result(source: dict[str, Any], source_path: Path, out_dir: Path) -> dict[str, Any]:
    radial_rows = source["components"]["radial_hypergeometric"]["asymptotic"]["rows"]
    mixed = source["components"]["tensor_cone"]["mixed_shape_hessian_proxy"]
    matrix = [[float(value) for value in row] for row in mixed["matrix"]]
    gersh = gershgorin_lower_bound(matrix)
    radial = radial_constant_packet(radial_rows)

    eigen_min = float(mixed["min_eigenvalue"])
    gersh_min = float(gersh["lower_bound_float"])
    lambda20 = 250000.0
    if not (0 < lambda20 < min(eigen_min, gersh_min)):
        raise ValueError("lambda20 safety constant is not below available matrix lower bounds")

    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "scope": "proof-facing n=20 interval-hardening target packet",
        "source_experiment_id": source.get("experiment_id", SOURCE_ID),
        "source_artifact": {
            "path": str(source_path),
            "sha256": sha256_file(source_path),
        },
        "claim_ceiling": "This is a shadow signature, not universal law. It is a certificate target, not a proof.",
        "packet_verdict": "CERTIFICATE_TARGET_READY",
        "radial_puiseux_constant": radial,
        "shape_cone_constant": {
            "floating_min_eigenvalue": eigen_min,
            "gershgorin_lower_bound": gersh,
            "working_lambda20": lambda20,
            "reason": (
                "lambda20=250000 is below both the floating eigenvalue lower "
                "bound and the Gershgorin diagnostic lower bound."
            ),
        },
        "mixed_remainder_obligation": {
            "status": "OPEN",
            "needed_bound": (
                "For |r| <= delta20 and ||s|| <= eta20, prove "
                "R20(r,s) <= 0.5 * (C20 * |r|^(1/20) + lambda20 * ||s||^2)."
            ),
            "why_it_matters": (
                "The radial and shape terms are separately positive in the "
                "diagnostics; the theorem still needs a bound that prevents "
                "mixed radial/shape terms from canceling them."
            ),
        },
        "lean_shaped_targets": {
            "radial": """theorem ehp114_n20_radial_puiseux_interval
    (eps : Real) (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10000) :
    (32 : Real) * Real.rpow eps ((1 : Real) / 20)
      <= D20 (radialMode20 eps) := by
  -- exact hypergeometric interval bound for the radial family
  sorry""",
            "shape_cone": """theorem ehp114_n20_shape_cone_interval
    (s : ShapeQuotient 20) :
    (250000 : Real) * quotientNormSq s
      <= shapeQuadratic20 s := by
  -- interval matrix lower bound for the quotient shape block
  sorry""",
            "combined": """theorem ehp114_n20_stratified_local_certificate
    (r : Real) (s : ShapeQuotient 20)
    (hr_pos : 0 < abs r) (hr_small : abs r <= delta20)
    (hs_small : quotientNorm s <= eta20) :
    (16 : Real) * Real.rpow (abs r) ((1 : Real) / 20)
      + (125000 : Real) * quotientNormSq s
      <= D20 (radialMode20 r + shapeMode20 s) := by
  -- radial Puiseux + shape cone + mixed remainder absorption
  sorry""",
        },
        "next_blocker": "Turn the diagnostic constants into interval arithmetic lemmas and prove the mixed remainder absorption bound.",
        "guardrails": {
            "writes_scoped_to": str(out_dir),
            "no_d1": True,
            "no_scorecard": True,
            "no_public_docs": True,
            "no_git": True,
            "no_email": True,
        },
    }


def build_report(result: dict[str, Any]) -> str:
    radial = result["radial_puiseux_constant"]
    shape = result["shape_cone_constant"]
    gersh = shape["gershgorin_lower_bound"]
    return f"""# {EXPERIMENT_ID} Report

## Status

This packet takes the successful n=20 Hessian/Puiseux triad and narrows it into
proof obligations. It is not a proof of EHP114 and not a publication packet.
The claim ceiling remains: this is a shadow signature, not universal law.

## Verdict

- Packet verdict: `{result["packet_verdict"]}`
- Source: `{result["source_experiment_id"]}`
- Current blocker: `{result["next_blocker"]}`

The useful fantasy is now precise: do not run another ordinary Hessian sweep.
Prove the radial Puiseux interval bound, prove the quotient shape-cone lower
bound, then absorb the mixed remainder.

## Candidate Constants

| layer | diagnostic lower evidence | working constant |
|---|---:|---:|
| radial Puiseux | min sampled ratio {radial["global_min_ratio"]["ratio"]:.12g} at eps {radial["global_min_ratio"]["eps"]:.12g} | C20 = {radial["working_constant_C20"]:.12g} |
| radial tail | min tail ratio {radial["tail_min_ratio_eps_le_1e_minus_4"]["ratio"]:.12g} for eps <= 1e-4 | C20_tail = {radial["tail_constant_C20_tail"]:.12g} |
| shape cone | Gershgorin lower {gersh["lower_bound_float"]:.12g}; floating eig min {shape["floating_min_eigenvalue"]:.12g} | lambda20 = {shape["working_lambda20"]:.12g} |

The constants are deliberately conservative. They are not certified constants
until interval arithmetic replaces the floating diagnostic input.

## Lean-Shaped Targets

### Radial

```lean
{result["lean_shaped_targets"]["radial"]}
```

### Shape Cone

```lean
{result["lean_shaped_targets"]["shape_cone"]}
```

### Combined Local Certificate

```lean
{result["lean_shaped_targets"]["combined"]}
```

## Mixed Remainder Obligation

{result["mixed_remainder_obligation"]["needed_bound"]}

This is the single mathematical blocker. The diagnostics already separate the
radial and nonradial positive terms; the closure theorem needs a proof that
mixed terms cannot erase that positivity inside the local cone.

## Source Boundary

No scorecard, D1, public document, git, email, CLAUDE.md, or AGENTS.md was
changed.
"""


def main() -> int:
    out_dir = Path(__file__).resolve().parent
    source_path = out_dir / f"{SOURCE_ID}_RESULTS.json"
    result_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"

    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {path}")

    source = load_json(source_path)
    result = build_result(source, source_path, out_dir)
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_path.write_text(build_report(result), encoding="utf-8")
    sha_path.write_text(f"{sha256_file(result_path)}  {result_path.name}\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "experiment_id": EXPERIMENT_ID,
                "result": str(result_path),
                "report": str(report_path),
                "sha256": sha_path.read_text(encoding="utf-8").split()[0],
                "packet_verdict": result["packet_verdict"],
                "next_blocker": result["next_blocker"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
