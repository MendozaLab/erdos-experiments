#!/usr/bin/env python3
"""n=14 shape-cone / remainder splice packet for EHP114.

This consumes the newly certified radial Puiseux lane and runs the existing
floating tensor-cone scaffold at n=14. It does not certify the shape cone. It
decides whether the next local-stability theorem target is well-posed and
records the exact mixed-remainder obligation.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


sys.dont_write_bytecode = True

EXPERIMENT_ID = "EXP-MATH-EHP114-N14-SHAPE-CONE-SPLICE-20260505-01"
RADIAL_COMPACT_ID = "EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01"
RADIAL_TAIL_ID = "EXP-MATH-EHP114-N14-RADIAL-TAIL-ARB-20260505-01"
RADIAL_C14 = 24.0
SHAPE_LAMBDA14 = 100000.0


def repo_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


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
    return {
        "method": "Gershgorin lower bound on floating mixed-shape Hessian proxy",
        "row_count": len(rows),
        "lower_bound_float": min(row["diagonal_minus_radius"] for row in rows),
        "all_rows_positive": all(row["diagonal_minus_radius"] > 0 for row in rows),
        "worst_row": min(rows, key=lambda row: row["diagonal_minus_radius"]),
        "rows": rows,
    }


def run_tensor_probe(root: Path) -> tuple[dict[str, Any], float]:
    module = load_module(
        "ehp114_tensor_cone_scaffold",
        root / "erdos-experiments" / "scripts" / "probe_ehp114_tensor_cone_scaffold.py",
    )
    started = time.monotonic()
    result = module.run_probe(
        n=14,
        eps_values=[0.04, 0.02, 0.01, 0.005],
        res=300,
        mixed_res=220,
        extent=3.0,
        result_dir=root / "erdos-experiments" / "results" / "erdos-114",
        include_mixed_matrix=True,
    )
    return result, time.monotonic() - started


def build_result() -> dict[str, Any]:
    root = repo_root_from_script()
    out_dir = Path(__file__).resolve().parent
    compact_path = out_dir / f"{RADIAL_COMPACT_ID}_RESULTS.json"
    tail_path = out_dir / f"{RADIAL_TAIL_ID}_RESULTS.json"
    compact = load_json(compact_path)
    tail = load_json(tail_path)
    tensor, elapsed = run_tensor_probe(root)

    mixed = tensor["mixed_shape_hessian_proxy"]
    eigvals = [float(value) for value in mixed["eigenvalues"]]
    matrix = [[float(value) for value in row] for row in mixed["matrix"]]
    gersh = gershgorin_lower_bound(matrix)
    min_eig = min(eigvals)
    shape_lambda_ok = 0.0 < SHAPE_LAMBDA14 < min(min_eig, gersh["lower_bound_float"])

    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "scope": "n=14 shape-cone/remainder splice target",
        "claim_ceiling": "This is a shadow signature, not universal law. Shape-cone values are diagnostic, not interval-certified.",
        "status": "SHAPE_SPLICE_TARGET_READY" if shape_lambda_ok else "SHAPE_SPLICE_NEEDS_REPAIR",
        "radial_certificates": {
            "compact": {
                "experiment_id": RADIAL_COMPACT_ID,
                "path": str(compact_path),
                "sha256": sha256_file(compact_path),
                "status": compact["status"],
                "min_margin": compact["min_margin"],
            },
            "tail": {
                "experiment_id": RADIAL_TAIL_ID,
                "path": str(tail_path),
                "sha256": sha256_file(tail_path),
                "status": tail["status"],
                "all_bins_pass": tail["all_bins_pass"],
            },
            "working_C14": RADIAL_C14,
            "spendable_radial_budget": RADIAL_C14 / 2.0,
        },
        "tensor_probe": {
            "runtime_seconds": elapsed,
            "degree": tensor["degree"],
            "parameters": tensor["parameters"],
            "summary": tensor["summary"],
            "mixed_min_eigenvalue_float": min_eig,
            "mixed_max_eigenvalue_float": max(eigvals),
            "mixed_positive_eigenvalues": sum(1 for value in eigvals if value > 0.0),
            "mixed_eigenvalue_count": len(eigvals),
            "gershgorin": gersh,
            "working_lambda14": SHAPE_LAMBDA14,
            "working_lambda14_below_diagnostic_bounds": shape_lambda_ok,
        },
        "mixed_remainder_obligation": {
            "status": "OPEN",
            "target": (
                "For 0 < eps <= 1e-1 and ||s|| <= eta14(eps), prove "
                "R14(eps,s) <= 12*eps^(1/14) + 0.5*lambda14*||s||^2."
            ),
            "why_it_matters": (
                "The radial lane now supplies 24*eps^(1/14). The splice target "
                "spends half that radial budget to absorb mixed radial/shape terms."
            ),
        },
        "lean_shaped_target": """theorem ehp114_n14_radial_shape_remainder_splice
    (eps : Real) (s : ShapeQuotient 14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hs : quotientNorm s <= eta14 eps) :
    (12 : Real) * Real.rpow eps ((1 : Real) / 14)
      + (50000 : Real) * quotientNormSq s
      <= D14 (radialMode14 eps + shapeMode14 s) := by
  -- radial certificates + interval shape cone + mixed-remainder absorption
  sorry""",
        "next_blocker": "Replace the floating n=14 tensor-cone matrix with an interval matrix certificate and prove mixed-remainder absorption.",
        "guardrails": {
            "no_scorecard_update": True,
            "no_d1_update": True,
            "no_public_claim": True,
            "no_lean_file_created": True,
            "no_email": True,
            "no_git": True,
        },
    }


def build_report(result: dict[str, Any]) -> str:
    tensor = result["tensor_probe"]
    gersh = tensor["gershgorin"]
    summary = tensor["summary"]
    shape_total = (
        summary["shape_hessian_like_count_slope_ge_1_5"]
        + summary["shape_singular_like_count_slope_lt_1_5"]
    )
    shape_positive = shape_total if summary["all_directions_positive_symmetric_deficit"] else "not-all"
    return f"""# EHP114 n=14 Shape-Cone Splice Packet

Experiment: `{EXPERIMENT_ID}`

## Meaning

The radial/Puiseux lane is now certified for the radial family. This packet
asks the next question: does the n=14 nonradial shape cone look strong enough
to formulate a local-stability splice theorem?

The answer is yes as a target, not as a proof. This is a shadow signature, not
universal law. The shape numbers below are floating diagnostics; the next step
is interval hardening.

## Verdict

- Status: `{result["status"]}`
- Radial compact status: `{result["radial_certificates"]["compact"]["status"]}`
- Radial tail status: `{result["radial_certificates"]["tail"]["status"]}`
- Shape lambda target: `{tensor["working_lambda14"]}`
- Shape lambda below diagnostic bounds: `{tensor["working_lambda14_below_diagnostic_bounds"]}`
- Current blocker: `{result["next_blocker"]}`

## Shape-Cone Diagnostic

| quantity | value |
|---|---:|
| quotient basis rank | {tensor["parameters"]["quotient_basis_rank"]} |
| expected rank `2n-3` | {tensor["parameters"]["expected_rank_2n_minus_3"]} |
| shape positive symmetric deficits | {shape_positive}/{shape_total} |
| shape slopes below Hessian threshold 1.5 | {summary["shape_singular_like_count_slope_lt_1_5"]}/{shape_total} |
| mixed positive eigenvalues | {tensor["mixed_positive_eigenvalues"]}/{tensor["mixed_eigenvalue_count"]} |
| floating min eigenvalue | {tensor["mixed_min_eigenvalue_float"]:.12g} |
| Gershgorin lower bound | {gersh["lower_bound_float"]:.12g} |
| worst Gershgorin row | {gersh["worst_row"]["index"]} |

## Mixed Remainder Obligation

{result["mixed_remainder_obligation"]["target"]}

This is now the live mathematical bottleneck. The radial singularity is no
longer the first obstruction; the obstruction is proving that mixed terms
cannot eat more than half the certified radial deficit.

## Lean-Shaped Target

```lean
{result["lean_shaped_target"]}
```

## Source Boundary

No scorecard, D1, public document, Lean file, git, or email state was changed.
"""


def main() -> int:
    out_dir = Path(__file__).resolve().parent
    result_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"

    if result_path.exists() and (not report_path.exists() or not sha_path.exists()):
        result = load_json(result_path)
    else:
        for path in (result_path, report_path, sha_path):
            if path.exists():
                raise SystemExit(f"Refusing to overwrite existing artifact: {path}")
        result = build_result()
        result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_path.write_text(build_report(result), encoding="utf-8")
    sha_path.write_text(f"{sha256_file(result_path)}  {result_path.name}\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "experiment_id": EXPERIMENT_ID,
                "status": result["status"],
                "report": str(report_path),
                "sha256": sha_path.read_text(encoding="utf-8").split()[0],
                "next_blocker": result["next_blocker"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
