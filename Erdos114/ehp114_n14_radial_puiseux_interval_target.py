#!/usr/bin/env python3
"""Build the n=14 radial Puiseux interval-hardening target packet.

This is the first concrete "middle kingdom" closure target after the Tao
constant-chase packet. It deliberately stays below proof language: it converts
the DOI-backed n=14 interval certificate plus the radial hypergeometric
calibration into a narrow theorem target for interval arithmetic.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import mpmath as mp


EXPERIMENT_ID = "EXP-MATH-EHP114-N14-RADIAL-PUISEUX-INTERVAL-TARGET-20260505-01"
CALIBRATION_ID = "EXP-MATH-EHP114-N14-RADIAL-HYPERGEOMETRIC-CALIBRATION-20260502-01"
INARI_ID = "EXP-MM-EHP-007-n14-inari"
TAO_PACKET_ID = "EXP-MATH-EHP114-TAO-MIDDLE-KINGDOM-20260505-01"
ZENODO_DOI = "10.5281/zenodo.19480329"


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


def radial_length_inside(n: int, a: mp.mpf) -> mp.mpf:
    """Exact radial length for p_a(z)=z^n-a, |a|<1."""
    p = mp.mpf(n - 1) / (2 * n)
    return 2 * mp.pi * mp.hyper([p, p], [1], a * a)


def radial_length_at_one(n: int) -> mp.mpf:
    """Exact radial length at a=1 from the gamma formula."""
    return (
        (mp.mpf(2) ** (mp.mpf(1) / n + 1))
        * mp.sqrt(mp.pi)
        * mp.gamma(mp.mpf(1) / (2 * n))
        / (2 * mp.gamma(mp.mpf(n + 1) / (2 * n)))
    )


def direct_radial_rows(n: int, eps_values: list[str]) -> list[dict[str, str]]:
    l_star = radial_length_at_one(n)
    exponent = mp.mpf(1) / n
    rows = []
    for eps_text in eps_values:
        eps = mp.mpf(eps_text)
        a = 1 - eps
        length = radial_length_inside(n, a)
        deficit = l_star - length
        rows.append(
            {
                "eps": eps_text,
                "a": mp.nstr(a, 40),
                "exact_length": mp.nstr(length, 60),
                "deficit_vs_exact_L_star": mp.nstr(deficit, 60),
                "deficit_over_eps_power_1_over_n": mp.nstr(deficit / (eps**exponent), 60),
            }
        )
    return rows


def min_ratio(rows: list[dict[str, Any]], tail_eps: mp.mpf | None = None) -> dict[str, Any]:
    candidates = []
    for row in rows:
        eps = mp.mpf(str(row["eps"]))
        if tail_eps is not None and eps > tail_eps:
            continue
        ratio = mp.mpf(str(row["deficit_over_eps_power_1_over_n"]))
        candidates.append((ratio, eps, row))
    if not candidates:
        raise ValueError("no rows in ratio window")
    ratio, eps, row = min(candidates, key=lambda item: item[0])
    return {
        "eps": mp.nstr(eps, 20),
        "ratio": mp.nstr(ratio, 30),
        "row": row,
    }


def calibration_admissible_rows(calibration: dict[str, Any]) -> list[dict[str, Any]]:
    rows = calibration.get("radial_rows") or calibration.get("asymptotic", {}).get("rows") or []
    return [row for row in rows if row.get("admissible_under_roots_in_unit_disk", True)]


def build_result() -> dict[str, Any]:
    mp.mp.dps = 90
    root = repo_root_from_script()
    out_dir = Path(__file__).resolve().parent
    results_dir = root / "erdos-experiments" / "results" / "erdos-114"
    calibration_path = results_dir / f"{CALIBRATION_ID}_RESULTS.json"
    inari_path = results_dir / f"{INARI_ID}_RESULTS.json"
    tao_path = out_dir / f"{TAO_PACKET_ID}_RESULTS.json"

    calibration = load_json(calibration_path)
    inari = load_json(inari_path)
    tao = load_json(tao_path)

    eps_values = ["0.1", "0.05", "0.02", "0.01", "0.001", "0.0001", "0.00001", "0.000001", "0.00000001"]
    direct_rows = direct_radial_rows(14, eps_values)
    source_rows = calibration_admissible_rows(calibration)

    direct_global = min_ratio(direct_rows)
    direct_tail = min_ratio(direct_rows, tail_eps=mp.mpf("0.0001"))
    source_global = min_ratio(source_rows)
    source_tail = min_ratio(source_rows, tail_eps=mp.mpf("0.0001"))

    exact_l_star = radial_length_at_one(14)
    inari_interval = {
        "l_star_lower": inari["l_star_lower"],
        "l_star_upper": inari["l_star_upper"],
        "exact_inside_interval": mp.mpf(str(inari["l_star_lower"])) <= exact_l_star <= mp.mpf(str(inari["l_star_upper"])),
    }

    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "scope": "fixed-n=14 radial Puiseux interval-hardening target",
        "packet_verdict": "CERTIFICATE_TARGET_READY",
        "claim_ceiling": "This is a shadow signature, not universal law. It is a certificate target, not a proof.",
        "sources": {
            "doi": ZENODO_DOI,
            "calibration": {
                "experiment_id": CALIBRATION_ID,
                "path": str(calibration_path),
                "sha256": sha256_file(calibration_path),
                "status": calibration.get("status"),
                "claim_ceiling": calibration.get("claim_ceiling"),
            },
            "inari_ieee1788": {
                "experiment_id": INARI_ID,
                "path": str(inari_path),
                "sha256": sha256_file(inari_path),
                "verdict": inari.get("verdict"),
                "rigor": inari.get("rigor"),
                "reduced_dim": inari.get("reduced_dim"),
                "bb_total_evals": inari.get("bb_total_evals"),
                "bb_proof_complete": inari.get("bb_proof_complete"),
                "outer_domain_safe": inari.get("outer_domain_safe"),
                "hessian_negative": inari.get("hessian_negative"),
                "total_time_secs": inari.get("total_time_secs"),
                "l_star_interval": inari_interval,
            },
            "tao_middle_kingdom_packet": {
                "experiment_id": TAO_PACKET_ID,
                "path": str(tao_path),
                "sha256": sha256_file(tao_path),
                "first_theorem_target": tao.get("first_theorem_target"),
            },
        },
        "radial_model": {
            "family": "p_a(z)=z^14-a",
            "exact_length_inside": "L_14(a)=2*pi*2F1(13/28,13/28;1;a^2) for |a|<1",
            "exact_length_at_boundary": "L_14(1)=2^(1/14+1)*sqrt(pi)*Gamma(1/28)/(2*Gamma(15/28))",
            "singularity_scale": "L_14(1)-L_14(1-eps) is Puiseux scale eps^(1/14), not quadratic Hessian scale.",
            "exact_L_star": mp.nstr(exact_l_star, 80),
        },
        "diagnostic_constants": {
            "direct_a_equals_1_minus_eps": {
                "rows": direct_rows,
                "global_min_ratio": direct_global,
                "tail_min_ratio_eps_le_1e_minus_4": direct_tail,
                "working_constant_C14": 24.0,
                "tail_constant_C14_tail": 26.0,
                "reason": (
                    "C14=24 sits below every sampled direct-radial ratio, including "
                    "the largest perturbation eps=0.1. C14_tail=26 is a sharper "
                    "tail-only target for eps <= 1e-4."
                ),
            },
            "calibration_radius_schedule": {
                "source_schedule": "radius = 1 - eps/sqrt(14), a = radius^14",
                "global_min_ratio": source_global,
                "tail_min_ratio_eps_le_1e_minus_4": source_tail,
                "working_constant_C14": 24.0,
                "tail_constant_C14_tail": 27.0,
                "reason": (
                    "The existing calibration schedule has larger sampled ratios "
                    "than direct a=1-eps; C14=24 stays deliberately conservative."
                ),
            },
        },
        "interval_hardening_strategy": [
            {
                "stage": "singular_tail",
                "domain": "0 < eps <= 1e-4",
                "target": "Use the Gauss 2F1 connection formula at z=1 to prove a lower Puiseux coefficient above 26 for n=14.",
                "claim_ceiling": "Fixed-n radial-family interval lemma only.",
            },
            {
                "stage": "compact_middle",
                "domain": "1e-4 <= eps <= 1e-1",
                "target": "Use interval subdivision on the exact hypergeometric/integral formula to certify the weaker C14=24 bound.",
                "claim_ceiling": "Numerical interval hardening of a one-parameter slice.",
            },
            {
                "stage": "radial_to_shape_splice",
                "domain": "n=14 local cone around z^14-1",
                "target": "Feed the radial deficit into the shape-cone and remainder packets; do not infer global EHP from radial data alone.",
                "claim_ceiling": "Bridge target between DOI finite proof and Tao constant chase.",
            },
        ],
        "lean_shaped_target": """theorem ehp114_n14_radial_puiseux_interval_direct
    (eps : Real) (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10000) :
    (24 : Real) * Real.rpow eps ((1 : Real) / 14)
      <= D14 (radialMode14 eps) := by
  -- fixed-n Gauss 2F1 connection formula + interval constants
  -- no global EHP conclusion follows from this lemma alone
  sorry""",
        "single_current_bottleneck": "Interval-hardening the radial hypergeometric/Puiseux singularity model.",
        "next_theorem_target_for_agent": "Prove ehp114_n14_radial_puiseux_interval_direct with an interval-certified 2F1 connection-formula lower bound.",
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
    direct = result["diagnostic_constants"]["direct_a_equals_1_minus_eps"]
    schedule = result["diagnostic_constants"]["calibration_radius_schedule"]
    inari = result["sources"]["inari_ieee1788"]
    return f"""# EHP114 n=14 Radial Puiseux Interval Target

Experiment: `{EXPERIMENT_ID}`

## Meaning

This packet turns the n=14 stress case into a narrow proof target. The point is
not to rerun the expensive branch-and-bound certificate. The point is to carve
out the first analytic bridge between the DOI-backed finite proof and Tao's
large-n proof: a fixed-degree radial Puiseux lower bound that can be interval
hardened.

The claim ceiling is unchanged: this is a shadow signature, not universal law.
It is not a proof of #114, not a global local-stability theorem, and not a
scorecard upgrade.

## Source Anchor

- DOI: `{ZENODO_DOI}`
- Finite certificate: `{inari["experiment_id"]}`
- Verdict: `{inari["verdict"]}`
- Rigor: `{inari["rigor"]}`
- Reduced dimension: `{inari["reduced_dim"]}`
- Interval evaluations: `{inari["bb_total_evals"]}`
- B&B complete: `{inari["bb_proof_complete"]}`
- Exact L* inside inari interval: `{inari["l_star_interval"]["exact_inside_interval"]}`

n=14 is the materially stringent finite anchor because it is the largest
byte-reconciled DOI-backed Rust/inari certificate in the current corpus. The
n=15/n=16 records remain useful hints, but their zero-evaluation provenance
needs reconciliation before they can carry the same weight.

## Candidate Constants

| radial parameterization | sampled lower evidence | working constant |
|---|---:|---:|
| direct `a=1-eps` | min ratio {direct["global_min_ratio"]["ratio"]} at eps {direct["global_min_ratio"]["eps"]} | C14 = {direct["working_constant_C14"]} |
| direct tail `eps<=1e-4` | min ratio {direct["tail_min_ratio_eps_le_1e_minus_4"]["ratio"]} at eps {direct["tail_min_ratio_eps_le_1e_minus_4"]["eps"]} | C14_tail = {direct["tail_constant_C14_tail"]} |
| calibration schedule `radius=1-eps/sqrt(14)` | min ratio {schedule["global_min_ratio"]["ratio"]} at eps {schedule["global_min_ratio"]["eps"]} | C14 = {schedule["working_constant_C14"]} |
| calibration tail `eps<=1e-4` | min ratio {schedule["tail_min_ratio_eps_le_1e_minus_4"]["ratio"]} at eps {schedule["tail_min_ratio_eps_le_1e_minus_4"]["eps"]} | C14_tail = {schedule["tail_constant_C14_tail"]} |

These are floating diagnostics. The usable next move is to replace them with
interval-certified constants. The deliberate target constant is conservative:
prove 24 first, then sharpen only if that proof is clean.

## Interval-Hardening Split

1. Singular tail: for `0 < eps <= 1e-4`, use the Gauss 2F1 connection formula
   at `z=1` to certify the Puiseux coefficient.
2. Compact middle: for `1e-4 <= eps <= 1e-1`, interval-subdivide the exact
   hypergeometric or integral formula and certify the weaker `C14=24` bound.
3. Radial-to-shape splice: feed this radial deficit into the shape cone and
   mixed-remainder packets. Do not infer the global theorem from the radial
   slice alone.

## Lean-Shaped Target

```lean
{result["lean_shaped_target"]}
```

## Single Bottleneck

{result["single_current_bottleneck"]}

The next agent should attack exactly this theorem target:
`{result["next_theorem_target_for_agent"]}`.

## Guardrails

No scorecard, D1, public document, Lean file, git, or email state was changed.
"""


def main() -> int:
    out_dir = Path(__file__).resolve().parent
    result_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"

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
                "result": str(result_path),
                "report": str(report_path),
                "sha256": sha_path.read_text(encoding="utf-8").split()[0],
                "packet_verdict": result["packet_verdict"],
                "next_theorem_target": result["next_theorem_target_for_agent"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
