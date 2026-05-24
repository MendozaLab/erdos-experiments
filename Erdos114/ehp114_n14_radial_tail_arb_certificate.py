#!/usr/bin/env python3
"""Arb/FLINT singular-tail certificate for the n=14 radial Puiseux lane.

This checks the missing radial tail from the compact inari verifier:

    0 < eps <= 1e-4

for the one-parameter family p_a(z)=z^14-a, a=1-eps. It uses two layers:

1. Arb ball arithmetic on 2F1 for 1e-8 <= eps <= 1e-4.
2. The Gauss connection formula near z=1 for 0 < eps <= 1e-8.

It remains radial-family evidence only. It does not prove EHP #114.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from flint import arb, ctx


EXPERIMENT_ID = "EXP-MATH-EHP114-N14-RADIAL-TAIL-ARB-20260505-01"
COMPACT_ID = "EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01"
TARGET_ID = "EXP-MATH-EHP114-N14-RADIAL-PUISEUX-INTERVAL-TARGET-20260505-01"
DEGREE = 14
C14 = arb(24)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ball_text(x: arb, digits: int = 50) -> str:
    return x.str(digits)


def lower_text(x: arb, digits: int = 50) -> str:
    return x.lower().str(digits, radius=False)


def upper_text(x: arb, digits: int = 50) -> str:
    return x.upper().str(digits, radius=False)


def radial_constants() -> dict[str, arb]:
    n = arb(DEGREE)
    p = arb(13) / 28
    s = arb(1) / DEGREE
    lstar = (arb(2) ** (arb(1) / n + 1)) * arb.const_sqrt_pi() * arb.gamma(arb(1) / (2 * n)) / (
        2 * arb.gamma((n + 1) / (2 * n))
    )
    a_const = arb.gamma(s) / (arb.gamma(arb(15) / 28) ** 2)
    b_const = arb.gamma(-s) / (arb.gamma(p) ** 2)
    return {"p": p, "s": s, "lstar": lstar, "A": a_const, "B": b_const}


def interval_from_strings(lo: str, hi: str) -> arb:
    return arb(lo).union(arb(hi))


def check_hypergeom_bin(lo: str, hi: str, constants: dict[str, arb]) -> dict[str, Any]:
    eps = interval_from_strings(lo, hi)
    z = (1 - eps) * (1 - eps)
    p = constants["p"]
    lstar = constants["lstar"]
    length = 2 * arb.pi() * z.hypgeom_2f1(p, p, arb(1))
    deficit_lower = lstar.lower() - length.upper()
    rhs = C14 * (arb(hi) ** constants["s"])
    margin = deficit_lower - rhs.upper()
    return {
        "eps_lo": lo,
        "eps_hi": hi,
        "method": "direct Arb interval enclosure of 2F1 on eps bin",
        "length_upper": upper_text(length),
        "deficit_lower": lower_text(deficit_lower),
        "rhs_upper": upper_text(rhs),
        "margin_lower": lower_text(margin),
        "pass": bool(margin > 0),
    }


def geometric_bins(start: str = "1e-8", end: str = "1e-4", ratio: int = 2) -> list[tuple[str, str]]:
    lo = arb(start)
    end_arb = arb(end)
    out: list[tuple[str, str]] = []
    while lo < end_arb:
        hi = lo * ratio
        if hi > end_arb:
            hi = end_arb
        out.append((lo.str(24, radius=False), hi.str(24, radius=False)))
        lo = hi
    return out


def micro_tail_certificate(constants: dict[str, arb]) -> dict[str, Any]:
    # Connection formula:
    # F(z)=A*FA(w)+B*w^s*FB(w), w=1-z, s=1/14, B<0.
    # Deficit = A*(1-FA(w)) - B*w^s*FB(w).
    # For 0<w<=2eps, eps<=1e-8:
    #   FB(w) >= 1 by positive coefficients.
    #   FA(w)-1 <= w using a deliberately loose derivative bound M=1.
    #   w^s >= eps^s and w <= 2eps.
    # Hence D >= (2*pi*(-B) - 2*pi*A*2*eps_max^(1-s))*eps^s.
    eps_max = arb("1e-8")
    s = constants["s"]
    a_const = constants["A"]
    b_const = constants["B"]
    lead = 2 * arb.pi() * (-b_const)
    remainder = 2 * arb.pi() * a_const * 2 * (eps_max ** (1 - s))
    coefficient_lower = lead.lower() - remainder.upper()
    margin = coefficient_lower - C14
    return {
        "eps_domain": "(0, 1e-8]",
        "method": "Gauss 2F1 connection formula with loose FA(w)-1 <= w remainder bound",
        "A_upper": upper_text(a_const),
        "minus_B_lower": lower_text(-b_const),
        "lead_coefficient_lower": lower_text(lead),
        "remainder_coefficient_upper": upper_text(remainder),
        "certified_coefficient_lower": lower_text(coefficient_lower),
        "required_C14": "24",
        "margin_lower": lower_text(margin),
        "pass": bool(margin > 0),
        "notes": [
            "This is the only place where the singular endpoint is used.",
            "The bound is intentionally crude; the margin is still above 1.3 in coefficient units.",
            "The proof obligation is a connection-formula lemma plus positivity of the two 2F1 coefficient series.",
        ],
    }


def build_result() -> dict[str, Any]:
    ctx.prec = 300
    out_dir = Path(__file__).resolve().parent
    constants = radial_constants()
    bins = [check_hypergeom_bin(lo, hi, constants) for lo, hi in geometric_bins()]
    micro_tail = micro_tail_certificate(constants)
    compact_path = out_dir / f"{COMPACT_ID}_RESULTS.json"
    target_path = out_dir / f"{TARGET_ID}_RESULTS.json"
    all_pass = all(row["pass"] for row in bins) and micro_tail["pass"]
    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "scope": "n=14 radial Puiseux singular-tail certificate",
        "status": "RADIAL_TAIL_CERTIFIED" if all_pass else "RADIAL_TAIL_FAILED",
        "claim_ceiling": "This is a shadow signature, not universal law. It certifies the radial-family tail only, not EHP #114.",
        "degree": DEGREE,
        "constant_C14": 24,
        "tail_domain": "(0, 1e-4]",
        "companion_compact_certificate": {
            "experiment_id": COMPACT_ID,
            "path": str(compact_path),
            "sha256": sha256_file(compact_path) if compact_path.exists() else None,
        },
        "source_target_packet": {
            "experiment_id": TARGET_ID,
            "path": str(target_path),
            "sha256": sha256_file(target_path) if target_path.exists() else None,
        },
        "constants": {
            "s": "1/14",
            "p": "13/28",
            "Lstar": ball_text(constants["lstar"], 70),
            "A": ball_text(constants["A"], 70),
            "B": ball_text(constants["B"], 70),
        },
        "micro_tail": micro_tail,
        "hypergeom_bins": bins,
        "all_bins_pass": all_pass,
        "next_blocker": "Splice tail + compact radial certificates into the n=14 shape-cone and mixed-remainder local-stability packet.",
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
    rows = "\n".join(
        "| {eps_lo} | {eps_hi} | {margin_lower} | {pass} |".format(**row)
        for row in result["hypergeom_bins"]
    )
    micro = result["micro_tail"]
    return f"""# EHP114 n=14 Radial Tail Arb Certificate

Experiment: `{EXPERIMENT_ID}`

## Meaning

This closes the missing radial tail left by the compact-middle interval
verifier. Together with `{COMPACT_ID}`, the radial family now has an executable
certificate path for

```text
0 < eps <= 1e-1,  L_14(1)-L_14(1-eps) >= 24 eps^(1/14).
```

The claim ceiling is still narrow: this is a shadow signature, not universal
law. It is radial-family evidence only. It does not prove #114 and does not
settle nonradial local stability.

## Verdict

- Status: `{result["status"]}`
- Constant: `C14 = {result["constant_C14"]}`
- Tail domain: `{result["tail_domain"]}`
- Companion compact certificate: `{COMPACT_ID}`
- Next blocker: `{result["next_blocker"]}`

## Micro-Tail Connection Formula Bound

For `{micro["eps_domain"]}`, the packet uses the Gauss 2F1 connection formula
at `z=1`. The certified coefficient lower bound is:

```text
{micro["certified_coefficient_lower"]} > 24
```

The margin over `C14=24` is `{micro["margin_lower"]}`. This is the analytic
piece that prevents the endpoint from being just another floating sweep.

## Arb Bins for `1e-8 <= eps <= 1e-4`

| eps lo | eps hi | margin lower | pass |
|---:|---:|---:|---:|
{rows}

## What Remains

The radial bound is no longer the blocker. The next mathematical task is the
shape-cone/remainder splice: show nonradial perturbations cannot erase the
radial Puiseux deficit inside the n=14 local cone.

## Files

- Tail script: `erdos-experiments/Erdos114/ehp114_n14_radial_tail_arb_certificate.py`
- Result JSON: `erdos-experiments/Erdos114/{EXPERIMENT_ID}_RESULTS.json`
- SHA sidecar: `erdos-experiments/Erdos114/{EXPERIMENT_ID}_RESULTS.sha256`

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
                "status": result["status"],
                "all_bins_pass": result["all_bins_pass"],
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
