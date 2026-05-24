#!/usr/bin/env python3
"""Exact radial calculations for EHP114.

This is a diagnostic calculation, not a proof of Erdos #114.

For the radial family

    p_a(z) = z^n - a

the lemniscate |p_a(z)| = 1 can be parameterized by z^n = a + exp(it).
Summing over the n branches gives the exact length

    L_n(a) = integral_0^(2 pi) |a + exp(it)|^(1/n - 1) dt.

For |a| < 1 this is the constant term of

    (1 + a exp(it))^beta (1 + a exp(-it))^beta,

where beta = (1/n - 1)/2. Hence

    L_n(a) = 2 pi * 2F1(p, p; 1; a^2), p = (n - 1)/(2n).

For a > 1, factor out a^(2 beta). At a = 1, use Gauss' gamma-formula
equivalent. This avoids the quadrature failure near the singular level set.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import mpmath as mp


EXPERIMENT_ID = "EXP-MATH-EHP114-RADIAL-HYPERGEOMETRIC-20260502-01"
DEFAULT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS = DEFAULT_ROOT / "erdos-experiments/results/erdos-114"
DEFAULT_SHORTCUT_RESULTS = DEFAULT_ROOT / "erdos-experiments/scripts/erdos-114"
DEFAULT_FOURIER_ARTIFACT = (
    DEFAULT_RESULTS / "EXP-MATH-EHP114-N15-FOURIER-HESSIAN-20260502-02_RESULTS.json"
)
DEFAULT_EPS = "0.02,0.01,0.005,0.0025,0.001,0.0005,0.0001"
DEFAULT_ASYMPTOTIC_EPS = (
    "0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,"
    "0.0002,0.0001,0.00005,0.00002,0.00001,0.000005,0.000001,0.0000001,0.00000001"
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def decimal(x: mp.mpf, digits: int = 50) -> str:
    return mp.nstr(x, digits)


def parse_mp_list(raw: str) -> list[mp.mpf]:
    values = [mp.mpf(part.strip()) for part in raw.split(",") if part.strip()]
    if not values:
        raise ValueError("empty numeric list")
    return values


def radial_length_limit(n: int) -> mp.mpf:
    alpha = mp.mpf(1) / n - 1
    return (
        mp.power(2, alpha + 2)
        * mp.sqrt(mp.pi)
        * mp.gamma((alpha + 1) / 2)
        / (2 * mp.gamma((alpha + 2) / 2))
    )


def radial_length_hypergeometric(n: int, a: mp.mpf) -> mp.mpf:
    beta = mp.mpf(1) / (2 * n) - mp.mpf(1) / 2
    p = -beta
    a = mp.mpf(a)
    if mp.almosteq(a, 1):
        return radial_length_limit(n)
    if a < 1:
        return 2 * mp.pi * mp.hyper([p, p], [1], a * a)
    return mp.power(a, 2 * beta) * 2 * mp.pi * mp.hyper([p, p], [1], 1 / (a * a))


def reference_artifact(root: Path, shortcut_root: Path, n: int) -> dict[str, Any] | None:
    candidates = [
        root / f"EXP-MM-EHP-007-n{n}-inari_RESULTS.json",
        shortcut_root / f"EXP-MM-EHP-007-n{n}-inari_RESULTS.json",
    ]
    for path in candidates:
        data = load_json(path)
        if data is not None:
            return {
                "path": str(path),
                "sha256": sha256_file(path),
                "data": data,
            }
    return None


def artifact_interval_compare(exact: mp.mpf, artifact: dict[str, Any] | None) -> dict[str, Any]:
    if artifact is None:
        return {
            "artifact_path": None,
            "artifact_sha256": None,
            "l_star_lower": None,
            "l_star_upper": None,
            "exact_inside_artifact_interval": None,
            "difference_to_midpoint": None,
        }
    data = artifact["data"]
    lower = data.get("l_star_lower")
    upper = data.get("l_star_upper")
    if lower is None or upper is None:
        return {
            "artifact_path": artifact["path"],
            "artifact_sha256": artifact["sha256"],
            "l_star_lower": lower,
            "l_star_upper": upper,
            "exact_inside_artifact_interval": None,
            "difference_to_midpoint": None,
            "artifact_verdict": data.get("verdict"),
        }
    lower_mp = mp.mpf(str(lower))
    upper_mp = mp.mpf(str(upper))
    midpoint = (lower_mp + upper_mp) / 2
    return {
        "artifact_path": artifact["path"],
        "artifact_sha256": artifact["sha256"],
        "l_star_lower": str(lower),
        "l_star_upper": str(upper),
        "exact_inside_artifact_interval": bool(lower_mp <= exact <= upper_mp),
        "difference_to_midpoint": decimal(exact - midpoint, 30),
        "artifact_verdict": data.get("verdict"),
    }


def linear_slope(xs: list[float], ys: list[float]) -> float:
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    denom = sum((x - x_mean) ** 2 for x in xs)
    return sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys)) / denom


def load_marching_radial_rows(path: Path) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    artifact = load_json(path)
    if artifact is None:
        return None, []
    rows = artifact.get("singularity_sanity", {}).get("rows", [])
    return artifact, rows


def find_marching_row(rows: list[dict[str, Any]], eps: mp.mpf, sign: int) -> dict[str, Any] | None:
    eps_float = float(eps)
    for row in rows:
        row_sign = -1 if float(row["radius"]) < 1 else 1
        if row_sign == sign and abs(float(row["eps"]) - eps_float) <= 1e-14:
            return row
    return None


def build_result(
    n: int,
    experiment_id: str,
    eps_values: list[mp.mpf],
    asymptotic_eps_values: list[mp.mpf],
    result_root: Path,
    shortcut_root: Path,
    fourier_artifact_path: Path,
) -> dict[str, Any]:
    l_star = radial_length_limit(n)
    fourier_artifact, marching_rows = load_marching_radial_rows(fourier_artifact_path)
    reference = reference_artifact(result_root, shortcut_root, n)

    reference_rows = []
    for ref_n in [14, 15, 16]:
        exact = radial_length_limit(ref_n)
        reference_rows.append(
            {
                "n": ref_n,
                "exact_radial_length_at_a_equals_1": decimal(exact, 50),
                "artifact_compare": artifact_interval_compare(
                    exact,
                    reference_artifact(result_root, shortcut_root, ref_n),
                ),
            }
        )

    radial_rows = []
    for eps in eps_values:
        for sign in [-1, 1]:
            radius = mp.mpf(1) + sign * eps / mp.sqrt(n)
            a = mp.power(radius, n)
            exact_length = radial_length_hypergeometric(n, a)
            marching = find_marching_row(marching_rows, eps, sign)
            marching_length = mp.mpf(str(marching["length"])) if marching else None
            radial_rows.append(
                {
                    "eps": decimal(eps, 20),
                    "sign": sign,
                    "radius": decimal(radius, 30),
                    "a": decimal(a, 30),
                    "admissible_under_roots_in_unit_disk": bool(radius <= 1),
                    "exact_length": decimal(exact_length, 50),
                    "deficit_vs_exact_L_star": decimal(l_star - exact_length, 50),
                    "deficit_over_eps_power_1_over_n": decimal(
                        (l_star - exact_length) / mp.power(eps, mp.mpf(1) / n),
                        40,
                    ),
                    "marching_length_from_fourier_artifact": (
                        decimal(marching_length, 30) if marching_length is not None else None
                    ),
                    "exact_minus_marching": (
                        decimal(exact_length - marching_length, 30)
                        if marching_length is not None
                        else None
                    ),
                }
            )

    asymptotic_rows = []
    for eps in asymptotic_eps_values:
        radius = mp.mpf(1) - eps / mp.sqrt(n)
        a = mp.power(radius, n)
        exact_length = radial_length_hypergeometric(n, a)
        deficit = l_star - exact_length
        asymptotic_rows.append(
            {
                "eps": decimal(eps, 20),
                "radius": decimal(radius, 30),
                "a": decimal(a, 30),
                "deficit_vs_exact_L_star": decimal(deficit, 50),
                "deficit_over_eps_power_1_over_n": decimal(
                    deficit / mp.power(eps, mp.mpf(1) / n),
                    40,
                ),
            }
        )

    slope_rows = []
    for window in [4, 5, 6, 8, 10, 12, 16]:
        if len(asymptotic_rows) < window:
            continue
        tail = asymptotic_rows[-window:]
        xs = [math.log(float(row["eps"])) for row in tail]
        ys = [math.log(float(row["deficit_vs_exact_L_star"])) for row in tail]
        slope_rows.append({"tail_window": window, "fitted_log_log_slope": linear_slope(xs, ys)})

    l0_marching = None
    if fourier_artifact is not None:
        l0_marching = fourier_artifact.get("reference", {}).get("L0_marching_squares")

    return {
        "experiment_id": experiment_id,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "degree": n,
        "status": "ANALYTIC_NUMERICAL_DIAGNOSTIC_NOT_PROOF",
        "claim_ceiling": (
            "Exact radial-family calculation. It explains the singular boundary layer "
            "near z^n - 1, but it does not prove global or local EHP114 maximality."
        ),
        "formula": {
            "family": "p_a(z) = z^n - a",
            "length_integral": "L_n(a) = integral_0^(2 pi) |a + exp(i t)|^(1/n - 1) dt",
            "hypergeometric_for_abs_a_lt_1": "L_n(a) = 2 pi * 2F1(p,p;1;a^2), p=(n-1)/(2n)",
            "hypergeometric_for_a_gt_1": "L_n(a) = a^(1/n - 1) * 2 pi * 2F1(p,p;1;a^(-2))",
            "a_equals_1_gamma_formula": (
                "L_n(1) = 2^(1/n + 1) * sqrt(pi) * Gamma(1/(2n)) "
                "/ (2 * Gamma((n+1)/(2n)))"
            ),
            "boundary_layer_prediction": (
                "L_n(1) - L_n(a) scales like |1-a|^(1/n), not like a quadratic Hessian."
            ),
        },
        "reference": {
            "exact_L_star_for_degree": decimal(l_star, 50),
            "artifact_compare": artifact_interval_compare(l_star, reference),
            "fourier_probe_artifact": str(fourier_artifact_path),
            "fourier_probe_artifact_sha256": (
                sha256_file(fourier_artifact_path) if fourier_artifact_path.exists() else None
            ),
            "fourier_probe_L0_marching_squares": l0_marching,
            "fourier_probe_L0_error_vs_exact_L_star": (
                decimal(mp.mpf(str(l0_marching)) - l_star, 40) if l0_marching is not None else None
            ),
        },
        "reference_rows": reference_rows,
        "radial_rows": radial_rows,
        "asymptotic": {
            "admissible_side": "radius < 1, roots remain inside the unit disk",
            "expected_slope": decimal(mp.mpf(1) / n, 30),
            "slope_rows": slope_rows,
            "rows": asymptotic_rows,
        },
        "interpretation": {
            "primary": (
                "The radial family around z^n - 1 is a singular boundary-layer problem. "
                "The deficit follows an eps^(1/n) scale, so a smooth finite Hessian at "
                "the singular point is the wrong local model."
            ),
            "what_survives_from_fourier_probe": (
                "The large positive finite-difference signal survives as a stratified "
                "deficit/unfolding signal. It should not be described as an ordinary "
                "smooth Hessian certificate."
            ),
            "next_mathematical_move": (
                "Replace the finite-Hessian framing with a Puiseux or hypergeometric "
                "singularity certificate for radial contractions, then build interval "
                "Fourier/tensor bounds for nonradial admissible root perturbations."
            ),
        },
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    ref = result["reference"]
    lines = [
        f"# {result['experiment_id']} Report",
        "",
        "## Status",
        "",
        result["claim_ceiling"],
        "",
        "This is a diagnostic calculation, not a proof of Erdos #114.",
        "",
        "## Exact Formula",
        "",
        "For `p_a(z) = z^n - a`, parameterize the lemniscate by `z^n = a + exp(i t)`.",
        "Summing the derivative length over all branches gives:",
        "",
        "`L_n(a) = integral_0^(2 pi) |a + exp(i t)|^(1/n - 1) dt`.",
        "",
        "For `|a| < 1`, the constant-term expansion gives:",
        "",
        "`L_n(a) = 2 pi * 2F1(p,p;1;a^2)`, where `p = (n - 1)/(2n)`.",
        "",
        "For `a = 1`, the Gauss/gamma limit gives the reference value directly.",
        "",
        "## Reference Check",
        "",
        f"- Degree: {result['degree']}",
        f"- Exact `L_n(1)`: `{ref['exact_L_star_for_degree']}`",
        f"- Artifact interval contains exact value: {ref['artifact_compare']['exact_inside_artifact_interval']}",
        f"- Fourier-probe marching `L0`: {ref['fourier_probe_L0_marching_squares']}",
        f"- Fourier-probe marching `L0 - exact L*`: `{ref['fourier_probe_L0_error_vs_exact_L_star']}`",
        "",
        "| n | exact L_n(1) | artifact contains exact? | artifact path |",
        "|---:|---:|---|---|",
    ]
    for row in result["reference_rows"]:
        cmp_row = row["artifact_compare"]
        lines.append(
            "| {n} | `{exact}` | {contains} | `{path}` |".format(
                n=row["n"],
                exact=row["exact_radial_length_at_a_equals_1"],
                contains=cmp_row["exact_inside_artifact_interval"],
                path=cmp_row["artifact_path"],
            )
        )

    lines.extend(
        [
            "",
            "## Radial Perturbation Table",
            "",
            "The inward rows are admissible under the roots-in-unit-disk constraint. The outward rows are diagnostic only.",
            "",
            "| eps | side | radius | admissible? | exact length | deficit vs exact L* | marching length | exact - marching |",
            "|---:|---|---:|---|---:|---:|---:|---:|",
        ]
    )
    for row in result["radial_rows"]:
        side = "in" if row["sign"] == -1 else "out"
        lines.append(
            "| {eps} | {side} | `{radius}` | {ok} | `{length}` | `{deficit}` | {marching} | {diff} |".format(
                eps=row["eps"],
                side=side,
                radius=row["radius"],
                ok=row["admissible_under_roots_in_unit_disk"],
                length=row["exact_length"],
                deficit=row["deficit_vs_exact_L_star"],
                marching=(
                    f"`{row['marching_length_from_fourier_artifact']}`"
                    if row["marching_length_from_fourier_artifact"] is not None
                    else ""
                ),
                diff=f"`{row['exact_minus_marching']}`" if row["exact_minus_marching"] is not None else "",
            )
        )

    lines.extend(
        [
            "",
            "## Boundary-Layer Scaling",
            "",
            "For admissible inward contractions, the deficit behaves like `eps^(1/n)`, not `eps^2`.",
            f"For n = {result['degree']}, the expected exponent is `{result['asymptotic']['expected_slope']}`.",
            "",
            "| tail window | fitted log-log slope |",
            "|---:|---:|",
        ]
    )
    for row in result["asymptotic"]["slope_rows"]:
        lines.append(f"| {row['tail_window']} | {row['fitted_log_log_slope']:.12f} |")

    lines.extend(
        [
            "",
            "Tail rows used for the scaling check:",
            "",
            "| eps | deficit vs exact L* | deficit / eps^(1/n) |",
            "|---:|---:|---:|",
        ]
    )
    for row in result["asymptotic"]["rows"][-8:]:
        lines.append(
            f"| {row['eps']} | `{row['deficit_vs_exact_L_star']}` | `{row['deficit_over_eps_power_1_over_n']}` |"
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            result["interpretation"]["primary"],
            "",
            result["interpretation"]["what_survives_from_fourier_probe"],
            "",
            result["interpretation"]["next_mathematical_move"],
            "",
            "## Source Boundary",
            "",
            f"- Fourier probe artifact: `{ref['fourier_probe_artifact']}`",
            f"- Fourier probe SHA-256: `{ref['fourier_probe_artifact_sha256']}`",
            f"- Reference artifact: `{ref['artifact_compare']['artifact_path']}`",
            f"- Reference artifact SHA-256: `{ref['artifact_compare']['artifact_sha256']}`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=15)
    parser.add_argument("--eps", default=DEFAULT_EPS)
    parser.add_argument("--asymptotic-eps", default=DEFAULT_ASYMPTOTIC_EPS)
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--shortcut-dir", type=Path, default=DEFAULT_SHORTCUT_RESULTS)
    parser.add_argument("--fourier-artifact", type=Path, default=DEFAULT_FOURIER_ARTIFACT)
    parser.add_argument("--experiment-id", default=EXPERIMENT_ID)
    args = parser.parse_args()

    mp.mp.dps = args.mp_dps
    eps_values = parse_mp_list(args.eps)
    asymptotic_eps_values = parse_mp_list(args.asymptotic_eps)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    result_path = args.out_dir / f"{args.experiment_id}_RESULTS.json"
    report_path = args.out_dir / f"{args.experiment_id}_REPORT.md"
    sha_path = args.out_dir / f"{args.experiment_id}_RESULTS.sha256"
    for path in [result_path, report_path, sha_path]:
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {path}")

    result = build_result(
        n=args.n,
        experiment_id=args.experiment_id,
        eps_values=eps_values,
        asymptotic_eps_values=asymptotic_eps_values,
        result_root=args.out_dir,
        shortcut_root=args.shortcut_dir,
        fourier_artifact_path=args.fourier_artifact,
    )
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(result, report_path)
    sha_path.write_text(f"{sha256_file(result_path)}  {result_path.name}\n", encoding="utf-8")
    print(json.dumps({"result": str(result_path), "report": str(report_path), "sha256": sha_path.read_text().split()[0]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
