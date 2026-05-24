#!/usr/bin/env python3
"""EHP114 Fourier/Hilbert Hessian probe.

This is a diagnostic, not a proof.

It tests whether the deficit functional

    D_n(p) = L(z^n - 1) - L(p)

has positive finite-difference curvature along Fourier-mode perturbations of
the regular root polygon. The goal is to replace raw coefficient-coordinate
checks with a Hilbert/Fourier basis that can later be upgraded to interval
arithmetic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


DEFAULT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS = DEFAULT_ROOT / "erdos-experiments/results/erdos-114"
DEFAULT_SHORTCUT_RESULTS = DEFAULT_ROOT / "erdos-experiments/scripts/erdos-114"
EXPERIMENT_ID = "EXP-MATH-EHP114-N15-FOURIER-HESSIAN-20260502-02"


@dataclass(frozen=True)
class BasisVector:
    label: str
    mode: int
    kind: str
    phase: str
    vector: np.ndarray


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json_or_none(path: Path) -> dict | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def eval_poly_grid(z: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    """Evaluate z^n + a_{n-2} z^{n-2} + ... + a_0 on a grid.

    coeffs are ascending: [a_0, a_1, ..., a_{n-2}]. The z^(n-1)
    coefficient is assumed zero by centering the roots.
    """
    n = len(coeffs) + 1
    acc = np.ones_like(z, dtype=np.complex128)
    acc = acc * z
    for k in range(n - 2, -1, -1):
        acc = acc * z + coeffs[k]
    return acc.real * acc.real + acc.imag * acc.imag - 1.0


def lemniscate_length(coeffs: np.ndarray, res: int, extent: float) -> float:
    step = 2.0 * extent / res
    xs = np.linspace(-extent, extent, res + 1)
    ys = np.linspace(-extent, extent, res + 1)
    x, y = np.meshgrid(xs, ys, indexing="ij")
    z = x + 1j * y
    f = eval_poly_grid(z, coeffs)

    fsw = f[:-1, :-1]
    fse = f[1:, :-1]
    fne = f[1:, 1:]
    fnw = f[:-1, 1:]

    case = (
        (fsw > 0).astype(np.uint8)
        | ((fse > 0).astype(np.uint8) << 1)
        | ((fne > 0).astype(np.uint8) << 2)
        | ((fnw > 0).astype(np.uint8) << 3)
    )
    mask = (case != 0) & (case != 15)
    if not np.any(mask):
        return 0.0

    idx = np.where(mask)
    c = case[idx]
    sw = fsw[idx]
    se = fse[idx]
    ne_v = fne[idx]
    nw = fnw[idx]
    x0 = xs[idx[0]]
    y0 = ys[idx[1]]

    def interp(fa: np.ndarray, fb: np.ndarray) -> np.ndarray:
        d = fa - fb
        out = np.full_like(fa, 0.5, dtype=np.float64)
        np.divide(fa, d, out=out, where=np.abs(d) >= 1e-30)
        return out

    sx = x0 + interp(sw, se) * step
    sy = y0
    ex = x0 + step
    ey = y0 + interp(se, ne_v) * step
    nx = x0 + interp(nw, ne_v) * step
    ny = y0 + step
    wx = x0
    wy = y0 + interp(sw, nw) * step

    def seg(ax: np.ndarray, ay: np.ndarray, bx: np.ndarray, by: np.ndarray) -> np.ndarray:
        return np.sqrt((ax - bx) ** 2 + (ay - by) ** 2)

    total = np.zeros(len(c), dtype=np.float64)
    avg = (sw + se + ne_v + nw) / 4.0

    for cases, a, b in [
        ([1, 14], (sx, sy), (wx, wy)),
        ([2, 13], (sx, sy), (ex, ey)),
        ([3, 12], (wx, wy), (ex, ey)),
        ([4, 11], (ex, ey), (nx, ny)),
        ([6, 9], (sx, sy), (nx, ny)),
        ([7, 8], (wx, wy), (nx, ny)),
    ]:
        m = np.isin(c, cases)
        if np.any(m):
            total[m] = seg(a[0][m], a[1][m], b[0][m], b[1][m])

    m5 = c == 5
    if np.any(m5):
        pos = avg[m5] > 0
        la = seg(sx[m5], sy[m5], wx[m5], wy[m5]) + seg(ex[m5], ey[m5], nx[m5], ny[m5])
        lb = seg(sx[m5], sy[m5], ex[m5], ey[m5]) + seg(wx[m5], wy[m5], nx[m5], ny[m5])
        total[m5] = np.where(pos, la, lb)

    m10 = c == 10
    if np.any(m10):
        pos = avg[m10] > 0
        la = seg(sx[m10], sy[m10], ex[m10], ey[m10]) + seg(wx[m10], wy[m10], nx[m10], ny[m10])
        lb = seg(sx[m10], sy[m10], wx[m10], wy[m10]) + seg(ex[m10], ey[m10], nx[m10], ny[m10])
        total[m10] = np.where(pos, la, lb)

    return float(total.sum())


def roots_of_unity(n: int) -> np.ndarray:
    j = np.arange(n, dtype=np.float64)
    return np.exp(2j * np.pi * j / n)


def centered_coeffs_from_roots(roots: np.ndarray) -> tuple[np.ndarray, float]:
    centered = roots - roots.mean()
    poly = np.poly(centered)
    z_n_minus_1_coeff = poly[1]
    coeffs = poly[2:][::-1].astype(np.complex128)
    return coeffs, float(abs(z_n_minus_1_coeff))


def candidate_directions(n: int) -> list[BasisVector]:
    omega = roots_of_unity(n)
    theta = 2.0 * np.pi * np.arange(n, dtype=np.float64) / n
    candidates: list[BasisVector] = []

    candidates.append(BasisVector("m0_radial", 0, "radial", "constant", omega.copy()))
    candidates.append(BasisVector("m0_tangent_rotation", 0, "tangent", "constant", 1j * omega))

    for m in range(1, n // 2 + 1):
        cos_m = np.cos(m * theta)
        sin_m = np.sin(m * theta)
        for phase, scalar in [("cos", cos_m), ("sin", sin_m)]:
            candidates.append(BasisVector(f"m{m}_{phase}_radial", m, "radial", phase, omega * scalar))
            candidates.append(BasisVector(f"m{m}_{phase}_tangent", m, "tangent", phase, 1j * omega * scalar))

    return candidates


def as_real_vector(v: np.ndarray) -> np.ndarray:
    return np.concatenate([v.real, v.imag])


def as_complex_vector(v: np.ndarray) -> np.ndarray:
    n2 = len(v) // 2
    return v[:n2] + 1j * v[n2:]


def orthonormal_fourier_basis(n: int, tol: float = 1e-10) -> list[BasisVector]:
    """Build a centered real Hilbert basis from Fourier root perturbations.

    Translation directions are projected out by subtracting the complex mean.
    The global rotation direction is omitted because it is a length symmetry.
    """
    basis: list[np.ndarray] = []
    meta: list[tuple[str, int, str, str]] = []
    for cand in candidate_directions(n):
        if cand.label == "m0_tangent_rotation":
            continue
        v = cand.vector - cand.vector.mean()
        rv = as_real_vector(v)
        for q in basis:
            rv = rv - float(np.dot(q, rv)) * q
        norm = float(np.linalg.norm(rv))
        if norm <= tol:
            continue
        basis.append(rv / norm)
        meta.append((cand.label, cand.mode, cand.kind, cand.phase))

    return [
        BasisVector(label=label, mode=mode, kind=kind, phase=phase, vector=as_complex_vector(rv))
        for rv, (label, mode, kind, phase) in zip(basis, meta)
    ]


def run_probe(
    n: int,
    eps_values: list[float],
    res: int,
    extent: float,
    shortcut_dir: Path,
) -> dict:
    omega = roots_of_unity(n)
    base_coeffs, base_center_error = centered_coeffs_from_roots(omega)
    l0 = lemniscate_length(base_coeffs, res=res, extent=extent)
    basis = orthonormal_fourier_basis(n)

    shortcut_path = shortcut_dir / f"EXP-MM-EHP-007-n{n}-inari_RESULTS.json"
    shortcut = load_json_or_none(shortcut_path)
    shortcut_l_star = shortcut.get("l_star_lower") if shortcut else None

    singular_radial_rows = []
    for eps in eps_values:
        radius_step = eps / math.sqrt(n)
        for sign in [-1.0, 1.0]:
            radius = 1.0 + sign * radius_step
            coeffs = np.zeros(n - 1, dtype=np.complex128)
            coeffs[0] = -(radius**n)
            singular_radial_rows.append(
                {
                    "eps": eps,
                    "radius": radius,
                    "length": lemniscate_length(coeffs, res=res, extent=extent),
                    "deficit_vs_L0_marching": l0 - lemniscate_length(coeffs, res=res, extent=extent),
                }
            )

    rows: list[dict] = []
    for b in basis:
        eps_rows = []
        for eps in eps_values:
            coeffs_p, center_error_p = centered_coeffs_from_roots(omega + eps * b.vector)
            coeffs_m, center_error_m = centered_coeffs_from_roots(omega - eps * b.vector)
            lp = lemniscate_length(coeffs_p, res=res, extent=extent)
            lm = lemniscate_length(coeffs_m, res=res, extent=extent)
            deficit_curvature = (2.0 * l0 - lp - lm) / (eps * eps)
            length_second_difference = -deficit_curvature
            eps_rows.append(
                {
                    "eps": eps,
                    "L_plus": lp,
                    "L_minus": lm,
                    "deficit_curvature": deficit_curvature,
                    "length_second_difference": length_second_difference,
                    "positive_deficit_curvature": deficit_curvature > 0.0,
                    "center_error_plus": center_error_p,
                    "center_error_minus": center_error_m,
                }
            )
        curvatures = [row["deficit_curvature"] for row in eps_rows]
        rows.append(
            {
                "label": b.label,
                "mode": b.mode,
                "kind": b.kind,
                "phase": b.phase,
                "eps_rows": eps_rows,
                "min_deficit_curvature": min(curvatures),
                "max_deficit_curvature": max(curvatures),
                "all_eps_positive": all(c > 0.0 for c in curvatures),
            }
        )

    positive_all = [r for r in rows if r["all_eps_positive"]]
    worst = min(rows, key=lambda r: r["min_deficit_curvature"])
    best = max(rows, key=lambda r: r["max_deficit_curvature"])

    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "status": "DIAGNOSTIC_NUMERICAL_PROBE_NOT_PROOF",
        "degree": n,
        "method": {
            "functional": "D_n(p) = L(z^n - 1) - L(p)",
            "coordinate_lift": "centered root perturbations in finite real Hilbert space",
            "basis": "Fourier radial/tangent modes, projected off translation; global rotation omitted",
            "curvature_estimator": "(2*L0 - L(+eps*v) - L(-eps*v)) / eps^2",
            "length_estimator": "floating-point marching squares",
            "rigorous": False,
        },
        "parameters": {
            "eps_values": eps_values,
            "grid_resolution": res,
            "extent": extent,
            "basis_rank": len(basis),
            "expected_reduced_dim_2n_minus_3": 2 * n - 3,
        },
        "reference": {
            "L0_marching_squares": l0,
            "L0_relative_error_vs_shortcut_lstar": (
                abs(l0 - shortcut_l_star) / shortcut_l_star if shortcut_l_star else None
            ),
            "base_centering_residual_coeff_abs": base_center_error,
            "shortcut_artifact": str(shortcut_path) if shortcut_path.exists() else None,
            "shortcut_artifact_sha256": sha256_file(shortcut_path) if shortcut_path.exists() else None,
            "shortcut_l_star_lower": shortcut.get("l_star_lower") if shortcut else None,
            "shortcut_l_star_upper": shortcut.get("l_star_upper") if shortcut else None,
            "shortcut_verdict_not_used_as_proof": shortcut.get("verdict") if shortcut else None,
        },
        "singularity_sanity": {
            "why": (
                "z^n - 1 has a multiple critical point at the origin lying on the "
                "lemniscate. Small radial/constant perturbations unfold this singular "
                "level set, so an ordinary smooth Hessian interpretation is unsafe."
            ),
            "radial_family": "p_R(z) = z^n - R^n, with R = 1 +/- eps/sqrt(n)",
            "rows": singular_radial_rows,
            "large_reference_error_flag": (
                abs(l0 - shortcut_l_star) / shortcut_l_star > 0.05 if shortcut_l_star else None
            ),
        },
        "basis_rows": rows,
        "summary": {
            "basis_rows": len(rows),
            "rows_positive_all_eps": len(positive_all),
            "rows_not_positive_all_eps": len(rows) - len(positive_all),
            "worst_label": worst["label"],
            "worst_min_deficit_curvature": worst["min_deficit_curvature"],
            "best_label": best["label"],
            "best_max_deficit_curvature": best["max_deficit_curvature"],
            "all_modes_positive_all_eps": len(positive_all) == len(rows),
        },
        "interpretation": {
            "primary": (
                "Positive deficit curvature means the regular root polygon locally "
                "beats the tested Fourier perturbation direction under this floating-point estimator."
            ),
            "caution": (
                "This is not an interval certificate and does not prove local maximality. "
                "Because z^n - 1 is a singular lemniscate, the measured effect is better "
                "read as a stratified deficit/unfolding signal than as an ordinary smooth "
                "Hessian certificate."
            ),
            "next_step": (
                "If all or most directions are positive, upgrade this to an interval "
                "or analytic Fourier-mode Hessian certificate and then bound tensor remainders."
            ),
        },
    }


def write_report(result: dict, report_path: Path) -> None:
    summary = result["summary"]
    ref = result["reference"]
    sanity = result["singularity_sanity"]
    lines = [
        f"# {result['experiment_id']} Report",
        "",
        "## Claim Tested",
        "",
        "Can the EHP114 deficit near `z^n - 1` be seen in Hilbert/Fourier coordinates",
        "instead of raw coefficient coordinates?",
        "",
        "This is a diagnostic numerical probe, not a proof.",
        "",
        "## Result",
        "",
        f"- Degree: {result['degree']}",
        f"- Basis rank: {result['parameters']['basis_rank']} / expected `2n-3` = {result['parameters']['expected_reduced_dim_2n_minus_3']}",
        f"- Grid resolution: {result['parameters']['grid_resolution']}",
        f"- eps values: {result['parameters']['eps_values']}",
        f"- L0 marching-squares estimate: {ref['L0_marching_squares']:.12f}",
        f"- L0 relative error vs shortcut L*: {ref['L0_relative_error_vs_shortcut_lstar']:.6f}",
        f"- Rows positive for all eps: {summary['rows_positive_all_eps']}/{summary['basis_rows']}",
        f"- Worst mode: `{summary['worst_label']}` with min curvature {summary['worst_min_deficit_curvature']:.6e}",
        f"- Best mode: `{summary['best_label']}` with max curvature {summary['best_max_deficit_curvature']:.6e}",
        f"- All tested modes positive: {summary['all_modes_positive_all_eps']}",
        "",
        "## Interpretation",
        "",
        result["interpretation"]["primary"],
        "",
        result["interpretation"]["caution"],
        "",
        result["interpretation"]["next_step"],
        "",
        "## Singularity Sanity Check",
        "",
        sanity["why"],
        "",
        f"- Radial family checked: `{sanity['radial_family']}`",
        f"- Large reference-error flag: {sanity['large_reference_error_flag']}",
        "",
        "| eps | radius | length | deficit vs L0 marching |",
        "|---:|---:|---:|---:|",
    ]
    for row in sanity["rows"]:
        lines.append(
            f"| {row['eps']:.6g} | {row['radius']:.12f} | {row['length']:.12f} | {row['deficit_vs_L0_marching']:.12f} |"
        )
    lines.extend(
        [
            "",
            "## Mode Table",
            "",
            "| label | mode | kind | phase | min curvature | max curvature | all eps positive? |",
            "|---|---:|---|---|---:|---:|---|",
        ]
    )
    for row in result["basis_rows"]:
        lines.append(
            "| {label} | {mode} | {kind} | {phase} | {mn:.6e} | {mx:.6e} | {ok} |".format(
                label=row["label"],
                mode=row["mode"],
                kind=row["kind"],
                phase=row["phase"],
                mn=row["min_deficit_curvature"],
                mx=row["max_deficit_curvature"],
                ok="yes" if row["all_eps_positive"] else "no",
            )
        )
    lines.extend(
        [
            "",
            "## Source Boundary",
            "",
            f"- Shortcut L* artifact consulted only for reference: `{ref['shortcut_artifact']}`",
            f"- Shortcut artifact SHA-256: `{ref['shortcut_artifact_sha256']}`",
            "",
            "The shortcut artifact verdict is not used as proof in this probe.",
            "",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=15)
    parser.add_argument("--res", type=int, default=360)
    parser.add_argument("--extent", type=float, default=3.0)
    parser.add_argument("--eps", default="0.02,0.01,0.005")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--shortcut-dir", type=Path, default=DEFAULT_SHORTCUT_RESULTS)
    args = parser.parse_args()

    eps_values = [float(x.strip()) for x in args.eps.split(",") if x.strip()]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    result_path = args.out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = args.out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = args.out_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"

    for path in (result_path, report_path, sha_path):
        if path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {path}")

    result = run_probe(
        n=args.n,
        eps_values=eps_values,
        res=args.res,
        extent=args.extent,
        shortcut_dir=args.shortcut_dir,
    )
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(result, report_path)
    sha_path.write_text(sha256_file(result_path) + "\n", encoding="utf-8")

    print(json.dumps(result["summary"], indent=2, sort_keys=True))
    print(f"Wrote {result_path}")
    print(f"Wrote {report_path}")
    print(f"Wrote {sha_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
