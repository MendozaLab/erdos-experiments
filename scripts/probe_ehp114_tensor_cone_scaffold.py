#!/usr/bin/env python3
"""EHP114 tensor-cone scaffold probe.

This is a diagnostic scaffold, not a proof.

Purpose:
  * build a normalized quotient tangent basis around the regular root polygon,
  * separate the singular radial direction from shape-changing directions,
  * test whether finite-difference deficit behavior is Hessian-like or
    Puiseux/stratified, and
  * emit a versioned artifact that can guide the interval/tensor program.

The key functional is

    D_n(x) = L(z^n - 1) - L(p_x).

Lengths for perturbed configurations are estimated by floating marching
squares; the reference L(z^n - 1) is evaluated by the exact gamma formula.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np


EXPERIMENT_ID = "EXP-MATH-EHP114-N10-TENSOR-CONE-SCAFFOLD-20260502-02"
DEFAULT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS = DEFAULT_ROOT / "erdos-experiments/results/erdos-114"
DEFAULT_EPS = "0.04,0.02,0.01,0.005"


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


def load_json_or_none(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def parse_float_list(raw: str) -> list[float]:
    values = [float(part.strip()) for part in raw.split(",") if part.strip()]
    if not values:
        raise ValueError("empty eps list")
    return values


def exact_l_star(n: int) -> mp.mpf:
    alpha = mp.mpf(1) / n - 1
    return (
        mp.power(2, alpha + 2)
        * mp.sqrt(mp.pi)
        * mp.gamma((alpha + 1) / 2)
        / (2 * mp.gamma((alpha + 2) / 2))
    )


def eval_poly_grid(z: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    """Evaluate z^n + a_(n-1) z^(n-1) + ... + a_0 on a grid.

    coeffs are ascending: [a_0, a_1, ..., a_(n-1)]. The leading coefficient is
    monic and omitted.
    """
    acc = np.ones_like(z, dtype=np.complex128)
    for k in range(len(coeffs) - 1, -1, -1):
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


def coeffs_from_roots(roots: np.ndarray) -> np.ndarray:
    return np.poly(roots)[1:][::-1].astype(np.complex128)


def as_real_vector(v: np.ndarray) -> np.ndarray:
    return np.concatenate([v.real, v.imag])


def as_complex_vector(v: np.ndarray) -> np.ndarray:
    half = len(v) // 2
    return v[:half] + 1j * v[half:]


def candidate_directions(n: int) -> list[BasisVector]:
    omega = roots_of_unity(n)
    theta = 2.0 * np.pi * np.arange(n, dtype=np.float64) / n
    candidates = [
        BasisVector("radial_singular_m0", 0, "radial", "constant", omega.copy()),
        BasisVector("rotation_symmetry_m0", 0, "tangent", "constant", 1j * omega),
    ]
    for m in range(1, n // 2 + 1):
        cos_m = np.cos(m * theta)
        sin_m = np.sin(m * theta)
        for phase, scalar in [("cos", cos_m), ("sin", sin_m)]:
            candidates.append(BasisVector(f"m{m}_{phase}_radial", m, "radial", phase, omega * scalar))
            candidates.append(BasisVector(f"m{m}_{phase}_tangent", m, "tangent", phase, 1j * omega * scalar))
    return candidates


def quotient_basis(n: int, tol: float = 1e-10) -> list[BasisVector]:
    """Orthonormal basis after translation and rotation quotient.

    The singular radial vector is kept and marked. Translation is removed by
    subtracting the complex mean from each candidate; global rotation is omitted.
    """
    basis: list[np.ndarray] = []
    meta: list[tuple[str, int, str, str]] = []
    for cand in candidate_directions(n):
        if cand.label == "rotation_symmetry_m0":
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


def fit_log_slope(eps_values: list[float], values: list[float]) -> float | None:
    positive = [(e, v) for e, v in zip(eps_values, values) if e > 0.0 and v > 0.0 and math.isfinite(v)]
    if len(positive) < 3:
        return None
    xs = [math.log(e) for e, _ in positive]
    ys = [math.log(v) for _, v in positive]
    xm = sum(xs) / len(xs)
    ym = sum(ys) / len(ys)
    denom = sum((x - xm) ** 2 for x in xs)
    if denom <= 0.0:
        return None
    return sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / denom


def direction_deficits(
    roots: np.ndarray,
    direction: np.ndarray,
    l_star: float,
    eps_values: list[float],
    res: int,
    extent: float,
) -> list[dict[str, Any]]:
    rows = []
    for eps in eps_values:
        roots_p = roots + eps * direction
        roots_m = roots - eps * direction
        lp = lemniscate_length(coeffs_from_roots(roots_p), res=res, extent=extent)
        lm = lemniscate_length(coeffs_from_roots(roots_m), res=res, extent=extent)
        dp = l_star - lp
        dm = l_star - lm
        symmetric_deficit = 0.5 * (dp + dm)
        rows.append(
            {
                "eps": eps,
                "L_plus": lp,
                "L_minus": lm,
                "D_plus": dp,
                "D_minus": dm,
                "symmetric_deficit": symmetric_deficit,
                "hessian_proxy_D_over_eps2": symmetric_deficit / (eps * eps),
                "root_radius_max_plus": float(np.max(np.abs(roots_p))),
                "root_radius_max_minus": float(np.max(np.abs(roots_m))),
                "roots_inside_unit_disk_plus": bool(np.max(np.abs(roots_p)) <= 1.0 + 1e-12),
                "roots_inside_unit_disk_minus": bool(np.max(np.abs(roots_m)) <= 1.0 + 1e-12),
            }
        )
    return rows


def mixed_hessian_proxy(
    roots: np.ndarray,
    basis: list[BasisVector],
    l_star: float,
    eps: float,
    res: int,
    extent: float,
) -> dict[str, Any]:
    dim = len(basis)
    matrix = np.zeros((dim, dim), dtype=np.float64)

    def deficit_for(vec: np.ndarray) -> float:
        length = lemniscate_length(coeffs_from_roots(roots + eps * vec), res=res, extent=extent)
        return l_star - length

    single = []
    for b in basis:
        d_plus = deficit_for(b.vector)
        d_minus = deficit_for(-b.vector)
        value = 0.5 * (d_plus + d_minus) / (eps * eps)
        single.append(value)
    for i, value in enumerate(single):
        matrix[i, i] = value

    for i in range(dim):
        vi = basis[i].vector
        for j in range(i + 1, dim):
            vj = basis[j].vector
            dpp = deficit_for(vi + vj)
            dpm = deficit_for(vi - vj)
            dmp = deficit_for(-vi + vj)
            dmm = deficit_for(-vi - vj)
            value = (dpp - dpm - dmp + dmm) / (4 * eps * eps)
            matrix[i, j] = value
            matrix[j, i] = value

    eigvals = np.linalg.eigvalsh(matrix)
    return {
        "eps": eps,
        "basis_labels": [b.label for b in basis],
        "matrix": matrix.tolist(),
        "eigenvalues": eigvals.tolist(),
        "min_eigenvalue": float(eigvals[0]),
        "max_eigenvalue": float(eigvals[-1]),
        "condition_proxy_abs_max_over_min": (
            float(np.max(np.abs(eigvals)) / abs(eigvals[0])) if abs(eigvals[0]) > 0 else None
        ),
    }


def reference_compare(n: int, result_dir: Path, exact_value: mp.mpf) -> dict[str, Any]:
    path = result_dir / f"EXP-MM-EHP-007-n{n}-inari_RESULTS.json"
    artifact = load_json_or_none(path)
    if artifact is None:
        return {"artifact_path": str(path), "found": False}
    lower = mp.mpf(str(artifact.get("l_star_lower")))
    upper = mp.mpf(str(artifact.get("l_star_upper")))
    return {
        "artifact_path": str(path),
        "artifact_sha256": sha256_file(path),
        "found": True,
        "verdict": artifact.get("verdict"),
        "rigor": artifact.get("rigor"),
        "reduced_dim": artifact.get("reduced_dim"),
        "bb_proof_complete": artifact.get("bb_proof_complete"),
        "l_star_lower": str(artifact.get("l_star_lower")),
        "l_star_upper": str(artifact.get("l_star_upper")),
        "exact_inside_interval": bool(lower <= exact_value <= upper),
        "exact_minus_midpoint": mp.nstr(exact_value - (lower + upper) / 2, 30),
    }


def run_probe(
    n: int,
    eps_values: list[float],
    res: int,
    mixed_res: int,
    extent: float,
    result_dir: Path,
    include_mixed_matrix: bool,
) -> dict[str, Any]:
    roots = roots_of_unity(n)
    l_star_mp = exact_l_star(n)
    l_star = float(l_star_mp)
    basis = quotient_basis(n)
    radial = [b for b in basis if b.label == "radial_singular_m0"]
    shape = [b for b in basis if b.label != "radial_singular_m0"]

    direction_rows = []
    for b in basis:
        rows = direction_deficits(roots, b.vector, l_star, eps_values, res=res, extent=extent)
        symmetric = [row["symmetric_deficit"] for row in rows]
        slope = fit_log_slope(eps_values, symmetric)
        direction_rows.append(
            {
                "label": b.label,
                "mode": b.mode,
                "kind": b.kind,
                "phase": b.phase,
                "cone_role": "singular_radial" if b.label == "radial_singular_m0" else "shape_candidate",
                "rows": rows,
                "fitted_loglog_slope_symmetric_deficit": slope,
                "smallest_eps_hessian_proxy": rows[-1]["hessian_proxy_D_over_eps2"],
                "all_symmetric_deficits_positive": all(row["symmetric_deficit"] > 0.0 for row in rows),
                "all_tested_plus_minus_inside_unit_disk": all(
                    row["roots_inside_unit_disk_plus"] and row["roots_inside_unit_disk_minus"]
                    for row in rows
                ),
            }
        )

    shape_slopes = [
        row["fitted_loglog_slope_symmetric_deficit"]
        for row in direction_rows
        if row["cone_role"] == "shape_candidate" and row["fitted_loglog_slope_symmetric_deficit"] is not None
    ]
    radial_slope = (
        direction_rows[0]["fitted_loglog_slope_symmetric_deficit"]
        if direction_rows and direction_rows[0]["label"] == "radial_singular_m0"
        else None
    )

    mixed = None
    if include_mixed_matrix:
        mixed = mixed_hessian_proxy(
            roots=roots,
            basis=shape,
            l_star=l_star,
            eps=eps_values[-1],
            res=mixed_res,
            extent=extent,
        )

    hessian_like_threshold = 1.5
    shape_hessian_like_count = sum(1 for slope in shape_slopes if slope >= hessian_like_threshold)
    shape_singular_like_count = len(shape_slopes) - shape_hessian_like_count

    return {
        "experiment_id": EXPERIMENT_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "status": "DIAGNOSTIC_TENSOR_CONE_SCAFFOLD_NOT_PROOF",
        "degree": n,
        "method": {
            "functional": "D_n(x) = L(z^n - 1) - L(p_x)",
            "reference_length": "exact gamma formula for L(z^n - 1)",
            "perturbed_length_estimator": "floating-point marching squares",
            "basis": "quotient root tangent basis: translation projected out, rotation omitted, singular radial retained",
            "shape_cone_probe": "all quotient directions except singular radial m0",
            "unit_disk_flags": (
                "diagnostic only; the EHP monic-polynomial problem does not impose "
                "a roots-in-unit-disk constraint"
            ),
            "rigorous": False,
        },
        "parameters": {
            "eps_values": eps_values,
            "grid_resolution": res,
            "mixed_matrix_grid_resolution": mixed_res if include_mixed_matrix else None,
            "extent": extent,
            "quotient_basis_rank": len(basis),
            "shape_basis_rank": len(shape),
            "radial_basis_rank": len(radial),
            "expected_rank_2n_minus_3": 2 * n - 3,
        },
        "reference": {
            "exact_l_star": mp.nstr(l_star_mp, 50),
            "artifact_compare": reference_compare(n, result_dir, l_star_mp),
        },
        "direction_rows": direction_rows,
        "mixed_shape_hessian_proxy": mixed,
        "summary": {
            "radial_slope": radial_slope,
            "shape_slope_min": min(shape_slopes) if shape_slopes else None,
            "shape_slope_max": max(shape_slopes) if shape_slopes else None,
            "shape_slope_mean": sum(shape_slopes) / len(shape_slopes) if shape_slopes else None,
            "shape_hessian_like_count_slope_ge_1_5": shape_hessian_like_count,
            "shape_singular_like_count_slope_lt_1_5": shape_singular_like_count,
            "all_directions_positive_symmetric_deficit": all(
                row["all_symmetric_deficits_positive"] for row in direction_rows
            ),
            "mixed_min_eigenvalue": mixed["min_eigenvalue"] if mixed else None,
            "mixed_max_eigenvalue": mixed["max_eigenvalue"] if mixed else None,
        },
        "interpretation": {
            "primary": (
                "If shape slopes are near 2, a smooth tensor/Hessian shape cone may be "
                "viable. If slopes are well below 2, the local proof must be stratified "
                "rather than an ordinary tensor Hessian."
            ),
            "claim_boundary": (
                "This scaffold uses floating marching-squares estimates for perturbed "
                "lengths and is not an interval certificate."
            ),
            "next_step": (
                "If n=10 shows stable positive deficit and a coherent scaling class, "
                "repeat for n=11-14, then interval-harden the first case whose structure "
                "is representative."
            ),
        },
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    summary = result["summary"]
    ref = result["reference"]
    lines = [
        f"# {result['experiment_id']} Report",
        "",
        "## Status",
        "",
        "Diagnostic tensor-cone scaffold. Not a proof.",
        "",
        "## Claim Tested",
        "",
        "Does the `n = 10` quotient tangent cone look like an ordinary smooth Hessian",
        "problem, or does it already show stratified/Puiseux behavior?",
        "",
        "## Reference",
        "",
        f"- Degree: {result['degree']}",
        f"- Exact L*: `{ref['exact_l_star']}`",
        f"- Reference artifact: `{ref['artifact_compare'].get('artifact_path')}`",
        f"- Artifact verdict: `{ref['artifact_compare'].get('verdict')}`",
        f"- Exact value inside artifact interval: {ref['artifact_compare'].get('exact_inside_interval')}",
        "",
        "## Basis",
        "",
        f"- Quotient basis rank: {result['parameters']['quotient_basis_rank']}",
        f"- Expected rank `2n - 3`: {result['parameters']['expected_rank_2n_minus_3']}",
        f"- Shape basis rank after removing singular radial: {result['parameters']['shape_basis_rank']}",
        f"- eps values: {result['parameters']['eps_values']}",
        f"- Grid resolution: {result['parameters']['grid_resolution']}",
        "",
        "## Summary",
        "",
        f"- Radial slope: {summary['radial_slope']:.6f}",
        f"- Shape slope min: {summary['shape_slope_min']:.6f}",
        f"- Shape slope mean: {summary['shape_slope_mean']:.6f}",
        f"- Shape slope max: {summary['shape_slope_max']:.6f}",
        f"- Shape directions with slope >= 1.5: {summary['shape_hessian_like_count_slope_ge_1_5']}",
        f"- Shape directions with slope < 1.5: {summary['shape_singular_like_count_slope_lt_1_5']}",
        f"- All directions have positive symmetric deficit: {summary['all_directions_positive_symmetric_deficit']}",
        "",
    ]
    if result["mixed_shape_hessian_proxy"] is not None:
        mixed = result["mixed_shape_hessian_proxy"]
        lines.extend(
            [
                "## Mixed Shape Hessian Proxy",
                "",
                f"- eps: {mixed['eps']}",
                f"- min eigenvalue: {mixed['min_eigenvalue']:.6e}",
                f"- max eigenvalue: {mixed['max_eigenvalue']:.6e}",
                f"- condition proxy abs(max)/abs(min): {mixed['condition_proxy_abs_max_over_min']:.6e}",
                "",
                "This matrix is a finite-difference proxy, not an interval Hessian.",
                "",
            ]
        )

    lines.extend(
        [
            "## Direction Table",
            "",
            "| label | role | mode | kind | phase | fitted slope | smallest-eps D/eps^2 | positive? |",
            "|---|---|---:|---|---|---:|---:|---|",
        ]
    )
    for row in result["direction_rows"]:
        slope = row["fitted_loglog_slope_symmetric_deficit"]
        lines.append(
            "| {label} | {role} | {mode} | {kind} | {phase} | {slope:.6f} | {proxy:.6e} | {ok} |".format(
                label=row["label"],
                role=row["cone_role"],
                mode=row["mode"],
                kind=row["kind"],
                phase=row["phase"],
                slope=slope if slope is not None else float("nan"),
                proxy=row["smallest_eps_hessian_proxy"],
                ok="yes" if row["all_symmetric_deficits_positive"] else "no",
            )
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            result["interpretation"]["primary"],
            "",
            result["interpretation"]["claim_boundary"],
            "",
            result["interpretation"]["next_step"],
            "",
            "## Source Boundary",
            "",
            f"- Reference artifact SHA-256: `{ref['artifact_compare'].get('artifact_sha256')}`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--eps", default=DEFAULT_EPS)
    parser.add_argument("--res", type=int, default=340)
    parser.add_argument("--mixed-res", type=int, default=240)
    parser.add_argument("--extent", type=float, default=3.0)
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--skip-mixed-matrix", action="store_true")
    args = parser.parse_args()

    if args.n != 10:
        raise SystemExit("This version is pinned to n=10; version the experiment id before changing n.")

    mp.mp.dps = args.mp_dps
    args.out_dir.mkdir(parents=True, exist_ok=True)
    result_path = args.out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = args.out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = args.out_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"
    for artifact_path in [result_path, report_path, sha_path]:
        if artifact_path.exists():
            raise SystemExit(f"Refusing to overwrite existing artifact: {artifact_path}")

    result = run_probe(
        n=args.n,
        eps_values=parse_float_list(args.eps),
        res=args.res,
        mixed_res=args.mixed_res,
        extent=args.extent,
        result_dir=args.out_dir,
        include_mixed_matrix=not args.skip_mixed_matrix,
    )
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_report(result, report_path)
    sha_path.write_text(f"{sha256_file(result_path)}  {result_path.name}\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "result": str(result_path),
                "report": str(report_path),
                "sha256": sha_path.read_text().split()[0],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
