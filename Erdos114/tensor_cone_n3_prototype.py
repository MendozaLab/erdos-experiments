#!/usr/bin/env python3
"""Tensor-cone scoping prototype for Erdős #114 (EHP) at n=3.

Experiment ID: EXP-MATH-EHP114-TENSOR-CONE-V1-20260502

This is a SCOPING artifact — see TENSOR_CONE_DESIGN_2026-05-02.md for the
honest scope statement. The script:

  1. Defines the candidate cone K_3 = K_(1)^{\otimes 3} ∩ {trace ≤ 1} on the
     real coefficient slice (x_0, y_0, x_1, y_1, x_2, y_2) ∈ R^6.

  2. Verifies that a* = (-1, 0, 0, 0, 0, 0) — corresponding to z^3 - 1 — sits
     at the cone boundary (both atomic 2×2 PSD blocks rank-deficient AND the
     trace constraint saturated).

  3. Tests three alternative polynomials inside the constrained region.
     For each, we (a) check cone membership, (b) compute L(p) numerically via
     adaptive boundary tracing of {|p|=1}, (c) confirm L(p) ≤ L(z^3 - 1).

  4. Emits EXP-MATH-EHP114-TENSOR-CONE-V1-20260502_RESULTS.json with cone
     definition, boundary verification, and three alternative-poly checks.

NOTE on numerical L(p): the v5 preprint uses certified IEEE-1788 interval
arithmetic for L. We use plain double-precision with a documented contour-
length quadrature and report the value as best-estimate, NOT certified. This
is consistent with the scoping label.

Dependencies: numpy, scipy. cvxpy is NOT available in this environment, so
the SDP feasibility check is performed as the primal eigenvalue check on each
atomic block plus the trace inequality (this is exact for the diagonal cone
defined in the design doc — no SDP solver is needed because the cone factors
into independent 2×2 blocks).
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from dataclasses import dataclass
from datetime import datetime, timezone

import numpy as np
from scipy.linalg import eigvalsh

# ---------------------------------------------------------------------------
# Constants — closed form L(z^n - 1) (Theorem in v5 preprint, §2.2)
# ---------------------------------------------------------------------------

def L_zn_minus_1(n: int) -> float:
    """Closed-form lemniscate length of z^n - 1 (v5 preprint Eq. 1)."""
    return (2.0 ** (1.0 / n)) * math.sqrt(math.pi) * math.gamma(1.0 / (2 * n)) / math.gamma(1.0 / (2 * n) + 0.5)


# ---------------------------------------------------------------------------
# Atomic block K_(1) — coordinate-PSD constraint per coefficient
# ---------------------------------------------------------------------------

def atomic_block_eigs(x: float, y: float) -> tuple[np.ndarray, np.ndarray]:
    """Return eigenvalues of the two 2×2 PSD blocks for one coefficient.

    Block_x = [[1, x], [x, 1]]  has eigenvalues {1-x, 1+x}
    Block_y = [[1, y], [y, 1]]  has eigenvalues {1-y, 1+y}

    Smallest eigenvalue ≥ 0 iff |x| ≤ 1 (resp. |y| ≤ 1).
    """
    Bx = np.array([[1.0, x], [x, 1.0]])
    By = np.array([[1.0, y], [y, 1.0]])
    return eigvalsh(Bx), eigvalsh(By)


def cone_membership(coeffs: list[complex], tol: float = 1e-12) -> dict:
    """Check membership of a coefficient vector in K_n.

    coeffs: list of (a_0, a_1, ..., a_{n-1})  — leading coefficient is implicit 1.

    Returns dict with per-block min-eigenvalue, trace, and overall verdict.
    """
    n = len(coeffs)
    block_min_eigs: list[float] = []
    trace = 0.0
    for k, ak in enumerate(coeffs):
        x, y = ak.real, ak.imag
        ex, ey = atomic_block_eigs(x, y)
        block_min_eigs.append(float(min(ex.min(), ey.min())))
        trace += x * x + y * y

    min_block_eig = min(block_min_eigs)
    in_cone_blocks = min_block_eig >= -tol
    in_cone_trace = trace <= 1.0 + tol

    boundary_blocks = abs(min_block_eig) <= tol
    boundary_trace = abs(trace - 1.0) <= tol

    return {
        "n": n,
        "block_min_eigenvalues": block_min_eigs,
        "min_block_eigenvalue": min_block_eig,
        "energy_trace": trace,
        "in_cone_per_block": in_cone_blocks,
        "in_cone_trace": in_cone_trace,
        "in_cone": in_cone_blocks and in_cone_trace,
        "on_boundary_blocks": boundary_blocks,
        "on_boundary_trace": boundary_trace,
        "on_boundary_overall": boundary_blocks and boundary_trace,
        "tol": tol,
    }


# ---------------------------------------------------------------------------
# Numerical L(p) — adaptive contour-length quadrature
# ---------------------------------------------------------------------------

def lemniscate_length_numerical(coeffs: list[complex], n_components_hint: int | None = None,
                                 grid_points: int = 4000, refine_passes: int = 2) -> dict:
    """Estimate L(p) = H^1({|p|=1}) numerically via co-area formula.

    Co-area formula for u(z) = |p(z)|:
        L({u=1}) = lim_{eps->0} (1 / (2 eps)) * Area({|u - 1| < eps})

    Equivalently, for a fine grid in z with spacing h covering a box B ⊃
    {|p| ≤ 1+eps}:
        L ≈ (h^2 / (2 eps)) * #{grid points where |1 - |p(z)|| < eps}
              * (correction for grad|p|)

    Better: the rigorous co-area form is
        ∫_R^2 |∇u|(z) δ(u(z) - 1) dA = L({u=1})
    so for a thin strip of half-width eps,
        L ≈ (1 / (2 eps)) * ∫_{|u-1|<eps} |∇u| dA
    Numerically:
        L ≈ (h^2 / (2 eps)) * Σ_{|u(z_i)-1|<eps} |∇u(z_i)|

    This is robust (no curve following), gives a deterministic estimate, and
    we report a bracket from two values of eps.

    SCOPING quality only — not certified.
    """
    n = len(coeffs)
    poly_coeffs = [1.0] + [coeffs[n - 1 - k] for k in range(n)]
    p = np.poly1d(poly_coeffs)
    dp = p.deriv()

    # Co-area: build a complex grid covering {|p| <= some bound} and integrate
    # |grad u| * indicator(|u-1| < eps) over the grid.
    # First find the bounding box: scan radii on rays out to where |p| > 1+eps.
    # For monic z^n + lower, |p(z)| ~ |z|^n at large |z|, so a generous box of
    # radius (1 + 1.5)^{1/n} * (n+1) is safe for our test polynomials with
    # |a_k| <= 1. Conservatively: box radius = 3.0.
    R = 3.0
    grid_n = 1200  # 1200x1200 grid → 1.44M points; ~few seconds in numpy
    xs = np.linspace(-R, R, grid_n)
    ys = np.linspace(-R, R, grid_n)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    Z = X + 1j * Y
    h = xs[1] - xs[0]

    # Vectorized polynomial evaluation
    pZ = np.zeros_like(Z, dtype=complex)
    dpZ = np.zeros_like(Z, dtype=complex)
    for c in poly_coeffs:
        pZ = pZ * Z + c
    # Derivative coeffs
    dp_coeffs = [poly_coeffs[i] * (n - i) for i in range(n)]  # leading n*z^{n-1} ...
    for c in dp_coeffs:
        dpZ = dpZ * Z + c

    abs_p = np.abs(pZ)
    # |grad u(z)| where u = |p(z)|: gradient of sqrt(p p_bar) in (x,y) gives
    # |grad u| = |p'(z)| (Cauchy-Riemann; standard).
    grad_u = np.abs(dpZ)

    # Sample a few epsilons and integrate.
    eps_values = [0.04, 0.02, 0.01]
    L_estimates = []
    for eps in eps_values:
        mask = np.abs(abs_p - 1.0) < eps
        # Riemann sum: L ≈ (1 / (2 eps)) * h^2 * Σ |grad u|
        L_eps = (h * h / (2.0 * eps)) * float(np.sum(grad_u[mask]))
        L_estimates.append({"eps": eps, "L_estimate": L_eps, "n_points_in_strip": int(mask.sum())})

    # Best estimate: smallest eps; bracket: range over eps values
    L_best = L_estimates[-1]["L_estimate"]
    L_min = min(e["L_estimate"] for e in L_estimates)
    L_max = max(e["L_estimate"] for e in L_estimates)

    return {
        "method": "co_area_grid_integration",
        "grid_size": grid_n,
        "box_radius": R,
        "grid_spacing": h,
        "L_estimate": L_best,
        "L_estimates_by_eps": L_estimates,
        "L_lower_bracket": L_min,
        "L_upper_bracket": L_max,
        "completion_fraction": 1.0,
        "note": (
            "Co-area formula: L = lim (1/2eps) ∫_{|u-1|<eps} |grad u| dA. "
            "Estimated via Riemann sum on a 1200x1200 grid in [-3,3]^2. "
            "SCOPING quality only; bracket from eps in {0.04, 0.02, 0.01}."
        ),
    }
    # End of co-area implementation. Old tangent-ODE code below preserved
    # for reference but unreachable.

    # Step 1: from each root, find one starting point on {|p|=1} by walking
    # outward on a generic ray (slightly skewed to avoid landing on critical
    # points exactly).
    seed_points: list[complex] = []
    for r0 in roots:
        for skew in [1.0 + 0.0j, np.exp(1j * 0.37), np.exp(1j * 1.13), np.exp(1j * 2.71)]:
            r_search = 1e-7
            prev = abs(p(r0 + r_search * skew))
            if prev >= 1.0:
                continue
            found_z = None
            for _ in range(2000):
                r_search *= 1.02
                z_test = r0 + r_search * skew
                cur = abs(p(z_test))
                if (prev - 1.0) * (cur - 1.0) < 0.0:
                    # Bisect on the segment
                    z_lo = r0 + (r_search / 1.02) * skew
                    z_hi = z_test
                    for _ in range(100):
                        z_mid = 0.5 * (z_lo + z_hi)
                        v = abs(p(z_mid)) - 1.0
                        if v > 0:
                            z_hi = z_mid
                        else:
                            z_lo = z_mid
                        if abs(z_hi - z_lo) < 1e-14:
                            break
                    found_z = 0.5 * (z_lo + z_hi)
                    break
                prev = cur
                if r_search > 100.0:
                    break
            if found_z is not None:
                seed_points.append(complex(found_z))
                break

    if not seed_points:
        return {
            "method": "tangent_following_ode",
            "L_estimate": None,
            "note": "No seed points found on {|p|=1}.",
            "completion_fraction": 0.0,
        }

    # Step 2: dedupe seed points (any two within 0.01 are on the same component
    # — we'll re-detect this when tracing returns through a previous seed).
    unique_seeds: list[complex] = []
    for s in seed_points:
        if all(abs(s - u) > 0.05 for u in unique_seeds):
            unique_seeds.append(s)

    # Step 3: from each unique seed, trace the closed curve via tangent ODE
    # dz/ds = i * p'(z) / |p'(z)|, with periodic Newton projection to keep
    # |p|=1.
    components_traced: list[dict] = []
    consumed_seeds: set[int] = set()

    for si, seed in enumerate(unique_seeds):
        if si in consumed_seeds:
            continue
        # RK4 trace around the closed curve. ds is adaptive: small fraction of
        # local curvature scale.
        z = complex(seed)
        z_start = complex(seed)
        path: list[complex] = [z]
        L_curve = 0.0
        max_steps = 200_000
        # Initial tangent direction
        ds = 1e-3
        # Compute curvature scale at start to set ds
        d1 = dp(z)
        if abs(d1) < 1e-12:
            # Critical point — skip (would need bifurcation handling)
            continue
        ds = min(0.01, 0.1 / max(1.0, abs(p.deriv(2)(z) / d1)))

        for step_i in range(max_steps):
            def tangent(zz: complex) -> complex:
                d = dp(zz)
                if abs(d) < 1e-14:
                    return 0j
                return 1j * d / abs(d)

            # RK4
            k1 = tangent(z)
            k2 = tangent(z + 0.5 * ds * k1)
            k3 = tangent(z + 0.5 * ds * k2)
            k4 = tangent(z + ds * k3)
            z_new = z + (ds / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            # Newton-project back onto {|p|=1}: solve |p(z + t*grad)|=1 for t
            # using grad |p|^2 / |grad |p|^2| direction; do up to 5 Newton steps.
            for _ in range(5):
                pz = p(z_new)
                dpz = dp(z_new)
                F = (pz * pz.conjugate()).real - 1.0
                if abs(F) < 1e-13:
                    break
                # gradient of |p|^2 in z direction is conj(p)*p' — projection
                # along its direction.
                g = pz.conjugate() * dpz
                if abs(g) < 1e-14:
                    break
                # Move along -F * g_hat / (2|g|) in real coords
                z_new -= F * g / (2.0 * abs(g) ** 2)
            L_curve += abs(z_new - z)
            z = z_new
            path.append(z)

            # Termination: have we returned to z_start? Need to have moved
            # away first (at least 5 steps) and now be within ds of start.
            if step_i > 50 and abs(z - z_start) < 1.5 * ds:
                # Add the closing segment
                L_curve += abs(z_start - z)
                break

            # Mark any seed within ds of z as consumed (multi-seed-on-same-curve)
            for sj, other_seed in enumerate(unique_seeds):
                if sj == si or sj in consumed_seeds:
                    continue
                if abs(z - other_seed) < 1.5 * ds:
                    consumed_seeds.add(sj)

        else:
            # Did not close — likely numerical drift on a noisy component
            components_traced.append({
                "seed_real": seed.real,
                "seed_imag": seed.imag,
                "n_steps": max_steps,
                "L_partial": L_curve,
                "closed": False,
            })
            continue

        consumed_seeds.add(si)
        components_traced.append({
            "seed_real": seed.real,
            "seed_imag": seed.imag,
            "n_steps": step_i + 1,
            "L_estimate": L_curve,
            "closed": True,
        })

    if not components_traced:
        return {
            "method": "tangent_following_ode",
            "L_estimate": None,
            "n_seeds": len(unique_seeds),
            "note": "All component traces failed to close.",
        }

    closed_comps = [c for c in components_traced if c.get("closed")]
    L_total = sum(c["L_estimate"] for c in closed_comps)

    return {
        "method": "tangent_following_ode_RK4",
        "n_components": len(closed_comps),
        "n_seeds_initial": len(seed_points),
        "n_seeds_unique": len(unique_seeds),
        "L_estimate": L_total,
        "components": components_traced,
        "completion_fraction": float(len(closed_comps)) / max(1, len(unique_seeds)),
        "note": (
            "Tangent-following ODE (dz/ds = i * p'(z)/|p'(z)|) with Newton "
            "projection. SCOPING quality only — not certified."
        ),
    }


def _trace_component(p: np.poly1d, dp: np.poly1d, centroid: complex,
                     grid_points: int, refine_passes: int) -> dict:
    """Trace one component of {|p|=1} via angular sweep around `centroid`.

    Returns dict with L_estimate (or None on failure) and estimates list.
    """
    thetas = np.linspace(0, 2 * np.pi, grid_points, endpoint=False)
    radii = np.zeros(grid_points)
    found = np.zeros(grid_points, dtype=bool)

    for i, theta in enumerate(thetas):
        direction = np.exp(1j * theta)
        # Find smallest r > 0 where |p(centroid + r d)| = 1, starting from r=0
        # where |p(centroid)| might be 0 (root) or non-zero.
        # We seek the FIRST crossing of |p|=1 going outward.
        r_search = 1e-6
        prev_val = abs(p(centroid + r_search * direction))
        if prev_val > 1.0:
            # Already outside contour — centroid is not enclosed; skip
            continue
        r_max = 5.0
        crossed = False
        step = 0.005
        while r_search < r_max:
            r_search += step
            cur_val = abs(p(centroid + r_search * direction))
            if (prev_val - 1.0) * (cur_val - 1.0) < 0.0:
                r_lo = r_search - step
                r_hi = r_search
                crossed = True
                break
            prev_val = cur_val
        if not crossed:
            continue
        # Bisect
        for _ in range(80):
            r_mid = 0.5 * (r_lo + r_hi)
            v = abs(p(centroid + r_mid * direction)) - 1.0
            if v > 0:
                r_hi = r_mid
            else:
                r_lo = r_mid
            if r_hi - r_lo < 1e-13:
                break
        radii[i] = 0.5 * (r_lo + r_hi)
        found[i] = True

    if not found.all():
        return {"L_estimate": None, "completion_fraction": float(found.mean())}

    dtheta = thetas[1] - thetas[0]
    dr_dtheta = np.gradient(radii, dtheta)
    integrand = np.sqrt(radii ** 2 + dr_dtheta ** 2)
    L_estimate = float(np.trapezoid(integrand, dx=dtheta))

    estimates = [L_estimate]
    # Refinement: 2× density with Newton refinement on each ray
    for _ in range(refine_passes):
        thetas_f = np.linspace(0, 2 * np.pi, grid_points * 2, endpoint=False)
        radii_f = np.interp(thetas_f, thetas, radii, period=2 * np.pi)
        for i, theta in enumerate(thetas_f):
            direction = np.exp(1j * theta)
            r0 = radii_f[i]
            for _ in range(20):
                z = centroid + r0 * direction
                pz = p(z)
                dpz = dp(z)
                F = (pz * pz.conjugate()).real - 1.0
                dFdr = 2.0 * (pz.conjugate() * dpz * direction).real
                if abs(dFdr) < 1e-14:
                    break
                step = F / dFdr
                r0 -= step
                if abs(step) < 1e-14:
                    break
            radii_f[i] = r0
        dtheta_f = thetas_f[1] - thetas_f[0]
        dr_dtheta_f = np.gradient(radii_f, dtheta_f)
        integrand_f = np.sqrt(radii_f ** 2 + dr_dtheta_f ** 2)
        L_estimate = float(np.trapezoid(integrand_f, dx=dtheta_f))
        estimates.append(L_estimate)
        thetas, radii, grid_points = thetas_f, radii_f, grid_points * 2

    return {
        "L_estimate": estimates[-1],
        "estimates": estimates,
        "completion_fraction": 1.0,
        "grid_points_final": grid_points,
        "n_refinement_passes": refine_passes,
    }


# ---------------------------------------------------------------------------
# Main scoping run
# ---------------------------------------------------------------------------

def run_n3_scoping() -> dict:
    n = 3
    L_target = L_zn_minus_1(n)

    # Boundary witness: z^3 - 1  → coeffs (a_0, a_1, a_2) = (-1, 0, 0)
    a_star = [complex(-1.0, 0.0), complex(0.0, 0.0), complex(0.0, 0.0)]
    boundary_check = cone_membership(a_star)

    # Cross-check the closed-form against numerical for sanity (NOT certifying)
    L_num_witness = lemniscate_length_numerical(a_star)

    # Three alternative polynomials. All chosen with E(p) < 1 to live strictly
    # inside K_3 (so cone says "no certificate"), giving us a clean test that
    # the cone passes them through and that L(p) does not exceed L_target.
    alternatives = [
        # (label, coeffs as (a_0, a_1, a_2))
        ("z^3 - 0.9",        [complex(-0.9, 0.0),  complex(0.0, 0.0), complex(0.0, 0.0)]),
        ("z^3 + 0.1z - 0.8", [complex(-0.8, 0.0),  complex(0.1, 0.0), complex(0.0, 0.0)]),
        ("z^3 + 0.05i z^2 - 0.9", [complex(-0.9, 0.0), complex(0.0, 0.0), complex(0.0, 0.05)]),
    ]

    alt_results = []
    for label, coeffs in alternatives:
        cone = cone_membership(coeffs)
        Lp = lemniscate_length_numerical(coeffs)
        L_est = Lp.get("L_estimate")
        within_bound = (L_est is not None) and (L_est <= L_target + 1e-6)
        alt_results.append({
            "label": label,
            "coeffs_real": [c.real for c in coeffs],
            "coeffs_imag": [c.imag for c in coeffs],
            "cone_membership": cone,
            "L_numerical": Lp,
            "L_target_zn_minus_1": L_target,
            "satisfies_L_le_L_target": within_bound,
            "L_gap_to_target": (L_target - L_est) if L_est is not None else None,
        })

    return {
        "experiment_id": "EXP-MATH-EHP114-TENSOR-CONE-V1-20260502",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "design_doc": "Math/erdos-experiments/Erdos114/TENSOR_CONE_DESIGN_2026-05-02.md",
        "scope": "SCOPING ARTIFACT — see honest scope statement in design doc",
        "n": n,
        "L_target_closed_form": L_target,
        "L_target_v5_preprint_certified_interval": [9.17972422234315, 9.17972422234317],
        "cone_definition": {
            "atomic_block_K_1": (
                "K_(1) = {(x,y) ∈ R^2 : [[1,x],[x,1]] ⪰ 0 AND [[1,y],[y,1]] ⪰ 0} "
                "= {|x|≤1, |y|≤1}"
            ),
            "tensor_structure": "K_n = K_(1)^{⊗ n} ∩ {trace ≤ 1}",
            "explicit_constraints": (
                "For each k = 0..n-1: |a_k| ≤ 1 (per-coordinate). "
                "Plus global energy: Σ |a_k|^2 ≤ 1."
            ),
            "psd_witness_z_n_minus_1": (
                "At a* = (-1,0,0,...,0): per-coord block at k=0 is rank-1 PSD "
                "([[1,-1],[-1,1]]); other blocks identity-rank-2; trace = 1 (saturated)."
            ),
        },
        "boundary_verification_z_n_minus_1": {
            "coefficient_vector": [c.real for c in a_star] + [c.imag for c in a_star],
            "cone_membership": boundary_check,
            "expected": "on_boundary_overall = True (both blocks AND trace saturated)",
            "L_numerical_for_witness": L_num_witness,
            "L_closed_form": L_target,
            "L_numerical_minus_closed_form": (
                L_num_witness["L_estimate"] - L_target if L_num_witness.get("L_estimate") else None
            ),
        },
        "alternative_polynomials": alt_results,
        "summary": {
            "boundary_check_passed": boundary_check["on_boundary_overall"],
            "n_alternatives_tested": len(alt_results),
            "n_alternatives_in_cone": sum(1 for r in alt_results if r["cone_membership"]["in_cone"]),
            "n_alternatives_satisfying_L_le_L_target": sum(1 for r in alt_results if r["satisfies_L_le_L_target"]),
            "v5_preprint_n3_margin_pct": 6.1,
            "this_cone_quantitative_margin": (
                "NOT extracted by this cone — the cone certifies coefficient-energy bound, "
                "not a quantitative L margin. See §6 of design doc for the analytic gap."
            ),
        },
        "honest_scope_statement": (
            "This is a prototype scoping artifact. The cone K_n here is candidate-defined "
            "and only verified computationally at n=3 with three alternative polynomials. "
            "To become a closure-relevant certificate, the cone definition must be (1) shown "
            "analytically to enclose all valid coefficient vectors at all n, (2) shown to "
            "have z^n-1 on its boundary at every n, (3) verifiable via interval-arithmetic "
            "SDP at every n in a tractable way. None achieved in this scoping run."
        ),
    }


def write_results(results: dict, results_dir: str) -> tuple[str, str]:
    os.makedirs(results_dir, exist_ok=True)
    base = "EXP-MATH-EHP114-TENSOR-CONE-V1-20260502"
    json_path = os.path.join(results_dir, base + "_RESULTS.json")
    if os.path.exists(json_path):
        # Per Math/CLAUDE.md Rule 3: never overwrite. Bump suffix.
        idx = 2
        while True:
            json_path = os.path.join(results_dir, f"{base}-r{idx:02d}_RESULTS.json")
            if not os.path.exists(json_path):
                break
            idx += 1
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    sha = hashlib.sha256()
    with open(json_path, "rb") as f:
        sha.update(f.read())
    sha_path = json_path[:-5] + ".sha256"
    with open(sha_path, "w") as f:
        f.write(sha.hexdigest() + "  " + os.path.basename(json_path) + "\n")
    return json_path, sha_path


def main():
    results = run_n3_scoping()
    # Print a compact summary to stdout for the running session
    s = results["summary"]
    print("=== EXP-MATH-EHP114-TENSOR-CONE-V1-20260502 — n=3 scoping ===")
    print(f"L(z^3 - 1) closed form        : {results['L_target_closed_form']:.12f}")
    print(f"L(z^3 - 1) v5 certified interval: {results['L_target_v5_preprint_certified_interval']}")
    print(f"Boundary check passed          : {s['boundary_check_passed']}")
    bv = results["boundary_verification_z_n_minus_1"]
    print(f"  block min eig (a*)           : {bv['cone_membership']['min_block_eigenvalue']:+.3e}")
    print(f"  energy trace (a*)            : {bv['cone_membership']['energy_trace']:.6f}")
    print(f"  L_num(z^3 - 1) - L_closed    : {bv['L_numerical_minus_closed_form']:+.3e}")
    print(f"# alternatives tested           : {s['n_alternatives_tested']}")
    print(f"# alternatives in cone          : {s['n_alternatives_in_cone']}")
    print(f"# alternatives w/ L ≤ L_target  : {s['n_alternatives_satisfying_L_le_L_target']}")
    for r in results["alternative_polynomials"]:
        ln = r["L_numerical"].get("L_estimate")
        gap = r["L_gap_to_target"]
        ic = r["cone_membership"]["in_cone"]
        print(f"  - {r['label']:32s}  in_cone={ic}  L≈{ln:.6f}  gap=+{gap:.4f}")

    results_dir = "/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114"
    json_path, sha_path = write_results(results, results_dir)
    print(f"\nWrote: {json_path}")
    print(f"Wrote: {sha_path}")


if __name__ == "__main__":
    main()
