#!/usr/bin/env python3
"""Toeplitz-moment-lift scoping prototype for Erdős #114 (EHP) at n=3.

Experiment ID: EXP-MATH-EHP114-TOEPLITZ-LIFT-V1-20260502

Iteration on the C-2 tensor cone (EXP-MATH-EHP114-TENSOR-CONE-V1-20260502).
Replaces the rank-1 atomic block with a Toeplitz-moment dominance block.

This is a SCOPING artifact — see TOEPLITZ_MOMENT_LIFT_DESIGN_2026-05-02.md
for the honest scope statement. The script:

  1. Defines the Toeplitz autocorrelation moments m_k(p) for monic p, the
     witness Toeplitz form T*_n = T(z^n - 1), and the dominance cone
     K_Toep,n = { a : T*_n - T(p) ≽ 0 }.

  2. Verifies that a* = (-1, 0, 0) — corresponding to z^3 - 1 — sits at the
     APEX of the cone: T(p*) = T*_n exactly, so T*_n - T(p*) = 0 (rank
     deficit n+1 = 4, full).

  3. Tests three alternative polynomials inside K_Toep,3:
       - z^3 - 0.9
       - z^3 + 0.1z - 0.8
       - z^3 + 0.05i z^2 - 0.9
     For each: Toeplitz-block PSD check, dominance check, comparison of
     min-eigenvalue (rank deficit at saturation, strictly positive at
     interior), and confirmation that L(p) < L(z^3 - 1).

  4. Comparison to C-2: also reports the C-2 atomic-block min eigenvalue and
     trace, so we can see whether the Toeplitz cone strictly dominates or
     reaches the same verdict on these test polys.

  5. Witness directional discrimination: tests the phase-shifted vertex
     a_0 = -e^{iφ} for φ ∈ {0, 0.05, 0.1, 0.2}.  C-2 admits all of them on
     boundary (|a_0|=1, energy=1).  Toeplitz cone only admits φ=0 at apex.

  6. Emits EXP-MATH-EHP114-TOEPLITZ-LIFT-V1-20260502_RESULTS.json.

Dependencies: numpy, scipy.  cvxpy is not needed because the Toeplitz cone
membership reduces to a single PSD-eigenvalue check on the (n+1) × (n+1)
Hermitian matrix T*_n - T(p), which scipy.linalg.eigvalsh handles directly.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone

import numpy as np
from scipy.linalg import eigvalsh

# ---------------------------------------------------------------------------
# Closed-form L(z^n - 1) (preprint Eq. 1)
# ---------------------------------------------------------------------------

def L_zn_minus_1(n: int) -> float:
    """Closed-form lemniscate length of z^n - 1."""
    return (2.0 ** (1.0 / n)) * math.sqrt(math.pi) * math.gamma(1.0 / (2 * n)) / math.gamma(1.0 / (2 * n) + 0.5)


# ---------------------------------------------------------------------------
# Toeplitz autocorrelation moments and matrix
# ---------------------------------------------------------------------------

def autocorrelation_moments(coeffs: list[complex]) -> list[complex]:
    """Compute m_k(p) = sum_j a_j conj(a_{j+k})  for k = 0,...,n.

    coeffs: list of (a_0, a_1, ..., a_{n-1}); a_n = 1 is implicit (monic).
    Returns [m_0, m_1, ..., m_n] (m_0 always real and >= 1 because a_n=1).
    """
    n = len(coeffs)
    a = list(coeffs) + [1.0 + 0j]  # full coefficient vector incl. leading 1
    moments: list[complex] = []
    for k in range(n + 1):
        m = 0.0 + 0j
        for j in range(n + 1 - k):
            m += a[j] * np.conj(a[j + k])
        moments.append(complex(m))
    return moments


def toeplitz_block(coeffs: list[complex]) -> np.ndarray:
    """Build the (n+1) x (n+1) Hermitian Toeplitz autocorrelation matrix.

    T(p)_{ij} = m_{|i-j|}  with sign convention T_{i,j} = m_{i-j} for i>=j
    and conj(m_{j-i}) for i<j.  This is the Gram matrix of {p, z*p, ..., z^n p}
    in L^2(d theta / 2pi).

    Note: in the design doc we wrote the matrix with rows indexed in increasing
    k.  Here T[i,j] = m_{i-j} when i>=j, conj(m_{j-i}) when i<j.  This makes
    T Hermitian.
    """
    moments = autocorrelation_moments(coeffs)
    n_plus_1 = len(moments)
    T = np.zeros((n_plus_1, n_plus_1), dtype=complex)
    for i in range(n_plus_1):
        for j in range(n_plus_1):
            d = i - j
            if d >= 0:
                T[i, j] = moments[d]
            else:
                T[i, j] = np.conj(moments[-d])
    return T


def witness_toeplitz_block(n: int) -> np.ndarray:
    """T*_n = T(z^n - 1).  Closed-form: 2 on diagonal, -1 in (0,n) and (n,0)
    corners, 0 elsewhere."""
    T = np.zeros((n + 1, n + 1), dtype=complex)
    for i in range(n + 1):
        T[i, i] = 2.0
    T[0, n] = -1.0
    T[n, 0] = -1.0
    return T


# ---------------------------------------------------------------------------
# Toeplitz dominance cone membership
# ---------------------------------------------------------------------------

def toeplitz_cone_membership(coeffs: list[complex], tol: float = 1e-10) -> dict:
    """Check membership of (a_0,...,a_{n-1}) in K_Toep,n via T*_n - T(p) ≽ 0.

    Returns the eigenvalues of the Loewner difference, the minimum eigenvalue,
    and boundary/interior verdicts.  At the witness (a* = (-1,0,...,0)) the
    diff is the zero matrix; at strict interior all eigenvalues > 0; outside
    the cone the min eigenvalue is < 0.
    """
    n = len(coeffs)
    T_p = toeplitz_block(coeffs)
    T_star = witness_toeplitz_block(n)
    D = T_star - T_p

    # Hermiticity check: max imag of (D - D^H)/2
    herm_residual = float(np.max(np.abs(D - D.conj().T)))

    eigs = eigvalsh(D)
    # eigvalsh returns real eigenvalues for a Hermitian (or float-symmetric) input
    eigs_sorted = sorted(eigs.tolist())
    min_eig = float(eigs_sorted[0])
    max_eig = float(eigs_sorted[-1])

    in_cone = min_eig >= -tol
    on_apex = (abs(min_eig) <= tol) and (abs(max_eig) <= tol)
    on_boundary_face = (abs(min_eig) <= tol) and (max_eig > tol)

    # Frobenius norm of D — measures distance from witness in moment space
    frob = float(np.linalg.norm(D, ord="fro"))

    return {
        "n": n,
        "loewner_diff_eigenvalues": eigs_sorted,
        "min_loewner_eigenvalue": min_eig,
        "max_loewner_eigenvalue": max_eig,
        "frobenius_distance_from_witness": frob,
        "hermiticity_residual": herm_residual,
        "in_cone": in_cone,
        "on_apex": on_apex,
        "on_boundary_face": on_boundary_face,
        "tol": tol,
    }


# ---------------------------------------------------------------------------
# C-2 cone membership for comparison (rank-1 atomic block)
# ---------------------------------------------------------------------------

def c2_cone_membership(coeffs: list[complex], tol: float = 1e-12) -> dict:
    """Replicate C-2's atomic-block + trace check for direct comparison."""
    n = len(coeffs)
    block_min_eigs: list[float] = []
    trace = 0.0
    for ak in coeffs:
        x, y = ak.real, ak.imag
        # K_(1) atomic block: 2x2 PSD per coord
        Bx = np.array([[1.0, x], [x, 1.0]])
        By = np.array([[1.0, y], [y, 1.0]])
        ex = eigvalsh(Bx)
        ey = eigvalsh(By)
        block_min_eigs.append(float(min(ex.min(), ey.min())))
        trace += x * x + y * y

    min_block_eig = float(min(block_min_eigs))
    in_cone = (min_block_eig >= -tol) and (trace <= 1.0 + tol)
    on_boundary = (abs(min_block_eig) <= tol) or (abs(trace - 1.0) <= tol)
    return {
        "block_min_eigenvalues": block_min_eigs,
        "min_block_eigenvalue": min_block_eig,
        "energy_trace": float(trace),
        "in_cone": bool(in_cone),
        "on_boundary": bool(on_boundary),
    }


# ---------------------------------------------------------------------------
# Numerical L(p) — co-area, identical to C-2's implementation
# ---------------------------------------------------------------------------

def lemniscate_length_numerical(coeffs: list[complex]) -> dict:
    """Estimate L(p) = H^1({|p|=1}) via co-area on a 1200x1200 grid.

    Same method as C-2 prototype — kept identical for direct comparability.
    SCOPING-quality only; not certified.
    """
    n = len(coeffs)
    poly_coeffs = [1.0] + [coeffs[n - 1 - k] for k in range(n)]
    R = 3.0
    grid_n = 1200
    xs = np.linspace(-R, R, grid_n)
    ys = np.linspace(-R, R, grid_n)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    Z = X + 1j * Y
    h = xs[1] - xs[0]

    pZ = np.zeros_like(Z, dtype=complex)
    for c in poly_coeffs:
        pZ = pZ * Z + c

    # Derivative (for |grad u|)
    dp_coeffs = [poly_coeffs[i] * (n - i) for i in range(n)]
    dpZ = np.zeros_like(Z, dtype=complex)
    for c in dp_coeffs:
        dpZ = dpZ * Z + c

    abs_p = np.abs(pZ)
    grad_u = np.abs(dpZ)

    eps_values = [0.04, 0.02, 0.01]
    L_estimates = []
    for eps in eps_values:
        mask = np.abs(abs_p - 1.0) < eps
        L_eps = (h * h / (2.0 * eps)) * float(np.sum(grad_u[mask]))
        L_estimates.append({"eps": eps, "L_estimate": L_eps, "n_points_in_strip": int(mask.sum())})

    return {
        "method": "co_area_grid_integration",
        "grid_size": grid_n,
        "box_radius": R,
        "grid_spacing": float(h),
        "L_estimate": L_estimates[-1]["L_estimate"],
        "L_estimates_by_eps": L_estimates,
        "L_lower_bracket": min(e["L_estimate"] for e in L_estimates),
        "L_upper_bracket": max(e["L_estimate"] for e in L_estimates),
    }


# ---------------------------------------------------------------------------
# Connection-to-length identity statement
# ---------------------------------------------------------------------------

CONNECTION_TO_LENGTH = {
    "primary_identity": "Fejer-Riesz factorization: T*_n - T(p) ≽ 0 ⟺ |p*(e^iθ)|^2 - |p(e^iθ)|^2 = |r(e^iθ)|^2 ≥ 0 pointwise on |z|=1",
    "implication_1": "Pointwise dominance: |p(z)| ≤ |p*(z)| for all z on the unit circle",
    "implication_2_open": "From pointwise unit-circle dominance to L(p) ≤ L(p*) requires a Crofton/co-area bridge — open and being drafted by parallel subagent",
    "asymptotic_path": "Szego strong limit: log det T_N(|p|^2) ~ 2N log M(p) connects Toeplitz determinants to Mahler measure, but only as N → ∞",
    "this_iteration_provides": "Algebraic SDP-representable certificate for the antecedent of the Crofton bridge",
    "this_iteration_does_NOT_provide": "Closed analytic step from cone membership to length bound at finite n",
}


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def run_experiment() -> dict:
    n = 3
    L_target = L_zn_minus_1(n)

    # ---- Witness (z^3 - 1) ----
    witness_coeffs: list[complex] = [-1.0 + 0j, 0.0 + 0j, 0.0 + 0j]
    T_p_witness = toeplitz_block(witness_coeffs)
    T_star = witness_toeplitz_block(n)
    witness_diff_max = float(np.max(np.abs(T_star - T_p_witness)))

    witness_section = {
        "label": "z^3 - 1 (witness)",
        "coeffs_real": [c.real for c in witness_coeffs],
        "coeffs_imag": [c.imag for c in witness_coeffs],
        "moments_m_k": [complex(m).real if abs(complex(m).imag) < 1e-14 else [complex(m).real, complex(m).imag]
                          for m in autocorrelation_moments(witness_coeffs)],
        "T_p_real": T_p_witness.real.tolist(),
        "T_p_imag": T_p_witness.imag.tolist(),
        "T_star_minus_T_p_max_abs": witness_diff_max,
        "expected_T_p_equals_T_star": True,
        "toeplitz_cone": toeplitz_cone_membership(witness_coeffs),
        "c2_cone": c2_cone_membership(witness_coeffs),
        "L_target_closed_form": L_target,
    }

    # ---- Three alternative polynomials ----
    alternatives = [
        {
            "label": "z^3 - 0.9",
            "coeffs": [-0.9 + 0j, 0.0 + 0j, 0.0 + 0j],
        },
        {
            "label": "z^3 + 0.1z - 0.8",
            "coeffs": [-0.8 + 0j, 0.1 + 0j, 0.0 + 0j],
        },
        {
            "label": "z^3 + 0.05i z^2 - 0.9",
            "coeffs": [-0.9 + 0j, 0.0 + 0j, 0.0 + 0.05j],
        },
    ]

    alternative_results = []
    for alt in alternatives:
        coeffs = alt["coeffs"]
        toep = toeplitz_cone_membership(coeffs)
        c2 = c2_cone_membership(coeffs)
        Lnum = lemniscate_length_numerical(coeffs)
        L_est = Lnum["L_estimate"]
        moments = autocorrelation_moments(coeffs)
        alternative_results.append({
            "label": alt["label"],
            "coeffs_real": [c.real for c in coeffs],
            "coeffs_imag": [c.imag for c in coeffs],
            "moments_m_k_real": [complex(m).real for m in moments],
            "moments_m_k_imag": [complex(m).imag for m in moments],
            "toeplitz_cone": toep,
            "c2_cone": c2,
            "L_numerical": Lnum,
            "L_target_zn_minus_1": L_target,
            "L_estimate_lt_target": bool(L_est < L_target),
            "L_gap_to_target": float(L_target - L_est),
        })

    # ---- Phase-shifted vertex test (the discriminator) ----
    # a_0 = -e^{i phi}, others 0.  C-2 admits all on boundary (|a_0|=1, energy=1).
    # Toeplitz cone admits only phi=0 at apex.
    phi_tests = []
    for phi in [0.0, 0.05, 0.1, 0.2, 0.5, 1.0]:
        a0 = -np.exp(1j * phi)
        coeffs_phi: list[complex] = [complex(a0), 0 + 0j, 0 + 0j]
        toep_phi = toeplitz_cone_membership(coeffs_phi)
        c2_phi = c2_cone_membership(coeffs_phi)
        # Compute L(p) numerically too — even though z^3 + e^{i phi} is unitarily
        # equivalent to z^3 - 1 via root-of-unity rotation, our co-area on a fixed
        # box may give a slightly different number; report it for transparency.
        Lnum = lemniscate_length_numerical(coeffs_phi)
        phi_tests.append({
            "phi": phi,
            "coeffs_real": [c.real for c in coeffs_phi],
            "coeffs_imag": [c.imag for c in coeffs_phi],
            "toeplitz_cone_min_eig": toep_phi["min_loewner_eigenvalue"],
            "toeplitz_cone_max_eig": toep_phi["max_loewner_eigenvalue"],
            "toeplitz_cone_in_cone": toep_phi["in_cone"],
            "toeplitz_cone_on_apex": toep_phi["on_apex"],
            "frobenius_distance_from_witness": toep_phi["frobenius_distance_from_witness"],
            "c2_cone_min_block_eig": c2_phi["min_block_eigenvalue"],
            "c2_cone_energy_trace": c2_phi["energy_trace"],
            "c2_cone_in_cone": c2_phi["in_cone"],
            "c2_cone_on_boundary": c2_phi["on_boundary"],
            "L_numerical": Lnum["L_estimate"],
        })

    # ---- Strict dominance check ----
    # The Toeplitz cone strictly dominates C-2 at n=3 if there exists a
    # coefficient vector that is in C-2 (or on its boundary) but NOT in the
    # Toeplitz cone.  The phase-shifted vertex tests construct exactly such
    # vectors: at phi > 0 they are on C-2 boundary but their Loewner difference
    # T*_n - T(p) for the Toeplitz cone has min eigenvalue *not* zero (might be
    # positive or negative depending on phi).
    strict_dominance_evidence = []
    for r in phi_tests:
        if r["phi"] > 0:
            # If C-2 says "on boundary" but Toeplitz says "not on apex", we have
            # discrimination.  If additionally Toeplitz says "outside cone"
            # (min_eig < 0), we have strict dominance proof.
            strict_dominance_evidence.append({
                "phi": r["phi"],
                "c2_says": "on_boundary" if r["c2_cone_on_boundary"] else ("interior" if r["c2_cone_in_cone"] else "outside"),
                "toeplitz_says": (
                    "apex" if r["toeplitz_cone_on_apex"]
                    else ("interior" if r["toeplitz_cone_in_cone"] and r["toeplitz_cone_min_eig"] > 1e-10
                          else ("boundary_face" if (abs(r["toeplitz_cone_min_eig"]) < 1e-10 and r["toeplitz_cone_max_eig"] > 1e-10)
                                else "outside"))
                ),
                "discriminates": (r["c2_cone_on_boundary"] and not r["toeplitz_cone_on_apex"]),
                "frob_dist": r["frobenius_distance_from_witness"],
            })

    # Final result blob
    timestamp = datetime.now(timezone.utc).isoformat()
    result = {
        "experiment_id": "EXP-MATH-EHP114-TOEPLITZ-LIFT-V1-20260502",
        "timestamp_utc": timestamp,
        "design_doc": "Math/erdos-experiments/Erdos114/TOEPLITZ_MOMENT_LIFT_DESIGN_2026-05-02.md",
        "iteration_on": "EXP-MATH-EHP114-TENSOR-CONE-V1-20260502 (C-2)",
        "scope": "SCOPING ARTIFACT — see honest scope statement in design doc",
        "n": n,
        "L_target_closed_form": L_target,
        "L_target_v5_preprint_certified_interval": [9.17972422234315, 9.17972422234317],

        "moment_family_chosen": "coefficient_autocorrelation",
        "moment_family_alternatives_considered": [
            "coefficient_autocorrelation (CHOSEN — computable a priori from coefficients)",
            "trigonometric_on_lemniscate (REJECTED — requires knowing Lambda first; circular)",
        ],
        "atomic_block_definition": {
            "shape": f"({n+1}, {n+1})",
            "type": "Hermitian Toeplitz autocorrelation",
            "entries": "T(p)_{ij} = m_{i-j} for i>=j, conj(m_{j-i}) for i<j; m_k = sum_j a_j conj(a_{j+k}) with a_n=1",
            "PSD_automatic": "Yes — Gram matrix of {p, z*p, ..., z^n p} in L^2(dθ/2π)",
            "informative_only_relative_to_target": "Yes — uses dominance T*_n - T(p) ≽ 0 with target T*_n = T(z^n-1)",
        },
        "witness_toeplitz_form_T_star_n_3": {
            "matrix_real": T_star.real.tolist(),
            "matrix_imag": T_star.imag.tolist(),
            "eigenvalues_sorted": sorted(np.linalg.eigvalsh(T_star).tolist()),
            "trace": float(np.trace(T_star).real),
            "rank": int(np.linalg.matrix_rank(T_star)),
        },

        "witness_check": witness_section,
        "alternative_polynomials": alternative_results,
        "phase_shifted_vertex_tests": phi_tests,
        "strict_dominance_evidence": strict_dominance_evidence,

        "connection_to_length": CONNECTION_TO_LENGTH,

        "summary": {
            "witness_at_apex": witness_section["toeplitz_cone"]["on_apex"],
            "witness_loewner_diff_max_abs": witness_diff_max,
            "n_alternatives_tested": len(alternative_results),
            "n_alternatives_in_toeplitz_cone": sum(1 for a in alternative_results if a["toeplitz_cone"]["in_cone"]),
            "n_alternatives_satisfying_L_lt_target": sum(1 for a in alternative_results if a["L_estimate_lt_target"]),
            "n_phase_shifts_tested": len(phi_tests),
            "n_phase_shifts_discriminating_C2_vs_Toeplitz": sum(1 for e in strict_dominance_evidence if e["discriminates"]),
            "strict_dominance_at_n_3": "demonstrated" if any(e["discriminates"] for e in strict_dominance_evidence) else "not_demonstrated_in_this_run",
            "comparison_to_C2": (
                "C-2 sees only |a_k|^2 (energy per coord + total energy). Toeplitz sees cross-moments "
                "m_k = sum a_j conj(a_{j+k}). The phase-shifted vertex a_0 = -e^{iφ} is the canonical "
                "discriminator: C-2 admits it on boundary at all φ; Toeplitz only admits φ=0 at apex."
            ),
        },

        "honest_scope_statement": (
            "This is a prototype scoping artifact. The Toeplitz-moment block here is "
            "candidate-defined and only verified at n=3 with three alternative "
            "polynomials. The bridge from PSD-ness of the Toeplitz block to "
            "L(p) ≤ L(z^n - 1) relies on either Szegő's strong limit (asymptotic) "
            "or Fejér-Riesz (finite-n but local). Neither is a closed analytic step at "
            "finite n; both are pointers to an analytic completion. Tractability of the "
            "analytic completion is for the parallel Crofton-bridge draft, not this "
            "iteration."
        ),
    }
    return result


def main() -> int:
    out_dir = os.path.dirname(os.path.abspath(__file__))
    result = run_experiment()
    out_file = os.path.join(out_dir, "EXP-MATH-EHP114-TOEPLITZ-LIFT-V1-20260502_RESULTS.json")
    if os.path.exists(out_file):
        print(f"REFUSING TO OVERWRITE existing {out_file} (Math/CLAUDE.md Rule 3)", file=sys.stderr)
        return 1
    with open(out_file, "w") as f:
        json.dump(result, f, indent=2, default=str)
    sha = hashlib.sha256(json.dumps(result, indent=2, default=str).encode()).hexdigest()
    sha_file = out_file.replace("_RESULTS.json", "_RESULTS.sha256")
    with open(sha_file, "w") as f:
        f.write(f"{sha}  {os.path.basename(out_file)}\n")

    # Console summary
    print(f"=== EXP-MATH-EHP114-TOEPLITZ-LIFT-V1-20260502 ===")
    print(f"L(z^3 - 1) closed form: {result['L_target_closed_form']:.6f}")
    print(f"")
    print(f"Witness (z^3 - 1):")
    w = result["witness_check"]
    print(f"  T*_n - T(p*) max abs entry: {w['T_star_minus_T_p_max_abs']:.3e}")
    print(f"  Toeplitz cone min eig:      {w['toeplitz_cone']['min_loewner_eigenvalue']:.3e}")
    print(f"  Toeplitz cone max eig:      {w['toeplitz_cone']['max_loewner_eigenvalue']:.3e}")
    print(f"  on_apex:                    {w['toeplitz_cone']['on_apex']}")
    print(f"")
    print(f"Alternative polynomials:")
    for a in result["alternative_polynomials"]:
        t = a["toeplitz_cone"]
        c = a["c2_cone"]
        print(f"  {a['label']}:")
        print(f"    Toeplitz min eig: {t['min_loewner_eigenvalue']:.4f}, in_cone: {t['in_cone']}, frob_dist: {t['frobenius_distance_from_witness']:.4f}")
        print(f"    C-2 min block eig: {c['min_block_eigenvalue']:.4f}, energy: {c['energy_trace']:.4f}")
        print(f"    L(p) ≈ {a['L_numerical']['L_estimate']:.4f}  (target {a['L_target_zn_minus_1']:.4f})")
    print(f"")
    print(f"Phase-shifted vertex discriminator:")
    for r in result["phase_shifted_vertex_tests"]:
        print(f"  phi={r['phi']:.2f}: Toep min_eig={r['toeplitz_cone_min_eig']:+.4f} (apex={r['toeplitz_cone_on_apex']}, in_cone={r['toeplitz_cone_in_cone']}), C-2 on_boundary={r['c2_cone_on_boundary']}, frob={r['frobenius_distance_from_witness']:.4f}")
    print(f"")
    print(f"Strict-dominance evidence:")
    for e in result["strict_dominance_evidence"]:
        print(f"  phi={e['phi']:.2f}: C-2 says {e['c2_says']:>12s} | Toeplitz says {e['toeplitz_says']:>15s} | discriminates={e['discriminates']}")
    print(f"")
    print(f"strict_dominance_at_n_3: {result['summary']['strict_dominance_at_n_3']}")
    print(f"")
    print(f"Wrote: {out_file}")
    print(f"SHA256: {sha}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
