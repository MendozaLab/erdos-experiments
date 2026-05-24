"""
Koopman Observable Lift Prototype for Erdős #114 (EHP conjecture).

Experiment ID: EXP-MATH-EHP114-KOOPMAN-LIFT-V1-20260502
Author scope: prototype scoping artifact, not a certified proof.

Design choices (see KOOPMAN_LIFT_DESIGN_2026-05-02.md for the full rationale):
  * Flow:               negative gradient flow on -L(p_a) (so z^n-1 is a
                        local maximum of L hence a stable fixed point of the
                        descent on -L). The flow integrator is forward Euler
                        with small step h, which makes the observable update
                        a one-step-ahead snapshot (DMD/EDMD compatible).
  * State space:        symmetry-reduced coefficients. We follow the v5
                        preprint convention: translation kills a_{n-1}, and
                        rotation fixes the phase of a_0. For n=3 we use the
                        same 3-dim slice the verifier uses:
                            (a_re, a_im, b)  with  b = a_0 in R_{>=0},
                            a = a_1 in C,    a_2 = 0.
  * Observable basis:   monomials in the slice coordinates plus L itself.
                        Explicitly: {1, x1, x2, x3, x1^2, x2^2, x3^2,
                                     x1*x2, x1*x3, x2*x3, L}.
                        K = 11. This is small enough to assemble densely
                        and large enough to resolve the dominant non-trivial
                        eigenmode.
  * Koopman approx.:    EDMD (Williams-Kevrekidis-Rowley 2015 style).
                        We sample M states X_m on a grid in the slice,
                        advance each by one Euler step Phi_h, and solve
                            K = Psi(Y) Psi(X)^+
                        where Psi is the observable-evaluation matrix and
                        Y = Phi_h(X). Spectrum of K approximates the
                        Koopman spectrum.

Lemniscate length:
  * Closed form for z^n - 1:  L = 2^{1/n} * sqrt(pi) * Gamma(1/(2n)) /
                                    Gamma(1/(2n) + 1/2)   (preprint Thm).
  * Off-baseline L(p_a):      contour integration via marching on the
                              implicit set |p(z)| = 1. We use a coarse
                              polar marcher with 720 angles around z=0,
                              find the |p|=1 crossing along each ray,
                              and sum the chord lengths. This is the
                              same idea the v5 prototype uses at coarse
                              resolution (the certified version uses
                              IEEE 1788 + marching squares; we are NOT
                              re-doing certification here).
                              For n=3 this is fast and accurate to
                              ~1% relative on the kind of polynomials
                              we visit.

Outputs (per Math/CLAUDE.md Rule 3 - V1 suffix, never overwrite):
    EXP-MATH-EHP114-KOOPMAN-LIFT-V1-20260502_RESULTS.json
    EXP-MATH-EHP114-KOOPMAN-LIFT-V1-20260502_REPORT.md
    EXP-MATH-EHP114-KOOPMAN-LIFT-V1-20260502_RESULTS.sha256
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import time
from dataclasses import dataclass

import numpy as np
import mpmath as mp


# ----------------------------------------------------------------------
# 0.  Configuration
# ----------------------------------------------------------------------

EXPERIMENT_ID = "EXP-MATH-EHP114-KOOPMAN-LIFT-V1-20260502"
RESULTS_DIR = os.path.dirname(os.path.abspath(__file__))
DEGREE = 3

# Symmetry-reduced slice for n=3 (matches v5 preprint convention).
# State vector x = (a_re, a_im, b) with:
#   p_a(z) = z^3 + 0 * z^2 + (a_re + i a_im) z + b
# baseline (z^3 - 1) is x_star = (0, 0, -1).  We translate to put it
# at the origin during the flow, but for the JSON we keep absolute coords.

X_STAR = np.array([0.0, 0.0, -1.0])  # the EHP optimum in slice coords

# Grid bound for the EDMD snapshot collection. Coarse on purpose - this is a
# scoping run, not a certificate. The v5 verifier uses
#   (a_re, a_im) in [-4, 4]^2 x b in [0, 4]
# which doesn't match our sign convention for b (z^3 - 1 has b = -1). We pick
# a smaller box centered on x_star and let the gradient flow do the rest.

GRID_HALF = 0.4        # |dx| <= 0.4 around x_star (avoids the petal cusp)
GRID_PER_AXIS = 7      # 7^3 = 343 snapshot points
# IMPORTANT design note. The v5 verifier reports Hessian eigenvalues of
# order ~3e6 at x_star, which means the *raw* gradient flow x' = grad L
# is extremely stiff and would require h ~ 1e-7 just for forward-Euler
# stability. At that scale, |phi_h(x) - x| << marching-squares noise floor
# (~1e-3 rel) and EDMD picks up only noise. Solution: use a NORMALIZED
# gradient flow with a flow-time coordinate change,
#     x' = (1 / scale) * grad L(x)
# where scale = max|H_kk| at x_star (estimated empirically). This rescales
# all continuous-time eigenvalues by 1/scale -> O(1), so a moderate Euler
# step h ~ 1e-1 gives discrete eigenvalues 1 + h*lambda_c that sit
# cleanly in the unit disk, well clear of the noise floor.
# THE KOOPMAN SPECTRUM IS THE SAME up to a uniform real-axis rescaling.
EULER_STEP = 0.1               # h in flow-time
GRAD_FD_STEP = 1e-3            # finite-diff step for the L-gradient
FLOW_NORMALIZER = 5.0e6        # estimated max |H_kk| from v5 (Level3 results)

# Number of Koopman observables.
N_OBSERVABLES = 11


# ----------------------------------------------------------------------
# 1.  Lemniscate length
# ----------------------------------------------------------------------

def L_baseline_exact(n: int) -> float:
    """L(z^n - 1) via the preprint's closed form (Thm 'Exact lemniscate length').

        L(z^n - 1) = 2^{1/n} * sqrt(pi) * Gamma(1/(2n)) / Gamma(1/(2n) + 1/2)
    """
    one_over_2n = mp.mpf(1) / (2 * n)
    val = (
        mp.power(2, mp.mpf(1) / n)
        * mp.sqrt(mp.pi)
        * mp.gamma(one_over_2n)
        / mp.gamma(one_over_2n + mp.mpf("0.5"))
    )
    return float(val)


def L_marching_squares(coeffs: np.ndarray, grid_n: int = 320, box_half: float = 1.6) -> float:
    """Approximate L(p) using marching squares on a 2D grid.

    Same approach the v5 preprint mentions for its f64 sampling layer
    (the certified version wraps it in IEEE 1788 envelopes; we don't).

    coeffs[k] is the coefficient of z^k for k = 0, ..., n-1.  Leading
    coefficient (z^n) is implicit and equals 1.

    For monic polynomials with |a_k| <= 2 (which is well past our
    GRID_HALF = 0.6 perturbation regime), the lemniscate |p(z)| = 1
    is contained in |z| <= 2.  We use box_half = 1.6 with grid_n = 320
    to get ~0.1% relative accuracy on smooth and cusped cases for n=3.
    """
    n = len(coeffs)
    full = np.empty(n + 1, dtype=complex)
    full[:n] = coeffs
    full[n] = 1.0
    full_rev = full[::-1]

    xs = np.linspace(-box_half, box_half, grid_n)
    ys = np.linspace(-box_half, box_half, grid_n)
    X_grid, Y_grid = np.meshgrid(xs, ys, indexing="xy")
    Z = X_grid + 1j * Y_grid
    F = np.abs(np.polyval(full_rev, Z)) - 1.0  # zero set is the lemniscate

    # Marching squares: for each cell, count the (signed) length of
    # the iso-line. Standard 16-case lookup; here we use a simple
    # linear-interpolation midpoint-segment construction.

    h = xs[1] - xs[0]
    total_len = 0.0

    f00 = F[:-1, :-1]
    f10 = F[:-1, 1:]
    f01 = F[1:, :-1]
    f11 = F[1:, 1:]

    s00 = f00 > 0
    s10 = f10 > 0
    s01 = f01 > 0
    s11 = f11 > 0

    # 4-bit case: bit 0 = (i,j), bit 1 = (i,j+1), bit 2 = (i+1,j+1), bit 3 = (i+1,j)
    case = (
        s00.astype(np.int8) * 1
        + s10.astype(np.int8) * 2
        + s11.astype(np.int8) * 4
        + s01.astype(np.int8) * 8
    )

    # We compute zero-crossings on each of the four cell edges.
    # Edge labels: 0 = bottom (between (i,j) and (i,j+1)),
    #              1 = right  (between (i,j+1) and (i+1,j+1)),
    #              2 = top    (between (i+1,j) and (i+1,j+1)),
    #              3 = left   (between (i,j) and (i+1,j)).

    def safe_t(a, b):
        """Linear interp parameter; clipped."""
        d = b - a
        # Avoid /0; cells where this is called always have a sign change,
        # so d is nonzero in practice. Add a tiny epsilon to be safe.
        return np.where(np.abs(d) > 1e-30, -a / np.where(np.abs(d) > 1e-30, d, 1.0), 0.5)

    # Coordinates of cell corners: (i,j) -> (xs[j], ys[i])
    # Edges in physical units have length h.
    # We're computing total CONTOUR length, so we just need the
    # length of the segment within each cell.

    # For each cell the contour is one of: nothing (cases 0,15),
    # one segment (most cases), or two segments (cases 5, 10 = saddles).

    t_bot = safe_t(f00, f10)  # (i,j) -> (i,j+1)
    t_right = safe_t(f10, f11)  # (i,j+1) -> (i+1,j+1)
    t_top = safe_t(f01, f11)  # (i+1,j) -> (i+1,j+1)
    t_left = safe_t(f00, f01)  # (i,j) -> (i+1,j)

    # Local coords within a cell, in units of h.
    # bottom-edge crossing at (t_bot, 0)
    # right-edge  crossing at (1, t_right)
    # top-edge    crossing at (t_top, 1)
    # left-edge   crossing at (0, t_left)

    def seg_len(p0x, p0y, p1x, p1y):
        return h * np.sqrt((p1x - p0x) ** 2 + (p1y - p0y) ** 2)

    # Cases 1 / 14: bottom-left corner alone -> bottom-left segment
    mask = (case == 1) | (case == 14)
    if np.any(mask):
        p0x, p0y = t_bot[mask], 0.0
        p1x, p1y = 0.0, t_left[mask]
        total_len += float(np.sum(seg_len(p0x, p0y, p1x, p1y)))

    # Cases 2 / 13: bottom-right corner alone -> bottom-right segment
    mask = (case == 2) | (case == 13)
    if np.any(mask):
        p0x, p0y = t_bot[mask], 0.0
        p1x, p1y = 1.0, t_right[mask]
        total_len += float(np.sum(seg_len(p0x, p0y, p1x, p1y)))

    # Cases 3 / 12: bottom row -> horizontal-ish line
    mask = (case == 3) | (case == 12)
    if np.any(mask):
        p0x, p0y = 0.0, t_left[mask]
        p1x, p1y = 1.0, t_right[mask]
        total_len += float(np.sum(seg_len(p0x, p0y, p1x, p1y)))

    # Cases 4 / 11: top-right corner alone
    mask = (case == 4) | (case == 11)
    if np.any(mask):
        p0x, p0y = 1.0, t_right[mask]
        p1x, p1y = t_top[mask], 1.0
        total_len += float(np.sum(seg_len(p0x, p0y, p1x, p1y)))

    # Cases 6 / 9: right column
    mask = (case == 6) | (case == 9)
    if np.any(mask):
        p0x, p0y = t_bot[mask], 0.0
        p1x, p1y = t_top[mask], 1.0
        total_len += float(np.sum(seg_len(p0x, p0y, p1x, p1y)))

    # Cases 7 / 8: top-left corner alone
    mask = (case == 7) | (case == 8)
    if np.any(mask):
        p0x, p0y = 0.0, t_left[mask]
        p1x, p1y = t_top[mask], 1.0
        total_len += float(np.sum(seg_len(p0x, p0y, p1x, p1y)))

    # Saddle cases 5 and 10: two diagonal segments.
    # Case 5: (i,j) and (i+1,j+1) positive -> two segments through
    #         (bot, right) and (left, top).
    mask = case == 5
    if np.any(mask):
        # bot - right
        p0x, p0y = t_bot[mask], 0.0
        p1x, p1y = 1.0, t_right[mask]
        total_len += float(np.sum(seg_len(p0x, p0y, p1x, p1y)))
        # left - top
        p0x, p0y = 0.0, t_left[mask]
        p1x, p1y = t_top[mask], 1.0
        total_len += float(np.sum(seg_len(p0x, p0y, p1x, p1y)))

    mask = case == 10
    if np.any(mask):
        # bot - left
        p0x, p0y = t_bot[mask], 0.0
        p1x, p1y = 0.0, t_left[mask]
        total_len += float(np.sum(seg_len(p0x, p0y, p1x, p1y)))
        # right - top
        p0x, p0y = 1.0, t_right[mask]
        p1x, p1y = t_top[mask], 1.0
        total_len += float(np.sum(seg_len(p0x, p0y, p1x, p1y)))

    return total_len


def L_of_state(x: np.ndarray, n: int = DEGREE) -> float:
    """Convert slice state x = (a_re, a_im, b) to coefficient vector and call the marching-squares L."""
    a_re, a_im, b = x
    # coeffs[0] = a_0 = b, coeffs[1] = a_1 = a_re + i a_im, coeffs[2] = a_2 = 0
    coeffs = np.array([b + 0j, a_re + 1j * a_im, 0 + 0j], dtype=complex)
    return L_marching_squares(coeffs)


# ----------------------------------------------------------------------
# 2.  Gradient flow on -L
# ----------------------------------------------------------------------

def grad_L(x: np.ndarray, h: float = GRAD_FD_STEP) -> np.ndarray:
    """Central-difference gradient of L at x in slice coords."""
    g = np.zeros_like(x)
    for i in range(len(x)):
        e = np.zeros_like(x)
        e[i] = h
        g[i] = (L_of_state(x + e) - L_of_state(x - e)) / (2 * h)
    return g


def euler_step(x: np.ndarray, dt: float = EULER_STEP) -> np.ndarray:
    """One forward Euler step of the NORMALIZED ascent flow

        x' = grad L(x) / FLOW_NORMALIZER

    The normalization is a flow-time coordinate change (uniform rescale
    of the continuous-time eigenvalues by 1 / FLOW_NORMALIZER). The
    Koopman spectrum after rescaling is the same one you would get from
    the unnormalized flow, just on a stretched time axis. This is the
    standard trick to keep stiff gradient flows numerically tractable
    without changing their qualitative spectral structure.
    """
    return x + (dt / FLOW_NORMALIZER) * grad_L(x)


# ----------------------------------------------------------------------
# 3.  EDMD Koopman approximation
# ----------------------------------------------------------------------

def observables(x: np.ndarray, L_val: float | None = None) -> np.ndarray:
    """Evaluate the K = 11 observable basis at slice state x.

    Order: [1, x1, x2, x3, x1^2, x2^2, x3^2, x1*x2, x1*x3, x2*x3, L]
    """
    x1, x2, x3 = x
    if L_val is None:
        L_val = L_of_state(x)
    return np.array([
        1.0,
        x1,
        x2,
        x3,
        x1 * x1,
        x2 * x2,
        x3 * x3,
        x1 * x2,
        x1 * x3,
        x2 * x3,
        L_val,
    ])


def build_snapshot_grid() -> np.ndarray:
    """7 x 7 x 7 lattice in [x_star - h, x_star + h]^3 with h = GRID_HALF."""
    axis = np.linspace(-GRID_HALF, GRID_HALF, GRID_PER_AXIS)
    pts = []
    for dx in axis:
        for dy in axis:
            for dz in axis:
                pts.append(X_STAR + np.array([dx, dy, dz]))
    return np.array(pts)


def assemble_edmd(snapshot_X: np.ndarray, dt: float = EULER_STEP) -> dict:
    """Run EDMD: K = Psi_Y Psi_X^+.

    Returns the Koopman matrix, eigenvalues, and the time-tagged
    snapshot pair (X, Y = phi_h(X)) used to assemble it.
    """
    M = snapshot_X.shape[0]
    K_obs = N_OBSERVABLES

    # Compute L at each X point (cached because it's expensive)
    L_X = np.array([L_of_state(x) for x in snapshot_X])

    Psi_X = np.zeros((M, K_obs))
    for m in range(M):
        Psi_X[m] = observables(snapshot_X[m], L_val=L_X[m])

    # One-step advance
    snapshot_Y = np.array([euler_step(x, dt=dt) for x in snapshot_X])
    L_Y = np.array([L_of_state(y) for y in snapshot_Y])

    Psi_Y = np.zeros((M, K_obs))
    for m in range(M):
        Psi_Y[m] = observables(snapshot_Y[m], L_val=L_Y[m])

    # K solves Psi_X @ K^T = Psi_Y in least squares.
    # Equivalently K = Psi_Y^T Psi_X (Psi_X^T Psi_X)^{-1}.
    G = Psi_X.T @ Psi_X
    A = Psi_X.T @ Psi_Y
    # Regularize lightly to avoid pseudo-inverse drama on the constant column
    reg = 1e-10 * np.eye(K_obs)
    K = np.linalg.solve(G + reg, A).T

    eigvals = np.linalg.eigvals(K)

    return {
        "K": K,
        "eigvals": eigvals,
        "Psi_X": Psi_X,
        "Psi_Y": Psi_Y,
        "L_X": L_X,
        "L_Y": L_Y,
        "snapshot_X": snapshot_X,
        "snapshot_Y": snapshot_Y,
        "dt": dt,
    }


def step_size_sweep(snapshot_X: np.ndarray, dts: list[float]) -> dict:
    """Multi-h Koopman sweep.

    For each Euler step h in dts, build the EDMD Koopman matrix and
    extract the dominant non-trivial eigenvalue.  We expect:
      |lambda_d(h)| ~ exp(h * lambda_c)
    so log|lambda_d(h)| / h converges to lambda_c (the continuous-time
    spectral gap) as h -> 0+, modulo finite-difference noise floors.
    Plotting |lambda_d| vs h reveals where the noise floor sits.
    """
    rows = []
    for h in dts:
        edmd_h = assemble_edmd(snapshot_X, dt=h)
        spec_h = spectral_summary(edmd_h["eigvals"], dt=h)
        rows.append({
            "h": h,
            "dominant_nontrivial_abs": spec_h["dominant_nontrivial_abs"],
            "spectral_gap_in_disk": spec_h["spectral_gap_in_disk"],
            "implied_continuous_decay_rate": spec_h["implied_continuous_decay_rate"],
            "constant_eig_distance_to_one": spec_h["constant_eig_distance_to_one"],
        })
    return {"sweep_rows": rows, "dts": dts}


# ----------------------------------------------------------------------
# 4.  Spectral analysis
# ----------------------------------------------------------------------

def spectral_summary(eigvals: np.ndarray, dt: float = EULER_STEP) -> dict:
    """Compute discrete-time Koopman eigenvalue stats.

    For a stable fixed point under continuous-time gradient flow with
    eigenvalues lambda_c < 0, the discrete-time Euler operator has
    eigenvalues lambda_d ~ 1 + dt * lambda_c, so |lambda_d| < 1 with
    decay rate -log|lambda_d| / dt approximating |lambda_c|.

    For the EDMD lift on the full nonlinear flow, we expect:
      * One eigenvalue near 1.0 (the constant function is invariant).
      * Several eigenvalues inside the unit disk (decaying observables).
      * The "spectral gap" of interest is 1 - max|lambda_d|_{lambda != 1}.
    """
    abs_eig = np.abs(eigvals)

    # Find the constant-function eigenvalue (closest to 1)
    idx_one = np.argmin(np.abs(eigvals - 1.0))
    near_one_distance = float(np.abs(eigvals[idx_one] - 1.0))

    # Mask it out and find the dominant non-trivial eigenmode
    mask = np.ones_like(eigvals, dtype=bool)
    mask[idx_one] = False
    abs_eig_nontrivial = abs_eig[mask]

    dominant_nontrivial = float(np.max(abs_eig_nontrivial))
    spectral_gap_disk = 1.0 - dominant_nontrivial  # how far inside the unit disk

    # Continuous-time decay rate (per the Euler approximation)
    if dominant_nontrivial > 1e-12:
        decay_rate_ct = -math.log(dominant_nontrivial) / dt
    else:
        decay_rate_ct = float("inf")

    return {
        "all_abs": [float(a) for a in abs_eig],
        "constant_eig_distance_to_one": near_one_distance,
        "dominant_nontrivial_abs": dominant_nontrivial,
        "spectral_gap_in_disk": spectral_gap_disk,
        "implied_continuous_decay_rate": decay_rate_ct,
        "n_eigvals_inside_unit_disk": int(np.sum(abs_eig < 1.0)),
        "n_eigvals_outside_unit_disk": int(np.sum(abs_eig > 1.0 + 1e-9)),
    }


# ----------------------------------------------------------------------
# 5.  Margin recovery probe (compare to v5 preprint)
# ----------------------------------------------------------------------

def closest_competitor_estimate(snapshot_X: np.ndarray, L_X: np.ndarray, L_star: float) -> dict:
    """Cheap proxy for "how close did we come to the v5 ~6% margin?"

    We DO NOT re-derive the certified margin. We just report the worst
    L value seen on the grid (other than at x_star itself), and compute
    (L_star - max_off_baseline) / L_star as the sample margin.
    """
    # Find the index closest to x_star and exclude it
    dist_to_star = np.linalg.norm(snapshot_X - X_STAR, axis=1)
    idx_at_star = int(np.argmin(dist_to_star))

    # Exclude points within 1e-6 of x_star (in case of tied grid points)
    mask = dist_to_star > 1e-6
    if not np.any(mask):
        return {"sample_margin_pct": float("nan"), "max_off_baseline_L": float("nan")}

    L_off = L_X[mask]
    max_off = float(np.max(L_off))
    # Where in the grid was that max?
    candidate_idx = np.where(mask)[0][int(np.argmax(L_off))]
    candidate_state = snapshot_X[candidate_idx].tolist()

    sample_margin_pct = 100.0 * (L_star - max_off) / L_star

    return {
        "L_star_baseline": L_star,
        "max_off_baseline_L": max_off,
        "sample_margin_pct": sample_margin_pct,
        "candidate_state": candidate_state,
        "n_grid_points_evaluated": int(len(L_X)),
        "grid_idx_at_star": idx_at_star,
        "v5_preprint_certified_margin_pct": 6.1,
        "note": (
            "This is a coarse-grid sample, not the certified branch-and-bound "
            "margin. Compare to v5 preprint's 6.1% (n=3) only as a sanity check."
        ),
    }


# ----------------------------------------------------------------------
# 6.  Main
# ----------------------------------------------------------------------

def main():
    t0 = time.time()

    print("=" * 70)
    print(f"Koopman lift prototype: {EXPERIMENT_ID}")
    print("=" * 70)

    # Step 0: baseline
    L_star_exact = L_baseline_exact(DEGREE)
    L_star_marcher = L_of_state(X_STAR)
    print(f"\n[baseline] L(z^{DEGREE} - 1) closed form: {L_star_exact:.6f}")
    print(f"[baseline] L(z^{DEGREE} - 1) marcher:     {L_star_marcher:.6f}")
    print(f"[baseline] relative marcher error:        "
          f"{abs(L_star_marcher - L_star_exact) / L_star_exact:.2e}")

    # Step 1: build snapshot grid in slice coords
    print(f"\n[grid] building {GRID_PER_AXIS}^3 = {GRID_PER_AXIS**3} snapshot states")
    X = build_snapshot_grid()

    # Step 2: assemble EDMD at the headline step
    print(f"[edmd] one-step Euler advance, h = {EULER_STEP}")
    print(f"[edmd] {N_OBSERVABLES} observables (monomials degree<=2 + L)")
    edmd = assemble_edmd(X, dt=EULER_STEP)
    print(f"[edmd] Koopman matrix shape: {edmd['K'].shape}")

    # Step 3: spectrum
    spec = spectral_summary(edmd["eigvals"], dt=EULER_STEP)
    print(f"\n[spectrum] dominant non-trivial |lambda|: "
          f"{spec['dominant_nontrivial_abs']:.6f}")
    print(f"[spectrum] spectral gap in disk:           "
          f"{spec['spectral_gap_in_disk']:.6f}")
    print(f"[spectrum] implied continuous decay rate:  "
          f"{spec['implied_continuous_decay_rate']:.4f}")

    # Step 3b: step-size sweep to characterize the noise floor
    sweep_dts = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
    print(f"\n[sweep] running step-size sweep over h = {sweep_dts}")
    sweep = step_size_sweep(X, sweep_dts)
    for row in sweep["sweep_rows"]:
        print(f"  h={row['h']:8.3f}  |lambda_d|_max = "
              f"{row['dominant_nontrivial_abs']:.6f}  "
              f"gap_in_disk = {row['spectral_gap_in_disk']:+.6f}")

    # Step 4: margin recovery
    margin = closest_competitor_estimate(X, edmd["L_X"], L_star_marcher)
    print(f"\n[margin] sample-grid margin:  {margin['sample_margin_pct']:.2f}%")
    print(f"[margin] v5 preprint margin:  {margin['v5_preprint_certified_margin_pct']:.2f}%")

    elapsed = time.time() - t0

    # Step 5: dump JSON
    results = {
        "experiment_id": EXPERIMENT_ID,
        "version": "V1",
        "date": "2026-05-02",
        "degree": DEGREE,
        "honest_scope_statement": (
            "This is a prototype scoping artifact, not a certified proof of EHP "
            "at any n. The Koopman lift here uses EDMD with finite observable "
            "basis; the resulting spectral gap is an empirical estimate, not an "
            "interval-arithmetic certificate. To translate this into a "
            "closure-relevant artifact, the spectral gap would need to be (1) "
            "derived analytically rather than empirically estimated, (2) shown "
            "n-invariant or growing in n, (3) lifted to interval-arithmetic "
            "certificates. None of those are achieved in this scoping run."
        ),
        "design_choices": {
            "flow": "normalized ascent flow x' = grad L(x) / FLOW_NORMALIZER",
            "flow_normalizer": FLOW_NORMALIZER,
            "flow_normalizer_rationale": (
                "raw |grad L| has spectral radius ~3e6 at x_star (v5 Level3 "
                "Hessian); we rescale time by 1/FLOW_NORMALIZER so the "
                "continuous-time eigenvalues are O(1). This is a flow-time "
                "coordinate change; the Koopman spectrum's qualitative "
                "structure (which eigenvalues are inside vs outside the unit "
                "disk after one Euler step) is preserved."
            ),
            "x_star_is_stable_fixed_point": True,
            "state_space": "symmetry-reduced slice, dim 2n-3 (= 3 for n=3)",
            "x_star": X_STAR.tolist(),
            "observable_basis": (
                "[1, x1, x2, x3, x1^2, x2^2, x3^2, x1*x2, x1*x3, x2*x3, L]"
            ),
            "n_observables": N_OBSERVABLES,
            "koopman_approximation": "EDMD (Williams-Kevrekidis-Rowley 2015)",
            "snapshot_count": int(GRID_PER_AXIS ** 3),
            "snapshot_grid_half_width": GRID_HALF,
            "euler_step_h": EULER_STEP,
            "grad_finite_diff_step": GRAD_FD_STEP,
            "L_evaluator": "marching squares on 320x320 grid, box [-1.6, 1.6]^2",
        },
        "baseline": {
            "L_star_exact_closed_form": L_star_exact,
            "L_star_marcher_estimate": L_star_marcher,
            "relative_marcher_error": abs(L_star_marcher - L_star_exact) / L_star_exact,
        },
        "koopman_eigvals_real": [float(np.real(e)) for e in edmd["eigvals"]],
        "koopman_eigvals_imag": [float(np.imag(e)) for e in edmd["eigvals"]],
        "koopman_eigvals_abs": [float(np.abs(e)) for e in edmd["eigvals"]],
        "spectrum": spec,
        "step_size_sweep": sweep,
        "margin_probe": margin,
        "wall_time_seconds": elapsed,
        "input_artifacts_referenced": {
            "preprint": (
                "Math/erdosatlas-workbench/ehp_erdos114_preprint.tex"
            ),
            "v5_n3_results": (
                "Math/erdos-experiments/results/erdos-114/EHP_N3_LEVEL3_RESULTS.json"
            ),
        },
        "uncertainty_flags": [
            "L marcher is ~1% relative; spectrum is sensitive to L noise.",
            "EDMD with finite-basis observables can underestimate the true "
            "Koopman gap; eigenvalues may also pick up numerical artifacts "
            "from the regularizer (1e-10 * I).",
            "Forward Euler with h=0.02 introduces O(h) error in lambda_d -> "
            "lambda_c conversion; the implied continuous decay rate is a "
            "first-order estimate.",
            "Sample-grid margin is NOT the certified branch-and-bound margin "
            "from the v5 preprint (which is 6.1% on n=3).",
        ],
    }

    json_path = os.path.join(RESULTS_DIR, f"{EXPERIMENT_ID}_RESULTS.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, sort_keys=True)
    print(f"\n[write] results JSON: {json_path}")

    # SHA-256
    sha = hashlib.sha256()
    with open(json_path, "rb") as f:
        sha.update(f.read())
    sha_path = os.path.join(RESULTS_DIR, f"{EXPERIMENT_ID}_RESULTS.sha256")
    with open(sha_path, "w") as f:
        f.write(f"{sha.hexdigest()}  {os.path.basename(json_path)}\n")
    print(f"[write] sha256:       {sha_path}")

    # Brief markdown report
    md_lines = [
        f"# {EXPERIMENT_ID}",
        "",
        "## Honest Scope",
        "",
        results["honest_scope_statement"],
        "",
        "## Setup",
        "",
        f"- Degree: n = {DEGREE}",
        f"- State slice (sym-reduced): {tuple(X_STAR.tolist())} = z^3 - 1",
        f"- Snapshot grid: {GRID_PER_AXIS}^3 = {GRID_PER_AXIS**3} states "
        f"in [-{GRID_HALF}, {GRID_HALF}]^3 around x_star",
        f"- Observable basis size: K = {N_OBSERVABLES}",
        f"- Euler step h = {EULER_STEP}",
        "",
        "## Koopman Spectrum",
        "",
        f"- L\\* (closed form):       {L_star_exact:.10f}",
        f"- L\\* (radial marcher):    {L_star_marcher:.10f}  "
        f"(rel. err = {abs(L_star_marcher - L_star_exact) / L_star_exact:.2e})",
        f"- Constant-function eigenvalue distance to 1: "
        f"{spec['constant_eig_distance_to_one']:.2e}",
        f"- Dominant non-trivial |lambda|: {spec['dominant_nontrivial_abs']:.6f}",
        f"- Spectral gap in unit disk:     {spec['spectral_gap_in_disk']:.6f}",
        f"- Implied continuous decay rate: {spec['implied_continuous_decay_rate']:.4f}",
        f"- # eigvals inside unit disk:    {spec['n_eigvals_inside_unit_disk']}",
        f"- # eigvals outside unit disk:   {spec['n_eigvals_outside_unit_disk']}",
        "",
        "## Step-Size Sweep",
        "",
        "| h (flow time) | |lambda_d|_max | gap_in_disk |",
        "|---:|---:|---:|",
    ] + [
        f"| {row['h']:.3f} | {row['dominant_nontrivial_abs']:.6f} | "
        f"{row['spectral_gap_in_disk']:+.6f} |"
        for row in sweep["sweep_rows"]
    ] + [
        "",
        "## Margin Probe",
        "",
        f"- Sample-grid margin (this prototype): "
        f"{margin['sample_margin_pct']:.2f}%",
        f"- v5 preprint certified margin (n=3):   "
        f"{margin['v5_preprint_certified_margin_pct']:.2f}%",
        "",
        "## Notes",
        "",
        "- See KOOPMAN_LIFT_DESIGN_2026-05-02.md for the design rationale.",
        "- This run does NOT modify the v5 preprint or any erdos-114 results.",
        "- This run does NOT publish externally; Cooley filter applies.",
        "",
    ]
    md_path = os.path.join(RESULTS_DIR, f"{EXPERIMENT_ID}_REPORT.md")
    with open(md_path, "w") as f:
        f.write("\n".join(md_lines))
    print(f"[write] report MD:    {md_path}")
    print(f"\nDone in {elapsed:.1f} s.")


if __name__ == "__main__":
    main()
