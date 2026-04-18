#!/usr/bin/env python3
"""
EXP-009C: Spectral Torus Phase Portrait

Compute FULL eigenvalue spectra for small W (dense matrix, exact), then:
  1. Map eigenvalues to torus coordinates (|λ|, arg(λ))
  2. Track how spectral density evolves with W → spectral flow
  3. Apply Koopman/DMD to the spectral DISTRIBUTION (not just λ_max)
  4. Detect spectral phase transition via level spacing statistics (GUE/GOE/Poisson)
  5. Extract the renormalization group beta function from the spectral flow

The key insight: the transfer matrix T(W) has REAL entries (it's a 0/1 matrix),
so eigenvalues are either real or come in complex conjugate pairs. The
spectral torus is actually S¹ (the unit circle in the complex plane,
after normalizing by λ_max). As W → ∞, the spectral distribution on S¹
encodes the thermodynamic limit.

This is the proper implementation of the torus/KvN idea.
"""

import json
import math
import numpy as np
from collections import defaultdict


# ============================================================
# DENSE TRANSFER MATRIX (for full eigenvalue spectrum)
# ============================================================

def enumerate_sidon_states(W):
    """Enumerate all Sidon subsets of {0,...,W-1}."""
    states = []

    def backtrack(start, current, diffs):
        mask = frozenset(current)
        states.append(tuple(sorted(current)))

        for x in range(start, W):
            new_diffs = []
            valid = True
            for a in current:
                d = abs(x - a)
                if d in diffs or d in new_diffs:
                    valid = False
                    break
                new_diffs.append(d)
            if valid:
                current.append(x)
                diffs_new = diffs | set(new_diffs)
                backtrack(x + 1, current, diffs_new)
                current.pop()

    backtrack(0, [], set())
    return states


def build_dense_transfer_matrix(W):
    """Build the full transfer matrix as a dense numpy array."""
    states = enumerate_sidon_states(W)
    state_set = {s: i for i, s in enumerate(states)}
    n = len(states)

    T = np.zeros((n, n))

    for i, state in enumerate(states):
        elements = set(state)
        diffs = set()
        for a in elements:
            for b in elements:
                if a < b:
                    diffs.add(b - a)

        # Shift right: each position p → p-1, drop those at 0
        shifted = tuple(sorted(p - 1 for p in elements if p > 0))
        shifted_set = set(shifted)
        shifted_diffs = set()
        for a in shifted_set:
            for b in shifted_set:
                if a < b:
                    shifted_diffs.add(b - a)

        # Option 1: don't add W-1
        if shifted in state_set:
            j = state_set[shifted]
            T[i, j] += 1.0

        # Option 2: add W-1 to the shifted state
        new_pos = W - 1
        valid = True
        for p in shifted_set:
            d = new_pos - p
            if d in shifted_diffs:
                valid = False
                break

        # Also check that new diffs among themselves don't collide
        if valid:
            new_diffs_list = [new_pos - p for p in shifted_set]
            if len(new_diffs_list) != len(set(new_diffs_list)):
                valid = False

        if valid:
            new_state = tuple(sorted(list(shifted) + [new_pos]))
            if new_state in state_set:
                j = state_set[new_state]
                T[i, j] += 1.0

    return T, states


def full_spectrum(W):
    """Compute full eigenvalue spectrum at window size W."""
    T, states = build_dense_transfer_matrix(W)
    eigenvalues = np.linalg.eigvals(T)
    # Sort by magnitude (largest first)
    idx = np.argsort(-np.abs(eigenvalues))
    eigenvalues = eigenvalues[idx]
    return eigenvalues, len(states)


# ============================================================
# SPECTRAL TORUS ANALYSIS
# ============================================================

def level_spacing_ratio(eigenvalues):
    """
    Compute the level spacing ratio <r> = <min(s_i, s_{i+1}) / max(s_i, s_{i+1})>.

    This distinguishes:
      - Poisson (integrable): <r> ≈ 0.386
      - GOE (real symmetric, time-reversal): <r> ≈ 0.530
      - GUE (complex hermitian): <r> ≈ 0.603

    For our real 0/1 matrix, GOE is the expected universality class
    if the system is chaotic; Poisson if integrable.
    """
    # Use absolute values (real eigenvalues sorted by magnitude)
    real_eigs = np.sort(np.abs(eigenvalues))[::-1]

    # Only use eigenvalues that are significantly nonzero
    threshold = 1e-10 * real_eigs[0] if len(real_eigs) > 0 else 0
    real_eigs = real_eigs[real_eigs > threshold]

    if len(real_eigs) < 4:
        return None, None

    spacings = np.diff(real_eigs)
    spacings = np.abs(spacings)

    if len(spacings) < 2:
        return None, None

    ratios = []
    for i in range(len(spacings) - 1):
        s1, s2 = spacings[i], spacings[i + 1]
        if max(s1, s2) > 1e-15:
            ratios.append(min(s1, s2) / max(s1, s2))

    if not ratios:
        return None, None

    r_mean = np.mean(ratios)

    # Classify
    if r_mean < 0.45:
        universality = "Poisson (integrable)"
    elif r_mean < 0.57:
        universality = "GOE (time-reversal symmetric chaos)"
    else:
        universality = "GUE (broken time-reversal)"

    return r_mean, universality


def spectral_entropy(eigenvalues):
    """Shannon entropy of the eigenvalue magnitude distribution."""
    mags = np.abs(eigenvalues)
    mags = mags[mags > 1e-15]
    if len(mags) == 0:
        return 0.0
    p = mags / np.sum(mags)
    return -np.sum(p * np.log(p + 1e-30))


def torus_coordinates(eigenvalues, lambda_max):
    """
    Map eigenvalues to torus coordinates.

    For each eigenvalue λ_i:
      - r_i = |λ_i| / λ_max  (radial coordinate, 0 to 1)
      - θ_i = arg(λ_i)        (angular coordinate, -π to π)

    Returns (r, theta) arrays.
    """
    r = np.abs(eigenvalues) / abs(lambda_max) if abs(lambda_max) > 0 else np.abs(eigenvalues)
    theta = np.angle(eigenvalues)
    return r, theta


def spectral_moments(eigenvalues, k_max=6):
    """Compute spectral moments: M_k = (1/N) Σ |λ_i|^k / λ_max^k."""
    lmax = np.max(np.abs(eigenvalues))
    if lmax < 1e-15:
        return {}

    normalized = np.abs(eigenvalues) / lmax
    moments = {}
    for k in range(1, k_max + 1):
        moments[f"M_{k}"] = np.mean(normalized ** k)
    return moments


# ============================================================
# KOOPMAN ON SPECTRAL FLOW
# ============================================================

def spectral_flow_koopman(spectra_by_W, n_bins=20):
    """
    The real KvN method: treat the spectral DISTRIBUTION at each W as a
    state vector in a function space, then find the Koopman operator.

    1. Bin the normalized eigenvalue magnitudes into a histogram (n_bins)
    2. Stack histograms: each row = one W, each column = one bin
    3. DMD on the histogram sequence → Koopman eigenvalues of the spectral flow
    4. Extrapolate the spectral distribution to W=∞
    """
    W_list = sorted(spectra_by_W.keys())
    if len(W_list) < 4:
        return None

    # Build histogram matrix
    hist_matrix = []
    for W in W_list:
        eigs = spectra_by_W[W]
        lmax = np.max(np.abs(eigs))
        if lmax < 1e-15:
            continue
        normalized = np.abs(eigs) / lmax
        hist, bin_edges = np.histogram(normalized, bins=n_bins, range=(0, 1), density=True)
        hist_matrix.append(hist)

    if len(hist_matrix) < 4:
        return None

    H = np.array(hist_matrix)  # shape: (n_W, n_bins)

    # DMD: find A such that H[i+1] ≈ A @ H[i]
    X = H[:-1].T  # (n_bins, n_W-1)
    Y = H[1:].T   # (n_bins, n_W-1)

    U, S, Vh = np.linalg.svd(X, full_matrices=False)
    rank = max(1, np.sum(S > 0.01 * S[0]))
    rank = min(rank, len(S))

    Ur = U[:, :rank]
    Sr = S[:rank]
    Vr = Vh[:rank, :]

    Atilde = Ur.T @ Y @ Vr.T @ np.diag(1.0 / Sr)
    koopman_eigs, koopman_vecs = np.linalg.eig(Atilde)

    # Sort by magnitude
    idx = np.argsort(-np.abs(koopman_eigs))
    koopman_eigs = koopman_eigs[idx]

    # DMD modes in full space
    modes = Y @ Vr.T @ np.diag(1.0 / Sr) @ koopman_vecs[:, idx]

    # Extrapolate: predict histogram at W=∞ (many steps forward)
    x0 = H[-1]  # last observed histogram
    b = np.linalg.lstsq(modes, x0, rcond=None)[0]

    # Predict 50 steps ahead
    x_future = modes @ (b * koopman_eigs ** 50)
    hist_inf = x_future.real

    # The spectral density at W=∞ tells us the thermodynamic limit
    # The fraction of eigenvalue weight at r=1 is the "condensate fraction"
    condensate_frac = hist_inf[-1] / (np.sum(np.abs(hist_inf)) + 1e-30) if len(hist_inf) > 0 else 0

    return {
        "koopman_eigenvalues": [{"real": e.real, "imag": e.imag, "abs": abs(e)} for e in koopman_eigs],
        "koopman_rank": int(rank),
        "n_W_used": len(hist_matrix),
        "extrapolated_histogram": hist_inf.tolist(),
        "bin_edges": bin_edges.tolist(),
        "condensate_fraction": float(condensate_frac),
        "dominant_koopman_abs": float(np.abs(koopman_eigs[0])),
    }


# ============================================================
# RENORMALIZATION GROUP BETA FUNCTION
# ============================================================

def rg_beta_function(spectra_by_W):
    """
    The RG beta function: how does each spectral observable change with scale?

    β(g) = dg/d(ln W) where g is a coupling constant (spectral observable).

    If β(g*) = 0, we have a fixed point → the system is critical at W=∞.
    The eigenvalues of dβ/dg at the fixed point give the critical exponents.
    """
    W_list = sorted(spectra_by_W.keys())
    if len(W_list) < 4:
        return None

    # Track observables as functions of W
    observables = defaultdict(list)

    for W in W_list:
        eigs = spectra_by_W[W]
        lmax = np.max(np.abs(eigs))

        # Gap ratio
        sorted_abs = np.sort(np.abs(eigs))[::-1]
        if len(sorted_abs) >= 2 and sorted_abs[0] > 1e-15:
            observables["gap_ratio"].append((W, sorted_abs[1] / sorted_abs[0]))

        # Spectral entropy (normalized)
        se = spectral_entropy(eigs)
        se_max = np.log(len(eigs)) if len(eigs) > 1 else 1
        observables["spectral_entropy_norm"].append((W, se / se_max))

        # Level spacing ratio
        r_mean, _ = level_spacing_ratio(eigs)
        if r_mean is not None:
            observables["level_spacing_ratio"].append((W, r_mean))

        # Moments
        moments = spectral_moments(eigs)
        for mk, mv in moments.items():
            observables[mk].append((W, mv))

        # Fraction of real eigenvalues
        n_real = np.sum(np.abs(eigs.imag) < 1e-10 * lmax)
        observables["real_fraction"].append((W, n_real / len(eigs)))

    # Compute beta functions
    betas = {}
    for name, data in observables.items():
        if len(data) < 4:
            continue

        Ws = np.array([d[0] for d in data], dtype=float)
        gs = np.array([d[1] for d in data])

        # β(g) = dg/d(ln W)
        log_W = np.log(Ws)
        # Numerical derivative
        dg = np.diff(gs)
        dlogW = np.diff(log_W)
        beta_values = dg / dlogW
        W_mid = np.exp((log_W[:-1] + log_W[1:]) / 2)

        # Is beta → 0? (fixed point)
        if len(beta_values) >= 3:
            last_3 = beta_values[-3:]
            trend = np.polyfit(np.arange(3), last_3, 1)
            approaching_zero = abs(last_3[-1]) < abs(beta_values[0]) * 0.5

            betas[name] = {
                "last_value": float(gs[-1]),
                "last_beta": float(beta_values[-1]),
                "beta_trend_slope": float(trend[0]),
                "approaching_fixed_point": approaching_zero,
                "values": [(float(w), float(g)) for w, g in data[-5:]],
            }

    return betas


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 80)
    print("EXP-009C: Spectral Torus Phase Portrait")
    print("  Full eigenvalue spectra → torus mapping → Koopman on spectral flow")
    print("=" * 80)

    max_W = 16  # Dense eigendecomposition limit in sandbox

    spectra = {}
    all_results = []

    print(f"\n{'W':>4} | {'#states':>8} | {'λ_max':>10} | {'λ₂/λ₁':>8} | {'<r>':>6} | {'class':>10} | {'S_norm':>7} | {'f_real':>6}")
    print("-" * 80)

    for W in range(3, max_W + 1):
        eigenvalues, n_states = full_spectrum(W)
        spectra[W] = eigenvalues

        lambda_max = np.abs(eigenvalues[0])
        lambda_2 = np.abs(eigenvalues[1]) if len(eigenvalues) > 1 else 0
        gap_ratio = lambda_2 / lambda_max if lambda_max > 0 else 0

        r_mean, universality = level_spacing_ratio(eigenvalues)
        se = spectral_entropy(eigenvalues)
        se_norm = se / np.log(max(2, len(eigenvalues)))

        n_real = np.sum(np.abs(eigenvalues.imag) < 1e-10 * lambda_max)
        f_real = n_real / len(eigenvalues)

        moments = spectral_moments(eigenvalues)

        r_torus, theta_torus = torus_coordinates(eigenvalues, lambda_max)

        r_str = f"{r_mean:.4f}" if r_mean is not None else "  N/A "
        class_str = universality[:10] if universality else "N/A"

        print(f"{W:4} | {n_states:8} | {lambda_max:10.6f} | {gap_ratio:8.4f} | {r_str} | {class_str:>10} | {se_norm:7.4f} | {f_real:6.3f}")

        all_results.append({
            "W": W,
            "n_states": n_states,
            "n_eigenvalues": len(eigenvalues),
            "lambda_max": float(lambda_max),
            "gap_ratio": float(gap_ratio),
            "level_spacing_ratio": float(r_mean) if r_mean is not None else None,
            "universality_class": universality,
            "spectral_entropy_normalized": float(se_norm),
            "real_fraction": float(f_real),
            "moments": {k: float(v) for k, v in moments.items()},
            "torus_radii_mean": float(np.mean(r_torus)),
            "torus_radii_std": float(np.std(r_torus)),
            "n_distinct_phases": int(np.sum(np.abs(np.diff(np.sort(theta_torus))) > 0.01)),
        })

    # ---- Spectral Flow Koopman ----
    print(f"\n{'='*60}")
    print("KOOPMAN ON SPECTRAL FLOW")
    print(f"{'='*60}")

    for n_bins in [10, 15, 20]:
        koop = spectral_flow_koopman(spectra, n_bins=n_bins)
        if koop:
            print(f"\n  Bins={n_bins}:")
            print(f"    Koopman rank: {koop['koopman_rank']}")
            dom = koop['koopman_eigenvalues'][0]
            print(f"    Dominant Koopman eigenvalue: |μ₁| = {dom['abs']:.6f}")
            if dom['abs'] < 1:
                print(f"    → Spectral distribution CONVERGES (contractive)")
                print(f"    → Rate: {-np.log(dom['abs']):.4f} per W step")
            else:
                print(f"    → Spectral distribution DIVERGES or oscillates")

            print(f"    Condensate fraction at W→∞: {koop['condensate_fraction']:.4f}")

            # Show top Koopman eigenvalues
            for i, ke in enumerate(koop['koopman_eigenvalues'][:5]):
                print(f"    μ_{i+1} = {ke['abs']:.4f}∠{np.degrees(np.arctan2(ke['imag'], ke['real'])):.1f}°")

    # ---- RG Beta Function ----
    print(f"\n{'='*60}")
    print("RENORMALIZATION GROUP BETA FUNCTIONS")
    print(f"{'='*60}")

    betas = rg_beta_function(spectra)
    if betas:
        for name, data in sorted(betas.items()):
            arrow = "→ 0 (FIXED POINT)" if data["approaching_fixed_point"] else "≠ 0 (flowing)"
            print(f"\n  {name}:")
            print(f"    Current value: g = {data['last_value']:.6f}")
            print(f"    β(g) = {data['last_beta']:.6f}  {arrow}")
            print(f"    β trend: {'↘ (decelerating)' if data['beta_trend_slope'] < 0 else '↗ (accelerating)'}")

    # ---- Phase Transition Detection ----
    print(f"\n{'='*60}")
    print("PHASE TRANSITION DETECTION")
    print(f"{'='*60}")

    # Track gap_ratio(W) — if it → 1, we have a phase transition
    gap_ratios = [(r["W"], r["gap_ratio"]) for r in all_results]
    if len(gap_ratios) >= 3:
        Ws = np.array([g[0] for g in gap_ratios], dtype=float)
        grs = np.array([g[1] for g in gap_ratios])

        # Fit: 1 - gap_ratio ~ W^{-ν}
        deficit = 1.0 - grs
        if np.all(deficit > 0):
            slope, intercept = np.polyfit(np.log(Ws), np.log(deficit), 1)
            nu = -slope
            print(f"\n  Gap closing: 1 - λ₂/λ₁ ~ W^(-{nu:.4f})")
            print(f"  → Spectral gap closes as W^(-{nu:.4f})")
            print(f"  → At W=50: gap_ratio ≈ {1.0 - np.exp(intercept) * 50**slope:.4f}")
            print(f"  → At W=100: gap_ratio ≈ {1.0 - np.exp(intercept) * 100**slope:.4f}")
            print(f"  → At W=∞: gap_ratio → 1.0 (continuous phase transition)")

        # Level spacing ratio trend
        lsr = [(r["W"], r["level_spacing_ratio"]) for r in all_results
               if r["level_spacing_ratio"] is not None]
        if len(lsr) >= 3:
            lsr_vals = [l[1] for l in lsr]
            print(f"\n  Level spacing ratio <r>:")
            print(f"    Trend: {lsr_vals[0]:.4f} → {lsr_vals[-1]:.4f}")
            if lsr_vals[-1] > 0.50:
                print(f"    → Approaching GOE (random matrix, chaotic)")
            elif lsr_vals[-1] < 0.42:
                print(f"    → Poisson-like (integrable)")
            else:
                print(f"    → Intermediate (crossover region)")

    # ---- Summary ----
    print(f"\n{'='*60}")
    print("TORUS PORTRAIT SUMMARY")
    print(f"{'='*60}")

    print(f"\n  Data: W = 3..{max_W} ({len(spectra)} spectra)")
    print(f"  States at W={max_W}: {all_results[-1]['n_states']}")
    print(f"  Eigenvalues at W={max_W}: {all_results[-1]['n_eigenvalues']}")

    # Key observables at largest W
    last = all_results[-1]
    print(f"\n  At W={max_W}:")
    print(f"    λ_max = {last['lambda_max']:.6f}")
    print(f"    λ₂/λ₁ = {last['gap_ratio']:.4f}")
    print(f"    <r> = {last['level_spacing_ratio']:.4f}" if last['level_spacing_ratio'] else "    <r> = N/A")
    print(f"    Class: {last['universality_class']}")
    print(f"    Entropy (norm): {last['spectral_entropy_normalized']:.4f}")
    print(f"    Real fraction: {last['real_fraction']:.3f}")

    print(f"\n  For the Rust crate (W=50+), run:")
    print(f"    cargo run --release -- -W 60 -k 20 --spectral-torus -o spectral_torus_60.json")

    # Save
    output = {
        "experiment_id": "EXP-009C",
        "title": "Spectral Torus Phase Portrait",
        "max_W": max_W,
        "per_W_data": all_results,
    }

    if betas:
        output["rg_beta_functions"] = {
            k: {kk: vv for kk, vv in v.items() if kk != "values"}
            for k, v in betas.items()
        }

    with open("EXP-009C_SPECTRAL_TORUS_RESULTS.json", "w") as f:
        json.dump(output, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else str(x))

    print(f"\nSaved to EXP-009C_SPECTRAL_TORUS_RESULTS.json")


if __name__ == "__main__":
    main()
