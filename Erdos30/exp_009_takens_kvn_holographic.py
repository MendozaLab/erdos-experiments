#!/usr/bin/env python3
"""
EXP-009: Takens Embedding + KvN Holographic Extrapolation

Idea (Mendoza): The eigenvalue spectrum {λ_i(W)} for each window size W
lives naturally on a torus (phases of complex eigenvalues). Instead of
brute-force power-law fitting of λ_max(W), embed the spectrum using
Takens' theorem to reconstruct the underlying dynamical attractor, then
use the attractor geometry to extrapolate λ_∞ = lim_{W→∞} λ_max(W).

The KvN (Koopman-von Neumann) formalism lifts the classical map
W → λ_max(W) to a linear operator on L²(torus). The eigenvalues of
the Koopman operator encode the asymptotic dynamics — the boundary
data (finite W) holographically encodes the bulk (W→∞).

Three-pronged attack:
  1. Takens embedding: embed the sequence {λ_max(W)} in R^d using
     delay coordinates. The attractor dimension reveals how many
     "degrees of freedom" the convergence has.
  2. Koopman/DMD: fit a linear operator to the delay-embedded data.
     The Koopman eigenvalues give the decay rates → extrapolate.
  3. Torus geometry: map eigenvalue phases to a torus, compute
     winding numbers and rotation numbers → detect quasiperiodic
     structure that constrains the limit.

If this works, we get λ_∞ from W≤47 data that would otherwise
require W=200+ by brute force. Geometry accelerates compute.
"""

import json
import math
import time

import numpy as np


# ============================================================
# DATA: C++ results from EXP-007C (W=3..47)
# ============================================================

# λ_max values from the C++ transfer matrix computation
SIDON_DATA = {
    3: 1.618034,
    4: 1.839287,
    5: 1.754878,
    6: 1.725498,
    7: 1.669696,
    8: 1.659010,
    9: 1.620138,
    10: 1.607803,
    11: 1.578709,
    12: 1.566497,
    13: 1.545082,
    14: 1.533478,
    15: 1.516586,
    16: 1.506160,
    17: 1.492412,
    18: 1.483108,
    19: 1.471639,
    20: 1.452216,
    21: 1.443527,
    22: 1.434271,
    23: 1.422800,
    24: 1.414916,
    25: 1.405792,
}

# Sum-free data from EXP-008 (known answer: density 1/2)
SUMFREE_DATA = {}  # Will be loaded from EXP-008 results


def load_sidon_full():
    """Load full Sidon data from EXP-007C JSON if available."""
    try:
        with open("EXP-MATH-ERDOS30-SIDON-007C_RESULTS.json") as f:
            data = json.load(f)
        result = {}
        for row in data.get("data", []):
            result[row["W"]] = row["lambda_max"]
        if len(result) > len(SIDON_DATA):
            return result
    except (FileNotFoundError, json.JSONDecodeError):
        pass
    return SIDON_DATA


def load_sumfree():
    """Load sum-free data from EXP-008 JSON."""
    try:
        with open("EXP-008_SUMFREE_CALIBRATION_RESULTS.json") as f:
            data = json.load(f)
        result = {}
        for row in data.get("data", []):
            result[row["W"]] = row["lambda_max"]
        return result
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


# ============================================================
# 1. TAKENS EMBEDDING
# ============================================================

def takens_embed(series, dim, tau=1):
    """
    Takens delay embedding of a 1D time series.

    Given x[0], x[1], ..., x[N-1], construct vectors:
      v[i] = (x[i], x[i+tau], x[i+2*tau], ..., x[i+(dim-1)*tau])

    Returns (N - (dim-1)*tau) x dim array.
    """
    N = len(series)
    M = N - (dim - 1) * tau
    if M <= 0:
        raise ValueError(f"Series too short ({N}) for dim={dim}, tau={tau}")

    embedded = np.zeros((M, dim))
    for d in range(dim):
        embedded[:, d] = series[d * tau: d * tau + M]

    return embedded


def estimate_embedding_dimension(series, max_dim=8, tau=1):
    """
    Estimate optimal embedding dimension using false nearest neighbors (FNN).

    Simple version: compute the fraction of "false neighbors" that disappear
    when embedding dimension increases by 1.
    """
    results = []

    for dim in range(1, max_dim + 1):
        try:
            emb = takens_embed(series, dim, tau)
        except ValueError:
            break

        if len(emb) < 3:
            break

        # Compute pairwise distances
        n = len(emb)
        # Find nearest neighbor for each point
        false_nn = 0
        total = 0

        if dim < max_dim:
            try:
                emb_next = takens_embed(series, dim + 1, tau)
            except ValueError:
                break

            for i in range(min(n, len(emb_next))):
                # Find nearest neighbor in dim-dimensional embedding
                dists = np.sqrt(np.sum((emb[:n] - emb[i]) ** 2, axis=1))
                dists[i] = np.inf  # exclude self
                nn_idx = np.argmin(dists)
                nn_dist = dists[nn_idx]

                if nn_dist < 1e-15:
                    continue

                # Check if this neighbor is "false" — distance changes significantly
                # in (dim+1)-dimensional embedding
                if nn_idx < len(emb_next) and i < len(emb_next):
                    dist_higher = np.sqrt(np.sum((emb_next[nn_idx] - emb_next[i]) ** 2))
                    ratio = abs(dist_higher - nn_dist) / nn_dist
                    if ratio > 2.0:  # threshold for false neighbor
                        false_nn += 1
                total += 1

        fnn_frac = false_nn / total if total > 0 else 0
        results.append({
            "dim": dim,
            "fnn_fraction": fnn_frac,
            "n_points": len(emb),
        })

    return results


# ============================================================
# 2. DYNAMIC MODE DECOMPOSITION (Koopman approximation)
# ============================================================

def dmd(X, Y, rank=None):
    """
    Dynamic Mode Decomposition: find best-fit linear operator A such that
    Y ≈ A @ X, where X and Y are delay-embedded snapshots.

    Returns eigenvalues (Koopman eigenvalues), modes, and amplitudes.

    This is the Koopman operator approximation — the linear lift of
    the nonlinear dynamics W → λ_max(W) into a Hilbert space.
    """
    U, S, Vh = np.linalg.svd(X.T, full_matrices=False)

    if rank is None:
        # Auto-select rank: keep singular values > 1% of max
        rank = max(1, np.sum(S > 0.01 * S[0]))

    rank = min(rank, len(S))

    Ur = U[:, :rank]
    Sr = S[:rank]
    Vr = Vh[:rank, :]

    # Reduced DMD operator
    Atilde = Ur.T @ Y.T @ Vr.T @ np.diag(1.0 / Sr)

    # Eigendecomposition of reduced operator
    eigvals, eigvecs_r = np.linalg.eig(Atilde)

    # DMD modes (project back to full space)
    modes = Y.T @ Vr.T @ np.diag(1.0 / Sr) @ eigvecs_r

    return eigvals, modes, Atilde


def koopman_extrapolate(series, embed_dim, n_future, tau=1):
    """
    Use Koopman/DMD to extrapolate a time series.

    1. Takens embed the series
    2. Fit Koopman operator via DMD
    3. Propagate forward n_future steps
    4. Extract the first coordinate (original observable)
    """
    emb = takens_embed(series, embed_dim, tau)

    # X = emb[:-1], Y = emb[1:]  (consecutive pairs)
    X = emb[:-1]
    Y = emb[1:]

    eigvals, modes, Atilde = dmd(X, Y)

    # Sort by magnitude (dominant modes first)
    idx = np.argsort(-np.abs(eigvals))
    eigvals = eigvals[idx]
    modes = modes[:, idx]

    # Initial condition: last embedded point
    x0 = emb[-1]

    # Solve for amplitudes: x0 = Φ @ b
    b = np.linalg.lstsq(modes, x0, rcond=None)[0]

    # Propagate
    predictions = []
    for k in range(1, n_future + 1):
        xk = modes @ (b * eigvals**k)
        predictions.append(xk[0].real)  # first coordinate = original observable

    return {
        "predictions": predictions,
        "koopman_eigenvalues": eigvals,
        "koopman_modes": modes,
        "amplitudes": b,
        "extrapolated_limit": predictions[-1] if predictions else None,
    }


# ============================================================
# 3. TORUS GEOMETRY / WINDING NUMBERS
# ============================================================

def torus_analysis(series):
    """
    Map the eigenvalue convergence to torus coordinates.

    The sequence λ_max(W) → λ_∞ can be decomposed:
      λ_max(W) = λ_∞ + Σ_k A_k × exp(i ω_k W) × R_k^W

    where (ω_k, R_k) are the "torus coordinates" — angular frequency
    and radial decay. Pure exponential decay has ω=0. Oscillatory
    convergence has ω≠0 (the sequence spirals in on the torus).

    The winding numbers ω_k / (2π) tell us about quasiperiodic structure.
    """
    # Compute successive differences
    diffs = np.diff(series)

    # Compute successive ratios of differences (convergence acceleration)
    ratios = []
    for i in range(len(diffs) - 1):
        if abs(diffs[i]) > 1e-15:
            ratios.append(diffs[i + 1] / diffs[i])
    ratios = np.array(ratios)

    # If ratios are roughly constant → geometric convergence
    # If ratios oscillate → quasiperiodic component
    ratio_mean = np.mean(ratios) if len(ratios) > 0 else 0
    ratio_std = np.std(ratios) if len(ratios) > 0 else 0

    # FFT of the differences to find oscillatory components
    if len(diffs) >= 4:
        fft_diffs = np.fft.fft(diffs)
        freqs = np.fft.fftfreq(len(diffs))
        power = np.abs(fft_diffs) ** 2

        # Find dominant frequencies (exclude DC)
        power_no_dc = power.copy()
        power_no_dc[0] = 0
        dominant_idx = np.argsort(-power_no_dc)[:3]
        dominant_freqs = freqs[dominant_idx]
        dominant_power = power[dominant_idx]
    else:
        dominant_freqs = np.array([])
        dominant_power = np.array([])

    # Winding numbers
    winding = dominant_freqs * 2 * np.pi if len(dominant_freqs) > 0 else np.array([])

    # Richardson-like extrapolation using the geometric structure
    # If λ(W) ≈ λ_∞ + c₁ r₁^W + c₂ r₂^W, and we know r₁, r₂,
    # we can solve for λ_∞ from 3 data points
    if len(series) >= 6 and abs(ratio_mean) < 1:
        # Use last 6 points for a 2-term Richardson
        tail = series[-6:]
        # Fit: λ(W) = a + b*r^W
        # Three equations, three unknowns
        # Use Aitken's Δ² method for acceleration
        aitken_estimates = []
        for i in range(len(series) - 2):
            s0, s1, s2 = series[i], series[i + 1], series[i + 2]
            denom = s2 - 2 * s1 + s0
            if abs(denom) > 1e-15:
                aitken = s0 - (s1 - s0) ** 2 / denom
                aitken_estimates.append(aitken)
    else:
        aitken_estimates = []

    return {
        "diff_ratios_mean": ratio_mean,
        "diff_ratios_std": ratio_std,
        "dominant_frequencies": dominant_freqs.tolist(),
        "dominant_power": dominant_power.tolist(),
        "winding_numbers": winding.tolist(),
        "aitken_estimates": aitken_estimates,
        "aitken_limit": np.mean(aitken_estimates[-5:]) if len(aitken_estimates) >= 5 else None,
        "convergence_type": "geometric" if ratio_std < 0.1 * abs(ratio_mean) else "quasiperiodic",
    }


# ============================================================
# 4. HOLOGRAPHIC BOUNDARY EXTRAPOLATION
# ============================================================

def holographic_extrapolate(series, Ws):
    """
    Holographic principle applied to finite-size scaling:

    The boundary (finite W data) encodes the bulk (W→∞) through
    the correlation structure. We use multiple independent methods
    to triangulate λ_∞:

    1. Power-law fit: λ - 1 = c W^{-α}
    2. Aitken Δ² acceleration
    3. Koopman/DMD extrapolation (3 embedding dims)
    4. Padé approximant (rational function fit)
    5. BST algorithm (Bulirsch-Stoer-like extrapolation)
    """
    results = {}
    series = np.array(series)
    Ws = np.array(Ws, dtype=float)

    # --- Method 1: Power-law fit ---
    excess = series - 1.0
    if np.all(excess > 0):
        log_W = np.log(Ws)
        log_ex = np.log(excess)
        A = np.vstack([np.ones_like(log_W), log_W]).T
        coef, _, _, _ = np.linalg.lstsq(A, log_ex, rcond=None)
        alpha = -coef[1]
        c = math.exp(coef[0])
        results["powerlaw"] = {
            "alpha": alpha,
            "c": c,
            "lambda_inf": 1.0,  # Power law assumes λ→1
            "formula": f"λ-1 = {c:.4f} × W^(-{alpha:.4f})",
        }

    # --- Method 2: Aitken Δ² ---
    torus = torus_analysis(series)
    results["aitken"] = {
        "lambda_inf": torus["aitken_limit"],
        "convergence_type": torus["convergence_type"],
        "last_5_estimates": torus["aitken_estimates"][-5:] if torus["aitken_estimates"] else [],
    }

    # --- Method 3: Koopman/DMD (multiple embedding dims) ---
    koopman_results = {}
    for dim in [3, 4, 5, 6]:
        if len(series) > dim + 2:
            try:
                koop = koopman_extrapolate(series, dim, n_future=50)
                koopman_results[f"dim_{dim}"] = {
                    "lambda_inf": koop["predictions"][-1],
                    "koopman_eigenvalues": [
                        {"real": e.real, "imag": e.imag, "abs": abs(e)}
                        for e in koop["koopman_eigenvalues"][:dim]
                    ],
                }
            except Exception as e:
                koopman_results[f"dim_{dim}"] = {"error": str(e)}
    results["koopman"] = koopman_results

    # --- Method 4: Padé approximant ---
    # Fit λ(1/W) as a rational function P(x)/Q(x) where x = 1/W
    # Then evaluate at x=0 (W=∞)
    x = 1.0 / Ws
    y = series

    # Fit [2/2] Padé: (a₀ + a₁x + a₂x²) / (1 + b₁x + b₂x²)
    # Rearrange: y(1 + b₁x + b₂x²) = a₀ + a₁x + a₂x²
    # y = a₀ + a₁x + a₂x² - yb₁x - yb₂x²
    # y = a₀ + (a₁ - yb₁)x + (a₂ - yb₂)x²
    if len(x) >= 5:
        # Use last N points for stability
        n_use = min(len(x), 15)
        x_fit = x[-n_use:]
        y_fit = y[-n_use:]

        # Design matrix: y = a₀ + a₁x + a₂x² - b₁xy - b₂x²y
        M = np.column_stack([
            np.ones_like(x_fit),
            x_fit,
            x_fit**2,
            -x_fit * y_fit,
            -x_fit**2 * y_fit,
        ])

        try:
            coefs, _, _, _ = np.linalg.lstsq(M, y_fit, rcond=None)
            a0, a1, a2, b1, b2 = coefs
            # λ(∞) = a₀ / 1 = a₀
            results["pade"] = {
                "lambda_inf": a0,
                "coefficients": {"a0": a0, "a1": a1, "a2": a2, "b1": b1, "b2": b2},
            }
        except Exception as e:
            results["pade"] = {"error": str(e)}

    # --- Method 5: BST (Bulirsch-Stoer) extrapolation ---
    # Use the sequence of λ_max(W) values and extrapolate to W=∞ (h=1/W→0)
    # using the BST rational interpolation tableau
    h = 1.0 / Ws
    n = len(h)
    T = np.zeros((n, n))
    T[:, 0] = series

    for j in range(1, n):
        for i in range(n - j):
            if j == 1:
                denom = (h[i] / h[i + j]) * (1 - T[i + 1, j - 1] / T[i, j - 1]) + (h[i] / h[i + j]) - 1
            else:
                ratio = (h[i] / h[i + j])
                diff = T[i + 1, j - 1] - T[i, j - 2] if j >= 2 else 0
                if abs(diff) < 1e-30:
                    T[i, j] = T[i + 1, j - 1]
                else:
                    T[i, j] = T[i + 1, j - 1] + (T[i + 1, j - 1] - T[i, j - 1]) / ((ratio) * (1 - (T[i + 1, j - 1] - T[i, j - 1]) / diff) - 1)

    # The best estimate is T[0, n-1], but convergence is from the diagonal
    bst_diagonal = [T[0, j] for j in range(min(n, 10)) if not np.isnan(T[0, j]) and abs(T[0, j]) < 100]
    results["bst"] = {
        "lambda_inf": bst_diagonal[-1] if bst_diagonal else None,
        "diagonal": bst_diagonal,
    }

    # --- Consensus ---
    estimates = []
    for method in ["powerlaw", "aitken", "pade", "bst"]:
        if method in results and results[method].get("lambda_inf") is not None:
            val = results[method]["lambda_inf"]
            if isinstance(val, (int, float)) and 0.5 < val < 3.0:
                estimates.append((method, val))

    # Add Koopman estimates
    for key, val in koopman_results.items():
        if "lambda_inf" in val and isinstance(val["lambda_inf"], (int, float)):
            if 0.5 < val["lambda_inf"] < 3.0:
                estimates.append((f"koopman_{key}", val["lambda_inf"]))

    if estimates:
        values = [v for _, v in estimates]
        results["consensus"] = {
            "methods": {m: v for m, v in estimates},
            "mean": np.mean(values),
            "std": np.std(values),
            "median": np.median(values),
            "range": [min(values), max(values)],
        }

    return results


# ============================================================
# 5. MAIN
# ============================================================

def main():
    print("=" * 80)
    print("EXP-009: Takens Embedding + KvN Holographic Extrapolation")
    print("         'Let geometry accelerate the compute'")
    print("=" * 80)

    # Load data
    sidon = load_sidon_full()
    sumfree = load_sumfree()

    Ws_sidon = sorted(sidon.keys())
    lambdas_sidon = np.array([sidon[w] for w in Ws_sidon])

    print(f"\n  Sidon data: W = {Ws_sidon[0]}..{Ws_sidon[-1]} ({len(Ws_sidon)} points)")

    # ---- SECTION 1: Takens embedding dimension ----
    print("\n" + "-" * 60)
    print("1. TAKENS EMBEDDING ANALYSIS")
    print("-" * 60)

    # Use the monotonically decreasing portion (skip early W where λ isn't monotone)
    # Actually use all data — the non-monotonicity IS part of the dynamics
    series = lambdas_sidon

    fnn = estimate_embedding_dimension(series, max_dim=7)
    print(f"\n  False Nearest Neighbors analysis:")
    for r in fnn:
        marker = " ← optimal" if r["fnn_fraction"] < 0.05 and r["dim"] > 1 else ""
        print(f"    dim={r['dim']}: FNN fraction = {r['fnn_fraction']:.3f} ({r['n_points']} pts){marker}")

    # Find optimal dimension
    optimal_dim = 3  # default
    for r in fnn:
        if r["dim"] >= 2 and r["fnn_fraction"] < 0.05:
            optimal_dim = r["dim"]
            break

    print(f"\n  Selected embedding dimension: {optimal_dim}")

    # ---- SECTION 2: Torus geometry ----
    print("\n" + "-" * 60)
    print("2. TORUS GEOMETRY / WINDING NUMBERS")
    print("-" * 60)

    torus = torus_analysis(series)
    print(f"\n  Convergence type: {torus['convergence_type']}")
    print(f"  Difference ratio: mean={torus['diff_ratios_mean']:.4f}, std={torus['diff_ratios_std']:.4f}")

    if torus["dominant_frequencies"]:
        print(f"  Dominant frequencies: {[f'{f:.4f}' for f in torus['dominant_frequencies'][:3]]}")
        print(f"  Winding numbers: {[f'{w:.4f}' for w in torus['winding_numbers'][:3]]}")

    if torus["aitken_limit"] is not None:
        print(f"\n  Aitken Δ² accelerated limit: λ_∞ ≈ {torus['aitken_limit']:.6f}")
        if torus["aitken_estimates"]:
            print(f"  Last 5 Aitken estimates: {[f'{a:.6f}' for a in torus['aitken_estimates'][-5:]]}")

    # ---- SECTION 3: Koopman/DMD ----
    print("\n" + "-" * 60)
    print("3. KOOPMAN/DMD EXTRAPOLATION")
    print("-" * 60)

    for dim in [3, 4, 5, 6]:
        if len(series) > dim + 2:
            try:
                koop = koopman_extrapolate(series, dim, n_future=100)
                koop_eigs = koop["koopman_eigenvalues"]
                print(f"\n  Embedding dim={dim}:")
                print(f"    Koopman eigenvalues: {[f'{abs(e):.4f}∠{np.angle(e)*180/np.pi:.1f}°' for e in koop_eigs[:dim]]}")
                print(f"    Predicted λ(W→∞): {koop['predictions'][-1]:.6f}")

                # Show convergence of predictions
                steps = [10, 25, 50, 100]
                preds_at = [koop["predictions"][min(s-1, len(koop["predictions"])-1)] for s in steps]
                print(f"    Predictions at +{steps}: {[f'{p:.6f}' for p in preds_at]}")

                # Check: are Koopman eigenvalues inside unit circle? (stability)
                stable = all(abs(e) <= 1.0 + 1e-10 for e in koop_eigs)
                print(f"    Koopman stable: {stable}")

            except Exception as e:
                print(f"\n  Embedding dim={dim}: ERROR — {e}")

    # ---- SECTION 4: Holographic consensus ----
    print("\n" + "-" * 60)
    print("4. HOLOGRAPHIC CONSENSUS (all methods)")
    print("-" * 60)

    holo = holographic_extrapolate(series, Ws_sidon)

    if "powerlaw" in holo:
        p = holo["powerlaw"]
        print(f"\n  Power-law: {p['formula']}")
        print(f"    → λ_∞ = {p['lambda_inf']:.6f} (assumes convergence to 1)")

    if "aitken" in holo:
        a = holo["aitken"]
        print(f"\n  Aitken Δ²: λ_∞ ≈ {a['lambda_inf']:.6f}" if a["lambda_inf"] else "\n  Aitken: insufficient data")
        print(f"    Convergence: {a['convergence_type']}")

    if "pade" in holo and "lambda_inf" in holo["pade"]:
        print(f"\n  Padé [2/2]: λ_∞ ≈ {holo['pade']['lambda_inf']:.6f}")

    if "bst" in holo and holo["bst"]["lambda_inf"] is not None:
        print(f"\n  BST rational: λ_∞ ≈ {holo['bst']['lambda_inf']:.6f}")
        if holo["bst"]["diagonal"]:
            print(f"    Diagonal convergence: {[f'{d:.4f}' for d in holo['bst']['diagonal'][-5:]]}")

    if "koopman" in holo:
        for key, val in holo["koopman"].items():
            if "lambda_inf" in val:
                print(f"\n  Koopman ({key}): λ_∞ ≈ {val['lambda_inf']:.6f}")

    if "consensus" in holo:
        c = holo["consensus"]
        print(f"\n  {'='*50}")
        print(f"  CONSENSUS: λ_∞ = {c['mean']:.6f} ± {c['std']:.6f}")
        print(f"  Median: {c['median']:.6f}")
        print(f"  Range: [{c['range'][0]:.6f}, {c['range'][1]:.6f}]")
        print(f"  Methods: {len(c['methods'])}")
        for m, v in sorted(c['methods'].items()):
            print(f"    {m:25s}: {v:.6f}")

    # ---- SECTION 5: Calibration against sum-free ----
    if sumfree:
        print("\n" + "-" * 60)
        print("5. CALIBRATION: Sum-Free (known answer: λ_∞ = 1)")
        print("-" * 60)

        Ws_sf = sorted(sumfree.keys())
        lambdas_sf = np.array([sumfree[w] for w in Ws_sf])

        holo_sf = holographic_extrapolate(lambdas_sf, Ws_sf)

        if "consensus" in holo_sf:
            c_sf = holo_sf["consensus"]
            print(f"\n  CONSENSUS: λ_∞ = {c_sf['mean']:.6f} ± {c_sf['std']:.6f}")
            print(f"  Known answer: 1.000000")
            print(f"  Error: {abs(c_sf['mean'] - 1.0):.6f}")

            for m, v in sorted(c_sf['methods'].items()):
                err = abs(v - 1.0)
                quality = "✓" if err < 0.05 else "✗"
                print(f"    {m:25s}: {v:.6f}  (err={err:.4f}) {quality}")

            # Which method is best on the calibration problem?
            best_method = min(c_sf['methods'].items(), key=lambda x: abs(x[1] - 1.0))
            print(f"\n  Best method on calibration: {best_method[0]} (err={abs(best_method[1]-1.0):.6f})")

            # Use that method's Sidon estimate
            if best_method[0] in holo.get("consensus", {}).get("methods", {}):
                best_sidon = holo["consensus"]["methods"][best_method[0]]
                print(f"  → That method's Sidon estimate: λ_∞ = {best_sidon:.6f}")

    # ---- SECTION 6: Implications for h(N) ----
    print("\n" + "-" * 60)
    print("6. IMPLICATIONS FOR h(N) = √N + correction")
    print("-" * 60)

    if "consensus" in holo:
        lambda_inf = holo["consensus"]["median"]

        if lambda_inf > 1.001:
            # λ_∞ > 1 means the Sidon lattice gas has positive free energy
            # at the thermodynamic limit → the constraint allows exponential
            # growth in the number of states per site
            fe = math.log(lambda_inf)
            print(f"\n  λ_∞ = {lambda_inf:.6f}  →  free energy/site = {fe:.6f}")
            print(f"  This means: the transfer matrix does NOT reach criticality")
            print(f"  Interpretation: the window-W formulation misses long-range")
            print(f"  correlations that enforce the √N constraint.")
            print(f"  The PMF needs a renormalization step (DMRG/tensor network)")
            print(f"  to capture these correlations.")
        elif lambda_inf < 0.999:
            print(f"\n  λ_∞ = {lambda_inf:.6f} < 1  →  the system is subcritical")
            print(f"  This would be surprising. Check data quality.")
        else:
            print(f"\n  λ_∞ ≈ 1.000  →  the system reaches criticality!")
            print(f"  The correction term is governed by how fast λ → 1.")

            if "powerlaw" in holo:
                alpha = holo["powerlaw"]["alpha"]
                correction_exp = 1.0 / (2 * alpha)
                print(f"  Power-law α = {alpha:.4f}")
                print(f"  → h(N) = √N + O(N^{{{correction_exp:.4f}}})")

                if correction_exp < 0.25:
                    print(f"  THIS BEATS LINDSTRÖM (N^0.25)!")
                elif correction_exp < 0.5:
                    print(f"  Better than trivial (N^0.5) but not SOTA.")
                else:
                    print(f"  Does not beat trivial bound.")

    # ---- Save results ----
    out = {
        "experiment_id": "EXP-009",
        "title": "Takens Embedding + KvN Holographic Extrapolation",
        "method": "Geometry-accelerated finite-size scaling",
        "sidon_data_range": f"W={Ws_sidon[0]}..{Ws_sidon[-1]}",
        "takens_optimal_dim": optimal_dim,
        "torus_analysis": {
            "convergence_type": torus["convergence_type"],
            "diff_ratio_mean": torus["diff_ratios_mean"],
            "aitken_limit": torus["aitken_limit"],
        },
        "holographic_consensus": holo.get("consensus"),
    }

    # Add calibration
    if sumfree and "consensus" in locals().get("holo_sf", {}):
        out["calibration_sumfree"] = {
            "known_answer": 1.0,
            "consensus_mean": holo_sf["consensus"]["mean"],
            "consensus_std": holo_sf["consensus"]["std"],
            "error": abs(holo_sf["consensus"]["mean"] - 1.0),
        }

    with open("EXP-009_TAKENS_KVN_RESULTS.json", "w") as f:
        json.dump(out, f, indent=2, default=lambda x: float(x) if isinstance(x, np.floating) else None)
    print(f"\nSaved to EXP-009_TAKENS_KVN_RESULTS.json")


if __name__ == "__main__":
    main()
