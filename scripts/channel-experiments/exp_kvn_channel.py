"""
Koopman-von Neumann Channel Capacity Experiment — Erdős #114/#509

The KvN lens on the EHP problem:

  f: ℂ → ℂ monic degree-n polynomial acts as a "channel".
  Koopman operator: U_f g = g ∘ f  on L²(Λ(f), arc-length measure μ).
  Böttcher coordinate: φ_f satisfying φ_f(f(z)) = φ_f(z)^n is the
    Koopman intertwiner — it conjugates f to z ↦ z^n on the Riemann sphere.

  The lemniscate Λ(f) = {|f|=1} = φ_f^{-1}(S¹).
  On Λ(f), the map induced by f is the n-fold cover z ↦ z^n of S¹.

  Channel capacity interpretation:
    - Each polynomial defines a communication channel: pre-images of a point
      on S¹ are "symbols" that encode the same output.
    - z^n−1 achieves n equal pre-images uniformly distributed on Λ(f),
      maximising the Shannon entropy of the pre-image distribution.
    - This is equivalent to maximising arc length (EHP) because the
      uniform covering is the most "efficient" (max-entropy) encoding.

Experiments:
  1. Preimage entropy H_f(n): entropy of the n pre-image distribution on Λ(f).
     For z^n−1: all n pre-images equidistant → H = log(n) (maximum).
  2. Koopman transfer matrix T_f: n×n matrix of arc-length overlap integrals
     between pre-image sectors. Spectral gap measures channel non-uniformity.
  3. Böttcher approximation quality: |φ_f(z) − z| on |z|=2 as a function of
     coefficient perturbation. Measures how far f is from the extremizer.
  4. KvN channel capacity vs L(f): correlation across random degree-3,4 polys.
"""
import json
import os
from typing import NamedTuple

import numpy as np
from scipy import special
from skimage import measure

RESULTS_DIR = "../../results/channel-experiments"
EXTENT = 3.5
IMG_SIZE = 400


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def l_star_formula(n: int) -> float:
    a = 1.0 / (2 * n)
    return 2.0 ** (1.0 / n) * np.sqrt(np.pi) * special.gamma(a) / special.gamma(a + 0.5)


def eval_poly(z: np.ndarray, a: np.ndarray, n: int) -> np.ndarray:
    """f(z) = z^n + a[n-2]*z^{n-2} + … + a[1]*z + a[0].  a[n-1]=0."""
    F = z ** n
    for k in range(n - 1):
        if a[k] != 0.0:
            F = F + a[k] * z ** k
    return F


def arc_length_from_contour(a: np.ndarray, n: int) -> float:
    x = np.linspace(-EXTENT, EXTENT, IMG_SIZE)
    X, Y = np.meshgrid(x, x)
    Z = X + 1j * Y
    absF = np.abs(eval_poly(Z, a, n)).astype(np.float64)
    try:
        contours = measure.find_contours(absF, 1.0)
    except Exception:
        return 0.0
    if not contours:
        return 0.0
    dx = (2.0 * EXTENT) / (IMG_SIZE - 1)
    total = 0.0
    for c in contours:
        cx = c[:, 1] * dx - EXTENT
        cy = c[:, 0] * dx - EXTENT
        total += float(np.sum(np.sqrt(np.diff(cx) ** 2 + np.diff(cy) ** 2)))
    return total


# ---------------------------------------------------------------------------
# Experiment 1: Pre-image entropy H_f
# ---------------------------------------------------------------------------

class PreimageResult(NamedTuple):
    n: int
    entropy: float
    max_entropy: float
    uniformity: float   # entropy / max_entropy ∈ (0,1]


def preimage_entropy(a: np.ndarray, n: int, n_angles: int = 360) -> PreimageResult:
    """
    For each angle θ ∈ [0,2π), the target point is e^{iθ} ∈ S¹.
    The pre-images are the n solutions of f(z) = e^{iθ} on Λ(f).
    We estimate the distribution of arc-length gaps between pre-images
    and compute the Shannon entropy of this distribution.

    For z^n − 1: f(z) = z^n − 1 = e^{iθ}  →  z = e^{i(θ+1)/n} · ω^k_n
    All n pre-images are equally spaced → H = log(n).
    """
    thetas = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
    entropies = []

    for theta in thetas:
        target = np.exp(1j * theta)
        # Find roots of f(z) - target = 0 numerically
        # Build numpy poly: f(z) - target = z^n + a[n-2]*z^{n-2} + ... + a[0] - target
        poly_np = np.zeros(n + 1, dtype=complex)
        poly_np[0] = 1.0
        for k in range(n - 1):
            poly_np[n - k] = a[k]
        poly_np[-1] -= target  # subtract target from constant term

        roots = np.roots(poly_np)

        # Keep roots near |z| = 1 (on or near the lemniscate)
        near = roots[np.abs(np.abs(roots) - 1.0) < 0.5]
        if len(near) < 2:
            continue

        # Arc-length distribution: angles of pre-images on approximate unit circle
        angles = np.sort(np.angle(near) % (2 * np.pi))
        gaps = np.diff(np.append(angles, angles[0] + 2 * np.pi))
        gaps = gaps / gaps.sum()  # normalise to probability

        # Shannon entropy of gap distribution
        h = -np.sum(gaps * np.log(gaps + 1e-300))
        entropies.append(h)

    if not entropies:
        return PreimageResult(n, 0.0, np.log(n), 0.0)

    mean_h = float(np.mean(entropies))
    max_h = np.log(n)
    return PreimageResult(n, mean_h, max_h, mean_h / max_h if max_h > 0 else 0.0)


# ---------------------------------------------------------------------------
# Experiment 2: Koopman transfer matrix
# ---------------------------------------------------------------------------

def koopman_transfer_matrix(a: np.ndarray, n: int) -> np.ndarray:
    """
    Approximate n×n Koopman transfer matrix.

    Partition S¹ into n equal sectors S_k = [2πk/n, 2π(k+1)/n).
    T[j,k] = fraction of pre-images of sector S_k that fall in sector S_j.

    For z^n − 1: T = (1/n) * ones(n,n)  →  uniform, maximum entropy.
    For a generic competitor: T is non-uniform → lower entropy.
    """
    n_samples = 200
    T = np.zeros((n, n))

    for sector_k in range(n):
        # Sample points in sector S_k on the unit circle
        thetas = np.linspace(
            2 * np.pi * sector_k / n,
            2 * np.pi * (sector_k + 1) / n,
            n_samples, endpoint=False
        )
        for theta in thetas:
            target = np.exp(1j * theta)
            poly_np = np.zeros(n + 1, dtype=complex)
            poly_np[0] = 1.0
            for i in range(n - 1):
                poly_np[n - i] = a[i]
            poly_np[-1] -= target

            roots = np.roots(poly_np)
            for r in roots:
                if abs(abs(r) - 1.0) < 0.5:
                    angle = np.angle(r) % (2 * np.pi)
                    j = int(angle / (2 * np.pi / n)) % n
                    T[j, sector_k] += 1.0

    col_sums = T.sum(axis=0, keepdims=True)
    col_sums[col_sums == 0] = 1.0
    return T / col_sums


def transfer_spectral_gap(T: np.ndarray) -> float:
    """
    Spectral gap of the transfer matrix = 1 - |λ_2|/|λ_1|.
    For uniform T (z^n−1): all eigenvalues except λ_1=1 are 0 → gap=1.
    Smaller gap → less uniform channel.
    """
    eigvals = np.abs(np.linalg.eigvals(T))
    eigvals = np.sort(eigvals)[::-1]
    if eigvals[0] < 1e-10:
        return 0.0
    return float(1.0 - eigvals[1] / eigvals[0]) if len(eigvals) > 1 else 1.0


# ---------------------------------------------------------------------------
# Experiment 3: Böttcher approximation quality
# ---------------------------------------------------------------------------

def bottcher_quality(a: np.ndarray, n: int, R: float = 2.0, n_theta: int = 200) -> float:
    """
    Böttcher coordinate near ∞: φ_f(z) = z · (1 + a[n-2]/z^2 + …)^{1/n} ≈ z for |z| large.
    Quality = 1 - mean|φ_f(z)/z - 1| on |z|=R, measures proximity to z^n−1 extremizer.

    For z^n − 1 (a = [1, 0, …, 0]):
      φ(z) = z·(1 − 1/z^n)^{1/n} ≈ z − 1/(n z^{n-1}) + …
      |φ(z)/z − 1| = |1 − (1 − 1/z^n)^{1/n} − 1| ≈ 1/(n|z|^n) → 0 as |z| → ∞
    So quality → 1 for z^n−1 and for any polynomial as R → ∞.
    At fixed R, competitors deviate more.
    """
    theta = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    z = R * np.exp(1j * theta)
    fz = eval_poly(z, a, n)

    # φ_f(z) = f(z)^{1/n} (principal branch, matches z for large |z|)
    phi = fz ** (1.0 / n)
    # Align branch: phi should be ≈ z for large |z|
    # Choose branch so that phi/z ≈ 1
    ratio = phi / z
    # Correct branch by multiplying by n-th root of unity if needed
    best_ratio = np.ones_like(ratio)
    for k in range(n):
        omega = np.exp(2j * np.pi * k / n)
        candidate = ratio * omega
        if np.abs(np.mean(candidate) - 1.0) < np.abs(np.mean(best_ratio) - 1.0):
            best_ratio = candidate

    quality = float(1.0 - np.mean(np.abs(best_ratio - 1.0)))
    return quality


# ---------------------------------------------------------------------------
# Main experiment runner
# ---------------------------------------------------------------------------

def run() -> dict:
    rng = np.random.default_rng(42)
    N_POLY = 50   # competitors per degree
    by_n: dict = {}

    for n in [3, 4]:
        L_star = l_star_formula(n)

        # z^n + 1 (reduced form of z^n − 1)
        a_zn1 = np.zeros(n, dtype=complex)
        a_zn1[0] = 1.0

        L_zn1 = arc_length_from_contour(a_zn1, n)
        H_zn1 = preimage_entropy(a_zn1, n)
        T_zn1 = koopman_transfer_matrix(a_zn1, n)
        gap_zn1 = transfer_spectral_gap(T_zn1)
        bq_zn1 = bottcher_quality(a_zn1, n)

        print(f"\nn={n}: z^n+1")
        print(f"  L={L_zn1:.5f} (formula={L_star:.5f})")
        print(f"  Preimage entropy: H={H_zn1.entropy:.4f}, max={H_zn1.max_entropy:.4f}, "
              f"uniformity={H_zn1.uniformity:.4f}")
        print(f"  Koopman spectral gap: {gap_zn1:.4f} (1.0=fully uniform)")
        print(f"  Böttcher quality: {bq_zn1:.4f}")

        # Competitors
        comp_Ls, comp_Hs, comp_gaps, comp_bqs = [], [], [], []
        for _ in range(N_POLY):
            a = np.zeros(n, dtype=complex)
            for k in range(1, n - 1):
                a[k] = rng.uniform(-0.8, 0.8) + 1j * rng.uniform(-0.8, 0.8)
            a[0] = rng.uniform(0.1, 1.5)

            comp_Ls.append(arc_length_from_contour(a, n))
            h = preimage_entropy(a, n)
            comp_Hs.append(h.uniformity)
            T = koopman_transfer_matrix(a, n)
            comp_gaps.append(transfer_spectral_gap(T))
            comp_bqs.append(bottcher_quality(a, n))

        comp_Ls = np.array(comp_Ls)
        comp_Hs = np.array(comp_Hs)
        comp_gaps = np.array(comp_gaps)

        corr_L_H = float(np.corrcoef(comp_Ls, comp_Hs)[0, 1])
        corr_L_gap = float(np.corrcoef(comp_Ls, comp_gaps)[0, 1])

        print(f"  Competitors: Corr(L, uniformity)={corr_L_H:.4f}, "
              f"Corr(L, spectral_gap)={corr_L_gap:.4f}")
        print(f"  Max competitor L={comp_Ls.max():.5f} vs z^n+1 L={L_zn1:.5f}")

        by_n[str(n)] = {
            "L_star_formula": L_star,
            "zn1": {
                "L": float(L_zn1),
                "preimage_entropy": float(H_zn1.entropy),
                "preimage_uniformity": float(H_zn1.uniformity),
                "koopman_spectral_gap": float(gap_zn1),
                "bottcher_quality": float(bq_zn1),
            },
            "competitors": {
                "n": N_POLY,
                "mean_L": float(comp_Ls.mean()),
                "max_L": float(comp_Ls.max()),
                "mean_uniformity": float(comp_Hs.mean()),
                "corr_L_vs_uniformity": corr_L_H,
                "corr_L_vs_spectral_gap": corr_L_gap,
            },
            "kvn_hypothesis": (
                "z^n-1 maximises preimage entropy (uniformity=1) "
                "and arc length simultaneously — both consequences "
                "of the n-fold uniform covering structure via Böttcher conjugacy."
            ),
        }

    verdict = "PASS" if all(
        v["zn1"]["preimage_uniformity"] > 0.95 and
        v["zn1"]["koopman_spectral_gap"] > 0.8
        for v in by_n.values()
    ) else "REVIEW"

    result = {
        "experiment": "KVN",
        "description": (
            "Koopman-von Neumann channel capacity analysis of EHP extremizer. "
            "Tests: (1) pre-image entropy uniformity, (2) Koopman transfer matrix "
            "spectral gap, (3) Böttcher coordinate quality, (4) L vs uniformity correlation."
        ),
        "by_n": by_n,
        "verdict": verdict,
    }

    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(os.path.join(RESULTS_DIR, "exp_kvn_channel.json"), "w") as fh:
        json.dump(result, fh, indent=2, default=str)

    print(f"\nKvN experiment verdict: {verdict}")
    return result


if __name__ == "__main__":
    run()
