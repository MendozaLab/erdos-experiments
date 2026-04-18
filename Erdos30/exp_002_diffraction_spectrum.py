#!/usr/bin/env python3
"""
EXP-MATH-ERDOS30-SIDON-002: Diffraction Spectrum of Singer Sets

Leg 4 experiment for the PHYS-QC-001 (Quasicrystal/Meyer Set) morphism.

Question: Do Singer perfect difference sets exhibit quasicrystalline
diffraction (pure-point spectrum), and does the continuous spectral
residual decay as q → ∞?

If yes: strongest evidence that near-extremal Sidon sets are model sets,
which would imply h(N) = √N + O_ε(N^ε) via model set diffraction theory.

If no: the quasicrystal morphism fails Leg 4, and O(N^{1/4}) may be tight.
"""

import json
import math
import time
from pathlib import Path

import numpy as np

# ─── Singer PDS Construction ─────────────────────────────────────────────

def singer_pds(q: int) -> list[int]:
    """Construct Singer perfect difference set for prime power q.
    Returns a Sidon set of size q+1 in Z_{q²+q+1}."""
    n = q * q + q + 1
    # For prime q: use primitive root of GF(q³) mod the ideal
    # Simplified: find a set S ⊂ Z_n with |S|=q+1 where all
    # differences are distinct mod n
    if q == 2:
        return [0, 1, 3]
    elif q == 3:
        return [0, 1, 3, 9]
    elif q == 4:
        return [0, 1, 4, 14, 16]
    elif q == 5:
        return [0, 1, 3, 8, 12, 18]
    elif q == 7:
        return [0, 1, 3, 13, 32, 36, 43, 52]
    elif q == 8:
        return [0, 1, 3, 7, 15, 31, 36, 54, 63]
    elif q == 9:
        return [0, 1, 3, 9, 27, 49, 56, 61, 77, 81]
    elif q == 11:
        return [0, 1, 3, 12, 20, 34, 38, 81, 88, 94, 104, 109]
    elif q == 13:
        return [0, 1, 3, 16, 23, 28, 42, 76, 82, 86, 119, 137, 154, 175]
    else:
        # Greedy Sidon construction as fallback for larger q
        return greedy_sidon(q * q + q + 1, q + 1)


def greedy_sidon(n: int, target_size: int) -> list[int]:
    """Greedy Sidon set construction in Z_n."""
    A = [0]
    sums = {0}
    for x in range(1, n):
        new_sums = set()
        valid = True
        for a in A:
            s = (x + a) % (2 * n)  # track all pairwise sums
            if s in sums or s in new_sums:
                valid = False
                break
            new_sums.add(s)
        if valid:
            A.append(x)
            sums.update(new_sums)
            if len(A) >= target_size:
                break
    return A


# ─── Diffraction Measure Computation ────────────────────────────────────

def compute_diffraction(A: list[int], n: int, num_freq: int = 2048) -> dict:
    """
    Compute the diffraction measure of set A ⊂ Z_n.

    The diffraction measure is: γ̂(ξ) = (1/|A|) |Σ_{a∈A} e^{2πiaξ/n}|²

    Returns dict with:
      - frequencies: array of ξ values
      - spectrum: γ̂(ξ) values
      - bragg_peaks: identified Bragg peak positions and heights
      - bragg_fraction: fraction of spectral mass in Bragg peaks
      - continuous_residual: 1 - bragg_fraction
    """
    A_arr = np.array(A, dtype=float)
    k = len(A)

    # Compute diffraction at integer frequencies (these are the natural Bragg positions for Z_n)
    bragg_spectrum = np.zeros(n)
    for m in range(n):
        # Exponential sum
        phases = 2 * np.pi * A_arr * m / n
        re = np.sum(np.cos(phases))
        im = np.sum(np.sin(phases))
        bragg_spectrum[m] = (re**2 + im**2) / k  # normalized

    # Total spectral mass (should equal k by Parseval)
    total_mass = np.sum(bragg_spectrum) / n

    # Identify Bragg peaks: values significantly above the mean
    mean_val = total_mass  # average value of γ̂ over all frequencies
    threshold = 3 * mean_val  # peak = 3× above mean

    peak_indices = np.where(bragg_spectrum > threshold * n / k)[0]
    peak_mass = np.sum(bragg_spectrum[peak_indices]) / n

    bragg_fraction = peak_mass / total_mass if total_mass > 0 else 0

    # Also compute the diffraction at high resolution (continuous part)
    xi = np.linspace(0, 1, num_freq, endpoint=False)
    continuous_spectrum = np.zeros(num_freq)
    for i, x in enumerate(xi):
        phases = 2 * np.pi * A_arr * x
        re = np.sum(np.cos(phases))
        im = np.sum(np.sin(phases))
        continuous_spectrum[i] = (re**2 + im**2) / k

    # L² norm of the continuous part (subtract Bragg peaks)
    # The continuous spectral density
    cont_l2 = np.mean(continuous_spectrum**2)
    bragg_l2 = np.mean(bragg_spectrum**2) / n  # normalize

    return {
        "bragg_spectrum": bragg_spectrum.tolist(),
        "bragg_fraction": float(bragg_fraction),
        "continuous_residual": float(1 - bragg_fraction),
        "total_mass": float(total_mass),
        "num_peaks": int(len(peak_indices)),
        "peak_indices": peak_indices.tolist()[:20],  # first 20
        "max_peak": float(np.max(bragg_spectrum)),
        "mean_spectrum": float(np.mean(bragg_spectrum)),
        "spectral_flatness_db": float(
            10 * np.log10(np.exp(np.mean(np.log(bragg_spectrum + 1e-15))) /
                          (np.mean(bragg_spectrum) + 1e-15))
        ),
    }


# ─── Chowla Cosine Maximum ──────────────────────────────────────────────

def chowla_cosine_max(A: list[int], n: int) -> float:
    """Compute max_t |Σ_{a∈A} cos(2πat/n)| for t = 1, ..., n-1."""
    A_arr = np.array(A, dtype=float)
    max_val = 0.0
    for t in range(1, n):
        val = abs(np.sum(np.cos(2 * np.pi * A_arr * t / n)))
        max_val = max(max_val, val)
    return float(max_val)


# ─── Deviation Process Entropy ───────────────────────────────────────────

def deviation_entropy(A: list[int], n: int) -> float:
    """Compute Shannon entropy of the discretized deviation process X(t) = |A∩[1,t]| - √t."""
    A_set = set(A)
    count = 0
    deviations = []
    for t in range(1, n + 1):
        if t in A_set:
            count += 1
        dev = count - math.sqrt(t)
        deviations.append(round(dev, 1))  # discretize to 0.1

    # Empirical distribution
    from collections import Counter
    counts = Counter(deviations)
    total = len(deviations)
    entropy = 0.0
    for c in counts.values():
        p = c / total
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy


# ─── Main Experiment ─────────────────────────────────────────────────────

def run_experiment():
    primes = [2, 3, 5, 7, 11, 13]
    # Also use greedy for larger sizes
    large_primes = [17, 19, 23, 29, 31, 37, 41, 43, 47]

    results = {
        "experiment_id": "EXP-MATH-ERDOS30-SIDON-002",
        "title": "Diffraction Spectrum of Singer Sets — Quasicrystal Leg 4",
        "morphism": "PHYS-QC-001 (Quasicrystal / Meyer Set)",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "singer_results": [],
        "greedy_results": [],
        "summary": {},
    }

    print("=" * 70)
    print("EXP-MATH-ERDOS30-SIDON-002: Diffraction Spectrum of Singer Sets")
    print("=" * 70)

    # ── Singer sets (exact construction) ──
    print("\n--- Singer Perfect Difference Sets ---")
    print(f"{'q':>4s} | {'N':>7s} | {'k':>4s} | {'Bragg%':>7s} | {'R(q)':>8s} | {'Peaks':>5s} | {'Chowla':>8s} | {'H(X)':>6s}")
    print("-" * 70)

    bragg_fractions = []
    residuals = []
    chowla_vals = []
    entropy_vals = []

    for q in primes:
        n = q * q + q + 1
        A = singer_pds(q)
        k = len(A)

        diff = compute_diffraction(A, n)
        chowla = chowla_cosine_max(A, n)
        entropy = deviation_entropy(A, n)

        bf = diff["bragg_fraction"]
        r = diff["continuous_residual"]
        bragg_fractions.append(bf)
        residuals.append(r)
        chowla_vals.append(chowla)
        entropy_vals.append(entropy)

        print(f"{q:4d} | {n:7d} | {k:4d} | {bf*100:6.2f}% | {r:8.5f} | {diff['num_peaks']:5d} | {chowla:8.4f} | {entropy:6.2f}")

        results["singer_results"].append({
            "q": q,
            "n": n,
            "k": k,
            "bragg_fraction": bf,
            "continuous_residual": r,
            "num_peaks": diff["num_peaks"],
            "chowla_cosine_max": chowla,
            "chowla_normalized": chowla / math.sqrt(k * math.log(n + 1)),
            "deviation_entropy": entropy,
            "spectral_flatness_db": diff["spectral_flatness_db"],
            "total_mass": diff["total_mass"],
        })

    # ── Greedy sets for comparison ──
    print("\n--- Greedy Sidon Sets (non-algebraic control) ---")
    print(f"{'q':>4s} | {'N':>7s} | {'k':>4s} | {'Bragg%':>7s} | {'R(q)':>8s} | {'Peaks':>5s} | {'Chowla':>8s} | {'H(X)':>6s}")
    print("-" * 70)

    for q in [2, 3, 5, 7, 11]:
        n = q * q + q + 1
        A = greedy_sidon(n, q + 1)
        k = len(A)

        diff = compute_diffraction(A, n)
        chowla = chowla_cosine_max(A, n)
        entropy = deviation_entropy(A, n)

        bf = diff["bragg_fraction"]
        r = diff["continuous_residual"]

        print(f"{q:4d} | {n:7d} | {k:4d} | {bf*100:6.2f}% | {r:8.5f} | {diff['num_peaks']:5d} | {chowla:8.4f} | {entropy:6.2f}")

        results["greedy_results"].append({
            "q": q,
            "n": n,
            "k": k,
            "bragg_fraction": bf,
            "continuous_residual": r,
            "num_peaks": diff["num_peaks"],
            "chowla_cosine_max": chowla,
            "chowla_normalized": chowla / math.sqrt(k * math.log(n + 1)),
            "deviation_entropy": entropy,
            "spectral_flatness_db": diff["spectral_flatness_db"],
        })

    # ── Analysis ──
    print("\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)

    # Trend in continuous residual
    if len(residuals) >= 3:
        # Fit R(q) = C * q^(-alpha)
        log_q = [math.log(primes[i]) for i in range(len(residuals)) if residuals[i] > 1e-10]
        log_r = [math.log(residuals[i]) for i in range(len(residuals)) if residuals[i] > 1e-10]

        if len(log_q) >= 2:
            # Linear regression
            n_pts = len(log_q)
            sx = sum(log_q)
            sy = sum(log_r)
            sxx = sum(x*x for x in log_q)
            sxy = sum(x*y for x, y in zip(log_q, log_r))

            denom = n_pts * sxx - sx * sx
            if abs(denom) > 1e-15:
                alpha = -(n_pts * sxy - sx * sy) / denom
                log_C = (sy + alpha * sx) / n_pts
                C = math.exp(log_C)

                # R²
                y_mean = sy / n_pts
                ss_tot = sum((y - y_mean)**2 for y in log_r)
                ss_res = sum((log_r[i] - (log_C - alpha * log_q[i]))**2 for i in range(n_pts))
                r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0

                print(f"\nContinuous residual decay: R(q) ≈ {C:.4f} × q^(-{alpha:.4f})")
                print(f"R² = {r_sq:.4f}")

                results["summary"]["residual_decay_exponent"] = alpha
                results["summary"]["residual_decay_coefficient"] = C
                results["summary"]["residual_decay_r_squared"] = r_sq

    # Chowla trend
    print(f"\nChowla cosine maxima (Singer): {[f'{v:.3f}' for v in chowla_vals]}")
    chowla_norms = [results["singer_results"][i]["chowla_normalized"] for i in range(len(primes))]
    print(f"Chowla normalized (÷√(k·ln(n))): {[f'{v:.3f}' for v in chowla_norms]}")

    if chowla_norms[-1] < chowla_norms[0]:
        print("→ TREND: Chowla normalized ratio DECREASING — Singer sets are anomalously flat")
    else:
        print("→ TREND: Chowla normalized ratio NOT decreasing — spectral flatness inconclusive")

    # Entropy trend
    print(f"\nDeviation entropy (Singer): {[f'{v:.2f}' for v in entropy_vals]}")

    # Verdict
    print("\n" + "=" * 70)
    print("LEG 4 VERDICT: PHYS-QC-001 (Quasicrystal / Meyer Set)")
    print("=" * 70)

    final_bf = bragg_fractions[-1] if bragg_fractions else 0
    if final_bf > 0.90:
        verdict = "PASS"
        msg = f"Bragg fraction {final_bf*100:.1f}% → consistent with pure-point diffraction"
    elif final_bf > 0.70:
        verdict = "PARTIAL"
        msg = f"Bragg fraction {final_bf*100:.1f}% → trending toward pure-point but not definitive"
    else:
        verdict = "FAIL"
        msg = f"Bragg fraction {final_bf*100:.1f}% → significant continuous spectral component"

    print(f"\nVerdict: {verdict}")
    print(f"Reason: {msg}")

    alpha_val = results["summary"].get("residual_decay_exponent", 0)
    if alpha_val > 0:
        print(f"Residual decay: R(q) ~ q^(-{alpha_val:.2f})")
        if alpha_val > 0.5:
            print("→ Fast decay: consistent with subpolynomial error term")
        else:
            print("→ Slow decay: consistent with polynomial error (O(N^{1/4}) may be tight)")

    results["summary"]["verdict"] = verdict
    results["summary"]["bragg_fractions"] = bragg_fractions
    results["summary"]["final_bragg_fraction"] = final_bf

    # Save
    out_path = Path(__file__).parent / "EXP-MATH-ERDOS30-SIDON-002_RESULTS.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")

    return results


if __name__ == "__main__":
    run_experiment()
