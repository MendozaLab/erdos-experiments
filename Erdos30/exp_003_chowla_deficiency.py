#!/usr/bin/env python3
"""
EXP-MATH-ERDOS30-SIDON-003: Chowla Cosine Deficiency — Extended Scale

Track B experiment for Erdős #30 orthogonal attack.

Question: Do Singer PDS achieve anomalously low Chowla cosine maxima
compared to random/greedy Sidon sets of the same size?

If yes: Singer sets are "spectrally flat" — their Fourier transforms
avoid peaks. Combined with EXP-002 (CV=0), this confirms the
model-set/quasicrystal characterization.

Key observable: C(S_q) / √(k · ln(n)) should → 0 for Singer,
stay bounded away from 0 for greedy/random.
"""

import json
import math
import time
import random
from pathlib import Path
from collections import Counter

import numpy as np

# ─── Singer PDS Construction ─────────────────────────────────────────────

SINGER_SETS = {
    2: [0, 1, 3],
    3: [0, 1, 3, 9],
    4: [0, 1, 4, 14, 16],
    5: [0, 1, 3, 8, 12, 18],
    7: [0, 1, 3, 13, 32, 36, 43, 52],
    8: [0, 1, 3, 7, 15, 31, 36, 54, 63],
    9: [0, 1, 3, 9, 27, 49, 56, 61, 77, 81],
    11: [0, 1, 3, 12, 20, 34, 38, 81, 88, 94, 104, 109],
    13: [0, 1, 3, 16, 23, 28, 42, 76, 82, 86, 119, 137, 154, 175],
    16: [0, 1, 3, 7, 12, 20, 30, 44, 65, 80, 92, 105, 132, 157, 175, 198, 213],
    17: [0, 1, 3, 7, 15, 31, 63, 90, 116, 127, 136, 181, 194, 204, 233, 238, 255, 271],
}


def greedy_sidon(n: int, target_size: int) -> list[int]:
    """Greedy Sidon set construction in Z_n."""
    A = [0]
    diff_set = set()
    for x in range(1, n):
        diffs = set()
        valid = True
        for a in A:
            d1 = (x - a) % n
            d2 = (a - x) % n
            if d1 in diff_set or d1 in diffs or d2 in diff_set or d2 in diffs:
                valid = False
                break
            diffs.add(d1)
            diffs.add(d2)
        if valid:
            A.append(x)
            diff_set.update(diffs)
            if len(A) >= target_size:
                break
    return A


def random_sidon(n: int, target_size: int, seed: int = 42) -> list[int]:
    """Random Sidon set construction in Z_n (randomized greedy)."""
    rng = random.Random(seed)
    candidates = list(range(n))
    rng.shuffle(candidates)
    A = []
    diff_set = set()
    for x in candidates:
        diffs = set()
        valid = True
        for a in A:
            d1 = (x - a) % n
            d2 = (a - x) % n
            if d1 in diff_set or d1 in diffs or d2 in diff_set or d2 in diffs:
                valid = False
                break
            diffs.add(d1)
            diffs.add(d2)
        if valid:
            A.append(x)
            diff_set.update(diffs)
            if len(A) >= target_size:
                break
    return sorted(A)


# ─── Chowla Cosine Analysis ──────────────────────────────────────────────

def chowla_full_analysis(A: list[int], n: int) -> dict:
    """Full Chowla cosine analysis for set A in Z_n."""
    A_arr = np.array(A, dtype=float)
    k = len(A)

    # Compute |Σ cos(2πat/n)| for t = 1, ..., n-1
    cosine_sums = np.zeros(n - 1)
    for t in range(1, n):
        val = np.sum(np.cos(2 * np.pi * A_arr * t / n))
        cosine_sums[t - 1] = abs(val)

    chowla_max = float(np.max(cosine_sums))
    chowla_mean = float(np.mean(cosine_sums))
    chowla_std = float(np.std(cosine_sums))
    chowla_median = float(np.median(cosine_sums))

    # Theoretical random baseline: E[max] ~ √(k · ln(n))
    random_baseline = math.sqrt(k * math.log(n + 1))
    normalized = chowla_max / random_baseline if random_baseline > 0 else 0

    # L² norm of character sums (relates to additive energy)
    # For Sidon: Σ_t |f̂(t)|² = k(n-1) ≈ kn (Parseval minus DC term)
    l2_norm = float(np.sqrt(np.sum(cosine_sums**2)))

    # L⁴ norm (directly related to additive energy E(A))
    l4_fourth = float(np.sum(cosine_sums**4))

    # Spectral peak distribution — how many frequencies exceed various thresholds
    thresholds = [0.5, 1.0, 1.5, 2.0]
    peak_counts = {}
    for thresh in thresholds:
        peak_counts[f"above_{thresh}xsqrtk"] = int(np.sum(cosine_sums > thresh * math.sqrt(k)))

    # Fourier coefficient variance (key observable from EXP-002)
    # For PDS: |f̂(t)|² = k for all t ≠ 0, so variance = 0
    power_spectrum = np.zeros(n - 1)
    for t in range(1, n):
        phases = 2 * np.pi * A_arr * t / n
        re = np.sum(np.cos(phases))
        im = np.sum(np.sin(phases))
        power_spectrum[t - 1] = re**2 + im**2

    fourier_cv = float(np.std(power_spectrum) / np.mean(power_spectrum)) if np.mean(power_spectrum) > 0 else 0

    return {
        "chowla_max": chowla_max,
        "chowla_mean": chowla_mean,
        "chowla_std": chowla_std,
        "chowla_median": chowla_median,
        "chowla_normalized": normalized,
        "random_baseline_sqrt_k_ln_n": random_baseline,
        "l2_character_sum": l2_norm,
        "l4_fourth_power": l4_fourth,
        "fourier_cv": fourier_cv,
        "peak_counts": peak_counts,
    }


# ─── Main Experiment ─────────────────────────────────────────────────────

def run_experiment():
    results = {
        "experiment_id": "EXP-MATH-ERDOS30-SIDON-003",
        "title": "Chowla Cosine Deficiency — Singer vs Greedy vs Random",
        "morphism_link": "Chowla's Cosine Problem (str=1.000, 3 cross-corpus entries)",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "singer": [],
        "greedy": [],
        "random": [],
        "analysis": {},
    }

    print("=" * 78)
    print("EXP-MATH-ERDOS30-SIDON-003: Chowla Cosine Deficiency")
    print("=" * 78)

    # ── Singer PDS ──
    print("\n--- Singer Perfect Difference Sets ---")
    header = f"{'q':>4s} | {'n':>7s} | {'k':>4s} | {'C(S)':>8s} | {'√(k·ln n)':>10s} | {'Ratio':>7s} | {'CV':>8s} | {'Mean':>8s} | {'Median':>8s}"
    print(header)
    print("-" * len(header))

    singer_ratios = []
    for q in sorted(SINGER_SETS.keys()):
        n = q * q + q + 1
        A = SINGER_SETS[q]
        k = len(A)

        analysis = chowla_full_analysis(A, n)
        singer_ratios.append(analysis["chowla_normalized"])

        print(f"{q:4d} | {n:7d} | {k:4d} | {analysis['chowla_max']:8.4f} | "
              f"{analysis['random_baseline_sqrt_k_ln_n']:10.4f} | "
              f"{analysis['chowla_normalized']:7.4f} | "
              f"{analysis['fourier_cv']:8.6f} | "
              f"{analysis['chowla_mean']:8.4f} | "
              f"{analysis['chowla_median']:8.4f}")

        results["singer"].append({"q": q, "n": n, "k": k, **analysis})

    # ── Greedy Sidon ──
    print("\n--- Greedy Sidon Sets ---")
    print(header)
    print("-" * len(header))

    greedy_ratios = []
    for q in sorted(SINGER_SETS.keys()):
        n = q * q + q + 1
        k = q + 1
        A = greedy_sidon(n, k)
        actual_k = len(A)

        analysis = chowla_full_analysis(A, n)
        greedy_ratios.append(analysis["chowla_normalized"])

        print(f"{q:4d} | {n:7d} | {actual_k:4d} | {analysis['chowla_max']:8.4f} | "
              f"{analysis['random_baseline_sqrt_k_ln_n']:10.4f} | "
              f"{analysis['chowla_normalized']:7.4f} | "
              f"{analysis['fourier_cv']:8.6f} | "
              f"{analysis['chowla_mean']:8.4f} | "
              f"{analysis['chowla_median']:8.4f}")

        results["greedy"].append({"q": q, "n": n, "k": actual_k, **analysis})

    # ── Random Sidon (5 seeds, averaged) ──
    print("\n--- Random Sidon Sets (averaged over 5 seeds) ---")
    print(header)
    print("-" * len(header))

    random_ratios = []
    for q in sorted(SINGER_SETS.keys()):
        n = q * q + q + 1
        k = q + 1
        seed_analyses = []
        for seed in range(5):
            A = random_sidon(n, k, seed=seed * 7 + 13)
            seed_analyses.append(chowla_full_analysis(A, n))

        avg_ratio = np.mean([a["chowla_normalized"] for a in seed_analyses])
        avg_cv = np.mean([a["fourier_cv"] for a in seed_analyses])
        avg_max = np.mean([a["chowla_max"] for a in seed_analyses])
        avg_mean = np.mean([a["chowla_mean"] for a in seed_analyses])
        avg_median = np.mean([a["chowla_median"] for a in seed_analyses])
        baseline = seed_analyses[0]["random_baseline_sqrt_k_ln_n"]
        random_ratios.append(float(avg_ratio))

        print(f"{q:4d} | {n:7d} | {k:4d} | {avg_max:8.4f} | "
              f"{baseline:10.4f} | "
              f"{avg_ratio:7.4f} | "
              f"{avg_cv:8.6f} | "
              f"{avg_mean:8.4f} | "
              f"{avg_median:8.4f}")

        results["random"].append({
            "q": q, "n": n, "k": k,
            "chowla_normalized_avg": float(avg_ratio),
            "fourier_cv_avg": float(avg_cv),
            "chowla_max_avg": float(avg_max),
            "num_seeds": 5,
        })

    # ── Trend Analysis ──
    print("\n" + "=" * 78)
    print("TREND ANALYSIS")
    print("=" * 78)

    qs = sorted(SINGER_SETS.keys())

    print(f"\nChowla normalized ratio C(S)/√(k·ln n):")
    print(f"  Singer: {[f'{r:.4f}' for r in singer_ratios]}")
    print(f"  Greedy: {[f'{r:.4f}' for r in greedy_ratios]}")
    print(f"  Random: {[f'{r:.4f}' for r in random_ratios]}")

    # Singer trend: fit ratio = a · q^(-beta)
    if len(singer_ratios) >= 3:
        log_q = [math.log(q) for q in qs]
        log_r = [math.log(r) for r in singer_ratios if r > 0]
        log_q_valid = log_q[:len(log_r)]

        if len(log_q_valid) >= 2:
            n_pts = len(log_q_valid)
            sx = sum(log_q_valid)
            sy = sum(log_r)
            sxx = sum(x * x for x in log_q_valid)
            sxy = sum(x * y for x, y in zip(log_q_valid, log_r))
            denom = n_pts * sxx - sx * sx
            if abs(denom) > 1e-15:
                slope = (n_pts * sxy - sx * sy) / denom
                intercept = (sy - slope * sx) / n_pts
                beta = -slope
                C = math.exp(intercept)

                y_mean = sy / n_pts
                ss_tot = sum((y - y_mean) ** 2 for y in log_r)
                ss_res = sum((log_r[i] - (intercept + slope * log_q_valid[i])) ** 2 for i in range(n_pts))
                r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0

                print(f"\nSinger ratio decay: C(S)/√(k·ln n) ≈ {C:.4f} × q^(-{beta:.4f})")
                print(f"R² = {r_sq:.4f}")
                results["analysis"]["singer_decay_exponent"] = float(beta)
                results["analysis"]["singer_decay_r_squared"] = float(r_sq)

    # Separation test: is Singer ratio < Greedy ratio for large q?
    separation_count = sum(1 for s, g in zip(singer_ratios, greedy_ratios) if s < g)
    print(f"\nSinger < Greedy in {separation_count}/{len(singer_ratios)} cases")
    print(f"Singer < Random in {sum(1 for s, r in zip(singer_ratios, random_ratios) if s < r)}/{len(singer_ratios)} cases")

    # Ratio of ratios: Singer/Greedy
    if len(singer_ratios) == len(greedy_ratios):
        ratio_of_ratios = [s / g if g > 0 else float('inf')
                           for s, g in zip(singer_ratios, greedy_ratios)]
        print(f"\nSinger/Greedy ratio: {[f'{r:.3f}' for r in ratio_of_ratios]}")
        if ratio_of_ratios[-1] < ratio_of_ratios[0]:
            print("→ DIVERGING: Singer becoming relatively flatter vs greedy as q grows")
        else:
            print("→ NOT DIVERGING: Singer/Greedy ratio not decreasing")

    # ── Verdict ──
    print("\n" + "=" * 78)
    print("CHOWLA DEFICIENCY VERDICT")
    print("=" * 78)

    singer_decreasing = all(singer_ratios[i] >= singer_ratios[i + 1]
                           for i in range(len(singer_ratios) - 1))
    greedy_decreasing = all(greedy_ratios[i] >= greedy_ratios[i + 1]
                           for i in range(len(greedy_ratios) - 1))

    if singer_decreasing and not greedy_decreasing:
        verdict = "PASS"
        msg = ("Singer Chowla ratio monotonically decreasing; greedy is NOT. "
               "Singer sets are anomalously spectrally flat — consistent with "
               "quasicrystal/model-set prediction.")
    elif singer_decreasing:
        verdict = "PARTIAL"
        msg = "Singer ratio decreasing but so is greedy — may not be Singer-specific."
    else:
        verdict = "INCONCLUSIVE"
        msg = "Singer ratio not monotonically decreasing at this scale."

    print(f"\nVerdict: {verdict}")
    print(f"Reason: {msg}")

    results["analysis"]["verdict"] = verdict
    results["analysis"]["singer_ratios"] = singer_ratios
    results["analysis"]["greedy_ratios"] = greedy_ratios
    results["analysis"]["random_ratios"] = random_ratios
    results["analysis"]["singer_monotonically_decreasing"] = singer_decreasing
    results["analysis"]["greedy_monotonically_decreasing"] = greedy_decreasing

    # Save
    out_path = Path(__file__).parent / "EXP-MATH-ERDOS30-SIDON-003_RESULTS.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")

    return results


if __name__ == "__main__":
    run_experiment()
