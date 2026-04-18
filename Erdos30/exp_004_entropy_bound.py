#!/usr/bin/env python3
"""
EXP-MATH-ERDOS30-SIDON-004: Entropy Bound on Sidon Deviations

Track C experiment for Erdős #30 orthogonal attack.

Question: Do Singer sets have LOWER deviation-process entropy than
greedy/random Sidon sets? If so, this supports the discrepancy-entropy
connection (morphism to Erdős Discrepancy Problem, str=1.000).

The deviation process: X_A(t) = |A ∩ [1,t]| - (k/n)·t
This measures how far A's counting function deviates from uniform density.

For Singer PDS: the deviation should be rigid (low entropy, low max deviation).
For random Sidon: higher entropy (more structural freedom).

Connection to Tao's entropy decrement: if near-extremal Sidon sets with
large deviations have high entropy, and algebraic Sidon sets are forced to
have low entropy, then large deviations are incompatible with near-extremality.
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
    """Random Sidon set construction in Z_n."""
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


# ─── Deviation Process Analysis ──────────────────────────────────────────

def deviation_analysis(A: list[int], n: int) -> dict:
    """Full deviation process analysis for set A in Z_n."""
    A_set = set(A)
    k = len(A)
    density = k / n

    # Deviation process: X(t) = |A ∩ [0,t]| - density · (t+1)
    deviations = []
    count = 0
    for t in range(n):
        if t in A_set:
            count += 1
        dev = count - density * (t + 1)
        deviations.append(dev)

    dev_arr = np.array(deviations)

    # Basic statistics
    max_dev = float(np.max(np.abs(dev_arr)))
    mean_dev = float(np.mean(dev_arr))
    std_dev = float(np.std(dev_arr))
    rms_dev = float(np.sqrt(np.mean(dev_arr**2)))

    # Shannon entropy of discretized deviation
    discretized = [round(d, 1) for d in deviations]
    counts = Counter(discretized)
    total = len(discretized)
    entropy = 0.0
    for c in counts.values():
        p = c / total
        if p > 0:
            entropy -= p * math.log2(p)

    # Normalized entropy: H(X) / log₂(n)
    max_entropy = math.log2(n)
    normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0

    # Discrepancy: max |X(t)| / √k (Lindström-type normalization)
    lindstrom_discrepancy = max_dev / math.sqrt(k) if k > 0 else 0

    # Autocorrelation of deviation process (structural rigidity indicator)
    if len(dev_arr) > 1:
        dev_centered = dev_arr - np.mean(dev_arr)
        autocorr_1 = float(np.correlate(dev_centered[:-1], dev_centered[1:])[0] /
                          (np.sum(dev_centered**2) + 1e-15))
    else:
        autocorr_1 = 0.0

    # Gap statistics: spacings between consecutive elements
    A_sorted = sorted(A)
    if len(A_sorted) > 1:
        gaps = [A_sorted[i + 1] - A_sorted[i] for i in range(len(A_sorted) - 1)]
        gap_arr = np.array(gaps, dtype=float)
        gap_mean = float(np.mean(gap_arr))
        gap_std = float(np.std(gap_arr))
        gap_cv = gap_std / gap_mean if gap_mean > 0 else 0

        # Gap entropy
        gap_counts = Counter(gaps)
        gap_total = len(gaps)
        gap_entropy = 0.0
        for c in gap_counts.values():
            p = c / gap_total
            if p > 0:
                gap_entropy -= p * math.log2(p)
    else:
        gap_cv = 0.0
        gap_entropy = 0.0

    # Kolmogorov-Smirnov-style statistic: max deviation from uniform CDF
    if k > 0:
        uniform_cdf = np.linspace(1 / k, 1.0, k)
        empirical_cdf = np.array(sorted(A)) / n
        ks_stat = float(np.max(np.abs(empirical_cdf - uniform_cdf)))
    else:
        ks_stat = 0.0

    return {
        "max_deviation": max_dev,
        "rms_deviation": rms_dev,
        "std_deviation": std_dev,
        "shannon_entropy": entropy,
        "normalized_entropy": normalized_entropy,
        "lindstrom_discrepancy": lindstrom_discrepancy,
        "autocorrelation_lag1": autocorr_1,
        "gap_cv": gap_cv,
        "gap_entropy": gap_entropy,
        "ks_statistic": ks_stat,
    }


# ─── Main Experiment ─────────────────────────────────────────────────────

def run_experiment():
    results = {
        "experiment_id": "EXP-MATH-ERDOS30-SIDON-004",
        "title": "Entropy Bound on Sidon Deviations — Discrepancy-Entropy Connection",
        "morphism_link": "Erdős Discrepancy Problem (str=1.000, solved by Tao 2015)",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "singer": [],
        "greedy": [],
        "random": [],
        "analysis": {},
    }

    print("=" * 85)
    print("EXP-MATH-ERDOS30-SIDON-004: Entropy Bound on Sidon Deviations")
    print("=" * 85)

    # ── Singer ──
    print("\n--- Singer Perfect Difference Sets ---")
    header = (f"{'q':>4s} | {'n':>7s} | {'k':>4s} | {'H(X)':>7s} | {'H_norm':>7s} | "
              f"{'maxDev':>8s} | {'L_disc':>7s} | {'AC(1)':>7s} | {'GapCV':>7s} | {'KS':>7s}")
    print(header)
    print("-" * len(header))

    singer_entropies = []
    singer_norm_ent = []
    singer_discrep = []

    for q in sorted(SINGER_SETS.keys()):
        n = q * q + q + 1
        A = SINGER_SETS[q]
        k = len(A)

        analysis = deviation_analysis(A, n)
        singer_entropies.append(analysis["shannon_entropy"])
        singer_norm_ent.append(analysis["normalized_entropy"])
        singer_discrep.append(analysis["lindstrom_discrepancy"])

        print(f"{q:4d} | {n:7d} | {k:4d} | {analysis['shannon_entropy']:7.3f} | "
              f"{analysis['normalized_entropy']:7.4f} | "
              f"{analysis['max_deviation']:8.3f} | "
              f"{analysis['lindstrom_discrepancy']:7.4f} | "
              f"{analysis['autocorrelation_lag1']:7.4f} | "
              f"{analysis['gap_cv']:7.4f} | "
              f"{analysis['ks_statistic']:7.4f}")

        results["singer"].append({"q": q, "n": n, "k": k, **analysis})

    # ── Greedy ──
    print("\n--- Greedy Sidon Sets ---")
    print(header)
    print("-" * len(header))

    greedy_entropies = []
    greedy_norm_ent = []
    greedy_discrep = []

    for q in sorted(SINGER_SETS.keys()):
        n = q * q + q + 1
        k = q + 1
        A = greedy_sidon(n, k)
        actual_k = len(A)

        analysis = deviation_analysis(A, n)
        greedy_entropies.append(analysis["shannon_entropy"])
        greedy_norm_ent.append(analysis["normalized_entropy"])
        greedy_discrep.append(analysis["lindstrom_discrepancy"])

        print(f"{q:4d} | {n:7d} | {actual_k:4d} | {analysis['shannon_entropy']:7.3f} | "
              f"{analysis['normalized_entropy']:7.4f} | "
              f"{analysis['max_deviation']:8.3f} | "
              f"{analysis['lindstrom_discrepancy']:7.4f} | "
              f"{analysis['autocorrelation_lag1']:7.4f} | "
              f"{analysis['gap_cv']:7.4f} | "
              f"{analysis['ks_statistic']:7.4f}")

        results["greedy"].append({"q": q, "n": n, "k": actual_k, **analysis})

    # ── Random (5 seeds) ──
    print("\n--- Random Sidon Sets (averaged over 5 seeds) ---")
    short_header = f"{'q':>4s} | {'n':>7s} | {'k':>4s} | {'H(X)':>7s} | {'H_norm':>7s} | {'maxDev':>8s} | {'L_disc':>7s}"
    print(short_header)
    print("-" * len(short_header))

    random_entropies = []
    random_norm_ent = []
    random_discrep = []

    for q in sorted(SINGER_SETS.keys()):
        n = q * q + q + 1
        k = q + 1
        seed_analyses = []
        for seed in range(5):
            A = random_sidon(n, k, seed=seed * 7 + 13)
            seed_analyses.append(deviation_analysis(A, n))

        avg_ent = float(np.mean([a["shannon_entropy"] for a in seed_analyses]))
        avg_norm = float(np.mean([a["normalized_entropy"] for a in seed_analyses]))
        avg_maxdev = float(np.mean([a["max_deviation"] for a in seed_analyses]))
        avg_discrep = float(np.mean([a["lindstrom_discrepancy"] for a in seed_analyses]))
        random_entropies.append(avg_ent)
        random_norm_ent.append(avg_norm)
        random_discrep.append(avg_discrep)

        print(f"{q:4d} | {n:7d} | {k:4d} | {avg_ent:7.3f} | {avg_norm:7.4f} | "
              f"{avg_maxdev:8.3f} | {avg_discrep:7.4f}")

        results["random"].append({
            "q": q, "n": n, "k": k,
            "shannon_entropy_avg": avg_ent,
            "normalized_entropy_avg": avg_norm,
            "max_deviation_avg": avg_maxdev,
            "lindstrom_discrepancy_avg": avg_discrep,
        })

    # ── Trend Analysis ──
    print("\n" + "=" * 85)
    print("TREND ANALYSIS")
    print("=" * 85)

    qs = sorted(SINGER_SETS.keys())

    print(f"\nNormalized entropy H(X)/log₂(n):")
    print(f"  Singer: {[f'{e:.4f}' for e in singer_norm_ent]}")
    print(f"  Greedy: {[f'{e:.4f}' for e in greedy_norm_ent]}")
    print(f"  Random: {[f'{e:.4f}' for e in random_norm_ent]}")

    # Entropy ratio: Singer/Greedy
    ent_ratios = [s / g if g > 0 else float('inf')
                  for s, g in zip(singer_norm_ent, greedy_norm_ent)]
    print(f"\nSinger/Greedy entropy ratio: {[f'{r:.3f}' for r in ent_ratios]}")

    # Discrepancy comparison
    print(f"\nLindström discrepancy (max|X|/√k):")
    print(f"  Singer: {[f'{d:.4f}' for d in singer_discrep]}")
    print(f"  Greedy: {[f'{d:.4f}' for d in greedy_discrep]}")
    print(f"  Random: {[f'{d:.4f}' for d in random_discrep]}")

    disc_ratios = [s / g if g > 0 else float('inf')
                   for s, g in zip(singer_discrep, greedy_discrep)]
    print(f"\nSinger/Greedy discrepancy ratio: {[f'{r:.3f}' for r in disc_ratios]}")

    # ── Verdict ──
    print("\n" + "=" * 85)
    print("ENTROPY-DISCREPANCY VERDICT")
    print("=" * 85)

    # Test 1: Singer entropy < Greedy entropy (consistently)?
    ent_lower = sum(1 for s, g in zip(singer_norm_ent, greedy_norm_ent) if s < g)
    # Test 2: Singer discrepancy < Greedy discrepancy?
    disc_lower = sum(1 for s, g in zip(singer_discrep, greedy_discrep) if s < g)
    # Test 3: Entropy ratio trending down?
    ent_ratio_decreasing = all(ent_ratios[i] >= ent_ratios[i + 1]
                               for i in range(len(ent_ratios) - 1)) if len(ent_ratios) > 1 else False

    print(f"\n  Singer entropy < Greedy: {ent_lower}/{len(singer_norm_ent)} cases")
    print(f"  Singer discrepancy < Greedy: {disc_lower}/{len(singer_discrep)} cases")
    print(f"  Entropy ratio monotonically decreasing: {ent_ratio_decreasing}")

    if ent_lower > len(singer_norm_ent) // 2 and disc_lower > len(singer_discrep) // 2:
        verdict = "PASS"
        msg = ("Singer sets have consistently lower entropy AND lower discrepancy "
               "than non-algebraic Sidon sets. This supports the hypothesis that "
               "near-extremal Sidon sets are structurally rigid (low entropy) and "
               "well-distributed (low discrepancy) — the entropy-discrepancy connection "
               "predicted by the Erdős Discrepancy morphism.")
    elif ent_lower > len(singer_norm_ent) // 2 or disc_lower > len(singer_discrep) // 2:
        verdict = "PARTIAL"
        msg = "One of entropy or discrepancy is lower for Singer, but not both."
    else:
        verdict = "FAIL"
        msg = "Singer sets do not show entropy/discrepancy advantage over greedy/random."

    print(f"\n  Verdict: {verdict}")
    print(f"  Reason: {msg}")

    results["analysis"]["verdict"] = verdict
    results["analysis"]["singer_norm_entropies"] = singer_norm_ent
    results["analysis"]["greedy_norm_entropies"] = greedy_norm_ent
    results["analysis"]["random_norm_entropies"] = random_norm_ent
    results["analysis"]["entropy_ratio_singer_greedy"] = ent_ratios
    results["analysis"]["discrepancy_ratio_singer_greedy"] = disc_ratios

    # Save
    out_path = Path(__file__).parent / "EXP-MATH-ERDOS30-SIDON-004_RESULTS.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")

    return results


if __name__ == "__main__":
    run_experiment()
