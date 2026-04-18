#!/usr/bin/env python3
"""
Spectral universality test: at Singer moduli N = q²+q+1,
compare the λ_max/k trajectory for:
  (a) Singer (PDS) prefixes
  (b) Greedy Sidon prefixes
  (c) Random Sidon sets (multiple trials)

Key question: does the post-peak descent of λ_max/k collapse
onto a universal curve, or is PDS behavior anomalous?
"""

import numpy as np
import math
import time
import json
from collections import Counter

# ── Singer set construction ──

SINGER = {
    2: ([0, 1, 3], 7),
    3: ([0, 1, 3, 9], 13),
    4: ([0, 1, 4, 14, 16], 21),
    5: ([0, 3, 5, 9, 23, 24], 31),
    7: ([0, 4, 26, 33, 36, 42, 44, 56], 57),
    8: ([0, 16, 38, 40, 43, 47, 53, 61, 72], 73),
    9: ([0, 10, 33, 42, 50, 70, 72, 76, 77, 88], 91),
    11: ([4, 6, 34, 40, 41, 44, 52, 61, 66, 85, 108, 124], 133),
    13: ([12, 14, 15, 23, 44, 61, 79, 112, 116, 122, 156, 170, 175, 182], 183),
}

def greedy_sidon(N, target_k=None):
    """Build a greedy Sidon set in Z_N."""
    A = [0]
    sums = {0}
    for x in range(1, N):
        new_sums = []
        conflict = False
        for a in A:
            s = x + a
            if s in sums:
                conflict = True
                break
            new_sums.append(s)
        if not conflict:
            s_xx = x + x
            if s_xx in sums:
                continue
            new_sums.append(s_xx)
            for s in new_sums:
                sums.add(s)
            A.append(x)
            if target_k and len(A) >= target_k:
                break
    return A

def random_sidon(N, target_k, max_attempts=50):
    """Build a random Sidon set of size target_k in Z_N via random insertion order."""
    import random
    for _ in range(max_attempts):
        candidates = list(range(N))
        random.shuffle(candidates)
        A = []
        sums = set()
        for x in candidates:
            new_sums = []
            conflict = False
            for a in A:
                s = x + a
                if s in sums:
                    conflict = True
                    break
                new_sums.append(s)
            if not conflict:
                s_xx = x + x
                if s_xx in sums:
                    continue
                new_sums.append(s_xx)
                for s in new_sums:
                    sums.add(s)
                A.append(x)
                if len(A) >= target_k:
                    break
        if len(A) >= target_k:
            return sorted(A[:target_k])
    return None

def compute_lambda_max(A, N):
    indicator = np.zeros(N)
    for a in A:
        indicator[a % N] = 1.0
    fft = np.fft.fft(indicator)
    spectrum = np.abs(fft) ** 2
    return float(np.max(spectrum[1:]))

def is_pds(A, N):
    """Check if A is a perfect difference set in Z_N."""
    k = len(A)
    if k * (k - 1) != N - 1:
        return False
    diffs = Counter()
    for i in range(k):
        for j in range(k):
            if i != j:
                diffs[(A[i] - A[j]) % N] += 1
    return all(diffs[d] == 1 for d in range(1, N))

# ── Main experiment ──

print("=" * 70)
print("SPECTRAL UNIVERSALITY TEST")
print("At Singer moduli: PDS vs non-PDS Sidon trajectories")
print("=" * 70)

results = {}

for q in [4, 5, 7, 8, 9, 11, 13]:
    A_singer, N = SINGER[q]
    k_singer = len(A_singer)
    sqrtN = math.sqrt(N)

    print(f"\n--- q={q}, N={N}, k_singer={k_singer}, √N={sqrtN:.2f}, k/√N={k_singer/sqrtN:.3f} ---")

    # Singer λ_max
    lmax_singer = compute_lambda_max(A_singer, N)
    print(f"  Singer: λ_max={lmax_singer:.4f}, k-1={k_singer-1}, ratio={lmax_singer/(k_singer-1):.6f}")

    # Greedy at same N, full set + prefixes
    A_greedy = greedy_sidon(N)
    k_greedy = len(A_greedy)
    print(f"  Greedy: k_max={k_greedy}, k/√N={k_greedy/sqrtN:.3f}")

    greedy_trajectory = []
    for k in range(3, k_greedy + 1):
        prefix = A_greedy[:k]
        lmax = compute_lambda_max(prefix, N)
        greedy_trajectory.append({
            'k': k,
            'k_over_sqrtN': k / sqrtN,
            'lambda_max': lmax,
            'lambda_max_over_k': lmax / k,
        })

    # Find greedy peak
    peak = max(greedy_trajectory, key=lambda r: r['lambda_max_over_k'])
    print(f"  Greedy peak: k={peak['k']} (k/√N={peak['k_over_sqrtN']:.3f}), λ_max/k={peak['lambda_max_over_k']:.3f}")

    final = greedy_trajectory[-1]
    print(f"  Greedy final: k={final['k']}, λ_max/k={final['lambda_max_over_k']:.3f}")

    # Random Sidon sets at k = k_singer (matched size)
    n_random = 20
    random_lmax = []
    for trial in range(n_random):
        A_rand = random_sidon(N, k_singer)
        if A_rand is not None:
            lmax = compute_lambda_max(A_rand, N)
            random_lmax.append(lmax)

    if random_lmax:
        mean_rand = np.mean(random_lmax)
        std_rand = np.std(random_lmax)
        min_rand = min(random_lmax)
        max_rand = max(random_lmax)
        print(f"  Random (k={k_singer}, n={len(random_lmax)} trials):")
        print(f"    λ_max: mean={mean_rand:.3f}, std={std_rand:.3f}, min={min_rand:.3f}, max={max_rand:.3f}")
        print(f"    λ_max/k: mean={mean_rand/k_singer:.3f}, Singer={lmax_singer/k_singer:.3f}")
        print(f"    Gap: random_mean/Singer = {mean_rand/lmax_singer:.3f}")

    # Also: random Sidon at various k values to trace trajectory
    random_trajectories = []
    test_ks = list(range(3, min(k_singer + 3, k_greedy + 1)))
    for k in test_ks:
        trial_lmax = []
        for _ in range(10):
            A_rand = random_sidon(N, k)
            if A_rand is not None:
                trial_lmax.append(compute_lambda_max(A_rand, N))
        if trial_lmax:
            random_trajectories.append({
                'k': k,
                'k_over_sqrtN': k / sqrtN,
                'lambda_max_mean': float(np.mean(trial_lmax)),
                'lambda_max_std': float(np.std(trial_lmax)),
                'lambda_max_over_k_mean': float(np.mean(trial_lmax)) / k,
                'n_trials': len(trial_lmax),
            })

    results[q] = {
        'N': N,
        'k_singer': k_singer,
        'sqrtN': sqrtN,
        'singer_lmax': lmax_singer,
        'greedy_trajectory': greedy_trajectory,
        'random_at_k_singer': random_lmax,
        'random_trajectories': random_trajectories,
    }

# ── Cross-q summary ──
print("\n" + "=" * 70)
print("CROSS-q SUMMARY: Singer vs Random at matched k")
print("=" * 70)
print(f"{'q':>3s} {'N':>5s} {'k':>3s} {'k/√N':>6s} {'Singer λ/k':>11s} {'Random λ/k':>11s} {'Gap':>6s}")
print("-" * 50)
for q in sorted(results.keys()):
    r = results[q]
    k = r['k_singer']
    s_ratio = r['singer_lmax'] / k
    if r['random_at_k_singer']:
        r_ratio = np.mean(r['random_at_k_singer']) / k
        gap = r_ratio / s_ratio
    else:
        r_ratio = float('nan')
        gap = float('nan')
    print(f"{q:3d} {r['N']:5d} {k:3d} {k/r['sqrtN']:6.3f} {s_ratio:11.4f} {r_ratio:11.4f} {gap:6.2f}x")

# ── Key test: at Singer moduli, do non-Singer sets EVER achieve λ_max/k close to 1? ──
print("\n" + "=" * 70)
print("KEY TEST: Closest any non-Singer set gets to λ_max = k-1")
print("=" * 70)
for q in sorted(results.keys()):
    r = results[q]
    k = r['k_singer']
    target = k - 1  # Singer floor
    if r['random_at_k_singer']:
        closest = min(r['random_at_k_singer'])
        print(f"  q={q}: Singer floor={target}, closest random={closest:.3f}, "
              f"ratio={closest/target:.3f}, gap={closest-target:.3f}")

# Save
serializable = {}
for q, r in results.items():
    serializable[str(q)] = {
        'N': r['N'],
        'k_singer': r['k_singer'],
        'sqrtN': r['sqrtN'],
        'singer_lmax': r['singer_lmax'],
        'greedy_trajectory': r['greedy_trajectory'],
        'random_at_k_singer': r['random_at_k_singer'],
        'random_trajectories': r['random_trajectories'],
    }
with open('/tmp/spectral_universality_results.json', 'w') as f:
    json.dump(serializable, f, indent=2)
print(f"\nResults saved to /tmp/spectral_universality_results.json")
