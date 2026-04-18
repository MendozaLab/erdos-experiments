#!/usr/bin/env python3
"""
Attack A test: λ_max = k-1 iff Perfect Difference Set (Singer).

At Singer moduli N = q²+q+1, compare:
  - Singer set: λ_max = k-1 (expected, verified)
  - Greedy set of SAME size k: λ_max > k-1 (expected)
  - All Sidon sets of same size (exhaustive for small N)

If λ_max = k-1 is achievable ONLY by PDS, this characterizes
spectral flatness as equivalent to the PDS condition.

Also test: at non-Singer N, does ANY Sidon set achieve λ_max = k-1?
"""

import numpy as np
import math
from itertools import combinations

def is_sidon(A):
    sums = set()
    for i in range(len(A)):
        for j in range(i, len(A)):
            s = A[i] + A[j]
            if s in sums:
                return False
            sums.add(s)
    return True

def compute_lambda_max(A, N):
    indicator = np.zeros(N)
    for a in A:
        indicator[a % N] = 1.0
    fft = np.fft.fft(indicator)
    spectrum = np.abs(fft) ** 2
    return float(np.max(spectrum[1:]))

def compute_full_spectrum(A, N):
    indicator = np.zeros(N)
    for a in A:
        indicator[a % N] = 1.0
    fft = np.fft.fft(indicator)
    return np.abs(fft) ** 2

def greedy_sidon(N, target_k=None):
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

# Singer sets
SINGER = {
    2: ([0, 1, 3], 7),
    3: ([0, 1, 3, 9], 13),
    4: ([0, 1, 4, 14, 16], 21),
    5: ([0, 3, 5, 9, 23, 24], 31),
}

def all_sidon_sets(N, k):
    """Enumerate all Sidon sets of size k in {0,...,N-1}. Only for small N."""
    count = 0
    results = []
    for combo in combinations(range(N), k):
        A = list(combo)
        if is_sidon(A):
            lmax = compute_lambda_max(A, N)
            results.append((A, lmax))
            count += 1
    return results

print("=" * 70)
print("ATTACK A: Uniqueness of Spectral Flatness")
print("=" * 70)

# ─── Test 1: Singer moduli, exhaustive for small N ───
print("\n─── Test 1: Exhaustive enumeration at Singer moduli ───")

for q in [2, 3]:
    A_singer, N = SINGER[q]
    k = q + 1

    print(f"\nq={q}: N={N}, k={k}")
    print(f"  Singer set: {A_singer}")
    lmax_singer = compute_lambda_max(A_singer, N)
    print(f"  Singer λ_max = {lmax_singer:.6f}, k-1 = {k-1}")

    # Enumerate ALL Sidon sets of size k
    all_sets = all_sidon_sets(N, k)
    print(f"  Total Sidon sets of size {k}: {len(all_sets)}")

    # Check which achieve λ_max = k-1
    flat_sets = [(A, lm) for A, lm in all_sets if abs(lm - (k-1)) < 0.01]
    nonflat_sets = [(A, lm) for A, lm in all_sets if abs(lm - (k-1)) >= 0.01]

    print(f"  Sets with λ_max = k-1 (flat): {len(flat_sets)}")
    for A, lm in flat_sets:
        print(f"    {A} → λ_max={lm:.4f}")

    print(f"  Sets with λ_max > k-1 (non-flat): {len(nonflat_sets)}")
    if nonflat_sets:
        # Sort by λ_max
        nonflat_sets.sort(key=lambda x: x[1])
        print(f"    min λ_max: {nonflat_sets[0][1]:.4f} (set: {nonflat_sets[0][0]})")
        print(f"    max λ_max: {nonflat_sets[-1][1]:.4f} (set: {nonflat_sets[-1][0]})")
        print(f"    mean λ_max: {np.mean([lm for _, lm in nonflat_sets]):.4f}")

    # Check if flat sets are PDS
    print(f"\n  Are all flat sets perfect difference sets?")
    for A, lm in flat_sets:
        from collections import Counter
        diffs = Counter()
        for i in range(len(A)):
            for j in range(len(A)):
                if i != j:
                    diffs[(A[i] - A[j]) % N] += 1
        is_pds = all(diffs[d] == 1 for d in range(1, N))
        print(f"    {A}: PDS = {is_pds}")

# ─── Test 2: Non-Singer N, check if any set achieves λ_max = k-1 ───
print("\n─── Test 2: Non-Singer N (no PDS exists) ───")

for N in [10, 11, 14, 15, 17, 19, 20, 23, 25]:
    # Find max k for Sidon sets in this N
    A_greedy = greedy_sidon(N)
    k = len(A_greedy)

    all_sets = all_sidon_sets(N, k)
    min_lmax = min(lm for _, lm in all_sets) if all_sets else float('inf')
    threshold = k - 1

    flat_count = sum(1 for _, lm in all_sets if abs(lm - threshold) < 0.01)

    is_singer_N = any(N == q*q+q+1 for q in range(2, 20))
    marker = "★ Singer N" if is_singer_N else ""

    print(f"  N={N:3d} k={k} total_sidon={len(all_sets):5d} "
          f"min_λ_max={min_lmax:.3f} k-1={threshold} "
          f"flat_sets={flat_count} {marker}")

# ─── Test 3: q=4 (N=21) with larger search ───
print("\n─── Test 3: q=4 (N=21, k=5) — larger exhaustive search ───")
q = 4
A_singer, N = SINGER[q]
k = q + 1

lmax_singer = compute_lambda_max(A_singer, N)
print(f"  Singer: {A_singer}, λ_max = {lmax_singer:.6f}, k-1 = {k-1}")

all_sets = all_sidon_sets(N, k)
print(f"  Total Sidon sets of size {k} in {{0,...,{N-1}}}: {len(all_sets)}")

flat_sets = [(A, lm) for A, lm in all_sets if abs(lm - (k-1)) < 0.01]
print(f"  Sets with λ_max = k-1: {len(flat_sets)}")
for A, lm in flat_sets:
    from collections import Counter
    diffs = Counter()
    for i in range(len(A)):
        for j in range(len(A)):
            if i != j:
                diffs[(A[i] - A[j]) % N] += 1
    is_pds = all(diffs[d] == 1 for d in range(1, N))
    print(f"    {A}: PDS={is_pds}, λ_max={lm:.6f}")

if all_sets:
    lmax_vals = [lm for _, lm in all_sets]
    print(f"\n  λ_max distribution:")
    print(f"    min={min(lmax_vals):.4f}, max={max(lmax_vals):.4f}")
    print(f"    mean={np.mean(lmax_vals):.4f}, std={np.std(lmax_vals):.4f}")
    # Histogram
    bins = np.linspace(min(lmax_vals) - 0.01, max(lmax_vals) + 0.01, 20)
    hist, edges = np.histogram(lmax_vals, bins=bins)
    print(f"    Histogram:")
    for i in range(len(hist)):
        if hist[i] > 0:
            bar = '#' * min(hist[i], 60)
            print(f"      [{edges[i]:.2f}-{edges[i+1]:.2f}): {hist[i]:4d} {bar}")

# ─── Test 4: q=5 (N=31, k=6) ───
print("\n─── Test 4: q=5 (N=31, k=6) ───")
q = 5
A_singer, N = SINGER[q]
k = q + 1

lmax_singer = compute_lambda_max(A_singer, N)
print(f"  Singer: {A_singer}, λ_max = {lmax_singer:.6f}, k-1 = {k-1}")

print(f"  Enumerating all Sidon sets of size {k} in {{0,...,{N-1}}}...")
all_sets = all_sidon_sets(N, k)
print(f"  Total: {len(all_sets)}")

flat_sets = [(A, lm) for A, lm in all_sets if abs(lm - (k-1)) < 0.01]
print(f"  Sets with λ_max = k-1: {len(flat_sets)}")
for A, lm in flat_sets[:20]:  # show first 20
    from collections import Counter
    diffs = Counter()
    for i in range(len(A)):
        for j in range(len(A)):
            if i != j:
                diffs[(A[i] - A[j]) % N] += 1
    is_pds = all(diffs[d] == 1 for d in range(1, N))
    print(f"    {A}: PDS={is_pds}")

if flat_sets:
    print(f"  All flat sets are PDS: {all(all(Counter((a-b)%N for a in A for b in A if a!=b)[d]==1 for d in range(1,N)) for A, _ in flat_sets)}")

if all_sets:
    lmax_vals = sorted(set(round(lm, 4) for _, lm in all_sets))
    print(f"\n  Distinct λ_max values: {len(lmax_vals)}")
    print(f"  Smallest 5: {lmax_vals[:5]}")
    print(f"  Largest 5: {lmax_vals[-5:]}")
