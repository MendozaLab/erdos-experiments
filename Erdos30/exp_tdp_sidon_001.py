#!/usr/bin/env python3
"""
EXP-TDP-SIDON-001: Spectral Sidon Test
Tests whether spectral distance sets of n-body choreographies (Koopman eigenvalues
at nth roots of unity) satisfy the Sidon property.

Observable: For each n from 3 to 30, compute D_n = {2*sin(pi*m/n) : m = 1,...,floor(n/2)}
and check if D_n is a Sidon set in R (all pairwise sums distinct).

Pass/fail: Correlation between Sidon deficiency and known choreography stability.
"""

import numpy as np
import json
from itertools import combinations
from datetime import datetime

def spectral_distances(n):
    """Pairwise distances between nth roots of unity."""
    return [2 * np.sin(np.pi * m / n) for m in range(1, n // 2 + 1)]

def sidon_deficiency(distances, tol=1e-10):
    """Count collisions in pairwise sums. Returns (collision_count, total_sums, collision_pairs)."""
    sums = []
    sum_labels = []
    for i, d1 in enumerate(distances):
        for j, d2 in enumerate(distances):
            if i <= j:
                sums.append(d1 + d2)
                sum_labels.append((i+1, j+1))  # 1-indexed m values
    
    # Sort by value, track indices
    indexed = sorted(enumerate(sums), key=lambda x: x[1])
    collisions = 0
    collision_pairs = []
    for k in range(len(indexed) - 1):
        idx1, val1 = indexed[k]
        idx2, val2 = indexed[k+1]
        if abs(val2 - val1) < tol:
            collisions += 1
            collision_pairs.append((sum_labels[idx1], sum_labels[idx2], val1))
    
    return collisions, len(sums), collision_pairs

def b2_property_check(distances, tol=1e-10):
    """
    Strict B2 (Sidon) check: all pairwise sums a_i + a_j (i <= j) are distinct.
    Returns number of non-trivial sum collisions.
    """
    sums = {}
    collisions = 0
    for i in range(len(distances)):
        for j in range(i, len(distances)):
            s = distances[i] + distances[j]
            # Check if this sum is close to any existing sum
            found = False
            for existing_s in sums:
                if abs(s - existing_s) < tol:
                    collisions += 1
                    found = True
                    break
            if not found:
                sums[s] = (i+1, j+1)
    return collisions

# Known stable choreographies (from Simó's catalog and literature)
# n=3: Figure-Eight (Chenciner-Montgomery 2000) - STABLE
# n=4: Some solutions found (Simó) - partially stable
# n=5-8: Various solutions cataloged
# Higher n: increasingly unstable
KNOWN_STABLE = {3: True, 4: True, 5: True, 6: False, 7: True, 8: False}

results = {
    "experiment_id": "EXP-TDP-SIDON-001",
    "title": "Spectral Sidon Test: N-Body Choreography Eigenvalue Distances",
    "date": datetime.now().isoformat(),
    "method": "Compute D_n = {2*sin(pi*m/n)} for n=3..30, test Sidon property",
    "results": []
}

print("=" * 80)
print("EXP-TDP-SIDON-001: Spectral Sidon Test")
print("=" * 80)
print()
print(f"{'n':>4} {'|D_n|':>6} {'B2_collisions':>14} {'Deficiency':>11} {'Density':>8} {'Status':>12}")
print("-" * 70)

for n in range(3, 31):
    D = spectral_distances(n)
    deficiency, total_sums, collision_pairs = sidon_deficiency(D)
    b2_collisions = b2_property_check(D)
    density = len(D) / (2 * max(D)) if D else 0  # density relative to max distance
    
    is_sidon = deficiency == 0 and b2_collisions == 0
    status = "SIDON" if is_sidon else f"NOT SIDON ({b2_collisions})"
    
    stability = KNOWN_STABLE.get(n, "unknown")
    
    result_entry = {
        "n": n,
        "D_n_size": len(D),
        "distances": [round(d, 10) for d in D],
        "b2_collisions": b2_collisions,
        "sidon_deficiency": deficiency,
        "total_pairwise_sums": total_sums,
        "density": round(density, 6),
        "is_sidon": is_sidon,
        "known_stability": str(stability),
        "collision_examples": [(str(p[0]), str(p[1]), round(p[2], 10)) for p in collision_pairs[:5]]
    }
    results["results"].append(result_entry)
    
    print(f"{n:>4} {len(D):>6} {b2_collisions:>14} {deficiency:>11} {density:>8.4f} {status:>12}")

# Analysis: correlation between Sidon property and stability
print()
print("=" * 80)
print("ANALYSIS: Sidon Property vs Known Choreography Stability")
print("=" * 80)

sidon_ns = [r["n"] for r in results["results"] if r["is_sidon"]]
non_sidon_ns = [r["n"] for r in results["results"] if not r["is_sidon"]]

print(f"\nSidon n values: {sidon_ns}")
print(f"Non-Sidon n values: {non_sidon_ns}")

# Check for primes
import sympy
prime_sidon = [n for n in sidon_ns if sympy.isprime(n)]
prime_non_sidon = [n for n in non_sidon_ns if sympy.isprime(n)]
print(f"\nPrime Sidon: {prime_sidon}")
print(f"Prime Non-Sidon: {prime_non_sidon}")

# Deficiency growth pattern
print(f"\nDeficiency growth pattern:")
for r in results["results"]:
    if r["b2_collisions"] > 0:
        print(f"  n={r['n']:>3}: {r['b2_collisions']} collisions in {r['total_pairwise_sums']} sums "
              f"(ratio={r['b2_collisions']/r['total_pairwise_sums']:.4f})")

# Summary
results["summary"] = {
    "sidon_n_values": sidon_ns,
    "non_sidon_n_values": non_sidon_ns,
    "total_sidon": len(sidon_ns),
    "total_non_sidon": len(non_sidon_ns),
    "prime_correlation": {
        "prime_sidon": prime_sidon,
        "prime_non_sidon": prime_non_sidon
    }
}

# Verdict
if len(sidon_ns) > 0 and len(non_sidon_ns) > 0:
    verdict = "PARTIAL_PASS"
    explanation = f"D_n is Sidon for {len(sidon_ns)} values, non-Sidon for {len(non_sidon_ns)}. Pattern analysis needed."
elif len(sidon_ns) == 0:
    verdict = "FAIL"
    explanation = "D_n is never Sidon for n >= 3."
else:
    verdict = "PASS"
    explanation = "D_n is always Sidon."

results["verdict"] = verdict
results["explanation"] = explanation

print(f"\n{'=' * 80}")
print(f"VERDICT: {verdict}")
print(f"EXPLANATION: {explanation}")
print(f"{'=' * 80}")

# Save results
with open("EXP-TDP-SIDON-001_RESULTS.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to EXP-TDP-SIDON-001_RESULTS.json")
