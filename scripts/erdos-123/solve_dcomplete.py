#!/usr/bin/env python3
"""Erdos Problem #123: d-completeness of {a^k * b^l * c^m} for coprime a,b,c."""

import json
from itertools import product as iterproduct
from datetime import datetime

def generate_abc_set(a, b, c, max_val):
    """Generate {a^k * b^l * c^m : k,l,m >= 0} up to max_val."""
    S = set()
    ak = 1
    while ak <= max_val:
        bl = 1
        while ak * bl <= max_val:
            cm = 1
            while ak * bl * cm <= max_val:
                S.add(ak * bl * cm)
                cm *= c
            bl *= b
        ak *= a
    return sorted(S)

def check_d_complete(S, N):
    """Check if every integer 1..N can be represented as sum of distinct elements of S."""
    # DP: reachable[i] = True if i can be made as sum of distinct elements
    reachable = [False] * (N + 1)
    reachable[0] = True
    for s in S:
        if s > N:
            break
        # Traverse backwards to ensure distinct elements
        for i in range(N, s - 1, -1):
            if reachable[i - s]:
                reachable[i] = True
    # Find largest non-representable
    largest_gap = 0
    for i in range(N, 0, -1):
        if not reachable[i]:
            largest_gap = i
            break
    all_representable = all(reachable[1:N+1])
    gaps = [i for i in range(1, N+1) if not reachable[i]]
    return all_representable, largest_gap, len(gaps), gaps[:20]

# Test triples
triples = [(2,3,5), (2,3,7), (2,5,7), (3,5,7), (2,3,11), (2,7,11),
           (3,5,11), (3,7,11), (5,7,11), (2,3,13)]
N = 10000
results = []

for a, b, c in triples:
    S = generate_abc_set(a, b, c, N)
    complete, largest_gap, num_gaps, sample_gaps = check_d_complete(S, N)
    result = {
        "triple": [a, b, c],
        "N": N,
        "set_size": len(S),
        "d_complete_up_to_N": complete,
        "largest_gap": largest_gap,
        "num_gaps": num_gaps,
        "sample_gaps": sample_gaps
    }
    results.append(result)
    status = "D-COMPLETE" if complete else f"GAPS ({num_gaps}, largest={largest_gap})"
    print(f"  ({a},{b},{c}): |S|={len(S)}, {status}")

output = {
    "experiment_id": "EXP-ERDOS123-DCOMPLETE-01",
    "problem_id": 123,
    "date": datetime.utcnow().isoformat() + "Z",
    "N": N,
    "results": results,
    "verdict": "PASS" if all(r["d_complete_up_to_N"] for r in results) else "PARTIAL"
}

# Save
import os
os.makedirs("/Users/kenbengoetxea/container-projects/apps/H2/Research-Hub/erdos/results", exist_ok=True)
with open("/Users/kenbengoetxea/container-projects/apps/H2/Research-Hub/erdos/results/EXP-ERDOS123-DCOMPLETE-01_RESULTS.json", "w") as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to EXP-ERDOS123-DCOMPLETE-01_RESULTS.json")
