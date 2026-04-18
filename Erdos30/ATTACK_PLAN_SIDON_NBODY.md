# Attack Plan: Erdős #30 (Sidon Sets) × N-Body Choreographies

**Author:** Kenneth A. Mendoza / MendozaLab
**Date:** 2026-04-15
**Status:** ATTACK PLAN — pre-experiment
**Existing asset:** Erdős #30 Lean proof COMPILED (0 sorries, lake build PASS)

---

## The Bridge

An n-body choreography with Z_n symmetry has Koopman eigenvalues at the nth roots of unity: λ_k = exp(2πik/n) for k = 0, ..., n-1.

The **spectral distance set** is:

    D_n = { |λ_j − λ_k| : 0 ≤ j < k ≤ n−1 }

which equals:

    D_n = { |exp(2πij/n) − exp(2πik/n)| } = { 2·sin(π(j−k)/n) : 1 ≤ j−k ≤ ⌊n/2⌋ }

**Question 1 (Sidon condition):** For which n is D_n a Sidon set — i.e., are all pairwise sums of spectral distances distinct?

**Question 2 (Stability):** Does the Sidon property of D_n correlate with linear stability of the n-body choreography?

**Question 3 (Capacity):** Does the Sidon bound |A| ≤ √N + O(N^{1/4}) constrain the maximum number of *stable* choreographies at a given n?

---

## Why This Is Novel

The literature has two separate streams with no cross-citation:

- **Stream A (Additive Combinatorics):** Sidon sets on Z_n, rainbow Sidon numbers, spectral gap properties on character groups (Bożejko-Pytlik, Blacker)
- **Stream B (Dynamical Systems):** Koopman operators for choreographies, Z_n equivariant spectral theory (Montaldi-Steckles 2013, Koopman symmetry papers)

Nobody has asked: "Is the spectral distance set of an n-body choreography a Sidon set?"

---

## Computational Attack (Phase 1 — Falsifiable)

### Experiment EXP-TDP-SIDON-001: Spectral Sidon Test

**Observable signature:** For each n from 3 to 20, compute D_n = {2·sin(π·m/n) : m = 1, ..., ⌊n/2⌋} and check:
1. Is D_n itself a Sidon set in ℝ? (All pairwise sums distinct?)
2. If not exactly Sidon, what is the "Sidon deficiency" — how many collisions exist?
3. Does deficiency correlate with known choreography stability (Simó's stability data)?

**Pass/fail:** 
- PASS if Sidon deficiency = 0 for n where stable choreographies exist, deficiency > 0 where they don't
- PARTIAL PASS if correlation exists but isn't perfect
- FAIL if no correlation between Sidon property and stability

**Implementation:** Pure Python, ~50 lines. The distance set is deterministic — no simulation needed.

```python
import numpy as np
from itertools import combinations

def spectral_distances(n):
    """Pairwise distances between nth roots of unity."""
    return [2 * np.sin(np.pi * m / n) for m in range(1, n // 2 + 1)]

def sidon_deficiency(distances, tol=1e-10):
    """Count collisions in pairwise sums."""
    sums = []
    for i, d1 in enumerate(distances):
        for j, d2 in enumerate(distances):
            if i <= j:
                sums.append(d1 + d2)
    sums.sort()
    collisions = 0
    for i in range(len(sums) - 1):
        if abs(sums[i+1] - sums[i]) < tol:
            collisions += 1
    return collisions

# Test n = 3 to 20
for n in range(3, 21):
    D = spectral_distances(n)
    deficiency = sidon_deficiency(D)
    print(f"n={n:2d}  |D|={len(D):2d}  Sidon deficiency={deficiency}  "
          f"{'SIDON' if deficiency == 0 else 'NOT SIDON'}")
```

### Experiment EXP-TDP-SIDON-002: Sidon Bound vs Choreography Count

**Observable:** For each n, compare:
- Sidon capacity bound: |A| ≤ √(2N) where N = max element of the "discretized" distance set
- Known count of distinct choreographic solutions (from Simó's catalog)

**Pass/fail:** If Sidon bound ≈ choreography count for small n (within factor of 2), that's a publishable correlation.

### Experiment EXP-TDP-SIDON-003: MDL Encoding of Spectral Distances

**Observable:** Encode D_n using the MDL framework from TDP:
- L_model = cost of specifying which distances are "active" (choreography-compatible)
- L_data = residual after fitting the Sidon structure

**Pass/fail:** If MDL-optimal subset of D_n matches the observed choreography spectrum, MDL ↔ Sidon ↔ stability triangle closes.

---

## Formal Attack (Phase 2 — Lean 4)

### Theorem Target: `spectral_sidon_of_choreography`

**Statement (informal):** If φ is a measure-preserving Z_n-choreography and K is its Koopman operator, then the spectral distance set D(K) satisfies the Sidon sum-count bound |D(K)|·(|D(K)|+1)/2 ≤ 2n.

**Proof strategy:**
1. **Import existing assets:**
   - `sidon_sum_count` from Sidon_SumCount_Fix.lean (COMPILED, 0 sorries)
   - `koopman_isometry` from Layer2_KoopmanOperator.lean (TDP, 0 sorries in that theorem)
   - `eigenvalue_on_circle` from Layer3a_SpectralUnitary.lean (TDP, 1 sorry)

2. **New lemma needed:** `spectral_distances_of_roots_of_unity` — the pairwise distances between nth roots of unity form a set of size ⌊n/2⌋ with explicit formula 2·sin(πm/n).

3. **Bridge lemma:** `sidon_sum_count` applied to the spectral distance set gives the choreography constraint.

**This would be the first formal proof connecting additive combinatorics to celestial mechanics.**

---

## Connection to Existing TDP Architecture

```
TDP Layer 2 (Koopman)          Erdős #30 (Sidon)
     │                              │
     │  eigenvalues = roots of      │  pairwise sums
     │  unity for Z_n               │  must be distinct
     │                              │
     └──────────┬───────────────────┘
                │
    Spectral Distance Set D_n
                │
         ┌──────┴──────┐
         │             │
    Sidon bound    Stability
    constrains     of n-body
    |D_n|          choreography
```

The TDP paper (Layer 3b) already proves Z₃ → periodicity via roots of unity. This attack extends it: roots of unity → Sidon structure → constraints on which n admit stable choreographies.

---

## What This Buys You

1. **A second Erdős-to-celestial bridge** (alongside #114 lemniscate)
2. **Cross-citation between two compiled Lean proofs** (#30 Sidon + TDP Koopman)
3. **A novel result in additive combinatorics ∩ dynamical systems** that nobody has published
4. **Strengthens the Erdős Atlas patent** — morphism #30↔choreography is a concrete example of the discovery engine
5. **Natural Mathstodon post** — "I found a Sidon set hiding inside the 3-body problem"
6. **Falsifiable in 50 lines of Python** before any formal work

---

## Priority

Run EXP-TDP-SIDON-001 first. It takes 30 seconds and tells you immediately whether the bridge holds or collapses. If the Sidon deficiency pattern correlates with anything meaningful, proceed to formal proof. If it's random, kill this line and move on.

Anti-Slop Rule compliance: Observable signature defined. Pass/fail criterion written. Run it.
