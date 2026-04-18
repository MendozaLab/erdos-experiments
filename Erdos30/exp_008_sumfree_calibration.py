#!/usr/bin/env python3
"""
EXP-008: Sum-Free Set Calibration — PMF Validation

Known answer: max sum-free subset of {1,...,N} has size ceil(N/2).
(Take all odd numbers, or all numbers > N/2.)

If the transfer matrix gives λ_max → 1 with the right scaling for
sum-free sets, the PMF machinery is validated. If not, the Sidon
transfer matrix formulation needs debugging.

Sum-free constraint: no a + b = c for a, b, c ∈ A.
(Weaker than Sidon: Sidon requires all pairwise SUMS distinct,
sum-free just requires no sum lands back in the set.)

Transfer matrix state: within window {n-W+1, ..., n}, which sites
are occupied? The constraint is that for occupied sites i < j,
if i + j ≤ n (within window), then site i+j must be unoccupied.
"""

import json
import math
import time
from collections import defaultdict

import numpy as np


def enumerate_sumfree_states(W):
    """Enumerate all sum-free subsets of {0, 1, ..., W-1}."""
    states = [frozenset()]
    state_index = {frozenset(): 0}

    def is_sumfree(positions):
        pos_set = set(positions)
        for i, a in enumerate(positions):
            for b in positions[i:]:
                if a + b in pos_set and a + b != a and a + b != b:
                    return False
                # Also check: a + b might equal a third element
                s = a + b
                if s in pos_set and s != a and s != b:
                    return False
        return True

    def backtrack(start, current):
        if current:
            fs = frozenset(current)
            idx = len(states)
            states.append(fs)
            state_index[fs] = idx

        for x in range(start, W):
            # Check if adding x violates sum-free property
            valid = True
            for a in current:
                # Check if a + x is in current
                if a + x < W and (a + x) in set(current):
                    valid = False
                    break
                # Check if x + a = some existing element
                # Check if there exist b in current with b + x = a (i.e. a - x = b)
                if a > x and (a - x) in set(current):
                    valid = False
                    break
                # Check if x - a is in current and x - a + a = x... no, that's trivial
                # The constraint is: no three elements a, b, c with a + b = c
                # Adding x: check if x = a + b for some a, b in current
                # or if a = x + b for some b in current (i.e. a - x in current)
                # or if b = x + a... already covered

            # More careful check
            if valid:
                cur_set = set(current)
                # Is x = a + b for some a, b in current?
                for a in current:
                    b = x - a
                    if b > 0 and b in cur_set and b != a:
                        valid = False
                        break
                    if b == a and current.count(a) >= 1:
                        # a + a = x requires two copies of a, but sets don't repeat
                        pass

            if valid:
                # Also check: x + a for each a in current — is it in current?
                for a in current:
                    s = x + a
                    if s in cur_set:
                        valid = False
                        break

            if valid:
                current.append(x)
                backtrack(x + 1, current)
                current.pop()

    backtrack(0, [])
    return states, state_index


def build_sumfree_transfer(W, states, state_index):
    """Build sparse transfer matrix for sum-free lattice gas."""
    n = len(states)
    adj = [[] for _ in range(n)]

    for i, state_i in enumerate(states):
        # Shift right: positions become {p-1 : p in state_i, p > 0}
        shifted = frozenset(p - 1 for p in state_i if p > 0)
        shifted_list = sorted(shifted)
        shifted_set = set(shifted)

        # Option 1: don't occupy W-1
        if shifted in state_index:
            j = state_index[shifted]
            adj[i].append(j)

        # Option 2: occupy W-1
        new_pos = W - 1
        valid = True

        # Check: new_pos + a must not be in shifted (but new_pos + a >= W, so outside window — OK)
        # Check: is new_pos = a + b for some a, b in shifted?
        for a in shifted_list:
            b = new_pos - a
            if b > 0 and b in shifted_set and b != a:
                valid = False
                break

        # Check: a + new_pos for a in shifted — these exceed W-1, so they're outside the window
        # But we need to check if any pair (a, new_pos) sums to another element in shifted
        if valid:
            for a in shifted_list:
                if a + new_pos in shifted_set:
                    # a + new_pos is in the window and occupied
                    valid = False
                    break

        if valid:
            new_state = shifted | {new_pos}
            if new_state in state_index:
                j = state_index[new_state]
                adj[i].append(j)

    return adj


def power_iteration(adj, n, max_iter=500, tol=1e-12):
    rng = np.random.RandomState(42)
    v = rng.random(n)
    v /= np.linalg.norm(v)

    lambda_old = 0
    for iteration in range(max_iter):
        w = np.zeros(n)
        for i in range(n):
            for j in adj[i]:
                w[j] += v[i]
        norm = np.linalg.norm(w)
        if norm < 1e-30:
            break
        v = w / norm
        if abs(norm - lambda_old) < tol * abs(norm):
            return norm, v, iteration + 1
        lambda_old = norm
    return lambda_old, v, max_iter


def main():
    print("=" * 80)
    print("EXP-008: Sum-Free Set Transfer Matrix — PMF Calibration")
    print("Known answer: max sum-free set in {1,...,N} has size ceil(N/2)")
    print("=" * 80)

    print(f"\n{'W':>4s} | {'#states':>8s} | {'#edges':>8s} | {'λ_max':>10s} | {'f.e./site':>10s}")
    print(f"{'-'*4}-+-{'-'*8}-+-{'-'*8}-+-{'-'*10}-+-{'-'*10}")

    results = []
    for W in range(3, 26):
        t0 = time.time()
        states, state_index = enumerate_sumfree_states(W)
        n_states = len(states)

        if n_states > 500000:
            print(f"  W={W}: {n_states} states — stopping")
            break

        adj = build_sumfree_transfer(W, states, state_index)
        n_edges = sum(len(a) for a in adj)

        lmax, _, iters = power_iteration(adj, n_states)
        fe = math.log(lmax) / W if lmax > 0 else 0

        elapsed = time.time() - t0
        print(f"{W:4d} | {n_states:8d} | {n_edges:8d} | {lmax:10.6f} | {fe:10.6f}")

        results.append({"W": W, "n_states": n_states, "lambda_max": lmax, "fe": fe})

    # Analysis
    print("\n" + "=" * 80)
    print("CALIBRATION CHECK")
    print("=" * 80)

    if len(results) >= 5:
        lambdas = [r["lambda_max"] for r in results]
        fes = [r["fe"] for r in results]

        print(f"\n  λ_max trend: {[f'{l:.4f}' for l in lambdas[-5:]]}")
        print(f"  f.e./site trend: {[f'{f:.6f}' for f in fes[-5:]]}")

        # For sum-free sets, the answer is N/2.
        # The density is 1/2, so the free energy should reflect this.
        # If the transfer matrix is correct, λ_max should → some value
        # that encodes the 1/2 density.

        # A sum-free set of density 1/2 in {1,...,N} exists (all odds).
        # The number of sum-free subsets of {1,...,N} of size N/2 should
        # grow exponentially → λ_max > 1.
        # The "critical density" is 1/2.

        last = results[-1]
        print(f"\n  At W={last['W']}: λ_max = {last['lambda_max']:.6f}")
        print(f"  For comparison, Sidon at W={last['W']} (from EXP-007C): check data")

        # The key test: does the sum-free transfer matrix have a DIFFERENT
        # convergence behavior than the Sidon one?
        if len(results) >= 3:
            Ws = np.array([r["W"] for r in results if r["W"] >= 5], dtype=float)
            lms = np.array([r["lambda_max"] for r in results if r["W"] >= 5], dtype=float)

            # Fit log(λ-1) vs log(W)
            excess = lms - 1.0
            if np.all(excess > 0):
                log_W = np.log(Ws)
                log_ex = np.log(excess)
                A = np.vstack([np.ones_like(log_W), log_W]).T
                coef, _, _, _ = np.linalg.lstsq(A, log_ex, rcond=None)
                alpha = -coef[1]
                c = math.exp(coef[0])
                print(f"\n  Fit: λ_max - 1 = {c:.4f} × W^(-{alpha:.4f})")
                print(f"  Sum-free α = {alpha:.4f}")
                print(f"  (Sidon α was 0.4416)")

                if alpha > 0.5:
                    print(f"  → Sum-free scaling is FASTER than Sidon — expected!")
                    print(f"    (sum-free answer is exact: density 1/2, not √N)")

    # Save
    out = {
        "experiment_id": "EXP-008",
        "title": "Sum-Free Calibration",
        "known_answer": "ceil(N/2)",
        "data": results,
    }
    with open("EXP-008_SUMFREE_CALIBRATION_RESULTS.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to EXP-008_SUMFREE_CALIBRATION_RESULTS.json")


if __name__ == "__main__":
    main()
