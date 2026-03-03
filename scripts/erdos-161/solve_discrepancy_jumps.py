#!/usr/bin/env python3
"""
Erdős Problem #161 — Discrepancy Jumps in Hypergraph Colorings ($500)

PROBLEM: For fixed t >= 2 and n, define F^{(t)}(n, alpha) as the smallest m
such that there EXISTS a 2-coloring of the t-element subsets of [n] where
for EVERY m-element subset X, both colors appear on at least alpha * C(|X|,t)
of the t-subsets of X.

QUESTION: Does F^{(t)}(n, alpha) increase continuously in alpha, or jump?

This script computes F^{(2)}(n, alpha) exactly for small n (t=2, graph edges)
by exhaustive enumeration of all 2-colorings of E(K_n).

Known (Conlon-Fox-Sudakov 2011): For t=3, there is a single jump at alpha=0.
The t=2 case (graph edges) should be tractable computationally.

Author: Kenneth A. Mendoza (ORCID: 0009-0000-9475-5938)
Date: 2026-03-03
"""

import itertools
import json
import math
import os
import sys
import time
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np


def comb(n: int, k: int) -> int:
    """Binomial coefficient C(n,k)."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def get_edges(n: int) -> List[Tuple[int, int]]:
    """Return all edges of K_n (pairs from {0,...,n-1})."""
    return list(itertools.combinations(range(n), 2))


def get_t_subsets(elements: List[int], t: int) -> List[Tuple[int, ...]]:
    """Return all t-element subsets of the given elements."""
    return list(itertools.combinations(elements, t))


def check_coloring_for_m(
    n: int,
    t: int,
    coloring: Dict[Tuple[int, ...], int],
    m: int,
    alpha: float,
) -> bool:
    """
    Check if a given 2-coloring satisfies: for EVERY m-element subset X of [n],
    both colors appear on at least alpha * C(|X|, t) of the t-subsets of X.

    Returns True if the coloring is "good" for this (m, alpha).
    """
    vertices = list(range(n))

    # Enumerate all m-element subsets of [n]
    for subset in itertools.combinations(vertices, m):
        subset_list = list(subset)
        # Get all t-subsets of this m-subset
        t_subs = get_t_subsets(subset_list, t)
        total = len(t_subs)
        if total == 0:
            continue

        threshold = alpha * total

        # Count colors
        color_count = [0, 0]
        for ts in t_subs:
            c = coloring[ts]
            color_count[c] += 1

        # Both colors must appear at least threshold times
        if color_count[0] < threshold or color_count[1] < threshold:
            return False

    return True


def compute_F_exact(n: int, t: int, alpha: float) -> int:
    """
    Compute F^{(t)}(n, alpha) exactly by enumerating all 2-colorings.

    F^{(t)}(n, alpha) = smallest m such that THERE EXISTS a 2-coloring where
    every m-subset X has both colors on >= alpha * C(|X|,t) t-subsets.

    Strategy: For each m from t+1 down to... we want the SMALLEST m.
    Actually, we want the smallest m such that a good coloring exists.

    Note: For m < t+1, the condition is vacuous (C(m,t) might be small).
    For m = n, we need a globally balanced coloring.

    The function F is non-decreasing in alpha (harder balance = need larger m
    to be able to guarantee it). Wait — re-read the problem:

    F^{(t)}(n, alpha) is the smallest m such that we CAN 2-color...
    so that for every m-element subset X, both colors have >= alpha fraction.

    So F is the smallest m for which such a coloring EXISTS.

    Actually, re-reading more carefully: as m increases, the constraint gets
    WEAKER (fewer subsets need to be checked? No — every m-subset, and there
    are fewer m-subsets for large m). Wait:

    For m close to n, there are few m-subsets (just C(n,m) of them), so it's
    easier to satisfy. For m = n, there's only one subset (all of [n]), so
    we just need one balanced coloring globally.

    For m small, there are MANY m-subsets, and each must be balanced — harder.

    So F^{(t)}(n, alpha) = smallest m where a good coloring exists.
    As alpha increases, we need more balance, so we might need larger m
    (where there are fewer constraints).

    Hmm, but the problem says F increases from something to something as
    alpha goes from 0 to 1/2.

    Let me re-read the erdosproblems.com statement:
    "the smallest m such that we can 2-color... so that for every m-element
     subset X contains at least alpha * C(|X|,t) of each color"

    So for alpha=0: any coloring works, any m works. F(n,0) = t+1 perhaps
    (or even smaller).

    For alpha close to 1/2: we need near-perfect balance on every m-subset.
    Only possible for large m (where there are few subsets to balance).
    So F(n, 1/2-epsilon) is close to n.

    The question is: does F jump?
    """
    if alpha <= 0:
        return t + 1  # Trivially satisfied

    # Get all t-subsets
    if t == 2:
        t_subsets = get_edges(n)
    else:
        t_subsets = get_t_subsets(list(range(n)), t)

    num_t_subsets = len(t_subsets)

    # Index t-subsets for fast lookup
    t_subset_index = {ts: i for i, ts in enumerate(t_subsets)}

    # Precompute: for each m, for each m-subset, which t-subsets are contained
    # This is expensive but necessary for exact computation

    # Try each m from small to large
    # F(n, alpha) = smallest m such that a good coloring exists
    for m in range(max(t + 1, 2), n + 1):
        # Check if there exists a coloring good for this m and alpha
        # Enumerate all 2-colorings of the t-subsets
        # (2^num_t_subsets colorings — feasible only for small n)

        if num_t_subsets > 25:
            print(f"  WARNING: {num_t_subsets} t-subsets, 2^{num_t_subsets} "
                  f"colorings — this may take very long for m={m}")
            if num_t_subsets > 30:
                print(f"  SKIPPING n={n}, t={t} — too many colorings")
                return -1  # Infeasible

        # Precompute m-subsets and their t-subset indices
        m_subsets = list(itertools.combinations(range(n), m))
        m_subset_t_indices = []
        for ms in m_subsets:
            ms_list = list(ms)
            t_subs = get_t_subsets(ms_list, t)
            indices = [t_subset_index[ts] for ts in t_subs]
            m_subset_t_indices.append(indices)

        found = False
        # Enumerate all 2-colorings as bitmasks
        for coloring_bits in range(2 ** num_t_subsets):
            # Check all m-subsets
            good = True
            for indices in m_subset_t_indices:
                total = len(indices)
                threshold = alpha * total
                count_1 = sum(1 for i in indices if (coloring_bits >> i) & 1)
                count_0 = total - count_1
                if count_0 < threshold or count_1 < threshold:
                    good = False
                    break

            if good:
                found = True
                break

        if found:
            return m

    return n + 1  # No valid m found (shouldn't happen for alpha < 1/2)


def compute_F_optimized(n: int, t: int, alpha: float) -> int:
    """
    Optimized version using numpy for faster enumeration.
    For t=2: edges of K_n.
    """
    if alpha <= 0:
        return t + 1

    if t == 2:
        edges = get_edges(n)
    else:
        edges = get_t_subsets(list(range(n)), t)

    num_edges = len(edges)
    edge_index = {e: i for i, e in enumerate(edges)}

    if num_edges > 28:
        print(f"  n={n}, t={t}: {num_edges} edges — too many for exhaustive "
              f"search (2^{num_edges}). Using sampling.")
        return compute_F_sampled(n, t, alpha, num_samples=500000)

    # For each m, precompute constraint matrix
    for m in range(max(t + 1, 2), n + 1):
        m_subsets = list(itertools.combinations(range(n), m))

        # For each m-subset, compute which edge indices are involved
        # and the threshold
        constraints = []
        for ms in m_subsets:
            ms_list = list(ms)
            t_subs = get_t_subsets(ms_list, t)
            indices = [edge_index[ts] for ts in t_subs]
            total = len(indices)
            threshold = alpha * total
            constraints.append((indices, total, threshold))

        # Enumerate all colorings
        found = False
        total_colorings = 2 ** num_edges

        # Use numpy for batch checking
        # Process in chunks to manage memory
        chunk_size = min(total_colorings, 2 ** 16)

        for chunk_start in range(0, total_colorings, chunk_size):
            chunk_end = min(chunk_start + chunk_size, total_colorings)
            actual_size = chunk_end - chunk_start

            # Generate coloring bits as array
            coloring_array = np.arange(chunk_start, chunk_end, dtype=np.int64)

            # For each edge position, extract bit
            bits = np.zeros((actual_size, num_edges), dtype=np.int8)
            for j in range(num_edges):
                bits[:, j] = (coloring_array >> j) & 1

            # Check all constraints for all colorings in chunk
            valid = np.ones(actual_size, dtype=bool)

            for indices, total, threshold in constraints:
                if not np.any(valid):
                    break
                # Sum color-1 for these edge indices
                count_1 = bits[valid][:, indices].sum(axis=1)
                count_0 = total - count_1
                # Both must be >= threshold
                constraint_ok = (count_0 >= threshold) & (count_1 >= threshold)
                # Update valid mask
                valid_indices = np.where(valid)[0]
                valid[valid_indices[~constraint_ok]] = False

            if np.any(valid):
                found = True
                # Find the first valid coloring for diagnostics
                first_valid = np.where(valid)[0][0] + chunk_start
                break

        if found:
            return m

    return n + 1


def compute_F_sampled(n: int, t: int, alpha: float,
                      num_samples: int = 500000) -> int:
    """
    For larger n: sample random 2-colorings and find the best achievable m.
    Returns an UPPER BOUND on F^{(t)}(n, alpha).
    """
    if alpha <= 0:
        return t + 1

    if t == 2:
        edges = get_edges(n)
    else:
        edges = get_t_subsets(list(range(n)), t)

    num_edges = len(edges)
    edge_index = {e: i for i, e in enumerate(edges)}

    # Precompute for all m
    all_m_data = {}
    for m in range(max(t + 1, 2), n + 1):
        m_subsets = list(itertools.combinations(range(n), m))
        constraints = []
        for ms in m_subsets:
            ms_list = list(ms)
            t_subs = get_t_subsets(ms_list, t)
            indices = [edge_index[ts] for ts in t_subs]
            total = len(indices)
            threshold = alpha * total
            constraints.append((np.array(indices), total, threshold))
        all_m_data[m] = constraints

    best_m = n + 1  # Start with worst case

    for _ in range(num_samples):
        # Random coloring
        coloring = np.random.randint(0, 2, size=num_edges, dtype=np.int8)

        # Find smallest m for which this coloring works
        for m in range(max(t + 1, 2), n + 1):
            good = True
            for indices, total, threshold in all_m_data[m]:
                count_1 = coloring[indices].sum()
                count_0 = total - count_1
                if count_0 < threshold or count_1 < threshold:
                    good = False
                    break
            if good:
                best_m = min(best_m, m)
                break

    return best_m


def compute_F_smart(n: int, t: int, alpha: float) -> Tuple[int, Optional[int]]:
    """
    Smart solver that uses exact computation for small cases
    and sampling for larger ones. Returns (F_value, coloring_bits_or_None).
    """
    if t == 2:
        edges = get_edges(n)
    else:
        edges = get_t_subsets(list(range(n)), t)

    num_edges = len(edges)

    if num_edges <= 21:
        # Exact computation feasible
        val = compute_F_optimized(n, t, alpha)
        return (val, None)
    else:
        # Use sampling — returns upper bound
        val = compute_F_sampled(n, t, alpha)
        return (val, None)


def run_t2_sweep():
    """
    Main computation: sweep F^{(2)}(n, alpha) for n=4..8 and many alpha values.
    """
    print("=" * 70)
    print("ERDŐS PROBLEM #161 — DISCREPANCY JUMPS")
    print("Computing F^{(2)}(n, α) for t=2 (graph edges)")
    print("=" * 70)

    # Alpha values to test — fine-grained to detect jumps
    alphas_coarse = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35,
                     0.40, 0.45, 0.49]
    alphas_fine = np.arange(0.0, 0.50, 0.01).tolist()

    # For exact computation: n=4,5,6 (feasible)
    # For sampling: n=7,8
    n_values = [4, 5, 6, 7, 8]

    results = {}

    for n in n_values:
        edges = get_edges(n)
        num_edges = len(edges)
        exact = num_edges <= 21
        method = "EXACT" if exact else "SAMPLED (upper bound)"

        print(f"\n{'─' * 60}")
        print(f"n = {n}, |E(K_n)| = {num_edges}, method = {method}")
        print(f"{'─' * 60}")

        # Use fine alphas for small n, coarse for large n
        if exact:
            alphas = alphas_fine
        else:
            alphas = alphas_coarse

        F_values = {}
        start_time = time.time()

        for alpha in alphas:
            alpha_key = round(alpha, 4)
            F_val, _ = compute_F_smart(n, 2, alpha)
            F_values[alpha_key] = F_val
            elapsed = time.time() - start_time
            print(f"  α={alpha_key:6.3f}  →  F^(2)({n}, α) = {F_val}"
                  f"  [{elapsed:.1f}s elapsed]")

        results[n] = {
            "n": n,
            "t": 2,
            "num_edges": num_edges,
            "method": method,
            "F_values": {str(k): v for k, v in F_values.items()},
        }

        # Detect jumps
        sorted_alphas = sorted(F_values.keys())
        jumps = []
        for i in range(1, len(sorted_alphas)):
            a_prev = sorted_alphas[i - 1]
            a_curr = sorted_alphas[i]
            f_prev = F_values[a_prev]
            f_curr = F_values[a_curr]
            if f_curr > f_prev:
                jumps.append({
                    "alpha_from": a_prev,
                    "alpha_to": a_curr,
                    "F_from": f_prev,
                    "F_to": f_curr,
                    "jump_size": f_curr - f_prev,
                })

        results[n]["jumps"] = jumps

        print(f"\n  JUMPS DETECTED for n={n}:")
        if jumps:
            for j in jumps:
                print(f"    α ∈ ({j['alpha_from']:.3f}, {j['alpha_to']:.3f}]: "
                      f"F jumps from {j['F_from']} to {j['F_to']} "
                      f"(Δ = {j['jump_size']})")
        else:
            print(f"    No jumps — F is constant at {F_values[sorted_alphas[0]]}")

    return results


def analyze_results(results: dict) -> dict:
    """Analyze results for patterns and jump structure."""
    analysis = {
        "problem": "Erdős #161",
        "description": "Discrepancy jumps in 2-colorings of t-uniform hypergraphs",
        "t": 2,
        "date": "2026-03-03",
        "findings": [],
    }

    for n, data in sorted(results.items()):
        F_values = {float(k): v for k, v in data["F_values"].items()}
        sorted_alphas = sorted(F_values.keys())

        # Distinct F values
        distinct_F = sorted(set(F_values.values()))

        # Alpha ranges for each F value
        F_ranges = defaultdict(list)
        for a in sorted_alphas:
            F_ranges[F_values[a]].append(a)

        plateaus = []
        for f_val in distinct_F:
            a_list = F_ranges[f_val]
            plateaus.append({
                "F_value": f_val,
                "alpha_min": min(a_list),
                "alpha_max": max(a_list),
                "alpha_range_width": max(a_list) - min(a_list),
            })

        finding = {
            "n": n,
            "num_edges": data["num_edges"],
            "method": data["method"],
            "distinct_F_values": distinct_F,
            "num_distinct_values": len(distinct_F),
            "plateaus": plateaus,
            "jumps": data["jumps"],
            "has_jumps": len(data["jumps"]) > 0,
            "num_jumps": len(data["jumps"]),
        }

        # Key observation: is F a step function?
        if len(distinct_F) <= 3 and len(data["jumps"]) <= 2:
            finding["pattern"] = "STEP_FUNCTION"
            finding["interpretation"] = (
                f"F^(2)({n}, α) appears to be a step function with "
                f"{len(data['jumps'])} jump(s)"
            )
        elif len(distinct_F) > 3:
            finding["pattern"] = "MULTI_STEP_OR_CONTINUOUS"
            finding["interpretation"] = (
                f"F^(2)({n}, α) takes {len(distinct_F)} distinct values — "
                f"either many jumps or approaching continuity"
            )

        analysis["findings"].append(finding)

    return analysis


def plot_results(results: dict, save_dir: str):
    """Generate plots of F^{(2)}(n, alpha) vs alpha."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plots")
        return

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    colors_map = {4: '#1f77b4', 5: '#ff7f0e', 6: '#2ca02c',
                  7: '#d62728', 8: '#9467bd'}

    for idx, (n, data) in enumerate(sorted(results.items())):
        if idx >= 5:
            break

        ax = axes[idx]
        F_values = {float(k): v for k, v in data["F_values"].items()}
        sorted_alphas = sorted(F_values.keys())
        F_list = [F_values[a] for a in sorted_alphas]

        ax.step(sorted_alphas, F_list, where='post',
                color=colors_map.get(n, 'black'), linewidth=2)
        ax.scatter(sorted_alphas, F_list, color=colors_map.get(n, 'black'),
                   s=20, zorder=5)

        # Mark jumps
        for jump in data.get("jumps", []):
            ax.axvline(x=jump["alpha_to"], color='red', linestyle='--',
                       alpha=0.5, linewidth=1)
            ax.annotate(
                f'Δ={jump["jump_size"]}',
                xy=(jump["alpha_to"], jump["F_to"]),
                xytext=(jump["alpha_to"] + 0.02, jump["F_to"] + 0.2),
                fontsize=8, color='red',
            )

        ax.set_xlabel('α', fontsize=12)
        ax.set_ylabel(f'F^(2)({n}, α)', fontsize=12)
        ax.set_title(f'n = {n} ({data["method"]})', fontsize=13)
        ax.set_xlim(-0.02, 0.52)
        ax.set_ylim(0, n + 1)
        ax.grid(True, alpha=0.3)

    # Combined plot on last subplot
    ax = axes[5] if len(axes) > 5 else axes[-1]
    for n, data in sorted(results.items()):
        F_values = {float(k): v for k, v in data["F_values"].items()}
        sorted_alphas = sorted(F_values.keys())
        F_list = [F_values[a] for a in sorted_alphas]

        label = f'n={n}'
        ax.step(sorted_alphas, F_list, where='post',
                color=colors_map.get(n, 'black'), linewidth=2, label=label)

    ax.set_xlabel('α', fontsize=12)
    ax.set_ylabel('F^(2)(n, α)', fontsize=12)
    ax.set_title('All n values combined', fontsize=13)
    ax.set_xlim(-0.02, 0.52)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.suptitle(
        'Erdős Problem #161: Discrepancy Jumps — F^{(2)}(n, α) vs α\n'
        'Does F increase continuously or does it have jumps?',
        fontsize=14, fontweight='bold',
    )
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    plot_path = os.path.join(save_dir, 'erdos161_discrepancy_jumps_t2.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {plot_path}")
    plt.close()


def main():
    # Output directories
    results_dir = "/Users/kenbengoetxea/container-projects/apps/H2/Research-Hub/erdos/results"
    os.makedirs(results_dir, exist_ok=True)

    # Also save solver to erdos-161 if possible
    solver_dir = "/Users/kenbengoetxea/container-projects/apps/H2/Math-Problems/erdos-161"
    os.makedirs(solver_dir, exist_ok=True)

    # Run computation
    print("\nStarting computation...")
    start = time.time()
    results = run_t2_sweep()
    elapsed = time.time() - start
    print(f"\nTotal computation time: {elapsed:.1f}s")

    # Analyze
    print("\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)
    analysis = analyze_results(results)

    for finding in analysis["findings"]:
        print(f"\nn = {finding['n']}:")
        print(f"  Distinct F values: {finding['distinct_F_values']}")
        print(f"  Number of jumps: {finding['num_jumps']}")
        print(f"  Pattern: {finding.get('pattern', 'N/A')}")
        print(f"  Interpretation: {finding.get('interpretation', 'N/A')}")
        if finding['plateaus']:
            print(f"  Plateaus:")
            for p in finding['plateaus']:
                print(f"    F={p['F_value']}: α ∈ [{p['alpha_min']:.3f}, "
                      f"{p['alpha_max']:.3f}] (width {p['alpha_range_width']:.3f})")

    # Save results
    output = {
        "experiment_id": "EXP-ERDOS161-DISCREPANCY-01",
        "problem": "Erdős #161 — Discrepancy Jumps",
        "prize": "$500",
        "date": "2026-03-03",
        "parameters": {
            "t": 2,
            "n_values": [4, 5, 6, 7, 8],
            "alpha_range": "[0, 0.49]",
        },
        "computation_time_seconds": round(elapsed, 1),
        "raw_results": {
            str(n): {
                "F_values": data["F_values"],
                "jumps": data["jumps"],
                "method": data["method"],
            }
            for n, data in results.items()
        },
        "analysis": analysis,
    }

    results_path = os.path.join(
        results_dir, "EXP-ERDOS161-DISCREPANCY-01_RESULTS.json"
    )
    with open(results_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults saved to: {results_path}")

    # Generate plots
    plot_results(results, results_dir)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY — ERDŐS PROBLEM #161")
    print("=" * 70)

    any_jumps = any(
        f["has_jumps"] for f in analysis["findings"]
    )
    all_jumps = all(
        f["has_jumps"] for f in analysis["findings"]
    )

    if all_jumps:
        print("\nRESULT: F^{(2)}(n, α) exhibits JUMPS for ALL tested n values.")
        print("This is a STEP FUNCTION in α, not continuous.")
        print("\nThis provides computational evidence that the answer to")
        print("Erdős's question is: F has JUMPS (at least for t=2).")
    elif any_jumps:
        print("\nRESULT: F^{(2)}(n, α) exhibits JUMPS for SOME n values.")
        print("The behavior may depend on n.")
    else:
        print("\nRESULT: No jumps detected — F appears continuous.")

    print("\nJump locations by n:")
    for finding in analysis["findings"]:
        n = finding["n"]
        if finding["has_jumps"]:
            for j in finding["jumps"]:
                print(f"  n={n}: jump at α ∈ ({j['alpha_from']:.3f}, "
                      f"{j['alpha_to']:.3f}], "
                      f"F: {j['F_from']} → {j['F_to']}")
        else:
            print(f"  n={n}: no jumps detected")


if __name__ == "__main__":
    main()
