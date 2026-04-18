#!/usr/bin/env python3
"""
EXP-MATH-ERDOS30-SIDON-007: Large-Window Transfer Matrix via Sparse Power Iteration

Push the Sidon transfer matrix to W=15-20+ using:
1. Backtracking enumeration of Sidon states (no brute-force 2^W)
2. Sparse matrix representation (most entries are 0)
3. Power iteration for λ_max (only need the largest eigenvalue)
4. Finite-size scaling to extract critical exponents

No external dependencies beyond numpy.
"""

import json
import math
import time
from pathlib import Path
from collections import defaultdict

import numpy as np


# ─── Efficient Sidon State Enumeration ────────────────────────────────────

def enumerate_sidon_states(W):
    """
    Enumerate all Sidon subsets of {0, 1, ..., W-1} using backtracking.
    A Sidon subset has all pairwise differences distinct.
    Returns list of frozensets.
    """
    states = [frozenset()]  # empty set is valid

    def backtrack(start, current, used_diffs):
        if current:
            states.append(frozenset(current))
        for x in range(start, W):
            new_diffs = []
            valid = True
            for a in current:
                d = x - a
                if d in used_diffs:
                    valid = False
                    break
                new_diffs.append(d)
            if valid:
                for d in new_diffs:
                    used_diffs.add(d)
                current.append(x)
                backtrack(x + 1, current, used_diffs)
                current.pop()
                for d in new_diffs:
                    used_diffs.remove(d)

    backtrack(0, [], set())
    return states


def get_diffs(positions):
    """Get all pairwise differences from a frozenset of positions."""
    diffs = set()
    pos_list = sorted(positions)
    for i in range(len(pos_list)):
        for j in range(i + 1, len(pos_list)):
            diffs.add(pos_list[j] - pos_list[i])
    return frozenset(diffs)


# ─── Sparse Transfer Matrix Construction ─────────────────────────────────

def build_sparse_transfer_matrix(W, states, state_index):
    """
    Build sparse transfer matrix as adjacency list.
    T[i] = list of (j, weight) pairs where transition i→j is valid.

    Transition rule: shift window right by 1.
    - Each position p becomes p-1; position 0 drops out
    - New position W-1 may or may not be occupied
    """
    # Precompute shifted states and their diff sets
    n_states = len(states)
    adj = [[] for _ in range(n_states)]  # adjacency list

    for i, state_i in enumerate(states):
        # Shift: positions become {p-1 : p in state_i, p > 0}
        shifted = frozenset(p - 1 for p in state_i if p > 0)
        shifted_diffs = get_diffs(shifted)

        # Option 1: don't occupy W-1
        if shifted in state_index:
            j = state_index[shifted]
            adj[i].append(j)

        # Option 2: occupy W-1
        new_with = shifted | {W - 1}
        valid = True
        for p in shifted:
            d = (W - 1) - p
            if d in shifted_diffs:
                valid = False
                break
        # Also check new differences among themselves (actually, only one new
        # element so new-new diffs don't exist, only new-old diffs)
        if valid:
            # Check that the combined set is Sidon
            # We already checked new element vs existing, but need to check
            # that the new diffs don't collide with each other
            new_diffs_list = [(W - 1) - p for p in shifted]
            if len(new_diffs_list) == len(set(new_diffs_list)):
                # No collision among new diffs
                combined_diffs = set(shifted_diffs) | set(new_diffs_list)
                if len(combined_diffs) == len(shifted_diffs) + len(new_diffs_list):
                    if new_with in state_index:
                        j = state_index[new_with]
                        adj[i].append(j)

    return adj


def power_iteration_sparse(adj, n, max_iter=500, tol=1e-12):
    """
    Power iteration on sparse adjacency list to find λ_max.
    adj[i] = list of j indices reachable from i (all weights = 1).
    """
    # Random initial vector
    rng = np.random.RandomState(42)
    v = rng.random(n)
    v /= np.linalg.norm(v)

    lambda_old = 0
    for iteration in range(max_iter):
        # Matrix-vector multiply via adjacency list
        w = np.zeros(n)
        for i in range(n):
            for j in adj[i]:
                w[j] += v[i]

        # Normalize
        norm = np.linalg.norm(w)
        if norm < 1e-30:
            break
        lambda_new = norm
        v = w / norm

        if abs(lambda_new - lambda_old) < tol * abs(lambda_new):
            return lambda_new, v, iteration + 1
        lambda_old = lambda_new

    return lambda_new, v, max_iter


def find_second_eigenvalue(adj, n, v1, lambda1, max_iter=500, tol=1e-10):
    """
    Find second-largest eigenvalue by deflation: project out v1 component.
    """
    rng = np.random.RandomState(123)
    v = rng.random(n)
    v -= np.dot(v, v1) * v1
    v /= np.linalg.norm(v)

    lambda_old = 0
    for iteration in range(max_iter):
        # Matrix-vector multiply
        w = np.zeros(n)
        for i in range(n):
            for j in adj[i]:
                w[j] += v[i]

        # Project out v1 component
        w -= np.dot(w, v1) * v1

        norm = np.linalg.norm(w)
        if norm < 1e-30:
            return 0, v, iteration + 1
        lambda_new = norm
        v = w / norm

        if abs(lambda_new - lambda_old) < tol * abs(lambda_new):
            return lambda_new, v, iteration + 1
        lambda_old = lambda_new

    return lambda_new, v, max_iter


# ─── SVD-based TT-rank (for small matrices) ──────────────────────────────

def estimate_tt_rank_sparse(adj, n, max_rank=200):
    """Estimate TT-rank from sparse matrix by building dense subblock if small enough."""
    if n > max_rank:
        return {"note": f"Matrix too large ({n}) for dense SVD, skipping TT-rank"}

    # Build dense matrix
    T = np.zeros((n, n))
    for i in range(n):
        for j in adj[i]:
            T[i, j] += 1.0

    S = np.linalg.svd(T, compute_uv=False)
    total = np.sum(S)
    if total == 0:
        return {"effective_rank_99": 0}

    cumsum = np.cumsum(S) / total
    rank_99 = int(np.searchsorted(cumsum, 0.99)) + 1

    return {
        "full_rank": int(np.sum(S > 1e-10)),
        "effective_rank_99": rank_99,
        "sv_decay_rate": float(S[1] / S[0]) if S[0] > 0 and len(S) > 1 else 0,
    }


# ─── Finite-Size Scaling ─────────────────────────────────────────────────

def finite_size_scaling(Ws, lambdas):
    """
    Fit λ_max(W) to various scaling forms:
    1. λ(W) = λ_∞ + a·W^(-ν)  (power-law correction)
    2. λ(W) = λ_∞ + a·exp(-b·W)  (exponential correction)
    3. λ(W) = λ_∞ + a/(W·log(W))  (logarithmic correction)

    The critical exponent ν determines the correction to h(N):
    h(N) = √N + O(N^{1/(2ν)})
    If ν → ∞, the Erdős conjecture follows.
    """
    Ws = np.array(Ws, dtype=float)
    lambdas = np.array(lambdas, dtype=float)
    n = len(Ws)

    results = {}

    # Method 1: Fit λ = a + b·W^(-c) using Richardson extrapolation
    # Use last 3 points to estimate λ_∞
    if n >= 3:
        # Richardson: λ_∞ ≈ (λ_n · λ_{n-2} - λ_{n-1}²) / (λ_n + λ_{n-2} - 2·λ_{n-1})
        l1, l2, l3 = lambdas[-3], lambdas[-2], lambdas[-1]
        denom = l1 + l3 - 2 * l2
        if abs(denom) > 1e-10:
            lambda_inf_richardson = (l1 * l3 - l2**2) / denom
            results["richardson_lambda_inf"] = float(lambda_inf_richardson)

    # Method 2: Fit log(λ - λ_∞) ~ -ν·log(W) for various λ_∞ candidates
    best_r2 = -1
    best_nu = 0
    best_lambda_inf = lambdas[-1] * 0.9  # initial guess

    for trial_lambda in np.linspace(0.8, lambdas[-1] - 0.001, 50):
        residuals = lambdas - trial_lambda
        if np.any(residuals <= 0):
            continue
        log_res = np.log(residuals)
        log_W = np.log(Ws)

        # Linear fit: log(residual) = a - ν·log(W)
        A = np.vstack([np.ones(n), log_W]).T
        coef, res, _, _ = np.linalg.lstsq(A, log_res, rcond=None)
        nu = -coef[1]

        y_pred = coef[0] + coef[1] * log_W
        ss_res = np.sum((log_res - y_pred)**2)
        ss_tot = np.sum((log_res - np.mean(log_res))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        if r2 > best_r2:
            best_r2 = r2
            best_nu = nu
            best_lambda_inf = trial_lambda

    results["power_law_lambda_inf"] = float(best_lambda_inf)
    results["power_law_nu"] = float(best_nu)
    results["power_law_r2"] = float(best_r2)

    # Method 3: Exponential fit: λ = a + b·exp(-c·W)
    # Use log(λ_n - λ_{n-1}) differences
    if n >= 4:
        diffs = np.diff(lambdas)
        log_neg_diffs = np.log(-diffs) if np.all(diffs < 0) else None
        if log_neg_diffs is not None:
            W_mid = (Ws[:-1] + Ws[1:]) / 2
            A = np.vstack([np.ones(len(W_mid)), W_mid]).T
            coef, _, _, _ = np.linalg.lstsq(A, log_neg_diffs, rcond=None)
            exp_rate = -coef[1]
            results["exponential_decay_rate"] = float(exp_rate)

    # Interpretation
    if best_nu > 0:
        implied_correction = f"h(N) = √N + O(N^{{1/{2*best_nu:.1f}}})"
        results["implied_correction"] = implied_correction
        if best_nu > 2:
            results["erdos_support"] = "STRONG — ν > 2 means correction < N^{1/4} (better than Lindström)"
        elif best_nu > 1:
            results["erdos_support"] = "MODERATE — correction ~ N^{1/2ν}"
        else:
            results["erdos_support"] = "WEAK — correction exponent is large"

    return results


# ─── Main Experiment ─────────────────────────────────────────────────────

def run_experiment():
    results = {
        "experiment_id": "EXP-MATH-ERDOS30-SIDON-007",
        "title": "Large-Window Transfer Matrix via Sparse Power Iteration",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "transfer_matrices": [],
        "scaling": {},
    }

    print("=" * 90)
    print("EXP-MATH-ERDOS30-SIDON-007: Sparse Transfer Matrix Push to W=20+")
    print("=" * 90)

    all_Ws = []
    all_lambdas = []
    all_gap_ratios = []

    # Start from W=3, push as far as we can
    for W in range(3, 24):
        t0 = time.time()

        print(f"\n  W={W}: enumerating Sidon states...", end="", flush=True)
        states = enumerate_sidon_states(W)
        n_states = len(states)
        state_index = {s: i for i, s in enumerate(states)}
        t_enum = time.time() - t0
        print(f" {n_states} states ({t_enum:.1f}s)", end="", flush=True)

        if n_states > 500000:
            print(f" — TOO LARGE, stopping")
            break

        # Build sparse transfer matrix
        t1 = time.time()
        adj = build_sparse_transfer_matrix(W, states, state_index)
        n_edges = sum(len(a) for a in adj)
        t_build = time.time() - t1
        print(f", {n_edges} edges ({t_build:.1f}s)", end="", flush=True)

        # Power iteration for λ_max
        t2 = time.time()
        lambda_max, v1, iters1 = power_iteration_sparse(adj, n_states)
        t_power = time.time() - t2
        print(f", λ₁={lambda_max:.6f} ({iters1} iters, {t_power:.1f}s)", end="", flush=True)

        # Second eigenvalue
        t3 = time.time()
        lambda_2, v2, iters2 = find_second_eigenvalue(adj, n_states, v1, lambda_max)
        t_second = time.time() - t3
        gap_ratio = lambda_2 / lambda_max if lambda_max > 0 else 0
        print(f", λ₂={lambda_2:.6f}, λ₂/λ₁={gap_ratio:.4f} ({t_second:.1f}s)")

        # Free energy
        fe = math.log(lambda_max) / W if lambda_max > 0 else 0

        # TT-rank (only for small matrices)
        tt_info = estimate_tt_rank_sparse(adj, n_states, max_rank=300)

        record = {
            "W": W,
            "n_states": n_states,
            "n_edges": n_edges,
            "lambda_max": float(lambda_max),
            "lambda_2": float(lambda_2),
            "gap_ratio": float(gap_ratio),
            "free_energy_per_site": float(fe),
            "power_iter_count": iters1,
            "time_total_s": time.time() - t0,
        }
        if "effective_rank_99" in tt_info:
            record["tt_rank_99"] = tt_info["effective_rank_99"]

        results["transfer_matrices"].append(record)
        all_Ws.append(W)
        all_lambdas.append(float(lambda_max))
        all_gap_ratios.append(float(gap_ratio))

        # Safety: if total time per step exceeds 120s, stop
        if time.time() - t0 > 120:
            print(f"  ⚠ Step took {time.time()-t0:.0f}s — stopping to avoid timeout")
            break

    # ── Summary Table ──
    print("\n" + "=" * 90)
    print("SUMMARY TABLE")
    print("=" * 90)
    print(f"  {'W':>4s} | {'#states':>8s} | {'#edges':>8s} | {'λ_max':>10s} | {'λ₂':>10s} | "
          f"{'λ₂/λ₁':>8s} | {'f.e./site':>10s}")
    print(f"  {'-'*4}-+-{'-'*8}-+-{'-'*8}-+-{'-'*10}-+-{'-'*10}-+-{'-'*8}-+-{'-'*10}")
    for r in results["transfer_matrices"]:
        tt_str = f", TT={r.get('tt_rank_99','?')}" if 'tt_rank_99' in r else ""
        print(f"  {r['W']:4d} | {r['n_states']:8d} | {r['n_edges']:8d} | "
              f"{r['lambda_max']:10.6f} | {r['lambda_2']:10.6f} | "
              f"{r['gap_ratio']:8.4f} | {r['free_energy_per_site']:10.6f}{tt_str}")

    # ── Finite-Size Scaling ──
    print("\n" + "=" * 90)
    print("FINITE-SIZE SCALING")
    print("=" * 90)

    if len(all_Ws) >= 5:
        scaling = finite_size_scaling(all_Ws, all_lambdas)
        results["scaling"] = scaling

        print(f"\n  Richardson extrapolation λ_∞: {scaling.get('richardson_lambda_inf', 'N/A')}")
        print(f"\n  Power-law fit: λ(W) = {scaling.get('power_law_lambda_inf', '?'):.4f} "
              f"+ a·W^(-{scaling.get('power_law_nu', '?'):.3f})")
        print(f"  R² = {scaling.get('power_law_r2', 0):.4f}")
        if "exponential_decay_rate" in scaling:
            print(f"  Exponential decay rate: {scaling['exponential_decay_rate']:.4f}")
        if "implied_correction" in scaling:
            print(f"\n  Implied: {scaling['implied_correction']}")
        if "erdos_support" in scaling:
            print(f"  Erdős support: {scaling['erdos_support']}")

    # ── Gap ratio extrapolation ──
    print(f"\n  Gap ratio λ₂/λ₁ trend:")
    print(f"    {[f'{g:.4f}' for g in all_gap_ratios]}")
    if len(all_gap_ratios) >= 3:
        # Linear extrapolation: at what W does λ₂/λ₁ → 1?
        Ws_arr = np.array(all_Ws[-5:], dtype=float)
        gaps_arr = np.array(all_gap_ratios[-5:], dtype=float)
        A = np.vstack([np.ones_like(Ws_arr), Ws_arr]).T
        coef, _, _, _ = np.linalg.lstsq(A, gaps_arr, rcond=None)
        if coef[1] > 0:
            W_critical = (1.0 - coef[0]) / coef[1]
            print(f"    Linear extrapolation: λ₂/λ₁ → 1 at W ≈ {W_critical:.0f}")
            results["scaling"]["gap_closure_W"] = float(W_critical)

    # ── State space growth analysis ──
    if len(all_Ws) >= 3:
        log_states = np.log([r["n_states"] for r in results["transfer_matrices"]])
        Ws_arr = np.array(all_Ws, dtype=float)
        A = np.vstack([np.ones_like(Ws_arr), Ws_arr]).T
        coef, _, _, _ = np.linalg.lstsq(A, log_states, rcond=None)
        growth_base = math.exp(coef[1])
        print(f"\n  State space growth: ~{growth_base:.3f}^W (base = {growth_base:.4f})")
        results["scaling"]["state_space_growth_base"] = float(growth_base)

    # Save
    out_path = Path(__file__).parent / "EXP-MATH-ERDOS30-SIDON-007_RESULTS.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")

    return results


if __name__ == "__main__":
    run_experiment()
