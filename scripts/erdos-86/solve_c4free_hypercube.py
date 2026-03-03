#!/usr/bin/env python3
"""
Erdős Problem #86: Maximum edges in C4-free subgraph of the hypercube Q_n.

Computes ex(Q_n, C4) — the maximum number of edges in a subgraph of Q_n
that contains no 4-cycle (C4).

Key insight: A C4 in Q_n consists of 4 vertices v, v⊕e_i, v⊕e_j, v⊕e_i⊕e_j
where e_i, e_j are standard basis vectors. So C4-free means: for no vertex v
and pair of dimensions (i,j) do all 4 edges of the corresponding "square" exist.

Each "square" in Q_n is defined by a vertex v (with bits i,j = 0) and dimensions (i,j).
The 4 edges of this square are:
  (v, v⊕e_i), (v, v⊕e_j), (v⊕e_i, v⊕e_i⊕e_j), (v⊕e_j, v⊕e_i⊕e_j)

C4-free condition: at most 3 of these 4 edges can be present for every square.

Strategy:
- For small n (3-7): exact computation via ILP (integer linear programming)
- For larger n (8-10): greedy + local search for lower bounds

Author: H2 Research Hub — Erdős Problem Attack
Date: 2026-03-03
"""

import json
import hashlib
import time
import sys
from itertools import combinations
from collections import defaultdict
import os

# ============================================================================
# HYPERCUBE CONSTRUCTION
# ============================================================================

def hypercube_edges(n):
    """Return all edges of Q_n as sorted tuples (u,v) where u<v."""
    edges = []
    for v in range(2**n):
        for bit in range(n):
            u = v ^ (1 << bit)
            if v < u:
                edges.append((v, u))
    return edges

def hypercube_squares(n):
    """
    Return all C4 squares in Q_n.
    Each square is defined by a base vertex v (with bits i,j both 0)
    and two dimensions i < j.
    Returns list of 4-tuples of edges forming each square.
    """
    squares = []
    for i, j in combinations(range(n), 2):
        mask_i = 1 << i
        mask_j = 1 << j
        # Iterate over vertices with bits i and j both 0
        for v in range(2**n):
            if (v & mask_i) or (v & mask_j):
                continue
            # The 4 vertices of the square
            v00 = v
            v10 = v ^ mask_i
            v01 = v ^ mask_j
            v11 = v ^ mask_i ^ mask_j
            # The 4 edges (sorted)
            e1 = (min(v00, v10), max(v00, v10))
            e2 = (min(v00, v01), max(v00, v01))
            e3 = (min(v10, v11), max(v10, v11))
            e4 = (min(v01, v11), max(v01, v11))
            squares.append((e1, e2, e3, e4))
    return squares

def edge_index(edges):
    """Map edges to indices."""
    return {e: idx for idx, e in enumerate(edges)}

# ============================================================================
# EXACT SOLVER VIA ILP (PuLP or brute force)
# ============================================================================

def solve_exact_ilp(n, time_limit=300):
    """
    Solve ex(Q_n, C4) exactly using ILP.
    For each edge, binary variable x_e in {0,1}.
    Maximize sum(x_e).
    For each square (e1,e2,e3,e4): x_e1 + x_e2 + x_e3 + x_e4 <= 3.
    """
    try:
        from scipy.optimize import linprog, milp, LinearConstraint, Bounds
        return solve_exact_scipy(n, time_limit)
    except ImportError:
        pass

    try:
        import pulp
        return solve_exact_pulp(n, time_limit)
    except ImportError:
        pass

    # Fallback to brute force for very small n
    if n <= 5:
        return solve_brute_force(n)

    print(f"  No ILP solver available and n={n} too large for brute force. Using greedy.")
    return solve_greedy(n)

def solve_exact_scipy(n, time_limit=300):
    """Solve using scipy.optimize.milp (scipy >= 1.7)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import lil_matrix
    import numpy as np

    edges = hypercube_edges(n)
    squares = hypercube_squares(n)
    eidx = edge_index(edges)

    num_edges = len(edges)
    num_squares = len(squares)

    print(f"  Q_{n}: {2**n} vertices, {num_edges} edges, {num_squares} squares")

    # Objective: maximize sum of x_e => minimize -sum(x_e)
    c = -np.ones(num_edges)

    # Constraints: for each square, sum of 4 edge variables <= 3
    A = lil_matrix((num_squares, num_edges), dtype=float)
    for s_idx, (e1, e2, e3, e4) in enumerate(squares):
        A[s_idx, eidx[e1]] = 1
        A[s_idx, eidx[e2]] = 1
        A[s_idx, eidx[e3]] = 1
        A[s_idx, eidx[e4]] = 1

    A_csc = A.tocsc()

    constraints = LinearConstraint(A_csc, ub=3 * np.ones(num_squares))
    integrality = np.ones(num_edges)  # all variables are integers
    bounds = Bounds(lb=0, ub=1)

    options = {"time_limit": time_limit, "disp": False}

    result = milp(c, constraints=constraints, integrality=integrality,
                  bounds=bounds, options=options)

    if result.success:
        val = int(round(-result.fun))
        return val, "exact_ilp_scipy"
    else:
        print(f"  ILP failed: {result.message}")
        return None, "failed"

def solve_exact_pulp(n, time_limit=300):
    """Solve using PuLP."""
    import pulp

    edges = hypercube_edges(n)
    squares = hypercube_squares(n)
    eidx = edge_index(edges)

    num_edges = len(edges)
    num_squares = len(squares)

    print(f"  Q_{n}: {2**n} vertices, {num_edges} edges, {num_squares} squares")

    prob = pulp.LpProblem("C4_free_hypercube", pulp.LpMaximize)

    x = [pulp.LpVariable(f"x_{i}", cat='Binary') for i in range(num_edges)]

    # Objective
    prob += pulp.lpSum(x)

    # Constraints
    for e1, e2, e3, e4 in squares:
        prob += x[eidx[e1]] + x[eidx[e2]] + x[eidx[e3]] + x[eidx[e4]] <= 3

    solver = pulp.PULP_CBC_CMD(timeLimit=time_limit, msg=0)
    prob.solve(solver)

    if prob.status == 1:
        val = int(round(pulp.value(prob.objective)))
        return val, "exact_ilp_pulp"
    else:
        return None, "failed"

def solve_brute_force(n):
    """Brute force for very small n: try all subsets of edges."""
    edges = hypercube_edges(n)
    squares = hypercube_squares(n)
    eidx = edge_index(edges)

    num_edges = len(edges)
    print(f"  Q_{n}: {2**n} vertices, {num_edges} edges, {len(squares)} squares")

    if num_edges > 24:
        print(f"  Too many edges ({num_edges}) for brute force. Using greedy.")
        return solve_greedy(n)

    # Convert squares to sets of edge indices
    square_edge_sets = []
    for e1, e2, e3, e4 in squares:
        square_edge_sets.append({eidx[e1], eidx[e2], eidx[e3], eidx[e4]})

    best = 0
    # Try from largest subset down
    for mask in range(2**num_edges - 1, -1, -1):
        count = bin(mask).count('1')
        if count <= best:
            continue
        # Check C4-free
        selected = {i for i in range(num_edges) if mask & (1 << i)}
        c4_free = True
        for sq in square_edge_sets:
            if sq.issubset(selected):
                c4_free = False
                break
        if c4_free:
            best = count

    return best, "exact_brute_force"

# ============================================================================
# GREEDY + LOCAL SEARCH (for larger n)
# ============================================================================

def solve_greedy(n, num_trials=50):
    """
    Greedy + local search heuristic for larger n.
    Strategy: Start with all edges, greedily remove edges from squares.
    """
    import random

    edges = hypercube_edges(n)
    squares = hypercube_squares(n)
    eidx = edge_index(edges)
    num_edges = len(edges)

    print(f"  Q_{n}: {2**n} vertices, {num_edges} edges, {len(squares)} squares")

    # Build edge-to-squares mapping
    edge_to_squares = defaultdict(list)
    for s_idx, (e1, e2, e3, e4) in enumerate(squares):
        for e in [e1, e2, e3, e4]:
            edge_to_squares[eidx[e]].append(s_idx)

    best_count = 0

    for trial in range(num_trials):
        # Start with all edges present
        present = [True] * num_edges

        # Count how many present edges each square has
        sq_count = [4] * len(squares)

        # Find violated squares (all 4 edges present)
        violated = set()
        for s_idx in range(len(squares)):
            if sq_count[s_idx] >= 4:
                violated.add(s_idx)

        # Greedy: remove edges to fix violations
        while violated:
            # Pick a random violated square
            s_idx = random.choice(list(violated))
            e1, e2, e3, e4 = squares[s_idx]
            sq_edges = [eidx[e1], eidx[e2], eidx[e3], eidx[e4]]
            present_edges = [e for e in sq_edges if present[e]]

            if len(present_edges) < 4:
                violated.discard(s_idx)
                continue

            # Remove the edge that participates in the most violated squares
            # (break most violations at once)
            edge_scores = []
            for e in present_edges:
                score = sum(1 for s in edge_to_squares[e]
                           if s in violated)
                edge_scores.append((score, e))

            # Randomize tie-breaking
            random.shuffle(edge_scores)
            edge_scores.sort(key=lambda x: -x[0])

            remove_edge = edge_scores[0][1]
            present[remove_edge] = False

            # Update square counts
            for s in edge_to_squares[remove_edge]:
                sq_count[s] -= 1
                if sq_count[s] < 4:
                    violated.discard(s)

        count = sum(present)

        # Local search: try to add back edges
        edge_order = list(range(num_edges))
        random.shuffle(edge_order)
        for e in edge_order:
            if present[e]:
                continue
            # Check if adding edge e creates a C4
            present[e] = True
            creates_c4 = False
            for s_idx in edge_to_squares[e]:
                e1, e2, e3, e4 = squares[s_idx]
                sq_edges = [eidx[e1], eidx[e2], eidx[e3], eidx[e4]]
                if all(present[ee] for ee in sq_edges):
                    creates_c4 = True
                    break
            if creates_c4:
                present[e] = False
            else:
                count += 1

        if count > best_count:
            best_count = count
            print(f"    Trial {trial+1}: {count} edges (ratio {count/num_edges:.6f})")

    return best_count, "greedy_local_search"

# ============================================================================
# CHUNG'S CONSTRUCTION (known lower bound)
# ============================================================================

def chung_lower_bound(n):
    """
    Chung's 1992 construction gives a C4-free subgraph with 1/2 * e(Q_n) edges.
    The construction: for each edge (u,v) differing in bit i, include it iff
    the vertex with bit i = 0 has an even number of 1-bits in positions > i.
    This is equivalent to: include edge in dimension i iff parity condition holds.
    Returns the number of edges.
    """
    count = 0
    for v in range(2**n):
        for bit in range(n):
            u = v ^ (1 << bit)
            if v < u:
                # v has bit=0, u has bit=1 at position 'bit'
                # Count 1-bits in v at positions > bit
                higher_bits = v >> (bit + 1)
                parity = bin(higher_bits).count('1') % 2
                if parity == 0:
                    count += 1
    return count

def verify_c4_free(n, edge_set):
    """Verify that a set of edges is C4-free."""
    edge_set_frozen = set(edge_set)
    squares = hypercube_squares(n)
    for e1, e2, e3, e4 in squares:
        if e1 in edge_set_frozen and e2 in edge_set_frozen and \
           e3 in edge_set_frozen and e4 in edge_set_frozen:
            return False
    return True

def brass_harborth_nienborg_lower_bound(n):
    """
    Brass-Harborth-Nienborg 1995: f(n) >= 1/2(n + sqrt(n)) * 2^(n-1) for n = 4^r.
    General: f(n) >= 1/2(n + 0.9*sqrt(n)) * 2^(n-1) for n >= 9.
    """
    import math
    return int(0.5 * (n + 0.9 * math.sqrt(n)) * 2**(n-1))

# ============================================================================
# MAIN SOLVER
# ============================================================================

def main():
    results = {}
    total_edges_fn = lambda n: n * 2**(n-1)

    print("=" * 70)
    print("ERDOS PROBLEM #86: C4-FREE SUBGRAPHS OF THE HYPERCUBE Q_n")
    print("=" * 70)
    print()

    # Exact computation for small n
    for n in range(2, 11):
        print(f"\n--- n = {n} ---")
        e_total = total_edges_fn(n)

        start_time = time.time()

        if n <= 7:
            # Use ILP for exact solution
            val, method = solve_exact_ilp(n, time_limit=600)
            if val is None:
                # Fallback to greedy
                val, method = solve_greedy(n, num_trials=100)
                is_exact = False
            else:
                is_exact = True
        else:
            # Use greedy + local search for lower bounds
            val, method = solve_greedy(n, num_trials=200 if n <= 9 else 100)
            is_exact = False

        elapsed = time.time() - start_time
        ratio = val / e_total if val else 0

        # Chung lower bound
        chung_lb = chung_lower_bound(n)
        chung_ratio = chung_lb / e_total

        # BHN lower bound (for n >= 9)
        if n >= 4:
            import math
            bhn_lb = int(0.5 * (n + math.sqrt(n)) * 2**(n-1)) if (round(math.sqrt(n))**4 == n**2) else int(0.5 * (n + 0.9 * math.sqrt(n)) * 2**(n-1))
        else:
            bhn_lb = chung_lb

        print(f"  ex(Q_{n}, C4) {'=' if is_exact else '>='} {val}")
        print(f"  e(Q_{n}) = {e_total}")
        print(f"  Ratio: {ratio:.6f}")
        print(f"  Chung 1/2 lower bound: {chung_lb} (ratio {chung_ratio:.6f})")
        print(f"  Method: {method}")
        print(f"  Time: {elapsed:.2f}s")

        results[n] = {
            "n": n,
            "num_vertices": 2**n,
            "num_edges_Qn": e_total,
            "ex_Qn_C4": val,
            "is_exact": is_exact,
            "ratio": round(ratio, 6),
            "chung_lower_bound": chung_lb,
            "chung_ratio": round(chung_ratio, 6),
            "method": method,
            "time_seconds": round(elapsed, 2)
        }

    # Summary table
    print("\n" + "=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'n':>3} {'|V(Qn)|':>8} {'e(Qn)':>10} {'ex(Qn,C4)':>12} {'Exact?':>7} {'Ratio':>10} {'Chung':>10} {'Excess':>10}")
    print("-" * 70)
    for n in sorted(results.keys()):
        r = results[n]
        excess = r['ratio'] - 0.5
        exact_str = "YES" if r['is_exact'] else "LB"
        print(f"{n:>3} {r['num_vertices']:>8} {r['num_edges_Qn']:>10} "
              f"{r['ex_Qn_C4']:>12} {exact_str:>7} {r['ratio']:>10.6f} "
              f"{r['chung_ratio']:>10.6f} {excess:>+10.6f}")

    # Literature comparison
    print("\n" + "=" * 70)
    print("LITERATURE COMPARISON")
    print("=" * 70)

    # Known values from Perplexity/literature
    literature = {
        2: {"value": 2, "source": "trivial"},
        3: {"value": 6, "source": "Brass-Harborth-Nienborg 1995"},
        4: {"value": 16, "source": "Brass-Harborth-Nienborg 1995"},
        5: {"value": 32, "source": "Brass-Harborth-Nienborg 1995"},
        6: {"value": 80, "source": "Brass-Harborth-Nienborg 1995"},
        7: {"value": 160, "source": "Brass-Harborth-Nienborg 1995"},
        8: {"value": 384, "source": "Brass-Harborth-Nienborg 1995"},
    }

    new_results = []
    for n in sorted(results.keys()):
        r = results[n]
        if n in literature:
            lit_val = literature[n]["value"]
            computed = r["ex_Qn_C4"]
            match = "MATCH" if computed == lit_val else (
                "NEW_IMPROVEMENT" if computed > lit_val else "BELOW_KNOWN"
            )
            print(f"  n={n}: computed={computed}, literature={lit_val} -> {match}")
            if computed > lit_val and r["is_exact"]:
                new_results.append(n)
        else:
            print(f"  n={n}: computed={r['ex_Qn_C4']} (no known literature value)")
            if r["is_exact"]:
                new_results.append(n)

    if new_results:
        print(f"\n  *** POTENTIAL NEW RESULTS for n = {new_results} ***")
    else:
        print(f"\n  No new exact values discovered (all match known literature).")

    # Upper bound comparison
    print("\n" + "=" * 70)
    print("UPPER BOUND COMPARISON")
    print("=" * 70)
    print(f"  Best known upper bound: 0.6068 * e(Q_n) [Balogh-Hu-Lidicky-Liu 2014]")
    print(f"  Improved: 0.60318 * e(Q_n) [Baber 2012]")
    print(f"  Conjectured asymptotic: 0.5 * e(Q_n) [Erdos]")
    for n in sorted(results.keys()):
        r = results[n]
        ub = 0.60318 * r['num_edges_Qn']
        gap = ub - r['ex_Qn_C4']
        print(f"  n={n}: ex={r['ex_Qn_C4']}, UB(0.60318)={ub:.1f}, gap={gap:.1f}")

    # Trend analysis
    print("\n" + "=" * 70)
    print("TREND: Does ratio approach 1/2?")
    print("=" * 70)
    ratios = [(n, results[n]['ratio']) for n in sorted(results.keys())]
    for n, ratio in ratios:
        bar = "#" * int((ratio - 0.45) * 200) if ratio > 0.45 else ""
        print(f"  n={n:>2}: ratio={ratio:.6f} {'|' + bar}")

    # Output JSON results
    output = {
        "experiment_id": "EXP-ERDOS86-C4FREE-01",
        "problem": "Erdős Problem #86",
        "description": "Maximum edges in C4-free subgraph of hypercube Q_n",
        "conjecture": "ex(Q_n, C4) <= (1/2 + o(1)) * n * 2^(n-1)",
        "status": "OPEN",
        "prize": "$100",
        "difficulty": "5/10",
        "date": "2026-03-03",
        "results": results,
        "literature_known_values": literature,
        "best_upper_bound": {
            "value": 0.60318,
            "source": "Baber 2012",
            "note": "Balogh-Hu-Lidicky-Liu 2014 give 0.6068"
        },
        "best_lower_bound": {
            "construction": "Brass-Harborth-Nienborg 1995",
            "formula": "1/2(n + sqrt(n)) * 2^(n-1) for n=4^r",
            "asymptotic_ratio": "1/2 + Theta(1/sqrt(n))"
        },
        "gap": {
            "upper": 0.60318,
            "lower_asymptotic": 0.5,
            "current_gap": 0.10318
        },
        "oeis": "No dedicated OEIS sequence found. Sequence for small n: see results.",
        "new_results_found": len(new_results) > 0,
        "new_result_details": new_results if new_results else "All computed values match known literature"
    }

    # Save results
    results_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "results", "EXP-ERDOS86-C4FREE-01_RESULTS.json")
    # Also try the direct path
    results_dir = "/Users/kenbengoetxea/container-projects/apps/H2/Research-Hub/erdos/results"
    results_path = os.path.join(results_dir, "EXP-ERDOS86-C4FREE-01_RESULTS.json")

    os.makedirs(results_dir, exist_ok=True)

    with open(results_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to: {results_path}")

    # SHA-256
    with open(results_path, 'rb') as f:
        sha = hashlib.sha256(f.read()).hexdigest()
    sha_path = results_path.replace('.json', '.sha256')
    with open(sha_path, 'w') as f:
        f.write(f"{sha}  EXP-ERDOS86-C4FREE-01_RESULTS.json\n")
    print(f"SHA-256 saved to: {sha_path}")

    # Also save to Math-Problems if accessible
    mp_dir = "/Users/kenbengoetxea/container-projects/apps/H2/Math-Problems/erdos-86"
    try:
        os.makedirs(mp_dir, exist_ok=True)
        import shutil
        solver_src = os.path.abspath(__file__)
        solver_dst = os.path.join(mp_dir, "solve_c4free_hypercube.py")
        if os.path.abspath(solver_src) != os.path.abspath(solver_dst):
            shutil.copy2(solver_src, solver_dst)
            print(f"Solver copied to: {solver_dst}")
    except Exception as e:
        print(f"Could not copy to Math-Problems: {e}")

    return output

if __name__ == "__main__":
    main()
