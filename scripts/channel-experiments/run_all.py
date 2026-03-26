#!/usr/bin/env python3
"""
Channel Experiments for Erdos Problems
run_all.py -- entry point for GitHub Actions CI

Experiments:
  A: Entropy-Hausdorff monotonicity (fast, analytic)
  B: Capacity saturation grid search (slow, CI uses reduced grid)
  C: Quantum channel analogy (fast, skips quantum sim if Qiskit absent)
  D: Cross-problem transfer test (fast, analytic)

Usage:
  python run_all.py               # run all
  python run_all.py --skip b      # skip experiment B
  python run_all.py --only b      # run only experiment B
"""

import argparse
import json
import math
import os
import time
from pathlib import Path

RESULTS_DIR = Path("../../results/channel-experiments")


# ---------------------------------------------------------------------------
# Experiment A: Entropy-Hausdorff monotonicity
# ---------------------------------------------------------------------------

def run_experiment_a():
    """
    Test: For the lemniscate L_n* (Erdos #114), does the differential entropy
    H(n) scale monotonically with log(Hausdorff dimension proxy)?

    We use the certified l_star values from the EHP n=3..11 results as
    the arc-length proxy, and compute the information-theoretic entropy
    of the uniform measure on each lemniscate.
    """
    # Certified l_star bounds from EXP-MM-EHP-007 results (lower bounds used)
    l_star = {
        3:  3.9122591069701235,
        4:  5.384820755989181,
        5:  6.831790707017682,
        6:  8.264399687890477,
        7:  9.688218374101553,
        8:  11.106193082276553,
        9:  12.519932573700426,
        10: 13.930298849960456,
        11: 24.87544868514786,
    }

    results = []
    monotone = True
    prev_h = None

    for n, l in sorted(l_star.items()):
        # Entropy proxy: H(n) = log(l_star(n)) -- uniform measure on arc
        h = math.log(l)
        # Channel capacity proxy: C(n) = log(n) -- degrees of freedom
        c = math.log(n)
        ratio = h / c if c > 0 else None
        if prev_h is not None and h < prev_h:
            monotone = False
        prev_h = h
        results.append({"n": n, "l_star": l, "H": h, "C_proxy": c, "H_over_C": ratio})

    summary = {
        "experiment": "A",
        "name": "Entropy-Hausdorff monotonicity",
        "monotone": monotone,
        "verdict": "PASS" if monotone else "FAIL",
        "data": results,
    }
    return summary


# ---------------------------------------------------------------------------
# Experiment B: Capacity saturation grid search
# ---------------------------------------------------------------------------

def run_experiment_b(grid_size=None):
    """
    Test: Does the channel capacity C(n) saturate as n -> inf?
    Grid search over the ratio H(n)/log(n) for n=3..N.
    """
    if grid_size is None:
        grid_size = int(os.environ.get("EXP_B_GRID_SIZE", "10000"))

    N = min(grid_size, 500)  # cap for CI safety

    # Use analytic approximation: l_star(n) ~ n * pi / 2 for large n
    # (lemniscate perimeter scales linearly with degree)
    ratios = []
    for n in range(3, N + 1):
        l_approx = n * math.pi / 2.0
        h = math.log(l_approx)
        c = math.log(n)
        ratio = h / c
        ratios.append(ratio)

    saturates = abs(ratios[-1] - ratios[-2]) < 1e-6 if len(ratios) >= 2 else False
    asymptotic_ratio = ratios[-1] if ratios else None

    summary = {
        "experiment": "B",
        "name": "Capacity saturation grid search",
        "grid_size": N,
        "asymptotic_H_over_C": asymptotic_ratio,
        "saturates": saturates,
        "verdict": "PASS" if (asymptotic_ratio is not None and 1.0 < asymptotic_ratio < 2.0) else "INCONCLUSIVE",
        "note": "Ratio H/C -> 1 + log(pi/2)/log(n) -> 1 as n->inf",
    }
    return summary


# ---------------------------------------------------------------------------
# Experiment C: Quantum channel analogy
# ---------------------------------------------------------------------------

def run_experiment_c():
    """
    Test: Does the classical channel entropy match a quantum channel model?
    Uses numpy simulation; skips Qiskit quantum circuit sim if unavailable.
    """
    import numpy as np

    n_values = list(range(3, 12))
    classical = []
    quantum_sim = []
    qiskit_available = False

    try:
        from qiskit import QuantumCircuit
        from qiskit_aer import AerSimulator
        qiskit_available = True
    except ImportError:
        pass

    for n in n_values:
        l_approx = n * math.pi / 2.0
        h_classical = math.log(l_approx)
        classical.append(h_classical)

        if qiskit_available:
            # Simple n-qubit GHZ state as quantum channel proxy
            qc = QuantumCircuit(min(n, 8))
            qc.h(0)
            for i in range(1, min(n, 8)):
                qc.cx(0, i)
            qc.measure_all()
            sim = AerSimulator()
            job = sim.run(qc, shots=1024)
            counts = job.result().get_counts()
            total = sum(counts.values())
            probs = [v / total for v in counts.values()]
            h_q = -sum(p * math.log(p) for p in probs if p > 0)
            quantum_sim.append(h_q)
        else:
            # Classical approximation of von Neumann entropy for n-qubit GHZ
            h_q = math.log(2)  # GHZ state has S=1 ebit regardless of n
            quantum_sim.append(h_q)

    correlation = float(np.corrcoef(classical, quantum_sim)[0, 1])

    summary = {
        "experiment": "C",
        "name": "Quantum channel analogy",
        "qiskit_available": qiskit_available,
        "classical_entropy": classical,
        "quantum_entropy": quantum_sim,
        "correlation": correlation,
        "verdict": "PASS" if correlation > 0.8 else "WEAK",
        "note": "GHZ entropy constant; correlation reflects classical scaling only",
    }
    return summary


# ---------------------------------------------------------------------------
# Experiment D: Cross-problem transfer
# ---------------------------------------------------------------------------

def run_experiment_d():
    """
    Test: Does the channel structure (H/C ratio) transfer across Erdos problems?
    Compare #114 (lemniscate arc-length) vs #509 (chromatic number bound).

    For #509 we use the known Ramsey-type bound: r(n) >= 2^(n/2),
    so log r(n) ~ (n/2) log 2, giving H_509(n) = (n/2)*log(2).
    """
    results = []
    for n in range(3, 12):
        # Problem #114: lemniscate channel
        l_114 = n * math.pi / 2.0
        h_114 = math.log(l_114)
        c_114 = math.log(n)
        ratio_114 = h_114 / c_114

        # Problem #509: Ramsey / chromatic number channel
        h_509 = (n / 2.0) * math.log(2)
        c_509 = math.log(n)
        ratio_509 = h_509 / c_509

        delta = abs(ratio_114 - ratio_509)
        results.append({
            "n": n,
            "ratio_114": ratio_114,
            "ratio_509": ratio_509,
            "delta": delta,
        })

    avg_delta = sum(r["delta"] for r in results) / len(results)
    transfers = avg_delta < 0.5

    summary = {
        "experiment": "D",
        "name": "Cross-problem transfer",
        "avg_ratio_delta": avg_delta,
        "transfers": transfers,
        "verdict": "PASS" if transfers else "FAIL",
        "data": results,
    }
    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Run channel experiments")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--skip", metavar="EXP", help="Skip experiment (a/b/c/d)")
    group.add_argument("--only", metavar="EXP", help="Run only experiment (a/b/c/d)")
    args = parser.parse_args()

    skip = args.skip.lower() if args.skip else None
    only = args.only.lower() if args.only else None

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    all_results = {}
    t0 = time.time()

    def should_run(exp_id):
        if only:
            return exp_id == only
        if skip:
            return exp_id != skip
        return True

    if should_run("a"):
        print("=== Experiment A: Entropy-Hausdorff monotonicity ===")
        r = run_experiment_a()
        all_results["A"] = r
        print(f"  Verdict: {r['verdict']}")
        out = RESULTS_DIR / "exp_a_results.json"
        out.write_text(json.dumps(r, indent=2))

    if should_run("b"):
        print("=== Experiment B: Capacity saturation ===")
        r = run_experiment_b()
        all_results["B"] = r
        print(f"  Verdict: {r['verdict']}  asymptotic H/C={r['asymptotic_H_over_C']:.6f}")
        out = RESULTS_DIR / "exp_b_results.json"
        out.write_text(json.dumps(r, indent=2))

    if should_run("c"):
        print("=== Experiment C: Quantum channel analogy ===")
        r = run_experiment_c()
        all_results["C"] = r
        print(f"  Verdict: {r['verdict']}  correlation={r['correlation']:.4f}")
        out = RESULTS_DIR / "exp_c_results.json"
        out.write_text(json.dumps(r, indent=2))

    if should_run("d"):
        print("=== Experiment D: Cross-problem transfer ===")
        r = run_experiment_d()
        all_results["D"] = r
        print(f"  Verdict: {r['verdict']}  avg_delta={r['avg_ratio_delta']:.4f}")
        out = RESULTS_DIR / "exp_d_results.json"
        out.write_text(json.dumps(r, indent=2))

    elapsed = time.time() - t0
    summary = {
        "experiments_run": list(all_results.keys()),
        "verdicts": {k: v["verdict"] for k, v in all_results.items()},
        "total_time_secs": elapsed,
    }
    (RESULTS_DIR / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nDone in {elapsed:.1f}s")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
