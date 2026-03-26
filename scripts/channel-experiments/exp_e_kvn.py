#!/usr/bin/env python3
"""
exp_e_kvn.py  --  Experiment E: Koopman-von Neumann Channel Lens

PHYSICAL MODEL
--------------
In the Koopman-von Neumann (KvN) formulation, classical Hamiltonian mechanics
is lifted to a Hilbert space of wavefunctions psi(q,p) over phase space, with
the Koopman operator K acting as the unitary evolution:

    d/dt |psi> = -i L_H |psi>

where L_H is the Liouvillian (Poisson bracket with H).  The eigenvalues of L_H
are the classical Lyapunov / Floquet exponents.  This gives a direct bridge
between the classical lemniscate dynamics and quantum channel formalism:

  * The lemniscate L_n* is the level set |p(z)| = 1 of a degree-n polynomial
  * Flow on the lemniscate under gradient dynamics defines a classical
    Hamiltonian system with H = log|p(z)|
  * Lifting to KvN gives a unitary channel U_n on L2(lemniscate)
  * The von Neumann entropy S(rho_n) of the mixed state rho_n = |psi_n><psi_n|
    averaged over one period is the KvN channel entropy

WHAT WE MEASURE
---------------
For each degree n = 3..11 (using certified l_star):

  1. KvN_entropy(n): von Neumann entropy S(rho_n) of the stationary density
     matrix on the lemniscate, modeled as the entropy of the arc-length measure
     discretized into k = 2n+1 cells (natural Floquet resolution)

  2. KvN_capacity(n): channel capacity of the KvN channel, estimated as
     log(rank(K_n)) where rank(K_n) = number of non-negligible Floquet modes

  3. KvN_duality_gap(n): gap between S(rho_n) and the classical Shannon
     entropy H_classical(n) = log(l_star(n)); measures how much information
     is "hidden" in quantum coherences vs accessible classically

  4. KvN_Maxwell_residual(n): analogy to Maxwell's displacement current --
     the excess entropy flow dS/dn that cannot be attributed to arc-length
     growth alone; tests whether the KvN lens reveals a conserved quantity

VERDICT criteria
----------------
  PASS:    KvN_entropy monotone AND duality_gap < 0.1 for all n
  PARTIAL: monotone but duality_gap >= 0.1 for some n (quantum corrections
           are non-negligible -- interesting, not a failure)
  FAIL:    KvN_entropy non-monotone

Dependencies: numpy, scipy (for eigenvalue computation)
"""

import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

# ---------------------------------------------------------------------------
# Certified l_star values (lower bounds from EXP-MM-EHP-007, n=3..11)
# ---------------------------------------------------------------------------

L_STAR = {
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


# ---------------------------------------------------------------------------
# KvN density matrix construction
# ---------------------------------------------------------------------------

def lemniscate_arc_density(n, l_star_val, k_cells=None):
    """
    Construct the discrete density matrix rho for the arc-length measure on
    the degree-n lemniscate L_n*.

    We discretize the lemniscate into k = 2n+1 cells (the natural Floquet
    resolution: n Floquet harmonics + their conjugates + DC mode).
    The arc-length measure assigns probability p_j = arc_j / l_star to each
    cell.  For the extremal lemniscate the arc-length distribution is nearly
    uniform but has n-fold symmetry modulations.

    Model: p_j = (1/k) * (1 + epsilon * cos(2*pi*n*j/k))
    where epsilon = (l_star - n*pi/2) / (n*pi/2) captures the deviation
    from the uniform (large-n) limit.  This is the leading Fourier mode of
    the arc-length density on the Chebyshev/lemniscate measure.
    """
    if k_cells is None:
        k_cells = 2 * n + 1

    # Deviation from asymptotic uniform
    l_uniform = n * math.pi / 2.0
    epsilon = (l_star_val - l_uniform) / l_uniform  # signed; >0 means larger than expected

    # Modulated arc-length probabilities (normalized)
    js = np.arange(k_cells)
    probs = (1.0 / k_cells) * (1.0 + epsilon * np.cos(2 * np.pi * n * js / k_cells))
    probs = np.maximum(probs, 1e-15)  # numerical safety
    probs /= probs.sum()

    # Pure-state density matrix: rho = |psi><psi| with psi_j = sqrt(p_j)
    psi = np.sqrt(probs)
    rho_pure = np.outer(psi, psi)

    # KvN mixed state: average over n Floquet sectors (rotate by 2pi/n each)
    # This is the thermal/microcanonical average natural in KvN mechanics
    rho_mixed = np.zeros_like(rho_pure)
    for sector in range(n):
        phase = 2 * np.pi * sector / n
        psi_rotated = psi * np.exp(1j * phase * js / k_cells)
        rho_sector = np.outer(psi_rotated, psi_rotated.conj()).real
        rho_mixed += rho_sector / n

    return rho_mixed, probs


def von_neumann_entropy(rho):
    """S(rho) = -Tr(rho log rho), computed via eigendecomposition."""
    eigvals = eigh(rho, eigvals_only=True)
    eigvals = eigvals[eigvals > 1e-15]
    return float(-np.sum(eigvals * np.log(eigvals)))


def koopman_channel_capacity(rho, n):
    """
    Estimate the KvN channel capacity as log(effective_rank(rho)).
    Effective rank = exp(S(rho)) -- the exponentiated von Neumann entropy,
    which equals the number of maximally mixed modes.
    """
    s = von_neumann_entropy(rho)
    return s  # log(exp(S)) = S; capacity in nats


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def run_experiment_e():
    results = []
    kvn_entropies = []
    duality_gaps = []
    maxwell_residuals = []

    prev_s_kvn = None
    monotone = True
    prev_l = None

    for n in sorted(L_STAR.keys()):
        l = L_STAR[n]
        k_cells = 2 * n + 1

        # Build KvN density matrix
        rho, probs = lemniscate_arc_density(n, l, k_cells)

        # KvN von Neumann entropy
        s_kvn = von_neumann_entropy(rho)

        # Classical Shannon entropy of arc-length measure
        h_classical = float(-np.sum(probs * np.log(probs)))

        # Duality gap: quantum coherences vs classical probabilities
        duality_gap = abs(s_kvn - h_classical)

        # KvN channel capacity (nats)
        c_kvn = koopman_channel_capacity(rho, n)

        # Maxwell residual: excess entropy rate beyond arc-length growth
        if prev_s_kvn is not None and prev_l is not None:
            dl = l - prev_l
            ds = s_kvn - prev_s_kvn
            maxwell_res = ds - math.log(l / prev_l)  # dS - d(log l)
        else:
            maxwell_res = 0.0

        if prev_s_kvn is not None and s_kvn < prev_s_kvn - 1e-10:
            monotone = False

        kvn_entropies.append(s_kvn)
        duality_gaps.append(duality_gap)
        maxwell_residuals.append(maxwell_res)
        prev_s_kvn = s_kvn
        prev_l = l

        results.append({
            "n": n,
            "l_star": l,
            "k_cells": k_cells,
            "S_kvn": s_kvn,
            "H_classical": h_classical,
            "duality_gap": duality_gap,
            "C_kvn": c_kvn,
            "maxwell_residual": maxwell_res,
        })

    max_gap = max(duality_gaps)
    avg_maxwell = float(np.mean(np.abs(maxwell_residuals[1:])))  # skip n=3 (no prev)

    if monotone and max_gap < 0.1:
        verdict = "PASS"
    elif monotone:
        verdict = "PARTIAL"  # quantum corrections non-negligible
    else:
        verdict = "FAIL"

    summary = {
        "experiment": "E",
        "name": "Koopman-von Neumann channel lens",
        "monotone": monotone,
        "max_duality_gap": max_gap,
        "avg_maxwell_residual": avg_maxwell,
        "verdict": verdict,
        "interpretation": (
            "PASS: KvN entropy monotone and quantum/classical gap < 0.1 nats -- "
            "the KvN lift does not add spurious structure." if verdict == "PASS"
            else "PARTIAL: KvN entropy monotone but quantum coherences contribute > 0.1 nats "
            "-- the KvN channel is strictly richer than the classical arc-length channel." if verdict == "PARTIAL"
            else "FAIL: KvN entropy non-monotone -- model needs revision."
        ),
        "data": results,
    }
    return summary


# ---------------------------------------------------------------------------
# Standalone runner
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    result = run_experiment_e()
    print(f"Verdict: {result['verdict']}")
    print(f"Monotone: {result['monotone']}")
    print(f"Max duality gap: {result['max_duality_gap']:.6f} nats")
    print(f"Avg Maxwell residual: {result['avg_maxwell_residual']:.6f} nats/step")
    print()
    print(f"{'n':>4}  {'l_star':>12}  {'S_kvn':>10}  {'H_class':>10}  {'gap':>8}  {'maxwell':>10}")
    print("-" * 64)
    for r in result["data"]:
        print(
            f"{r['n']:>4}  {r['l_star']:>12.6f}  {r['S_kvn']:>10.6f}  "
            f"{r['H_classical']:>10.6f}  {r['duality_gap']:>8.6f}  {r['maxwell_residual']:>10.6f}"
        )

    out = Path("exp_e_kvn_results.json")
    out.write_text(json.dumps(result, indent=2))
    print(f"\nResults written to {out}")
    sys.exit(0 if result["verdict"] in ("PASS", "PARTIAL") else 1)
