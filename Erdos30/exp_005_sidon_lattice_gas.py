#!/usr/bin/env python3
"""
EXP-MATH-ERDOS30-SIDON-005: Sidon Lattice Gas — The Collider Model

Erdős Collider framing: Sidon sets as molecular assemblies.
  - ATOMS: prime numbers (the irreducible building blocks)
  - MOLECULES: integers in {1,...,N}, each with a prime factorization
  - MATERIAL: a Sidon set A ⊂ {1,...,N} — a collection of molecules
  - BOND CONSTRAINT: all pairwise sums distinct (hard-core exclusion in difference space)
  - CLOSE-PACKING DENSITY: h(N)/√N → 1 (the Erdős conjecture)

This is a lattice gas on Z_n with a specific exclusion rule:
  site i is "occupied" iff i ∈ A
  constraint: for all pairs (i,j) and (k,l) in A×A, i+j = k+l → {i,j} = {k,l}

Materials science observables:
  1. Formation energy: cost of adding/removing one element
  2. Pair potential: interaction energy between elements via shared differences
  3. Phonon spectrum: Fourier transform of the occupation function (= EXP-002!)
  4. Chemical potential: prime factorization profile of elements near the surface
  5. Packing fraction: |A|/h_max(N) — how close to close-packing
  6. Defect formation energy: cost of the most favorable single-site vacancy/insertion
  7. Phase diagram: algebraic (Singer) vs amorphous (greedy) vs liquid (random)

The key insight: the decorated cospan A →^i M ←^o B where
  A = prime-factored elements ("reagents")
  M = the Sidon constraint ("reaction vessel" / encoding mechanism)
  B = the Fourier structure ("product")
  d = MDL cost of the encoding
"""

import json
import math
import time
import random
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np

# ─── Singer PDS (verified, q ≤ 13) ──────────────────────────────────────

SINGER_SETS = {
    2: [0, 1, 3],
    3: [0, 1, 3, 9],
    4: [0, 1, 4, 14, 16],
    5: [0, 1, 3, 8, 12, 18],
    7: [0, 1, 3, 13, 32, 36, 43, 52],
    8: [0, 1, 3, 7, 15, 31, 36, 54, 63],
    9: [0, 1, 3, 9, 27, 49, 56, 61, 77, 81],
    11: [0, 1, 3, 12, 20, 34, 38, 81, 88, 94, 104, 109],
    13: [0, 1, 3, 16, 23, 28, 42, 76, 82, 86, 119, 137, 154, 175],
}


def greedy_sidon(n, target_size):
    A = [0]
    diff_set = set()
    for x in range(1, n):
        diffs = set()
        valid = True
        for a in A:
            d1 = (x - a) % n
            d2 = (a - x) % n
            if d1 in diff_set or d1 in diffs or d2 in diff_set or d2 in diffs:
                valid = False
                break
            diffs.add(d1)
            diffs.add(d2)
        if valid:
            A.append(x)
            diff_set.update(diffs)
            if len(A) >= target_size:
                break
    return A


def random_sidon(n, target_size, seed=42):
    rng = random.Random(seed)
    candidates = list(range(n))
    rng.shuffle(candidates)
    A = []
    diff_set = set()
    for x in candidates:
        diffs = set()
        valid = True
        for a in A:
            d1 = (x - a) % n
            d2 = (a - x) % n
            if d1 in diff_set or d1 in diffs or d2 in diff_set or d2 in diffs:
                valid = False
                break
            diffs.add(d1)
            diffs.add(d2)
        if valid:
            A.append(x)
            diff_set.update(diffs)
            if len(A) >= target_size:
                break
    return sorted(A)


# ─── Prime Factorization (Atomic Structure) ──────────────────────────────

def prime_factorize(n):
    """Return prime factorization as {prime: exponent}."""
    if n <= 1:
        return {}
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def atomic_number(n):
    """Number of prime factors with multiplicity (Ω function)."""
    if n <= 1:
        return 0
    return sum(prime_factorize(n).values())


def atomic_weight(n):
    """Sum of prime factors with multiplicity."""
    if n <= 1:
        return 0
    return sum(p * e for p, e in prime_factorize(n).items())


def atomic_diversity(n):
    """Number of distinct prime factors (ω function)."""
    if n <= 1:
        return 0
    return len(prime_factorize(n))


def molecular_descriptor(n):
    """Full molecular descriptor for element n."""
    factors = prime_factorize(n)
    return {
        "value": n,
        "is_prime": len(factors) == 1 and list(factors.values())[0] == 1 if factors else False,
        "is_prime_power": len(factors) == 1 if factors else False,
        "omega": len(factors),                          # distinct prime factors
        "Omega": sum(factors.values()) if factors else 0,  # total prime factors w/ multiplicity
        "smallest_prime": min(factors.keys()) if factors else 0,
        "largest_prime": max(factors.keys()) if factors else 0,
        "radical": math.prod(factors.keys()) if factors else 1,
        "factorization": factors,
    }


# ─── Lattice Gas Thermodynamics ──────────────────────────────────────────

def difference_set(A, n):
    """Compute all pairwise differences mod n."""
    diffs = set()
    for i, a in enumerate(A):
        for b in A[i+1:]:
            diffs.add((a - b) % n)
            diffs.add((b - a) % n)
    return diffs


def formation_energy(A, n, site, occupied=True):
    """
    Energy cost of adding (occupied=True) or removing (occupied=False)
    a single site from the Sidon set.

    Energy = number of difference collisions created (add) or released (remove).
    For a valid Sidon set, adding any occupied element had energy 0 (no collisions).
    The "surface energy" is the minimum insertion energy over all unoccupied sites.
    """
    A_set = set(A)
    existing_diffs = difference_set(A, n)

    if occupied:
        # Cost of adding site to A
        if site in A_set:
            return 0  # already present
        new_diffs = set()
        collisions = 0
        for a in A:
            d1 = (site - a) % n
            d2 = (a - site) % n
            if d1 in existing_diffs:
                collisions += 1
            if d2 in existing_diffs:
                collisions += 1
            new_diffs.add(d1)
            new_diffs.add(d2)
        return collisions
    else:
        # Cost of removing site from A (= energy released)
        if site not in A_set:
            return 0
        A_without = [a for a in A if a != site]
        diffs_without = difference_set(A_without, n)
        released = len(existing_diffs) - len(diffs_without)
        return released


def surface_energy(A, n):
    """
    Minimum insertion energy over all unoccupied sites.
    If 0: the Sidon set is not maximal (can still grow).
    If > 0: the set is maximal — every insertion creates collisions.
    """
    A_set = set(A)
    existing_diffs = difference_set(A, n)
    min_energy = float('inf')
    min_site = -1

    # Sample unoccupied sites (for large n, sample rather than exhaustive)
    unoccupied = [x for x in range(n) if x not in A_set]
    if len(unoccupied) > 200:
        rng = random.Random(42)
        unoccupied = rng.sample(unoccupied, 200)

    for site in unoccupied:
        e = formation_energy(A, n, site, occupied=True)
        if e < min_energy:
            min_energy = e
            min_site = site

    return min_energy, min_site


def chemical_potential_profile(A, n):
    """
    Analyze the prime-factorization "chemistry" of the Sidon set elements.
    Returns aggregate molecular descriptors.
    """
    nonzero_elements = [a for a in A if a > 1]
    if not nonzero_elements:
        return {"note": "all elements ≤ 1"}

    descriptors = [molecular_descriptor(a) for a in nonzero_elements]

    primes_in_set = [d for d in descriptors if d["is_prime"]]
    prime_powers = [d for d in descriptors if d["is_prime_power"]]

    omegas = [d["omega"] for d in descriptors]
    Omegas = [d["Omega"] for d in descriptors]
    smallest_primes = [d["smallest_prime"] for d in descriptors if d["smallest_prime"] > 0]

    # "Chemical formula" — aggregate prime content
    total_prime_content = Counter()
    for d in descriptors:
        for p, e in d["factorization"].items():
            total_prime_content[p] += e

    return {
        "n_elements": len(nonzero_elements),
        "n_primes": len(primes_in_set),
        "prime_fraction": len(primes_in_set) / len(nonzero_elements),
        "n_prime_powers": len(prime_powers),
        "prime_power_fraction": len(prime_powers) / len(nonzero_elements),
        "mean_omega": float(np.mean(omegas)),
        "mean_Omega": float(np.mean(Omegas)),
        "mean_smallest_prime": float(np.mean(smallest_primes)) if smallest_primes else 0,
        "median_smallest_prime": float(np.median(smallest_primes)) if smallest_primes else 0,
        "top_5_primes": dict(total_prime_content.most_common(5)),
        "total_atomic_number": sum(Omegas),  # total "atoms" in the material
    }


def packing_fraction(A, n):
    """How close to theoretical close-packing. For Z_n, h_max ≈ √n."""
    return len(A) / math.sqrt(n)


def vacancy_spectrum(A, n, max_vacancies=50):
    """
    Compute the energy cost of each possible single-site vacancy.
    Low-cost vacancies = structural defects; high-cost = load-bearing elements.
    """
    energies = []
    for a in A:
        if a == 0:
            continue
        e = formation_energy(A, n, a, occupied=False)
        desc = molecular_descriptor(a)
        energies.append({
            "site": a,
            "vacancy_energy": e,
            "is_prime": desc["is_prime"],
            "omega": desc["omega"],
        })

    energies.sort(key=lambda x: x["vacancy_energy"])
    return energies[:max_vacancies]


def pair_correlation(A, n):
    """
    Pair correlation function g(r) — probability of finding an element
    at distance r from another element, normalized by uniform density.
    This is the materials science standard for characterizing short-range order.
    """
    k = len(A)
    density = k / n

    # Count pairs at each distance
    dist_counts = Counter()
    for i, a in enumerate(A):
        for b in A[i+1:]:
            d = min((a - b) % n, (b - a) % n)
            dist_counts[d] += 1

    # g(r) = (count at r) / (expected count at r for random placement)
    # Expected: k*(k-1)/2 * (1 if d appears) / n ≈ density * k/2
    g_r = {}
    expected_per_distance = k * (k - 1) / (2 * (n // 2))
    for r in range(1, n // 2 + 1):
        g_r[r] = dist_counts.get(r, 0) / expected_per_distance if expected_per_distance > 0 else 0

    # For a Sidon set: each distance appears AT MOST once → g(r) ∈ {0, 1/expected}
    # Sidon constraint = hard-core exclusion at distance 0 in difference-pair space

    g_values = list(g_r.values())
    return {
        "g_mean": float(np.mean(g_values)),
        "g_std": float(np.std(g_values)),
        "g_max": float(np.max(g_values)),
        "g_cv": float(np.std(g_values) / np.mean(g_values)) if np.mean(g_values) > 0 else 0,
        "occupied_distances": len(dist_counts),
        "total_distances": n // 2,
        "distance_filling": len(dist_counts) / (n // 2),
    }


# ─── Phase Diagram: Temperature Scan ─────────────────────────────────────

def thermal_occupation(A, n, temperature):
    """
    At temperature T, allow "defects" with probability proportional to
    exp(-E/T) where E is the number of difference collisions.

    Returns the expected packing fraction at temperature T.
    At T=0: only collision-free (Sidon) configurations survive.
    At T→∞: random occupation, packing fraction → density.
    """
    A_set = set(A)
    existing_diffs = difference_set(A, n)

    # For each unoccupied site, compute Boltzmann weight
    total_weight = 0
    expected_additions = 0

    unoccupied = [x for x in range(n) if x not in A_set]
    if len(unoccupied) > 300:
        rng = random.Random(42)
        unoccupied = rng.sample(unoccupied, 300)
        scale = (n - len(A)) / 300
    else:
        scale = 1.0

    for site in unoccupied:
        e = formation_energy(A, n, site, occupied=True)
        if temperature > 0:
            boltzmann = math.exp(-e / temperature)
        else:
            boltzmann = 1.0 if e == 0 else 0.0
        expected_additions += boltzmann

    expected_additions *= scale
    return (len(A) + expected_additions) / n


# ─── Main Experiment ─────────────────────────────────────────────────────

def run_experiment():
    results = {
        "experiment_id": "EXP-MATH-ERDOS30-SIDON-005",
        "title": "Sidon Lattice Gas — Collider Model (Sets as Molecules, Primes as Atoms)",
        "collider_cospan": "A(prime-factored elements) →[i] M(Sidon constraint) ←[o] B(Fourier structure)",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "singer": [],
        "greedy": [],
        "random": [],
        "phase_diagram": [],
        "analysis": {},
    }

    print("=" * 90)
    print("EXP-MATH-ERDOS30-SIDON-005: Sidon Lattice Gas (Collider Model)")
    print("Sets as molecules, primes as atoms, Sidon constraint as bond rule")
    print("=" * 90)

    # ── Per-construction analysis ──
    for label, get_set in [
        ("Singer", lambda q, n, k: SINGER_SETS[q]),
        ("Greedy", lambda q, n, k: greedy_sidon(n, k)),
        ("Random", lambda q, n, k: random_sidon(n, k, seed=42)),
    ]:
        print(f"\n{'─' * 40}")
        print(f"  {label} Construction")
        print(f"{'─' * 40}")

        construction_results = []
        for q in [3, 5, 7, 11, 13]:
            n = q * q + q + 1
            k = q + 1
            A = get_set(q, n, k)
            actual_k = len(A)

            print(f"\n  q={q}, n={n}, |A|={actual_k}")

            # 1. Packing fraction
            pf = packing_fraction(A, n)
            print(f"    Packing fraction: {pf:.4f}")

            # 2. Chemical potential profile
            chem = chemical_potential_profile(A, n)
            print(f"    Prime fraction: {chem['prime_fraction']:.3f}, "
                  f"Mean ω: {chem['mean_omega']:.2f}, "
                  f"Mean Ω: {chem['mean_Omega']:.2f}")

            # 3. Surface energy
            surf_e, surf_site = surface_energy(A, n)
            print(f"    Surface energy (min insertion cost): {surf_e}")
            if surf_e > 0:
                print(f"      → MAXIMAL: every insertion creates ≥{surf_e} collisions")
            else:
                print(f"      → NOT maximal: site {surf_site} can be added")

            # 4. Pair correlation
            pcorr = pair_correlation(A, n)
            print(f"    Pair correlation g(r): mean={pcorr['g_mean']:.3f}, "
                  f"CV={pcorr['g_cv']:.3f}, "
                  f"distance filling={pcorr['distance_filling']:.3f}")

            # 5. Vacancy spectrum (top 3 easiest vacancies)
            vacancies = vacancy_spectrum(A, n, max_vacancies=5)
            if vacancies:
                print(f"    Vacancy spectrum (easiest removals):")
                for v in vacancies[:3]:
                    print(f"      site {v['site']}: energy={v['vacancy_energy']}, "
                          f"prime={v['is_prime']}, ω={v['omega']}")

            record = {
                "q": q, "n": n, "k": actual_k,
                "packing_fraction": pf,
                "surface_energy": surf_e,
                "chemistry": chem,
                "pair_correlation": pcorr,
                "vacancy_top3": vacancies[:3] if vacancies else [],
            }
            construction_results.append(record)

        if label == "Singer":
            results["singer"] = construction_results
        elif label == "Greedy":
            results["greedy"] = construction_results
        else:
            results["random"] = construction_results

    # ── Phase Diagram: Temperature Scan ──
    print(f"\n{'=' * 90}")
    print("PHASE DIAGRAM: Temperature Scan")
    print("=" * 90)
    print("At T=0: only Sidon configurations. As T→∞: random occupation.")
    print("Phase transition = the temperature at which the Sidon constraint 'melts'.\n")

    for q in [5, 7, 11]:
        n = q * q + q + 1
        k = q + 1
        A_singer = SINGER_SETS[q]
        A_greedy = greedy_sidon(n, k)

        temps = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
        print(f"  q={q}, n={n}:")
        print(f"    {'T':>6s} | {'Singer ρ':>10s} | {'Greedy ρ':>10s}")
        print(f"    {'-'*6}-+-{'-'*10}-+-{'-'*10}")

        phase_q = {"q": q, "n": n, "temperatures": [], "singer_density": [], "greedy_density": []}
        for T in temps:
            rho_s = thermal_occupation(A_singer, n, T)
            rho_g = thermal_occupation(A_greedy, n, T)
            print(f"    {T:6.1f} | {rho_s:10.4f} | {rho_g:10.4f}")
            phase_q["temperatures"].append(T)
            phase_q["singer_density"].append(float(rho_s))
            phase_q["greedy_density"].append(float(rho_g))

        results["phase_diagram"].append(phase_q)

    # ── Synthesis ──
    print(f"\n{'=' * 90}")
    print("COLLIDER SYNTHESIS")
    print("=" * 90)

    # Compare Singer vs Greedy across all q
    print("\nComparative Table:")
    print(f"  {'q':>4s} | {'Type':>7s} | {'Pack':>6s} | {'Surf.E':>6s} | {'%Prime':>7s} | {'g(r)CV':>7s} | {'Dist.Fill':>9s}")
    print(f"  {'-'*4}-+-{'-'*7}-+-{'-'*6}-+-{'-'*6}-+-{'-'*7}-+-{'-'*7}-+-{'-'*9}")

    for i, q in enumerate([3, 5, 7, 11, 13]):
        for label, data in [("Singer", results["singer"]), ("Greedy", results["greedy"]), ("Random", results["random"])]:
            r = data[i]
            print(f"  {q:4d} | {label:>7s} | {r['packing_fraction']:6.3f} | "
                  f"{r['surface_energy']:6d} | {r['chemistry']['prime_fraction']:7.3f} | "
                  f"{r['pair_correlation']['g_cv']:7.3f} | "
                  f"{r['pair_correlation']['distance_filling']:9.3f}")

    # Key findings
    singer_pf = [r["packing_fraction"] for r in results["singer"]]
    greedy_pf = [r["packing_fraction"] for r in results["greedy"]]
    singer_surf = [r["surface_energy"] for r in results["singer"]]
    greedy_surf = [r["surface_energy"] for r in results["greedy"]]
    singer_prime = [r["chemistry"]["prime_fraction"] for r in results["singer"]]
    greedy_prime = [r["chemistry"]["prime_fraction"] for r in results["greedy"]]
    singer_gcv = [r["pair_correlation"]["g_cv"] for r in results["singer"]]
    greedy_gcv = [r["pair_correlation"]["g_cv"] for r in results["greedy"]]

    print(f"\n  Singer packing: {[f'{p:.3f}' for p in singer_pf]}")
    print(f"  Greedy packing: {[f'{p:.3f}' for p in greedy_pf]}")
    print(f"\n  Singer surface energy: {singer_surf}")
    print(f"  Greedy surface energy: {greedy_surf}")
    print(f"\n  Singer prime fraction: {[f'{p:.3f}' for p in singer_prime]}")
    print(f"  Greedy prime fraction: {[f'{p:.3f}' for p in greedy_prime]}")
    print(f"\n  Singer g(r) CV: {[f'{c:.3f}' for c in singer_gcv]}")
    print(f"  Greedy g(r) CV: {[f'{c:.3f}' for c in greedy_gcv]}")

    # Verdict
    pf_better = sum(1 for s, g in zip(singer_pf, greedy_pf) if s > g)
    gcv_lower = sum(1 for s, g in zip(singer_gcv, greedy_gcv) if s < g)

    print(f"\n  Singer packing > Greedy: {pf_better}/{len(singer_pf)}")
    print(f"  Singer g(r) CV < Greedy (more uniform): {gcv_lower}/{len(singer_gcv)}")

    results["analysis"] = {
        "singer_packing_fractions": singer_pf,
        "greedy_packing_fractions": greedy_pf,
        "singer_surface_energies": singer_surf,
        "greedy_surface_energies": greedy_surf,
        "singer_prime_fractions": singer_prime,
        "greedy_prime_fractions": greedy_prime,
        "singer_pair_corr_cv": singer_gcv,
        "greedy_pair_corr_cv": greedy_gcv,
        "singer_packing_better_count": pf_better,
        "singer_gcv_lower_count": gcv_lower,
    }

    # Save
    out_path = Path(__file__).parent / "EXP-MATH-ERDOS30-SIDON-005_RESULTS.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")

    return results


if __name__ == "__main__":
    run_experiment()
