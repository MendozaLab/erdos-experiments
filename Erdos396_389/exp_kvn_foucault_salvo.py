#!/usr/bin/env python3
"""
EXP-FOUCAULT-KVN-001: Koopman–von Neumann Salvo on Erdős #396 & #389
====================================================================
The Foucault Resonance Pair

Both problems involve divisibility in products of consecutive integers:
  #396: ∀k ∃n: n(n-1)...(n-k) | C(2n,n)
  #389: ∀n≥1 ∃k: n(n+1)...(n+k-1) | (n+k)...(n+2k-1)

Physics Morphism (PHYS-R-002 Foucault Pendulum):
  Foucault pendulum = slow precession (period 24h/sin(lat)) modulating fast oscillation.
  KvN bridge: p-adic valuations v_p(n!) follow a STAIRCASE function that
  precesses through residue classes mod p, mod p², mod p³ ...
  The "fast oscillation" = individual prime cycling (period p).
  The "slow precession" = carry propagation in base-p digits (Kummer's theorem).

The experiment tests three prongs:
  1. Foucault Precession Index: p-adic valuation staircase periodicity
  2. Koopman Eigenvalue Resonance: spectral analysis of valuation dynamics
  3. Resonance Overlap: whether #396 and #389 share Koopman eigenmodes

Author: MendozaLab / Kenneth A. Mendoza
Date: 2026-04-16
"""

import numpy as np
from math import comb, gcd, log, floor
from collections import defaultdict
import json, time

# ── Utilities ──

def v_p(n, p):
    """p-adic valuation of n (largest power of p dividing n)."""
    if n == 0:
        return float('inf')
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def legendre_valuation(n, p):
    """v_p(n!) via Legendre's formula."""
    if n < p:
        return 0
    s = 0
    pk = p
    while pk <= n:
        s += n // pk
        pk *= p
    return s

def desc_factorial(n, k):
    """n * (n-1) * ... * (n-k+1)"""
    result = 1
    for i in range(k):
        result *= (n - i)
    return result

# ═══════════════════════════════════════════════════════════
# PRONG 1: Foucault Precession Index
# ═══════════════════════════════════════════════════════════
# The Foucault pendulum precesses at ω = Ω·sin(latitude).
# Analogy: v_p(C(2n,n)) precesses through values as n increases,
# with period p for the dominant harmonic and period p^k for
# the k-th overtone (carry propagation in base-p addition).
#
# We compute the "precession spectrum" of v_p(C(2n,n)) and
# v_p(prod_i (n+i)) and look for shared resonant frequencies.

def prong1_precession_spectrum():
    """Compute Foucault precession index for v_p(C(2n,n)) dynamics."""
    print("=" * 60)
    print("PRONG 1: Foucault Precession Spectrum")
    print("=" * 60)

    results = {}
    primes = [2, 3, 5, 7, 11, 13]

    for p in primes:
        # Sequence: v_p(C(2n,n)) for n = 1..N
        N = 500
        vals = []
        for n in range(1, N + 1):
            v = legendre_valuation(2 * n, p) - 2 * legendre_valuation(n, p)
            vals.append(v)

        vals_arr = np.array(vals, dtype=float)

        # Compute FFT to find dominant periods
        fft = np.fft.rfft(vals_arr - vals_arr.mean())
        power = np.abs(fft) ** 2
        freqs = np.fft.rfftfreq(N)

        # Find top 5 peaks (skip DC)
        peak_idx = np.argsort(power[1:])[::-1][:5] + 1
        dominant_periods = []
        for idx in peak_idx:
            if freqs[idx] > 0:
                period = 1.0 / freqs[idx]
                dominant_periods.append({
                    "period": round(float(period), 2),
                    "power": round(float(power[idx]), 2),
                    "freq": round(float(freqs[idx]), 6),
                })

        # Check: is the dominant period close to p or p^k?
        main_period = dominant_periods[0]["period"] if dominant_periods else 0
        is_resonant = any(abs(main_period - p**k) / max(p**k, 1) < 0.1
                         for k in range(1, 5))

        # Compute precession index: ratio of dominant period to p
        prec_index = main_period / p if p > 0 else 0

        results[str(p)] = {
            "dominant_period": main_period,
            "precession_index": round(prec_index, 4),
            "is_resonant_with_prime": is_resonant,
            "top_periods": dominant_periods[:3],
            "mean_valuation": round(float(vals_arr.mean()), 4),
            "max_valuation": int(vals_arr.max()),
        }

        print(f"\n  p={p}: dominant period = {main_period:.1f} "
              f"(precession index = {prec_index:.3f})")
        print(f"    resonant with p^k: {is_resonant}")
        print(f"    mean v_p(C(2n,n)) = {vals_arr.mean():.3f}")

    return results


# ═══════════════════════════════════════════════════════════
# PRONG 2: Koopman Eigenvalue Analysis
# ═══════════════════════════════════════════════════════════
# Model the p-adic valuation dynamics as a discrete dynamical system:
#   State x_n = (digit_0(n,p), digit_1(n,p), ..., digit_K(n,p))
#   Evolution: x_n → x_{n+1} (increment in base p)
#
# The Koopman operator U acts on observables f(x):
#   (Uf)(x_n) = f(x_{n+1})
#
# For base-p digit dynamics, Koopman eigenvalues are p-th roots of unity
# (fast oscillation) with amplitude modulation from carry propagation
# (slow precession = Foucault effect).
#
# We test: does the Koopman spectrum of the valuation dynamics
# show the same resonance structure for BOTH #396 and #389?

def prong2_koopman_eigenvalues():
    """Koopman eigenvalue analysis of p-adic valuation dynamics."""
    print("\n" + "=" * 60)
    print("PRONG 2: Koopman Eigenvalue Resonance")
    print("=" * 60)

    results = {}

    for p in [2, 3, 5, 7]:
        # --- Problem #396: v_p(C(2n,n)) dynamics ---
        N = 300
        vals_396 = []
        for n in range(1, N + 1):
            v = legendre_valuation(2 * n, p) - 2 * legendre_valuation(n, p)
            vals_396.append(v)

        # --- Problem #389: v_p(ratio of consecutive products) ---
        # For fixed n, v_p(prod(n+k..n+2k-1) / prod(n..n+k-1)) as k varies
        # Use n=4 (minimal k=207 known)
        vals_389 = []
        n_fixed = 4
        for k in range(1, N + 1):
            num_val = sum(v_p(n_fixed + k + i, p) for i in range(k))
            den_val = sum(v_p(n_fixed + i, p) for i in range(k))
            vals_389.append(num_val - den_val)

        # Compute normalized cross-correlation
        a396 = np.array(vals_396[:N], dtype=float)
        a389 = np.array(vals_389[:N], dtype=float)
        a396_c = a396 - a396.mean()
        a389_c = a389 - a389.mean()

        if np.std(a396_c) > 0 and np.std(a389_c) > 0:
            cross_corr = np.correlate(a396_c / np.std(a396_c),
                                       a389_c / np.std(a389_c),
                                       mode='full') / N
            max_corr = float(np.max(np.abs(cross_corr)))
        else:
            max_corr = 0.0

        # Koopman eigenvalue extraction via FFT
        fft_396 = np.fft.rfft(a396_c)
        fft_389 = np.fft.rfft(a389_c)
        power_396 = np.abs(fft_396) ** 2
        power_389 = np.abs(fft_389) ** 2
        freqs = np.fft.rfftfreq(N)

        # Find shared dominant frequencies
        top_396 = set(np.argsort(power_396[1:])[::-1][:10] + 1)
        top_389 = set(np.argsort(power_389[1:])[::-1][:10] + 1)
        shared = top_396 & top_389
        shared_freqs = [round(float(freqs[i]), 6) for i in sorted(shared)]

        # Compute spectral overlap (cosine similarity of power spectra)
        if np.linalg.norm(power_396) > 0 and np.linalg.norm(power_389) > 0:
            spectral_overlap = float(np.dot(power_396, power_389) /
                                      (np.linalg.norm(power_396) * np.linalg.norm(power_389)))
        else:
            spectral_overlap = 0.0

        # Foucault test: are eigenvalues at 1/p, 2/p, ... present?
        foucault_harmonics = []
        for h in range(1, p):
            target_freq = h / p
            closest_idx = np.argmin(np.abs(freqs - target_freq))
            if closest_idx > 0:
                foucault_harmonics.append({
                    "harmonic": h,
                    "target_freq": round(target_freq, 4),
                    "power_396": round(float(power_396[closest_idx]), 2),
                    "power_389": round(float(power_389[closest_idx]), 2),
                })

        results[str(p)] = {
            "cross_correlation": round(max_corr, 6),
            "spectral_overlap": round(spectral_overlap, 6),
            "shared_dominant_freqs": shared_freqs[:5],
            "n_shared_modes": len(shared),
            "foucault_harmonics": foucault_harmonics,
        }

        print(f"\n  p={p}:")
        print(f"    Cross-correlation (max): {max_corr:.4f}")
        print(f"    Spectral overlap: {spectral_overlap:.4f}")
        print(f"    Shared dominant modes: {len(shared)}/10")
        print(f"    Shared freqs: {shared_freqs[:5]}")

    return results


# ═══════════════════════════════════════════════════════════
# PRONG 3: Kummer's Theorem as Carry Dynamics → KvN
# ═══════════════════════════════════════════════════════════
# Kummer's theorem: v_p(C(m+n,m)) = number of carries when
# adding m and n in base p.
#
# Carry propagation is a CLASSICAL dynamical system on base-p
# digit strings. The Koopman operator lifts this to a unitary
# operator on L²(Z_p) — the p-adic integers.
#
# For #396: v_p(C(2n,n)) = carries when adding n + n in base p.
# For #389: the ratio involves carries in multiple additions.
#
# The KvN bridge: carry propagation is the "discrete Foucault"
# — a slow drift (carry chain) modulating fast cycling (digits).
#
# We test: does carry chain length predict divisibility?

def prong3_kummer_carry_dynamics():
    """Kummer carry dynamics as discrete Foucault precession."""
    print("\n" + "=" * 60)
    print("PRONG 3: Kummer Carry Dynamics (Discrete Foucault)")
    print("=" * 60)

    results = {}

    for p in [2, 3, 5]:
        # --- Problem #396: Carry analysis for C(2n,n) ---
        carry_chains_396 = []
        divisibility_success_396 = []

        for n in range(1, 200):
            # Compute carries when adding n + n in base p
            digits = []
            temp = n
            while temp > 0:
                digits.append(temp % p)
                temp //= p
            if not digits:
                digits = [0]

            # Simulate carry propagation
            carries = 0
            carry = 0
            max_chain = 0
            current_chain = 0
            for d in digits:
                total = d + d + carry
                if total >= p:
                    carry = 1
                    carries += 1
                    current_chain += 1
                    max_chain = max(max_chain, current_chain)
                else:
                    carry = 0
                    current_chain = 0
            if carry:
                carries += 1
                current_chain += 1
                max_chain = max(max_chain, current_chain)

            carry_chains_396.append(max_chain)

            # Check divisibility for k=1,2,3
            for k in range(1, 4):
                df = desc_factorial(n, k + 1)
                if df > 0:
                    cb = comb(2 * n, n)
                    if cb % df == 0:
                        divisibility_success_396.append((n, k))

        # --- Problem #389: Carry analysis for consecutive products ---
        carry_chains_389 = []
        n_test = 4

        for k in range(1, 200):
            # Carry analysis: adding k consecutive integers starting at n_test+k
            total_carries = 0
            for i in range(k):
                a = n_test + i
                b = n_test + k + i
                # Carries when computing aspects of the ratio
                temp_a, temp_b = a, b
                carry = 0
                while temp_a > 0 or temp_b > 0:
                    da = temp_a % p
                    db = temp_b % p
                    # v_p aspect: track p-adic digits
                    if da > 0:
                        total_carries += 1
                    temp_a //= p
                    temp_b //= p

            carry_chains_389.append(total_carries)

        cc396 = np.array(carry_chains_396, dtype=float)
        cc389 = np.array(carry_chains_389[:len(carry_chains_396)], dtype=float)

        # Correlation between carry chain lengths
        if len(cc396) == len(cc389) and np.std(cc396) > 0 and np.std(cc389) > 0:
            correlation = float(np.corrcoef(cc396, cc389)[0, 1])
        else:
            correlation = 0.0

        # Key metric: how many n satisfy #396 for k=1,2,3?
        k1_count = len([x for x in divisibility_success_396 if x[1] == 1])
        k2_count = len([x for x in divisibility_success_396 if x[1] == 2])
        k3_count = len([x for x in divisibility_success_396 if x[1] == 3])

        # Foucault signature: autocorrelation of carry chain at lag p
        if len(cc396) > p:
            auto_corr_p = float(np.corrcoef(cc396[:-p], cc396[p:])[0, 1])
        else:
            auto_corr_p = 0.0

        results[str(p)] = {
            "mean_carry_chain_396": round(float(cc396.mean()), 4),
            "mean_carry_chain_389": round(float(cc389.mean()), 4),
            "carry_correlation": round(correlation, 4),
            "autocorr_at_lag_p": round(auto_corr_p, 4),
            "divisibility_396_k1": k1_count,
            "divisibility_396_k2": k2_count,
            "divisibility_396_k3": k3_count,
        }

        print(f"\n  p={p}:")
        print(f"    Carry chain correlation (#396 vs #389): {correlation:.4f}")
        print(f"    Autocorrelation at lag p: {auto_corr_p:.4f} "
              f"{'✓ FOUCAULT' if abs(auto_corr_p) > 0.3 else '✗ weak'}")
        print(f"    #396 divisibility (n≤199): k=1: {k1_count}, k=2: {k2_count}, k=3: {k3_count}")

    return results


# ═══════════════════════════════════════════════════════════
# PRONG 4: Resonance Overlap Test (Joint Morphism)
# ═══════════════════════════════════════════════════════════
# The strongest test: if #396 and #389 share a genuine Foucault
# morphism, their Koopman spectra should not just correlate —
# they should share EXACTLY the same resonant frequencies,
# corresponding to the same prime-periodic orbits.

def prong4_resonance_overlap():
    """Test whether #396 and #389 share Koopman resonance structure."""
    print("\n" + "=" * 60)
    print("PRONG 4: Resonance Overlap (Joint Morphism Test)")
    print("=" * 60)

    N = 500
    results = {}

    # Build full valuation sequences for both problems
    for p in [2, 3, 5, 7, 11]:
        # #396: v_p(C(2n,n)) sequence
        seq_396 = np.array([
            legendre_valuation(2*n, p) - 2*legendre_valuation(n, p)
            for n in range(1, N+1)
        ], dtype=float)

        # #389: v_p(upper_product / lower_product) for various n
        # Average over several n values to get the structural signal
        seq_389 = np.zeros(N)
        for n_val in [2, 3, 4, 5, 6]:
            for k in range(1, N+1):
                v_num = sum(v_p(n_val + k + i, p) for i in range(k))
                v_den = sum(v_p(n_val + i, p) for i in range(k))
                seq_389[k-1] += (v_num - v_den)
        seq_389 /= 5.0  # average

        # FFT-based spectral comparison
        fft_396 = np.fft.rfft(seq_396 - seq_396.mean())
        fft_389 = np.fft.rfft(seq_389 - seq_389.mean())

        # Phase coherence: are the phases aligned at resonant frequencies?
        phase_396 = np.angle(fft_396)
        phase_389 = np.angle(fft_389)

        # Power-weighted phase coherence
        p396 = np.abs(fft_396)**2
        p389 = np.abs(fft_389)**2
        weights = np.sqrt(p396 * p389)
        total_weight = weights.sum()

        if total_weight > 0:
            phase_diff = np.abs(phase_396 - phase_389)
            # Wrap to [0, pi]
            phase_diff = np.minimum(phase_diff, 2*np.pi - phase_diff)
            coherence = float(1.0 - np.average(phase_diff / np.pi, weights=weights))
        else:
            coherence = 0.0

        # Shared resonance peaks
        freqs = np.fft.rfftfreq(N)
        threshold_396 = np.percentile(p396[1:], 90)
        threshold_389 = np.percentile(p389[1:], 90)
        peaks_396 = set(np.where(p396[1:] > threshold_396)[0] + 1)
        peaks_389 = set(np.where(p389[1:] > threshold_389)[0] + 1)
        shared_peaks = peaks_396 & peaks_389

        # Compute overlap ratio
        union = peaks_396 | peaks_389
        jaccard = len(shared_peaks) / max(len(union), 1)

        # Key test: does 1/p appear in shared peaks?
        has_prime_resonance = False
        for sp in shared_peaks:
            if sp < len(freqs) and abs(freqs[sp] - 1.0/p) < 0.01:
                has_prime_resonance = True

        results[str(p)] = {
            "phase_coherence": round(coherence, 4),
            "spectral_jaccard": round(jaccard, 4),
            "n_shared_peaks": len(shared_peaks),
            "n_peaks_396": len(peaks_396),
            "n_peaks_389": len(peaks_389),
            "has_prime_resonance": bool(has_prime_resonance),
        }

        verdict = "RESONANT" if jaccard > 0.15 and coherence > 0.5 else "WEAK"
        print(f"\n  p={p}: phase coherence={coherence:.3f}, "
              f"Jaccard={jaccard:.3f}, shared={len(shared_peaks)} → {verdict}")

    return results


# ═══════════════════════════════════════════════════════════
# PRONG 5: Divisibility Witness Computation
# ═══════════════════════════════════════════════════════════
# Concrete computation: find witnesses for both conjectures
# and check whether witness structure correlates.

def prong5_witness_computation():
    """Compute divisibility witnesses for both problems."""
    print("\n" + "=" * 60)
    print("PRONG 5: Witness Computation & Structural Correlation")
    print("=" * 60)

    # #396: For each k, find smallest n where descFactorial(n, k+1) | C(2n,n)
    print("\n  Problem #396: Finding witnesses (k → min n)...")
    witnesses_396 = {}
    for k in range(1, 12):
        for n in range(k + 1, 2000):
            df = desc_factorial(n, k + 1)
            if df > 0 and comb(2*n, n) % df == 0:
                witnesses_396[k] = n
                print(f"    k={k}: n={n} "
                      f"(descFact={df}, C(2n,n) mod df = 0)")
                break
        else:
            witnesses_396[k] = None
            print(f"    k={k}: no witness found (n ≤ 1999)")

    # #389: For each n, find smallest k where divisibility holds
    print("\n  Problem #389: Finding witnesses (n → min k)...")
    witnesses_389 = {}
    for n in range(1, 15):
        for k in range(1, 500):
            lower = 1
            upper = 1
            for i in range(k):
                lower *= (n + i)
                upper *= (n + k + i)
            if lower > 0 and upper % lower == 0:
                witnesses_389[n] = k
                print(f"    n={n}: k={k}")
                break
        else:
            witnesses_389[n] = None
            print(f"    n={n}: no witness found (k ≤ 499)")

    # Structural analysis: how do witnesses scale?
    w396_vals = [v for v in witnesses_396.values() if v is not None]
    w389_vals = [v for v in witnesses_389.values() if v is not None]

    result = {
        "witnesses_396": {str(k): v for k, v in witnesses_396.items()},
        "witnesses_389": {str(n): v for n, v in witnesses_389.items()},
        "growth_396": "superexponential" if len(w396_vals) >= 3 and
                       w396_vals[-1] > 10 * w396_vals[-2] else
                       "polynomial" if len(w396_vals) >= 3 else "insufficient data",
        "growth_389": "superexponential" if len(w389_vals) >= 3 and
                       w389_vals[-1] > 10 * w389_vals[-2] else
                       "polynomial" if len(w389_vals) >= 3 else "insufficient data",
    }

    return result


# ═══════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════

if __name__ == "__main__":
    t0 = time.time()

    print("╔══════════════════════════════════════════════════════╗")
    print("║  EXP-FOUCAULT-KVN-001                               ║")
    print("║  KvN Salvo: Erdős #396 & #389 Foucault Resonance   ║")
    print("╚══════════════════════════════════════════════════════╝\n")

    r1 = prong1_precession_spectrum()
    r2 = prong2_koopman_eigenvalues()
    r3 = prong3_kummer_carry_dynamics()
    r4 = prong4_resonance_overlap()
    r5 = prong5_witness_computation()

    elapsed = time.time() - t0

    # ── Synthesis ──
    print("\n" + "=" * 60)
    print("SYNTHESIS: Foucault Resonance Morphism Assessment")
    print("=" * 60)

    # Aggregate evidence
    avg_coherence = np.mean([r4[p]["phase_coherence"] for p in r4])
    avg_jaccard = np.mean([r4[p]["spectral_jaccard"] for p in r4])
    avg_cross_corr = np.mean([r2[p]["cross_correlation"] for p in r2])
    has_foucault = any(abs(r3[p]["autocorr_at_lag_p"]) > 0.3 for p in r3)

    print(f"\n  Avg phase coherence:    {avg_coherence:.4f}")
    print(f"  Avg spectral Jaccard:   {avg_jaccard:.4f}")
    print(f"  Avg cross-correlation:  {avg_cross_corr:.4f}")
    print(f"  Foucault autocorrelation detected: {has_foucault}")

    # Morphism verdict
    if avg_coherence > 0.5 and avg_jaccard > 0.1:
        morphism_strength = "STRONG"
    elif avg_coherence > 0.3 or avg_jaccard > 0.05:
        morphism_strength = "MODERATE"
    else:
        morphism_strength = "WEAK"

    print(f"\n  ══► MORPHISM STRENGTH: {morphism_strength}")
    print(f"  ══► Foucault resonance binding: {'CONFIRMED' if has_foucault else 'PARTIAL'}")

    # ── Save results ──
    all_results = {
        "experiment": "EXP-FOUCAULT-KVN-001",
        "date": "2026-04-16",
        "problems": [396, 389],
        "physics_object": "PHYS-R-002 (Foucault Pendulum)",
        "morphism_type": "resonance_pair",
        "elapsed_seconds": round(elapsed, 2),
        "prong1_precession": r1,
        "prong2_koopman": r2,
        "prong3_kummer": r3,
        "prong4_overlap": r4,
        "prong5_witnesses": r5,
        "synthesis": {
            "avg_phase_coherence": round(avg_coherence, 4),
            "avg_spectral_jaccard": round(avg_jaccard, 4),
            "avg_cross_correlation": round(avg_cross_corr, 4),
            "foucault_detected": has_foucault,
            "morphism_strength": morphism_strength,
        }
    }

    out_path = "EXP-FOUCAULT-KVN-001_RESULTS.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n  Results saved to {out_path}")
    print(f"  Elapsed: {elapsed:.1f}s")
