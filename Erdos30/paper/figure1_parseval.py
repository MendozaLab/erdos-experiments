#!/usr/bin/env python3
"""
Figure 1: λ_max vs k at fixed N — the spectral trajectory chart.

For each N in {10000, 100000}:
  - Build the full greedy Sidon set A = [a_0, ..., a_{k_max}]
  - For each prefix A_k = A[0:k], compute λ_max via FFT
  - Plot λ_max alongside:
    (a) Singer floor: k-1
    (b) Parseval budget: Nk - k²  (max possible λ_max)
    (c) Fourth moment ceiling: sqrt(2Nk² - Nk - k⁴)
    (d) Fitted divergence curve: c·k^{2.19}

Check: does the divergence curve ever exceed a binding constraint?
"""

import numpy as np
import math
import time
import json

def greedy_sidon(N):
    """Build a Sidon set in {0,...,N-1} greedily."""
    A = [0]
    sums = {0}
    for x in range(1, N):
        new_sums = []
        conflict = False
        for a in A:
            s = x + a
            if s in sums:
                conflict = True
                break
            new_sums.append(s)
        if not conflict:
            s_xx = x + x
            if s_xx in sums:
                continue
            new_sums.append(s_xx)
            for s in new_sums:
                sums.add(s)
            A.append(x)
    return A

def compute_lambda_max_fft(A_prefix, N):
    """Compute λ_max = max_{t≠0} |f̂_A(t)|² via FFT."""
    indicator = np.zeros(N)
    for a in A_prefix:
        indicator[a % N] = 1.0
    fft = np.fft.fft(indicator)
    spectrum = np.abs(fft) ** 2
    return np.max(spectrum[1:])

def main():
    all_results = {}

    for N in [10000, 100000]:
        print(f"\n{'='*70}")
        print(f"N = {N:,}  (√N = {math.sqrt(N):.1f})")
        print(f"{'='*70}")

        t0 = time.time()
        A = greedy_sidon(N)
        k_max = len(A)
        print(f"  Greedy set built: k_max = {k_max}, k_max/√N = {k_max/math.sqrt(N):.4f}")
        print(f"  Construction time: {time.time()-t0:.2f}s")

        results = []

        # Compute λ_max for every prefix from k=3 to k_max
        # For large N, do every prefix (it's fast with FFT)
        step = 1 if N <= 20000 else 2  # every other for N=100K to save time
        ks_to_test = list(range(3, k_max + 1, step))
        if k_max not in ks_to_test:
            ks_to_test.append(k_max)

        print(f"  Testing {len(ks_to_test)} prefix sizes...", flush=True)
        t1 = time.time()

        for k in ks_to_test:
            prefix = A[:k]
            lmax = compute_lambda_max_fft(prefix, N)

            # Constraint curves
            singer_floor = k - 1
            parseval_budget = N * k - k * k  # Σ_{t≠0} λ_t = Nk - k²
            E2 = 2 * k * k - k  # Sidon additive energy
            fourth_moment_total = N * E2 - k**4  # Σ_{t≠0} λ_t²
            fourth_moment_ceil = math.sqrt(max(fourth_moment_total, 0))

            results.append({
                'k': k,
                'k_over_sqrtN': k / math.sqrt(N),
                'lambda_max': float(lmax),
                'lambda_max_over_k': float(lmax / k),
                'singer_floor': singer_floor,
                'parseval_budget': parseval_budget,
                'fourth_moment_ceil': float(fourth_moment_ceil),
                'lmax_frac_parseval': float(lmax / parseval_budget) if parseval_budget > 0 else 0,
                'lmax_frac_4th': float(lmax / fourth_moment_ceil) if fourth_moment_ceil > 0 else 0,
            })

        elapsed = time.time() - t1
        print(f"  Computed in {elapsed:.1f}s")

        # Print trajectory table
        print(f"\n  {'k':>5s} {'k/√N':>6s} {'λ_max':>10s} {'λ_max/k':>8s} "
              f"{'Singer':>8s} {'4th ceil':>10s} {'%Parse':>7s} {'%4th':>7s}")
        print(f"  {'-'*5} {'-'*6} {'-'*10} {'-'*8} {'-'*8} {'-'*10} {'-'*7} {'-'*7}")

        for r in results[::max(1, len(results)//25)]:  # print ~25 rows
            print(f"  {r['k']:5d} {r['k_over_sqrtN']:6.3f} {r['lambda_max']:10.2f} "
                  f"{r['lambda_max_over_k']:8.3f} {r['singer_floor']:8d} "
                  f"{r['fourth_moment_ceil']:10.1f} "
                  f"{r['lmax_frac_parseval']*100:6.3f}% "
                  f"{r['lmax_frac_4th']*100:6.3f}%")

        # Power law fit on λ_max/k vs k (excluding small k)
        fit_data = [(r['k'], r['lambda_max_over_k']) for r in results if r['k'] >= 5]
        if len(fit_data) >= 5:
            ks = np.array([d[0] for d in fit_data])
            lmks = np.array([d[1] for d in fit_data])
            alpha, beta = np.polyfit(np.log(ks), np.log(lmks), 1)
            print(f"\n  Power law fit (k≥5): λ_max/k ~ k^{alpha:.4f}")
            print(f"  → λ_max ~ k^{1+alpha:.4f}")
            print(f"  → At k=√N: λ_max ~ N^{(1+alpha)/2:.4f}")

        # Key analysis: does λ_max approach any constraint?
        print(f"\n  CONSTRAINT ANALYSIS:")
        final = results[-1]
        print(f"  At k_max={final['k']}:")
        print(f"    λ_max = {final['lambda_max']:.1f}")
        print(f"    Parseval budget = {final['parseval_budget']:,}")
        print(f"    λ_max/budget = {final['lmax_frac_parseval']*100:.4f}%")
        print(f"    4th moment ceiling = {final['fourth_moment_ceil']:.1f}")
        print(f"    λ_max/4th_ceil = {final['lmax_frac_4th']*100:.4f}%")
        print(f"    Singer floor = {final['singer_floor']}")
        print(f"    λ_max/Singer = {final['lambda_max']/final['singer_floor']:.2f}")

        # Extrapolate: what would λ_max be at k = √N?
        if len(fit_data) >= 5:
            k_target = math.sqrt(N)
            lmax_extrap = math.exp(beta) * k_target ** (1 + alpha)
            parseval_at_sqrtN = N * k_target - k_target**2
            E2_at_sqrtN = 2 * k_target**2 - k_target
            fourth_at_sqrtN = math.sqrt(N * E2_at_sqrtN - k_target**4)

            print(f"\n  EXTRAPOLATION TO k = √N = {k_target:.1f}:")
            print(f"    λ_max (extrapolated) = {lmax_extrap:.1f}")
            print(f"    Parseval budget = {parseval_at_sqrtN:.1f}")
            print(f"    λ_max/budget = {lmax_extrap/parseval_at_sqrtN*100:.4f}%")
            print(f"    4th moment ceiling = {fourth_at_sqrtN:.1f}")
            print(f"    λ_max/4th_ceil = {lmax_extrap/fourth_at_sqrtN*100:.4f}%")
            print(f"    Singer floor = {k_target-1:.1f}")
            print(f"    λ_max/Singer = {lmax_extrap/(k_target-1):.2f}")

            if lmax_extrap > parseval_at_sqrtN:
                print(f"    ★★★ PARSEVAL CONTRADICTION: extrapolated λ_max EXCEEDS budget!")
            elif lmax_extrap > fourth_at_sqrtN:
                print(f"    ★★ FOURTH MOMENT CONTRADICTION: extrapolated λ_max EXCEEDS √(Σλ²)!")
            else:
                print(f"    No contradiction from Parseval or 4th moment alone.")
                headroom_p = (parseval_at_sqrtN - lmax_extrap) / parseval_at_sqrtN
                headroom_4 = (fourth_at_sqrtN - lmax_extrap) / fourth_at_sqrtN
                print(f"    Parseval headroom: {headroom_p*100:.1f}%")
                print(f"    4th moment headroom: {headroom_4*100:.1f}%")

        all_results[N] = results

    # ─── Cross-N comparison ───
    print(f"\n{'='*70}")
    print("CROSS-N COMPARISON: scaling exponents")
    print(f"{'='*70}")

    for N, results in all_results.items():
        fit_data = [(r['k'], r['lambda_max']) for r in results if r['k'] >= 5]
        ks = np.array([d[0] for d in fit_data])
        lms = np.array([d[1] for d in fit_data])
        gamma, c = np.polyfit(np.log(ks), np.log(lms), 1)
        print(f"  N={N:>7,}: λ_max ~ k^{gamma:.4f} (range k={int(ks[0])}-{int(ks[-1])})")

    # ─── Attack vector assessment ───
    print(f"\n{'='*70}")
    print("ATTACK VECTOR ASSESSMENT")
    print(f"{'='*70}")

    print("""
  Attack A (Uniqueness of flatness):
    λ_max = k-1 iff Singer (PDS)
    Data: Singer has λ_max = k-1 exactly. Greedy has λ_max >> k.
    Status: The "if" is classical. The "only if" needs proof.
    Testable: check non-PDS Sidon sets at Singer moduli.

  Attack B (Parseval divergence contradiction):
    Need: extrapolated λ_max > Parseval budget at k ≈ √N
    Check the numbers above — is the contradiction real?

  Attack C (λ_max ≥ k-1 lower bound):
    Prove: every Sidon set has λ_max ≥ k-1.
    This follows from: Σ_{t≠0} λ_t = Nk - k², mean = (Nk-k²)/(N-1) ≈ k.
    By pigeonhole: max ≥ mean iff all non-negative.
    But λ_t ≥ 0 always (they're |f̂|²). So λ_max ≥ (Nk-k²)/(N-1).
    For k ≈ √N: λ_max ≥ (N^{3/2}-N)/(N-1) ≈ √N ≈ k.
    Actually: λ_max ≥ k·(N-k)/(N-1) ≈ k for k << N.
    So λ_max ≥ k-1 is nearly FREE from pigeonhole!
    The hard part: is λ_max = k-1 achievable only for PDS?
    """)

    # Verify pigeonhole lower bound
    print("  Pigeonhole verification:")
    for N, results in all_results.items():
        for r in results[-3:]:
            k = r['k']
            pigeonhole_lb = k * (N - k) / (N - 1)
            print(f"    N={N}, k={k}: λ_max={r['lambda_max']:.1f}, "
                  f"pigeonhole_LB={pigeonhole_lb:.1f}, "
                  f"ratio={r['lambda_max']/pigeonhole_lb:.3f}")

    # Save all results
    serializable = {str(N): results for N, results in all_results.items()}
    with open('/tmp/figure1_data.json', 'w') as f:
        json.dump(serializable, f, indent=2)
    print(f"\n  Results saved to /tmp/figure1_data.json")

if __name__ == '__main__':
    main()
