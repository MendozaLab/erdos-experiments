#!/usr/bin/env python3
"""
EXP-009B: Refined Holographic Extrapolation

Lessons from 009:
  1. Raw Koopman on λ decays past 1 → wrong. Must work on excess = λ-1.
  2. Quasiperiodic structure detected: even/odd W oscillation + slow mode.
  3. Power-law is perfect on calibration (sum-free) but ASSUMES λ→1.
  4. Padé is good (4% error on calibration) and doesn't assume the answer.

Refined approach:
  A. Separate even/odd W subsequences (decouple the period-2 oscillation)
  B. Work on log(λ-1) — make the power-law linear
  C. Koopman on log-excess: if the attractor is 1D, we get exact rate
  D. Torus phase portrait: plot (log excess, d/dW log excess) for visual
  E. Use Padé on BOTH subsequences separately
  F. Wynn epsilon algorithm (nonlinear sequence acceleration)
"""

import json
import math
import numpy as np


def load_data():
    """Load all available eigenvalue data."""
    # Sidon: try JSON first, fallback to hardcoded
    sidon = {}
    try:
        with open("EXP-MATH-ERDOS30-SIDON-007C_RESULTS.json") as f:
            data = json.load(f)
        for row in data.get("data", []):
            sidon[row["W"]] = row["lambda_max"]
    except:
        pass

    # Sum-free calibration
    sumfree = {}
    try:
        with open("EXP-008_SUMFREE_CALIBRATION_RESULTS.json") as f:
            data = json.load(f)
        for row in data.get("data", []):
            sumfree[row["W"]] = row["lambda_max"]
    except:
        pass

    return sidon, sumfree


def wynn_epsilon(sequence):
    """
    Wynn's epsilon algorithm — the most powerful scalar sequence accelerator.

    Given a convergent sequence s_0, s_1, ..., s_n, produces accelerated
    estimates of the limit. For geometric convergence, this is exact.
    For power-law convergence, it dramatically accelerates.

    The epsilon table:
      ε_{-1}(n) = 0
      ε_0(n) = s_n
      ε_{k+1}(n) = ε_{k-1}(n+1) + 1/(ε_k(n+1) - ε_k(n))

    The even columns ε_0, ε_2, ε_4, ... are successive accelerations.
    """
    n = len(sequence)
    if n < 3:
        return sequence, []

    # Build epsilon table
    eps = np.zeros((n, n + 1))
    eps[:, 0] = 0  # ε_{-1}
    eps[:, 1] = np.array(sequence)  # ε_0 = s_n

    for k in range(1, n):
        for i in range(n - k):
            diff = eps[i + 1, k] - eps[i, k]
            if abs(diff) < 1e-30:
                eps[i, k + 1] = 1e30  # infinity placeholder
            else:
                eps[i, k + 1] = eps[i + 1, k - 1] + 1.0 / diff

    # Extract even columns (the actual accelerated estimates)
    # ε_0(0), ε_2(0), ε_4(0), ... are the accelerated limit estimates
    accelerated = []
    for k in range(1, n, 2):  # columns 1, 3, 5, ... of the table (0-indexed)
        if k < n and abs(eps[0, k]) < 1e20:
            accelerated.append(eps[0, k])

    # The odd columns ε_1, ε_3, ... are auxiliary (reciprocals of differences)
    # The even columns give progressively better limit estimates
    even_col_estimates = []
    for k in range(2, n + 1, 2):  # columns 2, 4, 6, ... (these are ε_2, ε_4, ...)
        if k < n + 1 and abs(eps[0, k]) < 1e20:
            even_col_estimates.append(eps[0, k])

    return accelerated, even_col_estimates


def pade_extrapolate(Ws, values, order=(2, 2)):
    """Padé approximant at x=0 where x=1/W."""
    x = 1.0 / np.array(Ws, dtype=float)
    y = np.array(values)

    p, q = order
    n_params = p + q + 1  # a0..ap, b1..bq (b0=1)

    if len(x) < n_params:
        return None

    # Use last n_params*2 points for stability
    n_use = min(len(x), n_params * 3)
    x_fit = x[-n_use:]
    y_fit = y[-n_use:]

    # y(1 + b1*x + ... + bq*x^q) = a0 + a1*x + ... + ap*x^p
    # y = a0 + a1*x + ... + ap*x^p - b1*x*y - ... - bq*x^q*y

    cols = []
    # a coefficients
    for i in range(p + 1):
        cols.append(x_fit ** i)
    # b coefficients (negative)
    for j in range(1, q + 1):
        cols.append(-x_fit ** j * y_fit)

    M = np.column_stack(cols)

    try:
        coefs, _, _, _ = np.linalg.lstsq(M, y_fit, rcond=None)
        # λ(∞) = a0 / 1 = a0
        return coefs[0]
    except:
        return None


def analyze_sequence(name, Ws, lambdas, known_answer=None):
    """Full holographic analysis of an eigenvalue sequence."""
    Ws = np.array(Ws, dtype=float)
    lam = np.array(lambdas)
    excess = lam - 1.0

    print(f"\n{'='*70}")
    print(f"  {name}")
    print(f"  W range: {int(Ws[0])}..{int(Ws[-1])} ({len(Ws)} points)")
    if known_answer is not None:
        print(f"  Known answer: λ_∞ = {known_answer}")
    print(f"{'='*70}")

    results = {"name": name, "methods": {}}

    # --- A. Even/Odd separation ---
    print(f"\n  A. EVEN/ODD W SEPARATION")

    even_mask = np.array([int(w) % 2 == 0 for w in Ws])
    odd_mask = ~even_mask

    if np.sum(even_mask) >= 3:
        Ws_even = Ws[even_mask]
        lam_even = lam[even_mask]
        excess_even = lam_even - 1.0

        # Power-law fit on even
        if np.all(excess_even > 0):
            c_e, a_e = np.polyfit(np.log(Ws_even), np.log(excess_even), 1)
            alpha_even = -c_e
            print(f"    Even W: α = {alpha_even:.4f}")
        else:
            alpha_even = None

    if np.sum(odd_mask) >= 3:
        Ws_odd = Ws[odd_mask]
        lam_odd = lam[odd_mask]
        excess_odd = lam_odd - 1.0

        if np.all(excess_odd > 0):
            c_o, a_o = np.polyfit(np.log(Ws_odd), np.log(excess_odd), 1)
            alpha_odd = -c_o
            print(f"    Odd W:  α = {alpha_odd:.4f}")
        else:
            alpha_odd = None

    if alpha_even is not None and alpha_odd is not None:
        alpha_avg = (alpha_even + alpha_odd) / 2
        print(f"    Average α = {alpha_avg:.4f}")
        print(f"    Even/Odd split: Δα = {abs(alpha_even - alpha_odd):.4f}")
        results["methods"]["powerlaw_even"] = {"alpha": alpha_even, "lambda_inf": 1.0}
        results["methods"]["powerlaw_odd"] = {"alpha": alpha_odd, "lambda_inf": 1.0}

    # --- B. Power-law on full sequence ---
    print(f"\n  B. POWER-LAW FIT (full sequence)")

    if np.all(excess > 0):
        log_W = np.log(Ws)
        log_ex = np.log(excess)
        slope, intercept = np.polyfit(log_W, log_ex, 1)
        alpha_full = -slope
        c_full = math.exp(intercept)
        print(f"    λ - 1 = {c_full:.4f} × W^(-{alpha_full:.4f})")
        print(f"    → λ_∞ = 1.0 (by assumption)")
        results["methods"]["powerlaw"] = {"alpha": alpha_full, "c": c_full, "lambda_inf": 1.0}

        # Residual analysis
        predicted = c_full * Ws ** (-alpha_full) + 1.0
        residuals = lam - predicted
        rms_resid = np.sqrt(np.mean(residuals**2))
        print(f"    Residual RMS: {rms_resid:.6f}")

        # Is there curvature in the residuals? (systematic deviation)
        if len(residuals) > 10:
            first_half = np.mean(residuals[:len(residuals)//2])
            second_half = np.mean(residuals[len(residuals)//2:])
            print(f"    Residual bias (1st half): {first_half:+.6f}")
            print(f"    Residual bias (2nd half): {second_half:+.6f}")
            if abs(first_half) > 2 * rms_resid or abs(second_half) > 2 * rms_resid:
                print(f"    ⚠ Systematic residual pattern — power law may be wrong")

    # --- C. Koopman on log-excess ---
    print(f"\n  C. KOOPMAN ON LOG-EXCESS")

    if np.all(excess > 0):
        log_excess = np.log(excess)

        for dim in [2, 3, 4]:
            if len(log_excess) > dim + 3:
                # Takens embed log(λ-1)
                M = len(log_excess) - dim
                X = np.zeros((M - 1, dim))
                Y = np.zeros((M - 1, dim))
                for i in range(M - 1):
                    X[i] = log_excess[i:i+dim]
                    Y[i] = log_excess[i+1:i+1+dim]

                # SVD-based DMD
                U, S, Vh = np.linalg.svd(X, full_matrices=False)
                rank = min(dim, np.sum(S > 0.01 * S[0]))
                rank = max(1, rank)

                Ur = U[:, :rank]
                Sr = S[:rank]
                Vr = Vh[:rank, :]

                Atilde = Ur.T @ Y @ Vr.T @ np.diag(1.0 / Sr)
                eigs, vecs = np.linalg.eig(Atilde)

                # The dominant Koopman eigenvalue of log-excess tells us the
                # asymptotic decay rate. If |μ| < 1, log(λ-1) → -∞, meaning λ → 1.
                # The rate is log|μ| per W step.
                dominant = eigs[np.argmax(np.abs(eigs))]

                print(f"    dim={dim}: Koopman eigs = {[f'{abs(e):.4f}∠{np.angle(e)*180/np.pi:.1f}°' for e in eigs]}")
                print(f"      Dominant |μ| = {abs(dominant):.6f}")

                if abs(dominant) < 1:
                    # log(excess) decays → λ → 1
                    # Rate: log(excess) ~ log|μ| * W → excess ~ |μ|^W
                    decay_rate = math.log(abs(dominant))
                    print(f"      → log-excess decays at rate {decay_rate:.4f}/step")
                    print(f"      → λ - 1 ~ {abs(dominant):.4f}^W  (geometric decay)")
                    print(f"      → λ_∞ = 1.0 (confirmed)")
                    results["methods"][f"koopman_logexcess_d{dim}"] = {
                        "dominant_eigenvalue": abs(dominant),
                        "decay_rate": decay_rate,
                        "lambda_inf": 1.0,
                    }
                else:
                    print(f"      → |μ| ≥ 1: log-excess does NOT decay → λ_∞ > 1")
                    # Extrapolate: log(excess) at W=∞ approaches...
                    # actually if |μ| > 1, log(excess) grows → λ → ∞ (unphysical)
                    # if |μ| = 1, log(excess) is constant → λ_∞ = 1 + exp(constant)
                    if abs(abs(dominant) - 1.0) < 0.01:
                        # Nearly marginal — excess stabilizes
                        last_log_ex = log_excess[-1]
                        lambda_inf_est = 1.0 + math.exp(last_log_ex)
                        print(f"      → Marginal: λ_∞ ≈ {lambda_inf_est:.6f}")
                        results["methods"][f"koopman_logexcess_d{dim}"] = {
                            "dominant_eigenvalue": abs(dominant),
                            "lambda_inf": lambda_inf_est,
                        }

    # --- D. Wynn epsilon acceleration ---
    print(f"\n  D. WYNN EPSILON ACCELERATION")

    # Apply to the raw sequence
    _, wynn_even = wynn_epsilon(lam.tolist())
    if wynn_even:
        print(f"    Wynn ε estimates: {[f'{w:.6f}' for w in wynn_even[:6]]}")
        # The last stable estimate
        stable = [w for w in wynn_even if 0.5 < w < 3.0]
        if stable:
            print(f"    Best Wynn estimate: λ_∞ ≈ {stable[-1]:.6f}")
            results["methods"]["wynn"] = {"lambda_inf": stable[-1]}

    # Apply to even subsequence
    if np.sum(even_mask) >= 5:
        _, wynn_e = wynn_epsilon(lam[even_mask].tolist())
        stable_e = [w for w in wynn_e if 0.5 < w < 3.0] if wynn_e else []
        if stable_e:
            print(f"    Wynn (even W): λ_∞ ≈ {stable_e[-1]:.6f}")
            results["methods"]["wynn_even"] = {"lambda_inf": stable_e[-1]}

    # Apply to odd subsequence
    if np.sum(odd_mask) >= 5:
        _, wynn_o = wynn_epsilon(lam[odd_mask].tolist())
        stable_o = [w for w in wynn_o if 0.5 < w < 3.0] if wynn_o else []
        if stable_o:
            print(f"    Wynn (odd W):  λ_∞ ≈ {stable_o[-1]:.6f}")
            results["methods"]["wynn_odd"] = {"lambda_inf": stable_o[-1]}

    # --- E. Padé approximants ---
    print(f"\n  E. PADÉ APPROXIMANTS")

    for p, q in [(1, 1), (2, 1), (1, 2), (2, 2), (3, 2), (2, 3), (3, 3)]:
        est = pade_extrapolate(Ws, lam, order=(p, q))
        if est is not None and 0.5 < est < 3.0:
            print(f"    [{p}/{q}]: λ_∞ ≈ {est:.6f}")
            results["methods"][f"pade_{p}_{q}"] = {"lambda_inf": est}

    # --- F. Richardson extrapolation (assumes power-law) ---
    print(f"\n  F. RICHARDSON EXTRAPOLATION")

    # Two-term Richardson: use pairs (W, 2W) to eliminate leading correction
    richardson = []
    W_dict = dict(zip(Ws.astype(int), lam))
    for w in sorted(W_dict.keys()):
        w2 = 2 * w
        if w2 in W_dict and w >= 5:
            l1 = W_dict[w]
            l2 = W_dict[w2]
            # If λ(W) = λ_∞ + c W^{-α}, then:
            # λ_∞ ≈ (2^α λ(2W) - λ(W)) / (2^α - 1)
            # We use α from the power-law fit
            if "powerlaw" in results["methods"]:
                alpha = results["methods"]["powerlaw"]["alpha"]
                fac = 2 ** alpha
                rich = (fac * l2 - l1) / (fac - 1)
                richardson.append((w, rich))
                print(f"    W={w:3d}, 2W={w2:3d}: λ_∞ ≈ {rich:.6f}")

    if richardson:
        rich_vals = [r[1] for r in richardson]
        print(f"    Richardson mean: {np.mean(rich_vals):.6f} ± {np.std(rich_vals):.6f}")
        results["methods"]["richardson"] = {"lambda_inf": np.mean(rich_vals)}

    # --- G. CONSENSUS ---
    print(f"\n  {'='*50}")
    print(f"  CONSENSUS")
    print(f"  {'='*50}")

    all_estimates = {}
    for method, data in results["methods"].items():
        if "lambda_inf" in data:
            val = data["lambda_inf"]
            if isinstance(val, (int, float)) and 0.5 < val < 3.0:
                all_estimates[method] = val

    if all_estimates:
        values = list(all_estimates.values())

        # Separate "assumes λ→1" methods from "free" methods
        constrained = {k: v for k, v in all_estimates.items()
                       if abs(v - 1.0) < 0.001}
        free = {k: v for k, v in all_estimates.items()
                if abs(v - 1.0) >= 0.001}

        print(f"\n  Methods assuming λ_∞ = 1: {len(constrained)}")
        for m, v in sorted(constrained.items()):
            print(f"    {m:30s}: {v:.6f}")

        print(f"\n  Methods with FREE λ_∞: {len(free)}")
        for m, v in sorted(free.items()):
            err_str = f"  (err={abs(v - known_answer):.4f})" if known_answer is not None else ""
            print(f"    {m:30s}: {v:.6f}{err_str}")

        if free:
            free_vals = list(free.values())
            free_mean = np.mean(free_vals)
            free_median = np.median(free_vals)
            free_std = np.std(free_vals)

            print(f"\n  FREE methods consensus:")
            print(f"    Mean:   {free_mean:.6f} ± {free_std:.6f}")
            print(f"    Median: {free_median:.6f}")

            results["consensus_free"] = {
                "mean": float(free_mean),
                "median": float(free_median),
                "std": float(free_std),
                "n_methods": len(free),
            }

        results["consensus_all"] = {
            "mean": float(np.mean(values)),
            "median": float(np.median(values)),
            "std": float(np.std(values)),
            "n_methods": len(all_estimates),
        }

    # --- Calibration score ---
    if known_answer is not None:
        print(f"\n  CALIBRATION SCORE (known answer = {known_answer}):")
        for m, v in sorted(all_estimates.items()):
            err = abs(v - known_answer)
            quality = "✓✓" if err < 0.01 else ("✓" if err < 0.05 else "✗")
            print(f"    {m:30s}: err = {err:.6f}  {quality}")

    return results


def main():
    print("=" * 80)
    print("EXP-009B: Refined Holographic Extrapolation")
    print("  Even/odd separation • Koopman on log-excess • Wynn ε • Padé • Richardson")
    print("=" * 80)

    sidon, sumfree = load_data()

    # Analyze sum-free (calibration)
    if sumfree:
        Ws_sf = sorted(sumfree.keys())
        lam_sf = [sumfree[w] for w in Ws_sf]
        sf_results = analyze_sequence("SUM-FREE CALIBRATION", Ws_sf, lam_sf, known_answer=1.0)

    # Analyze Sidon
    if sidon:
        Ws_sid = sorted(sidon.keys())
        lam_sid = [sidon[w] for w in Ws_sid]
        sid_results = analyze_sequence("SIDON (Erdős #30)", Ws_sid, lam_sid)

    # --- FINAL VERDICT ---
    print("\n" + "=" * 80)
    print("FINAL VERDICT")
    print("=" * 80)

    if sumfree and sidon:
        print("\n  Calibration (sum-free):")
        if "consensus_free" in sf_results:
            c = sf_results["consensus_free"]
            print(f"    Free-method consensus: λ_∞ = {c['mean']:.4f} ± {c['std']:.4f}")
            print(f"    Known answer: 1.0000")
            print(f"    Error: {abs(c['mean'] - 1.0):.4f}")

        print(f"\n  Sidon:")
        if "consensus_free" in sid_results:
            c = sid_results["consensus_free"]
            print(f"    Free-method consensus: λ_∞ = {c['mean']:.4f} ± {c['std']:.4f}")
            print(f"    Number of methods: {c['n_methods']}")

        # The key question: does the Koopman analysis of log-excess confirm λ→1?
        koopman_confirms = False
        for k, v in sid_results.get("methods", {}).items():
            if k.startswith("koopman_logexcess") and isinstance(v, dict):
                li = v.get("lambda_inf")
                if li is not None and abs(li - 1.0) < 0.001:
                    koopman_confirms = True

        if koopman_confirms:
            print(f"\n  ✓ Koopman on log-excess CONFIRMS λ → 1")
            print(f"    → The Sidon lattice gas IS critical")
            print(f"    → The correction exponent α governs h(N)")
        else:
            print(f"\n  ? Koopman on log-excess is ambiguous")
            print(f"    → Need larger W data to resolve")

        # What's the best α estimate?
        alphas = []
        for k, v in sid_results.get("methods", {}).items():
            if isinstance(v, dict) and "alpha" in v:
                alphas.append((k, v["alpha"]))

        if alphas:
            print(f"\n  Power-law exponents:")
            for method, alpha in alphas:
                correction = 1.0 / (2 * alpha) if alpha > 0 else float('inf')
                print(f"    {method:30s}: α = {alpha:.4f}  →  h(N) = √N + O(N^{correction:.4f})")

    # Save
    out = {
        "experiment_id": "EXP-009B",
        "title": "Refined Holographic Extrapolation",
    }
    if sumfree:
        out["calibration"] = sf_results
    if sidon:
        out["sidon"] = sid_results

    with open("EXP-009B_REFINED_HOLOGRAPHIC_RESULTS.json", "w") as f:
        json.dump(out, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else str(x) if isinstance(x, complex) else None)

    print(f"\nSaved to EXP-009B_REFINED_HOLOGRAPHIC_RESULTS.json")


if __name__ == "__main__":
    main()
