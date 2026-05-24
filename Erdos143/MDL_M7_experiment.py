"""M7 numerical experiment: VN entropy of bipartite WS dilation state.

Build |psi>_AK with A = greedy near-WS in [1,N], K = dilations [1..K_max].
Amplitude matrix M[x,k] = (1 if k*x in [1,N] else 0) / sqrt(Z).
Partial-trace out K register: rho_A = M @ M^T.
Compute eigenvalues, VN entropy, effective rank.

Compare to:
  - log(|A|) classical max entropy baseline
  - log(K_max) bipartite rank upper bound
  - Shannon entropy of row marginal H(P_A) where P_A(x) = sum_k M[x,k]^2
  - Mutual information I(A:K) = 2*S(rho_A) for pure states

Vary (|A|, N, K_max) to look for scaling pattern.
"""
import math, random, json
import numpy as np

def build_ws(N_max, target, seed=137, k_max_check=40):
    rng = random.Random(seed)
    A = []
    attempts = 0
    while len(A) < target and attempts < 200000:
        attempts += 1
        u = rng.uniform(math.log(1.5), math.log(N_max))
        c = math.exp(u)
        ok = True
        for x in A:
            if not ok: break
            for k in range(1, k_max_check+1):
                if abs(k*x - c) < 1.0 or abs(k*c - x) < 1.0:
                    ok = False; break
        if ok: A.append(c)
    return sorted(A)

def m7_experiment(A_size, N_max, K_max, seed=137):
    A = build_ws(N_max, A_size, seed=seed, k_max_check=30)
    n = len(A)
    A_arr = np.array(A)

    # Bipartite amplitude matrix M[x_idx, k-1] = 1 if k*A[x_idx] <= N_max else 0
    M = np.zeros((n, K_max))
    for k in range(1, K_max + 1):
        M[:, k-1] = (k * A_arr <= N_max).astype(float)
    Z = M.sum()  # total number of valid (x,k) pairs
    if Z == 0:
        return None
    M = M / math.sqrt(Z)

    # SVD: eigenvalues of rho_A = squared singular values
    s = np.linalg.svd(M, compute_uv=False)
    eig = s ** 2
    eig = eig[eig > 1e-15]
    S_VN = -float(np.sum(eig * np.log(eig)))  # natural log; convert to bits if desired

    # Marginal P_A(x) = sum_k M[x,k]^2 = row sum of M*M
    row_pmf = (M * M).sum(axis=1)
    row_pmf = row_pmf[row_pmf > 1e-15]
    H_Pa = -float(np.sum(row_pmf * np.log(row_pmf)))

    # Effective rank: exp(S_VN) gives "effective dimension"
    eff_rank = math.exp(S_VN)

    # Mutual information for pure state: I(A:K) = 2 S(rho_A)
    I_AK = 2 * S_VN

    return {
        "A_size": n,
        "N_max": N_max,
        "K_max": K_max,
        "Z_total_pairs": int(Z * Z) if Z * math.sqrt(Z) > 1 else None,  # rough
        "S_VN_rho_A": S_VN,
        "S_VN_bits": S_VN / math.log(2),
        "log_A_size": math.log(n),
        "log_K_max": math.log(K_max),
        "VN_over_log_A": S_VN / math.log(n) if n > 1 else 0,
        "VN_over_log_K": S_VN / math.log(K_max) if K_max > 1 else 0,
        "H_marginal_PA": H_Pa,
        "H_marginal_over_VN": H_Pa / S_VN if S_VN > 0 else 0,
        "effective_rank": eff_rank,
        "true_rank": int(np.sum(s > 1e-10)),
        "mutual_info_AK": I_AK,
        "top_5_eigenvalues": [float(x) for x in eig[:5]],
    }

# Scan: vary |A| and N together so |A|/sqrt(N) stays moderate
print("=== M7 bipartite VN-entropy scan ===\n")
print(f"{'|A|':>5} {'N':>9} {'K':>4} {'S_VN(rho_A)':>11} {'log|A|':>8} {'VN/log|A|':>10} "
      f"{'H(P_A)':>8} {'H/VN':>6} {'eff_rank':>8} {'true_rank':>9}")

rows = []
for A_size, N_max, K_max in [(100, 1e4, 20), (100, 1e4, 50), (200, 1e5, 50),
                              (200, 1e5, 100), (300, 1e6, 50), (500, 1e6, 100),
                              (500, 1e7, 100), (800, 1e7, 100), (1000, 1e8, 100)]:
    r = m7_experiment(A_size, N_max, K_max, seed=137)
    if r is None: continue
    rows.append(r)
    print(f"{r['A_size']:>5d} {r['N_max']:>9.0e} {r['K_max']:>4d} "
          f"{r['S_VN_rho_A']:>11.4f} {r['log_A_size']:>8.4f} "
          f"{r['VN_over_log_A']:>10.4f} {r['H_marginal_PA']:>8.4f} "
          f"{r['H_marginal_over_VN']:>6.3f} "
          f"{r['effective_rank']:>8.1f} {r['true_rank']:>9d}")

# Save
with open("/tmp/erdos143_m7_experiment.json", "w") as f:
    json.dump({"experiment": "M7_bipartite_VN_entropy",
               "setup": "indicator alpha(x,k)=1 if k*x<=N; uniform; rho_A=Tr_K |psi><psi|",
               "rows": rows}, f, indent=2)

print("\n=== Interpretation ===")
print("VN/log|A| ratio < 1 would indicate WS forces sub-maximal entanglement entropy")
print("(i.e., the bipartite state is LOW RANK relative to set size — novel signal vs Shannon M5)")
print("VN/log|A| close to 1 means no novel signal — VN entropy saturates the classical bound")
print("H(P_A) > S_VN means marginal entropy exceeds bipartite entanglement (typical)")
print("\nSaved: /tmp/erdos143_m7_experiment.json")
