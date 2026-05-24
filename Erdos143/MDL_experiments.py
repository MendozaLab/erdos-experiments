"""Numerical experiments for M2, M3, M5 morphism candidates on near-WS sets."""
import math, random, json
from collections import Counter

def build_ws(N_max, target=1500, seed=137, k_max=40):
    rng = random.Random(seed)
    A = []
    attempts = 0
    while len(A) < target and attempts < 600000:
        attempts += 1
        u = rng.uniform(math.log(1.5), math.log(N_max))
        c = math.exp(u)
        ok = True
        for x in A:
            if not ok: break
            for k in range(1, k_max+1):
                if abs(k*x - c) < 1.0 or abs(k*c - x) < 1.0:
                    ok = False; break
        if ok: A.append(c)
    return sorted(A)

# Build A at N_max = 1e8, target |A|=1500
import time
t0 = time.time()
A = build_ws(1e8, target=1500, seed=137)
elapsed = time.time() - t0
print(f"WS set built: |A|={len(A)}, max(A)={A[-1]:.2e}, in {elapsed:.1f}s")

# ============================================================
# M2 — Multiplicative Sidon analogue test
# Predicted bound: |A ∩ [1,N]| ≪ √N (Sidon rate)
# Test: ratio |A ∩ [1,N]| / sqrt(N) across N decades
# ============================================================
print("\n=== M2 — Sidon-rate test (|A∩[1,N]| / sqrt(N)) ===")
m2_data = []
for N in [1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]:
    cnt = sum(1 for x in A if x <= N)
    ratio = cnt / math.sqrt(N)
    m2_data.append({"N": N, "count": cnt, "ratio_over_sqrtN": ratio})
    print(f"  N={N:>10.0f}  count={cnt:>5d}  |A|/√N = {ratio:.4f}")

# ============================================================
# M3 — Log-position autocorrelation test
# Predicted: low autocorrelation of {log x} under integer log-shifts
# Test: compute |sum_{x∈A} exp(i·t·log x)|² at integer t (lattice-mode strength)
# If low for all t > 0: aperiodic. If spikes at some t: hidden periodicity.
# ============================================================
print("\n=== M3 — Log-position autocorrelation (lattice modes) ===")
logs = [math.log(x) for x in A]
m3_data = []
for t in range(1, 11):  # integer log-shifts t=1..10
    re_s = sum(math.cos(t*l) for l in logs)
    im_s = sum(math.sin(t*l) for l in logs)
    power = (re_s**2 + im_s**2) / len(logs)
    m3_data.append({"t": t, "power_per_element": power})
    print(f"  t={t:>3d}   |Σ exp(i·t·log x)|²/|A| = {power:.4f}")
# baseline: for Poisson process, expected power ≈ 1 per element
print(f"  (Poisson baseline: ≈ 1.0 per element; lower = more aperiodic)")

# ============================================================
# M5 — Entropy decay test
# Predicted: 𝟙_A on dyadic windows has decreasing entropy rate
# Test: bin A into dyadic windows [2^k, 2^(k+1)), compute density;
#       compute Shannon entropy of binary indicator at each scale
# ============================================================
print("\n=== M5 — Dyadic-block entropy decay ===")
m5_data = []
for k in range(2, 27):  # 2^2 to 2^27
    lo, hi = 2**k, 2**(k+1)
    if hi > A[-1]: break
    # Discretize [lo, hi) into unit-integer bins (well-separated sets are nearly integer-spaced)
    bins = hi - lo  # number of integer slots
    occupied = sum(1 for x in A if lo <= x < hi)
    if bins <= 1 or occupied == 0:
        m5_data.append({"k": k, "density": 0.0, "binary_entropy_per_bin": 0.0})
        continue
    p = occupied / bins
    # Shannon entropy per bin (bits)
    if 0 < p < 1:
        H = -p * math.log2(p) - (1-p) * math.log2(1-p)
    else:
        H = 0.0
    m5_data.append({"k": k, "block": f"[2^{k}, 2^{k+1})", "bins": bins,
                    "occupied": occupied, "density": p,
                    "binary_entropy_per_bin": H, "total_block_entropy_bits": H * bins})
    print(f"  k={k:>2d}  block=[2^{k},2^{k+1})  bins={bins:>10d}  occ={occupied:>4d}  "
          f"density={p:.3e}  H_bin={H:.3e}  H_total≈{H*bins:.2f} bits")

# ============================================================
# Save
# ============================================================
out = {
    "ws_set": {"size": len(A), "max": A[-1], "build_seconds": elapsed},
    "M2_sidon_density": m2_data,
    "M3_log_autocorrelation": m3_data,
    "M5_entropy_decay": m5_data,
}
with open("/tmp/erdos143_mdl_experiments.json", "w") as f:
    json.dump(out, f, indent=2)
print(f"\nSaved /tmp/erdos143_mdl_experiments.json")
