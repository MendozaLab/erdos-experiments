#!/usr/bin/env python3
"""
Generate figures for the Sidon von Neumann entropy paper.

Figure 1: λ_max/k vs k/√N at fixed N (non-monotone trajectory)
Figure 2: Singer anomaly — spectral gap grows with q
Figure 3: Entropy deficit vs N for Singer and greedy
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
import math

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 13,
    'legend.fontsize': 10,
    'figure.figsize': (7, 5),
    'text.usetex': False,
    'font.family': 'serif',
})

# ── Load data ──

with open('/tmp/figure1_data.json') as f:
    fig1_data = json.load(f)

with open('/tmp/spectral_universality_results.json') as f:
    univ_data = json.load(f)

# ══════════════════════════════════════════════════════════════
# Figure 1: λ_max/k vs k/√N — the spectral trajectory
# ══════════════════════════════════════════════════════════════

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

for idx, N_str in enumerate(['10000', '100000']):
    ax = [ax1, ax2][idx]
    results = fig1_data[N_str]
    N = int(N_str)
    sqrtN = math.sqrt(N)

    ks = np.array([r['k'] for r in results])
    k_norm = ks / sqrtN
    lmax_over_k = np.array([r['lambda_max'] / r['k'] for r in results])
    lmax = np.array([r['lambda_max'] for r in results])
    singer_floor = np.ones_like(ks)  # λ_max/k = 1 for Singer
    parseval_budget_over_k = np.array([(N * k - k**2) / k for k in ks])  # (Nk-k²)/k = N-k
    fourth_ceil_over_k = np.array([
        math.sqrt(max(N * (2*k**2 - k) - k**4, 0)) / k for k in ks
    ])

    # Main trajectory
    ax.plot(k_norm, lmax_over_k, 'b-', linewidth=1.5, label=r'$\lambda_{\max}/k$ (greedy)', zorder=3)

    # Singer floor
    ax.axhline(y=1, color='red', linestyle='--', linewidth=1.2, label='Singer floor ($\\lambda_{\\max}/k = 1$)')

    # Fourth moment ceiling (scaled)
    ax.plot(k_norm, fourth_ceil_over_k, 'g-.', linewidth=1, alpha=0.7, label=r'$\sqrt{\Sigma\lambda_t^2}/k$ (4th moment)')

    # Mark the peak
    peak_idx = np.argmax(lmax_over_k)
    ax.plot(k_norm[peak_idx], lmax_over_k[peak_idx], 'ro', markersize=8, zorder=5)
    ax.annotate(f'peak: k/{int(sqrtN)}≈{k_norm[peak_idx]:.2f}',
                xy=(k_norm[peak_idx], lmax_over_k[peak_idx]),
                xytext=(k_norm[peak_idx] + 0.05, lmax_over_k[peak_idx] * 0.85),
                fontsize=9, arrowprops=dict(arrowstyle='->', color='red'),
                color='red')

    ax.set_xlabel(r'$k / \sqrt{N}$')
    ax.set_ylabel(r'$\lambda_{\max} / k$')
    ax.set_title(f'N = {N:,}  ($\\sqrt{{N}}$ = {int(sqrtN)})')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)

fig.suptitle('Figure 1: Non-monotone spectral trajectory at fixed N', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('/tmp/sidon_entropy_paper/figure1_trajectory.pdf', bbox_inches='tight', dpi=300)
plt.savefig('/tmp/sidon_entropy_paper/figure1_trajectory.png', bbox_inches='tight', dpi=150)
plt.close()
print("Figure 1 saved.")

# ══════════════════════════════════════════════════════════════
# Figure 2: Singer anomaly — spectral gap grows with q
# ══════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(7, 5))

qs = sorted([int(q) for q in univ_data.keys()])
singer_ratios = []
random_ratios_mean = []
random_ratios_std = []

for q in qs:
    r = univ_data[str(q)]
    k = r['k_singer']
    singer_ratios.append(r['singer_lmax'] / k)
    if r['random_at_k_singer']:
        rand_vals = np.array(r['random_at_k_singer']) / k
        random_ratios_mean.append(np.mean(rand_vals))
        random_ratios_std.append(np.std(rand_vals))
    else:
        random_ratios_mean.append(float('nan'))
        random_ratios_std.append(0)

qs_arr = np.array(qs)
singer_arr = np.array(singer_ratios)
random_arr = np.array(random_ratios_mean)
random_std = np.array(random_ratios_std)

ax.plot(qs_arr, singer_arr, 'rs-', markersize=8, linewidth=2, label='Singer (PDS)', zorder=5)
ax.errorbar(qs_arr, random_arr, yerr=random_std, fmt='bo-', markersize=6,
            capsize=4, linewidth=1.5, label='Random Sidon (mean ± σ, n=20)')

# Fill between to show gap
ax.fill_between(qs_arr, singer_arr, random_arr, alpha=0.15, color='purple',
                label='Spectral gap (growing)')

ax.set_xlabel('Prime power $q$')
ax.set_ylabel(r'$\lambda_{\max} / k$')
ax.set_title('Figure 2: Singer anomaly at Singer moduli $N = q^2+q+1$')
ax.legend(loc='upper left')
ax.set_xticks(qs)
ax.grid(True, alpha=0.3)

# Annotate the gap ratios
for i, q in enumerate(qs):
    if not np.isnan(random_arr[i]):
        gap = random_arr[i] / singer_arr[i]
        ax.annotate(f'{gap:.1f}×', xy=(q, (singer_arr[i] + random_arr[i])/2),
                    fontsize=8, ha='center', color='purple', fontweight='bold')

plt.tight_layout()
plt.savefig('/tmp/sidon_entropy_paper/figure2_singer_anomaly.pdf', bbox_inches='tight', dpi=300)
plt.savefig('/tmp/sidon_entropy_paper/figure2_singer_anomaly.png', bbox_inches='tight', dpi=150)
plt.close()
print("Figure 2 saved.")

# ══════════════════════════════════════════════════════════════
# Figure 3: Entropy deficit (Singer shrinking, greedy constant)
# ══════════════════════════════════════════════════════════════

# Recompute deficit data from the universality results
fig, ax = plt.subplots(figsize=(7, 5))

# Singer deficits
singer_Ns = []
singer_deficits = []
for q in qs:
    r = univ_data[str(q)]
    N = r['N']
    k = r['k_singer']
    # Singer: all λ_t = k-1 for t≠0, λ_0 = k²
    # p_0 = k²/(Nk) = k/N, p_t = (k-1)/(Nk) for t≠0
    p0 = k / N
    pt = (k - 1) / (N * k)
    S1 = -p0 * np.log2(p0) - (N-1) * pt * np.log2(pt)
    S_max = np.log2(N)
    singer_Ns.append(N)
    singer_deficits.append(S_max - S1)

# Greedy deficits (from figure1 data, last point of each)
greedy_Ns_deficits = []
greedy_deficits_vals = []

# Recompute from raw greedy data
def greedy_sidon(N):
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

for N in [1000, 2000, 5000, 10000, 20000]:
    A = greedy_sidon(N)
    k = len(A)
    indicator = np.zeros(N)
    for a in A:
        indicator[a % N] = 1.0
    fft = np.fft.fft(indicator)
    spectrum = np.abs(fft) ** 2
    probs = spectrum / np.sum(spectrum)
    S1 = -np.sum(probs[probs > 0] * np.log2(probs[probs > 0]))
    S_max = np.log2(N)
    greedy_Ns_deficits.append(N)
    greedy_deficits_vals.append(S_max - S1)
    print(f"  Greedy N={N}: k={k}, deficit={S_max-S1:.4f}")

ax.plot(singer_Ns, singer_deficits, 'rs-', markersize=8, linewidth=2,
        label='Singer (PDS)')
ax.plot(greedy_Ns_deficits, greedy_deficits_vals, 'bo-', markersize=6,
        linewidth=1.5, label='Greedy (sub-extremal)')

ax.set_xlabel('$N$')
ax.set_ylabel('Entropy deficit $\\log_2 N - S_1(\\rho_A)$')
ax.set_title('Figure 3: Entropy deficit — Singer vanishes, greedy constant')
ax.set_xscale('log')
ax.legend()
ax.grid(True, alpha=0.3)

# Annotate
ax.annotate('Singer: deficit → 0\n(spectral flatness)',
            xy=(singer_Ns[-1], singer_deficits[-1]),
            xytext=(singer_Ns[-1]*0.3, singer_deficits[-1] + 0.15),
            fontsize=9, arrowprops=dict(arrowstyle='->', color='red'),
            color='red')

ax.annotate('Greedy: deficit ≈ 0.5 bits\n(sub-extremality floor)',
            xy=(greedy_Ns_deficits[-1], greedy_deficits_vals[-1]),
            xytext=(greedy_Ns_deficits[-1]*0.15, greedy_deficits_vals[-1] - 0.12),
            fontsize=9, arrowprops=dict(arrowstyle='->', color='blue'),
            color='blue')

plt.tight_layout()
plt.savefig('/tmp/sidon_entropy_paper/figure3_deficit.pdf', bbox_inches='tight', dpi=300)
plt.savefig('/tmp/sidon_entropy_paper/figure3_deficit.png', bbox_inches='tight', dpi=150)
plt.close()
print("Figure 3 saved.")

print("\nAll figures generated in /tmp/sidon_entropy_paper/")
