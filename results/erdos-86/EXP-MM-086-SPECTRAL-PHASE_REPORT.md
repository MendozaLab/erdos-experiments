# EXP-MM-086-SPECTRAL-PHASE — Experiment Report

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-086-SPECTRAL-PHASE |
| Erdős Problem | #86 — C4-free subgraphs of the hypercube Q_n |
| Problem URL | https://www.erdosproblems.com/86 |
| Date | 2026-04-13 |
| Probe | T1-086 (Tier 1 — Euler-Lagrange) |
| Results file | `EXP-MM-086-SPECTRAL-PHASE_RESULTS.json` |
| SHA-256 | `ebe59d34b8ee54e1e19341146b2845378c4c3b217646008c133290613fd329ba` |
| Data integrity | REAL_COMPUTATION — all values from Monte Carlo with stated seed |

## Hypothesis

The critical edge-retention probability p_c (where the mean level spacing ratio ⟨r⟩ crosses the Poisson-to-GOE threshold at 0.4407) decreases with n, and the C4-free constraint shifts p_c upward by suppressing spectral entropy ("entropic cooling").

**Anti-Artifact classification (SOP 11):** This experiment measures a NUMERICAL OBSERVABLE (⟨r⟩ phase transition). The monotonic decrease of p_c(n) is an OBSERVATION, not a proof. The "entropic cooling" hypothesis is UNTESTED by this experiment — it requires a paired experiment (T2-086) comparing ⟨r⟩(p) with and without the C4-free constraint.

## Method

1. Construct the adjacency matrix of Q_n (hypercube on 2^n vertices, n·2^(n-1) edges)
2. For each edge retention probability p in [0.05, 0.95] with step 0.02 (46 points):
   - Generate 100 independent Bernoulli-thinned subgraphs (retain each edge with probability p)
   - Compute eigenvalues via `scipy.linalg.eigvalsh` (exact, symmetric real)
   - Compute the mean level spacing ratio ⟨r⟩ (Oganesyan-Huse diagnostic)
3. Locate the critical p_c where ⟨r⟩ crosses the GOE threshold (0.4407) via linear interpolation

**Observable:** Mean level spacing ratio ⟨r⟩ — a standard diagnostic in random matrix theory:
- Poisson ensemble (uncorrelated): ⟨r⟩ ≈ 0.386
- GOE (correlated, Wigner-Dyson): ⟨r⟩ ≈ 0.530
- Threshold: 0.4407 (midpoint)

**References:**
- Oganesyan & Huse, *Phys. Rev. B* 75, 155111 (2007)
- Atas et al., *Phys. Rev. Lett.* 110, 084101 (2013)

## Results

### Critical p_c values

| Graph | Vertices | Edges | p_c | Bracket | Peak ⟨r⟩ (at plateau) |
|---|---|---|---|---|---|
| Q_4 | 16 | 32 | 0.4878 | [0.47, 0.49] | ~0.475 |
| Q_5 | 32 | 80 | 0.3856 | [0.37, 0.39] | ~0.510 |
| Q_6 | 64 | 192 | 0.3169 | [0.31, 0.33] | ~0.519 |
| Q_7 | 128 | 448 | 0.2611 | [0.25, 0.27] | ~0.527 |
| Q_8 | 256 | 1024 | 0.2210 | [0.21, 0.23] | ~0.528 |

### Key observations

1. **p_c is monotonically decreasing with n.** The Poisson-to-GOE transition occurs at progressively lower edge retention as the hypercube grows. This is consistent with the increasing connectivity of Q_n — fewer random edges are needed to induce spectral rigidity in larger hypercubes.

2. **The GOE plateau sharpens with n.** For Q_4, the ⟨r⟩ curve is broad and noisy (std ~0.12-0.18), with a peak barely above 0.475. For Q_8, the transition is steep and the plateau is flat near 0.528 (std ~0.025-0.030), much closer to the theoretical GOE value of 0.530.

3. **Non-monotonic behavior at high p.** All curves show ⟨r⟩ declining after the plateau — the full Q_n adjacency matrix does NOT have GOE statistics. This is expected: Q_n is a highly structured graph (vertex-transitive, distance-regular), and its eigenvalues are known exactly ({n - 2k : k=0,...,n} with multiplicities C(n,k)). The "interesting" spectral regime is at intermediate p where disorder competes with hypercube structure.

4. **Convergence toward GOE.** The plateau ⟨r⟩ values (0.475, 0.510, 0.519, 0.527, 0.528) are converging toward the GOE value 0.530, suggesting that in the thermodynamic limit (n → ∞), the transition sharpens into a true phase transition.

### Scaling law (empirical fit)

A log-linear regression on {(n, p_c)} for n = 4, 5, 6, 7, 8:

p_c(n) ≈ C · n^α with approximate α ≈ -1.15

This is a NUMERICAL OBSERVATION, not a proven scaling law. Extrapolation beyond n=8 is speculation. The true scaling may involve log corrections or different functional forms.

### Comparison with Gemini agent results

A prior Gemini agent session reported p_c(Q_8) ≈ 0.192 ("~19.2%"). Our result is p_c(Q_8) = 0.221. The discrepancy is likely due to:
- Different trial counts (Gemini used fewer trials)
- Different threshold definitions or interpolation methods
- Possible coarser p-grid in the Gemini run

Our experiment uses 100 trials per p-value with p-step = 0.02, which provides reasonable statistical power for Q_8 (256×256 matrices).

### Comparison with prior MendozaLab results

EXP-ERDOS86-C4FREE-01 computed exact C4-free independent set sizes for n ≤ 5 and greedy approximations for n ≤ 10. That experiment addressed the combinatorial question (how large is f(n)?), while this experiment addresses the spectral question (at what edge density does the spectrum become rigid?). They probe different aspects of the same problem and are complementary, not redundant.

### Context: Minamoto 2026

Perplexity flagged arXiv:2603.29127 (Minamoto, 2026) as providing new lower bounds for Erdős #86. This paper should be retrieved and analyzed — it may affect the interpretation of p_c values if the new bounds constrain the C4-free density more tightly.

## Anti-Artifact Assessment (SOP 11 — Honest Accounting)

**What this experiment DOES establish:**
- The Poisson-to-GOE transition exists in random subgraphs of Q_n
- p_c decreases monotonically with n for 4 ≤ n ≤ 8
- The GOE plateau sharpens and converges toward 0.530 with increasing n
- The observable is well-defined and reproducible (seed=42, 100 trials)

**What this experiment does NOT establish:**
- Whether the C4-free constraint shifts p_c (that requires T2-086)
- Whether p_c has a closed-form scaling law (5 data points are insufficient)
- Whether the spectral transition has any bearing on the combinatorial problem f(n) ≥ ?
- Whether "entropic cooling" is a valid physical analogy or merely suggestive metaphor

**What would falsify the cooling hypothesis (T2-086):**
If ⟨r⟩(p) curves for C4-free subgraphs show the SAME p_c as unconstrained random subgraphs, the cooling hypothesis is dead. The C4-free constraint must demonstrably shift p_c upward (require more edges to reach GOE) for the analogy to hold.

## Probe Queue Update

- **T1-086** → ANSWERED (this experiment)
- **T2-086** (SEED, priority 35) → Next: design paired experiment comparing ⟨r⟩ with/without C4-free constraint
- **T3-086** (SEED, priority 25) → Depends on T2-086 outcome + physics bridge PHYS-R-002

## Compute

| Graph | Time |
|---|---|
| Q_4 | 0.2s |
| Q_5 | 0.4s |
| Q_6 | 1.2s |
| Q_7 | 4.1s |
| Q_8 | 17.7s |
| **Total** | **23.5s** |

Platform: macOS, numpy + scipy, seed=42.

## Next Steps

1. **T2-086** (highest-priority next probe): Run paired experiment — same setup but with C4-free constraint enforced during edge retention. Compare p_c values.
2. **Minamoto 2026**: Retrieve arXiv:2603.29127 and assess whether new bounds constrain the spectral interpretation.
3. **n=9,10 extension**: Q_9 (512 vertices) and Q_10 (1024 vertices) are computationally feasible (~2min and ~15min respectively). More data points would strengthen or weaken the empirical scaling.
4. **Formal spectral analysis**: The exact eigenvalue structure of Q_n is known. Compare the random-thinned spectrum against the Krawchouk polynomial structure to understand WHY p_c decreases.
