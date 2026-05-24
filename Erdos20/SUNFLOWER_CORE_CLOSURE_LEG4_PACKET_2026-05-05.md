# Sunflower Core-Closure Leg-4 Packet

**Date:** 2026-05-05
**Problem:** Erdos #20, sunflower conjecture
**Packet purpose:** Convert the existing PMF/lattice-gas measurements into a focused Maxwell-style core-closure Leg-4 experiment.
**Status:** Internal protocol packet. No public claim. No new theorem. No scorecard or D1 update.

## Executive Claim Ceiling

This packet retires the weak "publishable observation" framing for the April 17 sunflower runs. The measured values are useful as a Leg-4 signal stream, not as lower-bound news.

The right target is:

> Erdos #20 exhibits a candidate core-closure shadow: the shared core must carry enough information to keep the petals coherent. The Maxwell analogy is a narrow displacement-current-style closure term, not a proof of the sunflower conjecture.

The phrase to preserve is: **shadow signature, not universal law**.

Current A-axis remains **A0**. Safe verbs are still formalize, encode, compute, measure, observe, recover, and specify. Do not say prove, solve, advance, improve, resolve, or establish unless a later artifact contains a real theorem beyond encoding and is verified through the normal Lean/build pipeline.

## Inputs Read

- `EXP-MATH-ERDOS20-SUNFLOWER-001_REPORT_2026-04-17.md`
- `EXP-MATH-ERDOS20-SUNFLOWER-001_RESULTS.json`
- `EXP-MATH-ERDOS20-SUNFLOWER-002_REPORT.md`
- `EXP-MATH-ERDOS20-SUNFLOWER-002_RESULTS.json`
- `Q1_LITERATURE_GATE_2026-04-17.md`
- `Math-Problems/proofs/archive/Sunflower_MDL.lean`
- `SHADOW_DYNAMICS_HYSTERETIC_INFORMATION_EXHAUST_LIVING_ANALYSIS.md`

## What Core Closure Cost Means

A k-sunflower is not only "petals around a core." It is a coherence condition: k sets agree on a shared intersection C, and after removing C the petals are pairwise disjoint.

For fixed ground set [n], uniformity w, petal count k=3, and a family F of w-subsets, a candidate core C of size s carries the information needed to answer:

1. Which previously admitted sets contain C?
2. What are their petals A \ C?
3. Which petal pairs are already disjoint?
4. Would a new candidate set B containing C close a forbidden triple by making three petals pairwise disjoint?

The **core closure cost** is the minimal bookkeeping needed to keep those answers correct as the family grows. Operationally, it is the information carried by the shared core channel so that petals are not treated as independent fragments when they are actually constrained by the same C.

The aggregate version can be measured from the existing enumeration outputs:

```text
N(n,w)      = binomial(n,w)
D_n,w(m)   = number of sunflower-free families of size m
g_n,w(m)   = mean valid extensions from size m to m+1
p_safe(m)  = g_n,w(m) / (N(n,w) - m), when N(n,w) > m
I_close(m) = -log2 p_safe(m)
```

`I_close(m)` is the number of bits needed to distinguish a safe unused addition from a random unused addition at family size m. It is not yet a per-core cost; it is the aggregate closure pressure of the sunflower constraint.

The per-core version is the next instrumentation target:

```text
p_safe(C,s,m)  = valid additions through core C / candidate additions through core C
I_core(C,s,m)  = -log2 p_safe(C,s,m)
I_core(s,m)    = distribution of I_core(C,s,m) over all cores C of size s
```

That per-core distribution is the true core-closure observable. The current April data gives enough to define the aggregate channel and jamming points, but not enough to claim per-core Maxwell behavior.

## Narrow Maxwell Analogy

The analogy is only this:

Maxwell's displacement current is a closure term. It repairs a naive current law that would otherwise leave the field equations incoherent across a gap.

The sunflower core-closure term plays the same structural role in the packet. The naive picture separates petals, but the shared core has to carry state: which petals have already used this core, which petal supports are disjoint, and whether adding another petal closes a forbidden 3-sunflower.

This does **not** say:

- Maxwell's equations prove anything about sunflowers.
- The sunflower conjecture follows from thermodynamics.
- The PMF lattice gas improves known sunflower bounds.
- The April small-n values are new lower bounds.

It says only that the sunflower system has a candidate missing-bookkeeping term whose behavior can be tested against a thermodynamic-floor classifier. If the test passes, the result is Leg-4 morphism evidence. If it fails, the Maxwell analogy stays literary/internal.

## Literature Gate Constraint

The Q1 gate blocks any new-lower-bound framing.

The existing measurements live on the Erdos-Rado w-uniform axis, not the Erdos-Szemeredi power-set axis. That axis correction is useful, but it does not make the measured lower bounds novel. Abbott-Hansen-Sauer dominates the small-n measurements: the gate records the classical lower bound

```text
f(w,3) >= 10^(w/2 - O(log w)),
so c_3 >= sqrt(10) ~= 3.162.
```

The recovered April values are weaker:

```text
w=2: M = 6,  c_3 >= sqrt(6) ~= 2.449
w=3: M = 12 tentative, c_3 >= 12^(1/3) ~= 2.289
w=4: M >= 15 at n=7, c_3 >= 15^(1/4) ~= 1.968
```

Therefore:

- No public observation note from these numbers.
- No claim of improved Erdos-Rado lower bounds.
- No claim that M(infinity,3,3)=12 or M(infinity,4,3)=15 is established.
- The useful contribution is the Leg-4 question: does the core-closure observable show a thermodynamic-floor pattern when measured against known strong constructions and geometry-exhausting sweeps?

## Existing Measured Sunflower Quantities

From `EXP-MATH-ERDOS20-SUNFLOWER-001_RESULTS.json`:

```text
w=3, k=3, n=4..8
N(n,3): 4, 10, 20, 35, 56
Z_n:    16, 388, 33652, 2485795, 148790380
M(n):   4, 6, 10, 12, 12
jamming m*: n=4 -> 2, n=5 -> 3, n=6 -> 5, n=7 -> 6
```

From `EXP-MATH-ERDOS20-SUNFLOWER-002_RESULTS.json`:

```text
w=2: n=3..10, saturation observed at M=6 from n=6 onward
w=3: n=4..8, M=12 at n=7,8, tentative only
w=4: n=5..7, M=15 at n=7, n=8 aborted under exhaustive enumeration
w=4, n=7 jamming m*: 9
```

These are aggregate lattice-gas quantities. They support PMF-amenability and a visible jamming transition, but they do not execute Leg 4.

## Leg-4 Experiment Specification

### Question

When sunflower-free families approach the jamming boundary, does the information needed to close the shared core behave like a displacement-current-style closure term, or does it remain ordinary combinatorial geometry?

### Hypothesis

The core-closure observable may show a thermodynamic-floor shadow when geometry is exhausted. If real, the signal should appear as a stable closure-cost scaling regime, not merely as small-n saturation or a lower-bound calculation.

### Inputs

Required existing inputs:

- `EXP-MATH-ERDOS20-SUNFLOWER-001_RESULTS.json`
- `EXP-MATH-ERDOS20-SUNFLOWER-002_RESULTS.json`
- `Q1_LITERATURE_GATE_2026-04-17.md`
- `sunflower_transfer_matrix.cpp`

Required new computation inputs:

- A real transfer-matrix or symmetry-reduced enumerator. Exhaustive backtracking is acceptable only for calibration.
- Sweeps over `w >= 3` with enough n to see whether the signal stabilizes after small-n effects.
- Per-core instrumentation for core sizes `s = 0..w-1`.
- Abbott-Hansen-Sauer-style construction baselines where available, so the comparison is not against weak small-n families.

### Outputs

Every real run must emit the normal experiment contract:

```text
*_RESULTS.json
*_REPORT.md
*_RESULTS.sha256
```

The results JSON should include at minimum:

```text
experiment_id
run_date
n, w, k
N = binomial(n,w)
Z_n
density_of_states D(m)
growth_rates g(m)
safe_fraction p_safe(m)
aggregate I_close(m)
per_core I_core(C,s,m) summaries by core size s
jamming m*
thermodynamic_floor model parameters
Mendoza Limit comparison table
classification: PASS / FAIL / INCONCLUSIVE
claim_ceiling
```

### Thermodynamic Floor Comparison

Use the Mendoza Limit only as a classifier probe:

```text
E_L(I,T) = I * k_B * T * ln 2
M_L(I,T) = I * k_B * T * ln 2 / c^2
```

Here `I` is the measured closure information in bits, preferably `I_core(C,s,m)` near the jamming boundary. The physical conversion does not make the mathematics physical by itself. It only asks whether the measured closure cost has floor-like scaling after the ordinary geometry has been accounted for.

The normalized diagnostic should be reported as:

```text
R(n,w,s,m) = measured closure cost / Mendoza-floor cost
```

where the numerator must be defined before the run. Acceptable numerators:

1. Energy-equivalent cost of erasing `I_core` bits under an explicitly stated T.
2. A dimensionless information ratio comparing observed `I_core` to the floor-predicted bit budget.
3. A construction-normalized closure cost using Abbott-Hansen-Sauer baseline families.

Do not treat `R ~= O(1)` alone as success. Small systems can make false friends. Scaling decides the class.

### PASS Criteria

PASS means "Leg-4 morphism evidence for core closure," not "progress on the sunflower conjecture."

A run can PASS only if all of the following hold:

1. Per-core closure instrumentation exists; aggregate-only data is insufficient.
2. The literature baseline is explicit, and Abbott-Hansen-Sauer is not being beaten or rebranded.
3. At least one geometry-exhausting sweep shows stable late-regime behavior across three consecutive n-values or an equivalent symmetry-reduced progression.
4. The late-regime normalized ratio `R` is stable with coefficient of variation <= 0.25 across the accepted window.
5. The log-log slope of `R` against ambient state count `N(n,w)` has absolute value <= 0.15 across the accepted window.
6. The same qualitative closure-cost pattern appears in at least two core sizes or two w-values, unless a pre-registered reason says only one core size is relevant.

### FAIL Criteria

FAIL means the Maxwell/core-closure analogy is not supported by the measured Leg-4 channel.

Classify FAIL if any of the following hold:

1. The signal is explainable by ordinary ambient geometry: `I_core` tracks `N(n,w)`, family size, or density without a stable floor-class regime.
2. `R` drifts monotonically with absolute log-log slope > 0.25 across the late window.
3. The apparent signal disappears when compared against Abbott-Hansen-Sauer-style baseline constructions.
4. The result depends on small-n saturation only.
5. Per-core instrumentation contradicts the aggregate reading.

### INCONCLUSIVE Criteria

Classify INCONCLUSIVE if:

1. Only aggregate growth rates are available.
2. The run remains at the current April scale (`w <= 4`, `n <= 8` with missing w=4 n=8).
3. The transfer-matrix rewrite is not done and exhaustive enumeration blocks the sweep.
4. Literature status for exact small-w extremals is not settled.
5. The Mendoza-floor numerator is not defined before the run.

The current 2026-04-17 data is therefore **INCONCLUSIVE for Leg 4**. It is useful setup, not evidence of a power morphism.

## Lean Status and Formal Claim Ceiling

`Math-Problems/proofs/archive/Sunflower_MDL.lean` contains useful formal vocabulary:

- `IsSunflower`
- `sunflower_pigeonhole`
- `erdos_rado_sunflower`
- `sets_sharing_kernel_bound`

For this packet, that file was read as a definition and proof-pattern source only. No current `lake build` was run, and this packet does not quote it as a current COMPILED artifact.

No optional Lean scratch file was created. The next low-risk Lean task, if needed, is only definitional: encode `CoreClosureState` and `safeExtensionFraction` as structures around finite families, with no conjecture claim and no theorem beyond basic projections.

## Retired Framing

Retire:

> The April sunflower measurements are a publishable observation or a new lower-bound contribution.

Replace with:

> The April sunflower measurements expose the right observable for a Leg-4 experiment: core closure cost. The Maxwell analogy is worth testing only as a closure-term signature. The claim ceiling is shadow signature, not universal law.

## Next Artifact Target

Create a new experiment only after the transfer-matrix/per-core instrumentation is ready. The recommended artifact shape is:

```text
EXP-MATH-ERDOS20-CORE-CLOSURE-LEG4-001_RESULTS.json
EXP-MATH-ERDOS20-CORE-CLOSURE-LEG4-001_REPORT.md
EXP-MATH-ERDOS20-CORE-CLOSURE-LEG4-001_RESULTS.sha256
```

Minimum viable run:

1. Reproduce April aggregate values for w=3 as calibration.
2. Add per-core summaries for w=3, n=7 and n=8.
3. Extend w=4 beyond n=7 using a real transfer matrix or symmetry reduction.
4. Compute `I_core(C,s,m)` near jamming.
5. Compare against the Mendoza-floor classifier and the Abbott-Hansen-Sauer baseline.
6. Return PASS / FAIL / INCONCLUSIVE with the claim ceiling embedded in the report.
