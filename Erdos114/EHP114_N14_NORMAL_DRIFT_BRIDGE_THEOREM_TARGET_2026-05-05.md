# EHP114 n=14 Normal-Drift Bridge Theorem Target

Date: 2026-05-05

Required claim language: shadow signature, not universal law.

## Meaning

The accepted Taylor-collar row contract is now strong enough to carry the
normal-drift error budget on the known hard n=14 subcell `(6,4)`. The finite
part is no longer the main blocker:

```text
row_count = 15793
row_failure_count = 0
normal_drift_error_sum = 2.4736484172357147
normal_drift_budget = 2.5620009612530126
normal_drift_slack = 0.08835254401729786
normal_error / budget = 0.9655142424403749
```

What this means in plain math: if the row-wise Taylor-collar estimates are
accepted and if normal drift is the right comparison functional, then the
selected hard subcell has enough additive slack to stay under the exact-length
reserve budget.

What it does not mean: this is not an exact lemniscate-length certificate and
not a proof of Erdős #114. The analytic comparison theorem is still missing.

## Source Artifacts

Finite row contract:

```text
erdos-experiments/Erdos114/row_contract/EXP-MATH-EHP114-N14-TAYLOR-COLLAR-ROW-CONTRACT-20260505-01/
```

Source diagnostic:

```text
erdos-experiments/scripts/erdos-114/bridge-diagnostic-taylor-collar-local-bound6-output-z16/
```

Source result SHA-256:

```text
accc2124a77606187ec73fd000d98da517fbb3adc3a1b9ae5ffea5c45332356b
```

Row contract SHA-256:

```text
8417f5903d3696fe300961359e65e4771c1c94269d067e4dd81db7341c422e47
```

## Theorem Target

There are two layers. The first is finite and certificate-shaped. The second is
analytic and is the current bottleneck.

### Layer 1: Finite Normal-Drift Budget

Plain statement:

```text
For every retained z16 Taylor-collar row B in the `(6,4)` root-affine parameter
cell, the local 6 x 6 bound resolves regularity and provides a normal-drift
error contribution e_B. If the finite sum of these e_B is at most
2.5620009612530126, then the normal-drift bridge error for the cell is within
the available reserve.
```

Lean-shaped statement:

```lean
structure BridgeRow where
  normalDriftError : ℚ
  exactLengthError : ℚ
  regularityResolved : Prop

def FiniteNormalDriftContract {ι : Type} [Fintype ι]
    (rows : ι → BridgeRow) (budget : ℚ) : Prop :=
  (∀ i, (rows i).regularityResolved ∧
    (rows i).exactLengthError ≤ (rows i).normalDriftError) ∧
  (∑ i, (rows i).normalDriftError) ≤ budget

theorem finite_normal_drift_contract_implies_exact_error_budget
    {ι : Type} [Fintype ι] (rows : ι → BridgeRow) (budget : ℚ)
    (h : FiniteNormalDriftContract rows budget) :
    (∑ i, (rows i).exactLengthError) ≤ budget
```

This layer has a small Lean scratch file at:

```text
erdos-experiments/Erdos114/lean_scratch/Ehp114NormalDriftBridgeScratch.lean
```

It only checks the finite-sum implication pattern. It intentionally does not
encode complex analysis, exact lemniscate length, or the full EHP statement.

### Layer 2: Analytic Normal-Drift Bridge

Plain statement:

```text
For the n=14 root-affine family on the selected `(6,4)` subcell, the exact
lemniscate length discrepancy between the grid oracle and the true level set
is bounded above by the finite normal-drift row sum.
```

More explicit target:

```text
L_exact(C) <= L_ms(C) + sum_B normal_drift_error(B)
```

where `B` ranges over the retained Taylor-collar rows covering the relevant
level-set collar in the subcell `C`.

The current accepted row contract gives:

```text
sum_B normal_drift_error(B) <= 2.5620009612530126
```

The previous exact-length lift budget says the worst accepted subcell can
tolerate:

```text
L_exact(C) <= L_ms(C) + 2.5620009612569206
```

The margins are numerically aligned, but the comparison theorem has not been
proved. The residual slack between these two displayed budgets is tiny, so the
analytic theorem must be stated with conservative constants and no hidden
rounding assumptions.

## What Is Already Certified

The following are computationally certified by immutable local artifacts:

```text
Taylor-collar regularity unresolved cells = 0
local-bound subdivision = 6
row failures = 0
finite normal-drift row sum <= budget
hash check for row-contract RESULTS.json = OK
```

The following is Lean-checked as a structural scratch theorem:

```text
if row-wise exact error is bounded by row-wise normal drift, and the finite
normal-drift sum is within budget, then the finite exact-error sum is within
budget.
```

The following remains outside the current certificate:

```text
normal drift really bounds exact lemniscate length discrepancy
marching-squares length plus normal drift is a valid one-sided exact-length upper bound
the selected n=14 cell implies any general n theorem
```

## Current Bottleneck

The single bottleneck is now the analytic bridge:

```text
prove a coarea, tubular-neighborhood, or implicit-curve comparison theorem
that turns the Taylor-collar normal-drift row sum into an exact lemniscate
length upper error bound.
```

This is narrower than the earlier vague "exact length" blocker. The finite
data says the target constant is barely feasible; the proof still has to earn
the comparison.

## Next Attack

The next theorem another agent can attack without re-reading the corpus is:

```text
NormalDriftControlsExactLength14:
For the selected n=14 `(6,4)` root-affine subcell and retained Taylor-collar
cover, assume:
  1. |p'| is bounded below on each retained collar row,
  2. |p| and |p''| are bounded above on the local 6 x 6 subboxes,
  3. every exact level-set segment crossing the row lies inside the retained
     collar,
  4. endpoint and double-counting contributions are explicitly bounded.
Then:
  L_exact(C) <= L_ms(C) + sum_B normal_drift_error(B).
```

If this fails or the endpoint term exceeds the remaining `0.08835254401729786`
slack, the fallback should be the relative-length route or a direct validated
implicit-curve enclosure. Do not inflate the normal-drift claim to cover those
routes.

## Claim Ceiling

This packet records a local bridge theorem target and a finite-sum Lean scratch
check. It is not a public theorem, not an exact-length certificate, and not a
proof of Erdős #114.

Book-facing language should stay at:

```text
shadow signature, not universal law
```

