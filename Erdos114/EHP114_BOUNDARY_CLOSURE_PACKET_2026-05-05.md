# EHP114 Boundary Closure Packet

Date: 2026-05-05
Worker: A
Scope: Erdos #114 radial hypergeometric/Puiseux boundary term
Status: closure packet, not a proof

## Bottom line

The current #114 boundary-correction sprint has a precise bottleneck:
**interval-hardening the radial hypergeometric/Puiseux singularity model**.

The radial file `erdos-experiments/Erdos30/lean/EhpRadialPuiseux.lean`
compiles and has no `sorry`, but the meaningful analytic work is still imported
through named axioms. The useful progress is that the remaining gap can now be
stated narrowly: replace the generic hypergeometric/Puiseux axioms by a fixed
n = 14 interval-certified theorem for the boundary singularity of
`p_a(z) = z^14 - a` as `a -> 1-`.

This is a shadow signature, not universal law. In #114 the "shadow" is the
boundary correction term: the local model is Puiseux/hypergeometric, not a
smooth Hessian. That is a research target, not a problem-resolution claim.

## Source files read

- `erdos-experiments/results/erdos-114/SCAFFOLD_CLOSURE_ASSESSMENT_2026-05-02.md`
- `erdos-experiments/Erdos30/lean/EhpRadialPuiseux.lean`
- `erdos-experiments/results/erdos-114/EXP-MATH-EHP114-N14-RADIAL-HYPERGEOMETRIC-CALIBRATION-20260502-01_REPORT.md`
- `erdos-experiments/results/erdos-114/EXP-MATH-EHP114-N14-RADIAL-HYPERGEOMETRIC-CALIBRATION-20260502-01_RESULTS.json`
- `erdos-experiments/results/erdos-114/EXP-MATH-EHP114-RADIAL-HYPERGEOMETRIC-20260502-01_REPORT.md`
- `erdos-experiments/results/erdos-114/EXP-MATH-EHP114-RADIAL-HYPERGEOMETRIC-20260502-01_RESULTS.json`
- `Math-Problems/Erdos-Standard/erdos-114/docs/*`
- `SHADOW_DYNAMICS_HYSTERETIC_INFORMATION_EXHAUST_LIVING_ANALYSIS.md`
- Local supporting files in this directory:
  `RADIAL_PUISEUX_FIT_2026-05-02.md`,
  `CROFTON_COAREA_DRAFT_2026-05-02.md`,
  `ATHENA_SPIKE_PROTOCOL_2026-05-02.md`

## What is already known

### Compiled

`EhpRadialPuiseux.lean` currently contains compiled Lean theorem bodies for:

```lean
theorem ehp_radial_L_form_n14
theorem ehp_radial_L_form_n14_boundary
theorem ehp_radial_puiseux_n14
```

Local check run during this packet:

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos30
lake env lean lean/EhpRadialPuiseux.lean
```

Result: PASS, no output.

Meaning: Lean accepts the reduction from the imported analytic axioms to the
n = 14 Puiseux corollary. It does **not** mean the analytic boundary theorem has
been proved in Mathlib.

### Probed

The radial-hypergeometric reports give exact-family diagnostics for
`p_a(z) = z^n - a`:

```text
L_n(a) = integral_0^(2*pi) |a + exp(i*t)|^(1/n - 1) dt
L_n(a) = 2*pi * 2F1(p,p;1;a^2), p = (n-1)/(2n), for |a| < 1
```

For n = 14, the exact boundary value reported is:

```text
L_14(1) = 30.852910841548542295844126311342204648430086988176
```

The n = 14 tail-window slopes approach `1/14`, with the four-point tail fit:

```text
0.071426652466 vs expected 0.071428571428571...
```

For n = 15, the same probe reports:

```text
0.066664875637 vs expected 1/15 = 0.066666666666...
```

Meaning: the radial boundary layer has the right Puiseux exponent. This is a
diagnostic and analytic-model confirmation, not a proof of global maximality.

### Certified

The finite verification frontier remains the existing interval/canonical
surface, not this packet. The standard-doc audit states:

```text
n = 3..14: canonical finite-verification frontier, with n=13 caveat
n = 15..16: noncanonical shortcut artifacts exist, claim status remains
             NEEDS_RECONCILIATION / NOT_PUBLIC_PROOF
```

The canonical n = 14 artifact is stronger than the shortcut surface because it
is tied to a large branch-and-bound run:

```text
bb_total_evals: 855638016
rigor: ieee_1788_interval_arithmetic_inari
```

This packet does not upgrade that frontier.

### Axiomatized

The following six axioms appear in `EhpRadialPuiseux.lean`:

```lean
axiom coarea_lemniscate_length
axiom fejer_riesz_factorization
axiom gauss_2F1_connection_at_one
axiom hypergeometric_lemniscate_radial
axiom hypergeometric_lemniscate_radial_boundary
axiom gauss_2F1_puiseux_lower_bound
```

Only the last three are load-bearing for `ehp_radial_puiseux_n14` as currently
written. The first three are architecture context or prior scaffold residue.

## Exact missing theorem

Plain math target, fixed n = 14:

Let

```text
F_14(z) = 2F1(13/28, 13/28; 1; z)
L_14(a) = length({z in C : |z^14 - a| = 1})
```

Then prove, with interval-certified constants, that:

```text
L_14(a) = 2*pi*F_14(a^2)              for 0 <= a < 1
L_14(1) = 2*pi*F_14(1)
```

and there exists an explicit `C14 > 0` such that for every
`0 < epsilon <= 1/4`,

```text
L_14(1) - L_14(1 - epsilon) >= C14 * epsilon^(1/14).
```

This is the boundary-correction theorem. It is radial only. It does not prove
that every nonradial polynomial perturbation pays at least this deficit.

Lean-shaped fixed theorem:

```lean
namespace Erdos114
namespace Radial

theorem ehp_radial_puiseux_n14_interval_hardened :
    exists C14 : R, 0 < C14 and
      forall epsilon : R, 0 < epsilon -> epsilon <= (1 / 4 : R) ->
        ehp_lemniscate_length (fun z : C => z^14 - 1) -
        ehp_lemniscate_length
          (fun z : C => z^14 - ((1 - epsilon : R) : C)) >=
          C14 * Real.rpow epsilon (1 / 14) := by
  -- no hypergeometric/Puiseux axiom here;
  -- proof must call certified fixed-n lemmas below
```

Lean-shaped narrower hypergeometric obligation:

```lean
theorem gauss_2F1_puiseux_lower_bound_n14_interval :
    exists K14 : R, 0 < K14 and
      forall w : R, 0 < w -> w <= (1 / 2 : R) ->
        gauss_2F1 (13/28) (13/28) 1 1 -
        gauss_2F1 (13/28) (13/28) 1 (1 - w) >=
          K14 * Real.rpow w (1 / 14) := by
  -- interval-hardened proof target
```

Lean-shaped boundary-value obligation:

```lean
theorem hypergeometric_lemniscate_radial_boundary_n14_interval :
    ehp_lemniscate_length (fun z : C => z^14 - 1) =
      2 * Real.pi * gauss_2F1 (13/28) (13/28) 1 1 := by
  -- boundary extension / Abel-dominated-convergence target,
  -- or an equivalent certified interval bridge to the boundary value
```

## Axiom triage

| Axiom | Current role | Reducible? | Honest status |
|---|---|---|---|
| `coarea_lemniscate_length` | Not used by current n = 14 radial theorem | Yes, for this packet | Classical analysis context. Generic coarea may exist in Mathlib pieces, but this polynomial lemniscate specialization is not needed for the present Lean proof. |
| `fejer_riesz_factorization` | Not used by current n = 14 radial theorem | Yes, for this packet | Useful for a later shape/cone route, not for the radial Puiseux boundary term. |
| `gauss_2F1_connection_at_one` | Not used by current theorem bodies | Yes, if the Puiseux bound remains the imported theorem | The limit theorem alone is too weak for a lower bound; it becomes useful only as part of proving the Puiseux-rate theorem. |
| `hypergeometric_lemniscate_radial` | Used for the interior formula | Partly | Too broad. Reduce to fixed n = 14 and eventually fixed `a = 1 - epsilon` if needed. The bridge from Hausdorff length to `2F1` is genuinely missing from current Mathlib. |
| `hypergeometric_lemniscate_radial_boundary` | Used for `a = 1` | Partly | Reduce to fixed n = 14. The endpoint continuity/Abel theorem bridge is the first real analytic boundary obligation. |
| `gauss_2F1_puiseux_lower_bound` | Used for the rate | Yes, by narrowing | Generic statement is overkill. Replace with fixed parameters `(13/28,13/28;1)` and a certified explicit positive constant. This is the central Mathlib-missing theorem for this sprint. |

The important reduction is not "prove every classical theorem now." It is:
delete all generic ambition from the sprint target and prove the single fixed
radial boundary inequality needed to recover the n = 14 singular model.

## Closure route tried

The current Lean proof of `ehp_radial_puiseux_n14` has this dependency shape:

```text
ehp_radial_puiseux_n14
  -> ehp_radial_L_form_n14
       -> hypergeometric_lemniscate_radial
  -> ehp_radial_L_form_n14_boundary
       -> hypergeometric_lemniscate_radial_boundary
  -> gauss_2F1_puiseux_lower_bound
```

I checked whether `gauss_2F1_connection_at_one` can replace
`gauss_2F1_puiseux_lower_bound`. It cannot. The connection-at-one axiom gives
only:

```text
F(z) -> F(1)
```

as `z -> 1`. It does not give:

```text
F(1) - F(1 - w) >= K * w^(1/14).
```

That rate is exactly the missing Puiseux boundary term. The honest narrowed
target is therefore not the generic DLMF theorem in full generality; it is the
fixed n = 14 interval theorem:

```text
G14(w) =
  (F_14(1) - F_14(1 - w)) / w^(1/14)
```

has a certified positive lower bound on `0 < w <= 1/2`.

### Proposed interval hardening split

1. Small interval near the singular point:

   Use the Gauss connection formula expansion for
   `a = b = 13/28`, `c = 1`, so `c - a - b = 1/14`. The singular coefficient is
   positive in the deficit convention because the raw connection coefficient
   multiplying `(1-z)^(1/14)` has the opposite sign. Prove:

   ```text
   G14(w) >= K_small > 0 for 0 < w <= w0.
   ```

2. Middle interval:

   Use Arb/interval arithmetic on the continuous function `G14(w)` over
   `[w0, 1/2]` and record a lower endpoint:

   ```text
   G14(w) >= K_mid > 0 for w0 <= w <= 1/2.
   ```

3. Combine:

   ```text
   K14 = min(K_small, K_mid)
   ```

4. Translate from `w` to radial contraction:

   For `a = 1 - epsilon`,

   ```text
   w = 1 - (1 - epsilon)^2 = epsilon * (2 - epsilon).
   ```

   On `0 < epsilon <= 1/4`, `w >= epsilon` and `w <= 1/2`, so

   ```text
   w^(1/14) >= epsilon^(1/14).
   ```

This route reduces `gauss_2F1_puiseux_lower_bound` from a generic analytic axiom
to a fixed-parameter interval certificate. That is the honest next theorem.

## What this would and would not close

If the narrowed theorem lands, it closes the **radial singularity model** for
n = 14:

```text
radial contraction from z^14 - 1 pays explicit Puiseux deficit
```

It does not close #114 because the full conjecture also needs:

```text
all nonradial admissible perturbations pay positive deficit
mixed radial/shape terms cannot cancel the Puiseux deficit
the local singular cone overlaps the certified/global region
```

The finite certificate already says n = 14 is handled by interval verification.
The value of the radial theorem is different: it explains the boundary layer in
a reusable way and gives the proof architecture a singular term instead of an
ordinary Hessian fantasy.

## Tao-grade claim ceiling

### Safe internal claim

```text
The n = 14 radial boundary model has a compiled Lean reduction from named
analytic axioms, and the numerical/exact-family probes identify the missing
boundary term as a Puiseux 1/14 deficit. The next proof obligation is
interval-hardening the fixed radial 2F1 lower bound.
```

### Safe public or collaborator phrasing, after care

Only after an interval certificate exists for `G14(w)`:

```text
We have isolated and certified the radial Puiseux boundary correction for the
EHP lemniscate at n = 14. This does not solve the full conjecture, but it
replaces the smooth-Hessian local model with the correct singular boundary
term and gives a sharper target for nonradial remainder bounds.
```

### What would justify Tao outreach

Outreach would be technically justified by all of the following:

1. A reproducible interval report for `gauss_2F1_puiseux_lower_bound_n14_interval`,
   with explicit constants and SHA-bound artifact.
2. A Lean file that no longer imports the generic Puiseux-rate axiom for n = 14.
3. A short note explaining why this boundary term is compatible with Tao/Fryntov-
   Nazarov Stokes-dispersion language.
4. A separate paragraph saying exactly what remains open: nonradial Fourier/
   tensor modes and mixed cone terms.

The outreach should ask for feedback on the boundary-correction mechanism, not
announce a small-n closure.

### Embarrassing overclaim

Do not say any of the following:

```text
Full-problem closure has been reached.
n = 15 follows from Tao's large-degree theorem.
The Lean theorem proves EHP114.
The radial Puiseux theorem proves global maximality.
The shadow dynamic is a universal law.
The shortcut n = 15/n = 16 JSONs are public proof certificates.
```

The present claim ceiling is: radial singularity identified, Lean reduction
compiled from explicit axioms, interval-hardening target now precise.

## Immediate handoff target

Create the next artifact, preferably under the same local/internal scope:

```text
EXP-MATH-EHP114-N14-PUISEUX-INTERVAL-HARDENING-20260505-01_RESULTS.json
EXP-MATH-EHP114-N14-PUISEUX-INTERVAL-HARDENING-20260505-01_REPORT.md
EXP-MATH-EHP114-N14-PUISEUX-INTERVAL-HARDENING-20260505-01_RESULTS.sha256
```

Minimum fields:

```json
{
  "problem": "Erdos #114 / EHP lemniscate perimeter",
  "degree": 14,
  "target": "gauss_2F1_puiseux_lower_bound_n14_interval",
  "parameters": {"a": "13/28", "b": "13/28", "c": "1"},
  "domain_w": "(0, 1/2]",
  "split": {"small_w": "...", "middle_w": "..."},
  "lower_bound_K14": "...",
  "method": "connection_formula_plus_interval_arithmetic",
  "claim_ceiling": "fixed radial boundary certificate only; no global EHP114 proof"
}
```

No scratch Lean was created for this packet because the useful reduction is
mathematical and documentary: the current compiled file already shows the exact
generic axiom to replace.
