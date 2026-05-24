# EHP114 n=14 Radial Puiseux Interval Target

Experiment: `EXP-MATH-EHP114-N14-RADIAL-PUISEUX-INTERVAL-TARGET-20260505-01`

## Meaning

This packet turns the n=14 stress case into a narrow proof target. The point is
not to rerun the expensive branch-and-bound certificate. The point is to carve
out the first analytic bridge between the DOI-backed finite proof and Tao's
large-n proof: a fixed-degree radial Puiseux lower bound that can be interval
hardened.

The claim ceiling is unchanged: this is a shadow signature, not universal law.
It is not a proof of #114, not a global local-stability theorem, and not a
scorecard upgrade.

## Source Anchor

- DOI: `10.5281/zenodo.19480329`
- Finite certificate: `EXP-MM-EHP-007-n14-inari`
- Verdict: `EHP_N14_PROVEN`
- Rigor: `ieee_1788_interval_arithmetic_inari`
- Reduced dimension: `25`
- Interval evaluations: `855638016`
- B&B complete: `True`
- Exact L* inside inari interval: `True`

n=14 is the materially stringent finite anchor because it is the largest
byte-reconciled DOI-backed Rust/inari certificate in the current corpus. The
n=15/n=16 records remain useful hints, but their zero-evaluation provenance
needs reconciliation before they can carry the same weight.

## Candidate Constants

| radial parameterization | sampled lower evidence | working constant |
|---|---:|---:|
| direct `a=1-eps` | min ratio 26.2084041204575507279830524281 at eps 0.1 | C14 = 24.0 |
| direct tail `eps<=1e-4` | min ratio 26.5870905812431276381831686628 at eps 0.0001 | C14_tail = 26.0 |
| calibration schedule `radius=1-eps/sqrt(14)` | min ratio 28.8151359533827510919383170269 at eps 0.02 | C14 = 24.0 |
| calibration tail `eps<=1e-4` | min ratio 29.2112033305747017067704054676 at eps 0.0001 | C14_tail = 27.0 |

These are floating diagnostics. The usable next move is to replace them with
interval-certified constants. The deliberate target constant is conservative:
prove 24 first, then sharpen only if that proof is clean.

## Interval-Hardening Split

1. Singular tail: for `0 < eps <= 1e-4`, use the Gauss 2F1 connection formula
   at `z=1` to certify the Puiseux coefficient.
2. Compact middle: for `1e-4 <= eps <= 1e-1`, interval-subdivide the exact
   hypergeometric or integral formula and certify the weaker `C14=24` bound.
3. Radial-to-shape splice: feed this radial deficit into the shape cone and
   mixed-remainder packets. Do not infer the global theorem from the radial
   slice alone.

## Lean-Shaped Target

```lean
theorem ehp114_n14_radial_puiseux_interval_direct
    (eps : Real) (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10000) :
    (24 : Real) * Real.rpow eps ((1 : Real) / 14)
      <= D14 (radialMode14 eps) := by
  -- fixed-n Gauss 2F1 connection formula + interval constants
  -- no global EHP conclusion follows from this lemma alone
  sorry
```

## Single Bottleneck

Interval-hardening the radial hypergeometric/Puiseux singularity model.

The next agent should attack exactly this theorem target:
`Prove ehp114_n14_radial_puiseux_interval_direct with an interval-certified 2F1 connection-formula lower bound.`.

## Guardrails

No scorecard, D1, public document, Lean file, git, or email state was changed.
