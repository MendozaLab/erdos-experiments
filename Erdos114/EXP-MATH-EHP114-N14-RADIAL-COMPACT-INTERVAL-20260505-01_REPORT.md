# EHP114 n=14 Radial Compact Interval Certificate

Experiment: `EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01`

## Meaning

This is the first executable interval-hardening result from the n=14 radial/Puiseux lane. It certifies the compact-middle radial interval `1e-4 <= eps <= 1e-1` for the one-parameter family `p_a(z)=z^14-a` with the conservative constant `C14=24`.

The claim ceiling is still narrow: this is a shadow signature, not universal law. It certifies one radial-family compact interval only. It does not prove #114, does not close local stability, and does not touch the global nonradial cone.

## Verdict

- Status: `COMPACT_MIDDLE_CERTIFIED`
- Constant: `C14 = 24.0`
- Domain certified: `[1e-4, 1e-1]`
- Missing domain: `(0, 1e-4] singular tail; requires hypergeometric connection formula`
- Minimum margin: `0.4683310451937004`
- Runtime: `4.927s`

## Method

The exact radial formula is rewritten as

```text
L_14(1-eps) = 2 int_0^pi (eps^2 + 2(1-eps)(1-cos u))^(-13/28) du.
```

For each epsilon bin, the verifier computes an IEEE-1788 interval upper Riemann enclosure of `L_14(1-eps)`, then checks

```text
L*_lower - upper(L_14(1-eps)) >= 24 * eps_hi^(1/14).
```

This is intentionally one-sided: an upper bound on the radial length gives a lower bound on the radial deficit.

## Certified Bins

| eps bin | u steps | deficit lower | RHS upper | margin | pass |
|---|---:|---:|---:|---:|---:|
| [0.0001, 0.0002] | 250000 | 13.7048916009 | 13.0616817695 | 0.643209831406 | True |
| [0.0002, 0.0005] | 160000 | 14.4134873 | 13.9451562548 | 0.468331045194 | True |
| [0.0005, 0.001] | 120000 | 15.4122007721 | 14.6529655118 | 0.759235260251 | True |
| [0.001, 0.002] | 90000 | 16.1996391416 | 15.3967007875 | 0.802938354136 | True |
| [0.002, 0.005] | 70000 | 17.0138646885 | 16.4381128005 | 0.575751888029 | True |
| [0.005, 0.01] | 50000 | 18.149933009 | 17.272456152 | 0.877476857004 | True |
| [0.01, 0.02] | 40000 | 19.0333460088 | 18.1491479676 | 0.884198041218 | True |
| [0.02, 0.05] | 30000 | 19.8817097492 | 19.3767317844 | 0.504977964773 | True |
| [0.05, 0.1] | 24000 | 21.0802264557 | 20.3602295579 | 0.719996897882 | True |

## What Remains

The compact middle is no longer the blocker. The next blocker is exactly:

`Prove the singular tail 0 < eps <= 1e-4 by Gauss 2F1/Puiseux connection formula.`

That means the next theorem target is the singular-tail version of the same bound, not a bigger brute-force run.

## Files

- Rust verifier: `erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_radial_compact_interval.rs`
- Result JSON: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01_RESULTS.json`
- SHA sidecar: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01_RESULTS.sha256`

No scorecard, D1, public document, Lean file, git, or email state was changed.
