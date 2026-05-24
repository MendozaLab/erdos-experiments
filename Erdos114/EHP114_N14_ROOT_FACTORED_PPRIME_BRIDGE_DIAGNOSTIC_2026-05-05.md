# EHP114 n=14 Root-Factored `p'` Bridge Diagnostic

Date: 2026-05-05

## Meaning

The cheapest proposed bridge fix was to evaluate the derivative in root-factored
form,

```text
p'(z) = sum_i prod_{j != i} (z - r_j),
```

instead of using the recurrent derivative carried through the polynomial
product. This tested whether the bridge failure was mainly an interval
dependency artifact in `p'`.

The result is negative but useful: root-factored `p'` does not improve the
current interval bridge on the known hard `n=14`, `eps=0.1`, worst subcell
`(6,4)`. It worsens the unresolved regularity count.

This remains a shadow signature, not universal law. It is not an exact
lemniscate-length certificate, not a Lean theorem, and not a proof of Erdős
#114.

## Verified Local Artifact

```text
engine = erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_bridge_diagnostic.rs
output = erdos-experiments/scripts/erdos-114/bridge-diagnostic-root-factored-pprime-output-z16/
experiment_id = EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01
derivative_mode = compare
z_subdivision = 16
status = BRIDGE_REGULARITY_INTERVAL_UNRESOLVED
```

Hash verification:

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/bridge-diagnostic-root-factored-pprime-output-z16
shasum -a 256 -c *_RESULTS.sha256
```

Output:

```text
EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01_RESULTS.json: OK
```

## Result

```text
recurrent unresolved cells     = 13560
root-factored unresolved cells = 17569
regularity resolution delta    = -4009
normal error / budget          = 200.17607112416812
relative error / budget        = 566.1821693511315
minimum gradient lower         = 0.0
```

The root-factored derivative is therefore not the next bridge route. It
expands the interval enclosure enough to make regularity harder, not easier.

## Next Proof Move

Do not spend more effort on blind z-subdivision or root-factored `p'` alone.
The next serious bridge attempt should be one of:

1. Bernstein or affine/Taylor model arithmetic for `p'` on the active
   level-set boxes.
2. A validated level-set collar that avoids checking empty ambient space.
3. Critical-point exclusion for `p'`, if the critical-point isolator can be made
   cheaper than the Bernstein/Taylor path.

The immediate recommended successor is Bernstein or affine/Taylor arithmetic on
the same worst subcell before scaling to the 164-cell packet.

