# EHP114 n=14 Taylor-Collar Local-Bound Bridge Candidate

Date: 2026-05-05

## Meaning

The prior Taylor-collar diagnostic changed the bridge problem from
"regularity unresolved" to "constants too loose." The next constant-tightening
test subdivided each retained Taylor-collar box internally when bounding
`|p|`, `|p'|`, and `|p''|`.

The `6 x 6` local-bound run crosses the normal-drift budget on the known hard
`n=14`, `eps=0.1`, worst subcell `(6,4)`:

```text
regularity unresolved cells = 0
normal error / budget       = 0.9655142424403749
relative error / budget     = 2.730675196761203
```

The row-contract extractor has now materialized the retained-box obligations
as a machine-readable artifact:

```text
row_contract = erdos-experiments/Erdos114/row_contract/EXP-MATH-EHP114-N14-TAYLOR-COLLAR-ROW-CONTRACT-20260505-01/
status = ROW_CONTRACT_ACCEPTED_NOT_THEOREM
row_count = 15793
row_failure_count = 0
normal_drift_slack = 0.08835254401729786
row_contract_sha256 = 8417f5903d3696fe300961359e65e4771c1c94269d067e4dd81db7341c422e47
```

This is the first local bridge candidate under the normal-drift budget. It is
not yet an exact lemniscate-length certificate, because the relative-length
candidate remains over budget and the normal-drift bridge still needs to be
justified as the correct exact-length comparison theorem.

This remains a shadow signature, not universal law. It is not a Lean theorem
and not a proof of Erdős #114.

## Verified Local Artifact

```text
engine = erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_bridge_diagnostic.rs
output = erdos-experiments/scripts/erdos-114/bridge-diagnostic-taylor-collar-local-bound6-output-z16/
experiment_id = EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01
degree = 14
eps = 0.1
subcell = (6,4)
z_subdivision = 16
derivative_mode = recurrent
regularity_strategy = taylor_collar
local_bound_subdivision = 6
status = BRIDGE_CANDIDATE_ERROR_UNDER_BUDGET_NOT_THEOREM
```

Hash verification:

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/bridge-diagnostic-taylor-collar-local-bound6-output-z16
shasum -a 256 -c *_RESULTS.sha256
```

Output:

```text
EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01_RESULTS.json: OK
```

## Progression

```text
ambient z16:
  regularity unresolved = 13560
  normal error / budget = 58.40007122458721

taylor-collar z16:
  regularity unresolved = 0
  normal error / budget = 2.1882768069742764

taylor-collar z16 + local-bound 4:
  regularity unresolved = 0
  normal error / budget = 1.0531457929156913

taylor-collar z16 + local-bound 5:
  regularity unresolved = 0
  normal error / budget = 0.9993520731171657

taylor-collar z16 + local-bound 6:
  regularity unresolved = 0
  normal error / budget = 0.9655142424403749
```

## Next Proof Move

The next target is no longer "find regularity." It is to formalize the
Taylor-collar local-bound row contract:

```text
For every retained z16 collar box B in the (6,4) root-affine parameter cell,
the local 6 x 6 sub-bound encloses |p|, |p'|, and |p''|;
the Taylor lower bound gives |p'| > 0 on B;
the sum of normal-drift errors is <= 2.5620009612530126.
```

If this row contract is accepted, the successor task is to either tighten the
relative-length candidate or prove that the normal-drift bridge is the correct
comparison theorem for exact length.
