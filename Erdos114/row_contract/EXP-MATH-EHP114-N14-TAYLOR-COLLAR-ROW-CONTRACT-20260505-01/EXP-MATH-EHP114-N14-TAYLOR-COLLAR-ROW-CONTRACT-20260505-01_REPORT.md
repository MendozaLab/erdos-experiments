# EHP114 n=14 Taylor-Collar Row Contract

Experiment: `EXP-MATH-EHP114-N14-TAYLOR-COLLAR-ROW-CONTRACT-20260505-01`

## Verdict

- Status: `ROW_CONTRACT_ACCEPTED_NOT_THEOREM`
- Source result SHA-256: `accc2124a77606187ec73fd000d98da517fbb3adc3a1b9ae5ffea5c45332356b`
- Row count: `15793`
- Row failure count: `0`
- Normal-drift error sum: `2.4736484172357147`
- Normal-drift budget: `2.5620009612530126`
- Normal-drift slack: `0.08835254401729786`
- Normal error / budget: `0.9655142424403749`
- Relative error / budget: `2.730675196761203`
- Row contract SHA-256: `8417f5903d3696fe300961359e65e4771c1c94269d067e4dd81db7341c422e47`

## Claim Ceiling

This is a Taylor-collar normal-drift row contract only. It is not an exact
lemniscate-length certificate, not a Lean theorem, and not a proof of Erdős
#114.

## Lean-Shaped Target

```text
For every retained z16 collar row B in the (6,4) root-affine parameter cell, the local 6x6 sub-bound encloses |p|, |p'|, and |p''|; Taylor gives |p'| > 0 on B; and the finite sum of normal-drift row errors is <= 2.5620009612530126.
```

## Next Blocker

Turn this row contract into a theorem-grade checker or Lean-side finite-sum certificate, then justify that the normal-drift bridge is the correct exact-length comparison theorem.
