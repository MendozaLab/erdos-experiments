# EXP-MATH-ERDOS30-SYM2-FORMAL-CARRIER-SPEC-20260508-01

Status: `REVIEW_ONLY / SYM2_FORMAL_CARRIER_SPEC_READY`

## Decision

`FORMALIZE_CARRIER_DISCIPLINE_BEFORE_OPERATOR_BOUND`

The Sym2 result teaches where the operator should live. It separates unordered-pair structure from ordered-pair swap noise. That is useful, but it is carrier discipline only; it does not improve the #30 coefficient.

## Formal Scout Target

`Erdos30_Sym2PairSumOperator_SCOUT.lean`

The target should define the unordered carrier, off-diagonal carrier, unordered sum map, and a pair-sum operator matrix over the correct quotient. The first useful theorem is a trace/counting identity, not an asymptotic bound.

## Claim Ceiling

SCOUT only. No coefficient pressure and no public #30 claim.
