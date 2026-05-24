# EXP-MATH-ERDOS30-BFR-REPRESENTATION-CANONICAL-PORT-20260508-01

## Meaning

Eratosthenes moved the #30 representation-function bridge from the isolated scout file into the canonical BFR Lean file. This is the right kind of progress for the current lane: not a new bound, not a public result, but a reusable formal block in the file where the BFR closure work will happen.

The landed bridge says:

- `bfrRepFunction A s` counts ordered representations `a + b = s`.
- Under the Sidon condition, each fiber has at most the two swap-related ordered witnesses.
- Over the natural support range `{0, ..., 2N}`, the sum of all representation fibers is `|A|^2`.

That is the bookkeeping needed before interval occupancy can be made precise.

## Verification

- Command: `lake env lean lean/Erdos30_BFR.lean`
- Result: exit code `0`
- Canonical file SHA-256 after port: `4ba292e877fcf436afaa37ab8ccd1217a850e0d32ba7df572dc1fc5c4ce7d921`
- Source scout SHA-256: `eaddd7c50aa9849c4f94e7c9c52ecf80482d7b81ab4c3b4a570d639f708a1af4`

## Boundary

The BFR core bound remains a review target. This port only adds representation mechanics and support-range counting. It does not move the public Sidon coefficient, update D1, update curated morphisms, or touch #114 artifacts.

## Next Move

The next useful block is interval occupancy:

1. Define finite intervals inside `bfrSumRange N`.
2. Define `y_i` as the sum of `bfrRepFunction A s` over each interval.
3. Connect total occupancy to `bfr_sum_repFunction_eq_card_sq`.
4. Only after that, approach the discrepancy inequality.
