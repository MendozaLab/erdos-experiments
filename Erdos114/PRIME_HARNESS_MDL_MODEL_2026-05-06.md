# Prime-Harness MDL Model

Date: 2026-05-06
Status: internal research model

## Claim Ceiling

This note defines a finite encoding-cost diagnostic for Beurling-Nyman MDL
experiments. It is not a Riemann Hypothesis claim, not a zeta formalization, and
not a quantum-mechanical claim.

## Meaning

The failed composite-first projection run suggests a better finite question.
Primes may not behave like later residual columns after composites. They may be
the addressing harness that makes composite-indexed structure decodable.

Safe working hypothesis:

```text
prime-indexed structure acts as an arithmetic decoding basis;
composite-indexed columns are cheap only when addressed through that basis.
```

This is an MDL question about finite encodings. The measurable object is:

```text
residual tolerance reached per description bit under different index codes.
```

## Encodings

For dictionary indices `j <= N`, compare:

1. `flat_index`: address each selected `j` directly using `ceil(log2(N + 1))`
   bits.
2. `composite_index`: address only `j = 1` and composite `j > 1` by rank inside
   the composite bucket. Prime indices are not representable in this code.
3. `factorized_reusable_harness`: address `j` by its prime factorization,
   assuming the prime harness up to `N` is already shared by the dictionary.
4. `factorized_with_harness`: same factorization code, but charges once for the
   prime harness itself.

The two factorized versions are intentionally both reported. If the reusable
version wins but the charged version loses, the conclusion is not that primes
are free. It is that a reusable prime-address layer may be doing useful MDL
work only after its setup cost is amortized.

## Projection Supports

For each dictionary, grid, and size, the diagnostic fits:

- integer-prefix supports `1..m`;
- composite-prefix supports `1` plus the first composite indices;
- factor-cost-ordered supports, selecting the cheapest factorization addresses.

The same support can be scored under multiple encodings. This separates a pure
address-code comparison from a projection-subspace comparison.

## Metrics

Each row reports:

- residual norm and certified upper residual after the Lean quantization
  penalty;
- active coefficient count;
- address bits;
- coefficient bits;
- total description bits;
- information bits, `-log2(certified_upper_relative_residual)`;
- description bits per information bit.

Tolerance tables then ask which encoding reaches fixed residual thresholds with
the fewest description bits.

## Run 01 Result

Immutable experiment:

```text
EXP-MATH-RH-BN-PRIME-HARNESS-MDL-20260506-01
```

The first strict finite diagnostic did not find a broad prime-harness MDL
advantage.

```text
status: NO_PRIME_HARNESS_MDL_SIGNAL
same-support comparable rows: 144
factorized reusable wins: 12 / 144
factorized with-harness wins: 0 / 144
median reusable savings vs flat: -28.0 bits
median charged savings vs flat: -91.5 bits
```

Tolerance winners were mixed:

```text
flat_index:                    51
composite_index:               47
factorized_reusable_harness:   20
factorized_with_harness:        0
```

The interpretation is negative but useful. A naive factorization address code
is usually more expensive than flat finite indexing in this small regime. The
only remaining signal is local: factorized reusable addressing wins some
tolerance rows when the support is selected by factorization cost. That is not
enough for a paper claim. The next version should test amortized harness cost
across many targets or tasks, where a prime-address layer can be reused rather
than charged inside a single approximation row.

## Amortized Run 01 Result

Immutable experiment:

```text
EXP-MATH-RH-BN-PRIME-HARNESS-AMORTIZED-20260506-01
```

The amortized bundle test also did not produce a broad positive signal.

```text
status: NO_AMORTIZED_PRIME_HARNESS_SIGNAL
bundles: 18
shared-harness wins: 8 / 18
median shared-harness savings vs flat: -80.5 bits
median reusable savings before setup: -21.0 bits
```

There are real local wins, especially for the geometric dictionary and some
`N = 48` bundles, but the median bundle still loses. This points to the next
modeling correction: homogeneous bit cost may be the wrong cost surface. If the
N-body analogy is load-bearing, the comparison should charge direct indexing,
prime lookup, multiplication, exponentiation, and coefficient movement as
different compute channels rather than one undifferentiated bit budget.

## Safe Language

Safe:

- prime-harness MDL diagnostic;
- arithmetic addressing basis;
- factorization-coded finite dictionary;
- reusable harness versus charged harness;
- finite accessible-description experiment.

Unsafe:

- identifying primes with a physical measurement theorem;
- claiming a quantum result;
- claiming an RH consequence;
- treating the finite encoding result as dictionary invariant.
