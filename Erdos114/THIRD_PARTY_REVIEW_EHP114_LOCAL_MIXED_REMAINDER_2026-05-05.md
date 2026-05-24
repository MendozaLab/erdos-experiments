# Third-Party Review — EHP114 Local Mixed-Remainder Target

Date: 2026-05-05  
Reviewer: Perplexity via `/third-party-review`  
Scope: Internal theorem-target review only

## Query

Review the fixed-n `n=14` Erdős #114 local-stability theorem target:

```text
D14(eps,s) = R14(eps) + Q14(s) - M14(eps,s)
R14(eps) >= 24 eps^(1/14)
Q14(s)   >= 100000 ||s||^2
target: M14(eps,s) <= 12 eps^(1/14) + 50000 ||s||^2
```

The target domain is:

```text
0 < eps <= 1/10
RootsInClosedUnitDisk(radialMode14(eps)+shapeMode14(s))
||s|| <= min(eta14Boundary(eps), 0.008)
```

The Lean scratch scaffold compiles and proves only the algebraic splice under
three named assumptions: `RadialCertificate14`, `ShapeConeCertificate14`, and
`MixedRemainderAbsorption14`.

## Perplexity Assessment

Perplexity's overall stance:

> The theorem-target decomposition is mathematically sane and appropriate for a
> rigorous local stability proof at fixed n=14 around `p(z)=z^14-1`, leveraging
> the correct radial singular scale `eps^(1/14)`.

Main points:

1. The decomposition
   `D14 = R14 + Q14 - M14` is a clean split between radial Puiseux deficit,
   shape-cone energy, and mixed remainder.

2. The scratch Lean scaffold is correctly scoped: it formalizes the algebraic
   consequence without pretending to prove the analytic estimates.

3. The radial exponent `eps^(1/14)` is plausible as the correct boundary scale
   for the radial family near `z^14 - 1`.

4. The local cone condition is mathematically reasonable but still ad hoc until
   `eta14Boundary(eps)` is turned into an analytic or interval-certified
   topology/locality statement.

5. The weakest link is the mixed-remainder absorption theorem. The finite scout
   is evidence, not a replacement for a proof.

6. Perplexity ranked the routes as:
   - best: validated interval Taylor model over `(eps, s)` boxes;
   - second: Cauchy/local derivative bounds;
   - third: Jordan/spectral cone normal form, useful for articulation but risky
     as the proof's critical path.

Perplexity bottom line:

> Pursue route B for `M14` with topology-augmented cone scouting, targeting a
> Lean-formalized local theorem; ceiling at cone-confined stability to avoid
> global overclaim.

## Codex Caution

The review supports the direction, but two parts should not be imported as
authority:

1. The cited/source mapping in the Perplexity response is uneven. Treat it as a
   reviewer sanity check, not as a literature packet.

2. Its wording around the radial asymptotic mixes heuristic language too
   loosely. Our internal claim should remain the artifact-backed fixed-n
   statement:

```text
For n=14, the radial family has an interval-backed Puiseux lower-bound target
with exponent 1/14 and conservative constant 24.
```

Do not upgrade that to a general theorem without a separate citation and proof
audit.

## Decision

The third-party review agrees with the current next move:

```text
build an interval Taylor / validated remainder packet for M14 on the local cone.
```

The next artifact should not be another broad scout. It should be a box-wise
Taylor certificate that tries to bound:

```text
M14(eps,s) <= 12 eps^(1/14) + 50000 ||s||^2
```

under:

```text
0 < eps <= 1/10
||s|| <= min(eta14Boundary(eps), 0.008)
```

## Claim Ceiling

Safe language:

```text
An external reviewer agrees that the n=14 local-stability packet has isolated
the right bottleneck: mixed-remainder absorption on a local admissible cone.
```

Unsafe language, paraphrased:

- Do not claim the conjecture is settled.
- Do not claim the local cone proves global maximality.
- Do not cite Perplexity as mathematical authority.
- Do not promote the finite scout to theorem status.

## Perplexity Source Links

- Tao arXiv source cited by Perplexity: https://arxiv.org/abs/2512.12455
- Tao blog source cited by Perplexity: https://terrytao.wordpress.com/2025/12/15/the-maximal-length-of-the-erdos-herzog-piranian-lemniscate-in-high-degree/
- Fryntov/Nazarov arXiv source cited by Perplexity: https://arxiv.org/pdf/0808.0717
- LMS source cited by Perplexity: https://londmathsoc.onlinelibrary.wiley.com/doi/abs/10.1112/plms/pdw039
- Erdős Problems thread cited by Perplexity: https://www.erdosproblems.com/forum/thread/1045
- Huhtanen source cited by Perplexity: https://math.aalto.fi/~mhuhtane/indef.pdf
- ADS mirror cited by Perplexity: https://ui.adsabs.harvard.edu/abs/2008arXiv0808.0717F/abstract
- PMC source cited by Perplexity: https://pmc.ncbi.nlm.nih.gov/articles/PMC11872773/

