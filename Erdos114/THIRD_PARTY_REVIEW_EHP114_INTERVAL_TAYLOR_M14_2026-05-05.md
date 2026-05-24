# Third-Party Review — EHP114 Interval Taylor / M14 Split Verdict

Date: 2026-05-05  
Reviewer: Perplexity via `/third-party-review`  
Scope: Internal interval Taylor packet review only

## Query

Review the fixed-n `n=14` interval Taylor / M14 packet after a split verdict:

```text
axis endpoint budget: passes
uniform positive Taylor matrix at radial bases: fails
spectral diagnosis: radial-base shape softening detected
```

The question was whether this weakens the local stability route or simply
changes the theorem target.

## Perplexity Assessment

Perplexity's overall stance:

> The interpretation is mathematically sane and correctly identifies that
> radial contraction induces shape softening. This does not kill the local
> stability route, but it forces a revised theorem targeting radial reserve
> dominance over curvature losses.

Main points:

1. The split verdict is meaningful. Passing admissible axis endpoints with a
   positive margin while failing the ambient matrix test falsifies the naive
   claim that the boundary shape cone remains uniformly positive after radial
   contraction.

2. The negative spectral values are not automatically fatal. They suggest the
   shape quadratic form softens around radially contracted bases, so the radial
   Puiseux reserve must absorb that softening.

3. The next theorem should not be:

```text
the shape cone remains positive after radial contraction
```

It should be:

```text
the radial reserve dominates radial-base shape softening on the local cone
```

4. Perplexity recommended shifting from a fixed local cap to an
   epsilon-dependent cone radius, with a candidate scaling:

```text
eta14Boundary(eps) proportional to eps^(1/28)
```

This scaling is plausible because a quadratic shape term at radius
`eps^(1/28)` is at the radial Puiseux scale `eps^(1/14)`.

5. Perplexity ranked the next routes as:

- first: direct Taylor model for the total deficit, with radial reserve and
  negative curvature handled in one inequality;
- second: restrict `eta14Boundary(eps)` so reserve dominates softening;
- third: boundary-normal coordinates, useful later but not the immediate
  critical path.

## Codex Caution

Use the recommendation, not the citations, as the durable value of this review.
Some cited links are likely noisy or not the exact source for the fixed-n claim.
The review should not be treated as literature authority.

The useful mathematical correction is internal and concrete:

```text
replace fixed cap ||s|| <= 0.008 with a cone radius tied to eps,
then prove the total deficit bound directly or by interval Taylor boxes.
```

The candidate scaling to test next is:

```text
||s|| <= min(eta0 * eps^(1/28), admissibility radius)
```

where `eta0` should be chosen by interval search, not guessed.

## Decision

The next packet should be:

```text
EXP-MATH-EHP114-N14-EPS-SCALED-CONE-M14-SEARCH-20260505-01
```

Goal:

```text
Find a conservative eta0 such that all tested box/axis directions satisfy
D14(eps,s) >= 12 eps^(1/14)
for ||s|| <= eta0 eps^(1/28).
```

Then the theorem target becomes:

```lean
theorem ehp114_n14_eps_scaled_cone_deficit
    (eps : Real) (s : ShapeQuotient14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hadm : RootsInClosedUnitDisk14 eps s)
    (hcone : quotientNorm s <= eta0 * Real.rpow eps ((1 : Real) / 28)) :
    (12 : Real) * Real.rpow eps ((1 : Real) / 14)
      <= totalDeficit14 eps s := by
  -- direct total-deficit interval Taylor certificate target
  sorry
```

## Claim Ceiling

Safe internal language:

```text
External review agrees that the interval Taylor failure should be interpreted
as radial-base shape softening, and recommends replacing fixed-radius cone
language with an epsilon-scaled cone theorem.
```

Unsafe language, paraphrased:

- Do not claim the conjecture is settled.
- Do not claim local stability is proved.
- Do not claim the epsilon-scaled cone theorem before running it.
- Do not cite Perplexity as mathematical authority.

## Perplexity Source Links

- Tao arXiv source cited by Perplexity: https://arxiv.org/abs/2512.12455
- Tao blog source cited by Perplexity: https://terrytao.wordpress.com/2025/12/15/the-maximal-length-of-the-erdos-herzog-piranian-lemniscate-in-high-degree/
- Tao arXiv HTML source cited by Perplexity: https://arxiv.org/html/2512.12455v1
- Erdős Problems source cited by Perplexity: https://www.erdosproblems.com/forum/thread/114
- Erdős Problems source cited by Perplexity: https://www.erdosproblems.com/1038
- arXiv source cited by Perplexity: https://arxiv.org/abs/2312.13673
- Erdős Problems source cited by Perplexity: https://www.erdosproblems.com/114

