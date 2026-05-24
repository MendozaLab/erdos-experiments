# EHP114 -> RH MDL Prior-Art Pass

Date: 2026-05-06
Method: Perplexity prior-art query plus Codex verification against primary
sources
Scope: internal prior-art memo for `EHP114_RH_MDL_OVERLOOK_2026-05-06.md`
Status: PRIOR_ART / CLAIM_CEILING, not theorem progress

## Meaning

The Perplexity pass supports the current claim ceiling with one important
qualification. The EHP and Nyman-Beurling sides are heavily covered prior art.
The apparent gap is narrower:

> An explicit MDL or description-length cost curve for the Nyman-Beurling /
> Baez-Duarte approximation problem appears not to be established in the
> peer-reviewed literature checked here.

However, "dictionary approximation" and "compressed sensing" language is not
fully untouched: a recent preprint on Preprints.org explicitly frames the
Baez-Duarte problem as dictionary approximation / compressed sensing. Treat that
as gray-literature risk, not as a reliable theorem source until checked more
deeply.

## Disagreement With Perplexity Output

DISAGREEMENT: Perplexity mixed correct prior-art direction with unreliable
details. It described Tao's EHP result as unpublished/blog-like, but arXiv
currently lists Terence Tao's `arXiv:2512.12455`, submitted 2025-12-13 and
revised 2025-12-22. It also included unrelated citation links in its final
source list, including ML/Stiefel links unrelated to EHP, RH, or
Nyman-Beurling. I ignored those and verified the load-bearing claims against
primary pages.

## Exact Prior Art That Must Be Cited

### EHP Lemniscate Extremality

1. **Fryntov-Nazarov, 2008.** `arXiv:0808.0717`, "New estimates for the length
   of the Erdos-Herzog-Piranian lemniscate."
   - Exact relevance: monic degree-`n` polynomial lemniscates, local maximum at
     `z^n - 1`, and asymptotically sharp bound `|L| < 2n + o(n)`.
   - Claim impact: preempts any broad novelty claim around EHP arclength
     extremality itself.

2. **Tao, 2025.** `arXiv:2512.12455`, "The maximal length of the
   Erdos-Herzog-Piranian lemniscate in high degree."
   - Exact relevance: proves the EHP conjecture for all sufficiently large
     `n`, building on Fryntov-Nazarov.
   - Claim impact: our EHP side can only be a calibration theorem / small-n /
     method-translation lane, not a new high-degree theorem.

### Nyman-Beurling / Baez-Duarte RH Approximation

3. **Nyman, 1950** and **Beurling, 1955.**
   - Exact relevance: original closure criterion linking approximation by
     fractional-part dilates to zero-free regions / RH.
   - Verified Beurling page: PNAS 41, 312-314, DOI
     `10.1073/pnas.41.5.312`.

4. **Baez-Duarte, 2002/2003.** `arXiv:math/0202141`, "A strengthening of the
   Nyman-Beurling criterion for the Riemann Hypothesis."
   - Exact relevance: RH equivalent to approximating the indicator by the
     smaller integer-dilate subspace.
   - Claim impact: preempts any claim that integer dictionaries are a new RH
     bridge.

5. **Baez-Duarte, 2005.** `arXiv:math/0505453`, "A general strong
   Nyman-Beurling Criterion for the Riemann Hypothesis."
   - Exact relevance: general strong criterion using co-Poisson / Muntz
     transforms and dilation closures.
   - Claim impact: any broader kernel/dilation formulation must cite this.

6. **Burnol, 2001/2002.** `arXiv:math/0103058`, "A lower bound in an
   approximation problem involving the zeros of the Riemann zeta function."
   - Exact relevance: lower bounds in the Nyman-Beurling approximation problem
     involving zeta zeros.
   - Claim impact: this is the main preexisting quantitative-rate anchor. Any
     MDL curve must be positioned as a cost refinement, not as the first
     approximation-rate question.

7. **Alouges-Darses-Hillion, 2022.** "Polynomial approximations in a
   generalized Nyman-Beurling criterion", JTNB 34 (2022), 767-785.
   - Exact relevance: Nyman-Beurling is an RH-equivalent approximation problem;
     randomized `theta_k` structures; coefficient-control condition; Gram
     matrices; block-Hankel simplification.
   - Claim impact: strongest immediate prior art for our proposed experiment's
     dictionary choices, coefficient control, and Gram/Hankel diagnostics.

8. **Darses-Hillion, 2021.** "On probabilistic generalizations of the
   Nyman-Beurling criterion for the zeta function", Confluentes Mathematici 13
   (2021), 43-59.
   - Exact relevance: probabilistic Nyman-Beurling variants and coefficient
     control.
   - Claim impact: preempts a loose claim that randomized dictionaries are new.

### Li / Positivity Side

9. **Lagarias, 2004/2007.** `arXiv:math/0404394`, "Li Coefficients for
   Automorphic L-Functions."
   - Exact relevance: Li-style positivity coefficients are tied to Riemann
     hypothesis criteria and Weil's quadratic functional.
   - Claim impact: "Li coefficients as slack" is at most an interpretation
     unless a new cost functional is proved.

## Gray-Literature / Preprint Risk

The web pass found a recent Preprints.org manuscript:

`Spectral and Analytic Structure of the Nyman-Beurling-Baez-Duarte
Approximation`, manuscript `202506.0772`.

Relevant sections explicitly state that the Baez-Duarte problem admits a
dictionary-approximation / compressed-sensing reformulation and discusses LASSO,
coherence, Gram kernels, Li coefficients, and speculative spectral links.

Risk assessment:

- It is not peer-reviewed prior art at the same level as Burnol or
  Alouges-Darses-Hillion.
- It does preempt casual language such as "nobody has viewed Nyman-Beurling as
  dictionary approximation."
- It does **not**, on this pass, appear to define a predeclared MDL /
  description-length cost curve `K(epsilon)` charging theta bits, coefficient
  bits, support size, and quantization robustness.

Safe response: cite or footnote it as gray literature if this becomes an
external manuscript, then make the novelty claim specifically about MDL cost
and theorem-shaped asymptotic implications.

## Apparent Novelty Gap

The defensible novelty gap is not:

```text
Nyman-Beurling is an approximation problem.
```

That is old.

It is not:

```text
Nyman-Beurling has coefficient-control, Gram, Hankel, or randomized variants.
```

That is also prior art.

The remaining viable gap is:

```text
Define and analyze a computable MDL cost curve K(epsilon) for
Nyman-Beurling/Baez-Duarte approximation, where K charges for support,
dictionary-parameter precision, coefficient precision, and quantization
robustness; then relate asymptotic regimes of K to RH-equivalent or
zero-free-region statements.
```

## Safest Next Theorem-Shaped Claim

Do not claim an RH equivalence first. Start with a definition plus an
unconditional stability theorem.

Suggested target:

```text
Definition.
Fix a dictionary Theta_N and a computable bit-cost C_N for support,
theta values, and coefficients. Define

K_N(epsilon) =
  min C_N(a, theta)
  subject to ||chi_(0,1] - sum_j a_j rho_theta_j||_2 <= epsilon.

Theorem target.
For the frozen dictionary and cost model, coefficient quantization at b bits
changes the L2 residual by at most an explicit function of b, the Gram
condition number, and the unquantized coefficient norm.
```

Why this is the right first claim:

- It is not an RH proof claim.
- It directly controls the artifact risk in numerical MDL curves.
- It turns the experiment into a mathematical object: residual decay only
  counts if it survives quantization and conditioning penalties.
- It cites the real prior art cleanly while leaving a genuine contribution.

## Bibliography Links

- Tao, 2025: https://arxiv.org/abs/2512.12455
- Fryntov-Nazarov, 2008: https://arxiv.org/abs/0808.0717
- Beurling, 1955 DOI page: https://doi.org/10.1073/pnas.41.5.312
- Baez-Duarte strengthening, 2002: https://arxiv.org/abs/math/0202141
- Baez-Duarte general strong criterion, 2005:
  https://arxiv.org/abs/math/0505453
- Burnol lower bound, 2001/2002: https://arxiv.org/abs/math/0103058
- Alouges-Darses-Hillion, 2022:
  https://jtnb.centre-mersenne.org/articles/10.5802/jtnb.1227/
- Darses-Hillion, 2021: https://numdam.org/articles/10.5802/cml.71/
- Lagarias on Li coefficients, 2004/2007:
  https://arxiv.org/abs/math/0404394
- Gray-literature compressed-sensing preprint:
  https://www.preprints.org/manuscript/202506.0772

## Bottom Line

The MDL-overlook survives Perplexity, but the novelty claim must be sharpened:

Safe:

> Prior art already makes Nyman-Beurling an RH-equivalent approximation
> problem with serious coefficient-control and Gram/Hankel structure. What
> appears open is a frozen, computable MDL cost curve with quantization and
> conditioning penalties, and a theorem relating that curve to known
> zero-free-region or RH-equivalent bounds.

Unsafe:

> EHP114 gives a new route to RH.

Unsafe:

> Nyman-Beurling has not been viewed as dictionary approximation.

Unsafe:

> Li coefficients are already established as MDL slack.
