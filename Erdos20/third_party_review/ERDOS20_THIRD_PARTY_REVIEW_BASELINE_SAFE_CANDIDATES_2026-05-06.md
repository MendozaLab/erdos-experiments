# Erdős #20 Third-Party Review: Baseline + SafeCandidates Bridge

Date: 2026-05-06

Reviewer: Perplexity, via `/3rd` default mode

Scope: independent critique of the next #20 bridge run. This is not a proof of the sunflower conjecture, not a new lower bound, and not a scorecard upgrade.

Claim ceiling: A0. This remains a shadow signature, not universal law.

## Query

We asked whether the next #20 step should use floor-normalized local core-closure information,

```text
primary_quantity = I_core_local(s,m) / I_floor(s,m)
secondary_control = Abbott-Hansen-Sauer-style construction-normalized closure cost
```

and whether the enumerator should emit candidate-level rows that bind

```text
safeCandidates = candidateExtensions.filter (SafeExtension core family)
```

to real data.

## Perplexity Assessment

Perplexity supported the proposed bridge as a reasonable next diagnostic, but not as theorem progress. It agreed that `I_core_local / I_floor` is defensible as a Leg-4 information-floor quantity, because it asks whether the observed core-closure channel is large relative to a predeclared information denominator rather than merely measuring raw local cost.

The main critique was calibration. Perplexity argued that Abbott-Hansen-Sauer-style normalization should be promoted from secondary control to co-primary reporting, because it ties the experiment back to known construction/literature baselines. Raw `I_core_local` remains useful as the numerator, but is not enough by itself.

Minimum evidence threshold proposed by Perplexity:

- candidate-level `safeCandidates` rows, not only summary counts;
- reproduction or calibration against known sunflower constructions/bounds;
- correlation between the selected ratio and sunflower emergence in benchmark cases;
- a Lean-checked `SafeExtension` predicate before any theorem-language upgrade.

## Sources Cited By Perplexity

- https://en.wikipedia.org/wiki/Sunflower_(mathematics)
- https://homepages.cwi.nl/~lex/files/ESsunflowerconjectureCarla.pdf
- https://arxiv.org/abs/2505.03671
- https://gilkalai.wordpress.com/2016/05/17/polymath-10-emergency-post-5-the-erdos-szemeredi-sunflower-conjecture-is-now-proven/
- https://www.renyi.hu/~pach/publications/OddSunflowerReprint032124.pdf
- https://www.quantamagazine.org/mathematicians-begin-to-tame-wild-sunflower-problem-20191021/
- https://www.erdosproblems.com/forum/thread/856
- https://github.com/google-deepmind/formal-conjectures/issues/2284

## Actionable Update To The #20 Plan

Run the bridge as planned, but report two primary columns side-by-side:

```text
floor_ratio = I_core_local(s,m) / I_floor(s,m)
ahs_ratio = I_core_local(s,m) / I_AHS(s,m)
```

Treat `floor_ratio` as the Maxwell/Mendoza Leg-4 quantity and `ahs_ratio` as the literature/construction calibration. A run that only reports `I_core_local` remains precursor evidence.

## Claim Ceiling

The third-party review strengthens the design discipline but does not improve the mathematical status. The next #20 artifact should say:

> The bridge run binds candidate-level safe-extension data to floor and construction-normalized closure ratios. Status remains A0: shadow signature, not universal law.
