# Prima Facie Next-Target Scan

Date: 2026-05-01  
Status: local fallback scan; not an attack-queue rerank

## What Was Queried

The documented attack-queue script path
`Math-Problems/erdos-mdl-mapper/scripts/build_morphism_index.py` was not present
in this checkout. A repository search found `erdos_problems.db` and older
`Math-Problems` tooling, but not the canonical mapper script named in the local
policy.

Fallback used: lexical scan over
`erdosatlas-workbench/public/data/workbench/library/erdos_*.md` for terms tied
to the compiled recipe: Sidon, additive, representation, extremal, face,
field, entropy, and constraint.

## Best Prima Facie Transfers

| Rank | Problem | Why it matches the recipe | Action |
|---:|---|---|---|
| 1 | #153 Sidon set sumset gap variance | Same Sidon universe; turns the object from `A` to a derived face-like structure `A+A` and its gap observables. | Best next mathematical target after PT-1. |
| 2 | #1 distinct subset sums | Additive representation uniqueness; finite families and observable spread are natural. | Good later target, but not before the Sidon lane is tightened. |
| 3 | #712 hypergraph Turan density | Extremal faces and exposed structures are central, but the machinery is much heavier. | Watchlist only. |
| 4 | #128 triangle-free sparse halves | Extremal graph constructions and stability faces. | Watchlist only. |
| 5 | #713 rational Turan exponents | Extremal graph-growth faces and algebraic construction classes. | Watchlist only. |
| 6 | #120 Erdos similarity | Avoidance/constraint problem with density fields, but far from current finite Sidon machinery. | Do not climb yet. |

## Conclusion

Do not broaden the Collider claim. The next climb should stay on the Sidon
ridge:

1. Keep the compiled n=61..64 selected-surface and selected-surface transition
   certificates as the finite foothold.
2. Turn the certified transition partition into a symbolic target involving
   `+1` persistence, new selected witnesses, and tiny selected surfaces.
3. Define one derived observable on `A+A` for #153, probably a gap-variance or
   local-spacing functional.
4. Test whether the Sidon extremal face that splits under prefix/mass fields
   also separates under the #153 sumset-gap observable.

That would be a real ascent: from abstract exposed-face lemma, to exact Sidon
witness certificate, to a neighboring Erdős problem with the same mathematical
object family.
