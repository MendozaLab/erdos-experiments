# EHP114 Tao-Style Proof Technique Note

Date: 2026-05-06

Status: private explanatory sidecar. This note is not part of the EHP114 proof
package, not evidence, not a collaborator statement, and not a public claim
surface.

Sources skimmed:

- Terence Tao, "A proof of Roth's theorem":
  https://terrytao.wordpress.com/2014/04/24/a-proof-of-roths-theorem/
- Terence Tao, "The transference principle, and linear equations in primes":
  https://terrytao.wordpress.com/2010/06/05/254b-lecture-notes-7-the-transference-principle-and-linear-equations-in-primes/
- Terence Tao, "Give yourself an epsilon of room":
  https://terrytao.wordpress.com/2009/02/28/tricks-wiki-give-yourself-an-epsilon-of-room/comment-page-1/
- Terence Tao, "The correspondence principle and finitary ergodic theory":
  https://terrytao.wordpress.com/2008/08/30/the-correspondence-principle-and-finitary-ergodic-theory/

## What To Imitate

The useful Tao-style move is not verbal flash. It is proof architecture.

1. State the theorem in a version where the proof mechanism is visible.
   In the Roth note, the statement is not introduced only in its most familiar
   integer form; it is moved to a compact abelian group where the Fourier and
   averaging mechanism can be seen cleanly.

2. Name the obstruction before solving it.
   In the transference notes, the primes are not casually treated as dense.
   The zero-density obstruction is stated first, and then the proof imports a
   dense model/pseudorandom measure framework to make the transfer legitimate.

3. Import an outside field as a typed interface, not as metaphor.
   Ergodic theory, topology, Fourier analysis, sieve theory, and functional
   analysis enter through precise roles: compactness, pseudorandomness,
   norms, dense models, local-to-global transfer, or limiting arguments.

4. Give yourself controlled room.
   The epsilon-of-room technique deliberately proves a perturbed or lossy
   statement first, then passes to the limit only after continuity and
   dependency are controlled.

5. Keep the proof's locality visible.
   In the correspondence-principle note, the transfer back from an infinitary
   object works only because the conclusion is local/finitary enough to survive
   the limiting step.

## Translation For EHP114

For EHP114, the proof-facing style should be:

```text
state the local hard-cell theorem
name the obstruction
replace geometry bookkeeping with the right invariant
certify each local closure
then run a separate integration certificate
```

The current invariant is not the operator-grammar sidecar. The current invariant
is concrete:

```text
on |p| = 1, a critical point must satisfy p'(z) = 0
```

That is Tao-style in the relevant sense: compress the obstruction into the
right structure, then verify the transfer back to the original problem.

## Translation For RH/MDL

The RH/MDL lane should copy the transference discipline, not the rhetoric.

Safe pattern:

```text
frozen finite dictionary
explicit bit-cost model
quantization stability theorem
audited finite rows
only then ask whether an asymptotic object exists
```

Unsafe pattern:

```text
finite MDL bend
therefore RH signal
```

The current finite theorem target is exactly the right kind of modest object:

```text
||A q_b(a)-y|| <= ||A a-y|| + sigma_max(A) sqrt(k) Delta_b / 2
```

It states what is controlled and what is not controlled.

## Practical Rule

When importing an outside field into this project, require three things:

```text
typed role        = what the imported field does
loss accounting   = what the import costs
return map        = how the imported statement returns to the original problem
```

For EHP114, the return map is L31-style integration. For RH/MDL, the return map
is a finite audit and then, later, an explicitly stated asymptotic theorem if
one exists.

## Claim Boundary

This note can guide writing and theorem design. It cannot be cited as
certificate evidence. Any proof-facing use must be restated as a theorem,
implemented as an artifact, and verified by checksum.
