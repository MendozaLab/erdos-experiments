# In Search of Erdős #114

Tao's large-degree theorem, the small-degree certificate ladder, and the
carpenter's hammer

Assembled: 2026-05-07  
Artifact horizon: local EHP114 artifacts through 2026-05-06  
Status: working narrative for review and orientation, not a proof package

## Claim Boundary

This note is a story of the search path. It is not evidence by itself.

It is not a proof of Erdős #114. It is not a global n=14 certificate. It is not
an exact lemniscate-length certificate. It does not use the RH/MDL/operator
grammar sidecars as mathematical evidence.

The only proof-facing objects in this lane are theorem-shaped local artifacts,
source hashes, interval arithmetic reports, residual-chain integration checks,
and the eventual global atlas certificate.

## The Problem

Erdős, Herzog, and Piranian asked whether the longest lemniscate of a monic
degree-n polynomial is the symmetric one:

```text
p(z) = z^n - 1.
```

The lemniscate is the curve

```text
|p(z)| = 1.
```

The conjectural extremizer is simple to say: put the roots evenly around the
unit circle. The difficulty is that a polynomial lemniscate is not just a
polygon or a smooth perturbation of a circle. Its length is controlled by root
geometry, critical points, singular boundary layers, and how local branches
can be counted without double-counting a thick numerical tube.

That is why this problem is a good test of proof architecture. The statement
is one line. The certificate is not.

## The First Public Pass: A Glimpse Of The Landscape

The first public artifact was the Zenodo record
`10.5281/zenodo.19480329`, titled "Computational Verification of the
Erdős-Herzog-Piranian Conjecture for Degrees 3 <= n <= 14." It preserved the
initial finite computational surface: interval-arithmetic code, result JSON
files, checksums, and the closed-form benchmark for the symmetric curve.

That record mattered because it made the work auditable. It gave us a map of
the landscape: where the numerical frontier was, what the expected extremal
lengths were, and which degrees could be tested by a reproducible engine rather
than by prose.

But a landscape map is not the same thing as an orthodox proof packet. The work
after Zenodo has been a tightening pass: take the promising finite surface and
ask what a skeptical analyst would require before accepting the n=14 frontier
as a theorem-shaped certificate.

The answer was not "more pictures." The answer was a chain of local certificates.

## Tao's Large-Degree Endpoint

Tao's high-degree work changes the strategic shape of the problem. It supplies
the large-n horizon: for sufficiently large degree, the symmetric extremizer is
the right object. That means the full completion path is no longer an
undirected search across all degrees.

The architecture becomes:

```text
large n: external analytic theorem to audit and effectivize
small n: finite certificate ladder
middle range: whatever remains after the large-n threshold is made explicit
```

This is a major reduction in how to think. It does not finish the finite work.
The phrase "sufficiently large" still has to be made operational if it is going
to close the entire problem. But it tells us where our work belongs: not as a
replacement for Tao's asymptotic theorem, but as the small-degree certificate
side of the bridge.

## The Carpenter's Hammer

The early mental model was Hessian-shaped: if the symmetric polynomial is the
extremizer, maybe local positivity of a transported shape cone will certify
nearby perturbations.

That hammer broke in a useful way.

The interval Taylor packet showed a split verdict. Axis endpoint checks looked
friendly, but the uniform Taylor-matrix positivity check failed. Spectral
diagnostics then showed the failure was not merely a loose Gershgorin artifact:
negative directions appeared under the stencil. The naive theorem "radial
contraction preserves uniform shape positivity" had to be retired.

That was not wasted work. It taught the right local sentence:

```text
radial reserve must dominate radial-base shape softening.
```

The radial direction is singular. It has a Puiseux boundary layer, not an
ordinary quadratic Hessian. The later radial/Puiseux packet certified the
radial reserve for the family

```text
p_a(z) = z^14 - a
```

with the fixed lower bound

```text
L_14(1) - L_14(1 - eps) >= 24 eps^(1/14)
```

on the stated local range. That clarified the analytic geometry but did not
close the nonradial problem. The hammer was not discarded; it was reshaped.

The lesson was simple: do not force the curve into the coordinate system you
wish it had. Follow the invariant it actually gives you.

## The Bridge That Failed For The Right Reason

The next natural attempt was to use marching squares and normal drift. That
was useful for locating geometry, but it was the wrong proof anchor. A
marching-squares curve is a diagnostic drawing; the exact lemniscate is an
implicit curve.

The direct validated-length run proved regularity and ownership, but its
length bound was far above the cap because it counted a thick interval tube.
It was measuring the whole uncertainty strip, not one owned branch.

That failure was precise enough to be progress. It forced the branch-atlas
frame:

```text
one owned curve branch
one chart or exclusion certificate
one length contribution
no duplicate ownership
```

The slab validator then showed the length budget was not the main enemy. On
the hard n=14 subcell `(6,4)`, the accepted branch length was

```text
20.316451752723314 < 20.672796062619668.
```

The enemy was isolation: unresolved branch tubes still had to be excluded or
certified.

## Collars, Krawczyk, Bernstein, And What They Taught

Several reasonable local routes were then tested and retired as proof-facing
methods:

- axis-aligned collars could not reliably get wall sign separation;
- gradient-normal collars certified no branches in the pilot;
- third-order local collar bounds exposed the center-strip term as the blocker;
- one-dimensional Krawczyk under full root-affine uncertainty certified no
  enough of the remaining tubes;
- parameter slicing certified zero pieces in the pilot;
- bivariate Bernstein/Krawczyk produced zero exclusions and zero certifications
  on the sampled hard pieces, with hulls far too wide.

The important discipline is what these failures do and do not mean.

They do not mean the curve is too long. They do not mean the symmetric
candidate is wrong. They mean generic subdivision bookkeeping is the wrong
primitive for the remaining obstruction.

The invariant that survived was this:

```text
on |p| = 1, a critical point of the level curve must satisfy p'(z) = 0.
```

That compressed the local obstruction from "all possible branch geometry" to
"exclude derivative roots from the remaining candidate boxes."

## The Current Proof-Facing Pipeline

The local n=14 proof lane now runs as a source-to-certificate chain for each
root-affine cell:

```text
slab source
source smoke / subcell contract
branch and critical-region diagnostics
regular residual closure
p-prime candidate exclusion
p-prime root-location exclusion
sharp p-prime root-location exclusion
residual-chain integration
local cell certificate packet
```

The chain is intentionally redundant. It checks that:

- the source artifact belongs to the requested subcell;
- no hard-cell data is silently reused;
- ownership keys do not duplicate curve branches;
- candidate counts do not drift between stages;
- source hashes pass;
- every local residual bucket is either closed or named as a blocker;
- the local length upper bound stays below `20.672796062619668`.

For the original hard cell `(6,4)`, this became the L32 local hard-cell packet.
For the global n=14 atlas, it then had to be parameterized cell by cell. That
software-proof boundary became its own result: the pipeline was repaired so
that non-hard-cell runs must carry explicit subcell identity and fail on source
mismatch.

## Where The n=14 Gauntlet Stands

As of the current local artifact scan, the n=14 grid has:

```text
48 / 64 local cell certificate packets passing
75.0% of the local cell gauntlet
16 cells still open
```

The tightest passing cells are:

```text
CELL-02-04  total <= 20.469469503184456  margin 0.2033265594352116
CELL-04-03  total <= 20.222651196010045  margin 0.45014486660962305
CELL-01-00  total <= 19.945524819071935  margin 0.7272712435477331
CELL-03-05  total <= 19.877489662022747  margin 0.795306400596921
CELL-05-02  total <= 19.78988981798287   margin 0.8829062446367963
```

The controlled remaining-cell batch stopped honestly at:

```text
CELL-05-06
```

The length bound there is still comfortably under the cap:

```text
17.70719360051437 < 20.672796062619668.
```

The blocker is not length. The blocker is residual-chain closure:

```text
60 / 64 residual buckets closed
4 residual buckets still require a local closure route
```

That is the next proof-track target. We should not restart global sweeps, and
we should not patch with prose. The correct next move is to inspect the four
unclosed `CELL-05-06` residual buckets, identify whether they are regular
slice, p-prime, or source-filter leftovers, and add the narrowest certificate
that closes them without changing the accounting contract.

## What Still Has To Happen For A Full Completion

The full proof ladder still has four gates.

First, finish the n=14 local atlas: all 64 local cell packets must pass under
the same source/subcell/ownership/hash discipline.

Second, emit the global n=14 atlas certificate. A pile of passing local packets
is not enough; there must be a separate integration artifact proving that the
cells cover the intended parameter atlas with no gaps, no duplicate ownership,
no source-filter mismatch, and no candidate-count drift.

Third, build the finite-degree certificate index. Degrees 3 through 14, and
any lower-degree literature anchors such as n=2, must be separated by source:
literature theorem, local computation, checksum-backed artifact, or unresolved
gap.

Fourth, audit Tao's high-degree threshold. If the large-n theorem can be made
explicit at a practical cutoff, the finite computation range is finite and
known. If not, the main remaining work is an analytic cutoff-improvement lane.

Only after those gates pass does a full synthesis packet make sense.

## The Meaning Of The Search

The story is not that a clever numerical picture became a theorem.

The story is that almost every tempting shortcut was forced to leave an
artifact-backed reason for its retirement. Hessian positivity became radial
reserve versus shape softening. Marching squares became branch ownership.
Branch collars became derivative-root exclusion. Hard-cell success became a
parameterized local pipeline. The local pipeline became a 64-cell atlas
problem.

That is the carpenter's hammer in its useful form: not one magic tool, but the
discipline of replacing a broken tool with a sharper invariant each time the
artifact says where it failed.

For an orthodox reader, that is the strongest current selling point. The work
does not ask for nimbleness instead of rigor. It asks for rigor about the
discarded routes too.

## Primary Local Anchors

- Public first-pass record: https://zenodo.org/records/19480329
- Zenodo DOI: https://doi.org/10.5281/zenodo.19480329
- Tao high-degree arXiv paper: https://arxiv.org/abs/2512.12455
- Tao high-degree blog note:
  https://terrytao.wordpress.com/2025/12/15/the-maximal-length-of-the-erdos-herzog-piranian-lemniscate-in-high-degree/comment-page-1/
- Radial/Puiseux synthesis:
  `erdos-experiments/Erdos114/EHP114_N14_RADIAL_PUISEUX_CLOSURE_SYNTHESIS_2026-05-05.md`
- Method ledger:
  `erdos-experiments/Erdos114/EHP114_METHOD_LEDGER_2026-05-06.md`
- Local hard-cell packet:
  `erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-HARD-CELL-CERTIFICATE-PACKET-20260506-01/`
- Controlled remaining-cell batch:
  `erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CONTROLLED-REMAINING-CELL-BATCH-20260506-03/`

## Review Use

Use this note to orient a reader. Do not cite it as a certificate. Every
mathematical statement that matters has to point back to a source artifact,
checksum, theorem statement, or literature theorem.
