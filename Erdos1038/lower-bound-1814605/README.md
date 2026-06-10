# Erdős #1038 — finite-atom lower-bound tail certificate at M = 1.814605

Checked certificates extending the finite-atom dual-forcing tail bound for Erdős
problem #1038 (the infimum of |{U_μ > 0}|) from the published M = 1.8146
(Hua Xu, erdosproblems.com #1038 thread, 2026-05-09) to **M = 1.814605**, a
+5.0e-6 step, together with independent interval re-verifications of the
published 1.8146 package itself.

## What is claimed, and what is not

**Claimed (checked-certificate level):**

- The sweep certificate `certificates/w2_uniform4_M1814605.json` (4 positive
  lanes / 5 atoms, K = 560 blocks, same piecewise structure and conventions as
  the published package) certifies positive on all 560 blocks under directed
  interval arithmetic (IEEE-1788, outward rounding, fail-closed). Per-block
  certified lower bounds: `certify-output/certify_w2_uniform4_M1814605_full.json`
  (worst certified lower bound 7.213773443572941e-09, block 191).
- The forcing certificate
  `certificates/conservative_forcing_interval_certificate_v2_fullrange.json`
  covers the FULL range a ∈ [−1.708, −√2] at cap L = 1.836 (6504 leaves, worst
  certified bound 1.5704408086181232e-06) and passes the independent pure-Python
  interval re-verifier in this bundle (structure / domain / exact-partition /
  inequality, all PASS).
- Re-verification of the published package: Hua Xu's 1.8146 sweep certificate
  certifies positive 560/560 (`certify-output/certify_reference_560blocks.json`),
  and the published 926-leaf forcing certificate certifies positive 926/926 with
  zero sign disagreements (`certify-output/crosscheck_vs_926leaf.json`).

**Not claimed:**

- This is a checked computational certificate, not a written-up theorem and not
  a formal (Lean/ITP) proof.
- The composition consumes thread-shared foundations as-is: the structural
  reductions and normalization (supp μ ⊆ {−1} ∪ [0,1]), the duality-principle
  step, and the two-atom (−√2, 0] window, per the erdosproblems.com #1038
  thread (jspier, J_Koizumi_144, natso26, Hua Xu, Terence Tao). These are not
  re-derived or re-verified here.
- Nothing here addresses the supremum side or any other part of the problem,
  and no claim is made that 1.814605 is optimal for other certificate families.
  Within THIS family the searches indicate a ceiling — M = 1.814606 fails at
  −1.0e-6 and resists polishing — but those searches are local: evidence, not
  exhaustion.
- The forcing certificate does not depend on M: its constants are fixed by the
  cap L = 1.836, so it composes with any sweep value M ≤ 1.836. Minimality
  enters only via the cap (here 1.814605 ≤ 1.82 < 1.836).

## Reproduction (Python only)

Requires Python 3.10+, `numpy`, `scipy`, `mpmath` (the two interval
re-verifiers, steps 2 and 3a, need only the standard library + `mpmath`).

1. **Sweep certificate, independent float-precision check** (mirrors the
   official thread verifier; validated `|diff| = 0` against the published
   1.8146 baseline):

   ```bash
   python3 verifiers/checker.py certificates/w2_uniform4_M1814605.json
   # expect: PASS, worst margin 7.597309819029618e-07, block 196, 0 bad blocks
   ```

   Note — two different worst-block numbers, do not confuse them: this
   float-precision checker reports a worst *float margin* of 7.5973e-07 at
   block 196, whereas the directed-interval kernel (step 3 data) reports a
   worst *certified interval lower bound* of 7.2138e-09 at block 191. They are
   different quantities on different blocks — the interval lower bound is the
   conservative, outward-rounded, fail-closed figure and is the one quoted as
   the certified margin; the float margin is the looser nominal-arithmetic
   slack from the thread-style checker. Both are positive, so both legs pass.

2. **Forcing certificate, full interval re-verification** (pure-Python interval
   backend, independent of the generator):

   ```bash
   python3 verifiers/verify_forcing_interval_v2.py \
       certificates/conservative_forcing_interval_certificate_v2_fullrange.json
   # expect: STRUCTURE/DOMAIN/PARTITION/INEQUALITY all PASS; worst 1.570441e-06
   ```

3. **Certified-bound audit (data):** the two `certify-output/certify_*.json`
   files carry per-block certified lower bounds for the sweep legs (this
   bundle's certificate and the published one). They were produced by an
   internal directed-interval kernel whose source is not included in this
   bundle; the per-block bounds are published as auditable data, and the sweep
   certificate itself is independently re-checkable at float precision via
   step 1. A pure-Python interval re-verifier for the sweep leg is included —
   see step 3a below.

   Honest scope of the asymmetry (as originally staged): the forcing leg was
   fully interval-re-checkable from this bundle alone; the sweep leg's
   interval pass was not. Step 3a closes this.

   Two notes for skeptical readers. (a) LP solver tolerances do not enter the
   validity argument: the LP only finds the candidate; the certificate's frozen
   f64 values are the witness, and the interval pass certifies those frozen
   values directly. (b) The per-block certified values are conservative lower
   bounds — bisection stops refining a box once it clears zero — not minima,
   so they sit well below the float-precision point margins by design
   (worst certified interval bound 7.21e-9 vs worst float point margin 7.60e-7
   on this certificate is expected behavior, not a discrepancy).

3a. **Sweep certificate, full interval re-verification (pure Python).** The
   asymmetry noted above is now closed: the sweep leg is independently
   re-verifiable from this bundle alone, with nothing but Python and `mpmath`
   (standard library otherwise — no numpy/scipy, no compiled code):

   ```bash
   python3 verifiers/verify_sweep_interval.py certificates/w2_uniform4_M1814605.json \
       --compare certify-output/certify_w2_uniform4_M1814605_full.json
   # expect: 560/560 CERTIFIED_POSITIVE, exit 0; worst certified lower bound
   # ~7.2138e-09 at block 191; verdicts match the bundled certify-output
   # receipt block-for-block
   ```

   The verifier uses outward-rounded interval arithmetic (mpmath `iv`) with
   adaptive bisection, taking per box the better of the natural interval
   extension and a mean-value form, with fail-closed verdicts
   (CERTIFIED_POSITIVE / FAIL_TO_CERTIFY / REFUTED) and a `--self-test`
   negative control that corrupts one weight in memory and confirms the
   REFUTED path fires. Certificate f64 values are treated as exact (dyadic
   rationals), the same convention as the other verifiers. Its receipt is
   `certify-output/certify_w2_uniform4_M1814605_full_pyiv.json` — 560/560
   positive, worst certified lower bound 7.213773776423008e-09 at block 191
   vs the bundled receipt's 7.213773443572941e-09 on the same block
   (agreement to ~8 significant digits; identical 25,724-box subdivision
   tree; exact equality of bounds is not expected across engines). The
   internal compiled kernel that produced the original `certify-output/`
   receipts remains withheld, but it is no longer load-bearing for trust:
   both legs of the M = 1.814605 result — sweep and forcing — are now
   independently checkable from this bundle in pure Python. (The
   re-verification receipts for Hua Xu's published 1.8146 package remain
   data-only here; the published certificate itself lives in the #1038
   thread.)

4. **Hashes:** `shasum -c SHA256SUMS` from the bundle root.

## Credit

The certificate family, conventions, K = 560 block sweep, and the 1.8146
frontier package are Hua Xu's; the dual-measure forcing construction and the
lower-bound approach follow jspier and J_Koizumi_144 (with natso26's write-up
of the duality logic and optimizations); the problem normalization / support
reduction (supp μ ⊆ {−1} ∪ [0,1]) follows Terence Tao's thread notes.
This bundle extends the tail by +5.0e-6 and adds interval-level checking; the
heavy lifting upstream is theirs. catsflowers5544's tiling/permutation question
prompted the assignment-layer analysis (answer: tiling kills false certificates
near 1.82, but no single-block reassignment tested opens headroom; in the
computed assignment matrices the binding M-endpoint block already holds its
best slot).

## AI acknowledgment

The research direction, problem identification, alternative framings, and
error-catching are the author's work. AI tools were used inside an integrity
environment designed and maintained by the author (independent-checker gates,
parity gating against validated baselines, fail-closed interval verification,
no-overwrite artifact provenance). Within that environment: Claude (Anthropic,
Fable 5 model family, 2026-05 – 2026-06) for code co-development, LP-search
orchestration, and the interval kernel, alongside Rust/Python toolchains; the
candidate certificate at M = 1.814605 was surfaced by the Claude-orchestrated
search and subsequently verified by the independent checkers included here. All
results were independently verified by the author via the bundled verifiers,
and the certificates are independently verifiable from this bundle alone. No
AI system is an author.

— Kenneth A. Mendoza
