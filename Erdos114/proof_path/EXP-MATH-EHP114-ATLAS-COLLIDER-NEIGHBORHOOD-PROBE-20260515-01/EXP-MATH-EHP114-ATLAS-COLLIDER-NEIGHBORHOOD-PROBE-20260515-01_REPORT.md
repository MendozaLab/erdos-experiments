# EXP-MATH-EHP114-ATLAS-COLLIDER-NEIGHBORHOOD-PROBE-20260515-01

## Final Report (Atlas + Collider neighborhood probe of Erdős Problem #114)

### Neighbor counts confirmed

Queried the canonical workbench Atlas DB (`Math/erdosatlas-workbench/erdosatlas.db`, 17 MB, 5,496 problems, 54,203 morphisms). #114 has **0 curated** edges, **3 scout** edges, **13 dense** edges, and 123 cross-corpus 51xxx edges (excluded per plan as low-signal bulk-cosine). The 16 curated neighbors break down as **8 solved, 4 disproven, 4 open** (plan recon said 7/4/5 — the off-by-one is from #114 itself being mis-counted as a neighbor). Recon list verified: technique-transfer (#139, #720, #1120), structural (#72 at edge 0.79, #115), generalization (#477), and ten constraint-parallel (#116, #1114, #511, #1048, #1046, #671, #1047, #1150, #228, #525). None of the seven Collider `l4_*` baby-steps targets #114 directly; they are MDL-floor infrastructure for the existing finite-n=15 EHP bridge program. The 75 prior `EXP-MATH-EHP114-*` directories all attack #114 directly via Fryntov-Nazarov 2008 + Tao 2512.12455 + finite-n bridge, not via neighbor transfer.

### Top transfer candidate

**#116 (solved, planar measure of |{|p|<1}|)** — rank 1, plausibility **LOW**. Technique is subharmonic potential theory on u = log|p| with Cartan-type exceptional-set covering. Hypothesis: a coarea-formula bridge (length-of-level-set integrated against |∇u|) could convert the area estimate into a length estimate, but per Perplexity call 1 the required pointwise gradient control near {u=0} is exactly what the Fryntov-Nazarov 2008 Stokes/Cauchy/area-integral toolkit already supplies in the live 114 program. Transfer lands back in the existing toolbox — no shortcut. Ranks 2 and 3 (#1120 dual-extremization SPECULATIVE; #115 Bernstein-Walsh LOW) are even weaker.

### Perplexity calls used

**4 of 12 budget** (sonar-low, $0.005 each = $0.02 total). Calls covered: #116 transfer (call 1), #1046/#1047/#1048 disproof counterexamples (call 2), #1120 + #115 (call 3), Littlewood/BBMS #228/#525/#1150 (call 4). All four confirmed: techniques in the neighborhood do not transfer to arc-length on unconstrained monic polynomials.

### Intractability verdict

**INTRACTABILITY_CONFIRMED_BY_NEIGHBORHOOD_PROBE.** The Atlas-curated neighborhood of #114 offers no non-brute-force shortcut. Disproof counterexamples (#511, #1046–#1048) bound topology not length. Littlewood-family results (#228, #525, #1150) live on a fixed circle under coefficient constraints, orthogonal to arc-length on arbitrary monic. The strongest open structural analog #1120 is a dual problem with no known bridge to #114. Solved polynomial problems #115/#116/#1114 share the analytic primitives but require an extra coarea/Stokes step that is precisely the existing Fryntov-Nazarov + Tao machinery. Off-domain edges to #72/#139/#720 are spurious Atlas tags. The live attack remains the finite-n≤14 bridge + Fryntov-Nazarov + Tao asymptotic-improvement chain.

This verdict is a probe diagnostic only: the *neighborhood* is silent, not #114 itself.
