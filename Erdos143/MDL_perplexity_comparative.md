[perp] model=sonar-pro tokens=767/2804 cost=$0.05036
[perp] citations (6):
  - https://real-j.mtak.hu/5464/1/StudScientMath_22.pdf
  - https://www.erdosproblems.com/143
  - https://archytas2.birs.ca/files/publications/proceedings/birs_2006_proceedings.pdf
  - https://www.ams.org/notices/201307/201307-full-issue.pdf
  - https://ndl.ethernet.edu.et/bitstream/123456789/77431/1/80.pdf
  - https://www.jstor.org/stable/pdf/421141.pdf
2. M2 – MULTIPLICATIVE SIDON ANALOGUE (B₂[ℓ]-dilation)

**(A) PRIOR ART**

- Closest objects: multiplicative Sidon / B₂-sets and primitive sets. Key refs:
  - Erdős–Tenenbaum–Zarankiewicz-type multiplicative Sidon sets (e.g. Erdős 1969; see survey: P. Erdős & P. T. Bateman, J. Number Theory 1962; also Ruzsa & others on B₂ sets).
  - Primitive sets: Erdős 1935; Banks–Martin; Granville’s survey.
- For Erdős #143 explicitly, the standard writeups (Tao’s blog, Koukoulopoulos lecture notes, KLL 2025) do not seem to recast it as a “multiplicative Sidon in disguise.” They use GCD graphs / common multiples, not collision multiplicities of the form  
  `kx = y + O(1)` as a Sidon-type combinatorial object.
- I am not aware of a paper that directly states “well-separated sets for Erdős #143 behave like multiplicative Sidon sets with |A∩[1,N]| ≪ N^{1/2+o(1)}.”

**(B) GENUINENESS**

- Conceptually orthogonal: KLL works with GCD-graphs, sieves, and weighted counts; M2 is trying to import *additive* Sidon combinatorics, translated to the dilation relation.
- However, once you formalize “no close kx–y collisions” via inequalities and double-counting, you tend to recover almost exactly the same inequalities KLL use for bounding the number of edges in the GCD-graph. So at a technical level, many Sidon-style arguments end up being “fancy rephrasings” of KLL’s counting.
- Verdict: genuinely different *heuristic* picture, but likely not deeply distinct analytically unless one proves something truly Sidon-sharp (|A∩[1,N]| ≪ N^{1/2+o(1)}). Without that, it collapses back to KLL-style combinatorics.

**(C) NOVEL PREDICTION**

- Strong claim:  
  \[
  |A\cap[1,N]| \ll_\varepsilon N^{1/2+\varepsilon}
  \]
  or even Sidon-quality |A∩[1,N]| ≪ √N.
- This is much stronger than anything in KLL and would in particular force  
  \(\sum_{x\in A} 1/(x\log x)\) to converge.
- It also predicts sharp upper bounds on local density in *short intervals* (Sidon sets spread out evenly), which KLL does not address structurally.

**(D) TRACTABILITY**

- Tools that would be needed:
  - Strong incidence / energy bounds for the dilation relation `(x,y,k) with |kx−y|<1`.
  - A multiplicative analogue of the classical B₂-set Fourier-analytic machinery (e.g., Cilleruelo, Ruzsa, Green–Ruzsa style arguments).
  - Potential use of sum–product phenomena over ℝ or ℚ to show that too many near-collisions force unintended structure.
- Current technology does not give Sidon-level bounds for primitive-like multiplicative constraints; this seems very ambitious.

**(E) RANK**

- Rank: **3/5** – conceptually attractive but likely very hard and technically converges towards KLL’s counting without a genuinely new incidence or sum–product breakthrough.


3. M3 – LOG-POSITION APERIODIC POINT PROCESS

**(A) PRIOR ART**

- Log map appears in several places:
  - Logarithmic density / multiplicative equidistribution in probabilistic number theory.
  - Logarithmic spacing of primes, Poisson processes (Montgomery–Odlyzko).
- For Erdős #143, I do not see log(A)–as–point-process used explicitly in the literature; KLL work in the original scale.
- Related ideas: Furstenberg’s correspondence principle for multiplicative structures; ergodic approaches to multiplicative functions (Tao, Frantzikinakis–Host), but these concern values of arithmetic functions, not this particular dilation avoidance constraint.

**(B) GENUINENESS**

- This is genuinely orthogonal conceptually: it moves from arithmetic combinatorics to dynamical systems / point processes on ℝ with constraints under multiplication by log k.
- Not obviously equivalent to the GCD graph: instead of an intersection graph on integers, you look at orbits of x ↦ x + log k in the log coordinate.
- Still, any rigorous inequality must be pulled back to counting arguments in the original scale; so the analytic content may still recapitulate sieve bounds if not handled carefully.

**(C) NOVEL PREDICTION**

- Predicts that the point process {log x : x∈A} has:
  - Very low “self-correlation” under integer shifts log k.
  - Possibly zero entropy or strong rigidity under the ℤ^+–action, leading to strong upper bounds on local intensity.
- This could yield *weighted* bounds such as  
  \(\sum_{x\in A} f(\log x)/x\) for suitable test functions f, going beyond what KLL’s sieve naturally sees.

**(D) TRACTABILITY**

- Would need tools from:
  - Ergodic theory for ℤ^d-actions with rigidity/aperiodicity constraints.
  - Spectral methods for point processes (Palm measures, correlation functions).
- The big gap is to translate an L²-type spectral estimate into sharp arithmetic bounds on ∑1/(x log x). This is nontrivial and currently not standard.
- Feels conceptually deep and long-term rather than a near-term proof architecture.

**(E) RANK**

- Rank: **2/5** – clearly orthogonal and might encode the “no approximations kx≈y” in a genuinely new way, but technically formidable.


4. M4 – ANTI-BEATTY

**(A) PRIOR ART**

- Beatty sequences, densities, and harmonic sums: classical (Beatty, Rayleigh, Weyl; see e.g. Niven–Zuckerman–Montgomery).
- Erdős-type problems with Beatty sequences arise in additive combinatorics, but I see no explicit “Erdős #143 ↔ Beatty-exclusion” in the literature.
- KLL do not use Beatty or Sturmian descriptions.

**(B) GENUINENESS**

- The idea “A cannot contain a dense approximate arithmetic progression like a Beatty sequence” is basically a structured/anti-structured dichotomy: if A looked like a Beatty, its kx–y approximations would be too frequent.
- However, making this precise may reduce to bounding the number of near-collisions kx≈y, which is again KLL’s business.
- Unless one leverages fine Diophantine properties of slopes (α for ⌊nα+β⌋), the Anti-Beatty framing risks being a rephrasing of density constraints.

**(C) NOVEL PREDICTION**

- Predicts: for each irrational α>1 there is a finite obstruction that prevents A from containing infinitely many terms of ⌊nα+β⌋ at bounded distortion.
- This would say more about *structured subsequences* of A than KLL currently provide; e.g., “A cannot shadow any Beatty sequence too often.”
- However, these predictions are qualitative, not immediately giving a quantitative harmonic sum bound.

**(D) TRACTABILITY**

- Tools:
  - Diophantine approximation and uniform distribution of nα mod 1.
  - Classical estimates for harmonic sums over Beatty sequences.
  - A transference principle: show that any set with “too big” ∑1/(x log x) must correlate with some Beatty; then rule that out.
- The missing ingredient is such a transference principle; there is no off-the-shelf machinery connecting large harmonic mass to Beatty correlation in this multiplicative context.

**(E) RANK**

- Rank: **5/5** – interesting heuristic, but least concrete as a proof architecture and closest to being a rephrasing of “A has no dense pseudo-progressions.”


5. M5 – INFORMATION-THEORETIC ENTROPY DECAY

**(A) PRIOR ART**

- Entropy and forbidden patterns are standard in symbolic dynamics and combinatorics on words (subshifts of finite type, sofic shifts).
- In number theory, there are entropy methods for multiplicative functions (Tao–Teräväinen, Frantzikinakis–Host).
- For Erdős #143 specifically, I see no explicit entropy-rate treatment of 1_A(n) with pattern constraints coming from |kx−y|≥1.

**(B) GENUINENESS**

- Mapping A to a binary sequence with long-range forbidden patterns (depending on multiplicative relations) is quite different from KLL’s graph-theoretic sieve.
- If you can show that any process satisfying these constraints has entropy rate zero, you are asserting an extreme sparsity that goes beyond current sieve bounds.
- Genuinely orthogonal in philosophy; not obviously equivalent.

**(C) NOVEL PREDICTION**

- Predicts: the upper Banach density of A is 0, and more strongly, that any shift-invariant measure supported on such sequences has entropy 0.
- In particular, any “typical” realization would have counting function n^{o(1)}, which forces ∑1/(x log x) < ∞.
- This ties the problem to dynamical classification of low-entropy subshifts defined by multiplicative exclusions.

**(D) TRACTABILITY**

- Would need:
  - A symbolic encoding where the WS condition appears as finitely many forbidden patterns (or a well-controlled almost-local constraint).
  - Subadditive entropy bounds for this subshift, akin to classical arguments for SFTs but with nonlocal interactions.
  - Possibly an application of Ornstein–Weiss or Gurevich entropy bounds for countable-state shifts.
- The crux: WS is “global and multiplicative,” not a finite-radius constraint. Turning this into a useable entropy bound is highly nontrivial; no ready-made theory handles exactly this.

**(E) RANK**

- Rank: **4/5** – truly orthogonal, but the gap from a global multiplicative constraint to quantitative entropy decay is wide.


6. M6 – KOOPMAN / TRANSFER-OPERATOR SPECTRAL

**(A) PRIOR ART**

- Koopman operators and transfer operators are standard in ergodic theory and dynamical zeta functions; in number theory, see work on the Gauss map, continued fractions (Baladi, Mayer, etc.).
- For multiplicative problems (e.g., squarefree numbers), transfer operators appear in thermodynamic formalism.
- For Erdős #143 or primitive sets, there is no standard Koopman-operator framing in the literature; KLL do not use spectral theory of operators.

**(B) GENUINENESS**

- This is close in spirit to M3/M5 but more analytic: define an operator encoding dilations and examine its spectrum.
- Not obviously equivalent to KLL; KLL’s bounds come from combinatorial and sieve estimates rather than spectral gaps of a linear operator.
- That said, if one only recovers L²-energy bounds of the type KLL already prove, then this is merely a spectral rephrasing. The value would come from a *new* spectral gap not visible by combinatorial means.

**(C) NOVEL PREDICTION**

- Predicts that the spectral radius of the “dilation operator” on a suitable Hilbert/Banach space is strictly less than 1 in the relevant direction, implying exponential decay of certain correlations.
- Translated back, this would give exponential or at least power savings in counts of kx≈y, much stronger than KLL’s o(log n/√log log n).
- Could imply very strong weighted sum bounds, e.g. ∑_{x∈A} x^{-s} convergent for s near 1, which KLL do not touch.

**(D) TRACTABILITY**

- Needs:
  - A precise dynamical system whose invariant sets correspond to well-separated A.
  - A function space where the transfer operator for dilations is quasi-compact with a gap.
  - Thermodynamic formalism or Dolgopyat-type estimates to prove the gap.
- Building such a system is highly nontrivial; there is no canonical Markov partition or symbolic coding as in continued fractions. This is heavy-duty analysis.

**(E) RANK**

- Rank: **1/5** – as an “orthogonal MDL-morphism,” this is the most structurally different from KLL and, if made precise, could give qualitatively new predictions. Technically, it is also the most ambitious, but it best fits your goal of an orthogonal architecture that might force part (ii).
