[perp] model=sonar-pro tokens=722/1881 cost=$0.03638
[perp] citations (6):
  - https://arxiv.org/abs/2502.09539
  - https://arxiv.org/html/2502.09539v1
  - https://www.renyi.hu/~p_erdos/1975-33.pdf
  - https://www.renyi.hu/~p_erdos/1959-10.pdf
  - https://www.scribd.com/document/707285733/401-Foundations-of-Discrete-Mathematics-With-Algorithms-and-Programming-by-Sriraman-Sridharan-R-Balakrishnan-Z-lib-org
  - https://www.softouch.on.ca/kb/data/CRC%20Standard%20Mathematical%20Tables%20and%20Formulas%2033E.pdf
Q1. Prior art / “quantum entropy = sieve”

I am not aware of any work in analytic number theory that literally recasts a Selberg sieve bound (or sieve for primitive / dilation-constrained sets) as a von Neumann entropy inequality on an explicitly constructed bipartite quantum state, in the sense you propose.

What *does* exist are:

- Information-theoretic or entropy-based proofs of additive/combinatorial results, e.g.  
  – Madiman, Kontoyiannis, et al. on entropy and sumsets;  
  – Tao & Vu use Shannon entropy occasionally in additive combinatorics.  
  But these are classical Shannon entropy, not quantum.

- Connections between multiplicative structures and spectral gaps/graphs that *look* quantum-ish (because you deal with eigenvalues of adjacency or Laplacian matrices), but are not phrased as von Neumann entropy. The “GCD graphs” in Koukoulopoulos–Lamzouri–Lichtman (KLL) are squarely in that camp: spectral graph theory + Selberg sieve, but no quantum formalism.

- A few papers that use “quantum” language as an analogy to discuss sieving or random walks (e.g. quantum walks used to *compute* number-theoretic objects), but not to rederive sieve bounds as quantum-entropy inequalities.

I checked common intersections (“Selberg sieve” with “von Neumann entropy”, “Holevo bound”, “mutual information” etc.) on arXiv and in MathSciNet-style summaries; nothing like “Selberg sieve via strong subadditivity” or “primitive sets via quantum mutual information” appears. The entropy that does arise in sieve literature is typically Kolmogorov/Shannon-type (e.g. entropy of probability distributions on moduli or residue classes), not von Neumann entropy of a constructed bipartite state.

So, as far as current literature goes: your M7 setup is genuinely new as a *formal* recasting; it is not a recognized ingredient in the proof of Erdős #143, nor is there a prior Holevo-type reformulation of the Selberg sieve bounds in KLL.

---

Q2. Does SSA on your |ψ⟩ give a strictly stronger bound than KLL?

Let’s fix your construction:

- Basis states |x⟩ in ℋ_A for x ∈ A ∩ [1,N];
- |k⟩ in ℋ_K for k ∈ {1,…,K};
- amplitudes α(x,k) = 1 if kx ∈ [1,N], else 0.

The amplitude matrix M has entries M_{xk} = α(x,k)/√Z. The well-separated (WS) condition says: for x ≠ y in A and any integer k ≥ 1, |kx − y| ≥ 1. In your context this eliminates certain “collisions” between dilates, which shows up as a sparsity/orthogonality constraint on rows of M: different x’s cannot share too many k’s with the property that kx and ky both land in [1,N] with small spacing.

Two observations:

1. **ρ_A is essentially a (normalized) Gram matrix of the indicator vectors (α(x,·)) in ℓ²({1,…,K}).**  
   The eigenvalues of ρ_A encode how much these dilate-patterns “overlap”. The WS condition ensures overlaps are small in a combinatorial sense. KLL’s GCD-graph + Selberg sieve analysis is precisely a refined way of bounding how many x’s can share strong overlaps in arithmetic constraints; in matrix language, this is controlling the spectrum of a certain incidence/adjacency-like matrix.

2. **von Neumann entropy S(ρ_A) is just H(λ(ρ_A)), the Shannon entropy of the eigenvalues λ.**  
   Applying strong subadditivity (SSA) with a tripartite split A : K_low : K_high amounts to bounding S(ρ_A) via inequalities of the form
   \[
   S(ρ_A) \le S(ρ_{AK_\mathrm{low}}) - S(ρ_{K_\mathrm{low}}) + S(ρ_{AK_\mathrm{high}}) - S(ρ_{K_\mathrm{high}}),
   \]
   or equivalent rearrangements. But since ρ_{AK} is pure, all entropy expressions reduce to entropies of reduced states determined by the same spectrum as ρ_A (up to zeros and multiplicities).

When you push this through, the only quantitative input you get is:

- rank(ρ_A) ≤ min{|A∩[1,N]|, K};  
- eigenvalues are controlled by the ℓ²-norms and pairwise inner products of rows of M.

Those ℓ²-norm and overlap bounds are exactly the combinatorial / “second moment” input that the Selberg sieve is already exploiting in KLL — albeit in a classical way: they count how many x’s survive after sieving by congruences and approximate independence of local conditions, then use large-sieve/ Selberg bounds to limit correlations.

Thus:

- If you rigorously implement M7 + SSA, you end up translating KLL’s combinatorial and sieve estimates into statements like “ρ_A has at most R effective eigenvalues and none larger than p_max”, hence S(ρ_A) ≤ log R + small.  
- Converting this back to a bound on Σ_{x≤N, x∈A} 1/x uses only very rough inequalities (e.g. N_eff ≤ something like Σ 1/x) that do not sharpen the sieve.

In short: SSA and pinching do not source any *new* inequality specific to the WS arithmetic structure; they give generic entropy inequalities valid for all bipartite states with a given Gram matrix. The arithmetic content is entirely in bounding the Gram-matrix entries, which is precisely where KLL already invoke Selberg’s sieve. So M7 is—at least at the level of current techniques—a notational reformulation that repackages the same eigenvalue constraints, not an orthogonal route that yields a strictly stronger Σ 1/(x log x)-type bound.

---

Q3. Mutual information / Araki–Lieb vs KLL’s bound

For a pure bipartite state |ψ⟩_{AK}, the mutual information satisfies
\[
I(A:K) = 2 S(ρ_A).
\]
So bounding I(A:K) boils down to bounding S(ρ_A). You might hope that:

- WS severely restricts I(A:K), and
- some inequality like Araki–Lieb (|S_A − S_K| ≤ S_{AK}) or a refinment thereof would force I(A:K) to be very small, hence constrain the size/weight of A.

But in this construction:

- Araki–Lieb is trivial: S_{AK}=0, so S_A = S_K; it contributes nothing beyond “ρ_A and ρ_K have the same nonzero spectrum”.
- Generic bounds like S(ρ_A) ≤ log rank(ρ_A) ≤ log min{|A|,K} are far weaker than KLL’s analytic bound, which roughly says the harmonic weight of A up to N is ≪ log N / √(log log N). You would need rank(ρ_A) or the entropy of its eigenvalue distribution to be dramatically smaller than log N to beat this, but WS alone does not enforce such a sharp spectral collapse.

Crucially, I(A:K) is invariant under any unitary acting on ℋ_K that mixes the k-labels. The Selberg sieve, by contrast, is exquisitely sensitive to the *arithmetic* placement of the dilation indices (local conditions modulo primes, distribution among congruence classes). The mutual information washes out that structure; it only “sees” how many different k’s each x is linked to and how much overlap there is between x-patterns at a coarse level.

Therefore:

- There is no known Araki–Lieb-style or general quantum-information inequality that, when specialized to this |ψ⟩ built from the dilation condition, yields a bound on I(A:K) stronger than what KLL derive from Selberg sieve.
- Any such inequality would have to be *non-universal*, i.e. depend explicitly on the arithmetical pattern defining α(x,k). Universal quantum inequalities (SSA, subadditivity, Fannes–Audenaert, etc.) are too coarse; they don’t encode congruence or GCD structure.

So currently, the quantum-information viewpoint is best seen as a conceptual or notational reformulation that mirrors KLL’s Gram/spectral estimates. It offers a nice language (entropy, mutual information, rank) for thinking about the combinatorial structure of dilates, but with present tools it does not give a provably tighter bound than Σ_{x∈A,x<N} 1/x = o(log N / √(log log N)).
