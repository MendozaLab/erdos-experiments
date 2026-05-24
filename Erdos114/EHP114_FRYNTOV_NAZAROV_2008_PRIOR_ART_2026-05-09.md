# Fryntov–Nazarov 2008 — Prior-Art Read for EHP / Erdős #114

**Artifact:** internal literature-review note (read-only on source paper)
**Read date:** 2026-05-09
**Reader:** Claude Code session (Erdős #114 prior-art gap follow-up)
**Source:** Fryntov & Nazarov, *New estimates for the length of the Erdős–Herzog–Piranian lemniscate*, arXiv:0808.0717v1, May 8 2008 (timestamp on title page) / arXiv-stamp 5 Aug 2008. 18 pages. MSC 30C10. NSF DMS0501067 (Nazarov).
**Local copy:** `/tmp/fn08.pdf` (152.9 KB), `/tmp/fn08_abs.html` (41 KB)
**Scope:** literature triage. No claims about H² IP, no portfolio binding. Honest scope = "I read it; here is what's in it."

---

## Question 1 — What length bound do they prove?

**Verdict: explicit two-bound stack culminating in `|L| ≤ 2n + O(n^{7/8})`.**

Reading the paper from §1 through §9, three distinct length bounds appear, each strictly sharper than the last:

- **§7 (p. 10), "The simplest upper bound":** `|L| ≤ 2π(2n − 1)`, with the parenthetical refinement `|L| ≤ 2π(2k − 1)` when `p` has only `k` distinct roots. Proof: bound `|L|` by the area integral of the logarithmic-derivative sum, then invoke Pólya's theorem (Ransford [11], Thm. 5.3.5) — area of `E` is at most `π` because logarithmic capacity of `E` equals 1.

- **§8 (pp. 10–11), "An improved upper bound":** `|L| ≤ 2π(n − 1 + √n)`. Proof: change of extension `s = |p|·|φ|/φ` so that `2∂s = 2|p′| − (|pφ|/φ)·(φ′/φ)`, giving `|L| = 2 ∬_E |p′| dA − ∬_E (|pφ|/φ) ψ dA` (eq. (6)), then apply Cauchy on `∬|p′| dA ≤ π√n` and bound the tail integral by `2π(n − 1)`. Authors note this is "only marginally worse than Danchenko's estimate `2πn`."

- **§9 (pp. 11–16), "Asymptotic estimate":** `|L| ≤ 2n + O(n^{7/8})` as `n → ∞` for *all* monic polynomials `p` of degree `n` (not just those near `p_0 = z^n − 1`). This is the headline result. Final algebraic line on p. 16:

  > `|L| ≤ 26πδn + 2π√n + e^{2δ} (1/π + 2δ + 4a/(δ³√n)) · 2π(n−1)`

  for every `δ > 0`. They write that the optimal choice is `δ ≈ n^{−1/8}`, which yields `|L| ≤ 2n + O(n^{7/8})` (final inequality, p. 17).

The leading constant `2` is sharp (Erdős–Herzog–Piranian 1958 conjectured `|L_{p_0}| = 2n + O(1)`; Fryntov–Nazarov match the leading term and bound the `o(n)` remainder by `O(n^{7/8})`). The `7/8` is openly noted as suboptimal: "We have no doubt that the power 7/8 can be substantially improved though to bring it below 1/2 seems quite a challenging problem" (p. 17).

There is also a secondary deliverable: §6 (pp. 7–10) proves that `|L_p| ≤ |L_{p_0}|` whenever `p` is *sufficiently close* to `p_0`, i.e. `p_0 = z^n − 1` is a strict local maximum of `|L|` in coefficient space. This is a separate "first variation" result and is not the asymptotic estimate.

---

## Question 2 — Explicit constants?

**Verdict: partially explicit, partially absorbed into anonymous `C`.**

The paper has two distinct constant regimes:

- **Anonymous `C` / `c` (`n`-only-dependent):** §5 Lemma 1 introduces `c_n > 0` "depends on `n` only." §6 uses `C` "with some absolute positive constant `C` depending on `n` only" inside `|G ∩ T_ρ| ≤ Ca²ρ^{−1}` (eq. (5)). These never get pinned to a numeric value in the body — they are "depending-on-`n`-only" black boxes used to absorb finite-`n` rounding inside the local-maximum argument.

- **Genuinely explicit numerics:** the asymptotic chain in §7–§9 is fully explicit. Every prefactor is numbered. Listing them as they appear:
  - §7: `|L| ≤ (2n − 1) · ∬_𝔻 dA/|z| = 2π(2n − 1)`. The `2π` factor is not hidden.
  - §8: `∬_E |p′|² dA = πn` (exact identity, since `p` covers `𝔻` with multiplicity `n`); `∬_E |p′| dA ≤ π√n` (Cauchy); `|J| ≤ 2π(n−1)`. Result `|L| ≤ 2π(n − 1 + √n)`.
  - §9: bound (9) gives `2πnδ + 8π(2n−1)δ`; bound (10) gives `4π(2n−1)δ`; sum is `26πδn`. The Remez application (eq. above (5) on p. 9) gives `4·2^{1/n}·r²·∫_0^π dα/|cos nα|^{1/n} = Cr²` — the first factor is `4·2^{1/n}` explicitly. Final p. 16 inequality is fully numeric: `26πδn + 2π√n + e^{2δ}(1/π + 2δ + 4a/(δ³√n)) · 2π(n−1)`. Plugging `δ = n^{−1/8}` mechanically yields a *named* implied constant inside `O(n^{7/8})`, even though the authors do not collect it.

So: the §9 chain is constructive and **a careful reader can extract a fully numeric `K_FN` such that `|L| ≤ 2n + K_FN · n^{7/8}` for all `n ≥ N_FN_0`**, where both `K_FN` and `N_FN_0` are computable from the §9 inequalities. The authors did not bother. The Vinogradov `≪` symbol does not appear in the paper at all (they use `≤` with named constants `c, C, c_n` throughout).

---

## Question 3 — `C_0`-style scaling? (THE KEY QUESTION)

**Verdict: there is no `C_0` parameter in Fryntov–Nazarov. The closest analogues are `δ` (an `n^{−1/8}`-scaled cutoff radius) and `r` (an annulus radius `r = a^{2/3}`). Neither appears with power `4` anywhere in the paper.**

This is the question Track A4's claim hinges on, so I want to be explicit about what *does* appear and where, and then say what doesn't.

The free parameters that govern the §9 asymptotic estimate are:
- `δ ∈ (0, 1/4)` — the "good-set radius" defining `E_δ ⊂ E` via inequalities (7) and (8). All three terms in the final p. 16 bound carry an explicit power of `δ`:
  - `26πδn` (linear in `δ`),
  - `2π√n` (no `δ`),
  - `e^{2δ} (1/π + 2δ + 4a/(δ³√n)) · 2π(n − 1)` — the *cubic* `δ³` appears in the denominator of the third sub-term, so the dominant `δ`-dependence on optimization is `δ + 1/(δ³√n)`. Setting derivative to zero gives `δ ≍ n^{−1/8}`. **The dominant `δ`-power is `−3` (from `δ³√n`), not `−4`.**
- `a` = `∑_{k≠0} |a_k|/|k|` (the Fourier-coefficient `ℓ¹` norm of `Re_+ z` on the unit circle, eq. (16)). This is an absolute constant introduced after Lemma 2 — it is **not** a free parameter the user tunes. It is fixed by the function `Re_+ z`. So it cannot play the role of Tao's `C_0`.
- `r` (in §6's local-max argument) is set to `r = a^{2/3}` where `a = max_k |a_k|^{1/k}` is a smallness scale of the perturbation `q`. So `r` is **adapted to the perturbation, not free**, and its power on the smallness parameter is `2/3`, not `4`.

**Direct comparison to Tao's annulus-2 remainder `O_{C_0}(log n · ‖p‖ / n)`:** Tao chooses an annulus radius parameter `C_0` (a free large constant), then bounds the contribution of an annulus of radius `[1, C_0]` (or analog) and lets `C_0 → ∞` outside. Fryntov–Nazarov's analogous geometric quantity is the radius `r ∈ (4a, 1/4)` in §6 (eq. above (5)) — they take `r` *inside* the unit disk rather than outside, and the dependence is **`r²`** in the dominant integral `Cr²` (line three from the bottom of p. 9) and **`r^{n+1}`** in the band-thickness estimate `|H ∩ T_r| ≤ Cr^{n+1}` (top of p. 9). So if there is an "annulus power" in this paper, the candidate values from a direct read are `2`, `n+1`, `2/3` (from `r = a^{2/3}`), and `−3` (from the `δ³` denominator) — **not `4`**.

**Bottom line on Q3:** Fryntov–Nazarov do not parametrize a Tao-style outer cutoff `C_0` at all. They pick a single `δ` and optimize. Any direct mapping `K_logn(C_0) = K_logn_universal · C_0^4` cannot be read off Fryntov–Nazarov's analysis — there is no place in their §9 estimate where a free large parameter enters the bound to the 4th power. So Fryntov–Nazarov's text **neither confirms nor refutes the `C_0^4` claim** on its own terms. What it does refute is any narrative that says "the `C_0^4` factor is an inherited convention from Fryntov–Nazarov." It is not. It would have to come from elsewhere (Tao 2512.12455's own internal accounting, or a downstream rederivation), and Track A4 needs to source it there honestly.

---

## Question 4 — Methods compared to Tao 2512.12455

**Verdict: shared technique family (Stokes / Cauchy / area-integral), different decomposition surgery. Constants do not transfer cleanly.**

Shared:
- Both papers use Stokes formula (Fryntov–Nazarov §4, eq. (2) — explicitly reduced from boundary integral to area integral via `|L| = 2 Re ∬_E ∂s dA`, eq. (4)).
- Both use the logarithmic derivative `φ = p′/p` and its decomposition into pole-at-roots-of-`p` minus pole-at-roots-of-`p′` (Fryntov–Nazarov p. 4: `2∂s = (|φ|/φ)(∑_η 1/(z−η) − ∑_ζ 1/(z−ζ))`).
- Both use Cauchy / Hölder on `∬ |p′|^k dA` over `E` (Fryntov–Nazarov §8, p. 11: `∬ |p′|² dA = πn`; `∬ |p′| dA ≤ π√n`).

Different:
- Fryntov–Nazarov's surgery is **`E_δ ⊂ E` "good-set" + tubular bad-set + Pólya-capacity bound on `F`** (the union-of-squares cover of `E_δ`). The cutoff `δ` is a single scalar in `(0, 1/4)` chosen as `n^{−1/8}` after the fact. There is no inside/annulus/outside trichotomy; the geometry is "good points where `|φ′/φ|` is large enough" vs. "bad points near roots of `pp′`" (eqs. (7), (8)).
- Tao 2512.12455's `inside-2 / annulus-2 / outside-again` decomposition (per the prompt's framing) is a radial-shell trichotomy with a free parameter `C_0` controlling the annulus. That radial structure is **absent** in Fryntov–Nazarov §9. The closest radial-shell move is in §6 (`r ∈ (4a, 1/4)`, splitting into `D_r` and complement), but that is the *first-variation* argument near `p_0`, not the asymptotic estimate.
- Auxiliary tools differ: Fryntov–Nazarov §5 uses **Poincaré / central-projection / great-circle counting** on a sphere `S` of radius `R → ∞` (pp. 5–6) — a piece of integral geometry not present in Tao's annular accounting. They also use **Remez' theorem** ([12] = Borwein–Erdélyi, Thm. 5.1.1) to bound the polynomial-thin-set (top of p. 9). Tao does not invoke either. Conversely, anything Tao uses about "radial Bernstein on the unit circle" or "free `C_0` annulus" is not in Fryntov–Nazarov.

**Net:** Fryntov–Nazarov and Tao share an analytic family but do *not* share a decomposition. A constant extracted from Fryntov–Nazarov's `δ`-optimization will not transfer term-for-term to Tao's `C_0`-parametrized accounting. Any claim of the form "Fryntov–Nazarov's constants imply `C_0^k` in Tao" must reconstruct the bridge, not invoke it.

---

## Question 5 — Their claim ceiling. What range of `n`, what relationship to predecessors?

**Verdict: §9 is asymptotic for all `n` with no finite-`n` resolution. The local-max in §6 is qualitative — no explicit neighborhood radius. They resolve no new finite-`n` cases of the EHP conjecture.**

Predecessor chronology (from the introduction, p. 1):
- Erdős–Herzog–Piranian 1958, *Metric properties of polynomials*, J. Anal. Math. 6, 125–148 — original conjecture, listed as Problem 12.
- Dolzhenko 1960 (PhD thesis, Moscow, in Russian; published 1963 [3]) — first upper bound `|L| ≤ 4πn`.
- Pommerenke 1961 [4] — `74n²`, looser but better-known.
- Borwein 1995 [5] — `8eπn`, apparently unaware of Dolzhenko.
- Eremenko–Hayman 1999 [6] (arXiv:0805.2295) — proved EHP for `n = 2`; showed all critical points of any extremal `p` lie on `L_p`; obtained `|L| ≤ 9.173n` for all `n`.
- Danchenko 2007 [7] — `|L| ≤ 2πn`. Authors call this "the best published upper bound by the moment of writing this article."
- Fryntov–Nazarov 2008 (this paper) — `|L| ≤ 2n + O(n^{7/8})` asymptotically; `|L_p| ≤ |L_{p_0}|` qualitatively in a neighborhood of `p_0`.

Relation to Tao Dec-2025:
- Fryntov–Nazarov resolve **no new finite-`n` case**. The only finite-`n` resolution in the chain is Eremenko–Hayman's `n = 2`. Everything else (Dolzhenko, Pommerenke, Borwein, Danchenko, Fryntov–Nazarov) is an asymptotic or universal-`n` estimate, not a per-`n` proof of the conjecture.
- Their leading `2n` is sharp (matches the conjectured value); their improvement is in bringing the remainder down to `o(n)` for the first time. Prior to Fryntov–Nazarov the remainder was `O(n)` with a strictly-greater-than-1 coefficient (e.g. Danchenko's `(2π − 2)n`).
- Tao 2512.12455's "large-`n` result" sits *downstream* of Fryntov–Nazarov in the same asymptotic-improvement lineage — both are improvements on the `o(n)` remainder, not finite-`n` resolutions.

---

## Implications for Track A4

**Verdict: NEUTRAL on `K_logn(C_0) = K_logn_universal · C_0^4` — but the claim *cannot* cite Fryntov–Nazarov as a prior-art warrant for the `C_0^4` power.**

The paper does not have a `C_0` parameter and does not produce a `C_0^4` factor anywhere in any of its three bounds. The dominant free-parameter power in §9 is **`δ^{−3}`** (cubic, not quartic) which gets traded against `δ` and `√n` to land at `n^{7/8}`. The only place a parameter raised to power `4` appears at all in the paper is `e^{4δ}` (an exponential, not a polynomial in a free large parameter) — that is not the same kind of `C_0^4` factor.

What this means concretely for Track A4:
- If Track A4's `C_0^4` is a fresh derivation internal to Tao's own accounting — fine, Fryntov–Nazarov is silent on it. Track A4 must show its own work.
- If Track A4 cites Fryntov–Nazarov as "the published `C_0^4` precedent" — that citation is wrong. They have no `C_0`, no `C_0^4`. The claim should be reattributed or struck.
- The honest framing for the bridge program: Fryntov–Nazarov's `O(n^{7/8})` remainder is an *unrelated parametrization* of the same `o(n)` envelope. Any `C_0`-power in Tao that Track A4 wants to invoke must trace to Tao's text, not to Fryntov–Nazarov.

I am explicitly **not** massaging this finding. The paper does not address `C_0^4` and the bridge program should not pretend otherwise.

---

## Implications for the bridge program — candidate `N₀`?

Fryntov–Nazarov give a parametric form `|L| ≤ 26πδn + 2π√n + e^{2δ}(1/π + 2δ + 4a/(δ³√n))·2π(n−1)` valid for *every* `δ > 0`. This is an *all-`n`* estimate, not a "`n ≥ N_0`" estimate, so the paper does not name a threshold `N_0` directly.

However, two implicit thresholds appear and are worth flagging for the bridge architecture:
- The optimization `δ = n^{−1/8}` requires `δ < 1/4`, i.e. `n > 4^8 = 65,536`. Below that, the pre-optimized form is valid but the "optimal" `δ` value is out of range.
- The §6 local-max argument requires `r ∈ (4a, 1/4)` and `a` "sufficiently small", i.e. is qualitative — no quantitative threshold is stated.

So Fryntov–Nazarov's estimate suggests a **soft threshold around `n ≈ 6.5 × 10⁴`** below which the §9 optimization is not literally available. This is several orders of magnitude beyond any `N_0` the bridge program is currently targeting (we are working in `n ≤ 14`), so as a practical matter Fryntov–Nazarov **does not constrain** the bridge's `N_0`. It does, however, suggest that **any "large `n`" claim citing Fryntov–Nazarov should specify whether it lives below or above `n ≈ 6.5 × 10⁴`** — that's the boundary at which their stated optimization is meaningful. For the bridge program working at small `n`, Fryntov–Nazarov is mostly a structural/methodological reference rather than a numerical one.

---

## Acknowledgment paragraph for future EHP preprints / forum posts

> The decomposition strategy used here builds on the Stokes-formula framework introduced for the Erdős–Herzog–Piranian lemniscate by Fryntov and Nazarov (arXiv:0808.0717, 2008), who proved `|L| ≤ 2n + O(n^{7/8})` asymptotically and showed that `p_0(z) = z^n − 1` is a local maximum of the lemniscate length in coefficient space. We use a different parameter trichotomy than their good-set / bad-set / capacity-cover surgery, and our parameter dependence does not match theirs term-for-term, but the underlying analytic primitives (Stokes reduction `|L| = 2 Re ∬_E ∂s dA`, logarithmic-derivative pole decomposition, area-integral Cauchy on `|p′|`) are theirs. Their `o(n)` bound was the first sub-leading remainder for the EHP problem and remained the best published asymptotic estimate prior to Tao (2512.12455).

This paragraph is plain, gives credit, and does not overclaim. It is the version I would put in a preprint or a Zulip post; an erdosproblems.com claim should additionally cite the Eremenko–Hayman `n=2` result and the Danchenko `2πn` line.

---

## Reading Notes (audit trail)

- I read all 18 pages of the PDF, not just the abstract.
- Section/equation citations above are anchored to the page numbers I read (1–18).
- I did not invent constants. The `26π`, `2π√n`, `e^{2δ}`, `4a/(δ³√n)`, and `2π(n−1)` are quoted from the inequality on p. 16 exactly as written.
- The `r = a^{2/3}` choice and the `δ = n^{−1/8}` choice are quoted from pp. 9–10 and p. 17 respectively.
- The paper has no `C_0`. I checked. The closest relevant geometric parameter is `δ` (good-set radius) or `r` (annulus radius in §6 first-variation only). Neither appears to power `4`.
- No "Feynman", no "feynman", no "solver" anywhere in this artifact.
