# Zulip Post — #Formal conjectures

**Topic:** `Erdős Problem 114 (EHP conjecture)`

**Post:**

I opened a PR for Erdős Problem 114 — the EHP lemniscate conjecture: google-deepmind/formal-conjectures#3712

Two theorems: `erdos_114` for the full open conjecture (all n), and `erdos_114_small_n` for the computationally certified range 3 ≤ n ≤ 14. The small-n result uses IEEE 1788-rigorous interval arithmetic (branch-and-bound over the parameter space of monic polynomials). Results are on Zenodo with SHA-256 checksums: doi.org/10.5281/zenodo.19322367

This is separate from #3422 — Moritz confirmed a standalone PR works. Used Claude + custom solver for the Lean 4 code. Happy to take feedback on the formalization choices, especially around `arcLength` via Hausdorff measure vs. a path-integral definition.
