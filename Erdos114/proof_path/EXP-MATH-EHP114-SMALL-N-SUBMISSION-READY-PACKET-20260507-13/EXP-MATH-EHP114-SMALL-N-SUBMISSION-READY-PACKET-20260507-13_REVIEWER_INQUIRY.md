# Reviewer inquiry — what we are asking for and not asking for

This packet asks reviewers two things and explicitly does not ask a third.

## Asking

1. (erdosproblems.com) Is a finite-frontier note of this form and length the right way to record the small-`n` side on Problem 114's page, given that the all-degree statement remains open and is not being claimed here?

2. (formal-conjectures) Does the proposed shape — keep `Erdos114.erdos_114` open; add `Erdos114.erdos_114_finite_le_14` as a separately named theorem; record each per-degree dependency as an explicit named axiom — match upstream conventions, or should certificates be packaged differently?

## Not asking

We are not asking reviewers to certify or accept the all-degree Erdős #114 conjecture, to mark the main theorem closed, or to extract a threshold from Tao's sufficiently-large-`n` theorem (arXiv:2512.12455). The finite variant and the all-degree conjecture are deliberately separate.

## What an external reviewer can verify in one sitting

- Each row of the certificate table links to a public JSON and a public SHA-256 sidecar.
- The SHA-256 sidecars match the row hashes printed in the table.
- The Eremenko–Hayman pin (`d = 2` extremal case) is in the abstract and the proof remark after Lemma 5 of arXiv:0805.2295.
- The MacLane historical citation resolves at doi:10.1307/mmj/1028989918.
- The n=13 row's re-run is reproducible: `cargo build --release --bin ehp_general_ieee1788` then `./target/release/ehp_general_ieee1788 --outdir <path> 13 13` against the patched binary at HEAD of `triage-recovered-2026-05-03` (commit `dae62b8` or later) reproduces `bb_total_evals = 197,132,288`, `EHP_N13_PROVEN`, SHA `c06c633b…`. Wall-clock varies with hardware (53 minutes on a Modal 32-CPU worker).

## What an external reviewer cannot do from this packet alone

Verify the soundness of the Rust + `inari` interval-arithmetic implementation. The IEEE-1788 standard is taken for granted; the `inari` crate's correctness is an external dependency. A future iteration of this work should either include a Lean formalization of the per-degree certificates or rely on a separately certified interval-arithmetic toolchain.

## Tooling disclosure

This inquiry was prepared with AI-assisted tooling. The questions are the author's responsibility.
