# Red-Team Round — 2026-05-07-r01

**Preprint:** `EXP-MATH-EHP114-SMALL-N-SUBMISSION-READY-PACKET-20260507-12_ERDOSPROBLEMS_POST.md`
**Initialized:** 2026-05-07T12:27:21

This folder is the self-contained audit trail for one round of pre-publication external red-team review. The trail is content-addressable (each subfolder has a `*.sha256` integrity file) and append-only (edits are tracked, not overwritten).

## Reading order for an auditor

Read in numerical order — the trail is designed to be reconstructable cold.

1. **`00_INPUT_SNAPSHOT/`** — verbatim copy of the working artifacts at the moment the round opened. `INPUT_MANIFEST.json` lists every file with its sha256. This is the immutable baseline of "what went in."

2. **`01_EXTERNAL_PROMPT/REDTEAM_PROMPT.md`** — the exact adversarial prompt that was handed to external agents. Self-contained — a reader who has never seen the project can reconstruct what was asked. `PROMPT.sha256` confirms the prompt was not altered after the fact.

3. **`02_EXTERNAL_CRITIQUE/`** — verbatim external responses, one file per source agent (OpenAI, Gemini, external Claude, etc.). `CRITIQUE_MANIFEST.json` gives source attribution, model name, timestamp, optional session URL, and per-file sha256.

4. **`03_EDIT_LOG/`** — append-only chronological log of every edit applied to the working artifacts in response to the critique. Each row has timestamp, plain-English description, and an optional diff snippet under `diffs/`. `EDIT_LOG.sha256` covers the log + every diff file.

5. **`04_FROZEN/`** — the final approved-for-publication artifacts. `FROZEN_MANIFEST.json` records each file's sha256, the approval timestamp, and the inventor name. `FROZEN.lock` is a sealed JSON record covering the manifest hash plus every upstream subfolder's sha256, so any post-freeze tamper is detectable.

## Verification

```bash
python ~/.claude/skills/prepub-redteam/prepub_redteam.py verify <preprint-path>
```

Exits 0 if every file in `04_FROZEN/` matches its hash in `FROZEN.lock` AND the manifest itself is unaltered. Exits 1 on tamper. The Zenodo / arXiv / mendozalab.io deposit scripts call `verify` as a hard pre-flight gate.

## Round reopen / new-round policy

If a critique surfaces an issue late (after `freeze`), open a new round:

```bash
python ~/.claude/skills/prepub-redteam/prepub_redteam.py init <preprint-path> --new-round
```

Each round is immutable after freeze. New rounds get incremented `r02`, `r03`, ... — earlier rounds are preserved verbatim.

INTERNAL — H² research-program private organism.
