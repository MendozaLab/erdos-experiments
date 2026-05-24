# n = 58 Branch Certificate

Date: 2026-04-30
Problem: Erdos #30
Status: DERIVED_N58_BRANCH_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Source

Source exact packet:

```text
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30
```

Derived packet:

```text
EXP-MM-030-PMF-N58-BRANCH-CERTIFICATE-2026-04-30
```

Verification script:

```text
erdos-experiments/Erdos30/scripts/n58_branch_certificate.py
```

SHA256 verification returned `OK`.

## Question

Is `n = 58` still translated-embedding behavior, or is it the first true
difference-skeleton branch row?

## Answer

Both, with an important split.

`n = 58` is the first local row with two difference skeletons, but the new
skeleton is not the Pareto face. The Pareto face remains on the original
translated-embedding skeleton.

## Skeleton Roles

| skeleton | witnesses | persists from n=57 | prefix | mass | joint | Pareto |
|---:|---|---|---|---|---|---|
| 0 | `[0,1,4,5,6,7,8,9]` | `[0,1,4,5,6,7]` | `[8,9]` | `[7]` | `[7]` | `[7,9]` |
| 1 | `[2,3]` | `[]` | `[2]` | `[]` | `[]` | `[]` |

So the second skeleton enters through prefix witness `2`, not through the
mass/joint/Pareto exposed handoff.

## Pareto Candidates

`n = 58` Pareto candidate `7`:

```text
[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]
```

`n = 58` Pareto candidate `9`:

```text
[3, 5, 17, 24, 32, 35, 48, 52, 57, 58]
```

Certified relation:

```text
index 9 = index 7 + 1
```

This relation holds as an integer translation and also mod `57` and mod `58`.

Persistence from `n = 57`:

```text
n=58 index 7 = n=57 index 5
n=58 index 9 = n=57 index 5 + 1
n=58 index 9 = n=57 index 3 + 2
```

## Three-Stage Mechanism

The local mechanism is now:

```text
n = 56:
  single exposed embedding

n = 57:
  translated-embedding handoff inside one Singer-modulus skeleton

n = 58:
  first local skeleton branch appears, but Pareto remains on translated chain
```

That is more precise than "57 is a skeleton split." It is not. The skeleton
split appears at `58`, and even there it enters through the prefix side before
it becomes a Pareto branch.

## Checks

All checks passed:

- all `n = 58` exported witnesses are interval Sidon;
- `n = 58` has two difference skeletons;
- Pareto candidates `7` and `9` share a skeleton;
- `index 9 = index 7 + 1`;
- `index 7` persists from `n = 57` index `5`;
- `index 9 = n = 57` index `5 + 1`;
- a new skeleton branch exists at `n = 58`;
- the new branch is prefix-side, not Pareto-side.

## Meaning

This stabilizes the finite proof-candidate object:

> the exact face first shows field-sensitive translated embedding at `57`; one
> step later, at `58`, a second difference skeleton appears, but the exposed
> Pareto face still rides the original translated chain.

This is now a good finite package for Lean certificate work.

## Claim Boundary

Safe:

- `n = 58 first local skeleton branch`
- `new branch is prefix-side, not Pareto-side`
- `Pareto remains on the translated embedding chain`
- `three-stage local finite mechanism`

Unsafe:

- `n = 58 proves a phase transition`
- `PMF proves Sidon`
- `physics solves Erdos #30`
- `SOTA theorem result`
