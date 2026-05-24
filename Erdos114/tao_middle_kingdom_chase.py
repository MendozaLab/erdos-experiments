"""Build a bounded attack packet for the EHP #114 middle range.

This script fetches the official arXiv source for Tao's high-degree EHP paper,
extracts the dependency landmarks that matter for constant-chasing, and writes
the experiment-contract artifacts:

  EXP-MATH-EHP114-TAO-MIDDLE-KINGDOM-20260505-01_RESULTS.json
  EXP-MATH-EHP114-TAO-MIDDLE-KINGDOM-20260505-01_REPORT.md
  EXP-MATH-EHP114-TAO-MIDDLE-KINGDOM-20260505-01_RESULTS.sha256

It does not claim a proof. It creates the control plane for shrinking the
finite gap between the DOI-backed n <= 14 certificate and Tao's unspecified
large-n threshold.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import tarfile
import urllib.request
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable


EXPERIMENT_ID = "EXP-MATH-EHP114-TAO-MIDDLE-KINGDOM-20260505-01"
ARXIV_ID = "2512.12455"
ARXIV_EPRINT_URL = f"https://arxiv.org/e-print/{ARXIV_ID}"
ARXIV_ABS_URL = f"https://arxiv.org/abs/{ARXIV_ID}"
ZENODO_DOI = "10.5281/zenodo.19480329"


@dataclass
class Landmark:
    label: str
    kind: str
    line_start: int
    line_end: int
    role: str
    extraction: str
    middle_gap_use: str


def repo_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def output_dir() -> Path:
    return Path(__file__).resolve().parent


def fetch_arxiv_source() -> bytes:
    req = urllib.request.Request(
        ARXIV_EPRINT_URL,
        headers={"User-Agent": "MendozaLab-EHP114-audit/1.0"},
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        return response.read()


def extract_lemniscate_tex(payload: bytes) -> str:
    # arXiv e-print payload is a gzip-compressed tar archive for this record.
    data = gzip.decompress(payload)
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:") as archive:
        member = archive.extractfile("lemniscate.tex")
        if member is None:
            raise RuntimeError("lemniscate.tex not found in arXiv source archive")
        return member.read().decode("utf-8", errors="replace")


def find_env(lines: list[str], label: str) -> tuple[int, int, str, str]:
    label_pat = re.compile(r"\\label\{" + re.escape(label) + r"\}")
    label_line = None
    for idx, line in enumerate(lines):
        if label_pat.search(line):
            label_line = idx
            break
    if label_line is None:
        raise RuntimeError(f"label not found: {label}")

    start = label_line
    while start >= 0 and "\\begin{" not in lines[start]:
        start -= 1
    if start < 0:
        start = label_line

    begin_match = re.search(r"\\begin\{([^}]+)\}(?:\[([^\]]*)\])?", lines[start])
    env = begin_match.group(1) if begin_match else "unknown"
    title = begin_match.group(2) if begin_match and begin_match.group(2) else ""

    end = label_line
    end_token = f"\\end{{{env}}}"
    while end < len(lines) and end_token not in lines[end]:
        end += 1
    if end >= len(lines):
        end = label_line

    return start + 1, end + 1, env, title


def snippet(lines: list[str], start: int, end: int, max_chars: int = 420) -> str:
    text = " ".join(line.strip() for line in lines[start - 1 : end])
    text = re.sub(r"\s+", " ", text)
    return text[:max_chars]


def build_landmarks(lines: list[str]) -> list[Landmark]:
    specs = [
        (
            "main-thm",
            "Tao high-degree ceiling",
            "States the asymptotic theorem: EHP holds for sufficiently large n.",
            "Defines the upper end of the finite middle range, but no practical N0.",
        ),
        (
            "erem",
            "Extremizer normalization",
            "Reduces the problem to normalized maximizers with connected lemniscates.",
            "Required if the finite certificate is to match Tao's proof variables.",
        ),
        (
            "stokes",
            "Area-to-length functional",
            "Represents lemniscate length through the X1..X5 error budget.",
            "This is the natural bridge for an interval-compressed certificate.",
        ),
        (
            "out",
            "p0 reference value",
            "Gives the exact length formula and asymptotics for z^n - 1.",
            "This is where DOI L* intervals and Tao's asymptotic reference meet.",
        ),
        (
            "geomcontrol",
            "First geometry control",
            "K controls near-roundness, dispersion, and derivative decay.",
            "First constant-chase choke point: every hidden constant enters later.",
        ),
        (
            "x2-lem",
            "Small-field error term",
            "Bounds X2 by mu times dispersion.",
            "Quantifying the 'mu sufficiently small' choice is a tractable first chase.",
        ),
        (
            "x3-lem",
            "Large-psi error term",
            "Bounds X3 by an explicit expression involving n, dispersion, mu.",
            "Pairs with X2 to optimize mu explicitly instead of asymptotically.",
        ),
        (
            "lemni",
            "Total size bound",
            "Upgrades dispersion and origin repulsion to O(1).",
            "Separates global shape control from final local deficit.",
        ),
        (
            "inside",
            "Inner annulus estimate",
            "Controls the small-radius part of the lemniscate.",
            "Ancestor of the radial Puiseux interval target.",
        ),
        (
            "annulus",
            "Intermediate region estimate",
            "Controls the main bulk annulus using zero-counting and local bounds.",
            "Middle-region term must be made quantitative for any finite splice.",
        ),
        (
            "outside",
            "Outer tip estimate",
            "Produces the 4 log 2 boundary contribution.",
            "This is where endpoint/tip constants enter the final comparison.",
        ),
        (
            "ets",
            "Critical-point collapse",
            "Improves dispersion to o(1).",
            "This is the first final-section gate toward uniqueness.",
        ),
        (
            "pots",
            "Total-size collapse",
            "Improves total size to o(1).",
            "This makes the final split local around p0.",
        ),
        (
            "inside-2",
            "Final inner deficit",
            "Adds a negative deficit term in the inner region.",
            "Closest analytic cousin of the radial Puiseux certificate.",
        ),
        (
            "annulus-2",
            "Final middle deficit",
            "Adds a negative deficit term over the intermediate annulus.",
            "The shape-cone interval bound should target this term.",
        ),
        (
            "outside-again",
            "Final outer deficit",
            "Compares the outer annulus directly against p0.",
            "Remainder absorption must keep this error below the inner/middle gain.",
        ),
    ]
    out: list[Landmark] = []
    for label, role, extraction, use in specs:
        start, end, env, title = find_env(lines, label)
        out.append(
            Landmark(
                label=label,
                kind=f"{env}{('[' + title + ']') if title else ''}",
                line_start=start,
                line_end=end,
                role=role,
                extraction=extraction,
                middle_gap_use=use,
            )
        )
    return out


def find_line_numbers(lines: list[str], needles: Iterable[str]) -> dict[str, list[int]]:
    result: dict[str, list[int]] = {}
    for needle in needles:
        result[needle] = [i + 1 for i, line in enumerate(lines) if needle in line]
    return result


def write_report(results: dict, landmarks: list[Landmark], tex_lines: list[str]) -> str:
    md: list[str] = []
    md.append("# EHP114 Tao Middle-Kingdom Closure Packet")
    md.append("")
    md.append(f"Experiment: `{EXPERIMENT_ID}`")
    md.append("")
    md.append("## Meaning")
    md.append("")
    md.append(
        "This packet turns the gap between the DOI-backed n <= 14 certificate "
        "and Tao's sufficiently-large-n theorem into a concrete attack surface."
    )
    md.append("")
    md.append(
        "The target is not brute-force computation past n = 20. The target is "
        "to make Tao's asymptotic proof quantitative enough to splice with "
        "Rust/inari interval certificates, while replacing global B&B with "
        "radial, shape-cone, and remainder interval lemmas."
    )
    md.append("")
    md.append("Preserve the book-facing phrase: shadow signature, not universal law.")
    md.append("")
    md.append("## Source")
    md.append("")
    md.append(f"- arXiv source: `{ARXIV_EPRINT_URL}`")
    md.append(f"- arXiv abstract: `{ARXIV_ABS_URL}`")
    md.append(f"- source SHA-256: `{results['source_sha256']}`")
    md.append(f"- local DOI anchor: `{ZENODO_DOI}`")
    md.append("")
    md.append("## What Tao Gives")
    md.append("")
    md.append(
        "Tao proves EHP for sufficiently large degree. The paper states that all "
        "implied constants are effectively computable, so the remaining problem "
        "is finite, but the paper does not optimize the resulting numerical bound."
    )
    md.append("")
    md.append("The practical consequence is a middle range:")
    md.append("")
    md.append("```text")
    md.append("n = 1,2        known / classical")
    md.append("n = 3..14      DOI-backed Rust/inari IEEE-1788 certificate")
    md.append("n = 15..20     plausible finite-certificate extension, not yet safe")
    md.append("n = 21..N0-1   middle kingdom")
    md.append("n >= N0        Tao, after explicit constant extraction")
    md.append("```")
    md.append("")
    md.append("## Dependency Landmarks")
    md.append("")
    md.append("| label | source lines | role | middle-gap use |")
    md.append("|---|---:|---|---|")
    for item in landmarks:
        md.append(
            f"| `{item.label}` | {item.line_start}-{item.line_end} | "
            f"{item.role} | {item.middle_gap_use} |"
        )
    md.append("")
    md.append("## The Attack Order")
    md.append("")
    md.append("1. **Constant chase Tao's final proof.**")
    md.append(
        "   Start at `inside-2`, `annulus-2`, and `outside-again`, then walk "
        "backward through `pots`, `ets`, `inside`, `annulus`, `outside`, and "
        "`geomcontrol`. The goal is not a beautiful bound; the first goal is any "
        "explicit N0."
    )
    md.append("")
    md.append("2. **Interval-harden the radial term first.**")
    md.append(
        "   This corresponds to the inner-deficit lane. It should reuse the "
        "n = 14 DOI certificate as calibration and target a fixed-n theorem for "
        "the radial hypergeometric/Puiseux deficit."
    )
    md.append("")
    md.append("3. **Then harden the shape cone.**")
    md.append(
        "   This corresponds to the intermediate-region deficit. The theorem "
        "should prove that nonradial perturbations cannot restore the lost "
        "length after the radial term is controlled."
    )
    md.append("")
    md.append("4. **Then absorb the mixed remainder.**")
    md.append(
        "   This corresponds to keeping the outer-region error and local mixed "
        "terms below the inner plus middle gains."
    )
    md.append("")
    md.append("5. **Only then extend finite n.**")
    md.append(
        "   Reconcile the n = 15 and n = 16 zero-eval artifacts before claiming "
        "them. Treat n = 20 as a diagnostic endpoint until a new versioned "
        "Rust/inari certificate exists."
    )
    md.append("")
    md.append("## First Concrete Theorem Target")
    md.append("")
    md.append("```lean")
    md.append("theorem ehp114_fixed_n_radial_puiseux_interval")
    md.append("    (n : Nat) (eps : Real)")
    md.append("    (hn : n = 14)")
    md.append("    (hpos : 0 < eps)")
    md.append("    (hsmall : eps <= (1 : Real) / 10000) :")
    md.append("    C14 * Real.rpow eps ((1 : Real) / 14)")
    md.append("      <= D14 (radialMode14 eps) := by")
    md.append("  -- hypergeometric connection formula + interval constants")
    md.append("  sorry")
    md.append("```")
    md.append("")
    md.append("This is intentionally fixed-n. If n = 14 closes, general n becomes a")
    md.append("parameterization problem. If n = 14 does not close, n = 20 is premature.")
    md.append("")
    md.append("## Claim Ceiling")
    md.append("")
    md.append("Safe: `finite certified computation n <= 14`, `Tao proves sufficiently large n`,")
    md.append("and `we now have a concrete middle-range closure map`.")
    md.append("")
    md.append("Unsafe: saying the middle range is closed, saying n = 20 is certified, or")
    md.append("claiming a Lean/formal proof of the analytic bridge.")
    md.append("")
    md.append("## Verification")
    md.append("")
    md.append("- arXiv source fetched and parsed.")
    md.append("- `lemniscate.tex` landmarks extracted by label.")
    md.append("- No scorecard, D1, public deployment, or Lean status was mutated.")
    md.append("")
    return "\n".join(md) + "\n"


def main() -> None:
    out_dir = output_dir()
    payload = fetch_arxiv_source()
    tex = extract_lemniscate_tex(payload)
    tex_lines = tex.splitlines()
    landmarks = build_landmarks(tex_lines)

    needle_lines = find_line_numbers(
        tex_lines,
        [
            "sufficiently large",
            "effectively computable",
            "no attempt to optimize",
            "small enough",
            "large enough",
            "absolute constant",
        ],
    )

    results = {
        "experiment_id": EXPERIMENT_ID,
        "date": "2026-05-05",
        "scope": "Tao large-n constant-chase map plus finite-certificate splice plan",
        "claim_ceiling": (
            "This is a closure control packet, not a proof. It identifies the "
            "middle-range bottlenecks between n <= 14 Rust/inari certificates "
            "and Tao's sufficiently-large-n theorem."
        ),
        "arxiv_id": ARXIV_ID,
        "arxiv_abs_url": ARXIV_ABS_URL,
        "arxiv_eprint_url": ARXIV_EPRINT_URL,
        "source_sha256": hashlib.sha256(tex.encode("utf-8")).hexdigest(),
        "source_line_count": len(tex_lines),
        "needle_lines": needle_lines,
        "doi_anchor": {
            "doi": ZENODO_DOI,
            "certified_range": "n=3..14",
            "rigor": "ieee_1788_interval_arithmetic_inari",
        },
        "landmarks": [asdict(item) for item in landmarks],
        "attack_order": [
            "constant_chase_final_section_inside2_annulus2_outside_again",
            "radial_puiseux_interval_hardening_at_n14",
            "shape_cone_interval_hardening",
            "mixed_remainder_absorption",
            "reconcile_n15_n16_zero_eval_records",
            "versioned_rust_inari_extension_feasibility_n17_n20",
        ],
        "first_theorem_target": {
            "name": "ehp114_fixed_n_radial_puiseux_interval",
            "scope": "fixed n=14 radial deficit lower bound",
            "reason": "n=14 is the DOI-backed material stress case; if n=14 fails, n=20 is premature.",
        },
        "status": "CONTROL_PACKET_READY",
    }

    results_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.json"
    report_path = out_dir / f"{EXPERIMENT_ID}_REPORT.md"
    sha_path = out_dir / f"{EXPERIMENT_ID}_RESULTS.sha256"

    results_json = json.dumps(results, indent=2, sort_keys=True)
    results_path.write_text(results_json + "\n")
    report_path.write_text(write_report(results, landmarks, tex_lines))
    sha_path.write_text(hashlib.sha256(results_json.encode("utf-8")).hexdigest() + f"  {results_path.name}\n")

    print(results_path)
    print(report_path)
    print(sha_path)


if __name__ == "__main__":
    main()

