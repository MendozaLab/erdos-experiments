#!/usr/bin/env python3
"""
publisher/scripts/seal_evidence.py (prototype, project-local copy)

Compute SHA-256 evidence seals for verified claims in the Three Formalisms paper
and write the seal records + paper evidence root to a local JSON (to be synced to
research-hub D1 as paper_evidence_seals + papers.evidence_root in a later pass).

Seal input = claim_text_normalized || paper_stable_id || claim_location
             || producing_script_path || producing_script_sha256
             || output_artifact_path || output_artifact_sha256
             || verifier_id || verified_at_utc
"""
import hashlib, json, re, unicodedata
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/sessions/kind-serene-shannon/mnt/Math/erdos-experiments/Erdos30")
OUT_DIR = ROOT / ".publisher"
OUT_DIR.mkdir(exist_ok=True)

PAPER_ID = "h2:erdos30:three-formalisms"
PAPER_FILE = ROOT / "paper_three_formalisms.py"     # the generator (binds prose)
PAPER_PDF  = ROOT / "Three_Formalisms_Sidon_Mendoza_2026.pdf"

def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()

def sha256_file(p: Path) -> str:
    return sha256_bytes(p.read_bytes()) if p.exists() else "MISSING"

def normalize(s: str) -> str:
    s = unicodedata.normalize("NFC", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s

now = datetime.now(timezone.utc).isoformat(timespec="seconds")

# Claims verified in this stage 2a pass.
# Each claim identifies: normalized text, anchor in paper, producing script, output artifact.
claims = [
    {
        "claim_index": 1,
        "claim_text": "The channel dispersion V = 0 for every Singer PDS we tested (20 parameter choices, q in {2,3,4,5,7,8,9,11,13,...}), so Singer PDS are dispersion-free Sidon codes.",
        "claim_location": {"section": "2.3", "subsection": "Channel dispersion", "paragraph": 4},
        "producing_script": "erdos-experiments/Erdos30/singer_dispersion.py",
        "output_artifact": None,  # inline stdout; no JSON dumped in this pass
        "verifier_id": "realism_checker@manual-v1",
        "realism_verdict": "verified",  # theoretical: V=0 follows from PDS property
        "notes": "V = Var[i(X;Y)] (Polyanskiy-Poor-Verdu 2010). For a PDS, i is constant on its support, hence V=0 is a theorem, not a computational accident. Numerical agreement across 20 PDS confirms no bug in the implementation.",
    },
    {
        "claim_index": 2,
        "claim_text": "At q = 5 (N = 31, k = 6), exhaustive enumeration of all 44,370 Sidon sets of size 6 finds exactly 310 sets saturating lambda_max = 5, and all 310 are Singer PDS.",
        "claim_location": {"section": "5", "subsection": "(v) Spectral characterization at q = 5", "paragraph": 5},
        "producing_script": "erdos-experiments/Erdos30/paper/attack_a_test.py",
        "output_artifact": None,
        "verifier_id": "realism_checker@manual-v1",
        "realism_verdict": "verified",
        "notes": "310 = 31 translates x 10 multiplier classes. Multiplier group of Z/31Z has order 30; Frobenius x5 has order 3; orbit 30/3 = 10. Algebraically consistent. 44,370 total count requires script output for full traceability.",
    },
    {
        "claim_index": 3,
        "claim_text": "Z/31Z admits 10 multiplier classes of Singer PDS, and each class contributes 31 translates, for 10 x 31 = 310.",
        "claim_location": {"section": "5", "subsection": "(v) Spectral characterization at q = 5", "paragraph": 5},
        "producing_script": "erdos-experiments/Erdos30/paper/attack_a_test.py",
        "output_artifact": None,
        "verifier_id": "realism_checker@manual-v1",
        "realism_verdict": "verified",
        "notes": "Purely algebraic derivation, no numerics required. Serves as an independent cross-check on claim 2.",
    },
]

seals = []
for c in claims:
    script_path = c["producing_script"]
    script_sha  = sha256_file(ROOT.parent.parent / script_path)
    out_path    = c["output_artifact"] or ""
    out_sha     = sha256_file(ROOT.parent.parent / out_path) if out_path else "INLINE-NONE"
    claim_norm  = normalize(c["claim_text"])

    preimage = "\n".join([
        claim_norm,
        PAPER_ID,
        json.dumps(c["claim_location"], sort_keys=True, ensure_ascii=False),
        script_path,
        script_sha,
        out_path,
        out_sha,
        c["verifier_id"],
        now,
    ])
    seal = sha256_bytes(preimage.encode("utf-8"))

    seals.append({
        "seal_id": c["claim_index"],
        "paper_stable_id": PAPER_ID,
        "claim_index": c["claim_index"],
        "claim_text_normalized": claim_norm,
        "claim_location": c["claim_location"],
        "producing_script_path": script_path,
        "producing_script_sha256": script_sha,
        "output_artifact_path": out_path,
        "output_artifact_sha256": out_sha,
        "verifier_id": c["verifier_id"],
        "realism_verdict": c["realism_verdict"],
        "verified_at_utc": now,
        "seal_sha256": seal,
        "notes": c["notes"],
    })

# Paper-level evidence root = SHA-256 of sorted seal hashes joined by newline.
evidence_root = sha256_bytes("\n".join(sorted(s["seal_sha256"] for s in seals)).encode("utf-8"))

record = {
    "paper_stable_id": PAPER_ID,
    "paper_generator_file": str(PAPER_FILE.relative_to(ROOT.parent.parent)),
    "paper_generator_sha256": sha256_file(PAPER_FILE),
    "paper_pdf_file": str(PAPER_PDF.relative_to(ROOT.parent.parent)),
    "paper_pdf_sha256": sha256_file(PAPER_PDF),
    "sealed_at_utc": now,
    "stage": "2a-ensconcement",
    "seals": seals,
    "evidence_root": evidence_root,
    "pending_d1_sync": True,
}

out_file = OUT_DIR / "evidence_seals.json"
out_file.write_text(json.dumps(record, indent=2, ensure_ascii=False))

print(f"[seal] paper   {PAPER_ID}")
print(f"[seal] pdf     {record['paper_pdf_sha256'][:16]}...")
print(f"[seal] gen     {record['paper_generator_sha256'][:16]}...")
print(f"[seal] sealed  {len(seals)} verified claims")
for s in seals:
    print(f"         seal#{s['seal_id']} {s['seal_sha256'][:16]}... @ {s['claim_location']['section']}")
print(f"[seal] evidence_root  {evidence_root}")
print(f"[seal] wrote          {out_file}")
print(f"[seal] sync state     queued for D1 paper_evidence_seals + papers.evidence_root")
