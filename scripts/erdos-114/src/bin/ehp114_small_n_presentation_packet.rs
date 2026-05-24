//! EHP #114 small-n submission-ready packet for erdosproblems/formal-conjectures.
//!
//! This emitter packages the finite n<15 theorem for review while keeping the
//! all-degree Tao bridge conditional. It does not recompute the finite interval
//! certificates and does not promote any full EHP114 claim.

#![recursion_limit = "256"]

use serde::Serialize;
use serde_json::{json, Value};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const EXPERIMENT_ID: &str = "EXP-MATH-EHP114-SMALL-N-SUBMISSION-READY-PACKET-20260507-11";
const FINITE_PACKET_ID: &str = "EXP-MATH-EHP114-FINITE-N-LESS-15-PROOF-PACKET-20260507-01";
const ZENODO_DOI: &str = "10.5281/zenodo.19480329";
const ZENODO_URL: &str = "https://zenodo.org/records/19480329";
const GITHUB_REPO_URL: &str = "https://github.com/MendozaLab/erdos-experiments";
const ERDOS114_URL: &str = "https://www.erdosproblems.com/114";
const FORMAL_CONJECTURES_URL: &str = "https://github.com/google-deepmind/formal-conjectures";
const TAO_ARXIV_URL: &str = "https://arxiv.org/abs/2512.12455";
const EREMENKO_HAYMAN_ARXIV_URL: &str = "https://arxiv.org/abs/0805.2295";
const EREMENKO_HAYMAN_PDF_URL: &str = "https://www.math.purdue.edu/~eremenko/dvi/erdos23.pdf";
const MACLANE_DOI_URL: &str = "https://doi.org/10.1307/mmj/1028989918";
const MACLANE_CITATION: &str =
    "G. R. MacLane, On a conjecture of Erdos, Herzog, and Piranian, Michigan Math. J. 2 (1953/54), 147-148, doi:10.1307/mmj/1028989918.";
const TAO_NON_BRIDGING_SENTENCE: &str =
    "No explicit threshold is extracted from Tao's sufficiently-large-n theorem in this finite packet, and no degree n >= 15 is claimed here.";
const GITHUB_RESULTS_DIR_URL: &str =
    "https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114";
const GITHUB_RAW_RESULTS_DIR_URL: &str =
    "https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114";
const ZENODO_FILE_URL_PREFIX: &str = "https://zenodo.org/records/19480329/files";

const THEOREM_TITLE: &str =
    "Finite verification of the Erdos-Herzog-Piranian lemniscate conjecture for n < 15";
const THEOREM_STATEMENT: &str = "For every monic complex polynomial p of degree n with 1 <= n < 15, the lemniscate length is at most the length for p(z)=z^n-1.";
const PUBLIC_ABSTRACT: &str = "I have assembled a finite-degree verification of the Erdos-Herzog-Piranian lemniscate length conjecture for 1 <= n < 15. The dependency split is: n=1 is elementary, n=2 is the Eremenko-Hayman degree-two case, and 3 <= n <= 14 are covered by reproducible interval certificates with row-level hashes and public certificate links. This is meant as a finite-frontier result, complementary to Tao's sufficiently-large-n theorem, not as a standalone proof of the all-degree conjecture.";
const FLIP_SIDE_SENTENCE: &str = "Tao proves the conjecture for all sufficiently large degrees. This packet verifies the opposite finite frontier through degree fourteen. The remaining bridge question is whether Tao's threshold can be made effective at or below n=15, or whether a finite middle range remains.";
const BRIDGE_APPENDIX_QUESTION: &str = "Is Tao's sufficiently-large threshold effective in a range that begins at or before n=15? If so, the present finite packet and Tao's theorem would combine into a full proof of Erdos #114.";
const REVIEWER_INQUIRY_LANGUAGE: &str = "I am not claiming the all-degree Erdos #114 statement here. The finite side covers n <= 14; Tao's theorem covers all sufficiently large n. Is this the right way to record the finite side and phrase the remaining bridge question: extract Tao's threshold, and if it is above 15, certify the intervening finite degrees?";
const CLAIM_CEILING: &str = "Submission preflight packet only. It prepares erdosproblems.com and Google DeepMind formal-conjectures wording for the finite n<15 result and the conditional Tao-threshold bridge; it does not assert the all-degree Erdos #114 statement and does not authorize n=15 or higher computation.";
const AI_TOOLING_DISCLOSURE: &str = "Tooling disclosure: this packet and draft wording were prepared with AI-assisted tooling. The mathematical claim rests only on the cited literature rows and SHA-checked certificate artifacts; no AI output is used as proof evidence. The author takes responsibility for the mathematical claims, certificate selection, and submission wording.";

#[derive(Clone, Debug)]
struct Config {
    outdir: PathBuf,
    finite_packet: PathBuf,
}

#[derive(Clone, Debug, Serialize)]
struct ShaCheck {
    role: String,
    path: String,
    sha_path: String,
    expected_sha256: String,
    actual_sha256: String,
    sha_status: String,
}

fn parse_args() -> Result<Config, String> {
    let mut cfg = Config {
        outdir: PathBuf::from(format!("../../Erdos114/proof_path/{EXPERIMENT_ID}")),
        finite_packet: PathBuf::from(format!(
            "../../Erdos114/proof_path/{FINITE_PACKET_ID}/{FINITE_PACKET_ID}_RESULTS.json"
        )),
    };

    let args: Vec<String> = env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--outdir" => {
                cfg.outdir = next_path(&args, i, "--outdir")?;
                i += 2;
            }
            "--finite-packet" => {
                cfg.finite_packet = next_path(&args, i, "--finite-packet")?;
                i += 2;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
    }
    Ok(cfg)
}

fn next_path(args: &[String], i: usize, flag: &str) -> Result<PathBuf, String> {
    if i + 1 >= args.len() {
        return Err(format!("{flag} requires a path"));
    }
    Ok(PathBuf::from(&args[i + 1]))
}

fn read_json(path: &Path) -> Result<Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    serde_json::from_str(&text)
        .map_err(|err| format!("failed to parse {} as JSON: {err}", path.display()))
}

fn bool_field(v: &Value, key: &str) -> Result<bool, String> {
    v.get(key)
        .and_then(Value::as_bool)
        .ok_or_else(|| format!("missing bool field {key}"))
}

fn usize_field(v: &Value, key: &str) -> Result<usize, String> {
    let n = v
        .get(key)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("missing integer field {key}"))?;
    Ok(n as usize)
}

fn str_field<'a>(v: &'a Value, key: &str) -> Result<&'a str, String> {
    v.get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("missing string field {key}"))
}

fn sha256_file(path: &Path) -> Result<String, String> {
    let bytes = fs::read(path)
        .map_err(|err| format!("failed to read {} for sha256: {err}", path.display()))?;
    Ok(sha256::digest(bytes))
}

fn sha_path_for(path: &Path) -> Result<PathBuf, String> {
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| format!("invalid file name in {}", path.display()))?;
    if !file_name.ends_with("_RESULTS.json") {
        return Err(format!(
            "result JSON must end in _RESULTS.json: {}",
            path.display()
        ));
    }
    Ok(path.with_file_name(file_name.replace("_RESULTS.json", "_RESULTS.sha256")))
}

fn verify_result_sha(role: &str, path: &Path) -> Result<ShaCheck, String> {
    let sha_path = sha_path_for(path)?;
    let sha_text = fs::read_to_string(&sha_path)
        .map_err(|err| format!("failed to read {}: {err}", sha_path.display()))?;
    let expected = sha_text
        .split_whitespace()
        .next()
        .ok_or_else(|| format!("empty sha file {}", sha_path.display()))?
        .to_string();
    let actual = sha256_file(path)?;
    let sha_status = if expected == actual { "PASS" } else { "FAIL" };
    Ok(ShaCheck {
        role: role.to_string(),
        path: path.display().to_string(),
        sha_path: sha_path.display().to_string(),
        expected_sha256: expected,
        actual_sha256: actual,
        sha_status: sha_status.to_string(),
    })
}

fn unix_timestamp_string() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs()
        .to_string()
}

fn artifact_paths(
    outdir: &Path,
) -> (
    PathBuf,
    PathBuf,
    PathBuf,
    PathBuf,
    PathBuf,
    PathBuf,
    PathBuf,
    PathBuf,
    PathBuf,
) {
    (
        outdir.join(format!("{EXPERIMENT_ID}_RESULTS.json")),
        outdir.join(format!("{EXPERIMENT_ID}_REPORT.md")),
        outdir.join(format!("{EXPERIMENT_ID}_README.md")),
        outdir.join(format!("{EXPERIMENT_ID}_ERDOSPROBLEMS_POST.md")),
        outdir.join(format!("{EXPERIMENT_ID}_FORMAL_CONJECTURES_PACKET.md")),
        outdir.join(format!("{EXPERIMENT_ID}_FORMAL_CONJECTURES_ISSUE_DRAFT.md")),
        outdir.join(format!("{EXPERIMENT_ID}_REVIEWER_INQUIRY.md")),
        outdir.join(format!("{EXPERIMENT_ID}_SHA_MANIFEST.md")),
        outdir.join(format!("{EXPERIMENT_ID}_RESULTS.sha256")),
    )
}

fn contains_forbidden_all_degree_claim(text: &str) -> bool {
    let lowered = text.to_ascii_lowercase();
    let forbidden = [
        "we prove erdos #114",
        "we prove the erdos-herzog-piranian conjecture",
        "we have proved erdos #114",
        "solves erdos #114",
        "solution of erdos #114",
        "all-degree proof",
        "completed erdos #114",
    ];
    forbidden.iter().any(|needle| lowered.contains(needle))
}

fn public_result_url(file_name: &str) -> String {
    format!("{GITHUB_RESULTS_DIR_URL}/{file_name}")
}

fn public_raw_result_url(file_name: &str) -> String {
    format!("{GITHUB_RAW_RESULTS_DIR_URL}/{file_name}")
}

fn zenodo_file_url(file_name: &str) -> String {
    format!("{ZENODO_FILE_URL_PREFIX}/{file_name}")
}

fn degree_rows(finite_packet: &Value) -> Result<&Vec<Value>, String> {
    finite_packet
        .get("degree_rows")
        .and_then(Value::as_array)
        .ok_or_else(|| "finite packet missing degree_rows array".to_string())
}

fn public_certificate_rows(finite_packet: &Value) -> Result<Vec<Value>, String> {
    let mut rows = Vec::new();
    for row in degree_rows(finite_packet)? {
        let degree = row
            .get("degree")
            .and_then(Value::as_u64)
            .ok_or_else(|| "degree row missing degree".to_string())?;
        if degree == 1 {
            rows.push(json!({
                "degree": 1,
                "source_kind": "DIRECT_ANALYTIC",
                "lemma_name": row["lemma_name"].clone(),
                "status": row["certificate_status"].clone(),
                "public_source": "Monic linear case: {|z+a|=1} is a translated unit circle of length 2*pi.",
                "row_hash": null,
                "result_json_url": null,
                "sha256_url": null,
                "reviewer_note": "No machine artifact is required for n=1."
            }));
        } else if degree == 2 {
            rows.push(json!({
                "degree": 2,
                "source_kind": "LITERATURE",
                "lemma_name": row["lemma_name"].clone(),
                "status": row["certificate_status"].clone(),
                "historical_citation": MACLANE_CITATION,
                "accessible_citation": "Alexandre Eremenko and Walter K. Hayman, On the length of lemniscates, Michigan Math. J. 46 (1999), 409-415; arXiv:0805.2295.",
                "source_urls": {
                    "maclane_doi": MACLANE_DOI_URL,
                    "arxiv": EREMENKO_HAYMAN_ARXIV_URL,
                    "purdue_pdf": EREMENKO_HAYMAN_PDF_URL
                },
                "source_pin": "MacLane is cited as the earlier degree-two literature row. The accessible Eremenko-Hayman source states in the abstract that for d=2 the extremal level set is Bernoulli's lemniscate. In the proof text, the remark after Lemma 5 says that z^2+1 is extremal for d=2; z^2+1 and z^2-1 are equivalent for length by rotation.",
                "row_hash": null,
                "result_json_url": null,
                "sha256_url": null,
                "reviewer_note": "Pinned literature row; no local computational artifact is used for n=2."
            }));
        } else {
            let experiment_id = row
                .get("experiment_id")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("degree {degree} missing experiment_id"))?;
            let result_file = format!("{experiment_id}_RESULTS.json");
            let sha_file = format!("{experiment_id}_RESULTS.sha256");
            rows.push(json!({
                "degree": degree,
                "source_kind": "DOI_BACKED_RUST_INARI_IEEE1788_INTERVAL_CERTIFICATE",
                "lemma_name": row["lemma_name"].clone(),
                "status": row["certificate_status"].clone(),
                "experiment_id": experiment_id,
                "verdict": row["verdict"].clone(),
                "rigor": row["rigor"].clone(),
                "row_hash": row["artifact_sha256"].clone(),
                "artifact_sha_status": row["artifact_sha_status"].clone(),
                "result_json_url": {
                    "zenodo_file": zenodo_file_url(&result_file),
                    "github_blob": public_result_url(&result_file),
                    "github_raw": public_raw_result_url(&result_file)
                },
                "sha256_url": {
                    "zenodo_file": zenodo_file_url(&sha_file),
                    "github_blob": public_result_url(&sha_file),
                    "github_raw": public_raw_result_url(&sha_file)
                },
                "l_star_lower": row["l_star_lower"].clone(),
                "l_star_upper": row["l_star_upper"].clone(),
                "bb_total_evals": row["bb_total_evals"].clone(),
                "bb_level_count": row["bb_level_count"].clone(),
                "route_exception": row["route_exception"].clone(),
                "reviewer_note": "Row-level public result JSON and SHA-256 sidecar are linked from Zenodo and GitHub."
            }));
        }
    }
    Ok(rows)
}

fn reviewer_grade_certificate_table(rows: &[Value]) -> String {
    let mut out = String::from(
        "| n | dependency | public source | SHA / pin | status |\n\
|---:|---|---|---|---|\n",
    );
    for row in rows {
        let degree = row["degree"].as_u64().unwrap_or_default();
        if degree == 1 {
            out.push_str(
                "| 1 | direct analytic | translated unit circle | not applicable | covered |\n",
            );
        } else if degree == 2 {
            out.push_str(&format!(
                "| 2 | MacLane / Eremenko-Hayman literature row | [MacLane DOI]({}) / [Eremenko-Hayman arXiv]({}) / [PDF]({}) | {}; Eremenko-Hayman accessible pin: Michigan Math. J. 46 (1999), 409-415; degree-two extremal case pinned to Bernoulli lemniscate | covered |\n",
                MACLANE_DOI_URL, EREMENKO_HAYMAN_ARXIV_URL, EREMENKO_HAYMAN_PDF_URL, MACLANE_CITATION
            ));
        } else {
            let experiment_id = row["experiment_id"].as_str().unwrap_or("UNKNOWN");
            let result_file = format!("{experiment_id}_RESULTS.json");
            let sha_file = format!("{experiment_id}_RESULTS.sha256");
            let hash = row["row_hash"].as_str().unwrap_or("UNKNOWN");
            let route_note = if row["route_exception"].is_null() {
                ""
            } else {
                " route exception flagged"
            };
            out.push_str(&format!(
                "| {} | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON]({}) / [raw]({}) / [Zenodo file]({}) | `{}`; [sha sidecar]({}) | {}{} |\n",
                degree,
                public_result_url(&result_file),
                public_raw_result_url(&result_file),
                zenodo_file_url(&result_file),
                hash,
                public_result_url(&sha_file),
                row["status"].as_str().unwrap_or("UNKNOWN"),
                route_note
            ));
        }
    }
    out
}

fn portable_readme_markdown(result: &Value) -> String {
    let table = result["reviewer_grade_certificate_table"]
        .as_str()
        .unwrap_or("");
    format!(
        "# EHP114 finite n < 15 submission preflight packet\n\n\
## Finite theorem\n\n\
{}\n\n\
## Public references\n\n\
- Zenodo certificate record: [{}]({})\n\
- Code and results repository: [{}]({})\n\
- Erdős Problems #114 page: [{}]({})\n\
- Google DeepMind formal-conjectures: [{}]({})\n\
- Tao large-degree paper: [{}]({})\n\n\
## Dependency split\n\n\
- `n=1`: direct analytic linear case.\n\
- `n=2`: MacLane / Eremenko-Hayman degree-two literature row. MacLane is the historical citation ([doi:10.1307/mmj/1028989918]({})); Eremenko-Hayman is the accessible source pin at Michigan Math. J. 46 (1999), 409-415 / arXiv:0805.2295.\n\
- `3 <= n <= 14`: DOI-backed Rust/inari IEEE-1788 interval-certificate rows with row-level result links and SHA-256 hashes.\n\
- Tao bridge status: {}\n\n\
## Reviewer-grade certificate table\n\n\
{}\n\
## Local preflight status\n\n\
- Finite packet status: `{}`\n\
- Covered degrees: `{}`\n\
- Missing degrees: `{}`\n\
- Source SHA failures: `{}`\n\
- Full EHP114 claim: `{}`\n\
- External submission authorized: `{}`\n\n\
## Claim ceiling\n\n\
{}\n\n\
## Tooling disclosure\n\n\
{}\n",
        THEOREM_STATEMENT,
        ZENODO_DOI,
        ZENODO_URL,
        GITHUB_REPO_URL,
        GITHUB_REPO_URL,
        ERDOS114_URL,
        ERDOS114_URL,
        FORMAL_CONJECTURES_URL,
        FORMAL_CONJECTURES_URL,
        TAO_ARXIV_URL,
        TAO_ARXIV_URL,
        MACLANE_DOI_URL,
        TAO_NON_BRIDGING_SENTENCE,
        table,
        result["finite_packet_status"].as_str().unwrap_or("UNKNOWN"),
        result["covered_degree_count"],
        result["missing_degree_count"],
        result["finite_source_sha_fail_count"],
        result["full_ehp114_claim"],
        result["external_submission_authorized"],
        CLAIM_CEILING,
        AI_TOOLING_DISCLOSURE
    )
}

fn sha_manifest_markdown(result: &Value) -> String {
    let table = result["reviewer_grade_certificate_table"]
        .as_str()
        .unwrap_or("");
    format!(
        "# SHA manifest\n\n\
This manifest is the portable checksum summary for the finite `n < 15` \
submission preflight packet. It intentionally does not include internal \
Tao-threshold engineering artifacts.\n\n\
## Checked finite packet source\n\n\
| Role | Source | SHA status | SHA-256 |\n\
|---|---|---|---|\n\
| finite n < 15 preflight result JSON | local preflight source, not a public citation | `{}` | `{}` |\n\n\
## Public references\n\n\
- Zenodo certificate record: [{}]({})\n\
- Code and results repository: [{}]({})\n\n\
## Row-level certificate manifest\n\n\
{}\n\
## Counts\n\n\
- Covered degrees: `{}`\n\
- Missing degrees: `{}`\n\
- Source SHA failures: `{}`\n\
- Full EHP114 claim: `{}`\n\n\
## Tooling disclosure\n\n\
{}\n",
        result["finite_packet_sha_status"].as_str().unwrap_or("UNKNOWN"),
        result["finite_packet_sha256"].as_str().unwrap_or("UNKNOWN"),
        ZENODO_DOI,
        ZENODO_URL,
        GITHUB_REPO_URL,
        GITHUB_REPO_URL,
        table,
        result["covered_degree_count"],
        result["missing_degree_count"],
        result["finite_source_sha_fail_count"],
        result["full_ehp114_claim"],
        AI_TOOLING_DISCLOSURE
    )
}

fn erdosproblems_post_markdown(result: &Value) -> String {
    let table = result["reviewer_grade_certificate_table"]
        .as_str()
        .unwrap_or("");
    format!(
        "# Draft erdosproblems.com post for Problem 114\n\n\
I would like to record a finite-degree verification related to the \
Erdos-Herzog-Piranian lemniscate conjecture.\n\n\
**Finite statement.** {}\n\n\
The dependency split is deliberately finite:\n\n\
1. `n=1`: elementary, since a monic linear polynomial has a translated unit \
circle as its unit lemniscate.\n\
2. `n=2`: the MacLane / Eremenko-Hayman degree-two case. Historical citation: \
{} ([DOI]({})). Accessible source pin: Eremenko-Hayman, \"On the length of lemniscates\", \
Michigan Math. J. 46 (1999), 409-415; arXiv:0805.2295. Their abstract states \
the `d=2` extremal case, and the proof remark after Lemma 5 identifies \
`z^2+1` as extremal for `d=2`; this is length-equivalent to `z^2-1` by \
rotation.\n\
3. `3 <= n <= 14`: reproducible Rust/inari IEEE-1788 interval-certificate rows, archived under DOI \
`{}`, with SHA checks in the accompanying certificate manifest.\n\n\
## Certificate table\n\n\
{}\n\
The preflight packet records `14` covered positive degrees below `15`, `0` \
missing degree rows, and `0` source SHA failures.\n\n\
This is not meant to assert the full conjecture. Rather, it is the small-degree \
frontier complementary to Tao's sufficiently-large-`n` theorem. The remaining \
bridge question is whether Tao's threshold can be made effective at or below \
`n=15`; if not, the intervening finite degrees would still need separate \
certification. {}\n\n\
I would appreciate guidance on whether this is the right form in which to \
record the finite side on the problem page, and whether the formal-conjectures \
entry should name the theorem as a finite variant rather than changing the \
status of the main conjecture.\n\n\
Certificate record: [{}]({})\n\
Code/results repository: [{}]({})\n\
Claim ceiling: finite `n < 15` result only; conditional Tao-threshold bridge \
only.\n\n\
{}\n",
        THEOREM_STATEMENT,
        MACLANE_CITATION,
        MACLANE_DOI_URL,
        ZENODO_DOI,
        table,
        TAO_NON_BRIDGING_SENTENCE,
        ZENODO_DOI,
        ZENODO_URL,
        GITHUB_REPO_URL,
        GITHUB_REPO_URL,
        AI_TOOLING_DISCLOSURE
    )
}

fn formal_conjectures_packet_markdown(result: &Value) -> String {
    let table = result["reviewer_grade_certificate_table"]
        .as_str()
        .unwrap_or("");
    format!(
        "# Formal-conjectures packet for Erdos Problem 114\n\n\
## Intended PR title\n\n\
`Erdos114: add finite n < 15 certificate theorem`\n\n\
## Intended status change\n\n\
Do **not** mark the main theorem `erdos_114` as closed. Keep the main statement \
open unless and until the Tao threshold bridge is explicit.\n\n\
## Proposed formal-conjectures shape\n\n\
- Keep namespace: `Erdos114`.\n\
- Keep the main theorem as `@[category research open, AMS 30] theorem erdos_114 ... := by sorry`.\n\
- Add or retain the finite variant `erdos_114_finite_lt_15`.\n\
- Mark the finite variant with the repository's resolved-research category only \
if the certificate axioms are accepted as explicit external dependencies.\n\
- Keep every machine-certified row as an explicit axiom or external-certificate \
lemma; do not hide computational certification inside a proof term.\n\n\
## Suggested docstring wording\n\n\
\"Finite certified range for Erdos Problem 114. For `1 <= n < 15`, the \
Erdos-Herzog-Piranian lemniscate inequality holds, using the direct `n=1` case, \
the MacLane / Eremenko-Hayman degree-two case, and DOI-backed Rust/inari \
IEEE-1788 interval certificates for `3 <= n <= 14`. The historical `n=2` \
source is MacLane, and the accessible source pin is Alexandre Eremenko and \
Walter K. Hayman, \"On the length of lemniscates\", Michigan Math. J. 46 \
(1999), 409-415; arXiv:0805.2295. This finite statement is separate from the \
open all-degree conjecture and from Tao's sufficiently-large-`n` theorem. \
{}\"\n\n\
## PR body\n\n\
This PR records a finite variant of Erdos Problem 114 rather than changing the \
status of the main conjecture. The finite theorem covers exactly `1 <= n < 15`.\n\
The proof dependencies are typed explicitly:\n\n\
- direct analytic lemma for `n=1`;\n\
- literature-only MacLane / Eremenko-Hayman input for `n=2`, with MacLane as \
historical citation [{}]({}) and Eremenko-Hayman pinned to Michigan Math. J. 46 (1999), \
409-415 / arXiv:0805.2295;\n\
- DOI-backed Rust/inari IEEE-1788 interval certificates for `3 <= n <= 14`, \
with row-level result links and SHA-256 sidecars;\n\
- Tao's sufficiently-large-`n` theorem is mentioned only as context for the \
remaining bridge question. {}\n\n\
## Reviewer-grade certificate table\n\n\
{}\n\
The public certificate record is [{}]({}); the code/results repository is \
[{}]({}).  Local preflight checks report missing degree rows `{}`, source SHA \
failures `{}`, and `full_ehp114_claim = false`.\n\n\
## Review question\n\n\
Is this the right formal-conjectures shape: keep `erdos_114` open, add the \
finite `n < 15` variant, and leave the Tao-threshold bridge outside this PR until \
an explicit threshold is available?\n\n\
## Claim ceiling\n\n\
{}\n\n\
## Tooling disclosure\n\n\
{}\n",
        TAO_NON_BRIDGING_SENTENCE,
        MACLANE_CITATION,
        MACLANE_DOI_URL,
        TAO_NON_BRIDGING_SENTENCE,
        table,
        ZENODO_DOI,
        ZENODO_URL,
        GITHUB_REPO_URL,
        GITHUB_REPO_URL,
        result["missing_degree_count"],
        result["finite_source_sha_fail_count"],
        CLAIM_CEILING,
        AI_TOOLING_DISCLOSURE
    )
}

fn formal_conjectures_issue_draft_markdown(result: &Value) -> String {
    format!(
        "# Draft GitHub Issue: Formal-conjectures shape for Erdos 114 finite variant\n\n\
## Title\n\n\
Question: preferred shape for an Erdos 114 finite `n < 15` variant?\n\n\
## Body\n\n\
I am preparing a small, finite-range entry for Erdős Problem 114 and would like \
maintainer guidance before opening a PR.\n\n\
The intended theorem is not the all-degree conjecture. It is the finite variant:\n\n\
`For every monic complex polynomial p of degree n with 1 <= n < 15, the \
lemniscate length is at most the length for p(z)=z^n-1.`\n\n\
Proposed shape:\n\n\
- keep the main `erdos_114` theorem open;\n\
- add a named finite theorem, e.g. `erdos_114_finite_lt_15`;\n\
- keep the dependencies explicit as external axioms/certificate lemmas rather \
than hiding the computational part in proof terms;\n\
- use the direct `n=1` row, the MacLane / Eremenko-Hayman `n=2` literature row, \
and DOI-backed Rust/inari IEEE-1788 certificate rows for `3 <= n <= 14`.\n\n\
The public certificate record for `3 <= n <= 14` is [{}]({}). Row-level result \
JSON files and SHA-256 sidecars are available in [{}]({}). The `n=2` literature \
row is historically MacLane ([{}]({})), with accessible source pin Eremenko-Hayman, \
*On the length of lemniscates*, Michigan Math. J. 46 (1999), 409-415; \
arXiv:0805.2295.\n\n\
Important non-claim: {}\n\n\
{}\n\n\
Question for maintainers: would this be welcome as a finite variant alongside \
the open all-degree statement, or would you prefer the computational certificate \
dependencies to be represented differently?\n\n\
Local preflight status: missing degree rows `{}`, dependency SHA failures `{}`, \
full all-degree claim `false`.\n",
        ZENODO_DOI,
        ZENODO_URL,
        GITHUB_REPO_URL,
        GITHUB_REPO_URL,
        MACLANE_CITATION,
        MACLANE_DOI_URL,
        TAO_NON_BRIDGING_SENTENCE,
        AI_TOOLING_DISCLOSURE,
        result["missing_degree_count"],
        result["finite_source_sha_fail_count"],
    )
}

fn reviewer_inquiry_markdown(_result: &Value) -> String {
    format!(
        "# External Review Inquiry Language\n\n\
## Safe Short Version\n\n\
{}\n\n\
## What We Are Asking\n\n\
Is this the right form for recording the small side? The finite packet handles \
the small-degree frontier through `n=14`; Tao's theorem handles sufficiently \
large degrees. The remaining issue is whether the sufficiently-large threshold \
can be made explicit and whether it begins at or before `15`.\n\n\
## Source Pins Added After Third-Party Review\n\n\
- `n=2`: historical citation: {} ({}) Accessible pin: Alexandre Eremenko and Walter \
K. Hayman, \"On the length of lemniscates\", Michigan Math. J. 46 (1999), \
409-415; arXiv:0805.2295. The degree-two case is pinned to the Bernoulli \
lemniscate extremal statement.\n\
- `3 <= n <= 14`: each Rust/inari IEEE-1788 row has a public result JSON link \
and SHA-256 sidecar in the certificate table.\n\
- Tao bridge: {}\n\n\
## What We Are Not Claiming\n\n\
- We are not claiming the all-degree Erdos-Herzog-Piranian result in this packet.\n\
- We are not using the internal Tao-threshold constant-chase artifacts as \
submitted proof evidence.\n\
- We are not authorizing `n=15` or higher computation from this presentation packet.\n\n\
## Current Bridge Status\n\n\
- Bridge status: `conditional question only; no internal threshold artifacts included`\n\
- Candidate `N0`: `not asserted in this packet`\n\
- `n=15+` computation authorization: `false`\n\
- Claim ceiling: {}\n\n\
## Tooling disclosure\n\n\
{}\n",
        REVIEWER_INQUIRY_LANGUAGE,
        MACLANE_CITATION,
        MACLANE_DOI_URL,
        TAO_NON_BRIDGING_SENTENCE,
        CLAIM_CEILING,
        AI_TOOLING_DISCLOSURE
    )
}

fn report_markdown(result: &Value) -> String {
    format!(
        "# EHP114 small-n submission preflight packet\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
Status: `{}`.\n\n\
This packet freezes the portable finite theorem surface and gives two carefully \
scoped external drafts: one for erdosproblems.com and one for Google DeepMind \
formal-conjectures. The style is theorem-first, obstruction-explicit, and \
conditional about the Tao threshold bridge.\n\n\
## Presentation Title\n\n\
{}\n\n\
## Presentation Theorem\n\n\
{}\n\n\
## Dependency Audit\n\n\
- Finite packet SHA status: `{}`\n\
- Finite theorem ready: `{}`\n\
- Conditional bridge only: `{}`\n\n\
## Reviewer Hardening\n\n\
- `n=2` citation: `{}`\n\
- Row-level certificate rows: `{}`\n\
\n\
## Claim Ceiling\n\n\
{}\n\n\
## Tooling Disclosure\n\n\
{}\n",
        EXPERIMENT_ID,
        result["status"].as_str().unwrap_or("UNKNOWN"),
        THEOREM_TITLE,
        THEOREM_STATEMENT,
        result["finite_packet_sha_status"]
            .as_str()
            .unwrap_or("UNKNOWN"),
        result["finite_theorem_ready"],
        result["conditional_bridge_only"],
        result["n2_literature_citation"]
            .as_str()
            .unwrap_or("UNKNOWN"),
        result["reviewer_grade_certificate_rows"]
            .as_array()
            .map(|rows| rows.len())
            .unwrap_or(0),
        CLAIM_CEILING,
        AI_TOOLING_DISCLOSURE
    )
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let cfg = parse_args()?;
    let (
        result_path,
        report_path,
        readme_path,
        erdosproblems_path,
        formal_conjectures_path,
        formal_conjectures_issue_path,
        inquiry_path,
        sha_manifest_path,
        sha_path,
    ) = artifact_paths(&cfg.outdir);
    for path in [
        &result_path,
        &report_path,
        &readme_path,
        &erdosproblems_path,
        &formal_conjectures_path,
        &formal_conjectures_issue_path,
        &inquiry_path,
        &sha_manifest_path,
        &sha_path,
    ] {
        if path.exists() {
            return Err(format!(
                "refusing to overwrite existing artifact: {}",
                path.display()
            ));
        }
    }
    fs::create_dir_all(&cfg.outdir)
        .map_err(|err| format!("failed to create {}: {err}", cfg.outdir.display()))?;

    let finite_packet = read_json(&cfg.finite_packet)?;
    let public_rows = public_certificate_rows(&finite_packet)?;
    let reviewer_table = reviewer_grade_certificate_table(&public_rows);

    let finite_sha = verify_result_sha("finite_n_less15_packet", &cfg.finite_packet)?;
    let sha_fail_count = if finite_sha.sha_status == "PASS" {
        0
    } else {
        1
    };

    let finite_theorem_ready = str_field(&finite_packet, "status")?
        == "FINITE_N_LESS_15_PROOF_PACKET_PASS_NOT_FULL_EHP_PROOF"
        && bool_field(&finite_packet, "all_n_less_15_pass")?
        && usize_field(&finite_packet, "degree_min")? == 1
        && usize_field(&finite_packet, "degree_max_exclusive")? == 15
        && usize_field(&finite_packet, "covered_degree_count")? == 14
        && usize_field(&finite_packet, "missing_degree_count")? == 0
        && usize_field(&finite_packet, "source_sha_fail_count")? == 0
        && !bool_field(&finite_packet, "full_ehp114_claim")?;

    let bridge_is_conditional = true;

    let generated_language = [
        THEOREM_TITLE,
        THEOREM_STATEMENT,
        PUBLIC_ABSTRACT,
        FLIP_SIDE_SENTENCE,
        BRIDGE_APPENDIX_QUESTION,
        REVIEWER_INQUIRY_LANGUAGE,
        CLAIM_CEILING,
    ]
    .join("\n");
    let forbidden_language_detected = contains_forbidden_all_degree_claim(&generated_language);

    let status = if sha_fail_count == 0
        && finite_theorem_ready
        && bridge_is_conditional
        && !forbidden_language_detected
    {
        "SMALL_N_PRESENTATION_PACKET_READY_CONDITIONAL_BRIDGE_ONLY"
    } else if !finite_theorem_ready {
        "SMALL_N_PRESENTATION_PACKET_BLOCKED_FINITE_PACKET"
    } else if !bridge_is_conditional {
        "SMALL_N_PRESENTATION_PACKET_BLOCKED_BRIDGE_STATUS"
    } else if forbidden_language_detected {
        "SMALL_N_PRESENTATION_PACKET_BLOCKED_CLAIM_LANGUAGE"
    } else {
        "SMALL_N_PRESENTATION_PACKET_BLOCKED_SHA"
    };

    let first_failed_condition = if sha_fail_count > 0 {
        "dependency SHA verification failed".to_string()
    } else if !finite_theorem_ready {
        "finite n<15 packet is not ready".to_string()
    } else if !bridge_is_conditional {
        "bridge decision is not conditional or authorizes n=15+ computation".to_string()
    } else if forbidden_language_detected {
        "generated language contains a forbidden all-degree claim".to_string()
    } else {
        "none".to_string()
    };

    let result = json!({
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "presentation_title": THEOREM_TITLE,
        "theorem_statement": THEOREM_STATEMENT,
        "audiences": ["erdosproblems.com", "Google DeepMind formal-conjectures"],
        "public_reference_urls": {
            "zenodo_record": ZENODO_URL,
            "github_repository": GITHUB_REPO_URL,
            "erdosproblems_114": ERDOS114_URL,
            "formal_conjectures": FORMAL_CONJECTURES_URL,
            "tao_large_n_arxiv": TAO_ARXIV_URL,
            "maclane_doi": MACLANE_DOI_URL,
            "eremenko_hayman_arxiv": EREMENKO_HAYMAN_ARXIV_URL,
            "eremenko_hayman_pdf": EREMENKO_HAYMAN_PDF_URL
        },
        "tao_style_rules_applied": [
            "theorem-first statement",
            "typed dependency split",
            "obstruction named before bridge",
            "conditional bridge phrasing",
            "no all-degree status change"
        ],
        "public_abstract": PUBLIC_ABSTRACT,
        "flip_side_of_tao_framing": FLIP_SIDE_SENTENCE,
        "bridge_appendix_title": "Toward a Small-to-Large Degree Completion",
        "bridge_appendix_question": BRIDGE_APPENDIX_QUESTION,
        "reviewer_inquiry_language": REVIEWER_INQUIRY_LANGUAGE,
        "ai_tooling_disclosure": AI_TOOLING_DISCLOSURE,
        "finite_packet_id": FINITE_PACKET_ID,
        "finite_packet_status": finite_packet["status"].clone(),
        "finite_packet_sha_status": finite_sha.sha_status,
        "finite_packet_sha256": finite_sha.actual_sha256,
        "reviewer_grade_certificate_rows": public_rows,
        "reviewer_grade_certificate_table": reviewer_table,
        "n2_literature_citation": "Historical citation: G. R. MacLane, On a conjecture of Erdos, Herzog, and Piranian, Michigan Math. J. 2 (1953/54), 147-148, doi:10.1307/mmj/1028989918. Accessible source pin: Alexandre Eremenko and Walter K. Hayman, On the length of lemniscates, Michigan Math. J. 46 (1999), 409-415; arXiv:0805.2295.",
        "n2_literature_pin": "MacLane is cited as the earlier degree-two literature row. The accessible Eremenko-Hayman source states in the abstract that for d=2 the extremal level set is Bernoulli's lemniscate; the proof remark after Lemma 5 says z^2+1 is extremal for d=2. This is length-equivalent to z^2-1 by rotation.",
        "tao_non_bridging_sentence": TAO_NON_BRIDGING_SENTENCE,
        "degree_min": finite_packet["degree_min"].clone(),
        "degree_max_exclusive": finite_packet["degree_max_exclusive"].clone(),
        "covered_degree_count": finite_packet["covered_degree_count"].clone(),
        "missing_degree_count": finite_packet["missing_degree_count"].clone(),
        "finite_source_sha_fail_count": finite_packet["source_sha_fail_count"].clone(),
        "finite_theorem_ready": finite_theorem_ready,
        "tao_threshold_engineering_artifacts_included": false,
        "tao_threshold_status": "NOT_INCLUDED_IN_SUBMISSION_PACKET",
        "candidate_N0": null,
        "n15_or_higher_computation_authorized": false,
        "conditional_bridge_only": bridge_is_conditional,
        "submit_as_finite_result_not_full_solution": true,
        "external_submission_authorized": false,
        "bridge_work_internal_until_cleaner": true,
        "internal_threshold_artifacts_submitted_as_evidence": false,
        "full_ehp114_claim": false,
        "forbidden_language_detected": forbidden_language_detected,
        "dependency_sha_checks": [finite_sha],
        "dependency_sha_fail_count": sha_fail_count,
        "first_failed_condition": first_failed_condition,
        "claim_ceiling": CLAIM_CEILING
    });

    let result_text = serde_json::to_string_pretty(&result)
        .map_err(|err| format!("failed to serialize result JSON: {err}"))?;
    fs::write(&result_path, result_text + "\n")
        .map_err(|err| format!("failed to write {}: {err}", result_path.display()))?;
    fs::write(&report_path, report_markdown(&result))
        .map_err(|err| format!("failed to write {}: {err}", report_path.display()))?;
    fs::write(&readme_path, portable_readme_markdown(&result))
        .map_err(|err| format!("failed to write {}: {err}", readme_path.display()))?;
    fs::write(&erdosproblems_path, erdosproblems_post_markdown(&result))
        .map_err(|err| format!("failed to write {}: {err}", erdosproblems_path.display()))?;
    fs::write(
        &formal_conjectures_path,
        formal_conjectures_packet_markdown(&result),
    )
    .map_err(|err| {
        format!(
            "failed to write {}: {err}",
            formal_conjectures_path.display()
        )
    })?;
    fs::write(
        &formal_conjectures_issue_path,
        formal_conjectures_issue_draft_markdown(&result),
    )
    .map_err(|err| {
        format!(
            "failed to write {}: {err}",
            formal_conjectures_issue_path.display()
        )
    })?;
    fs::write(&inquiry_path, reviewer_inquiry_markdown(&result))
        .map_err(|err| format!("failed to write {}: {err}", inquiry_path.display()))?;
    fs::write(&sha_manifest_path, sha_manifest_markdown(&result))
        .map_err(|err| format!("failed to write {}: {err}", sha_manifest_path.display()))?;

    let digest = sha256_file(&result_path)?;
    let result_file = result_path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("RESULTS.json");
    fs::write(&sha_path, format!("{digest}  {result_file}\n"))
        .map_err(|err| format!("failed to write {}: {err}", sha_path.display()))?;

    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "experiment_id": EXPERIMENT_ID,
            "status": result["status"].clone(),
            "finite_theorem_ready": result["finite_theorem_ready"].clone(),
            "conditional_bridge_only": result["conditional_bridge_only"].clone(),
            "full_ehp114_claim": result["full_ehp114_claim"].clone(),
            "sha256": digest,
            "outdir": cfg.outdir.display().to_string()
        }))
        .unwrap()
    );

    Ok(())
}

fn main() {
    if let Err(err) = run() {
        eprintln!("error: {err}");
        std::process::exit(1);
    }
}
