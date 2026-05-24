//! EHP #114 finite proof packet for positive degrees n < 15.
//!
//! This is a packet assembler/auditor. It does not recompute the interval
//! proof. It re-reads the finite-degree index, the n=14 atlas audit, and the
//! canonical Rust/inari result files for degrees 3 through 14, then emits a
//! theorem-shaped packet for the finite range 1 <= n < 15.

use serde::Serialize;
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const EXPERIMENT_ID: &str = "EXP-MATH-EHP114-FINITE-N-LESS-15-PROOF-PACKET-20260507-01";
const L34A_ID: &str = "EXP-MATH-EHP114-N14-ATLAS-THEOREM-PACKET-AUDIT-20260507-01";
const L34B_ID: &str = "EXP-MATH-EHP114-FINITE-DEGREE-CERTIFICATE-INDEX-20260507-01";
const DEGREE_MIN: usize = 1;
const DEGREE_MAX_EXCLUSIVE: usize = 15;
const DOI_RECORD: &str = "10.5281/zenodo.19480329";
const DOI_URL: &str = "https://zenodo.org/records/19480329";
const ERDOS114_URL: &str = "https://www.erdosproblems.com/history/114";

#[derive(Clone, Debug)]
struct Config {
    outdir: PathBuf,
    finite_index: PathBuf,
    n14_atlas_audit: PathBuf,
    inari_root: PathBuf,
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
        finite_index: PathBuf::from(format!(
            "../../Erdos114/proof_path/{L34B_ID}/{L34B_ID}_RESULTS.json"
        )),
        n14_atlas_audit: PathBuf::from(format!(
            "../../Erdos114/proof_path/{L34A_ID}/{L34A_ID}_RESULTS.json"
        )),
        inari_root: PathBuf::from("../../results/erdos-114"),
    };

    let args: Vec<String> = env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--outdir" => {
                cfg.outdir = next_path(&args, i, "--outdir")?;
                i += 2;
            }
            "--finite-index" => {
                cfg.finite_index = next_path(&args, i, "--finite-index")?;
                i += 2;
            }
            "--n14-atlas-audit" => {
                cfg.n14_atlas_audit = next_path(&args, i, "--n14-atlas-audit")?;
                i += 2;
            }
            "--inari-root" => {
                cfg.inari_root = next_path(&args, i, "--inari-root")?;
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

fn str_field<'a>(v: &'a Value, key: &str) -> Result<&'a str, String> {
    v.get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("missing string field {key}"))
}

fn usize_field(v: &Value, key: &str) -> Result<usize, String> {
    let n = v
        .get(key)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("missing integer field {key}"))?;
    Ok(n as usize)
}

fn f64_field(v: &Value, key: &str) -> Result<f64, String> {
    v.get(key)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("missing numeric field {key}"))
}

fn bool_field(v: &Value, key: &str) -> Result<bool, String> {
    v.get(key)
        .and_then(Value::as_bool)
        .ok_or_else(|| format!("missing bool field {key}"))
}

fn unix_timestamp_string() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs()
        .to_string()
}

fn artifact_paths(outdir: &Path) -> (PathBuf, PathBuf, PathBuf, PathBuf) {
    (
        outdir.join(format!("{EXPERIMENT_ID}_RESULTS.json")),
        outdir.join(format!("{EXPERIMENT_ID}_REPORT.md")),
        outdir.join(format!("{EXPERIMENT_ID}_THEOREM_PACKET.md")),
        outdir.join(format!("{EXPERIMENT_ID}_RESULTS.sha256")),
    )
}

fn inari_result_path(root: &Path, degree: usize) -> PathBuf {
    root.join(format!("EXP-MM-EHP-007-n{degree}-inari_RESULTS.json"))
}

fn finite_index_row<'a>(finite_index: &'a Value, degree: usize) -> Result<&'a Value, String> {
    let rows = finite_index
        .get("rows")
        .and_then(Value::as_array)
        .ok_or_else(|| "finite index missing rows array".to_string())?;
    rows.iter()
        .find(|row| row.get("degree").and_then(Value::as_u64) == Some(degree as u64))
        .ok_or_else(|| format!("finite index missing degree row n={degree}"))
}

fn build_inari_degree_row(
    finite_index: &Value,
    inari_root: &Path,
    degree: usize,
    sha_checks: &mut Vec<ShaCheck>,
    failures: &mut Vec<Value>,
) -> Result<Value, String> {
    let index_row = finite_index_row(finite_index, degree)?;
    let result_path = inari_result_path(inari_root, degree);
    let result = read_json(&result_path)?;
    let sha = verify_result_sha(&format!("inari_degree_{degree}"), &result_path)?;
    sha_checks.push(sha.clone());

    let route_exception = if degree == 13 {
        let exception = index_row
            .get("route_exception")
            .and_then(Value::as_str)
            .ok_or_else(|| "n=13 finite index row missing route_exception".to_string())?;
        Some(exception.to_string())
    } else {
        None
    };

    let bb_level_count = result
        .get("bb_levels")
        .and_then(Value::as_array)
        .map(Vec::len)
        .unwrap_or(0);
    let bb_total_evals = usize_field(&result, "bb_total_evals")?;
    let rigor = str_field(&result, "rigor")?.to_string();
    let mut degree_failures = Vec::new();
    if sha.sha_status != "PASS" {
        degree_failures.push("sha mismatch");
    }
    if str_field(index_row, "certificate_status")? != "LOCAL_INTERVAL_CERTIFIED" {
        degree_failures.push("finite index row is not LOCAL_INTERVAL_CERTIFIED");
    }
    if str_field(index_row, "artifact_sha_status")? != "PASS" {
        degree_failures.push("finite index row artifact_sha_status is not PASS");
    }
    if str_field(&result, "experiment")? != format!("EXP-MM-EHP-007-n{degree}-inari") {
        degree_failures.push("experiment id mismatch");
    }
    if usize_field(&result, "degree")? != degree {
        degree_failures.push("degree field mismatch");
    }
    if rigor != "ieee_1788_interval_arithmetic_inari" {
        degree_failures.push("unexpected rigor field");
    }
    if !bool_field(&result, "bb_proof_complete")? {
        degree_failures.push("bb_proof_complete false");
    }
    if !bool_field(&result, "outer_domain_safe")? {
        degree_failures.push("outer_domain_safe false");
    }
    if !bool_field(&result, "hessian_negative")? {
        degree_failures.push("hessian_negative false");
    }
    if degree == 13 {
        if route_exception.is_none() {
            degree_failures.push("n=13 route_exception missing");
        }
    } else {
        if bb_total_evals == 0 {
            degree_failures.push("zero B&B evaluations outside n=13 exception");
        }
        if bb_level_count == 0 {
            degree_failures.push("empty B&B level list outside n=13 exception");
        }
    }

    if !degree_failures.is_empty() {
        failures.push(json!({
            "degree": degree,
            "failures": degree_failures
        }));
    }

    Ok(json!({
        "degree": degree,
        "coverage_status": "COVERED",
        "certificate_status": "LOCAL_INTERVAL_CERTIFIED",
        "source_type": "RUST_INARI_IEEE1788",
        "lemma_name": format!("MachineCertificateDegree{degree}"),
        "experiment_id": str_field(&result, "experiment")?,
        "verdict": str_field(&result, "verdict")?,
        "rigor": rigor,
        "artifact_path": result_path.display().to_string(),
        "artifact_sha_status": sha.sha_status,
        "artifact_sha256": sha.actual_sha256,
        "external_record_doi": DOI_RECORD,
        "external_record_url": DOI_URL,
        "reduced_dim": usize_field(&result, "reduced_dim")?,
        "bb_proof_complete": bool_field(&result, "bb_proof_complete")?,
        "bb_total_evals": bb_total_evals,
        "bb_level_count": bb_level_count,
        "outer_domain_safe": bool_field(&result, "outer_domain_safe")?,
        "hessian_negative": bool_field(&result, "hessian_negative")?,
        "l_star_lower": f64_field(&result, "l_star_lower")?,
        "l_star_upper": f64_field(&result, "l_star_upper")?,
        "l_star_interval_width": f64_field(&result, "l_star_upper")? - f64_field(&result, "l_star_lower")?,
        "route_exception": route_exception,
        "n14_atlas_cross_reference": if degree == 14 {
            json!({
                "artifact_id": L34A_ID,
                "role": "additional local stress-case atlas audit; not a replacement for the inari degree row"
            })
        } else {
            Value::Null
        },
        "claim_ceiling": "Finite-degree local interval certificate row only."
    }))
}

fn theorem_packet_markdown(result: &Value) -> String {
    format!(
        "# FiniteEHPBelow15 Theorem Packet\n\n\
## Theorem Statement\n\n\
`FiniteEHPBelow15`: For every monic complex polynomial `p` of degree `n` with \
`1 <= n < 15`, the length of the lemniscate `{{ z in C : |p(z)| = 1 }}` is at \
most the corresponding length for `p(z)=z^n-1`, under the dependency model \
listed in this packet.\n\n\
## Proof Decomposition\n\n\
1. `n=1`: direct analytic lemma. A monic linear polynomial has a translated \
unit-circle lemniscate of length `2*pi`, matching `z-1`.\n\
2. `n=2`: literature lemma. The Eremenko-Hayman row is cited as recorded by \
Erdos Problems #114.\n\
3. `3 <= n <= 14`: machine-certificate lemma. Each degree is covered by a \
Rust/inari IEEE-1788 interval certificate row from DOI `{}` with local SHA \
verification.\n\
4. Finite union lemma. The rows above exhaust every integer degree satisfying \
`1 <= n < 15`.\n\n\
## Audit Summary\n\n\
- Status: `{}`\n\
- Covered degree count: `{}`\n\
- Missing degree count: `{}`\n\
- Source SHA fail count: `{}`\n\
- n=14 atlas cross-check: `{}`\n\
- n=13 route exception retained: `{}`\n\n\
## Claim Ceiling\n\n\
This packet proves only the finite range `n < 15` under the stated literature \
and interval-certificate dependencies. It does not prove Erdos #114 for all \
degrees and does not bridge Tao's opaque high-degree threshold.\n",
        DOI_RECORD,
        result["status"].as_str().unwrap_or("UNKNOWN"),
        result["covered_degree_count"],
        result["missing_degree_count"],
        result["source_sha_fail_count"],
        result["n14_atlas_cross_check_status"]
            .as_str()
            .unwrap_or("UNKNOWN"),
        result["n13_route_exception_present"]
    )
}

fn report_markdown(result: &Value) -> String {
    format!(
        "# EHP114 Finite Proof Packet For n < 15\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
Status: `{}`.\n\n\
This packet assembles the finite proof frontier for positive degrees below 15. \
It adds the direct `n=1` analytic row, imports `n=2` as a literature row, and \
revalidates the DOI-backed Rust/inari rows for `n=3..14` by SHA.\n\n\
## Counts\n\n\
- Degree range: `{} <= n < {}`\n\
- Covered degrees: `{}`\n\
- Missing degrees: `{}`\n\
- Trivial analytic rows: `{}`\n\
- Literature-only rows: `{}`\n\
- Local interval-certified rows: `{}`\n\
- Source SHA failures: `{}`\n\n\
## Important Caveat\n\n\
The `n=13` row remains accepted as part of the DOI-backed local interval surface \
but carries a route exception because its stored branch-and-bound profile has \
zero evaluations and an empty level list.\n\n\
## Claim Ceiling\n\n\
{}\n",
        EXPERIMENT_ID,
        result["status"].as_str().unwrap_or("UNKNOWN"),
        result["degree_min"],
        result["degree_max_exclusive"],
        result["covered_degree_count"],
        result["missing_degree_count"],
        result["trivial_analytic_degree_count"],
        result["literature_only_degree_count"],
        result["local_interval_certified_degree_count"],
        result["source_sha_fail_count"],
        result["claim_ceiling"]
            .as_str()
            .unwrap_or("finite packet only")
    )
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let cfg = parse_args()?;
    let (result_path, report_path, theorem_path, sha_path) = artifact_paths(&cfg.outdir);
    for path in [&result_path, &report_path, &theorem_path, &sha_path] {
        if path.exists() {
            return Err(format!(
                "refusing to overwrite existing artifact: {}",
                path.display()
            ));
        }
    }
    fs::create_dir_all(&cfg.outdir)
        .map_err(|err| format!("failed to create {}: {err}", cfg.outdir.display()))?;

    let finite_index = read_json(&cfg.finite_index)?;
    let n14_atlas = read_json(&cfg.n14_atlas_audit)?;
    let finite_index_sha = verify_result_sha("finite_degree_certificate_index", &cfg.finite_index)?;
    let n14_atlas_sha = verify_result_sha("n14_atlas_theorem_packet_audit", &cfg.n14_atlas_audit)?;

    let mut rows = Vec::new();
    let mut source_sha_checks = vec![finite_index_sha.clone(), n14_atlas_sha.clone()];
    let mut row_failures = Vec::new();

    rows.push(json!({
        "degree": 1,
        "coverage_status": "COVERED",
        "certificate_status": "TRIVIAL_ANALYTIC",
        "source_type": "DIRECT_ANALYTIC",
        "lemma_name": "LinearMonicLemniscateCircle",
        "statement": "For p(z)=z+a monic linear, {|p(z)|=1} is a translated unit circle of length 2*pi, matching z-1.",
        "artifact_sha_status": "NOT_APPLICABLE",
        "claim_ceiling": "Direct analytic finite-degree row only."
    }));

    let n2_index_row = finite_index_row(&finite_index, 2)?;
    if str_field(n2_index_row, "certificate_status")? != "LITERATURE_ONLY" {
        row_failures.push(json!({
            "degree": 2,
            "failures": ["finite index n=2 row is not LITERATURE_ONLY"]
        }));
    }
    rows.push(json!({
        "degree": 2,
        "coverage_status": "COVERED",
        "certificate_status": "LITERATURE_ONLY",
        "source_type": "LITERATURE",
        "lemma_name": "EremenkoHaymanDegreeTwo",
        "literature_reference": "Eremenko-Hayman n=2 row as recorded by Erdos Problems #114.",
        "literature_status_url": ERDOS114_URL,
        "artifact_sha_status": "LITERATURE_ONLY",
        "claim_ceiling": "Literature finite-degree row only."
    }));

    for degree in 3..=14 {
        let row = build_inari_degree_row(
            &finite_index,
            &cfg.inari_root,
            degree,
            &mut source_sha_checks,
            &mut row_failures,
        )?;
        rows.push(row);
    }

    let covered_degrees: BTreeSet<usize> = rows
        .iter()
        .filter_map(|row| {
            row.get("degree")
                .and_then(Value::as_u64)
                .map(|n| n as usize)
        })
        .collect();
    let expected_degrees: BTreeSet<usize> = (DEGREE_MIN..DEGREE_MAX_EXCLUSIVE).collect();
    let missing_degrees: Vec<usize> = expected_degrees
        .difference(&covered_degrees)
        .copied()
        .collect();
    let extra_degrees: Vec<usize> = covered_degrees
        .difference(&expected_degrees)
        .copied()
        .collect();

    let source_sha_fail_count = source_sha_checks
        .iter()
        .filter(|check| check.sha_status != "PASS")
        .count();
    let trivial_analytic_count = rows
        .iter()
        .filter(|row| row["certificate_status"] == "TRIVIAL_ANALYTIC")
        .count();
    let literature_only_count = rows
        .iter()
        .filter(|row| row["certificate_status"] == "LITERATURE_ONLY")
        .count();
    let local_interval_count = rows
        .iter()
        .filter(|row| row["certificate_status"] == "LOCAL_INTERVAL_CERTIFIED")
        .count();
    let n13_route_exception_present = rows.iter().any(|row| {
        row.get("degree").and_then(Value::as_u64) == Some(13)
            && !row.get("route_exception").unwrap_or(&Value::Null).is_null()
    });

    let finite_index_status_pass =
        str_field(&finite_index, "status")? == "FINITE_DEGREE_INDEX_PASS_NOT_FULL_PROOF";
    let n14_atlas_status_pass =
        str_field(&n14_atlas, "status")? == "N14_ATLAS_THEOREM_PACKET_AUDIT_PASS_NOT_FULL_PROOF";
    let n14_atlas_cross_check_status = if n14_atlas_status_pass
        && usize_field(&n14_atlas, "unique_certified_cell_count")? == 64
        && usize_field(&n14_atlas, "selected_child_packet_sha_fail_count")? == 0
    {
        "PASS"
    } else {
        "FAIL"
    };

    if !finite_index_status_pass {
        row_failures.push(json!({
            "dependency": L34B_ID,
            "failures": ["finite index status did not pass"]
        }));
    }
    if !n14_atlas_status_pass {
        row_failures.push(json!({
            "dependency": L34A_ID,
            "failures": ["n=14 atlas audit status did not pass"]
        }));
    }
    if !n13_route_exception_present {
        row_failures.push(json!({
            "degree": 13,
            "failures": ["required n=13 route_exception is absent"]
        }));
    }

    let all_n_less_15_pass = missing_degrees.is_empty()
        && extra_degrees.is_empty()
        && source_sha_fail_count == 0
        && row_failures.is_empty()
        && trivial_analytic_count == 1
        && literature_only_count == 1
        && local_interval_count == 12
        && n14_atlas_cross_check_status == "PASS";
    let status = if all_n_less_15_pass {
        "FINITE_N_LESS_15_PROOF_PACKET_PASS_NOT_FULL_EHP_PROOF"
    } else {
        "FINITE_N_LESS_15_PROOF_PACKET_FAIL"
    };
    let first_failed_condition = if source_sha_fail_count > 0 {
        "at least one source SHA check failed".to_string()
    } else if !missing_degrees.is_empty() {
        format!("missing degree rows: {missing_degrees:?}")
    } else if !extra_degrees.is_empty() {
        format!("unexpected extra degree rows: {extra_degrees:?}")
    } else if !row_failures.is_empty() {
        "at least one row failed packet validation".to_string()
    } else if n14_atlas_cross_check_status != "PASS" {
        "n=14 atlas cross-check failed".to_string()
    } else {
        "none".to_string()
    };

    let result = json!({
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "theorem_name": "FiniteEHPBelow15",
        "theorem_statement": "For every monic complex polynomial p of degree n with 1 <= n < 15, the lemniscate length is at most the corresponding length for p(z)=z^n-1, under the cited literature and interval-certificate dependencies.",
        "degree_min": DEGREE_MIN,
        "degree_max_exclusive": DEGREE_MAX_EXCLUSIVE,
        "covered_degree_count": covered_degrees.len(),
        "covered_degrees": covered_degrees.iter().copied().collect::<Vec<_>>(),
        "missing_degree_count": missing_degrees.len(),
        "missing_degrees": missing_degrees,
        "extra_degrees": extra_degrees,
        "trivial_analytic_degree_count": trivial_analytic_count,
        "literature_only_degree_count": literature_only_count,
        "local_interval_certified_degree_count": local_interval_count,
        "degree_rows": rows,
        "dependency_checks": [
            {
                "artifact_id": L34B_ID,
                "path": cfg.finite_index.display().to_string(),
                "status": finite_index["status"].clone(),
                "sha_status": finite_index_sha.sha_status
            },
            {
                "artifact_id": L34A_ID,
                "path": cfg.n14_atlas_audit.display().to_string(),
                "status": n14_atlas["status"].clone(),
                "sha_status": n14_atlas_sha.sha_status
            }
        ],
        "source_sha_checks": source_sha_checks,
        "source_sha_fail_count": source_sha_fail_count,
        "row_failure_count": row_failures.len(),
        "row_failures": row_failures,
        "n13_route_exception_present": n13_route_exception_present,
        "n14_atlas_cross_check_status": n14_atlas_cross_check_status,
        "n14_atlas_tightest_margin_cell": n14_atlas["tightest_margin_cell"].clone(),
        "external_record_doi": DOI_RECORD,
        "external_record_url": DOI_URL,
        "literature_status_url": ERDOS114_URL,
        "all_n_less_15_pass": all_n_less_15_pass,
        "full_ehp114_claim": false,
        "tao_threshold_bridge_status": "NOT_BRIDGED_IN_THIS_PACKET",
        "first_failed_condition": first_failed_condition,
        "claim_ceiling": "This packet proves only the finite range n < 15 under the stated literature and interval-certificate dependencies. It does not prove Erdos #114 for all n and does not bridge Tao's opaque high-degree threshold."
    });

    let report = report_markdown(&result);
    let theorem_packet = theorem_packet_markdown(&result);
    let result_text = serde_json::to_string_pretty(&result)
        .map_err(|err| format!("failed to serialize result JSON: {err}"))?;
    fs::write(&result_path, result_text + "\n")
        .map_err(|err| format!("failed to write {}: {err}", result_path.display()))?;
    fs::write(&report_path, report)
        .map_err(|err| format!("failed to write {}: {err}", report_path.display()))?;
    fs::write(&theorem_path, theorem_packet)
        .map_err(|err| format!("failed to write {}: {err}", theorem_path.display()))?;
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
            "all_n_less_15_pass": result["all_n_less_15_pass"].clone(),
            "covered_degree_count": result["covered_degree_count"].clone(),
            "missing_degree_count": result["missing_degree_count"].clone(),
            "source_sha_fail_count": result["source_sha_fail_count"].clone(),
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
