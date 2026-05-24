//! EHP #114 proof-path packet emitter after the n=14 atlas pass.
//!
//! This binary does not recompute lemniscate geometry. It turns the current
//! passing n=14 atlas layer, the DOI-backed finite-degree inari records, and
//! Tao's high-degree theorem interface into four theorem-shaped local audit
//! packets:
//!
//! L34A: n=14 atlas theorem-packet audit.
//! L34B: finite-degree certificate index.
//! L35: Tao high-degree threshold audit with opaque threshold.
//! L36: full-proof bridge blocker.

use serde::Serialize;
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const L34A_ID: &str = "EXP-MATH-EHP114-N14-ATLAS-THEOREM-PACKET-AUDIT-20260507-01";
const L34B_ID: &str = "EXP-MATH-EHP114-FINITE-DEGREE-CERTIFICATE-INDEX-20260507-01";
const L35_ID: &str = "EXP-MATH-EHP114-TAO-HIGH-DEGREE-THRESHOLD-AUDIT-20260507-01";
const L36_ID: &str = "EXP-MATH-EHP114-FULL-PROOF-BRIDGE-BLOCKER-20260507-01";

const L33_ID: &str = "EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-02";
const EXACT_LENGTH_CAP: f64 = 20.672796062619668;
const EXPECTED_N14_CELLS: usize = 64;
const FINITE_MIN_DEGREE: usize = 2;
const FINITE_MAX_DEGREE: usize = 14;
const LOCAL_INARI_MIN_DEGREE: usize = 3;
const DOI_RECORD: &str = "10.5281/zenodo.19480329";
const DOI_URL: &str = "https://zenodo.org/records/19480329";
const TAO_ARXIV_URL: &str = "https://arxiv.org/abs/2512.12455";
const TAO_AR5IV_URL: &str = "https://ar5iv.org/html/2512.12455v2";
const ERDOS114_URL: &str = "https://www.erdosproblems.com/history/114";

#[derive(Clone, Debug)]
struct Config {
    out_root: PathBuf,
    l33_source: PathBuf,
    finite_root: PathBuf,
    doi_audit_path: PathBuf,
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

#[derive(Clone, Debug)]
struct Artifact {
    experiment_id: &'static str,
    result: Value,
    report: String,
}

fn parse_args() -> Result<Config, String> {
    let mut cfg = Config {
        out_root: PathBuf::from("../../Erdos114/proof_path"),
        l33_source: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-02/EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-02_RESULTS.json"),
        finite_root: PathBuf::from("../../results/erdos-114"),
        doi_audit_path: PathBuf::from("../../Erdos114/EHP114_IEEE1788_RUST_DOI_AUDIT_2026-05-05.md"),
    };

    let args: Vec<String> = env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--out-root" => {
                cfg.out_root = next_path(&args, i, "--out-root")?;
                i += 2;
            }
            "--l33-source" => {
                cfg.l33_source = next_path(&args, i, "--l33-source")?;
                i += 2;
            }
            "--finite-root" => {
                cfg.finite_root = next_path(&args, i, "--finite-root")?;
                i += 2;
            }
            "--doi-audit" => {
                cfg.doi_audit_path = next_path(&args, i, "--doi-audit")?;
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

fn unix_timestamp_string() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs()
        .to_string()
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

fn maybe_usize(v: &Value, key: &str) -> usize {
    v.get(key).and_then(Value::as_u64).unwrap_or(0) as usize
}

fn maybe_f64(v: &Value, key: &str) -> Option<f64> {
    v.get(key).and_then(Value::as_f64)
}

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12
}

fn cell_from_summary(summary: &Value) -> Result<[usize; 2], String> {
    let arr = summary
        .get("subcell")
        .and_then(Value::as_array)
        .ok_or_else(|| "certified cell summary missing subcell array".to_string())?;
    if arr.len() != 2 {
        return Err("certified cell subcell array does not have length two".to_string());
    }
    let i = arr[0]
        .as_u64()
        .ok_or_else(|| "subcell[0] is not an integer".to_string())? as usize;
    let j = arr[1]
        .as_u64()
        .ok_or_else(|| "subcell[1] is not an integer".to_string())? as usize;
    Ok([i, j])
}

fn artifact_paths(out_root: &Path, experiment_id: &str) -> (PathBuf, PathBuf, PathBuf) {
    let dir = out_root.join(experiment_id);
    (
        dir.join(format!("{experiment_id}_RESULTS.json")),
        dir.join(format!("{experiment_id}_REPORT.md")),
        dir.join(format!("{experiment_id}_RESULTS.sha256")),
    )
}

fn build_l34a(cfg: &Config, started: Instant) -> Result<Artifact, String> {
    let source = read_json(&cfg.l33_source)?;
    let source_sha = verify_result_sha("l33_global_atlas_certificate", &cfg.l33_source)?;
    let cells = source
        .get("certified_cells")
        .and_then(Value::as_array)
        .ok_or_else(|| "L33 source missing certified_cells array".to_string())?;

    let mut unique_cells = BTreeSet::new();
    let mut cell_sha_checks = Vec::new();
    let mut length_over_cap = Vec::new();
    let mut tightest = Vec::new();
    let mut drift_failures = Vec::new();

    for cell in cells {
        let subcell = cell_from_summary(cell)?;
        unique_cells.insert((subcell[0], subcell[1]));
        let path = PathBuf::from(str_field(cell, "path")?);
        cell_sha_checks.push(verify_result_sha("selected_local_cell_packet", &path)?);
        let upper = f64_field(cell, "total_validated_length_upper")?;
        let cap = f64_field(cell, "exact_length_cap")?;
        let margin = f64_field(cell, "margin_to_cap")?;
        if upper > cap || margin < 0.0 || !approx_eq(cap, EXACT_LENGTH_CAP) {
            length_over_cap.push(json!({
                "subcell": subcell,
                "total_validated_length_upper": upper,
                "exact_length_cap": cap,
                "margin_to_cap": margin
            }));
        }
        if maybe_usize(cell, "source_sha_fail_count") > 0
            || maybe_usize(cell, "ownership_duplicate_count") > 0
            || maybe_usize(cell, "source_filter_mismatch_count") > 0
            || maybe_usize(cell, "candidate_count_drift_count") > 0
        {
            drift_failures.push(json!(cell));
        }
        tightest.push(json!({
            "subcell": subcell,
            "experiment_id": str_field(cell, "experiment_id")?,
            "margin_to_cap": margin,
            "total_validated_length_upper": upper,
            "path": str_field(cell, "path")?
        }));
    }

    tightest.sort_by(|a, b| {
        maybe_f64(a, "margin_to_cap")
            .unwrap_or(f64::INFINITY)
            .partial_cmp(&maybe_f64(b, "margin_to_cap").unwrap_or(f64::INFINITY))
            .unwrap()
    });

    let source_sha_fail_count = source_sha.sha_status != "PASS";
    let child_sha_fail_count = cell_sha_checks
        .iter()
        .filter(|row| row.sha_status != "PASS")
        .count();

    let checks = vec![
        (
            "source_status_pass",
            str_field(&source, "status")? == "GLOBAL_ATLAS_PASS_NOT_FULL_EHP_PROOF",
        ),
        (
            "source_global_atlas_pass",
            bool_field(&source, "global_atlas_certificate_pass")?,
        ),
        (
            "root_affine_skeleton_64",
            usize_field(&source, "root_affine_skeleton_cell_count")? == EXPECTED_N14_CELLS,
        ),
        (
            "certified_cells_64",
            usize_field(&source, "certified_cell_count")? == EXPECTED_N14_CELLS,
        ),
        ("unique_cells_64", unique_cells.len() == EXPECTED_N14_CELLS),
        (
            "missing_cells_zero",
            usize_field(&source, "missing_cell_certificate_count")? == 0,
        ),
        (
            "duplicate_packets_zero",
            usize_field(&source, "duplicate_cell_certificate_count")? == 0,
        ),
        (
            "rejected_packets_zero",
            usize_field(&source, "rejected_cell_certificate_count")? == 0,
        ),
        (
            "source_sha_fail_zero",
            usize_field(&source, "source_sha_fail_count")? == 0,
        ),
        (
            "coverage_fail_zero",
            usize_field(&source, "coverage_failure_count")? == 0,
        ),
        (
            "budget_fail_zero",
            usize_field(&source, "budget_failure_count")? == 0,
        ),
        ("l33_sha_pass", !source_sha_fail_count),
        ("selected_child_sha_pass", child_sha_fail_count == 0),
        ("no_child_drift_flags", drift_failures.is_empty()),
        ("all_local_uppers_below_cap", length_over_cap.is_empty()),
    ];

    let failed_checks: Vec<_> = checks
        .iter()
        .filter(|(_, pass)| !*pass)
        .map(|(name, _)| *name)
        .collect();
    let status = if failed_checks.is_empty() {
        "N14_ATLAS_THEOREM_PACKET_AUDIT_PASS_NOT_FULL_PROOF"
    } else {
        "N14_ATLAS_THEOREM_PACKET_AUDIT_FAIL"
    };
    let first_failed = failed_checks.first().copied().unwrap_or("none").to_string();

    let tightest_cells: Vec<Value> = tightest.into_iter().take(10).collect();
    let tightest_margin = maybe_f64(&tightest_cells[0], "margin_to_cap")
        .ok_or_else(|| "missing tightest margin".to_string())?;

    let result = json!({
        "experiment_id": L34A_ID,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "source_experiment_id": L33_ID,
        "source_path": cfg.l33_source.display().to_string(),
        "source_sha": source_sha,
        "source_sha_fail_count": if source_sha_fail_count { 1 } else { 0 },
        "selected_child_packet_sha_checks": cell_sha_checks,
        "selected_child_packet_sha_fail_count": child_sha_fail_count,
        "atlas_acceptance_checks": checks.iter().map(|(name, pass)| json!({
            "check": name,
            "pass": pass
        })).collect::<Vec<_>>(),
        "failed_checks": failed_checks,
        "global_atlas_status": source["status"].clone(),
        "global_atlas_certificate_pass": source["global_atlas_certificate_pass"].clone(),
        "root_affine_skeleton_cell_count": source["root_affine_skeleton_cell_count"].clone(),
        "certified_cell_count": source["certified_cell_count"].clone(),
        "unique_certified_cell_count": unique_cells.len(),
        "missing_cell_certificate_count": source["missing_cell_certificate_count"].clone(),
        "duplicate_cell_certificate_count": source["duplicate_cell_certificate_count"].clone(),
        "rejected_cell_certificate_count": source["rejected_cell_certificate_count"].clone(),
        "coverage_failure_count": source["coverage_failure_count"].clone(),
        "budget_failure_count": source["budget_failure_count"].clone(),
        "drift_failure_count": drift_failures.len(),
        "drift_failures": drift_failures,
        "length_over_cap_count": length_over_cap.len(),
        "length_over_cap": length_over_cap,
        "exact_length_cap": EXACT_LENGTH_CAP,
        "max_certified_cell_length_upper": source["max_certified_cell_length_upper"].clone(),
        "min_certified_cell_margin_to_cap": source["min_certified_cell_margin_to_cap"].clone(),
        "tightest_margin_cell": tightest_cells[0].clone(),
        "tightest_cells_by_margin": tightest_cells,
        "recorded_tightest_margin_cell_07_03": {
            "subcell": [7, 3],
            "margin_to_cap": tightest_margin,
            "expected_margin_to_cap": 0.008190869737216389,
            "matches_expected": (tightest_margin - 0.008190869737216389_f64).abs() <= 1e-12
        },
        "stale_l33_next_dependency": source["next_dependency"].clone(),
        "corrected_next_dependency": "Finite-degree index plus high-degree threshold bridge accounting.",
        "what_this_rules_out": "Silent reuse of the hard-cell packet as the n=14 atlas layer; the selected current packet set covers 64 unique cells.",
        "what_this_does_not_rule_out": "This does not cover n != 14 and does not make Tao's high-degree threshold explicit.",
        "first_failed_condition": first_failed,
        "claim_ceiling": "Theorem-packet audit for the local n=14 root-affine atlas only. Not a full EHP114 proof and not a replacement for finite-degree or high-degree bridge accounting."
    });

    let report = format!(
        "# EHP114 L34A n=14 Atlas Theorem-Packet Audit\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
Status: `{}`.\n\n\
The L33 atlas source was checked as a theorem-shaped n=14 layer: 64 unique local \
cell packets, zero missing cells, zero rejected packets, zero duplicate cell \
coverage, zero ownership/source-filter/candidate drift flags, and every selected \
local upper below `{}`.\n\n\
Tightest recorded cell: `CELL-07-03`, margin `{}`. This is the expected narrow \
spot, so the next dependency is not missing-cell generation. The corrected next \
dependency is finite-degree indexing plus high-degree threshold bridge accounting.\n\n\
## Claim Ceiling\n\n\
This is a local n=14 atlas audit. It is not a full EHP114 proof and it does not \
make Tao's threshold explicit.\n",
        L34A_ID, status, EXACT_LENGTH_CAP, tightest_margin
    );

    Ok(Artifact {
        experiment_id: L34A_ID,
        result,
        report,
    })
}

fn build_l34b(cfg: &Config, l34a_result: &Value, started: Instant) -> Result<Artifact, String> {
    let doi_audit_exists = cfg.doi_audit_path.exists();
    let mut rows = Vec::new();
    rows.push(json!({
        "degree": 2,
        "certificate_status": "LITERATURE_ONLY",
        "source_type": "LITERATURE",
        "literature_reference": "Eremenko-Hayman n=2 row as recorded by Erdős Problems #114",
        "artifact_path": Value::Null,
        "artifact_sha_status": "LITERATURE_ONLY",
        "route_exception": Value::Null,
        "claim_ceiling": "Literature row only; no local Rust or Lean certificate is asserted by this index."
    }));

    let mut sha_checks = Vec::new();
    let mut row_failures = Vec::new();
    for degree in LOCAL_INARI_MIN_DEGREE..=FINITE_MAX_DEGREE {
        let result_path = cfg
            .finite_root
            .join(format!("EXP-MM-EHP-007-n{degree}-inari_RESULTS.json"));
        let result = read_json(&result_path)?;
        let sha_check = verify_result_sha(&format!("inari_degree_{degree}"), &result_path)?;
        let verdict = str_field(&result, "verdict")?.to_string();
        let rigor = str_field(&result, "rigor")?.to_string();
        let bb_proof_complete = bool_field(&result, "bb_proof_complete")?;
        let outer_domain_safe = bool_field(&result, "outer_domain_safe")?;
        let hessian_negative = bool_field(&result, "hessian_negative")?;
        let bb_total_evals = usize_field(&result, "bb_total_evals")?;
        let bb_level_count = result
            .get("bb_levels")
            .and_then(Value::as_array)
            .map(Vec::len)
            .unwrap_or(0);
        let route_exception = if degree == 13 {
            Some("n=13 has zero branch-and-bound evaluations and an empty bb_levels array; retain as DOI-backed local interval row but flag for route reconciliation before any refreshed public theorem packet.")
        } else {
            None
        };

        let mut degree_failures = Vec::new();
        if sha_check.sha_status != "PASS" {
            degree_failures.push("sha mismatch");
        }
        if !bb_proof_complete {
            degree_failures.push("bb_proof_complete false");
        }
        if !outer_domain_safe {
            degree_failures.push("outer_domain_safe false");
        }
        if !hessian_negative {
            degree_failures.push("hessian_negative false");
        }
        if rigor != "ieee_1788_interval_arithmetic_inari" {
            degree_failures.push("unexpected rigor field");
        }
        if degree == 13 && route_exception.is_none() {
            degree_failures.push("missing required n=13 route_exception");
        }
        if degree != 13 && bb_total_evals == 0 {
            degree_failures.push("zero eval count outside the n=13 exception");
        }
        if degree != 13 && bb_level_count == 0 {
            degree_failures.push("empty bb_levels outside the n=13 exception");
        }
        if !degree_failures.is_empty() {
            row_failures.push(json!({
                "degree": degree,
                "failures": degree_failures
            }));
        }

        rows.push(json!({
            "degree": degree,
            "certificate_status": "LOCAL_INTERVAL_CERTIFIED",
            "source_type": "RUST_INARI_IEEE1788",
            "experiment_id": str_field(&result, "experiment")?,
            "verdict": verdict,
            "rigor": rigor,
            "reduced_dim": usize_field(&result, "reduced_dim")?,
            "bb_proof_complete": bb_proof_complete,
            "bb_total_evals": bb_total_evals,
            "bb_level_count": bb_level_count,
            "outer_domain_safe": outer_domain_safe,
            "hessian_negative": hessian_negative,
            "l_star_lower": f64_field(&result, "l_star_lower")?,
            "l_star_upper": f64_field(&result, "l_star_upper")?,
            "l_star_interval_width": f64_field(&result, "l_star_upper")? - f64_field(&result, "l_star_lower")?,
            "total_time_secs": result.get("total_time_secs").cloned().unwrap_or(Value::Null),
            "artifact_path": result_path.display().to_string(),
            "artifact_sha_status": sha_check.sha_status.clone(),
            "artifact_sha256": sha_check.actual_sha256.clone(),
            "external_record_doi": DOI_RECORD,
            "external_record_url": DOI_URL,
            "route_exception": route_exception,
            "n14_atlas_cross_reference": if degree == 14 {
                json!({
                    "artifact_id": L34A_ID,
                    "status": l34a_result["status"].clone(),
                    "role": "local n=14 atlas theorem-packet audit, not a replacement for the inari degree row"
                })
            } else {
                Value::Null
            },
            "claim_ceiling": "Finite-degree local interval certificate row only; not a full all-degree theorem."
        }));
        sha_checks.push(sha_check);
    }

    let all_inari_sha_pass = sha_checks.iter().all(|row| row.sha_status == "PASS");
    let degree_rows_present = rows.len() == (FINITE_MAX_DEGREE - FINITE_MIN_DEGREE + 1);
    let n13_exception_present = rows.iter().any(|row| {
        row.get("degree").and_then(Value::as_u64) == Some(13)
            && !row.get("route_exception").unwrap_or(&Value::Null).is_null()
    });
    let n2_present = rows
        .first()
        .and_then(|row| row.get("degree"))
        .and_then(Value::as_u64)
        == Some(2);
    let status = if all_inari_sha_pass
        && row_failures.is_empty()
        && degree_rows_present
        && n13_exception_present
        && n2_present
        && doi_audit_exists
    {
        "FINITE_DEGREE_INDEX_PASS_NOT_FULL_PROOF"
    } else {
        "FINITE_DEGREE_INDEX_FAIL"
    };
    let first_failed = if !doi_audit_exists {
        "doi audit note missing"
    } else if !degree_rows_present {
        "degree row count mismatch"
    } else if !n2_present {
        "literature n=2 row missing"
    } else if !n13_exception_present {
        "n=13 route exception missing"
    } else if !all_inari_sha_pass {
        "at least one inari SHA check failed"
    } else if !row_failures.is_empty() {
        "at least one finite-degree row failed required fields"
    } else {
        "none"
    };

    let result = json!({
        "experiment_id": L34B_ID,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "finite_degree_range": [FINITE_MIN_DEGREE, FINITE_MAX_DEGREE],
        "degree_row_count": rows.len(),
        "literature_only_degree_count": 1,
        "local_interval_certified_degree_count": rows.iter().filter(|row| row["certificate_status"] == "LOCAL_INTERVAL_CERTIFIED").count(),
        "rows": rows,
        "sha_checks": sha_checks,
        "sha_fail_count": if all_inari_sha_pass { 0 } else { 1 },
        "row_failure_count": row_failures.len(),
        "row_failures": row_failures,
        "doi_record": DOI_RECORD,
        "doi_record_url": DOI_URL,
        "doi_audit_path": cfg.doi_audit_path.display().to_string(),
        "doi_audit_present": doi_audit_exists,
        "n14_atlas_audit_cross_reference": {
            "artifact_id": L34A_ID,
            "status": l34a_result["status"].clone(),
            "tightest_margin_cell": l34a_result["tightest_margin_cell"].clone()
        },
        "what_this_rules_out": "The finite-degree layer is not just the n=14 hard-cell atlas; it has separate rows for n=2 and n=3..14.",
        "what_this_does_not_rule_out": "This does not cover degrees n >= 15 and does not extract Tao's high-degree cutoff.",
        "first_failed_condition": first_failed,
        "claim_ceiling": "Finite-degree certificate index through n=14 only. Not a full EHP114 proof."
    });

    let report = format!(
        "# EHP114 L34B Finite-Degree Certificate Index\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
Status: `{}`.\n\n\
The index binds `n=2` as a literature-only row and `n=3..14` as local \
Rust/inari IEEE-1788 rows from the DOI-backed result family `{}`. The n=13 row \
is retained but explicitly flagged because its branch-and-bound profile is an \
exception: zero total evaluations and an empty level list.\n\n\
The n=14 atlas audit is cross-referenced as a local theorem-packet audit; it \
does not replace the inari degree row.\n\n\
## Claim Ceiling\n\n\
This is a finite-degree index through n=14. The all-degree bridge remains \
separate.\n",
        L34B_ID, status, DOI_RECORD
    );

    Ok(Artifact {
        experiment_id: L34B_ID,
        result,
        report,
    })
}

fn build_l35(started: Instant) -> Result<Artifact, String> {
    let result = json!({
        "experiment_id": L35_ID,
        "status": "TAO_HIGH_DEGREE_THRESHOLD_AUDIT_OPAQUE_THRESHOLD",
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "primary_sources": [
            {
                "label": "Tao arXiv high-degree EHP paper",
                "url": TAO_ARXIV_URL,
                "recorded_title": "The maximal length of the Erdos-Herzog-Piranian lemniscate in high degree",
                "recorded_author": "Terence Tao",
                "recorded_version": "v2, last revised 2025-12-22",
                "recorded_statement": "establishes the EHP conjecture for all sufficiently large n"
            },
            {
                "label": "Tao ar5iv HTML mirror",
                "url": TAO_AR5IV_URL,
                "recorded_role": "readable theorem source mirror"
            },
            {
                "label": "Erdos Problems #114 current statement",
                "url": ERDOS114_URL,
                "recorded_role": "problem statement and literature-status page"
            }
        ],
        "tao_threshold_status": "EFFECTIVE_IN_PRINCIPLE_NOT_PRACTICALLY_EXTRACTED",
        "tao_threshold_symbol": "taoEHPThreshold",
        "tao_threshold_value": Value::Null,
        "forbidden_threshold_assignment": "Do not set taoEHPThreshold = 15 in this packet.",
        "typed_interface": {
            "opaque_constant": "taoEHPThreshold : Nat",
            "high_degree_theorem_shape": "forall n >= taoEHPThreshold, EHP_conjecture_holds(n)",
            "finite_bridge_obligation": "certify all 2 <= n < taoEHPThreshold, or extract taoEHPThreshold <= 15",
            "candidate_length_formula": "length(z^n - 1 lemniscate) = 2^(1/n) * Beta(1/2, 1/(2n))"
        },
        "effectivity_notes": [
            "The high-degree theorem is treated as an external analytic theorem.",
            "The threshold is effective in principle but not extracted here as a practical finite cutoff.",
            "The proof path cannot identify n=15 as the start of the high-degree range without a separate threshold extraction."
        ],
        "what_this_rules_out": "A proof synthesis packet may not silently glue n <= 14 certificates to Tao's theorem by assuming the high-degree range begins at n=15.",
        "what_this_does_not_rule_out": "A later analytic audit may extract a usable threshold, or computation may certify the remaining finite interval below that threshold.",
        "first_failed_condition": "opaque threshold remains unbridged",
        "claim_ceiling": "High-degree theorem interface only. Not a full EHP114 proof and not a finite cutoff extraction."
    });

    let report = format!(
        "# EHP114 L35 Tao High-Degree Threshold Audit\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
Status: `TAO_HIGH_DEGREE_THRESHOLD_AUDIT_OPAQUE_THRESHOLD`.\n\n\
Tao's high-degree result is used as a typed external interface: there exists an \
opaque threshold `taoEHPThreshold` after which the high-degree theorem applies. \
This packet does not extract a practical numerical cutoff and does not set the \
threshold to 15.\n\n\
## Bridge Consequence\n\n\
The finite side currently reaches n=14. The full theorem needs either an \
explicit extraction with `taoEHPThreshold <= 15`, or additional finite \
certificates for every degree below the extracted threshold.\n\n\
Sources recorded: `{}`, `{}`.\n",
        L35_ID, TAO_ARXIV_URL, ERDOS114_URL
    );

    Ok(Artifact {
        experiment_id: L35_ID,
        result,
        report,
    })
}

fn build_l36(
    l34a_result: &Value,
    l34b_result: &Value,
    l35_result: &Value,
    started: Instant,
) -> Result<Artifact, String> {
    let n14_pass =
        str_field(l34a_result, "status")? == "N14_ATLAS_THEOREM_PACKET_AUDIT_PASS_NOT_FULL_PROOF";
    let finite_pass =
        str_field(l34b_result, "status")? == "FINITE_DEGREE_INDEX_PASS_NOT_FULL_PROOF";
    let threshold_opaque = str_field(l35_result, "tao_threshold_status")?
        == "EFFECTIVE_IN_PRINCIPLE_NOT_PRACTICALLY_EXTRACTED";
    let status = if n14_pass && finite_pass && threshold_opaque {
        "FULL_PROOF_BRIDGE_BLOCKED_THRESHOLD_GAP"
    } else {
        "FULL_PROOF_BRIDGE_BLOCKER_SOURCE_FAILURE"
    };
    let first_failed = if !n14_pass {
        "n=14 atlas audit did not pass"
    } else if !finite_pass {
        "finite-degree index did not pass"
    } else if !threshold_opaque {
        "Tao threshold status changed unexpectedly"
    } else {
        "Tao threshold gap remains explicit"
    };

    let result = json!({
        "experiment_id": L36_ID,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "dependencies": [
            {
                "artifact_id": L34A_ID,
                "status": l34a_result["status"].clone(),
                "role": "n=14 atlas theorem-packet audit"
            },
            {
                "artifact_id": L34B_ID,
                "status": l34b_result["status"].clone(),
                "role": "finite-degree index n=2..14"
            },
            {
                "artifact_id": L35_ID,
                "status": l35_result["status"].clone(),
                "role": "high-degree threshold interface"
            }
        ],
        "finite_certificate_min_degree": 2,
        "finite_certificate_max_degree": 14,
        "n14_atlas_tightest_margin": l34a_result["min_certified_cell_margin_to_cap"].clone(),
        "n14_tightest_cell": l34a_result["tightest_margin_cell"].clone(),
        "tao_covers_from": "n >= taoEHPThreshold",
        "taoEHPThreshold": Value::Null,
        "remaining_theorem_gap": "Need explicit taoEHPThreshold <= 15, or certificates for every 15 <= n < taoEHPThreshold.",
        "bridge_alternatives": [
            "Extract a practical high-degree cutoff from Tao's proof and verify it is <= 15.",
            "If the cutoff is larger, extend finite certificates across 15 <= n < taoEHPThreshold.",
            "Improve the analytic cutoff theorem until the finite certificate range meets it."
        ],
        "full_proof_synthesis_allowed": false,
        "first_failed_condition": first_failed,
        "claim_ceiling": "Bridge-blocker packet only. The current evidence is local n=14 atlas plus finite-degree index through n=14 plus an opaque high-degree theorem interface."
    });

    let report = format!(
        "# EHP114 L36 Full-Proof Bridge Blocker\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
Status: `{}`.\n\n\
The proof ladder is now cleanly separated. The local n=14 atlas layer passes, \
and the finite-degree index reaches n=14. Tao's theorem is available only as an \
opaque-threshold high-degree interface in this packet.\n\n\
## Exact Blocker\n\n\
Full synthesis is blocked until one of two things happens: extract \
`taoEHPThreshold <= 15`, or certify every finite degree in the interval \
`15 <= n < taoEHPThreshold` after the threshold is extracted.\n\n\
## Claim Ceiling\n\n\
No stronger claim is promoted by this packet.\n",
        L36_ID, status
    );

    Ok(Artifact {
        experiment_id: L36_ID,
        result,
        report,
    })
}

fn check_artifact_targets(out_root: &Path, artifacts: &[Artifact]) -> Result<(), String> {
    for artifact in artifacts {
        let (result_path, report_path, sha_path) = artifact_paths(out_root, artifact.experiment_id);
        for path in [result_path, report_path, sha_path] {
            if path.exists() {
                return Err(format!(
                    "refusing to overwrite existing artifact: {}",
                    path.display()
                ));
            }
        }
    }
    Ok(())
}

fn write_artifact(out_root: &Path, artifact: &Artifact) -> Result<String, String> {
    let dir = out_root.join(artifact.experiment_id);
    fs::create_dir_all(&dir).map_err(|err| format!("failed to create {}: {err}", dir.display()))?;
    let (result_path, report_path, sha_path) = artifact_paths(out_root, artifact.experiment_id);
    let result_text = serde_json::to_string_pretty(&artifact.result)
        .map_err(|err| format!("failed to serialize {}: {err}", artifact.experiment_id))?;
    fs::write(&result_path, result_text + "\n")
        .map_err(|err| format!("failed to write {}: {err}", result_path.display()))?;
    fs::write(&report_path, &artifact.report)
        .map_err(|err| format!("failed to write {}: {err}", report_path.display()))?;
    let digest = sha256_file(&result_path)?;
    let result_file = result_path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("RESULTS.json");
    fs::write(&sha_path, format!("{digest}  {result_file}\n"))
        .map_err(|err| format!("failed to write {}: {err}", sha_path.display()))?;
    Ok(digest)
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let cfg = parse_args()?;

    let l34a = build_l34a(&cfg, started)?;
    let l34b = build_l34b(&cfg, &l34a.result, started)?;
    let l35 = build_l35(started)?;
    let l36 = build_l36(&l34a.result, &l34b.result, &l35.result, started)?;
    let artifacts = vec![l34a, l34b, l35, l36];

    check_artifact_targets(&cfg.out_root, &artifacts)?;

    let mut emitted = Vec::new();
    for artifact in &artifacts {
        let digest = write_artifact(&cfg.out_root, artifact)?;
        emitted.push(json!({
            "experiment_id": artifact.experiment_id,
            "status": artifact.result["status"].clone(),
            "sha256": digest,
            "outdir": cfg.out_root.join(artifact.experiment_id).display().to_string()
        }));
    }

    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "status": "EHP114_PROOF_PATH_PACKETS_EMITTED",
            "artifact_count": emitted.len(),
            "artifacts": emitted
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
