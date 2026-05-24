//! EHP #114 n=14 residual-universe accounting certificate.
//!
//! This audit binds an upstream branch-isolation residual universe to the
//! downstream L21/L24/L27/L28/L29/L31 chain. It does not close geometry. Its
//! only job is to prove that a variable-size source-backed residual universe is
//! fully accounted for before a local cell packet accepts a non-64 residual
//! count.

use ehp_n3_poc::ehp114_n14_cell::{ensure_source_subcell_matches, CellSpec};
use serde::Serialize;
use serde_json::{json, Value};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const DEFAULT_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-RESIDUAL-UNIVERSE-ACCOUNTING-CERT-CELL-05-06-20260506-01";
const DEGREE: usize = 14;
const EPS: f64 = 0.1;

#[derive(Clone, Debug)]
struct Config {
    branch: PathBuf,
    normal: PathBuf,
    third: PathBuf,
    global: PathBuf,
    l24: PathBuf,
    l27: PathBuf,
    l28: PathBuf,
    l29: PathBuf,
    l31: PathBuf,
    out_dir: PathBuf,
    experiment_id: String,
    cell: CellSpec,
}

#[derive(Clone, Debug, Serialize)]
struct SourceSha {
    experiment_id: String,
    role: String,
    path: String,
    sha_status: String,
    expected_sha256: String,
    actual_sha256: String,
}

fn parse_args() -> Result<Config, String> {
    let cell = CellSpec::new(5, 6)?;
    let base = PathBuf::from("../../Erdos114/validated_length");
    let mut cfg = Config {
        branch: base.join("EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-CELL-05-06-20260506-01/EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-CELL-05-06-20260506-01_RESULTS.json"),
        normal: base.join("EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-CELL-05-06-20260506-01/EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-CELL-05-06-20260506-01_RESULTS.json"),
        third: base.join("EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-05-06-20260506-01/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-05-06-20260506-01_RESULTS.json"),
        global: base.join("EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-CELL-05-06-20260506-01/EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-CELL-05-06-20260506-01_RESULTS.json"),
        l24: base.join("EXP-MATH-EHP114-N14-VARIABLE-REGULAR-RESIDUAL-CLOSURE-CELL-05-06-20260506-01/EXP-MATH-EHP114-N14-VARIABLE-REGULAR-RESIDUAL-CLOSURE-CELL-05-06-20260506-01_RESULTS.json"),
        l27: base.join("EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-CELL-05-06-20260506-01/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-CELL-05-06-20260506-01_RESULTS.json"),
        l28: base.join("EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-CELL-05-06-20260506-01/EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-CELL-05-06-20260506-01_RESULTS.json"),
        l29: base.join("EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-CELL-05-06-20260506-01/EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-CELL-05-06-20260506-01_RESULTS.json"),
        l31: base.join("EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-CELL-05-06-20260506-01/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-CELL-05-06-20260506-01_RESULTS.json"),
        out_dir: base.join(DEFAULT_EXPERIMENT_ID),
        experiment_id: DEFAULT_EXPERIMENT_ID.to_string(),
        cell,
    };

    let args: Vec<String> = env::args().collect();
    let mut sub_i = cfg.cell.sub_i;
    let mut sub_j = cfg.cell.sub_j;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--branch" => {
                cfg.branch = next_path(&args, i, "--branch")?;
                i += 2;
            }
            "--normal" => {
                cfg.normal = next_path(&args, i, "--normal")?;
                i += 2;
            }
            "--third" => {
                cfg.third = next_path(&args, i, "--third")?;
                i += 2;
            }
            "--global" => {
                cfg.global = next_path(&args, i, "--global")?;
                i += 2;
            }
            "--l24" => {
                cfg.l24 = next_path(&args, i, "--l24")?;
                i += 2;
            }
            "--l27" => {
                cfg.l27 = next_path(&args, i, "--l27")?;
                i += 2;
            }
            "--l28" => {
                cfg.l28 = next_path(&args, i, "--l28")?;
                i += 2;
            }
            "--l29" => {
                cfg.l29 = next_path(&args, i, "--l29")?;
                i += 2;
            }
            "--l31" => {
                cfg.l31 = next_path(&args, i, "--l31")?;
                i += 2;
            }
            "--outdir" => {
                cfg.out_dir = next_path(&args, i, "--outdir")?;
                i += 2;
            }
            "--experiment-id" => {
                cfg.experiment_id = next_string(&args, i, "--experiment-id")?;
                i += 2;
            }
            "--sub-i" => {
                sub_i = next_string(&args, i, "--sub-i")?
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --sub-i: {err}"))?;
                i += 2;
            }
            "--sub-j" => {
                sub_j = next_string(&args, i, "--sub-j")?
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --sub-j: {err}"))?;
                i += 2;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
    }
    cfg.cell = CellSpec::new(sub_i, sub_j)?;
    Ok(cfg)
}

fn next_path(args: &[String], i: usize, flag: &str) -> Result<PathBuf, String> {
    if i + 1 >= args.len() {
        return Err(format!("{flag} requires a path"));
    }
    Ok(PathBuf::from(&args[i + 1]))
}

fn next_string(args: &[String], i: usize, flag: &str) -> Result<String, String> {
    if i + 1 >= args.len() {
        return Err(format!("{flag} requires a value"));
    }
    Ok(args[i + 1].clone())
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
        .ok_or_else(|| format!("invalid source path {}", path.display()))?;
    if !file_name.ends_with("_RESULTS.json") {
        return Err(format!(
            "source path must end in _RESULTS.json: {}",
            path.display()
        ));
    }
    Ok(path.with_file_name(file_name.replace("_RESULTS.json", "_RESULTS.sha256")))
}

fn verify_source_sha(role: &str, experiment_id: &str, path: &Path) -> Result<SourceSha, String> {
    let sha_path = sha_path_for(path)?;
    let sha_text = fs::read_to_string(&sha_path)
        .map_err(|err| format!("failed to read {}: {err}", sha_path.display()))?;
    let expected = sha_text
        .split_whitespace()
        .next()
        .ok_or_else(|| format!("empty sha file {}", sha_path.display()))?
        .to_string();
    let actual = sha256_file(path)?;
    let sha_status = if expected == actual { "PASS" } else { "FAIL" }.to_string();
    Ok(SourceSha {
        experiment_id: experiment_id.to_string(),
        role: role.to_string(),
        path: path.display().to_string(),
        sha_status,
        expected_sha256: expected,
        actual_sha256: actual,
    })
}

fn expect_str<'a>(v: &'a Value, key: &str) -> Result<&'a str, String> {
    v.get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("missing string field {key}"))
}

fn expect_usize(v: &Value, key: &str) -> Result<usize, String> {
    let n = v
        .get(key)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("missing integer field {key}"))?;
    usize::try_from(n).map_err(|_| format!("integer field {key} is too large"))
}

fn expect_f64(v: &Value, key: &str) -> Result<f64, String> {
    v.get(key)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("missing numeric field {key}"))
}

fn array_len(v: &Value, key: &str) -> Result<usize, String> {
    Ok(v.get(key)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("missing array field {key}"))?
        .len())
}

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12
}

fn unix_timestamp_string() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs()
        .to_string()
}

fn status_for(
    source_sha_fail_count: usize,
    coverage_failures: &[String],
    total: f64,
    cap: f64,
) -> &'static str {
    if source_sha_fail_count > 0 {
        "RESIDUAL_UNIVERSE_ACCOUNTING_FAIL_SOURCE_SHA"
    } else if !coverage_failures.is_empty() {
        "RESIDUAL_UNIVERSE_ACCOUNTING_FAIL_COVERAGE"
    } else if total > cap {
        "RESIDUAL_UNIVERSE_ACCOUNTING_FAIL_BUDGET"
    } else {
        "RESIDUAL_UNIVERSE_ACCOUNTING_PASS_NOT_GLOBAL_PROOF"
    }
}

fn write_report(result: &Value, path: &Path) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Residual-Universe Accounting Certificate\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Source-backed leaf universe count: `{}`\n\
- Upstream closed count: `{}`\n\
- Downstream residual count: `{}`\n\
- Downstream closed count: `{}`\n\
- Coverage equation pass: `{}`\n\
- Source SHA fail count: `{}`\n\
- Ownership duplicate count: `{}`\n\
- Source-filter mismatch count: `{}`\n\
- Candidate-count drift count: `{}`\n\
- Total validated length upper: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
This certificate proves only accounting over the emitted source artifacts. It \
binds the upstream branch-isolation exclusions to the downstream residual \
chain so that a variable residual universe can be accepted without pretending \
that the residual count must always be 64.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        result["status"].as_str().unwrap_or("UNKNOWN"),
        result["source_backed_leaf_universe_count"],
        result["upstream_closed_count"],
        result["downstream_residual_count"],
        result["downstream_closed_count"],
        result["coverage_equation_pass"],
        result["source_sha_fail_count"],
        result["ownership_duplicate_count"],
        result["source_filter_mismatch_count"],
        result["candidate_count_drift_count"],
        result["total_validated_length_upper"],
        result["exact_length_cap"],
        result["margin_to_cap"],
        result["first_failed_condition"]
            .as_str()
            .unwrap_or("unknown"),
        result["claim_ceiling"]
            .as_str()
            .unwrap_or("local accounting certificate only"),
    );
    fs::write(path, report).map_err(|err| format!("failed to write {}: {err}", path.display()))
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let cfg = parse_args()?;
    fs::create_dir_all(&cfg.out_dir)
        .map_err(|err| format!("failed to create {}: {err}", cfg.out_dir.display()))?;

    let result_path = cfg
        .out_dir
        .join(format!("{}_RESULTS.json", cfg.experiment_id));
    let report_path = cfg.out_dir.join(format!("{}_REPORT.md", cfg.experiment_id));
    let sha_path = cfg
        .out_dir
        .join(format!("{}_RESULTS.sha256", cfg.experiment_id));
    for path in [&result_path, &report_path, &sha_path] {
        if path.exists() {
            return Err(format!(
                "refusing to overwrite existing artifact: {}",
                path.display()
            ));
        }
    }

    let branch = read_json(&cfg.branch)?;
    let normal = read_json(&cfg.normal)?;
    let third = read_json(&cfg.third)?;
    let global = read_json(&cfg.global)?;
    let l24 = read_json(&cfg.l24)?;
    let l27 = read_json(&cfg.l27)?;
    let l28 = read_json(&cfg.l28)?;
    let l29 = read_json(&cfg.l29)?;
    let l31 = read_json(&cfg.l31)?;

    let source_subcell_contracts = json!({
        "branch": ensure_source_subcell_matches(&branch, cfg.cell)?,
        "normal": ensure_source_subcell_matches(&normal, cfg.cell)?,
        "third": ensure_source_subcell_matches(&third, cfg.cell)?,
        "global": ensure_source_subcell_matches(&global, cfg.cell)?,
        "l24": ensure_source_subcell_matches(&l24, cfg.cell)?,
        "l27": ensure_source_subcell_matches(&l27, cfg.cell)?,
        "l28": ensure_source_subcell_matches(&l28, cfg.cell)?,
        "l29": ensure_source_subcell_matches(&l29, cfg.cell)?,
        "l31": ensure_source_subcell_matches(&l31, cfg.cell)?,
    });

    let source_sha_status = vec![
        verify_source_sha(
            "branch_isolation",
            expect_str(&branch, "experiment_id")?,
            &cfg.branch,
        )?,
        verify_source_sha(
            "normal_collar",
            expect_str(&normal, "experiment_id")?,
            &cfg.normal,
        )?,
        verify_source_sha(
            "third_order_collar",
            expect_str(&third, "experiment_id")?,
            &cfg.third,
        )?,
        verify_source_sha(
            "global_critical_target",
            expect_str(&global, "experiment_id")?,
            &cfg.global,
        )?,
        verify_source_sha(
            "l24_regular_residual",
            expect_str(&l24, "experiment_id")?,
            &cfg.l24,
        )?,
        verify_source_sha(
            "l27_pprime_p8",
            expect_str(&l27, "experiment_id")?,
            &cfg.l27,
        )?,
        verify_source_sha(
            "l28_pprime_root_location",
            expect_str(&l28, "experiment_id")?,
            &cfg.l28,
        )?,
        verify_source_sha(
            "l29_sharp_pprime_root_location",
            expect_str(&l29, "experiment_id")?,
            &cfg.l29,
        )?,
        verify_source_sha(
            "l31_residual_chain",
            expect_str(&l31, "experiment_id")?,
            &cfg.l31,
        )?,
    ];
    let source_sha_fail_count = source_sha_status
        .iter()
        .filter(|row| row.sha_status != "PASS")
        .count();

    let mut coverage_failures = Vec::new();
    let upstream_excluded_count = array_len(&branch, "excluded_regions")?;
    let upstream_certified_count = array_len(&branch, "certified_branches")?;
    let upstream_closed_count = upstream_excluded_count + upstream_certified_count;
    let branch_remaining_count = expect_usize(&branch, "remaining_unresolved_branch_count")?;
    let normal_processed = expect_usize(&normal, "processed_region_count")?;
    let third_processed = expect_usize(&third, "processed_region_count")?;
    let global_source = expect_usize(&global, "source_region_count")?;
    let global_processed = expect_usize(&global, "processed_region_count")?;
    let global_unprocessed = expect_usize(&global, "source_unprocessed_region_count")?;
    let l24_processed = expect_usize(&l24, "processed_regular_region_count")?;
    let l27_processed = expect_usize(&l27, "processed_critical_candidate_count")?;
    let l31_expected = expect_usize(&l31, "expected_residual_chain_total_count")?;
    let l31_closed = expect_usize(&l31, "residual_chain_total_closed_count")?;
    let l31_regular = expect_usize(&l31, "regular_regions_closed_by_l24")?;
    let l31_critical = expect_usize(&l31, "final_critical_closed_count")?;
    let source_backed_leaf_universe_count = upstream_closed_count + branch_remaining_count;

    if branch_remaining_count != normal_processed {
        coverage_failures.push(format!(
            "branch remaining count {branch_remaining_count} does not equal normal processed count {normal_processed}"
        ));
    }
    if normal_processed != third_processed {
        coverage_failures.push(format!(
            "normal processed count {normal_processed} does not equal third-order processed count {third_processed}"
        ));
    }
    if third_processed != global_source || global_source != global_processed {
        coverage_failures.push(format!(
            "third/global source count mismatch: third {third_processed}, global source {global_source}, global processed {global_processed}"
        ));
    }
    if global_unprocessed != 0 {
        coverage_failures.push(format!(
            "global source_unprocessed_region_count is {global_unprocessed}, expected 0"
        ));
    }
    if global_processed != l24_processed + l27_processed {
        coverage_failures.push(format!(
            "global processed count {global_processed} does not equal L24 regular {l24_processed} + L27 critical {l27_processed}"
        ));
    }
    if l31_expected != l24_processed + l27_processed {
        coverage_failures.push(format!(
            "L31 expected count {l31_expected} does not equal L24 regular {l24_processed} + L27 critical {l27_processed}"
        ));
    }
    if l31_closed != l31_expected {
        coverage_failures.push(format!(
            "L31 closed count {l31_closed} does not equal expected count {l31_expected}"
        ));
    }
    if l31_closed != l31_regular + l31_critical {
        coverage_failures.push(format!(
            "L31 closed count {l31_closed} does not equal regular {l31_regular} + critical {l31_critical}"
        ));
    }
    if l31_closed != branch_remaining_count {
        coverage_failures.push(format!(
            "downstream closed count {l31_closed} does not equal branch remaining residual count {branch_remaining_count}"
        ));
    }
    if expect_usize(&l31, "source_sha_fail_count")? != 0 {
        coverage_failures.push("L31 source_sha_fail_count is not zero".to_string());
    }
    if expect_usize(&l31, "ownership_duplicate_count")? != 0 {
        coverage_failures.push("L31 ownership_duplicate_count is not zero".to_string());
    }
    if expect_usize(&l31, "source_filter_mismatch_count")? != 0 {
        coverage_failures.push("L31 source_filter_mismatch_count is not zero".to_string());
    }
    if expect_usize(&l31, "candidate_count_drift_count")? != 0 {
        coverage_failures.push("L31 candidate_count_drift_count is not zero".to_string());
    }
    if expect_str(&l31, "status")? != "RESIDUAL_CHAIN_INTEGRATION_PASS_NOT_GLOBAL_PROOF" {
        coverage_failures.push("L31 status is not pass".to_string());
    }
    if expect_usize(&branch, "ownership_duplicate_count")? != 0 {
        coverage_failures.push("branch ownership_duplicate_count is not zero".to_string());
    }

    let cap = expect_f64(&l31, "exact_length_cap")?;
    let total = expect_f64(&l31, "total_validated_length_upper")?;
    let margin = cap - total;
    for (label, source) in [
        ("branch", &branch),
        ("normal", &normal),
        ("third", &third),
        ("global", &global),
        ("l24", &l24),
        ("l27", &l27),
        ("l28", &l28),
        ("l29", &l29),
    ] {
        if !approx_eq(expect_f64(source, "exact_length_cap")?, cap) {
            coverage_failures.push(format!("{label} exact_length_cap drifts from L31"));
        }
        let source_total = expect_f64(source, "total_validated_length_upper")?;
        if source_total > total && !approx_eq(source_total, total) && label != "branch" {
            coverage_failures.push(format!(
                "{label} total_validated_length_upper {source_total} exceeds L31 total {total}"
            ));
        }
    }

    let coverage_equation_pass = source_backed_leaf_universe_count
        == upstream_closed_count + l31_closed
        && l31_closed == branch_remaining_count
        && coverage_failures.is_empty();
    let status = status_for(source_sha_fail_count, &coverage_failures, total, cap);
    let first_failed_condition = if source_sha_fail_count > 0 {
        "at least one source checksum failed".to_string()
    } else if let Some(first) = coverage_failures.first() {
        first.clone()
    } else if total > cap {
        "total_validated_length_upper exceeds exact_length_cap".to_string()
    } else {
        "none".to_string()
    };

    let result = json!({
        "experiment_id": cfg.experiment_id,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "degree": DEGREE,
        "eps": EPS,
        "subcell": cfg.cell,
        "cell_tag": cfg.cell.tag(),
        "source_subcell_contracts": source_subcell_contracts,
        "source_experiment_ids": {
            "branch": branch["experiment_id"].clone(),
            "normal": normal["experiment_id"].clone(),
            "third": third["experiment_id"].clone(),
            "global": global["experiment_id"].clone(),
            "l24": l24["experiment_id"].clone(),
            "l27": l27["experiment_id"].clone(),
            "l28": l28["experiment_id"].clone(),
            "l29": l29["experiment_id"].clone(),
            "l31": l31["experiment_id"].clone()
        },
        "source_paths": {
            "branch": cfg.branch.display().to_string(),
            "normal": cfg.normal.display().to_string(),
            "third": cfg.third.display().to_string(),
            "global": cfg.global.display().to_string(),
            "l24": cfg.l24.display().to_string(),
            "l27": cfg.l27.display().to_string(),
            "l28": cfg.l28.display().to_string(),
            "l29": cfg.l29.display().to_string(),
            "l31": cfg.l31.display().to_string()
        },
        "source_sha_status": source_sha_status,
        "source_sha_fail_count": source_sha_fail_count,
        "source_backed_leaf_universe_count": source_backed_leaf_universe_count,
        "upstream_excluded_count": upstream_excluded_count,
        "upstream_certified_count": upstream_certified_count,
        "upstream_closed_count": upstream_closed_count,
        "downstream_residual_count": branch_remaining_count,
        "downstream_expected_count": l31_expected,
        "downstream_closed_count": l31_closed,
        "global_source_region_count": global_source,
        "global_processed_region_count": global_processed,
        "global_source_unprocessed_region_count": global_unprocessed,
        "l24_regular_region_count": l24_processed,
        "l27_critical_candidate_count": l27_processed,
        "coverage_equation": "source_backed_leaf_universe_count = upstream_closed_count + downstream_closed_count",
        "coverage_equation_pass": coverage_equation_pass,
        "coverage_failure_count": coverage_failures.len(),
        "coverage_failures": coverage_failures,
        "ownership_duplicate_count": expect_usize(&l31, "ownership_duplicate_count")? + expect_usize(&branch, "ownership_duplicate_count")?,
        "source_filter_mismatch_count": expect_usize(&l31, "source_filter_mismatch_count")?,
        "candidate_count_drift_count": expect_usize(&l31, "candidate_count_drift_count")?,
        "source_accepted_length_upper": total,
        "total_validated_length_upper": total,
        "exact_length_cap": cap,
        "margin_to_cap": margin,
        "first_failed_condition": first_failed_condition,
        "packet_acceptance_policy": "A non-hard-cell packet may accept a non-64 L31 residual count only when this source-backed residual-universe accounting certificate passes and the downstream closed count equals the source-backed downstream residual count.",
        "claim_ceiling": "Local n=14 source-backed residual-universe accounting certificate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate."
    });

    let result_text = serde_json::to_string_pretty(&result)
        .map_err(|err| format!("failed to serialize result JSON: {err}"))?;
    fs::write(&result_path, result_text + "\n")
        .map_err(|err| format!("failed to write {}: {err}", result_path.display()))?;
    write_report(&result, &report_path)?;
    let digest = sha256_file(&result_path)?;
    fs::write(
        &sha_path,
        format!(
            "{}  {}\n",
            digest,
            result_path
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("RESULTS.json")
        ),
    )
    .map_err(|err| format!("failed to write {}: {err}", sha_path.display()))?;

    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "experiment_id": result["experiment_id"],
            "status": result["status"],
            "source_backed_leaf_universe_count": result["source_backed_leaf_universe_count"],
            "upstream_closed_count": result["upstream_closed_count"],
            "downstream_closed_count": result["downstream_closed_count"],
            "coverage_equation_pass": result["coverage_equation_pass"],
            "source_sha_fail_count": result["source_sha_fail_count"],
            "margin_to_cap": result["margin_to_cap"],
            "first_failed_condition": result["first_failed_condition"],
            "sha256": digest,
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
