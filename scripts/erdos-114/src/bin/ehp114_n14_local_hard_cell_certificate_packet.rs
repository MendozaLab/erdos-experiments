#![recursion_limit = "256"]

//! EHP #114 n=14 local hard-cell certificate packet.
//!
//! This is a theorem-facing assembler over the already verified hard-cell
//! residual chain. It does not recompute geometry and does not promote a global
//! n=14 or EHP #114 claim. Its job is to bind the local length budget, L31
//! integration result, source checksums, and claim ceiling in one auditable
//! packet before any scaling work.

use ehp_n3_poc::ehp114_n14_cell::{ensure_source_subcell_matches, CellSpec};
use serde::Serialize;
use serde_json::{json, Value};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const DEFAULT_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-LOCAL-HARD-CELL-CERTIFICATE-PACKET-20260506-01";
const L31_ID: &str = "EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01";
const DEGREE: usize = 14;
const EPS: f64 = 0.1;
const EXACT_LENGTH_CAP: f64 = 20.672796062619668;
const ACCEPTED_LENGTH_UPPER: f64 = 20.316451752723314;
const EXPECTED_RESIDUAL_CLOSURE_COUNT: usize = 64;

#[derive(Clone, Debug)]
struct Config {
    l31: PathBuf,
    residual_universe: Option<PathBuf>,
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
    let mut cfg = Config {
        l31: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01_RESULTS.json"),
        residual_universe: None,
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-HARD-CELL-CERTIFICATE-PACKET-20260506-01"),
        experiment_id: DEFAULT_EXPERIMENT_ID.to_string(),
        cell: CellSpec::hard_cell(),
    };

    let args: Vec<String> = env::args().collect();
    let mut sub_i = cfg.cell.sub_i;
    let mut sub_j = cfg.cell.sub_j;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--l31" => {
                cfg.l31 = next_path(&args, i, "--l31")?;
                i += 2;
            }
            "--residual-universe" => {
                cfg.residual_universe = Some(next_path(&args, i, "--residual-universe")?);
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

fn child_sha_rows_from_l31(l31: &Value) -> Result<Vec<SourceSha>, String> {
    let rows = l31
        .get("source_sha_status")
        .and_then(Value::as_array)
        .ok_or_else(|| "missing L31 source_sha_status array".to_string())?;
    let mut out = Vec::new();
    for row in rows {
        let path = expect_str(row, "path")?;
        let experiment_id = expect_str(row, "experiment_id")?;
        let role = match experiment_id {
            "EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01" => {
                "l24_regular_monotone"
            }
            "EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01" => {
                "l27_pprime_p8"
            }
            "EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01" => {
                "l28_pprime_root_location"
            }
            "EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01" => {
                "l29_sharp_pprime_root_location"
            }
            _ => "l31_child_source",
        };
        out.push(verify_source_sha(role, experiment_id, Path::new(path))?);
    }
    Ok(out)
}

fn validate_residual_universe(
    residual_universe: Option<&Value>,
    l31: &Value,
    cell: CellSpec,
) -> Result<Vec<String>, String> {
    let mut failures = Vec::new();
    let Some(universe) = residual_universe else {
        failures.push(
            "L31 residual count is not 64 and no residual-universe accounting certificate was supplied"
                .to_string(),
        );
        return Ok(failures);
    };
    if cell.is_hard_cell() {
        failures
            .push("hard-cell packet may not replace the fixed 64 residual universe".to_string());
    }
    if expect_str(universe, "status")? != "RESIDUAL_UNIVERSE_ACCOUNTING_PASS_NOT_GLOBAL_PROOF" {
        failures.push("residual-universe accounting certificate did not pass".to_string());
    }
    if !universe
        .get("coverage_equation_pass")
        .and_then(Value::as_bool)
        .unwrap_or(false)
    {
        failures.push("residual-universe coverage equation did not pass".to_string());
    }
    if expect_usize(universe, "source_sha_fail_count")? != 0 {
        failures.push("residual-universe source_sha_fail_count is not zero".to_string());
    }
    if expect_usize(universe, "ownership_duplicate_count")? != 0 {
        failures.push("residual-universe ownership_duplicate_count is not zero".to_string());
    }
    if expect_usize(universe, "source_filter_mismatch_count")? != 0 {
        failures.push("residual-universe source_filter_mismatch_count is not zero".to_string());
    }
    if expect_usize(universe, "candidate_count_drift_count")? != 0 {
        failures.push("residual-universe candidate_count_drift_count is not zero".to_string());
    }
    if expect_usize(universe, "downstream_closed_count")?
        != expect_usize(l31, "residual_chain_total_closed_count")?
    {
        failures.push(
            "residual-universe downstream_closed_count does not equal L31 residual_chain_total_closed_count"
                .to_string(),
        );
    }
    if expect_usize(universe, "downstream_expected_count")?
        != expect_usize(l31, "expected_residual_chain_total_count")?
    {
        failures.push(
            "residual-universe downstream_expected_count does not equal L31 expected_residual_chain_total_count"
                .to_string(),
        );
    }
    if !approx_eq(
        expect_f64(universe, "total_validated_length_upper")?,
        expect_f64(l31, "total_validated_length_upper")?,
    ) {
        failures.push(
            "residual-universe total_validated_length_upper does not equal L31 total".to_string(),
        );
    }
    Ok(failures)
}

fn validate_l31(
    l31: &Value,
    residual_universe: Option<&Value>,
    cell: CellSpec,
) -> Result<Vec<String>, String> {
    let mut failures = Vec::new();
    if cell.is_hard_cell() && expect_str(l31, "experiment_id")? != L31_ID {
        failures.push("L31 experiment_id mismatch".to_string());
    }
    if expect_str(l31, "status")? != "RESIDUAL_CHAIN_INTEGRATION_PASS_NOT_GLOBAL_PROOF" {
        failures.push("L31 status is not the expected local pass status".to_string());
    }
    if expect_usize(l31, "source_sha_fail_count")? != 0 {
        failures.push("L31 reports source_sha_fail_count > 0".to_string());
    }
    if expect_usize(l31, "ownership_duplicate_count")? != 0 {
        failures.push("L31 reports ownership duplicates".to_string());
    }
    if expect_usize(l31, "source_filter_mismatch_count")? != 0 {
        failures.push("L31 reports source-filter mismatches".to_string());
    }
    if expect_usize(l31, "candidate_count_drift_count")? != 0 {
        failures.push("L31 reports candidate-count drift".to_string());
    }
    let l31_closed = expect_usize(l31, "residual_chain_total_closed_count")?;
    let l31_expected = expect_usize(l31, "expected_residual_chain_total_count")?;
    if l31_closed == EXPECTED_RESIDUAL_CLOSURE_COUNT
        && l31_expected == EXPECTED_RESIDUAL_CLOSURE_COUNT
    {
        if residual_universe.is_some() {
            failures.push(
                "residual-universe accounting certificate supplied but fixed 64 residual count already holds"
                    .to_string(),
            );
        }
    } else {
        failures.extend(validate_residual_universe(residual_universe, l31, cell)?);
    }
    if cell.is_hard_cell() {
        if !approx_eq(
            expect_f64(l31, "total_validated_length_upper")?,
            ACCEPTED_LENGTH_UPPER,
        ) {
            failures
                .push("L31 total_validated_length_upper drifted from local baseline".to_string());
        }
        if !approx_eq(
            expect_f64(l31, "source_accepted_length_upper")?,
            ACCEPTED_LENGTH_UPPER,
        ) {
            failures
                .push("L31 source_accepted_length_upper drifted from local baseline".to_string());
        }
    }
    if !approx_eq(expect_f64(l31, "exact_length_cap")?, EXACT_LENGTH_CAP) {
        failures.push("L31 exact_length_cap drifted".to_string());
    }
    if expect_f64(l31, "total_validated_length_upper")? > expect_f64(l31, "exact_length_cap")? {
        failures.push("L31 local length upper exceeds cap".to_string());
    }
    Ok(failures)
}

fn status_for(source_sha_fail_count: usize, l31_failures: &[String], margin: f64) -> &'static str {
    if source_sha_fail_count > 0 {
        "LOCAL_HARD_CELL_FAIL_SOURCE_SHA"
    } else if !l31_failures.is_empty() {
        "LOCAL_HARD_CELL_FAIL_RESIDUAL_CHAIN"
    } else if margin < 0.0 {
        "LOCAL_HARD_CELL_FAIL_BUDGET"
    } else {
        "LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF"
    }
}

fn first_failed_condition(
    source_sha_fail_count: usize,
    l31_failures: &[String],
    margin: f64,
) -> String {
    if source_sha_fail_count > 0 {
        "at least one source checksum failed".to_string()
    } else if let Some(first) = l31_failures.first() {
        first.clone()
    } else if margin < 0.0 {
        "local total_validated_length_upper exceeds exact_length_cap".to_string()
    } else {
        "none".to_string()
    }
}

fn write_report(result: &Value, path: &Path) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Local Hard-Cell Certificate Packet\n\n\
Experiment: `{}`\n\n\
## Local Theorem Statement\n\n\
For the requested n=14 root-affine subcell at `eps = 0.1`, the certified \
local hard-cell length upper assembled from the accepted slab length and the \
L24/L27/L28/L29 residual-chain closure is at most `{}`. This is below the \
local exact cap `{}` by margin `{}`.\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Local certificate pass: `{}`\n\
- Source SHA fail count: `{}`\n\
- L31 residual-chain failures: `{}`\n\
- Residual closures verified by L31: `{}` of `{}`\n\
- Source-backed leaf universe count: `{}`\n\
- Upstream closed count: `{}`\n\
- Downstream closed count: `{}`\n\
- Coverage equation pass: `{}`\n\
- Ownership duplicate count: `{}`\n\
- Source-filter mismatch count: `{}`\n\
- Candidate-count drift count: `{}`\n\
- First failed condition: `{}`\n\n\
## Dependency Table\n\n\
| Role | Experiment |\n\
|---|---|\n\
| L24 regular monotone closure | `EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01` |\n\
| L27 p-prime partition closure | `EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01` |\n\
| L28 p-prime root-location closure | `EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01` |\n\
| L29 sharp p-prime root-location closure | `EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01` |\n\
| L31 residual-chain integration | `EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01` |\n\n\
## Interpretation\n\n\
This packet is a local theorem-facing wrapper around verified artifacts. It \
does not recompute interval geometry. It checks that the L31 integration result \
and all child source checksums are stable, then states the exact local claim \
that follows. The next proof-facing dependency is global n=14 coverage over all \
root-affine cells using the same residual-chain discipline.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"].as_str().unwrap_or(DEFAULT_EXPERIMENT_ID),
        result["total_validated_length_upper"],
        result["exact_length_cap"],
        result["margin_to_cap"],
        result["status"].as_str().unwrap_or("UNKNOWN"),
        result["local_hard_cell_certificate_pass"],
        result["source_sha_fail_count"],
        result["l31_failure_count"],
        result["residual_chain_total_closed_count"],
        result["expected_residual_chain_total_count"],
        result["source_backed_leaf_universe_count"],
        result["upstream_closed_count"],
        result["downstream_closed_count"],
        result["coverage_equation_pass"],
        result["ownership_duplicate_count"],
        result["source_filter_mismatch_count"],
        result["candidate_count_drift_count"],
        result["first_failed_condition"].as_str().unwrap_or("unknown"),
        result["claim_ceiling"].as_str().unwrap_or("local certificate packet only"),
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

    let l31 = read_json(&cfg.l31)?;
    let l31_source_subcell_contract = ensure_source_subcell_matches(&l31, cfg.cell)?;
    let residual_universe = match &cfg.residual_universe {
        Some(path) => Some(read_json(path)?),
        None => None,
    };
    let residual_universe_source_subcell_contract = match &residual_universe {
        Some(value) => Some(ensure_source_subcell_matches(value, cfg.cell)?),
        None => None,
    };
    let l31_actual_id = expect_str(&l31, "experiment_id")?.to_string();
    let mut source_sha_status = Vec::new();
    source_sha_status.push(verify_source_sha(
        "l31_residual_chain_integration",
        &l31_actual_id,
        &cfg.l31,
    )?);
    if let (Some(path), Some(value)) = (&cfg.residual_universe, &residual_universe) {
        source_sha_status.push(verify_source_sha(
            "residual_universe_accounting",
            expect_str(value, "experiment_id")?,
            path,
        )?);
    }
    source_sha_status.extend(child_sha_rows_from_l31(&l31)?);
    let source_sha_fail_count = source_sha_status
        .iter()
        .filter(|row| row.sha_status != "PASS")
        .count();

    let l31_failures = validate_l31(&l31, residual_universe.as_ref(), cfg.cell)?;
    let total = expect_f64(&l31, "total_validated_length_upper")?;
    let cap = expect_f64(&l31, "exact_length_cap")?;
    let margin = cap - total;
    let status = status_for(source_sha_fail_count, &l31_failures, margin);
    let first_failed = first_failed_condition(source_sha_fail_count, &l31_failures, margin);

    let result = json!({
        "experiment_id": cfg.experiment_id,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "degree": DEGREE,
        "eps": EPS,
        "subcell": cfg.cell,
        "cell_tag": cfg.cell.tag(),
        "source_subcell_contracts": {
            "l31": l31_source_subcell_contract,
            "residual_universe": residual_universe_source_subcell_contract
        },
        "theorem_statement": format!("For local n=14 root-affine subcell {} at eps=0.1, the assembled validated length upper is <= {}, hence below the local exact cap {}, provided the cited source artifacts are accepted.", cfg.cell.tag(), total, cap),
        "theorem_scope": format!("local certificate for n=14 subcell {} only", cfg.cell.tag()),
        "source_experiment_ids": {
            "l31_residual_chain_integration": l31_actual_id,
            "residual_universe_accounting": residual_universe
                .as_ref()
                .map(|value| value["experiment_id"].clone())
                .unwrap_or(Value::Null),
            "l24_regular_monotone": l31["source_experiment_ids"]["l24_regular_monotone"].clone(),
            "l27_pprime_p8": l31["source_experiment_ids"]["l27_pprime_p8"].clone(),
            "l28_pprime_root_location": l31["source_experiment_ids"]["l28_pprime_root_location"].clone(),
            "l29_sharp_pprime_root_location": l31["source_experiment_ids"]["l29_sharp_pprime_root_location"].clone()
        },
        "source_paths": {
            "l31": cfg.l31.display().to_string(),
            "residual_universe": cfg.residual_universe
                .as_ref()
                .map(|path| Value::String(path.display().to_string()))
                .unwrap_or(Value::Null),
            "l24": l31["source_paths"]["l24"].clone(),
            "l27": l31["source_paths"]["l27"].clone(),
            "l28": l31["source_paths"]["l28"].clone(),
            "l29": l31["source_paths"]["l29"].clone()
        },
        "source_sha_status": source_sha_status,
        "source_sha_fail_count": source_sha_fail_count,
        "source_residual_chain_status": l31["status"].clone(),
        "l31_failure_count": l31_failures.len(),
        "l31_failures": l31_failures,
        "proof_obligations": [
            "source checksum verification",
            "residual-chain pass status",
            "regular-plus-critical closure count equals 64 or source-backed residual-universe accounting certificate passes",
            "no ownership duplicates",
            "no source-filter mismatches",
            "no candidate-count drift",
            "local accepted length upper remains below exact cap"
        ],
        "closure_counts": l31["closure_counts"].clone(),
        "regular_regions_closed_by_l24": l31["regular_regions_closed_by_l24"].clone(),
        "critical_regions_closed_by_l27_partition": l31["critical_regions_closed_by_l27_partition"].clone(),
        "critical_regions_closed_by_l28_root_location": l31["critical_regions_closed_by_l28_root_location"].clone(),
        "critical_regions_closed_by_l29_sharp_root_location": l31["critical_regions_closed_by_l29_sharp_root_location"].clone(),
        "final_critical_closed_count": l31["final_critical_closed_count"].clone(),
        "residual_chain_total_closed_count": l31["residual_chain_total_closed_count"].clone(),
        "expected_residual_chain_total_count": l31["expected_residual_chain_total_count"].clone(),
        "residual_universe_accounting_status": residual_universe
            .as_ref()
            .map(|value| value["status"].clone())
            .unwrap_or(Value::Null),
        "source_backed_leaf_universe_count": residual_universe
            .as_ref()
            .map(|value| value["source_backed_leaf_universe_count"].clone())
            .unwrap_or(Value::Null),
        "upstream_closed_count": residual_universe
            .as_ref()
            .map(|value| value["upstream_closed_count"].clone())
            .unwrap_or(Value::Null),
        "downstream_closed_count": residual_universe
            .as_ref()
            .map(|value| value["downstream_closed_count"].clone())
            .unwrap_or(Value::Null),
        "coverage_equation_pass": residual_universe
            .as_ref()
            .map(|value| value["coverage_equation_pass"].clone())
            .unwrap_or(Value::Null),
        "ownership_duplicate_count": l31["ownership_duplicate_count"].clone(),
        "source_filter_mismatch_count": l31["source_filter_mismatch_count"].clone(),
        "candidate_count_drift_count": l31["candidate_count_drift_count"].clone(),
        "budget_pass": margin >= 0.0,
        "source_accepted_length_upper": l31["source_accepted_length_upper"].clone(),
        "total_validated_length_upper": total,
        "exact_length_cap": cap,
        "margin_to_cap": margin,
        "local_hard_cell_certificate_pass": status == "LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF",
        "first_failed_condition": first_failed,
        "next_dependency": "Global n=14 atlas certificate over all root-affine cells; this local packet is necessary but not sufficient for global n=14 or full EHP114.",
        "claim_ceiling": "Local n=14 cell certificate packet only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. Operator-grammar and MDL sidecars are not evidence for this packet."
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
            "local_hard_cell_certificate_pass": result["local_hard_cell_certificate_pass"],
            "source_sha_fail_count": result["source_sha_fail_count"],
            "residual_chain_total_closed_count": result["residual_chain_total_closed_count"],
            "total_validated_length_upper": result["total_validated_length_upper"],
            "exact_length_cap": result["exact_length_cap"],
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
