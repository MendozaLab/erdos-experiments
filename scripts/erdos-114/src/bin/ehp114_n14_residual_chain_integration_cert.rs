//! EHP #114 n=14 residual-chain integration certificate.
//!
//! This is a bookkeeping certificate over existing local hard-cell artifacts.
//! It does not recompute geometry or promote length. It verifies that the
//! L24/L27/L28/L29 residual closures cover the intended source-filter chain
//! without count drift, ownership duplication, or source-filter gaps.

use ehp_n3_poc::ehp114_n14_cell::{ensure_source_subcell_matches, CellSpec};
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const DEFAULT_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01";
const DEGREE: usize = 14;
const EPS: f64 = 0.1;

const L24_ID: &str = "EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01";
const L27_ID: &str =
    "EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01";
const L28_ID: &str = "EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01";
const L29_ID: &str = "EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01";

#[derive(Clone, Debug)]
struct Config {
    l24: PathBuf,
    l27: PathBuf,
    l28: PathBuf,
    l29: PathBuf,
    out_dir: PathBuf,
    experiment_id: String,
    cell: CellSpec,
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize)]
struct RegionKey {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
}

#[derive(Clone, Debug, Serialize)]
struct SourceSha {
    experiment_id: String,
    path: String,
    sha_status: String,
    expected_sha256: String,
    actual_sha256: String,
}

fn parse_args() -> Result<Config, String> {
    let mut cfg = Config {
        l24: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01_RESULTS.json"),
        l27: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01_RESULTS.json"),
        l28: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01_RESULTS.json"),
        l29: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01_RESULTS.json"),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01"),
        experiment_id: DEFAULT_EXPERIMENT_ID.to_string(),
        cell: CellSpec::hard_cell(),
    };

    let args: Vec<String> = env::args().collect();
    let mut sub_i = cfg.cell.sub_i;
    let mut sub_j = cfg.cell.sub_j;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
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

fn verify_source_sha(experiment_id: &str, path: &Path) -> Result<SourceSha, String> {
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

fn region_key(row: &Value) -> Result<RegionKey, String> {
    Ok(RegionKey {
        source_index: expect_usize(row, "source_index")?,
        source_ownership_key: expect_str(row, "source_ownership_key")?.to_string(),
        split_path: expect_str(row, "split_path")?.to_string(),
    })
}

fn array_field<'a>(v: &'a Value, key: &str) -> Result<&'a Vec<Value>, String> {
    v.get(key)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("missing array field {key}"))
}

fn keys_with_status(
    v: &Value,
    array_key: &str,
    status_key: &str,
    status: &str,
) -> Result<BTreeSet<RegionKey>, String> {
    let mut out = BTreeSet::new();
    for row in array_field(v, array_key)? {
        if expect_str(row, status_key)? == status {
            out.insert(region_key(row)?);
        }
    }
    Ok(out)
}

fn keys_with_statuses(
    v: &Value,
    array_key: &str,
    status_key: &str,
    statuses: &[&str],
) -> Result<BTreeSet<RegionKey>, String> {
    let mut out = BTreeSet::new();
    for row in array_field(v, array_key)? {
        let row_status = expect_str(row, status_key)?;
        if statuses.iter().any(|status| *status == row_status) {
            out.insert(region_key(row)?);
        }
    }
    Ok(out)
}

fn all_keys(v: &Value, array_key: &str) -> Result<BTreeSet<RegionKey>, String> {
    let mut out = BTreeSet::new();
    for row in array_field(v, array_key)? {
        out.insert(region_key(row)?);
    }
    Ok(out)
}

fn duplicate_count(sets: &[&BTreeSet<RegionKey>]) -> usize {
    let mut seen = BTreeSet::new();
    let mut duplicates = 0;
    for set in sets {
        for key in set.iter() {
            if !seen.insert(key.clone()) {
                duplicates += 1;
            }
        }
    }
    duplicates
}

fn missing_keys(expected: &BTreeSet<RegionKey>, actual: &BTreeSet<RegionKey>) -> Vec<RegionKey> {
    expected.difference(actual).cloned().collect()
}

fn extra_keys(expected: &BTreeSet<RegionKey>, actual: &BTreeSet<RegionKey>) -> Vec<RegionKey> {
    actual.difference(expected).cloned().collect()
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
    sha_fail_count: usize,
    count_drift_count: usize,
    source_filter_mismatch_count: usize,
    ownership_duplicate_count: usize,
    total: f64,
    cap: f64,
) -> &'static str {
    if sha_fail_count > 0 {
        "RESIDUAL_CHAIN_FAIL_SOURCE_SHA"
    } else if count_drift_count > 0 {
        "RESIDUAL_CHAIN_FAIL_COUNT_DRIFT"
    } else if source_filter_mismatch_count > 0 {
        "RESIDUAL_CHAIN_FAIL_SOURCE_FILTER"
    } else if ownership_duplicate_count > 0 {
        "RESIDUAL_CHAIN_FAIL_OWNERSHIP"
    } else if total > cap {
        "RESIDUAL_CHAIN_FAIL_BUDGET"
    } else {
        "RESIDUAL_CHAIN_INTEGRATION_PASS_NOT_GLOBAL_PROOF"
    }
}

fn first_failed_condition(
    sha_fail_count: usize,
    count_drifts: &[String],
    source_mismatches: &[String],
    ownership_duplicate_count: usize,
    total: f64,
    cap: f64,
) -> String {
    if sha_fail_count > 0 {
        "at least one source checksum failed".to_string()
    } else if let Some(first) = count_drifts.first() {
        first.clone()
    } else if let Some(first) = source_mismatches.first() {
        first.clone()
    } else if ownership_duplicate_count > 0 {
        format!("{ownership_duplicate_count} duplicate ownership keys in final closure sets")
    } else if total > cap {
        format!("total_validated_length_upper {total} exceeds exact_length_cap {cap}")
    } else {
        "none".to_string()
    }
}

fn write_report(result: &Value, path: &Path) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Residual-Chain Integration Certificate\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Regular regions closed by L24: `{}`\n\
- Critical regions closed by L27 partition closure: `{}`\n\
- Critical regions closed by L28 root-location: `{}`\n\
- Critical regions closed by L29 sharp root-location: `{}`\n\
- Final critical closed count: `{}` of `{}`\n\
- Total residual closure count: `{}`\n\
- Ownership duplicate count: `{}`\n\
- Source-filter mismatch count: `{}`\n\
- Candidate-count drift count: `{}`\n\
- Total validated length upper: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
This certificate does not recompute the analytic inequalities. It verifies that \
the already emitted L24, L27, L28, and L29 local certificates compose into a \
single residual-chain closure for the n=14 hard cell `(6,4)`. The pass condition \
is bookkeeping integrity: source checksums, source filters, candidate counts, \
ownership uniqueness, and length-cap consistency must all agree.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        result["status"].as_str().unwrap_or("UNKNOWN"),
        result["regular_regions_closed_by_l24"],
        result["critical_regions_closed_by_l27_partition"],
        result["critical_regions_closed_by_l28_root_location"],
        result["critical_regions_closed_by_l29_sharp_root_location"],
        result["final_critical_closed_count"],
        result["critical_candidate_count"],
        result["residual_chain_total_closed_count"],
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
            .unwrap_or("local diagnostic only"),
    );
    fs::write(path, report).map_err(|err| format!("failed to write {}: {err}", path.display()))
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let cfg = parse_args()?;
    fs::create_dir_all(&cfg.out_dir)
        .map_err(|err| format!("failed to create {}: {err}", cfg.out_dir.display()))?;

    let l24 = read_json(&cfg.l24)?;
    let l27 = read_json(&cfg.l27)?;
    let l28 = read_json(&cfg.l28)?;
    let l29 = read_json(&cfg.l29)?;
    let l24_source_subcell_contract = ensure_source_subcell_matches(&l24, cfg.cell)?;
    let l27_source_subcell_contract = ensure_source_subcell_matches(&l27, cfg.cell)?;
    let l28_source_subcell_contract = ensure_source_subcell_matches(&l28, cfg.cell)?;
    let l29_source_subcell_contract = ensure_source_subcell_matches(&l29, cfg.cell)?;
    let l24_actual_id = expect_str(&l24, "experiment_id")?.to_string();
    let l27_actual_id = expect_str(&l27, "experiment_id")?.to_string();
    let l28_actual_id = expect_str(&l28, "experiment_id")?.to_string();
    let l29_actual_id = expect_str(&l29, "experiment_id")?.to_string();

    let source_sha_status = vec![
        verify_source_sha(&l24_actual_id, &cfg.l24)?,
        verify_source_sha(&l27_actual_id, &cfg.l27)?,
        verify_source_sha(&l28_actual_id, &cfg.l28)?,
        verify_source_sha(&l29_actual_id, &cfg.l29)?,
    ];
    let sha_fail_count = source_sha_status
        .iter()
        .filter(|row| row.sha_status != "PASS")
        .count();

    let mut count_drifts = Vec::new();
    let mut source_mismatches = Vec::new();

    if cfg.cell.is_hard_cell() && l24_actual_id != L24_ID {
        count_drifts.push("L24 experiment_id mismatch".to_string());
    }
    if cfg.cell.is_hard_cell() && l27_actual_id != L27_ID {
        count_drifts.push("L27 experiment_id mismatch".to_string());
    }
    if cfg.cell.is_hard_cell() && l28_actual_id != L28_ID {
        count_drifts.push("L28 experiment_id mismatch".to_string());
    }
    if cfg.cell.is_hard_cell() && l29_actual_id != L29_ID {
        count_drifts.push("L29 experiment_id mismatch".to_string());
    }
    if expect_str(&l28, "source_experiment_id")? != l27_actual_id {
        source_mismatches.push("L28 source_experiment_id does not point to L27".to_string());
    }
    if expect_str(&l29, "source_experiment_id")? != l28_actual_id {
        source_mismatches.push("L29 source_experiment_id does not point to L28".to_string());
    }

    let l24_closed = keys_with_statuses(
        &l24,
        "certificates",
        "status",
        &["EXCLUDED_MONOTONE_X_CHART", "CERTIFIED_REGULAR_X_CHART"],
    )?;
    let l27_all = all_keys(&l27, "region_summaries")?;
    let l27_closed = keys_with_status(
        &l27,
        "region_summaries",
        "region_status",
        "AFFINE_PARTITION_CLOSED_REGION",
    )?;
    let l27_still = keys_with_status(
        &l27,
        "region_summaries",
        "region_status",
        "STILL_CRITICAL_REGION",
    )?;
    let l28_all = all_keys(&l28, "region_summaries")?;
    let l28_closed = keys_with_status(
        &l28,
        "region_summaries",
        "region_status",
        "PPRIME_ROOT_EXCLUDED",
    )?;
    let l28_still = keys_with_status(
        &l28,
        "region_summaries",
        "region_status",
        "STILL_PPRIME_ROOT_NEAR",
    )?;
    let l29_all = all_keys(&l29, "region_summaries")?;
    let l29_closed = keys_with_status(
        &l29,
        "region_summaries",
        "region_status",
        "SHARP_PPRIME_ROOT_EXCLUDED",
    )?;
    let l29_still = keys_with_status(
        &l29,
        "region_summaries",
        "region_status",
        "STILL_SHARP_PPRIME_ROOT_NEAR",
    )?;

    let l24_processed = expect_usize(&l24, "processed_regular_region_count")?;
    let l24_source_regular = expect_usize(&l24, "source_regular_region_count")?;
    let l24_certified = expect_usize(&l24, "regular_regions_certified")?;
    let l24_excluded = expect_usize(&l24, "regular_regions_excluded")?;
    if l24_processed != l24_source_regular {
        count_drifts.push(format!(
            "L24 processed_regular_region_count {l24_processed} does not equal source_regular_region_count {l24_source_regular}"
        ));
    }
    if l24_certified + l24_excluded != l24_closed.len() {
        count_drifts.push(
            "L24 certified+excluded regular count does not equal closed certificate set size"
                .to_string(),
        );
    }
    if l24_closed.len() != l24_processed {
        count_drifts.push(format!(
            "L24 closed {} of {} processed regular regions",
            l24_closed.len(),
            l24_processed
        ));
    }
    if expect_usize(&l27, "processed_critical_candidate_count")? != l27_all.len() {
        count_drifts.push(
            "L27 processed_critical_candidate_count does not equal region_summaries length"
                .to_string(),
        );
    }
    if expect_usize(&l27, "affine_partition_closed_count")? != l27_closed.len() {
        count_drifts
            .push("L27 affine_partition_closed_count does not equal closed set size".to_string());
    }
    if expect_usize(&l27, "still_critical_candidate_count")? != l27_still.len() {
        count_drifts.push(
            "L27 still_critical_candidate_count does not equal still-critical set size".to_string(),
        );
    }
    if expect_usize(&l28, "processed_remaining_candidate_count")? != l28_all.len() {
        count_drifts.push(
            "L28 processed_remaining_candidate_count does not equal region_summaries length"
                .to_string(),
        );
    }
    if expect_usize(&l28, "pprime_root_excluded_count")? != l28_closed.len() {
        count_drifts
            .push("L28 pprime_root_excluded_count does not equal excluded set size".to_string());
    }
    if expect_usize(&l28, "still_unresolved_count")? != l28_still.len() {
        count_drifts
            .push("L28 still_unresolved_count does not equal still-near set size".to_string());
    }
    if expect_usize(&l29, "processed_remaining_candidate_count")? != l29_all.len() {
        count_drifts.push(
            "L29 processed_remaining_candidate_count does not equal region_summaries length"
                .to_string(),
        );
    }
    if expect_usize(&l29, "sharp_pprime_root_excluded_count")? != l29_closed.len() {
        count_drifts.push(
            "L29 sharp_pprime_root_excluded_count does not equal excluded set size".to_string(),
        );
    }
    if expect_usize(&l29, "still_unresolved_count")? != 0 || !l29_still.is_empty() {
        count_drifts.push("L29 still has unresolved sharp p-prime root-near regions".to_string());
    }
    if expect_usize(&l29, "unresolved_leaf_count")? != 0 {
        count_drifts.push("L29 unresolved_leaf_count is not zero".to_string());
    }

    if l28_all != l27_still {
        source_mismatches.push(format!(
            "L28 region set does not equal L27 still-critical set (missing {}, extra {})",
            missing_keys(&l27_still, &l28_all).len(),
            extra_keys(&l27_still, &l28_all).len()
        ));
    }
    if l29_all != l28_still {
        source_mismatches.push(format!(
            "L29 region set does not equal L28 still-near set (missing {}, extra {})",
            missing_keys(&l28_still, &l29_all).len(),
            extra_keys(&l28_still, &l29_all).len()
        ));
    }

    let mut final_critical = BTreeSet::new();
    final_critical.extend(l27_closed.iter().cloned());
    final_critical.extend(l28_closed.iter().cloned());
    final_critical.extend(l29_closed.iter().cloned());
    if final_critical != l27_all {
        source_mismatches.push(format!(
            "final critical closure set does not equal L27 processed set (missing {}, extra {})",
            missing_keys(&l27_all, &final_critical).len(),
            extra_keys(&l27_all, &final_critical).len()
        ));
    }

    let ownership_duplicate_count =
        duplicate_count(&[&l24_closed, &l27_closed, &l28_closed, &l29_closed]);
    let residual_chain_total_closed_count = l24_closed.len() + final_critical.len();
    let expected_residual_chain_total_count = l24_processed + l27_all.len();
    if residual_chain_total_closed_count != expected_residual_chain_total_count {
        count_drifts.push(format!(
            "residual_chain_total_closed_count {residual_chain_total_closed_count} does not equal expected {expected_residual_chain_total_count}"
        ));
    }

    let expected_cap = expect_f64(&l24, "exact_length_cap")?;
    let expected_total = expect_f64(&l24, "total_validated_length_upper")?;
    for (label, source) in [("L24", &l24), ("L27", &l27), ("L28", &l28), ("L29", &l29)] {
        let cap = expect_f64(source, "exact_length_cap")?;
        let total = expect_f64(source, "total_validated_length_upper")?;
        if !approx_eq(cap, expected_cap) {
            count_drifts.push(format!("{label} exact_length_cap drifted: {cap}"));
        }
        if label == "L24" && !approx_eq(total, expected_total) {
            count_drifts.push(format!(
                "{label} total_validated_length_upper drifted: {total}"
            ));
        } else if label != "L24" && total > expected_total && !approx_eq(total, expected_total) {
            count_drifts.push(format!(
                "{label} total_validated_length_upper exceeds L24 repaired total: {total}"
            ));
        }
        if let Some(n) = source
            .get("ownership_duplicate_count")
            .and_then(Value::as_u64)
        {
            if n != 0 {
                count_drifts.push(format!("{label} reports ownership_duplicate_count {n}"));
            }
        }
    }

    let mut closure_counts = BTreeMap::new();
    closure_counts.insert(
        "L24_regular_certified_or_excluded_regions",
        l24_closed.len(),
    );
    closure_counts.insert("L27_partition_closed_critical_regions", l27_closed.len());
    closure_counts.insert("L28_pprime_root_exclusions", l28_closed.len());
    closure_counts.insert("L29_sharp_pprime_root_exclusions", l29_closed.len());

    let candidate_count_drift_count = count_drifts.len();
    let source_filter_mismatch_count = source_mismatches.len();
    let total = expected_total;
    let cap = expected_cap;
    let margin = cap - total;
    let status = status_for(
        sha_fail_count,
        candidate_count_drift_count,
        source_filter_mismatch_count,
        ownership_duplicate_count,
        total,
        cap,
    );
    let first_failed = first_failed_condition(
        sha_fail_count,
        &count_drifts,
        &source_mismatches,
        ownership_duplicate_count,
        total,
        cap,
    );

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
            "l24": l24_source_subcell_contract,
            "l27": l27_source_subcell_contract,
            "l28": l28_source_subcell_contract,
            "l29": l29_source_subcell_contract
        },
        "source_experiment_ids": {
            "l24_regular_monotone": l24_actual_id,
            "l27_pprime_p8": l27_actual_id,
            "l28_pprime_root_location": l28_actual_id,
            "l29_sharp_pprime_root_location": l29_actual_id
        },
        "source_paths": {
            "l24": cfg.l24.display().to_string(),
            "l27": cfg.l27.display().to_string(),
            "l28": cfg.l28.display().to_string(),
            "l29": cfg.l29.display().to_string()
        },
        "source_sha_status": source_sha_status,
        "source_sha_fail_count": sha_fail_count,
        "regular_region_count": l24_processed,
        "critical_candidate_count": l27_all.len(),
        "regular_regions_closed_by_l24": l24_closed.len(),
        "critical_regions_closed_by_l27_partition": l27_closed.len(),
        "critical_regions_closed_by_l28_root_location": l28_closed.len(),
        "critical_regions_closed_by_l29_sharp_root_location": l29_closed.len(),
        "final_critical_closed_count": final_critical.len(),
        "residual_chain_total_closed_count": residual_chain_total_closed_count,
        "expected_residual_chain_total_count": expected_residual_chain_total_count,
        "closure_counts": closure_counts,
        "ownership_duplicate_count": ownership_duplicate_count,
        "source_filter_mismatch_count": source_filter_mismatch_count,
        "source_filter_mismatches": source_mismatches,
        "candidate_count_drift_count": candidate_count_drift_count,
        "candidate_count_drifts": count_drifts,
        "budget_pass": total <= cap,
        "length_accounting_policy": "L24 is the authoritative repaired regular-slice length total. L27/L28/L29 are critical-candidate filters and may carry a lower pre-repair source total, but may not exceed L24 or drift the cap.",
        "source_accepted_length_upper": total,
        "total_validated_length_upper": total,
        "exact_length_cap": cap,
        "margin_to_cap": margin,
        "first_failed_condition": first_failed,
        "proof_obligation": "Integration certificate over L24, L27, L28, and L29 residual closures; integrated length is inherited from the repaired L24 regular-slice certificate.",
        "claim_ceiling": "Local n=14 hard-cell residual-chain integration certificate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law."
    });

    let result_path = cfg
        .out_dir
        .join(format!("{}_RESULTS.json", cfg.experiment_id));
    let report_path = cfg.out_dir.join(format!("{}_REPORT.md", cfg.experiment_id));
    let sha_path = cfg
        .out_dir
        .join(format!("{}_RESULTS.sha256", cfg.experiment_id));
    let result_text = serde_json::to_string_pretty(&result)
        .map_err(|err| format!("failed to serialize result JSON: {err}"))?;
    fs::write(&result_path, result_text)
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
            "regular_regions_closed_by_l24": result["regular_regions_closed_by_l24"],
            "final_critical_closed_count": result["final_critical_closed_count"],
            "residual_chain_total_closed_count": result["residual_chain_total_closed_count"],
            "ownership_duplicate_count": result["ownership_duplicate_count"],
            "source_filter_mismatch_count": result["source_filter_mismatch_count"],
            "candidate_count_drift_count": result["candidate_count_drift_count"],
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
