//! EHP #114 n=14 CELL-00-00 repaired slab-source overlay.
//!
//! This binary materializes the L33F replacement accounting as a source-like
//! artifact. It does not recompute geometry; it binds an existing slab source
//! and a verified high-slope repair artifact into explicit machine-readable
//! branch replacements so downstream diagnostics do not rely on prose.

use ehp_n3_poc::ehp114_n14_cell::{ensure_source_subcell_matches, CellSpec};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const DEFAULT_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-CELL-00-00-REPAIRED-SLAB-SOURCE-OVERLAY-20260506-01";
const DEFAULT_SOURCE: &str = "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01_RESULTS.json";
const DEFAULT_REPAIR: &str = "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-00-HIGH-SLOPE-SLAB-REPAIR-20260506-01/EXP-MATH-EHP114-N14-CELL-00-00-HIGH-SLOPE-SLAB-REPAIR-20260506-01_RESULTS.json";
const DEFAULT_OUTDIR: &str =
    "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-00-REPAIRED-SLAB-SOURCE-OVERLAY-20260506-01";
const EXACT_LENGTH_CAP: f64 = 20.672796062619668;

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    repair: PathBuf,
    out_dir: PathBuf,
    experiment_id: String,
    cell: CellSpec,
}

#[derive(Clone, Debug)]
struct Replacement {
    ownership_key: String,
    original_length_upper: f64,
    repaired_length_upper: f64,
    repaired_segment_count: usize,
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from(DEFAULT_SOURCE),
        repair: PathBuf::from(DEFAULT_REPAIR),
        out_dir: PathBuf::from(DEFAULT_OUTDIR),
        experiment_id: DEFAULT_EXPERIMENT_ID.to_string(),
        cell: CellSpec::new(0, 0)?,
    };
    let mut sub_i = cfg.cell.sub_i;
    let mut sub_j = cfg.cell.sub_j;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--source" => {
                cfg.source = PathBuf::from(next_string(&args, i, "--source")?);
                i += 2;
            }
            "--repair" => {
                cfg.repair = PathBuf::from(next_string(&args, i, "--repair")?);
                i += 2;
            }
            "--outdir" => {
                cfg.out_dir = PathBuf::from(next_string(&args, i, "--outdir")?);
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

fn next_string(args: &[String], i: usize, flag: &str) -> Result<String, String> {
    if i + 1 >= args.len() {
        return Err(format!("{flag} requires a value"));
    }
    Ok(args[i + 1].clone())
}

fn unix_timestamp_string() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs()
        .to_string()
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

fn f64_field(value: &Value, key: &str) -> Result<f64, String> {
    value
        .get(key)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("missing numeric field {key}"))
}

fn usize_field(value: &Value, key: &str) -> usize {
    value.get(key).and_then(Value::as_u64).unwrap_or(0) as usize
}

fn str_field(value: &Value, key: &str) -> String {
    value
        .get(key)
        .and_then(Value::as_str)
        .unwrap_or("unknown")
        .to_string()
}

fn source_branches_by_key(source: &Value) -> Result<BTreeMap<String, Value>, String> {
    let branches = source
        .get("branches")
        .and_then(Value::as_array)
        .ok_or_else(|| "source missing branches array".to_string())?;
    let mut out = BTreeMap::new();
    for branch in branches {
        let key = str_field(branch, "ownership_key");
        if out.insert(key.clone(), branch.clone()).is_some() {
            return Err(format!("duplicate source branch ownership key {key}"));
        }
    }
    Ok(out)
}

fn parse_replacements(repair: &Value) -> Result<Vec<Replacement>, String> {
    let repairs = repair
        .get("repairs")
        .and_then(Value::as_array)
        .ok_or_else(|| "repair artifact missing repairs array".to_string())?;
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for item in repairs {
        if str_field(item, "status") != "TARGET_REPAIR_CERTIFIED" {
            return Err(format!(
                "repair target {} is not certified",
                item["source_branch"]["ownership_key"]
            ));
        }
        if usize_field(item, "unresolved_group_count") != 0 {
            return Err(format!(
                "repair target {} still has unresolved groups",
                item["source_branch"]["ownership_key"]
            ));
        }
        let key = str_field(&item["source_branch"], "ownership_key");
        if !seen.insert(key.clone()) {
            return Err(format!("duplicate repair replacement key {key}"));
        }
        out.push(Replacement {
            ownership_key: key,
            original_length_upper: f64_field(item, "original_length_upper")?,
            repaired_length_upper: f64_field(item, "repaired_length_upper")?,
            repaired_segment_count: usize_field(item, "repaired_segment_count"),
        });
    }
    Ok(out)
}

fn updated_branches(source: &Value, replacements: &[Replacement]) -> Result<Vec<Value>, String> {
    let replacement_by_key: BTreeMap<String, &Replacement> = replacements
        .iter()
        .map(|replacement| (replacement.ownership_key.clone(), replacement))
        .collect();
    let branches = source
        .get("branches")
        .and_then(Value::as_array)
        .ok_or_else(|| "source missing branches array".to_string())?;
    let mut out = Vec::with_capacity(branches.len());
    for branch in branches {
        let key = str_field(branch, "ownership_key");
        let mut next = branch.clone();
        if let Some(replacement) = replacement_by_key.get(&key) {
            let original = f64_field(branch, "length_upper")?;
            if (original - replacement.original_length_upper).abs() > 1.0e-9 {
                return Err(format!(
                    "replacement original length mismatch for {key}: source {original}, repair {}",
                    replacement.original_length_upper
                ));
            }
            next["length_upper"] = json!(replacement.repaired_length_upper);
            next["overlay_repaired"] = json!(true);
            next["overlay_original_length_upper"] = json!(replacement.original_length_upper);
            next["overlay_repaired_length_upper"] = json!(replacement.repaired_length_upper);
            next["overlay_repaired_segment_count"] = json!(replacement.repaired_segment_count);
        }
        out.push(next);
    }
    Ok(out)
}

fn write_report(result: &Value, path: &Path) -> Result<(), String> {
    let report_cell = result["cell_tag"].as_str().unwrap_or("unknown");
    let report = format!(
        "# EHP114 n=14 {report_cell} Repaired Slab-Source Overlay\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Source total upper: `{}`\n\
- Repaired total upper: `{}`\n\
- Exact cap: `{}`\n\
- Repaired margin to cap: `{}`\n\
- Replacement count: `{}`\n\
- Replacement duplicate count: `{}`\n\
- Remaining unresolved branches: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
This artifact materializes L33F as explicit source replacement accounting. It \
does not close the remaining unresolved branch chain and does not promote a \
global n=14 claim. It is a safer input contract for downstream local pipeline \
diagnostics than prose replacement arithmetic.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        result["status"].as_str().unwrap_or("UNKNOWN"),
        result["source_total_validated_length_upper"],
        result["total_validated_length_upper"],
        result["length_budget"]["exact_length_cap"],
        result["margin_to_cap"],
        result["replacement_count"],
        result["replacement_duplicate_count"],
        result["unresolved_branch_count"],
        result["first_failed_condition"]
            .as_str()
            .unwrap_or("unknown"),
        result["claim_ceiling"]
            .as_str()
            .unwrap_or("local overlay only"),
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

    let source = read_json(&cfg.source)?;
    let repair = read_json(&cfg.repair)?;
    let source_contract = ensure_source_subcell_matches(&source, cfg.cell)?;
    let repair_contract = ensure_source_subcell_matches(&repair, cfg.cell)?;
    let replacements = parse_replacements(&repair)?;
    let source_by_key = source_branches_by_key(&source)?;
    let mut missing = Vec::new();
    for replacement in &replacements {
        if !source_by_key.contains_key(&replacement.ownership_key) {
            missing.push(replacement.ownership_key.clone());
        }
    }
    if !missing.is_empty() {
        return Err(format!(
            "repair replacements missing source keys: {}",
            missing.join(",")
        ));
    }

    let updated_branches = updated_branches(&source, &replacements)?;
    let source_total = f64_field(&source, "total_validated_length_upper")?;
    let source_sum = f64_field(&source, "sum_branch_length_upper")?;
    let original_sum: f64 = replacements
        .iter()
        .map(|replacement| replacement.original_length_upper)
        .sum();
    let repaired_sum: f64 = replacements
        .iter()
        .map(|replacement| replacement.repaired_length_upper)
        .sum();
    let length_delta = repaired_sum - original_sum;
    let repaired_sum_branch_length_upper = source_sum + length_delta;
    let repaired_total = source_total + length_delta;
    let margin_to_cap = EXACT_LENGTH_CAP - repaired_total;
    let replacement_duplicate_count = replacements.len()
        - replacements
            .iter()
            .map(|replacement| replacement.ownership_key.clone())
            .collect::<BTreeSet<_>>()
            .len();
    let source_unresolved = usize_field(&source, "unresolved_branch_count");
    let status = if replacement_duplicate_count > 0 {
        "REPAIRED_SLAB_SOURCE_OVERLAY_FAIL_DUPLICATE_REPLACEMENT"
    } else if repaired_total > EXACT_LENGTH_CAP {
        "REPAIRED_SLAB_SOURCE_OVERLAY_FAIL_BUDGET"
    } else {
        "REPAIRED_SLAB_SOURCE_OVERLAY_READY_FOR_DOWNSTREAM_DIAGNOSTIC"
    };
    let first_failed_condition = if replacement_duplicate_count > 0 {
        format!("{replacement_duplicate_count} duplicate replacement keys")
    } else if repaired_total > EXACT_LENGTH_CAP {
        format!(
            "repaired total {} exceeds cap {} by {}",
            repaired_total,
            EXACT_LENGTH_CAP,
            repaired_total - EXACT_LENGTH_CAP
        )
    } else {
        "none".to_string()
    };

    let replacement_json: Vec<Value> = replacements
        .iter()
        .map(|replacement| {
            json!({
                "ownership_key": replacement.ownership_key,
                "original_length_upper": replacement.original_length_upper,
                "repaired_length_upper": replacement.repaired_length_upper,
                "length_delta": replacement.repaired_length_upper - replacement.original_length_upper,
                "repaired_segment_count": replacement.repaired_segment_count,
            })
        })
        .collect();
    let cell_tag = cfg.cell.tag();
    let what_this_rules_out = format!(
        "The largest high-slope accepted branch lengths no longer force the {cell_tag} source above cap once the repair artifact is explicitly applied."
    );
    let claim_ceiling = format!(
        "{cell_tag} repaired slab-source overlay only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate."
    );

    let result = json!({
        "experiment_id": cfg.experiment_id,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "status": status,
        "source_results_path": cfg.source,
        "repair_results_path": cfg.repair,
        "source_subcell_contract": source_contract,
        "repair_subcell_contract": repair_contract,
        "source_experiment_id": source["experiment_id"],
        "repair_experiment_id": repair["experiment_id"],
        "subcell": cfg.cell,
        "cell_tag": cell_tag,
        "length_budget": {
            "exact_length_cap": EXACT_LENGTH_CAP,
            "source_total_validated_length_upper": source_total,
            "source_sum_branch_length_upper": source_sum,
            "replacement_original_length_sum": original_sum,
            "replacement_repaired_length_sum": repaired_sum,
            "replacement_length_delta": length_delta,
            "repaired_sum_branch_length_upper": repaired_sum_branch_length_upper,
        },
        "source_total_validated_length_upper": source_total,
        "sum_branch_length_upper": repaired_sum_branch_length_upper,
        "endpoint_overlap_tax": source["endpoint_overlap_tax"],
        "total_validated_length_upper": repaired_total,
        "margin_to_cap": margin_to_cap,
        "replacement_count": replacements.len(),
        "replacement_duplicate_count": replacement_duplicate_count,
        "branch_repair_overrides": replacement_json,
        "slab_branch_count": usize_field(&source, "slab_branch_count"),
        "candidate_cell_count": usize_field(&source, "candidate_cell_count"),
        "excluded_cell_count": usize_field(&source, "excluded_cell_count"),
        "unresolved_branch_count": source_unresolved,
        "ownership_duplicate_count": usize_field(&source, "ownership_duplicate_count"),
        "branches": updated_branches,
        "unresolved_branches": source["unresolved_branches"].clone(),
        "first_failed_condition": first_failed_condition,
        "what_this_rules_out": what_this_rules_out,
        "what_this_does_not_rule_out": "This overlay does not certify the remaining unresolved branch chain, does not cover other cells, and does not prove global n=14 coverage.",
        "next_dependency": "Run the per-cell source-chain smoke and downstream L21 branch-isolation/collar pipeline against this overlay, with source/subcell matching enforced.",
        "claim_ceiling": claim_ceiling
    });

    fs::write(
        &result_path,
        serde_json::to_string_pretty(&result).unwrap() + "\n",
    )
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
            "replacement_count": result["replacement_count"],
            "source_total_validated_length_upper": result["source_total_validated_length_upper"],
            "total_validated_length_upper": result["total_validated_length_upper"],
            "exact_length_cap": result["length_budget"]["exact_length_cap"],
            "margin_to_cap": result["margin_to_cap"],
            "unresolved_branch_count": result["unresolved_branch_count"],
            "first_failed_condition": result["first_failed_condition"],
            "sha256": digest
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
