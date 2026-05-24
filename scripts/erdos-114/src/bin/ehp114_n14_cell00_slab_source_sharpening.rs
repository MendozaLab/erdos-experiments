//! EHP #114 n=14 CELL-00-00 slab-source sharpening diagnostic.
//!
//! This binary does not certify a new local cell. It audits the first
//! non-hard-cell slab source and identifies whether the over-budget bound is
//! concentrated in a small number of high-slope/low-denominator branches or is
//! a global source-quality problem.

use ehp_n3_poc::ehp114_n14_cell::{ensure_source_subcell_matches, CellSpec};
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const DEFAULT_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-CELL-00-00-SLAB-SOURCE-SHARPENING-20260506-01";
const DEFAULT_SOURCE: &str = "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01_RESULTS.json";
const DEFAULT_BASELINE: &str = "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03_RESULTS.json";
const DEFAULT_OUTDIR: &str =
    "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-00-SLAB-SOURCE-SHARPENING-20260506-01";

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    baseline: PathBuf,
    out_dir: PathBuf,
    experiment_id: String,
    cell: CellSpec,
}

#[derive(Clone, Debug)]
struct BranchMetric {
    ix: usize,
    group_index: usize,
    ownership_key: String,
    length_upper: f64,
    slope_abs_upper: f64,
    denominator_abs_lower: f64,
    numerator_abs_upper: f64,
    y_height: f64,
    y_cell_span: usize,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
}

#[derive(Clone, Debug)]
struct UnresolvedMetric {
    ix: usize,
    group_index: usize,
    ownership_key: String,
    reason: String,
    fx_abs_lower: f64,
    fy_abs_lower: f64,
    fx_sign: String,
    fy_sign: String,
    candidate_cell_count: usize,
    y_height: f64,
    y_cell_span: usize,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
}

#[derive(Clone, Debug, Default, Serialize)]
struct SlabAggregate {
    ix: usize,
    branch_count: usize,
    unresolved_count: usize,
    length_sum: f64,
    max_branch_length: f64,
    max_slope_abs_upper: f64,
    min_denominator_abs_lower: f64,
    unresolved_candidate_cell_count: usize,
}

#[derive(Clone, Debug, Serialize)]
struct TopBranch {
    ix: usize,
    group_index: usize,
    ownership_key: String,
    length_upper: f64,
    slope_abs_upper: f64,
    denominator_abs_lower: f64,
    numerator_abs_upper: f64,
    y_height: f64,
    y_cell_span: usize,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
}

#[derive(Clone, Debug, Serialize)]
struct TopUnresolved {
    ix: usize,
    group_index: usize,
    ownership_key: String,
    reason: String,
    fx_abs_lower: f64,
    fy_abs_lower: f64,
    fx_sign: String,
    fy_sign: String,
    candidate_cell_count: usize,
    y_height: f64,
    y_cell_span: usize,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
}

#[derive(Clone, Debug, Serialize)]
struct SourceStats {
    experiment_id: String,
    status: String,
    cell_tag: String,
    total_validated_length_upper: f64,
    exact_length_cap: f64,
    margin_to_cap: f64,
    branch_count: usize,
    unresolved_branch_count: usize,
    ownership_duplicate_count: usize,
    candidate_cell_count: usize,
    excluded_cell_count: usize,
    sum_branch_length_upper: f64,
    length_quantiles: BTreeMap<String, f64>,
    slope_quantiles: BTreeMap<String, f64>,
    denominator_quantiles: BTreeMap<String, f64>,
    y_height_quantiles: BTreeMap<String, f64>,
    high_slope_counts: BTreeMap<String, usize>,
    low_denominator_counts: BTreeMap<String, usize>,
    top_branch_length_sums: BTreeMap<String, f64>,
    top_branch_length_fraction_of_total: BTreeMap<String, f64>,
    top_branches_by_length: Vec<TopBranch>,
    top_slabs_by_length: Vec<SlabAggregate>,
    top_slabs_by_unresolved_count: Vec<SlabAggregate>,
    unresolved_reason_counts: BTreeMap<String, usize>,
    top_unresolved_by_candidate_cells: Vec<TopUnresolved>,
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from(DEFAULT_SOURCE),
        baseline: PathBuf::from(DEFAULT_BASELINE),
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
            "--baseline" => {
                cfg.baseline = PathBuf::from(next_string(&args, i, "--baseline")?);
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

fn str_field(value: &Value, key: &str, fallback: &str) -> String {
    value
        .get(key)
        .and_then(Value::as_str)
        .unwrap_or(fallback)
        .to_string()
}

fn usize_field(value: &Value, key: &str) -> usize {
    value.get(key).and_then(Value::as_u64).unwrap_or(0) as usize
}

fn f64_field(value: &Value, key: &str) -> Result<f64, String> {
    value
        .get(key)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("missing numeric field {key}"))
}

fn exact_cap(value: &Value) -> Result<f64, String> {
    value
        .get("length_budget")
        .and_then(|budget| budget.get("exact_length_cap"))
        .and_then(Value::as_f64)
        .ok_or_else(|| "missing length_budget.exact_length_cap".to_string())
}

fn pair_field(value: &Value, key: &str) -> Result<[f64; 2], String> {
    let arr = value
        .get(key)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("missing pair field {key}"))?;
    if arr.len() != 2 {
        return Err(format!("{key} does not have length 2"));
    }
    let a = arr[0]
        .as_f64()
        .ok_or_else(|| format!("{key}[0] is not numeric"))?;
    let b = arr[1]
        .as_f64()
        .ok_or_else(|| format!("{key}[1] is not numeric"))?;
    Ok([a, b])
}

fn y_cell_span(value: &Value) -> usize {
    let Some(arr) = value.get("y_cells").and_then(Value::as_array) else {
        return 0;
    };
    if arr.len() != 2 {
        return 0;
    }
    let a = arr[0].as_u64().unwrap_or(0) as usize;
    let b = arr[1].as_u64().unwrap_or(0) as usize;
    b.saturating_sub(a).saturating_add(1)
}

fn branch_metric(value: &Value) -> Result<BranchMetric, String> {
    let y = pair_field(value, "y_interval")?;
    Ok(BranchMetric {
        ix: usize_field(value, "ix"),
        group_index: usize_field(value, "group_index"),
        ownership_key: str_field(value, "ownership_key", "unknown"),
        length_upper: f64_field(value, "length_upper")?,
        slope_abs_upper: f64_field(value, "slope_abs_upper")?,
        denominator_abs_lower: f64_field(value, "denominator_abs_lower")?,
        numerator_abs_upper: f64_field(value, "numerator_abs_upper")?,
        y_height: (y[1] - y[0]).abs(),
        y_cell_span: y_cell_span(value),
        x_interval: pair_field(value, "x_interval")?,
        y_interval: y,
    })
}

fn unresolved_metric(value: &Value) -> Result<UnresolvedMetric, String> {
    let y = pair_field(value, "y_interval")?;
    Ok(UnresolvedMetric {
        ix: usize_field(value, "ix"),
        group_index: usize_field(value, "group_index"),
        ownership_key: str_field(value, "ownership_key", "unknown"),
        reason: str_field(value, "reason", "unknown"),
        fx_abs_lower: f64_field(value, "fx_abs_lower").unwrap_or(0.0),
        fy_abs_lower: f64_field(value, "fy_abs_lower").unwrap_or(0.0),
        fx_sign: str_field(value, "fx_sign", "unknown"),
        fy_sign: str_field(value, "fy_sign", "unknown"),
        candidate_cell_count: usize_field(value, "candidate_cell_count"),
        y_height: (y[1] - y[0]).abs(),
        y_cell_span: y_cell_span(value),
        x_interval: pair_field(value, "x_interval")?,
        y_interval: y,
    })
}

fn quantiles(values: &[f64]) -> BTreeMap<String, f64> {
    let mut map = BTreeMap::new();
    if values.is_empty() {
        return map;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    for (label, q) in [
        ("p50", 0.50),
        ("p75", 0.75),
        ("p90", 0.90),
        ("p95", 0.95),
        ("p99", 0.99),
        ("max", 1.00),
    ] {
        let idx = ((sorted.len() - 1) as f64 * q).round() as usize;
        map.insert(label.to_string(), sorted[idx]);
    }
    map
}

fn count_above(values: &[f64], threshold: f64) -> usize {
    values.iter().filter(|v| **v >= threshold).count()
}

fn count_below(values: &[f64], threshold: f64) -> usize {
    values.iter().filter(|v| **v <= threshold).count()
}

fn top_sum(sorted_desc: &[BranchMetric], n: usize) -> f64 {
    sorted_desc
        .iter()
        .take(n.min(sorted_desc.len()))
        .map(|b| b.length_upper)
        .sum()
}

fn top_branch(branch: &BranchMetric) -> TopBranch {
    TopBranch {
        ix: branch.ix,
        group_index: branch.group_index,
        ownership_key: branch.ownership_key.clone(),
        length_upper: branch.length_upper,
        slope_abs_upper: branch.slope_abs_upper,
        denominator_abs_lower: branch.denominator_abs_lower,
        numerator_abs_upper: branch.numerator_abs_upper,
        y_height: branch.y_height,
        y_cell_span: branch.y_cell_span,
        x_interval: branch.x_interval,
        y_interval: branch.y_interval,
    }
}

fn top_unresolved(item: &UnresolvedMetric) -> TopUnresolved {
    TopUnresolved {
        ix: item.ix,
        group_index: item.group_index,
        ownership_key: item.ownership_key.clone(),
        reason: item.reason.clone(),
        fx_abs_lower: item.fx_abs_lower,
        fy_abs_lower: item.fy_abs_lower,
        fx_sign: item.fx_sign.clone(),
        fy_sign: item.fy_sign.clone(),
        candidate_cell_count: item.candidate_cell_count,
        y_height: item.y_height,
        y_cell_span: item.y_cell_span,
        x_interval: item.x_interval,
        y_interval: item.y_interval,
    }
}

fn compute_source_stats(value: &Value) -> Result<SourceStats, String> {
    let branches_json = value
        .get("branches")
        .and_then(Value::as_array)
        .ok_or_else(|| "source missing branches array".to_string())?;
    let unresolved_json = value
        .get("unresolved_branches")
        .and_then(Value::as_array)
        .ok_or_else(|| "source missing unresolved_branches array".to_string())?;
    let mut branches = Vec::with_capacity(branches_json.len());
    for item in branches_json {
        branches.push(branch_metric(item)?);
    }
    let mut unresolved = Vec::with_capacity(unresolved_json.len());
    for item in unresolved_json {
        unresolved.push(unresolved_metric(item)?);
    }

    let mut by_length = branches.clone();
    by_length.sort_by(|a, b| b.length_upper.total_cmp(&a.length_upper));
    let branch_total: f64 = branches.iter().map(|b| b.length_upper).sum();
    let lengths: Vec<f64> = branches.iter().map(|b| b.length_upper).collect();
    let slopes: Vec<f64> = branches.iter().map(|b| b.slope_abs_upper).collect();
    let denoms: Vec<f64> = branches.iter().map(|b| b.denominator_abs_lower).collect();
    let heights: Vec<f64> = branches.iter().map(|b| b.y_height).collect();
    let mut top_sums = BTreeMap::new();
    let mut top_fractions = BTreeMap::new();
    for n in [1usize, 5, 10, 25, 50, 100] {
        let sum = top_sum(&by_length, n);
        top_sums.insert(format!("top_{n}"), sum);
        top_fractions.insert(
            format!("top_{n}"),
            if branch_total > 0.0 {
                sum / branch_total
            } else {
                0.0
            },
        );
    }

    let mut high_slope_counts = BTreeMap::new();
    for threshold in [100.0, 500.0, 1000.0, 5000.0] {
        high_slope_counts.insert(format!("ge_{threshold}"), count_above(&slopes, threshold));
    }
    let mut low_denominator_counts = BTreeMap::new();
    for threshold in [0.001, 0.01, 0.05, 0.1] {
        low_denominator_counts.insert(format!("le_{threshold}"), count_below(&denoms, threshold));
    }

    let mut slabs = BTreeMap::<usize, SlabAggregate>::new();
    for branch in &branches {
        let entry = slabs.entry(branch.ix).or_insert_with(|| SlabAggregate {
            ix: branch.ix,
            min_denominator_abs_lower: f64::INFINITY,
            ..SlabAggregate::default()
        });
        entry.branch_count += 1;
        entry.length_sum += branch.length_upper;
        entry.max_branch_length = entry.max_branch_length.max(branch.length_upper);
        entry.max_slope_abs_upper = entry.max_slope_abs_upper.max(branch.slope_abs_upper);
        entry.min_denominator_abs_lower = entry
            .min_denominator_abs_lower
            .min(branch.denominator_abs_lower);
    }
    let mut unresolved_reason_counts = BTreeMap::<String, usize>::new();
    for item in &unresolved {
        *unresolved_reason_counts
            .entry(item.reason.clone())
            .or_default() += 1;
        let entry = slabs.entry(item.ix).or_insert_with(|| SlabAggregate {
            ix: item.ix,
            min_denominator_abs_lower: f64::INFINITY,
            ..SlabAggregate::default()
        });
        entry.unresolved_count += 1;
        entry.unresolved_candidate_cell_count += item.candidate_cell_count;
    }
    let mut slab_values: Vec<SlabAggregate> = slabs
        .into_values()
        .map(|mut slab| {
            if !slab.min_denominator_abs_lower.is_finite() {
                slab.min_denominator_abs_lower = 0.0;
            }
            slab
        })
        .collect();
    slab_values.sort_by(|a, b| b.length_sum.total_cmp(&a.length_sum));
    let top_slabs_by_length = slab_values.iter().take(20).cloned().collect();
    slab_values.sort_by(|a, b| {
        b.unresolved_count
            .cmp(&a.unresolved_count)
            .then_with(|| b.length_sum.total_cmp(&a.length_sum))
    });
    let top_slabs_by_unresolved_count = slab_values.iter().take(20).cloned().collect();

    let mut unresolved_by_cells = unresolved.clone();
    unresolved_by_cells.sort_by(|a, b| {
        b.candidate_cell_count
            .cmp(&a.candidate_cell_count)
            .then_with(|| b.y_height.total_cmp(&a.y_height))
    });

    Ok(SourceStats {
        experiment_id: str_field(value, "experiment_id", "unknown"),
        status: str_field(value, "status", "unknown"),
        cell_tag: str_field(value, "cell_tag", "unknown"),
        total_validated_length_upper: f64_field(value, "total_validated_length_upper")?,
        exact_length_cap: exact_cap(value)?,
        margin_to_cap: f64_field(value, "margin_to_cap")?,
        branch_count: branches.len(),
        unresolved_branch_count: unresolved.len(),
        ownership_duplicate_count: usize_field(value, "ownership_duplicate_count"),
        candidate_cell_count: usize_field(value, "candidate_cell_count"),
        excluded_cell_count: usize_field(value, "excluded_cell_count"),
        sum_branch_length_upper: f64_field(value, "sum_branch_length_upper")?,
        length_quantiles: quantiles(&lengths),
        slope_quantiles: quantiles(&slopes),
        denominator_quantiles: quantiles(&denoms),
        y_height_quantiles: quantiles(&heights),
        high_slope_counts,
        low_denominator_counts,
        top_branch_length_sums: top_sums,
        top_branch_length_fraction_of_total: top_fractions,
        top_branches_by_length: by_length.iter().take(20).map(top_branch).collect(),
        top_slabs_by_length,
        top_slabs_by_unresolved_count,
        unresolved_reason_counts,
        top_unresolved_by_candidate_cells: unresolved_by_cells
            .iter()
            .take(20)
            .map(top_unresolved)
            .collect(),
    })
}

fn write_report(result: &Value, path: &Path) -> Result<(), String> {
    let report_cell = result["source_stats"]["cell_tag"]
        .as_str()
        .unwrap_or("unknown");
    let report = format!(
        "# EHP114 n=14 {report_cell} Slab Source Sharpening Diagnostic\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Source cell: `{}`\n\
- Source total upper: `{}`\n\
- Exact cap: `{}`\n\
- Margin to cap: `{}`\n\
- Unresolved branches: `{}`\n\
- Worst branch length: `{}`\n\
- Worst branch slope: `{}`\n\
- Top-10 branch length sum: `{}`\n\
- Source-chain recommendation: `{}`\n\n\
## Interpretation\n\n\
This diagnostic audits the first non-hard-cell slab source before downstream \
local-pipeline promotion. The key question is whether the over-budget bound is \
localized to a few high-slope/low-denominator slabs or spread globally. It does \
not certify a cell and does not promote length.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        result["status"].as_str().unwrap_or("UNKNOWN"),
        report_cell,
        result["source_stats"]["total_validated_length_upper"],
        result["source_stats"]["exact_length_cap"],
        result["source_stats"]["margin_to_cap"],
        result["source_stats"]["unresolved_branch_count"],
        result["source_stats"]["top_branches_by_length"][0]["length_upper"],
        result["source_stats"]["top_branches_by_length"][0]["slope_abs_upper"],
        result["source_stats"]["top_branch_length_sums"]["top_10"],
        result["recommended_next_route"]
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
    let baseline = read_json(&cfg.baseline)?;
    let source_contract = ensure_source_subcell_matches(&source, cfg.cell)?;
    let source_stats = compute_source_stats(&source)?;
    let baseline_stats = compute_source_stats(&baseline)?;
    let excess_over_cap =
        (source_stats.total_validated_length_upper - source_stats.exact_length_cap).max(0.0);
    let top_10_sum = source_stats
        .top_branch_length_sums
        .get("top_10")
        .copied()
        .unwrap_or(0.0);
    let worst_branch = source_stats.top_branches_by_length.first().cloned();
    let low_denom_count = *source_stats
        .low_denominator_counts
        .get("le_0.01")
        .unwrap_or(&0);
    let high_slope_count = *source_stats.high_slope_counts.get("ge_1000").unwrap_or(&0);
    let over_budget = source_stats.total_validated_length_upper > source_stats.exact_length_cap;
    let unresolved = source_stats.unresolved_branch_count > 0;
    let concentrated_enough_to_repair = top_10_sum >= excess_over_cap;
    let recommended_next_route = if over_budget && concentrated_enough_to_repair {
        "targeted_high_slope_slab_repair_before_downstream_pipeline"
    } else if over_budget && high_slope_count > 0 {
        "affine_or_rotated_chart_source_bound_for_high_slope_slabs"
    } else if unresolved {
        "root_isolation_sharpening_before_length_promotion"
    } else {
        "inspect_cap_or_decomposition_before_downstream_pipeline"
    };
    let status = if over_budget {
        "SLAB_SOURCE_SHARPENING_DIAGNOSTIC_FAIL_OVER_BUDGET"
    } else if unresolved {
        "SLAB_SOURCE_SHARPENING_DIAGNOSTIC_FAIL_UNRESOLVED"
    } else {
        "SLAB_SOURCE_SHARPENING_DIAGNOSTIC_READY_FOR_DOWNSTREAM"
    };
    let first_failed_condition = if over_budget {
        format!(
            "source total upper {} exceeds cap {} by {}",
            source_stats.total_validated_length_upper,
            source_stats.exact_length_cap,
            excess_over_cap
        )
    } else if unresolved {
        format!(
            "source has {} unresolved branches",
            source_stats.unresolved_branch_count
        )
    } else {
        "none".to_string()
    };

    let cell_tag = cfg.cell.tag();
    let what_this_rules_out = format!(
        "The {cell_tag} source absence blocker is gone, but the current z32 slab source cannot be used as a theorem-packet input if its own bound is over cap."
    );
    let what_this_does_not_rule_out = format!(
        "This does not rule out {cell_tag} mathematically. It only identifies source-bound sharpness, especially high-slope or low-denominator branch contributions, as the next repair target."
    );
    let claim_ceiling = format!(
        "{cell_tag} slab source-sharpening diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate."
    );

    let result = json!({
        "experiment_id": cfg.experiment_id,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "status": status,
        "source_results_path": cfg.source,
        "baseline_results_path": cfg.baseline,
        "source_subcell_contract": source_contract,
        "source_stats": source_stats,
        "baseline_stats": baseline_stats,
        "comparison_to_hard_cell": {
            "total_length_delta": source_stats.total_validated_length_upper - baseline_stats.total_validated_length_upper,
            "margin_delta": source_stats.margin_to_cap - baseline_stats.margin_to_cap,
            "branch_count_delta": source_stats.branch_count as i64 - baseline_stats.branch_count as i64,
            "unresolved_branch_count_delta": source_stats.unresolved_branch_count as i64 - baseline_stats.unresolved_branch_count as i64,
            "worst_branch_length_delta": source_stats.top_branches_by_length.first().map(|b| b.length_upper).unwrap_or(0.0)
                - baseline_stats.top_branches_by_length.first().map(|b| b.length_upper).unwrap_or(0.0),
            "top_10_length_sum_delta": source_stats.top_branch_length_sums.get("top_10").copied().unwrap_or(0.0)
                - baseline_stats.top_branch_length_sums.get("top_10").copied().unwrap_or(0.0),
        },
        "overage_analysis": {
            "over_budget": over_budget,
            "excess_over_cap": excess_over_cap,
            "top_10_branch_sum": top_10_sum,
            "top_10_sum_exceeds_excess_over_cap": concentrated_enough_to_repair,
            "low_denominator_branch_count_le_0_01": low_denom_count,
            "high_slope_branch_count_ge_1000": high_slope_count,
            "worst_branch": worst_branch,
        },
        "recommended_next_route": recommended_next_route,
        "first_failed_condition": first_failed_condition,
        "what_this_rules_out": what_this_rules_out,
        "what_this_does_not_rule_out": what_this_does_not_rule_out,
        "next_dependency": "Run a targeted high-slope slab repair on the top length-contributing slabs before downstream L21-L32 promotion.",
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
            "source_cell": result["source_stats"]["cell_tag"],
            "source_total": result["source_stats"]["total_validated_length_upper"],
            "cap": result["source_stats"]["exact_length_cap"],
            "margin_to_cap": result["source_stats"]["margin_to_cap"],
            "unresolved_branch_count": result["source_stats"]["unresolved_branch_count"],
            "top_10_branch_sum": result["overage_analysis"]["top_10_branch_sum"],
            "recommended_next_route": result["recommended_next_route"],
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
