//! EHP #114 n=14 global atlas certificate gate.
//!
//! This is an assembler/auditor over theorem-shaped local cell packets. It
//! does not recompute geometry. It enumerates the 8x8 root-affine atlas cells,
//! verifies source hashes, imports any local theorem packets it can find, and
//! fails explicitly if a cell is missing, duplicated, drifting, or over budget.
//!
//! Directory discovery intentionally selects the latest theorem-shaped packet
//! per cell. Older packet versions remain immutable evidence, but they are not
//! duplicate coverage in the current atlas layer.

use serde::Serialize;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const DEFAULT_EXPERIMENT_ID: &str = "EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-02";
const ROOT_AFFINE_SOURCE_ID: &str = "EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01";
const LOCAL_PACKET_STATUS: &str = "LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF";
const DEGREE: usize = 14;
const EPS: f64 = 0.1;
const SUBDIVISION: usize = 8;
const EXPECTED_CELL_COUNT: usize = SUBDIVISION * SUBDIVISION;
const EXACT_LENGTH_CAP: f64 = 20.672796062619668;

#[derive(Clone, Debug)]
struct Config {
    experiment_id: String,
    root_affine_source: PathBuf,
    cell_packet_dir: PathBuf,
    explicit_cell_packets: Vec<PathBuf>,
    out_dir: PathBuf,
}

#[derive(Clone, Debug, Serialize)]
struct SourceSha {
    role: String,
    experiment_id: String,
    path: String,
    sha_status: String,
    expected_sha256: String,
    actual_sha256: String,
}

#[derive(Clone, Debug, Serialize)]
struct CellSummary {
    subcell: [usize; 2],
    experiment_id: String,
    status: String,
    path: String,
    sha_status: String,
    total_validated_length_upper: f64,
    exact_length_cap: f64,
    margin_to_cap: f64,
    residual_chain_total_closed_count: usize,
    expected_residual_chain_total_count: usize,
    residual_universe_accounting_status: Option<String>,
    source_backed_leaf_universe_count: Option<usize>,
    upstream_closed_count: Option<usize>,
    downstream_closed_count: Option<usize>,
    coverage_equation_pass: Option<bool>,
    source_sha_fail_count: usize,
    ownership_duplicate_count: usize,
    source_filter_mismatch_count: usize,
    candidate_count_drift_count: usize,
    accepted: bool,
    rejection_reason: String,
}

fn parse_args() -> Result<Config, String> {
    let mut cfg = Config {
        experiment_id: DEFAULT_EXPERIMENT_ID.to_string(),
        root_affine_source: PathBuf::from("../../Erdos114/exact_length_lift/EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01_RESULTS.json"),
        cell_packet_dir: PathBuf::from("../../Erdos114/validated_length"),
        explicit_cell_packets: Vec::new(),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-02"),
    };

    let args: Vec<String> = env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--experiment-id" => {
                cfg.experiment_id = next_string(&args, i, "--experiment-id")?;
                i += 2;
            }
            "--root-affine-source" => {
                cfg.root_affine_source = next_path(&args, i, "--root-affine-source")?;
                i += 2;
            }
            "--cell-packet-dir" => {
                cfg.cell_packet_dir = next_path(&args, i, "--cell-packet-dir")?;
                i += 2;
            }
            "--cell-packet" => {
                cfg.explicit_cell_packets
                    .push(next_path(&args, i, "--cell-packet")?);
                i += 2;
            }
            "--outdir" => {
                cfg.out_dir = next_path(&args, i, "--outdir")?;
                i += 2;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
    }
    Ok(cfg)
}

fn next_string(args: &[String], i: usize, flag: &str) -> Result<String, String> {
    if i + 1 >= args.len() {
        return Err(format!("{flag} requires a value"));
    }
    Ok(args[i + 1].clone())
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
        role: role.to_string(),
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

fn value_usize(v: &Value, key: &str) -> usize {
    v.get(key).and_then(Value::as_u64).unwrap_or(0) as usize
}

fn optional_usize(v: &Value, key: &str) -> Option<usize> {
    v.get(key).and_then(Value::as_u64).map(|n| n as usize)
}

fn optional_bool(v: &Value, key: &str) -> Option<bool> {
    v.get(key).and_then(Value::as_bool)
}

fn optional_string(v: &Value, key: &str) -> Option<String> {
    v.get(key).and_then(Value::as_str).map(str::to_string)
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

fn expected_cells() -> BTreeSet<(usize, usize)> {
    let mut cells = BTreeSet::new();
    for i in 0..SUBDIVISION {
        for j in 0..SUBDIVISION {
            cells.insert((i, j));
        }
    }
    cells
}

fn row_subcell(row: &Value) -> Result<(usize, usize), String> {
    let sub_i = expect_usize(row, "sub_i")?;
    let sub_j = expect_usize(row, "sub_j")?;
    Ok((sub_i, sub_j))
}

fn validate_root_affine_source(
    source: &Value,
) -> Result<(Vec<String>, BTreeSet<(usize, usize)>), String> {
    let mut failures = Vec::new();
    if expect_str(source, "experiment_id")? != ROOT_AFFINE_SOURCE_ID {
        failures.push("root-affine source experiment_id mismatch".to_string());
    }
    if expect_str(source, "status")? != "BUDGET_ONLY_NOT_EXACT_LENGTH_CERTIFICATE" {
        failures.push("root-affine source has unexpected status".to_string());
    }
    let rows = source
        .get("rows")
        .and_then(Value::as_array)
        .ok_or_else(|| "root-affine source missing rows array".to_string())?;
    if rows.len() != EXPECTED_CELL_COUNT {
        failures.push(format!(
            "root-affine source row count is {}, expected {}",
            rows.len(),
            EXPECTED_CELL_COUNT
        ));
    }
    if value_usize(source, "subcell_count") != EXPECTED_CELL_COUNT {
        failures.push("root-affine source subcell_count is not 64".to_string());
    }

    let mut seen = BTreeSet::new();
    let mut duplicate_count = 0usize;
    let mut out_of_range_count = 0usize;
    for row in rows {
        let (i, j) = row_subcell(row)?;
        if i >= SUBDIVISION || j >= SUBDIVISION {
            out_of_range_count += 1;
        }
        if !seen.insert((i, j)) {
            duplicate_count += 1;
        }
    }
    if duplicate_count > 0 {
        failures.push(format!(
            "root-affine source has {duplicate_count} duplicate cell rows"
        ));
    }
    if out_of_range_count > 0 {
        failures.push(format!(
            "root-affine source has {out_of_range_count} out-of-range cell rows"
        ));
    }
    let expected = expected_cells();
    let missing: Vec<String> = expected
        .difference(&seen)
        .map(|(i, j)| format!("({i},{j})"))
        .collect();
    if !missing.is_empty() {
        failures.push(format!(
            "root-affine source is missing cells: {}",
            missing.join(", ")
        ));
    }
    Ok((failures, seen))
}

fn discover_cell_packets(cfg: &Config) -> Result<Vec<PathBuf>, String> {
    if !cfg.explicit_cell_packets.is_empty() {
        return Ok(cfg.explicit_cell_packets.clone());
    }
    let mut out = Vec::new();
    let entries = fs::read_dir(&cfg.cell_packet_dir)
        .map_err(|err| format!("failed to read {}: {err}", cfg.cell_packet_dir.display()))?;
    for entry in entries {
        let entry = entry.map_err(|err| format!("failed to read directory entry: {err}"))?;
        let path = entry.path();
        if !path.is_dir() {
            continue;
        }
        let name = path
            .file_name()
            .and_then(|v| v.to_str())
            .unwrap_or_default()
            .to_string();
        let is_hard_cell_packet =
            name.starts_with("EXP-MATH-EHP114-N14-LOCAL-HARD-CELL-CERTIFICATE-PACKET");
        let is_parameterized_cell_packet =
            name.starts_with("EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-");
        if !is_hard_cell_packet && !is_parameterized_cell_packet {
            continue;
        }
        let result = path.join(format!("{name}_RESULTS.json"));
        if result.exists() {
            out.push(result);
        }
    }
    out.sort();
    Ok(out)
}

fn parse_subcell(v: &Value) -> Result<[usize; 2], String> {
    let subcell = v
        .get("subcell")
        .ok_or_else(|| "local packet missing subcell field".to_string())?;
    if let Some(arr) = subcell.as_array() {
        if arr.len() != 2 {
            return Err("local packet subcell array must have length 2".to_string());
        }
        let i = arr[0]
            .as_u64()
            .ok_or_else(|| "local packet subcell[0] is not an integer".to_string())?
            as usize;
        let j = arr[1]
            .as_u64()
            .ok_or_else(|| "local packet subcell[1] is not an integer".to_string())?
            as usize;
        return Ok([i, j]);
    }
    if let Some(obj) = subcell.as_object() {
        let i = obj
            .get("sub_i")
            .and_then(Value::as_u64)
            .ok_or_else(|| "local packet subcell.sub_i is not an integer".to_string())?
            as usize;
        let j = obj
            .get("sub_j")
            .and_then(Value::as_u64)
            .ok_or_else(|| "local packet subcell.sub_j is not an integer".to_string())?
            as usize;
        return Ok([i, j]);
    }
    Err("local packet subcell must be either [i,j] or {sub_i,sub_j}".to_string())
}

fn summarize_cell_packet(path: &Path) -> Result<CellSummary, String> {
    let packet = read_json(path)?;
    let experiment_id = expect_str(&packet, "experiment_id")?.to_string();
    let status = expect_str(&packet, "status")?.to_string();
    let subcell = parse_subcell(&packet)?;
    let sha = verify_source_sha("local_cell_packet", &experiment_id, path)?;
    let total = expect_f64(&packet, "total_validated_length_upper")?;
    let cap = expect_f64(&packet, "exact_length_cap")?;
    let margin = expect_f64(&packet, "margin_to_cap")?;
    let residual_count = expect_usize(&packet, "residual_chain_total_closed_count")?;
    let expected_residual_count = expect_usize(&packet, "expected_residual_chain_total_count")?;
    let residual_universe_accounting_status =
        optional_string(&packet, "residual_universe_accounting_status");
    let source_backed_leaf_universe_count =
        optional_usize(&packet, "source_backed_leaf_universe_count");
    let upstream_closed_count = optional_usize(&packet, "upstream_closed_count");
    let downstream_closed_count = optional_usize(&packet, "downstream_closed_count");
    let coverage_equation_pass = optional_bool(&packet, "coverage_equation_pass");
    let source_sha_fail_count = expect_usize(&packet, "source_sha_fail_count")?;
    let ownership_duplicate_count = expect_usize(&packet, "ownership_duplicate_count")?;
    let source_filter_mismatch_count = expect_usize(&packet, "source_filter_mismatch_count")?;
    let candidate_count_drift_count = expect_usize(&packet, "candidate_count_drift_count")?;

    let mut reasons = Vec::new();
    if sha.sha_status != "PASS" {
        reasons.push("checksum failed");
    }
    if status != LOCAL_PACKET_STATUS {
        reasons.push("status is not local pass");
    }
    if subcell[0] >= SUBDIVISION || subcell[1] >= SUBDIVISION {
        reasons.push("subcell outside 8x8 atlas");
    }
    if !approx_eq(cap, EXACT_LENGTH_CAP) {
        reasons.push("cap drift");
    }
    if margin < 0.0 || total > cap {
        reasons.push("budget failure");
    }
    if source_sha_fail_count > 0 {
        reasons.push("child source checksum failure");
    }
    if ownership_duplicate_count > 0 {
        reasons.push("ownership duplicates");
    }
    if source_filter_mismatch_count > 0 {
        reasons.push("source-filter mismatch");
    }
    if candidate_count_drift_count > 0 {
        reasons.push("candidate-count drift");
    }
    if residual_count != expected_residual_count {
        reasons.push("local residual closed count does not equal expected count");
    }
    if expected_residual_count != 64 {
        if residual_universe_accounting_status.as_deref()
            != Some("RESIDUAL_UNIVERSE_ACCOUNTING_PASS_NOT_GLOBAL_PROOF")
        {
            reasons.push("non-64 residual universe lacks passing accounting cert");
        }
        if coverage_equation_pass != Some(true) {
            reasons.push("non-64 residual universe coverage equation did not pass");
        }
        if downstream_closed_count != Some(residual_count) {
            reasons.push("non-64 downstream closed count does not match local residual closure");
        }
        match (
            source_backed_leaf_universe_count,
            upstream_closed_count,
            downstream_closed_count,
        ) {
            (Some(total_leaves), Some(upstream), Some(downstream))
                if total_leaves == upstream + downstream => {}
            _ => reasons.push("non-64 source-backed leaf equation is incomplete"),
        }
    }

    let accepted = reasons.is_empty();
    let rejection_reason = if accepted {
        "none".to_string()
    } else {
        reasons.join("; ")
    };

    Ok(CellSummary {
        subcell,
        experiment_id,
        status,
        path: path.display().to_string(),
        sha_status: sha.sha_status,
        total_validated_length_upper: total,
        exact_length_cap: cap,
        margin_to_cap: margin,
        residual_chain_total_closed_count: residual_count,
        expected_residual_chain_total_count: expected_residual_count,
        residual_universe_accounting_status,
        source_backed_leaf_universe_count,
        upstream_closed_count,
        downstream_closed_count,
        coverage_equation_pass,
        source_sha_fail_count,
        ownership_duplicate_count,
        source_filter_mismatch_count,
        candidate_count_drift_count,
        accepted,
        rejection_reason,
    })
}

fn status_for(
    source_sha_fail_count: usize,
    coverage_failures: &[String],
    duplicate_cell_certificate_count: usize,
    rejected_cell_certificate_count: usize,
    missing_cell_certificate_count: usize,
    budget_failure_count: usize,
) -> &'static str {
    if source_sha_fail_count > 0 {
        "GLOBAL_ATLAS_FAIL_SOURCE_SHA"
    } else if !coverage_failures.is_empty() {
        "GLOBAL_ATLAS_FAIL_COVERAGE"
    } else if duplicate_cell_certificate_count > 0 {
        "GLOBAL_ATLAS_FAIL_DUPLICATE_COVERAGE"
    } else if rejected_cell_certificate_count > 0 {
        "GLOBAL_ATLAS_FAIL_CELL_CERTIFICATE"
    } else if missing_cell_certificate_count > 0 {
        "GLOBAL_ATLAS_FAIL_MISSING_CELL_CERTIFICATES"
    } else if budget_failure_count > 0 {
        "GLOBAL_ATLAS_FAIL_BUDGET"
    } else {
        "GLOBAL_ATLAS_PASS_NOT_FULL_EHP_PROOF"
    }
}

fn first_failed_condition(
    source_sha_fail_count: usize,
    coverage_failures: &[String],
    duplicate_cell_certificate_count: usize,
    rejected_cell_certificate_count: usize,
    missing_cells: &[[usize; 2]],
    budget_failure_count: usize,
) -> String {
    if source_sha_fail_count > 0 {
        "at least one source checksum failed".to_string()
    } else if let Some(first) = coverage_failures.first() {
        first.clone()
    } else if duplicate_cell_certificate_count > 0 {
        "duplicate local cell certificate for at least one subcell".to_string()
    } else if rejected_cell_certificate_count > 0 {
        "at least one discovered local cell certificate failed local acceptance checks".to_string()
    } else if let Some(first) = missing_cells.first() {
        format!(
            "missing theorem-grade local cell certificate for subcell ({},{})",
            first[0], first[1]
        )
    } else if budget_failure_count > 0 {
        "at least one local cell certificate exceeds the exact cap".to_string()
    } else {
        "none".to_string()
    }
}

fn write_report(result: &Value, path: &Path) -> Result<(), String> {
    let certified_cells = result["certified_cell_count"].as_u64().unwrap_or(0);
    let missing_cells = result["missing_cell_certificate_count"]
        .as_u64()
        .unwrap_or(0);
    let status = result["status"].as_str().unwrap_or("UNKNOWN");
    let first_failed = result["first_failed_condition"]
        .as_str()
        .unwrap_or("unknown");
    let report = format!(
        "# EHP114 n=14 Global Atlas Certificate Gate\n\n\
Experiment: `{}`\n\n\
## Theorem-Facing Obligation\n\n\
L33 asks whether the n=14 root-affine atlas has theorem-grade local certificates \
for every one of its `{}` cells. Each accepted cell packet must pass checksum, \
coverage, ownership, source-filter, candidate-count, and length-budget checks. \
The global certificate can pass only if all cells are present and accepted.\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Global atlas certificate pass: `{}`\n\
- Root-affine skeleton cell count: `{}`\n\
- Certified theorem-grade cells: `{}`\n\
- Missing theorem-grade cells: `{}`\n\
- Discovered local packet files: `{}`\n\
- Selected current packet files: `{}`\n\
- Superseded packet files ignored: `{}`\n\
- Rejected discovered cell packets: `{}`\n\
- Duplicate local cell certificates: `{}`\n\
- Source SHA fail count: `{}`\n\
- Coverage failure count: `{}`\n\
- Max certified cell length upper: `{}`\n\
- Exact cap: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
This run is an atlas-level integration certificate. It does not recompute the \
local geometry; it checks that the root-affine skeleton enumerates exactly 64 \
cells and that the selected latest local packet for each cell passes checksum, \
coverage, ownership, source-filter, candidate-drift, and length-budget checks. \
Older packet versions are preserved as immutable evidence but ignored as \
superseded current coverage.\n\n\
## Next Dependency\n\n\
If this gate passes, the next dependency is an independent theorem-packet audit \
over the atlas statement and then finite-degree/high-degree bridge accounting. \
Do not promote a full EHP114 claim from this local n=14 atlas gate alone.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        EXPECTED_CELL_COUNT,
        status,
        result["global_atlas_certificate_pass"],
        result["root_affine_skeleton_cell_count"],
        certified_cells,
        missing_cells,
        result["discovered_cell_packet_count"],
        result["selected_cell_packet_count"],
        result["superseded_cell_packet_count"],
        result["rejected_cell_certificate_count"],
        result["duplicate_cell_certificate_count"],
        result["source_sha_fail_count"],
        result["coverage_failure_count"],
        result["max_certified_cell_length_upper"],
        result["exact_length_cap"],
        first_failed,
        result["claim_ceiling"]
            .as_str()
            .unwrap_or("local integration gate only"),
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

    let root_source = read_json(&cfg.root_affine_source)?;
    let root_source_sha = verify_source_sha(
        "root_affine_coverage_skeleton",
        ROOT_AFFINE_SOURCE_ID,
        &cfg.root_affine_source,
    )?;
    let (coverage_failures, root_affine_cells) = validate_root_affine_source(&root_source)?;

    let packet_paths = discover_cell_packets(&cfg)?;
    let discovered_cell_packet_count = packet_paths.len();
    let mut source_sha_status = vec![root_source_sha];
    let mut cell_summaries = Vec::new();
    for path in packet_paths {
        let summary = summarize_cell_packet(&path)?;
        source_sha_status.push(verify_source_sha(
            "local_cell_packet",
            &summary.experiment_id,
            Path::new(&summary.path),
        )?);
        cell_summaries.push(summary);
    }

    let source_sha_fail_count = source_sha_status
        .iter()
        .filter(|row| row.sha_status != "PASS")
        .count();

    let mut discovered_by_cell: BTreeMap<(usize, usize), Vec<CellSummary>> = BTreeMap::new();
    for summary in cell_summaries {
        discovered_by_cell
            .entry((summary.subcell[0], summary.subcell[1]))
            .or_default()
            .push(summary);
    }

    let mut superseded_cell_packets = Vec::new();
    let mut selected_by_cell: BTreeMap<(usize, usize), Vec<CellSummary>> = BTreeMap::new();
    for (cell, mut summaries) in discovered_by_cell {
        summaries.sort_by(|a, b| {
            a.experiment_id
                .cmp(&b.experiment_id)
                .then_with(|| a.path.cmp(&b.path))
        });
        if let Some(selected) = summaries.pop() {
            for superseded in summaries {
                superseded_cell_packets.push(json!(superseded));
            }
            selected_by_cell.entry(cell).or_default().push(selected);
        }
    }
    let superseded_cell_packet_count = superseded_cell_packets.len();

    let expected = expected_cells();
    let mut certified_cells = Vec::new();
    let mut rejected_cells = Vec::new();
    let mut duplicate_cells = Vec::new();
    let mut budget_failure_count = 0usize;
    for ((i, j), summaries) in &selected_by_cell {
        if summaries.len() > 1 {
            duplicate_cells.push(json!({
                "subcell": [i, j],
                "certificate_count": summaries.len(),
                "paths": summaries.iter().map(|s| s.path.clone()).collect::<Vec<_>>()
            }));
        }
        for summary in summaries {
            if !summary.accepted {
                rejected_cells.push(json!(summary));
            }
            if summary.total_validated_length_upper > summary.exact_length_cap
                || summary.margin_to_cap < 0.0
            {
                budget_failure_count += 1;
            }
        }
        if summaries.len() == 1 && summaries[0].accepted && expected.contains(&(*i, *j)) {
            certified_cells.push(summaries[0].clone());
        }
    }

    let certified_set: BTreeSet<(usize, usize)> = certified_cells
        .iter()
        .map(|summary| (summary.subcell[0], summary.subcell[1]))
        .collect();
    let missing_cells: Vec<[usize; 2]> = expected
        .difference(&certified_set)
        .map(|(i, j)| [*i, *j])
        .collect();
    let max_certified_cell_length_upper = certified_cells
        .iter()
        .map(|summary| summary.total_validated_length_upper)
        .fold(0.0_f64, f64::max);
    let min_certified_cell_margin = certified_cells
        .iter()
        .map(|summary| summary.margin_to_cap)
        .fold(f64::INFINITY, f64::min);

    let status = status_for(
        source_sha_fail_count,
        &coverage_failures,
        duplicate_cells.len(),
        rejected_cells.len(),
        missing_cells.len(),
        budget_failure_count,
    );
    let first_failed = first_failed_condition(
        source_sha_fail_count,
        &coverage_failures,
        duplicate_cells.len(),
        rejected_cells.len(),
        &missing_cells,
        budget_failure_count,
    );

    let result = json!({
        "experiment_id": cfg.experiment_id,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "degree": DEGREE,
        "eps": EPS,
        "subdivision": SUBDIVISION,
        "expected_cell_count": EXPECTED_CELL_COUNT,
        "root_affine_source_experiment_id": ROOT_AFFINE_SOURCE_ID,
        "root_affine_source_path": cfg.root_affine_source.display().to_string(),
        "root_affine_skeleton_cell_count": root_affine_cells.len(),
        "root_affine_skeleton_status": root_source["status"].clone(),
        "root_affine_skeleton_claim_ceiling": root_source["claim_ceiling"].clone(),
        "source_sha_status": source_sha_status,
        "source_sha_fail_count": source_sha_fail_count,
        "discovered_cell_packet_count": discovered_cell_packet_count,
        "selected_cell_packet_count": selected_by_cell.len(),
        "superseded_cell_packet_count": superseded_cell_packet_count,
        "superseded_cell_packets": superseded_cell_packets,
        "coverage_failures": coverage_failures,
        "coverage_failure_count": coverage_failures.len(),
        "certified_cells": certified_cells,
        "certified_cell_count": certified_set.len(),
        "missing_cells": missing_cells,
        "missing_cell_certificate_count": missing_cells.len(),
        "duplicate_cell_certificates": duplicate_cells,
        "duplicate_cell_certificate_count": duplicate_cells.len(),
        "rejected_cell_certificates": rejected_cells,
        "rejected_cell_certificate_count": rejected_cells.len(),
        "budget_failure_count": budget_failure_count,
        "max_certified_cell_length_upper": max_certified_cell_length_upper,
        "min_certified_cell_margin_to_cap": if certified_set.is_empty() { Value::Null } else { json!(min_certified_cell_margin) },
        "exact_length_cap": EXACT_LENGTH_CAP,
        "global_atlas_certificate_pass": status == "GLOBAL_ATLAS_PASS_NOT_FULL_EHP_PROOF",
        "first_failed_condition": first_failed,
        "proof_obligations": [
            "root-affine skeleton enumerates exactly 64 cells",
            "every cell has one theorem-grade local packet",
            "every local packet checksum passes",
            "no duplicate cell certificates",
            "no ownership duplicates inside any local packet",
            "no source-filter mismatches inside any local packet",
            "no candidate-count drift inside any local packet",
            "each local validated length upper is below the exact cap"
        ],
        "what_this_rules_out": "The current local hard-cell packet alone is not sufficient for a global n=14 atlas certificate.",
        "what_this_does_not_rule_out": "Missing cells may be certifiable by generating theorem-grade residual-chain packets for those cells; this result is not a geometric counterexample.",
        "next_dependency": "Generate theorem-shaped local residual-chain packets for the missing n=14 root-affine cells, then rerun this L33 atlas gate.",
        "claim_ceiling": "L33 global atlas integration gate only. This run is not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. It may pass only after all 64 local theorem-grade cell packets are present and accepted."
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
            "experiment_id": cfg.experiment_id,
            "status": result["status"],
            "global_atlas_certificate_pass": result["global_atlas_certificate_pass"],
            "root_affine_skeleton_cell_count": result["root_affine_skeleton_cell_count"],
            "certified_cell_count": result["certified_cell_count"],
            "missing_cell_certificate_count": result["missing_cell_certificate_count"],
            "source_sha_fail_count": result["source_sha_fail_count"],
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
