//! EHP #114 n=14 parameterized local-pipeline smoke gate.
//!
//! This is a software-proof rigor gate. It checks that the local residual-chain
//! binaries expose explicit subcell arguments and no longer contain fixed
//! hard-cell constants. It does not generate missing-cell certificates.

use ehp_n3_poc::ehp114_n14_cell::{
    ensure_source_subcell_matches, source_subcell_from_value, CellSpec,
};
use serde::Serialize;
use serde_json::{json, Value};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const DEFAULT_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-PARAMETERIZED-LOCAL-PIPELINE-SMOKE-20260506-01";
const L33A_ID: &str = "EXP-MATH-EHP114-N14-MISSING-CELL-LOCAL-PACKETS-20260506-02";
const DEFAULT_L33A: &str = "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-MISSING-CELL-LOCAL-PACKETS-20260506-02/EXP-MATH-EHP114-N14-MISSING-CELL-LOCAL-PACKETS-20260506-02_RESULTS.json";
const DEFAULT_OUTDIR: &str = "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-PARAMETERIZED-LOCAL-PIPELINE-SMOKE-20260506-01";
const DEFAULT_SRC_BIN_DIR: &str = "src/bin";
const HARD_CELL_L18_SOURCE: &str = "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01_RESULTS.json";

const REQUIRED_BINS: &[&str] = &[
    "ehp114_n14_global_critical_point_exclusion_target.rs",
    "ehp114_n14_regular_residual_decomposition.rs",
    "ehp114_n14_critical_candidate_affine_gradient.rs",
    "ehp114_n14_pprime_root_location.rs",
    "ehp114_n14_sharp_pprime_root_location.rs",
    "ehp114_n14_residual_chain_integration_cert.rs",
    "ehp114_n14_local_hard_cell_certificate_packet.rs",
];

#[derive(Clone, Debug)]
struct Config {
    l33a: PathBuf,
    out_dir: PathBuf,
    src_bin_dir: PathBuf,
    experiment_id: String,
    smoke_cell: CellSpec,
}

#[derive(Clone, Debug, Serialize)]
struct BinaryAudit {
    binary: String,
    path: String,
    source_exists: bool,
    has_sub_i_cli: bool,
    has_sub_j_cli: bool,
    has_source_subcell_contract: bool,
    hardcoded_subcell_evidence: Vec<String>,
    blocking: bool,
}

fn parse_args() -> Result<Config, String> {
    let mut cfg = Config {
        l33a: PathBuf::from(DEFAULT_L33A),
        out_dir: PathBuf::from(DEFAULT_OUTDIR),
        src_bin_dir: PathBuf::from(DEFAULT_SRC_BIN_DIR),
        experiment_id: DEFAULT_EXPERIMENT_ID.to_string(),
        smoke_cell: CellSpec::new(0, 0)?,
    };
    let args: Vec<String> = env::args().collect();
    let mut sub_i = cfg.smoke_cell.sub_i;
    let mut sub_j = cfg.smoke_cell.sub_j;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--l33a" => {
                cfg.l33a = PathBuf::from(next_string(&args, i, "--l33a")?);
                i += 2;
            }
            "--outdir" => {
                cfg.out_dir = PathBuf::from(next_string(&args, i, "--outdir")?);
                i += 2;
            }
            "--src-bin-dir" => {
                cfg.src_bin_dir = PathBuf::from(next_string(&args, i, "--src-bin-dir")?);
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
    cfg.smoke_cell = CellSpec::new(sub_i, sub_j)?;
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

fn audit_binary(src_bin_dir: &Path, binary: &str) -> BinaryAudit {
    let path = src_bin_dir.join(binary);
    let Ok(text) = fs::read_to_string(&path) else {
        return BinaryAudit {
            binary: binary.to_string(),
            path: path.display().to_string(),
            source_exists: false,
            has_sub_i_cli: false,
            has_sub_j_cli: false,
            has_source_subcell_contract: false,
            hardcoded_subcell_evidence: vec!["source file missing".to_string()],
            blocking: true,
        };
    };
    let hardcoded_subcell_evidence: Vec<String> = text
        .lines()
        .enumerate()
        .filter_map(|(idx, line)| {
            let trimmed = line.trim();
            if trimmed.starts_with("const SUB_I:")
                || trimmed.starts_with("const SUB_J:")
                || trimmed.starts_with("const SUBCELL:")
            {
                Some(format!("{}:{}", idx + 1, trimmed))
            } else {
                None
            }
        })
        .collect();
    let has_sub_i_cli = text.contains("--sub-i");
    let has_sub_j_cli = text.contains("--sub-j");
    let has_source_subcell_contract = text.contains("ensure_source_subcell_matches");
    let blocking = !hardcoded_subcell_evidence.is_empty()
        || !has_sub_i_cli
        || !has_sub_j_cli
        || !has_source_subcell_contract;
    BinaryAudit {
        binary: binary.to_string(),
        path: path.display().to_string(),
        source_exists: true,
        has_sub_i_cli,
        has_sub_j_cli,
        has_source_subcell_contract,
        hardcoded_subcell_evidence,
        blocking,
    }
}

fn first_missing_cell(l33a: &Value) -> Result<CellSpec, String> {
    let cells = l33a
        .get("missing_cells")
        .and_then(Value::as_array)
        .ok_or_else(|| "L33A source missing missing_cells array".to_string())?;
    let first = cells
        .first()
        .and_then(Value::as_array)
        .ok_or_else(|| "L33A missing_cells[0] is not a cell array".to_string())?;
    if first.len() != 2 {
        return Err("L33A missing_cells[0] does not have length 2".to_string());
    }
    CellSpec::new(
        first[0]
            .as_u64()
            .ok_or_else(|| "L33A missing_cells[0][0] is not an integer".to_string())?
            as usize,
        first[1]
            .as_u64()
            .ok_or_else(|| "L33A missing_cells[0][1] is not an integer".to_string())?
            as usize,
    )
}

fn write_report(result: &Value, path: &Path) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Parameterized Local Pipeline Smoke\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Smoke subcell: `{}`\n\
- Required binaries: `{}`\n\
- Hardcoded blocker count: `{}`\n\
- Source contract status: `{}`\n\
- First failed condition: `{}`\n\n\
## Meaning\n\n\
This is a software-proof rigor gate. It checks that the seven local residual-chain \
binaries now take explicit subcell arguments and contain no fixed hard-cell constants. \
It does not certify any missing cell. If the first non-hard-cell source artifact is \
absent, that is reported as the blocker rather than silently reusing hard-cell data.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        result["status"].as_str().unwrap_or("UNKNOWN"),
        result["smoke_subcell"]["cell_tag"]
            .as_str()
            .unwrap_or("unknown"),
        result["required_binary_count"],
        result["hardcoded_blocker_count"],
        result["source_contract_status"]
            .as_str()
            .unwrap_or("UNKNOWN"),
        result["first_failed_condition"]
            .as_str()
            .unwrap_or("unknown"),
        result["claim_ceiling"]
            .as_str()
            .unwrap_or("local smoke gate only"),
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

    let l33a = read_json(&cfg.l33a)?;
    let l33a_experiment_id = l33a
        .get("experiment_id")
        .and_then(Value::as_str)
        .unwrap_or("UNKNOWN");
    if l33a_experiment_id != L33A_ID {
        return Err(format!(
            "unexpected L33A source experiment_id: {l33a_experiment_id}"
        ));
    }
    let l33a_first_missing = first_missing_cell(&l33a)?;
    let smoke_cell = cfg.smoke_cell;
    let audits: Vec<BinaryAudit> = REQUIRED_BINS
        .iter()
        .map(|binary| audit_binary(&cfg.src_bin_dir, binary))
        .collect();
    let hardcoded_blocker_count = audits.iter().filter(|row| row.blocking).count();

    let hard_source_path = PathBuf::from(HARD_CELL_L18_SOURCE);
    let hard_source_exists = hard_source_path.exists();
    let (source_contract_status, source_declared_subcell) = if hard_source_exists {
        let source = read_json(&hard_source_path)?;
        let declared = source_subcell_from_value(&source)?.map(|cell| cell.tag());
        let status = match ensure_source_subcell_matches(&source, smoke_cell) {
            Ok(status) => status,
            Err(err) if err.starts_with("SOURCE_SUBCELL_MISMATCH") => {
                "SOURCE_SUBCELL_MISMATCH".to_string()
            }
            Err(err) if err.starts_with("SOURCE_SUBCELL_UNVERIFIED") => {
                "SOURCE_SUBCELL_UNVERIFIED".to_string()
            }
            Err(err) => err,
        };
        (status, declared)
    } else {
        ("SOURCE_ARTIFACT_MISSING".to_string(), None)
    };

    let expected_per_cell_l18 = format!(
        "EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-{}-20260506-01_RESULTS.json",
        smoke_cell.tag()
    );
    let per_cell_source_exists = false;
    let status = if hardcoded_blocker_count > 0 {
        "PARAMETERIZED_LOCAL_PIPELINE_FAIL_HARDCODED_SUBCELL"
    } else if !per_cell_source_exists {
        "PARAMETERIZED_LOCAL_PIPELINE_BLOCKED_SOURCE_ARTIFACT_MISSING"
    } else if source_contract_status == "SOURCE_SUBCELL_MISMATCH" {
        "PARAMETERIZED_LOCAL_PIPELINE_FAIL_SOURCE_SUBCELL_MISMATCH"
    } else {
        "PARAMETERIZED_LOCAL_PIPELINE_READY"
    };
    let first_failed_condition = if hardcoded_blocker_count > 0 {
        "at least one proof-facing local binary still has a hard-coded subcell or lacks the subcell/source contract CLI".to_string()
    } else if !per_cell_source_exists {
        format!(
            "per-cell upstream source artifact is missing for {}; expected first source like {}",
            smoke_cell.tag(),
            expected_per_cell_l18
        )
    } else if source_contract_status == "SOURCE_SUBCELL_MISMATCH" {
        format!("hard-cell source cannot be reused for {}", smoke_cell.tag())
    } else {
        "none".to_string()
    };

    let result = json!({
        "experiment_id": cfg.experiment_id,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "source_l33a_experiment_id": l33a_experiment_id,
        "source_l33a_path": cfg.l33a,
        "missing_cell_count_from_l33a": l33a["missing_cell_count"].clone(),
        "l33a_first_missing_cell": {
            "sub_i": l33a_first_missing.sub_i,
            "sub_j": l33a_first_missing.sub_j,
            "cell_tag": l33a_first_missing.tag()
        },
        "smoke_subcell": {
            "sub_i": smoke_cell.sub_i,
            "sub_j": smoke_cell.sub_j,
            "cell_tag": smoke_cell.tag()
        },
        "required_binary_count": REQUIRED_BINS.len(),
        "hardcoded_blocker_count": hardcoded_blocker_count,
        "binary_audit": audits,
        "source_contract_status": source_contract_status,
        "hard_cell_source_probe": {
            "path": hard_source_path,
            "exists": hard_source_exists,
            "declared_subcell": source_declared_subcell,
            "expected_failure_for_nonhard_smoke": "SOURCE_SUBCELL_MISMATCH or SOURCE_SUBCELL_UNVERIFIED"
        },
        "per_cell_source_probe": {
            "exists": per_cell_source_exists,
            "expected_first_source": expected_per_cell_l18,
            "reason": "L33B parameterizes the local pipeline only; it intentionally does not generate the upstream per-cell L18 source artifacts."
        },
        "first_failed_condition": first_failed_condition,
        "what_this_rules_out": "The seven local residual-chain binaries are no longer allowed to hide fixed subcell constants when the smoke audit reports zero blockers.",
        "what_this_does_not_rule_out": "This does not certify any of the 63 missing cells and does not prove that per-cell upstream source artifacts exist.",
        "claim_ceiling": "Parameterized local-pipeline smoke gate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate."
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
            "smoke_subcell": result["smoke_subcell"],
            "hardcoded_blocker_count": result["hardcoded_blocker_count"],
            "source_contract_status": result["source_contract_status"],
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
