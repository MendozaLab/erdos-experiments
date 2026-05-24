//! EHP #114 n=14 per-cell source-generation smoke gate.
//!
//! This is a software-proof rigor gate for L33C. It checks whether the direct
//! upstream source chain for a non-hard cell can start honestly. It refuses to
//! substitute hard-cell artifacts for `CELL-00-00`.

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
    "EXP-MATH-EHP114-N14-PER-CELL-SOURCE-GENERATION-SMOKE-20260506-01";
const DEFAULT_OUTDIR: &str =
    "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-PER-CELL-SOURCE-GENERATION-SMOKE-20260506-01";
const DEFAULT_SRC_BIN_DIR: &str = "src/bin";

const UPSTREAM_BINS: &[&str] = &[
    "ehp114_n14_branch_isolation_collar_atlas.rs",
    "ehp114_n14_normal_collar_critical_exclusion_pilot.rs",
    "ehp114_n14_third_order_collar_remainder_pilot.rs",
    "ehp114_n14_global_critical_point_exclusion_target.rs",
];

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    source_was_explicit: bool,
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
    has_experiment_id_cli: bool,
    has_source_cli: bool,
    has_outdir_cli: bool,
    has_source_subcell_contract: bool,
    hardcoded_subcell_evidence: Vec<String>,
    blocking: bool,
}

fn expected_slab_source(cell: CellSpec) -> PathBuf {
    let id = format!(
        "EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-{}-20260506-01",
        cell.tag()
    );
    PathBuf::from(format!(
        "../../Erdos114/validated_length/{id}/{id}_RESULTS.json"
    ))
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut smoke_cell = CellSpec::new(0, 0)?;
    let mut source: Option<PathBuf> = None;
    let mut cfg = Config {
        source: expected_slab_source(smoke_cell),
        source_was_explicit: false,
        out_dir: PathBuf::from(DEFAULT_OUTDIR),
        src_bin_dir: PathBuf::from(DEFAULT_SRC_BIN_DIR),
        experiment_id: DEFAULT_EXPERIMENT_ID.to_string(),
        smoke_cell,
    };
    let mut sub_i = smoke_cell.sub_i;
    let mut sub_j = smoke_cell.sub_j;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--source" => {
                source = Some(PathBuf::from(next_string(&args, i, "--source")?));
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
    smoke_cell = CellSpec::new(sub_i, sub_j)?;
    cfg.smoke_cell = smoke_cell;
    cfg.source_was_explicit = source.is_some();
    cfg.source = source.unwrap_or_else(|| expected_slab_source(smoke_cell));
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
            has_experiment_id_cli: false,
            has_source_cli: false,
            has_outdir_cli: false,
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
    let has_experiment_id_cli = text.contains("--experiment-id");
    let has_source_cli = text.contains("--source");
    let has_outdir_cli = text.contains("--outdir");
    let has_source_subcell_contract = text.contains("ensure_source_subcell_matches");
    let blocking = !hardcoded_subcell_evidence.is_empty()
        || !has_sub_i_cli
        || !has_sub_j_cli
        || !has_experiment_id_cli
        || !has_source_cli
        || !has_outdir_cli
        || !has_source_subcell_contract;
    BinaryAudit {
        binary: binary.to_string(),
        path: path.display().to_string(),
        source_exists: true,
        has_sub_i_cli,
        has_sub_j_cli,
        has_experiment_id_cli,
        has_source_cli,
        has_outdir_cli,
        has_source_subcell_contract,
        hardcoded_subcell_evidence,
        blocking,
    }
}

fn source_contract_status(path: &Path, cell: CellSpec) -> Result<(String, Option<String>), String> {
    if !path.exists() {
        return Ok(("SOURCE_ARTIFACT_MISSING".to_string(), None));
    }
    let source = read_json(path)?;
    let declared = source_subcell_from_value(&source)?.map(|source_cell| source_cell.tag());
    let status = match ensure_source_subcell_matches(&source, cell) {
        Ok(status) => status,
        Err(err) if err.starts_with("SOURCE_SUBCELL_MISMATCH") => {
            "SOURCE_SUBCELL_MISMATCH".to_string()
        }
        Err(err) if err.starts_with("SOURCE_SUBCELL_UNVERIFIED") => {
            "SOURCE_SUBCELL_UNVERIFIED".to_string()
        }
        Err(err) => err,
    };
    Ok((status, declared))
}

fn write_report(result: &Value, path: &Path) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Per-Cell Source Generation Smoke\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Smoke subcell: `{}`\n\
- Upstream binaries audited: `{}`\n\
- Hardcoded blocker count: `{}`\n\
- Source contract status: `{}`\n\
- First missing source artifact: `{}`\n\
- First failed condition: `{}`\n\n\
## Meaning\n\n\
This L33C smoke gate asks whether the first missing cell can start the real \
source chain. It does not generate a certificate and does not rename hard-cell \
data as cell-local evidence. If the `CELL-00-00` slab source is absent or a \
source declares the wrong subcell, the proof pipeline stops there.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        result["status"].as_str().unwrap_or("UNKNOWN"),
        result["smoke_subcell"]["cell_tag"]
            .as_str()
            .unwrap_or("unknown"),
        result["upstream_binary_count"],
        result["hardcoded_blocker_count"],
        result["source_contract_status"]
            .as_str()
            .unwrap_or("UNKNOWN"),
        result["first_missing_source_artifact"]
            .as_str()
            .unwrap_or("none"),
        result["first_failed_condition"].as_str().unwrap_or("none"),
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

    let audits: Vec<BinaryAudit> = UPSTREAM_BINS
        .iter()
        .map(|binary| audit_binary(&cfg.src_bin_dir, binary))
        .collect();
    let hardcoded_blocker_count = audits.iter().filter(|row| row.blocking).count();
    let expected_source = expected_slab_source(cfg.smoke_cell);
    let expected_source_exists = expected_source.exists();
    let actual_source_exists = cfg.source.exists();
    let (source_contract_status, source_declared_subcell) =
        source_contract_status(&cfg.source, cfg.smoke_cell)?;

    let first_missing_source_artifact = if !expected_source_exists {
        Some(expected_source.display().to_string())
    } else if !actual_source_exists {
        Some(cfg.source.display().to_string())
    } else {
        None
    };
    let source_contract_failed = source_contract_status == "SOURCE_SUBCELL_MISMATCH"
        || source_contract_status == "SOURCE_SUBCELL_UNVERIFIED";
    let status = if hardcoded_blocker_count > 0 {
        "PER_CELL_SOURCE_BLOCKED_UPSTREAM_PARAMETERIZATION"
    } else if source_contract_failed {
        "PER_CELL_SOURCE_FAIL_SOURCE_SUBCELL_MISMATCH"
    } else if first_missing_source_artifact.is_some() {
        "PER_CELL_SOURCE_BLOCKED_SOURCE_ARTIFACT_MISSING"
    } else {
        "PER_CELL_SOURCE_CHAIN_READY"
    };
    let first_failed_condition = if hardcoded_blocker_count > 0 {
        "at least one upstream binary lacks the cell/source CLI contract or still has fixed subcell constants".to_string()
    } else if source_contract_failed {
        format!(
            "requested {} but source contract returned {}",
            cfg.smoke_cell.tag(),
            source_contract_status
        )
    } else if let Some(path) = &first_missing_source_artifact {
        format!(
            "required first slab source artifact is missing for {}: {}",
            cfg.smoke_cell.tag(),
            path
        )
    } else {
        "none".to_string()
    };

    let result = json!({
        "experiment_id": cfg.experiment_id,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "smoke_subcell": {
            "sub_i": cfg.smoke_cell.sub_i,
            "sub_j": cfg.smoke_cell.sub_j,
            "cell_tag": cfg.smoke_cell.tag()
        },
        "upstream_binary_count": UPSTREAM_BINS.len(),
        "upstream_binary_audit": audits,
        "hardcoded_blocker_count": hardcoded_blocker_count,
        "expected_first_source_artifact": expected_source,
        "expected_first_source_exists": expected_source_exists,
        "source_results_path": cfg.source,
        "source_was_explicit": cfg.source_was_explicit,
        "source_exists": actual_source_exists,
        "source_contract_status": source_contract_status,
        "source_declared_subcell": source_declared_subcell,
        "first_missing_source_artifact": first_missing_source_artifact,
        "first_failed_condition": first_failed_condition,
        "what_this_rules_out": "The direct L33C chain cannot silently begin from a hard-cell or metadata-free source when a non-hard cell is requested.",
        "what_this_does_not_rule_out": "This does not rule out CELL-00-00 or any other missing cell mathematically; it only shows that the required per-cell source artifact has not yet been generated.",
        "next_dependency": "Parameterize/generate the slab validated-length source artifact for CELL-00-00, then rerun the direct chain smoke.",
        "claim_ceiling": "Per-cell source-generation smoke gate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate."
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
            "smoke_subcell": result["smoke_subcell"],
            "hardcoded_blocker_count": result["hardcoded_blocker_count"],
            "source_contract_status": result["source_contract_status"],
            "first_missing_source_artifact": result["first_missing_source_artifact"],
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
