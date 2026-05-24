//! EHP #114 n=14 missing-cell local packet generation gate.
//!
//! This is the L33A bridge from the global atlas gate to actual per-cell
//! theorem packets. It consumes the L33 result, emits the missing-cell worklist,
//! and audits whether the local residual-chain pipeline is parameterized enough
//! to run honestly beyond the hard cell `(6,4)`.
//!
//! It intentionally does not synthesize local packets while the pipeline is
//! still hard-coded to `(6,4)`.

use serde::Serialize;
use serde_json::{json, Value};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const EXPERIMENT_ID: &str = "EXP-MATH-EHP114-N14-MISSING-CELL-LOCAL-PACKETS-20260506-02";
const L33_ID: &str = "EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-01";
const EXPECTED_MISSING_CELL_COUNT: usize = 63;

#[derive(Clone, Debug)]
struct Config {
    l33: PathBuf,
    src_bin_dir: PathBuf,
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
struct PipelineBinaryAudit {
    binary_name: String,
    path: String,
    role: String,
    hardcoded_subcell: bool,
    has_subcell_cli: bool,
    blocking: bool,
    evidence: Vec<String>,
}

fn parse_args() -> Result<Config, String> {
    let mut cfg = Config {
        l33: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-01/EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-01_RESULTS.json"),
        src_bin_dir: PathBuf::from("src/bin"),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-MISSING-CELL-LOCAL-PACKETS-20260506-02"),
    };

    let args: Vec<String> = env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--l33" => {
                cfg.l33 = next_path(&args, i, "--l33")?;
                i += 2;
            }
            "--src-bin-dir" => {
                cfg.src_bin_dir = next_path(&args, i, "--src-bin-dir")?;
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

fn unix_timestamp_string() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs()
        .to_string()
}

fn missing_cells(l33: &Value) -> Result<Vec<[usize; 2]>, String> {
    let arr = l33
        .get("missing_cells")
        .and_then(Value::as_array)
        .ok_or_else(|| "L33 missing_cells is not an array".to_string())?;
    let mut cells = Vec::with_capacity(arr.len());
    for row in arr {
        let pair = row
            .as_array()
            .ok_or_else(|| "L33 missing_cells entry is not an array".to_string())?;
        if pair.len() != 2 {
            return Err("L33 missing_cells entry is not length 2".to_string());
        }
        let i = pair[0]
            .as_u64()
            .ok_or_else(|| "L33 missing cell i is not integer".to_string())?
            as usize;
        let j = pair[1]
            .as_u64()
            .ok_or_else(|| "L33 missing cell j is not integer".to_string())?
            as usize;
        cells.push([i, j]);
    }
    Ok(cells)
}

fn audit_binary(
    src_bin_dir: &Path,
    name: &str,
    role: &str,
    required: bool,
) -> Result<PipelineBinaryAudit, String> {
    let path = src_bin_dir.join(name);
    let text = fs::read_to_string(&path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    let mut evidence = Vec::new();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with("const SUB_I:")
            || trimmed.starts_with("const SUB_J:")
            || trimmed.starts_with("const SUBCELL:")
        {
            evidence.push(trimmed.to_string());
        }
    }
    let hardcoded_subcell = !evidence.is_empty();
    let has_subcell_cli =
        text.contains("--sub-i") || text.contains("--sub-j") || text.contains("--subcell");
    let blocking = required && hardcoded_subcell && !has_subcell_cli;
    Ok(PipelineBinaryAudit {
        binary_name: name.trim_end_matches(".rs").to_string(),
        path: path.display().to_string(),
        role: role.to_string(),
        hardcoded_subcell,
        has_subcell_cli,
        blocking,
        evidence,
    })
}

fn audit_pipeline(src_bin_dir: &Path) -> Result<Vec<PipelineBinaryAudit>, String> {
    let binaries = [
        (
            "ehp114_n14_global_critical_point_exclusion_target.rs",
            "L21 source classifier: regular versus critical residual regions",
            true,
        ),
        (
            "ehp114_n14_regular_residual_decomposition.rs",
            "L24 regular monotone residual closure",
            true,
        ),
        (
            "ehp114_n14_critical_candidate_affine_gradient.rs",
            "L27 direct p-prime critical-candidate closure",
            true,
        ),
        (
            "ehp114_n14_pprime_root_location.rs",
            "L28 p-prime root-location closure",
            true,
        ),
        (
            "ehp114_n14_sharp_pprime_root_location.rs",
            "L29 sharp Taylor/Rouche p-prime root-location closure",
            true,
        ),
        (
            "ehp114_n14_residual_chain_integration_cert.rs",
            "L31 residual-chain integration packet",
            true,
        ),
        (
            "ehp114_n14_local_hard_cell_certificate_packet.rs",
            "L32 local theorem packet",
            true,
        ),
        (
            "ehp114_n14_global_atlas_certificate.rs",
            "L33 atlas gate",
            false,
        ),
    ];
    binaries
        .iter()
        .map(|(name, role, required)| audit_binary(src_bin_dir, name, role, *required))
        .collect()
}

fn status_for(
    source_sha_fail_count: usize,
    missing_count: usize,
    blocking_count: usize,
) -> &'static str {
    if source_sha_fail_count > 0 {
        "MISSING_CELL_LOCAL_PACKETS_FAIL_SOURCE_SHA"
    } else if missing_count == 0 {
        "MISSING_CELL_LOCAL_PACKETS_NO_MISSING_CELLS"
    } else if blocking_count > 0 {
        "MISSING_CELL_LOCAL_PACKETS_BLOCKED_PIPELINE_PARAMETERIZATION"
    } else {
        "MISSING_CELL_LOCAL_PACKETS_READY_TO_RUN"
    }
}

fn write_report(result: &Value, path: &Path) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Missing-Cell Local Packet Gate\n\n\
Experiment: `{}`\n\n\
## Meaning\n\n\
L33 showed that the global n=14 atlas skeleton has all 64 cells, but only one \
theorem-grade local packet exists. This L33A gate turns that into a concrete \
worklist and checks whether the local residual-chain pipeline can be run on \
cells other than `(6,4)` without code changes.\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Missing cells from L33: `{}`\n\
- Source SHA fail count: `{}`\n\
- Required binaries audited: `{}`\n\
- Blocking hard-coded binaries: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
The current software cannot honestly generate the 63 missing local theorem \
packets because the proof-facing residual-chain binaries still carry hard-coded \
`(6,4)` subcell constants. The next repair is parameterization, not another \
global run: add `--sub-i`, `--sub-j`, per-cell experiment IDs, and per-cell \
source paths to the required L21/L24/L27/L28/L29/L31/L32 binaries, then rerun \
this gate.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"].as_str().unwrap_or(EXPERIMENT_ID),
        result["status"].as_str().unwrap_or("UNKNOWN"),
        result["missing_cell_count"],
        result["source_sha_fail_count"],
        result["required_binary_count"],
        result["blocking_hardcoded_binary_count"],
        result["first_failed_condition"]
            .as_str()
            .unwrap_or("unknown"),
        result["claim_ceiling"]
            .as_str()
            .unwrap_or("local packet generation gate only"),
    );
    fs::write(path, report).map_err(|err| format!("failed to write {}: {err}", path.display()))
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let cfg = parse_args()?;
    fs::create_dir_all(&cfg.out_dir)
        .map_err(|err| format!("failed to create {}: {err}", cfg.out_dir.display()))?;

    let result_path = cfg.out_dir.join(format!("{EXPERIMENT_ID}_RESULTS.json"));
    let report_path = cfg.out_dir.join(format!("{EXPERIMENT_ID}_REPORT.md"));
    let sha_path = cfg.out_dir.join(format!("{EXPERIMENT_ID}_RESULTS.sha256"));
    for path in [&result_path, &report_path, &sha_path] {
        if path.exists() {
            return Err(format!(
                "refusing to overwrite existing artifact: {}",
                path.display()
            ));
        }
    }

    let l33 = read_json(&cfg.l33)?;
    let l33_sha = verify_source_sha("l33_global_atlas_gate", L33_ID, &cfg.l33)?;
    let source_sha_fail_count = usize::from(l33_sha.sha_status != "PASS");

    let mut validation_failures = Vec::new();
    if expect_str(&l33, "experiment_id")? != L33_ID {
        validation_failures.push("L33 experiment_id mismatch".to_string());
    }
    if expect_str(&l33, "status")? != "GLOBAL_ATLAS_FAIL_MISSING_CELL_CERTIFICATES" {
        validation_failures.push("L33 is not in the expected missing-cell status".to_string());
    }
    if expect_usize(&l33, "root_affine_skeleton_cell_count")? != 64 {
        validation_failures.push("L33 root-affine skeleton cell count is not 64".to_string());
    }
    if expect_usize(&l33, "certified_cell_count")? != 1 {
        validation_failures.push("L33 certified cell count is not 1".to_string());
    }

    let cells = missing_cells(&l33)?;
    if cells.len() != EXPECTED_MISSING_CELL_COUNT {
        validation_failures.push(format!(
            "L33 missing cell count is {}, expected {}",
            cells.len(),
            EXPECTED_MISSING_CELL_COUNT
        ));
    }
    let pipeline_audit = audit_pipeline(&cfg.src_bin_dir)?;
    let required_binary_count = pipeline_audit
        .iter()
        .filter(|row| row.role != "L33 atlas gate")
        .count();
    let blocking_hardcoded_binary_count = pipeline_audit.iter().filter(|row| row.blocking).count();
    let status = if !validation_failures.is_empty() {
        "MISSING_CELL_LOCAL_PACKETS_FAIL_L33_CONTRACT"
    } else {
        status_for(
            source_sha_fail_count,
            cells.len(),
            blocking_hardcoded_binary_count,
        )
    };
    let first_failed = if source_sha_fail_count > 0 {
        "L33 checksum failed".to_string()
    } else if let Some(first) = validation_failures.first() {
        first.clone()
    } else if blocking_hardcoded_binary_count > 0 {
        "required local residual-chain binaries are still hard-coded to subcell (6,4)".to_string()
    } else if cells.is_empty() {
        "none; L33 has no missing cells".to_string()
    } else {
        "none; pipeline is parameterized enough to start missing-cell packet generation".to_string()
    };

    let command_template = "./target/release/<parameterized-binary> --sub-i <i> --sub-j <j> --source <per-cell-source> --outdir <immutable-per-cell-artifact-dir>";
    let result = json!({
        "experiment_id": EXPERIMENT_ID,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "source_l33_experiment_id": L33_ID,
        "source_l33_path": cfg.l33.display().to_string(),
        "source_sha_status": [l33_sha],
        "source_sha_fail_count": source_sha_fail_count,
        "l33_validation_failures": validation_failures,
        "missing_cell_count": cells.len(),
        "missing_cells": cells,
        "expected_missing_cell_count": EXPECTED_MISSING_CELL_COUNT,
        "certified_cell_count_from_l33": l33["certified_cell_count"].clone(),
        "root_affine_skeleton_cell_count_from_l33": l33["root_affine_skeleton_cell_count"].clone(),
        "pipeline_audit": pipeline_audit,
        "required_binary_count": required_binary_count,
        "blocking_hardcoded_binary_count": blocking_hardcoded_binary_count,
        "packet_generation_started": false,
        "generated_local_packet_count": 0,
        "command_template_after_parameterization": command_template,
        "required_parameterization": [
            "add --sub-i and --sub-j to L21/L24/L27/L28/L29/L31/L32 binaries",
            "derive u0/u1 subcell intervals from the requested indices, not constants",
            "derive experiment IDs and output directories from the requested cell",
            "preserve per-cell source SHA verification",
            "make L31/L32 accept per-cell source artifact paths and expected subcell",
            "rerun L33 after generated packets exist"
        ],
        "first_failed_condition": first_failed,
        "what_this_rules_out": "It rules out starting a 63-cell batch by simply reusing hard-coded hard-cell binaries under different output names.",
        "what_this_does_not_rule_out": "It does not rule out the missing cells; it only shows the software pipeline must be parameterized before proof-facing batch generation can begin.",
        "next_dependency": "Parameterize the required residual-chain binaries, then rerun this L33A gate until it returns READY_TO_RUN.",
        "claim_ceiling": "Missing-cell local-packet generation gate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate."
    });

    let text = serde_json::to_string_pretty(&result)
        .map_err(|err| format!("failed to serialize result JSON: {err}"))?;
    fs::write(&result_path, text + "\n")
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
            "experiment_id": EXPERIMENT_ID,
            "status": result["status"],
            "missing_cell_count": result["missing_cell_count"],
            "blocking_hardcoded_binary_count": result["blocking_hardcoded_binary_count"],
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
