//! EHP #114 bridge packets from the finite n<15 certificate to Tao's large-n theorem.
//!
//! This emitter does not compute n=15 and does not extract Tao's threshold.
//! It builds the next bridge-control layer: a threshold extraction skeleton,
//! an effective-threshold checker with opaque rows, and a bridge decision
//! packet that blocks full synthesis until the threshold is numeric.

use serde_json::{json, Value};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const L37_ID: &str = "EXP-MATH-EHP114-TAO-THRESHOLD-EXTRACTION-SKELETON-20260507-02";
const L38_ID: &str = "EXP-MATH-EHP114-TAO-EFFECTIVE-THRESHOLD-CHECKER-20260507-02";
const L39_ID: &str = "EXP-MATH-EHP114-BRIDGE-DECISION-PACKET-20260507-02";

const MIDDLE_ID: &str = "EXP-MATH-EHP114-TAO-MIDDLE-KINGDOM-20260505-01";
const FINITE_ID: &str = "EXP-MATH-EHP114-FINITE-N-LESS-15-PROOF-PACKET-20260507-01";
const EXPECTED_TAO_SOURCE_SHA: &str =
    "ab792bd8fe806a985e44f692ccb8305cf386677509846979bb8ff35a276814a6";
const TAO_ARXIV_URL: &str = "https://arxiv.org/abs/2512.12455";
const ERDOS114_URL: &str = "https://www.erdosproblems.com/114";
const FINITE_FRONTIER_EXCLUSIVE: usize = 15;
const PRACTICAL_BOUND_DEFAULT: usize = 20;

#[derive(Clone, Debug)]
struct Config {
    out_root: PathBuf,
    middle_packet: PathBuf,
    finite_packet: PathBuf,
    practical_bound: usize,
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
        middle_packet: PathBuf::from(format!("../../Erdos114/{MIDDLE_ID}_RESULTS.json")),
        finite_packet: PathBuf::from(format!(
            "../../Erdos114/proof_path/{FINITE_ID}/{FINITE_ID}_RESULTS.json"
        )),
        practical_bound: PRACTICAL_BOUND_DEFAULT,
    };
    let args: Vec<String> = env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--out-root" => {
                cfg.out_root = next_path(&args, i, "--out-root")?;
                i += 2;
            }
            "--middle-packet" => {
                cfg.middle_packet = next_path(&args, i, "--middle-packet")?;
                i += 2;
            }
            "--finite-packet" => {
                cfg.finite_packet = next_path(&args, i, "--finite-packet")?;
                i += 2;
            }
            "--practical-bound" => {
                cfg.practical_bound = next_usize(&args, i, "--practical-bound")?;
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

fn next_usize(args: &[String], i: usize, flag: &str) -> Result<usize, String> {
    if i + 1 >= args.len() {
        return Err(format!("{flag} requires a value"));
    }
    args[i + 1]
        .parse::<usize>()
        .map_err(|err| format!("failed to parse {flag}: {err}"))
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

fn verify_result_sha(role: &str, path: &Path) -> Result<Value, String> {
    let sha_path = sha_path_for(path)?;
    let sha_text = fs::read_to_string(&sha_path)
        .map_err(|err| format!("failed to read {}: {err}", sha_path.display()))?;
    let expected = sha_text
        .split_whitespace()
        .next()
        .ok_or_else(|| format!("empty sha file {}", sha_path.display()))?
        .to_string();
    let actual = sha256_file(path)?;
    Ok(json!({
        "role": role,
        "path": path.display().to_string(),
        "sha_path": sha_path.display().to_string(),
        "expected_sha256": expected,
        "actual_sha256": actual,
        "sha_status": if expected == actual { "PASS" } else { "FAIL" }
    }))
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

fn bool_field(v: &Value, key: &str) -> Result<bool, String> {
    v.get(key)
        .and_then(Value::as_bool)
        .ok_or_else(|| format!("missing bool field {key}"))
}

fn value_usize(v: &Value) -> Option<usize> {
    v.as_u64().map(|n| n as usize)
}

fn artifact_paths(out_root: &Path, experiment_id: &str) -> (PathBuf, PathBuf, PathBuf) {
    let dir = out_root.join(experiment_id);
    (
        dir.join(format!("{experiment_id}_RESULTS.json")),
        dir.join(format!("{experiment_id}_REPORT.md")),
        dir.join(format!("{experiment_id}_RESULTS.sha256")),
    )
}

fn classify_dependency(phrase: &str, line: usize) -> (String, String, String) {
    let label = match line {
        1..=220 => "main_thm_and_normalization",
        221..=624 => "stokes_and_reference_value",
        625..=999 => "pre_geometry_constants",
        1000..=1240 => "geomcontrol",
        1241..=1300 => "x2_lem",
        1301..=1345 => "x3_lem",
        1346..=1390 => "lemni",
        1391..=1450 => "inside",
        1451..=1530 => "annulus_or_outside",
        1531..=1680 => "ets",
        1681..=1735 => "pots",
        1736..=1760 => "inside_2",
        1761..=1790 => "annulus_2",
        1791..=1965 => "outside_again",
        _ => "out_of_range",
    };
    let dependency_kind = match phrase {
        "sufficiently large" | "large enough" => "large_n_gate",
        "small enough" => "small_parameter_gate",
        "absolute constant" => "constant_dependency",
        "effectively computable" => "effectivity_notice",
        "no attempt to optimize" => "non_optimized_constant_notice",
        _ => "source_needle",
    };
    let priority = if matches!(
        label,
        "inside_2" | "annulus_2" | "outside_again" | "pots" | "ets"
    ) {
        "HIGH"
    } else if matches!(
        label,
        "geomcontrol" | "x2_lem" | "x3_lem" | "inside" | "annulus_or_outside"
    ) {
        "MEDIUM"
    } else {
        "LOW"
    };
    (
        label.to_string(),
        dependency_kind.to_string(),
        priority.to_string(),
    )
}

fn needle_dependency_rows(middle: &Value) -> Result<Vec<Value>, String> {
    let needle_lines = middle
        .get("needle_lines")
        .and_then(Value::as_object)
        .ok_or_else(|| "middle packet missing needle_lines object".to_string())?;
    let mut rows = Vec::new();
    for (phrase, lines) in needle_lines {
        let arr = lines
            .as_array()
            .ok_or_else(|| format!("needle_lines.{phrase} is not an array"))?;
        for line_value in arr {
            let line = value_usize(line_value)
                .ok_or_else(|| format!("needle line for {phrase} is not an integer"))?;
            let (label, kind, priority) = classify_dependency(phrase, line);
            rows.push(json!({
                "dependency_id": format!("needle-{}-{}", phrase.replace(' ', "-"), line),
                "source": "tao_arxiv_source",
                "source_line": line,
                "phrase": phrase,
                "classified_label": label,
                "dependency_kind": kind,
                "priority": priority,
                "classification_status": if label == "out_of_range" { "UNCLASSIFIED" } else { "CLASSIFIED" },
                "inequality_needed": "not_extracted_in_skeleton",
                "constants_required": "not_extracted_in_skeleton",
                "extraction_status": "OPAQUE_CONSTANT_DEPENDENCY"
            }));
        }
    }
    rows.sort_by_key(|row| row["source_line"].as_u64().unwrap_or(0));
    Ok(rows)
}

fn final_deficit_rows(middle: &Value) -> Result<Vec<Value>, String> {
    let landmarks = middle
        .get("landmarks")
        .and_then(Value::as_array)
        .ok_or_else(|| "middle packet missing landmarks array".to_string())?;
    let wanted = ["inside-2", "annulus-2", "outside-again"];
    let mut rows = Vec::new();
    for label in wanted {
        let item = landmarks
            .iter()
            .find(|row| row["label"].as_str() == Some(label))
            .ok_or_else(|| format!("missing final deficit landmark: {label}"))?;
        rows.push(json!({
            "dependency_id": format!("final-deficit-{label}"),
            "source": "tao_control_landmark",
            "source_label": label,
            "source_lines": [item["line_start"].clone(), item["line_end"].clone()],
            "role": item["role"].clone(),
            "middle_gap_use": item["middle_gap_use"].clone(),
            "priority": "CRITICAL",
            "classification_status": "CLASSIFIED",
            "inequality_needed": "explicit deficit bound strong enough to dominate accumulated errors",
            "constants_required": "not_extracted_in_skeleton",
            "extraction_status": "OPAQUE_CONSTANT_DEPENDENCY"
        }));
    }
    Ok(rows)
}

fn build_l37(cfg: &Config, started: Instant) -> Result<Artifact, String> {
    let middle = read_json(&cfg.middle_packet)?;
    let finite = read_json(&cfg.finite_packet)?;
    let middle_sha = verify_result_sha("tao_threshold_control_packet", &cfg.middle_packet)?;
    let finite_sha = verify_result_sha("finite_n_less15_packet", &cfg.finite_packet)?;
    let source_sha = str_field(&middle, "source_sha256")?;
    let source_sha_status = if source_sha == EXPECTED_TAO_SOURCE_SHA {
        "PASS_PINNED_SOURCE_SHA_MATCH"
    } else {
        "FAIL_PINNED_SOURCE_SHA_MISMATCH"
    };

    let needle_rows = needle_dependency_rows(&middle)?;
    let final_rows = final_deficit_rows(&middle)?;
    let all_rows = needle_rows
        .iter()
        .chain(final_rows.iter())
        .cloned()
        .collect::<Vec<_>>();
    let unclassified_count = all_rows
        .iter()
        .filter(|row| row["classification_status"] != "CLASSIFIED")
        .count();
    let middle_landmark_count = middle
        .get("landmarks")
        .and_then(Value::as_array)
        .map(Vec::len)
        .unwrap_or(0);
    let status = if str_field(&middle, "status")? != "CONTROL_PACKET_READY" {
        "THRESHOLD_EXTRACTION_BLOCKED_SOURCE_PARSE"
    } else if source_sha_status != "PASS_PINNED_SOURCE_SHA_MATCH" {
        "THRESHOLD_EXTRACTION_BLOCKED_SOURCE_PARSE"
    } else if unclassified_count > 0 {
        "THRESHOLD_EXTRACTION_BLOCKED_UNLABELED_DEPENDENCY"
    } else {
        "THRESHOLD_EXTRACTION_SKELETON_READY"
    };

    let result = json!({
        "experiment_id": L37_ID,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "source_artifacts": [
            {
                "artifact_id": MIDDLE_ID,
                "path": cfg.middle_packet.display().to_string(),
                "status": middle["status"].clone(),
                "sha_check": middle_sha
            },
            {
                "artifact_id": FINITE_ID,
                "path": cfg.finite_packet.display().to_string(),
                "status": finite["status"].clone(),
                "sha_check": finite_sha
            }
        ],
        "arxiv_source": {
            "url": TAO_ARXIV_URL,
            "source_sha256": source_sha,
            "expected_source_sha256": EXPECTED_TAO_SOURCE_SHA,
            "source_sha_status": source_sha_status,
            "source_line_count": middle["source_line_count"].clone()
        },
        "erdos114_problem_url": ERDOS114_URL,
        "finite_n_less15_status": finite["status"].clone(),
        "finite_n_less15_pass": finite["all_n_less_15_pass"].clone(),
        "tao_control_landmark_count": middle_landmark_count,
        "dependency_count": all_rows.len(),
        "needle_dependency_count": needle_rows.len(),
        "final_deficit_dependency_count": final_rows.len(),
        "classified_dependency_count": all_rows.len() - unclassified_count,
        "unclassified_dependency_count": unclassified_count,
        "dependency_rows": all_rows,
        "represented_final_deficit_labels": ["inside-2", "annulus-2", "outside-again"],
        "threshold_extraction_scope": "skeleton only; no explicit constants or N_i values extracted",
        "first_failed_condition": if status == "THRESHOLD_EXTRACTION_SKELETON_READY" { "none" } else { "source parse or dependency classification blocker" },
        "claim_ceiling": "Threshold extraction skeleton only. This packet classifies Tao-side dependencies but does not extract a numerical threshold and does not prove the full EHP114 statement."
    });

    let report = format!(
        "# EHP114 L37 Tao Threshold Extraction Skeleton\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
Status: `{}`.\n\n\
This packet consumes the pinned Tao threshold-control source map and the finite \
`n<15` proof packet. It classifies every recorded Tao source occurrence of \
`sufficiently large`, `large enough`, `small enough`, `absolute constant`, and \
the final deficit landmarks `inside-2`, `annulus-2`, and `outside-again`.\n\n\
## Meaning\n\n\
The bridge is now organized as a dependency table. No numerical threshold is \
claimed here; the next packet must turn these rows into explicit inequalities \
or name the first opaque constant blocker.\n",
        L37_ID, status
    );
    Ok(Artifact {
        experiment_id: L37_ID,
        result,
        report,
    })
}

fn build_l38(l37: &Value, started: Instant) -> Result<Artifact, String> {
    let dependency_rows = l37
        .get("dependency_rows")
        .and_then(Value::as_array)
        .ok_or_else(|| "L37 result missing dependency_rows".to_string())?;
    let mut checker_rows = Vec::new();
    for dep in dependency_rows {
        checker_rows.push(json!({
            "dependency_id": dep["dependency_id"].clone(),
            "source_line": dep.get("source_line").cloned().unwrap_or(Value::Null),
            "source_label": dep.get("source_label").cloned().unwrap_or(Value::Null),
            "classified_label": dep.get("classified_label").cloned().unwrap_or(Value::Null),
            "priority": dep["priority"].clone(),
            "inequality_needed": dep["inequality_needed"].clone(),
            "constants_required": dep["constants_required"].clone(),
            "extracted_N_i": Value::Null,
            "explicit_status": "OPAQUE",
            "blocker": "No explicit numerical constants extracted from Tao proof in this checker pass."
        }));
    }
    let dependency_count = checker_rows.len();
    let explicit_dependency_count = checker_rows
        .iter()
        .filter(|row| !row["extracted_N_i"].is_null())
        .count();
    let opaque_dependency_count = dependency_count - explicit_dependency_count;
    let status = if opaque_dependency_count == 0 {
        "TAO_THRESHOLD_EXTRACTED_OVERLAPS_SMALL_N"
    } else {
        "TAO_THRESHOLD_REMAINS_OPAQUE"
    };
    let first_blocker = checker_rows
        .iter()
        .find(|row| row["explicit_status"] == "OPAQUE")
        .map(|row| {
            row["dependency_id"]
                .as_str()
                .unwrap_or("unknown")
                .to_string()
        })
        .unwrap_or_else(|| "none".to_string());
    let result = json!({
        "experiment_id": L38_ID,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "source_experiment_id": L37_ID,
        "source_status": l37["status"].clone(),
        "dependency_count": dependency_count,
        "explicit_dependency_count": explicit_dependency_count,
        "opaque_dependency_count": opaque_dependency_count,
        "candidate_N0": Value::Null,
        "candidate_N0_rule": "candidate_N0 may be emitted only when every dependency row has an explicit extracted_N_i",
        "threshold_rows": checker_rows,
        "first_blocker": first_blocker,
        "first_failed_condition": if status == "TAO_THRESHOLD_REMAINS_OPAQUE" { "Tao threshold remains opaque; no candidate_N0 emitted" } else { "none" },
        "claim_ceiling": "Effective-threshold checker skeleton only. No numerical Tao threshold is extracted and no finite gap is authorized."
    });
    let report = format!(
        "# EHP114 L38 Effective Threshold Inequality Checker\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
Status: `{}`.\n\n\
The checker has `{}` dependency rows and `{}` opaque rows. Because at least one \
dependency is opaque, `candidate_N0` is deliberately null. The packet therefore \
does not authorize any `n=15` computation or full synthesis.\n",
        L38_ID, status, dependency_count, opaque_dependency_count
    );
    Ok(Artifact {
        experiment_id: L38_ID,
        result,
        report,
    })
}

fn build_l39(
    cfg: &Config,
    finite: &Value,
    l38: &Value,
    started: Instant,
) -> Result<Artifact, String> {
    let finite_pass = bool_field(finite, "all_n_less_15_pass")?;
    let threshold_status = str_field(l38, "status")?;
    let candidate_n0 = l38.get("candidate_N0").cloned().unwrap_or(Value::Null);
    let status = if !finite_pass {
        "BRIDGE_DECISION_BLOCKED_FINITE_PACKET"
    } else if threshold_status == "TAO_THRESHOLD_EXTRACTED_OVERLAPS_SMALL_N" {
        "BRIDGE_DECISION_FULL_SYNTHESIS_ALLOWED"
    } else if threshold_status == "TAO_THRESHOLD_EXTRACTED_FINITE_GAP" {
        "BRIDGE_DECISION_FINITE_GAP_REQUIRES_CERTIFICATES"
    } else {
        "BRIDGE_DECISION_OPAQUE_THRESHOLD_ANALYTIC_TIGHTENING_REQUIRED"
    };
    let finite_gap_degrees = if let Some(n0) = candidate_n0.as_u64() {
        if (FINITE_FRONTIER_EXCLUSIVE as u64) < n0 && n0 <= cfg.practical_bound as u64 {
            ((FINITE_FRONTIER_EXCLUSIVE as u64)..n0)
                .map(|n| json!(n))
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };
    let n15_cluster_run_authorized = status == "BRIDGE_DECISION_FINITE_GAP_REQUIRES_CERTIFICATES"
        && !finite_gap_degrees.is_empty();
    let result = json!({
        "experiment_id": L39_ID,
        "status": status,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "finite_n_less15_artifact_id": FINITE_ID,
        "finite_n_less15_status": finite["status"].clone(),
        "finite_n_less15_pass": finite_pass,
        "tao_threshold_checker_artifact_id": L38_ID,
        "tao_threshold_status": threshold_status,
        "candidate_N0": candidate_n0,
        "practical_bound": cfg.practical_bound,
        "finite_frontier_exclusive": FINITE_FRONTIER_EXCLUSIVE,
        "finite_gap_degrees_to_certify": finite_gap_degrees,
        "n15_or_higher_computation_authorized": n15_cluster_run_authorized,
        "full_proof_synthesis_allowed": status == "BRIDGE_DECISION_FULL_SYNTHESIS_ALLOWED",
        "next_action": if status == "BRIDGE_DECISION_OPAQUE_THRESHOLD_ANALYTIC_TIGHTENING_REQUIRED" {
            "Start analytic constant-tightening on inside-2, annulus-2, outside-again, then back-propagate through pots, ets, inside, annulus, outside, and geomcontrol."
        } else if status == "BRIDGE_DECISION_FINITE_GAP_REQUIRES_CERTIFICATES" {
            "Certify the listed finite gap degrees before full synthesis."
        } else if status == "BRIDGE_DECISION_FULL_SYNTHESIS_ALLOWED" {
            "Emit L40 full proof synthesis candidate."
        } else {
            "Repair finite n<15 packet before bridge work."
        },
        "l40_full_synthesis_emitted": false,
        "l40_reason": "No L40 packet is emitted because Tao threshold remains opaque unless L38 later extracts an overlapping threshold.",
        "first_failed_condition": if threshold_status == "TAO_THRESHOLD_REMAINS_OPAQUE" { "Tao threshold remains opaque" } else { "none" },
        "claim_ceiling": "Bridge decision packet only. It combines the finite n<15 packet with the Tao threshold checker, but does not prove the full all-degree EHP114 statement."
    });
    let report = format!(
        "# EHP114 L39 Bridge Decision Packet\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
Status: `{}`.\n\n\
The finite side passes for `1 <= n < 15`, but the Tao threshold checker remains \
opaque. Therefore no `n=15` or higher computation is authorized by this packet, \
and no full proof synthesis packet is emitted.\n\n\
## Next Action\n\n\
Start analytic constant-tightening on `inside-2`, `annulus-2`, and \
`outside-again`, then back-propagate through `pots`, `ets`, `inside`, \
`annulus`, `outside`, and `geomcontrol`.\n",
        L39_ID, status
    );
    Ok(Artifact {
        experiment_id: L39_ID,
        result,
        report,
    })
}

fn check_targets(out_root: &Path, artifacts: &[Artifact]) -> Result<(), String> {
    for artifact in artifacts {
        let (results, report, sha) = artifact_paths(out_root, artifact.experiment_id);
        for path in [results, report, sha] {
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
    let l37 = build_l37(&cfg, started)?;
    let l38 = build_l38(&l37.result, started)?;
    let finite = read_json(&cfg.finite_packet)?;
    let l39 = build_l39(&cfg, &finite, &l38.result, started)?;
    let artifacts = vec![l37, l38, l39];
    check_targets(&cfg.out_root, &artifacts)?;

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
            "status": "EHP114_TAO_THRESHOLD_BRIDGE_PACKETS_EMITTED",
            "artifact_count": emitted.len(),
            "l40_emitted": false,
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
