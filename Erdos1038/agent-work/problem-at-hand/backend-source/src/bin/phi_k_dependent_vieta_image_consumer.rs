// Fail-closed contract harness for the dependent Vieta image consumer scaffold.
//
// Scope: SCAFFOLD ONLY. This binary parses a key=value input file describing a
// synthetic dependent-Vieta certificate and emits a JSON status. It is not a
// theorem proof, not a real interval certificate, and not a substitute for a
// Rust/Inari directed interval audit of the Vieta image.

use inari::{interval, Interval};
use std::collections::BTreeMap;
use std::env;
use std::fs;

const ROOT_BOX_DIAMETER_HI_THRESHOLD: f64 = 1.0e-1;
const SCALED_VIETA_IMAGE_DIAMETER_HI_THRESHOLD: f64 = 1.0e-2;
const FIXED_CLOUD_DEVIATION_HI_THRESHOLD: f64 = 1.0e-2;

fn iv(lo: f64, hi: f64) -> Interval {
    interval!(lo, hi).expect("valid interval")
}

fn next_down(x: f64) -> f64 {
    if x == 0.0 {
        -f64::MIN_POSITIVE
    } else if x.is_sign_positive() {
        f64::from_bits(x.to_bits() - 1)
    } else {
        f64::from_bits(x.to_bits() + 1)
    }
}

fn next_up(x: f64) -> f64 {
    if x == 0.0 {
        f64::MIN_POSITIVE
    } else if x.is_sign_positive() {
        f64::from_bits(x.to_bits() + 1)
    } else {
        f64::from_bits(x.to_bits() - 1)
    }
}

fn point_enclosure(x: f64) -> Interval {
    iv(next_down(x), next_up(x))
}

fn json_escape(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

fn get_usize(
    fields: &BTreeMap<String, String>,
    key: &str,
    errors: &mut Vec<String>,
) -> Option<usize> {
    match fields.get(key) {
        Some(raw) => match raw.parse::<usize>() {
            Ok(value) => Some(value),
            Err(_) => {
                errors.push(format!("{} is not a usize", key));
                None
            }
        },
        None => {
            errors.push(format!("missing {}", key));
            None
        }
    }
}

fn get_f64(fields: &BTreeMap<String, String>, key: &str, errors: &mut Vec<String>) -> Option<f64> {
    match fields.get(key) {
        Some(raw) => match raw.parse::<f64>() {
            Ok(value) if value.is_finite() => Some(value),
            Ok(_) => {
                errors.push(format!("{} is not finite", key));
                None
            }
            Err(_) => {
                errors.push(format!("{} is not an f64", key));
                None
            }
        },
        None => {
            errors.push(format!("missing {}", key));
            None
        }
    }
}

fn get_bool(
    fields: &BTreeMap<String, String>,
    key: &str,
    errors: &mut Vec<String>,
) -> Option<bool> {
    match fields.get(key).map(String::as_str) {
        Some("true") | Some("1") => Some(true),
        Some("false") | Some("0") => Some(false),
        Some(_) => {
            errors.push(format!("{} is not a bool", key));
            None
        }
        None => {
            errors.push(format!("missing {}", key));
            None
        }
    }
}

fn parse_input(raw: &str) -> (BTreeMap<String, String>, Vec<String>) {
    let mut fields = BTreeMap::new();
    let mut errors = Vec::new();
    for (line_index, raw_line) in raw.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some((key, value)) = line.split_once('=') else {
            errors.push(format!("line {} is not key=value", line_index + 1));
            continue;
        };
        fields.insert(key.trim().to_string(), value.trim().to_string());
    }
    (fields, errors)
}

fn string_list_json(values: &[String]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| format!("\"{}\"", json_escape(value)))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn status_json(name: &str, status: &str, detail: &str) -> String {
    format!(
        "\"{}\":{{\"status\":\"{}\",\"detail\":\"{}\"}}",
        name,
        status,
        json_escape(detail)
    )
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("usage: phi_k_dependent_vieta_image_consumer <input.txt>");
        std::process::exit(2);
    }

    let raw = fs::read_to_string(&args[1]).expect("read input");
    let (fields, mut errors) = parse_input(&raw);

    let claimed_root_count = get_usize(&fields, "claimed_root_count", &mut errors);
    let multiplicity_sum = get_usize(&fields, "multiplicity_sum", &mut errors);
    let scaling_convention_match = get_bool(&fields, "scaling_convention_match", &mut errors);
    let dependency_relation_recorded =
        get_bool(&fields, "dependency_relation_recorded", &mut errors);

    let root_box_diameter_hi = get_f64(&fields, "root_box_diameter_hi", &mut errors);
    let scaled_vieta_image_diameter_hi =
        get_f64(&fields, "scaled_vieta_image_diameter_hi", &mut errors);
    let vieta_image_dimension = get_usize(&fields, "vieta_image_dimension", &mut errors);
    let vieta_image_independent_count =
        get_usize(&fields, "vieta_image_independent_count", &mut errors);
    let fixed_cloud_anchor_count = get_usize(&fields, "fixed_cloud_anchor_count", &mut errors);
    let fixed_cloud_deviation_hi = get_f64(&fields, "fixed_cloud_deviation_hi", &mut errors);

    let structure_pass = errors.is_empty()
        && claimed_root_count.is_some()
        && multiplicity_sum.is_some()
        && claimed_root_count == multiplicity_sum
        && claimed_root_count.unwrap_or(0) > 0
        && scaling_convention_match == Some(true)
        && dependency_relation_recorded == Some(true);

    let root_box_iv = root_box_diameter_hi.map(point_enclosure);
    let scaled_vieta_image_iv = scaled_vieta_image_diameter_hi.map(point_enclosure);
    let fixed_cloud_deviation_iv = fixed_cloud_deviation_hi.map(point_enclosure);

    let (v1_status, v1_detail) = if !structure_pass {
        (
            "BLOCKED_PRE_AUDIT",
            "Structural pre-audit did not pass; root-box enclosure is not consumed.",
        )
    } else if root_box_iv.map(|x| x.sup() > 0.0 && x.sup() < ROOT_BOX_DIAMETER_HI_THRESHOLD)
        == Some(true)
    {
        (
            "PASS_CONTRACT_FIXTURE",
            "Root-box diameter fixture is below threshold.",
        )
    } else {
        (
            "FAIL_CONTRACT_FIXTURE",
            "Root-box diameter fixture is missing, nonpositive, or above threshold.",
        )
    };

    let (v2_status, v2_detail) = if !structure_pass {
        (
            "BLOCKED_PRE_AUDIT",
            "Structural pre-audit did not pass; scaled Vieta image is not consumed.",
        )
    } else if scaled_vieta_image_iv
        .map(|x| x.sup() > 0.0 && x.sup() < SCALED_VIETA_IMAGE_DIAMETER_HI_THRESHOLD)
        == Some(true)
        && vieta_image_dimension.is_some()
        && vieta_image_independent_count == vieta_image_dimension
    {
        (
            "PASS_CONTRACT_FIXTURE",
            "Scaled Vieta image diameter and independence witness counts pass the fixture contract.",
        )
    } else {
        (
            "FAIL_CONTRACT_FIXTURE",
            "Scaled Vieta image diameter is above threshold or independence witness count does not match dimension.",
        )
    };

    let v1_pass = v1_status == "PASS_CONTRACT_FIXTURE";
    let v2_pass = v2_status == "PASS_CONTRACT_FIXTURE";
    let (v3_status, v3_detail) = if !structure_pass || !v1_pass || !v2_pass {
        (
            "BLOCKED_PRECONDITION",
            "V3 is blocked until structural pre-audit, V1, and V2 pass.",
        )
    } else if fixed_cloud_anchor_count.map(|n| n > 0) == Some(true)
        && fixed_cloud_deviation_iv
            .map(|x| x.sup() >= 0.0 && x.sup() < FIXED_CLOUD_DEVIATION_HI_THRESHOLD)
            == Some(true)
    {
        (
            "PASS_CONTRACT_FIXTURE",
            "Fixed-cloud anchor count is positive and deviation fixture is within threshold.",
        )
    } else {
        (
            "FAIL_CONTRACT_FIXTURE",
            "Fixed-cloud anchor count is not positive or deviation fixture is above threshold.",
        )
    };

    let overall_status = if !structure_pass {
        "STRUCTURAL_PRE_AUDIT_FAIL_CLOSED"
    } else if v1_status.starts_with("FAIL")
        || v2_status.starts_with("FAIL")
        || v3_status.starts_with("FAIL")
    {
        "FAIL_CLOSED_CONTRACT_FIXTURE_REJECTED"
    } else if v1_pass && v2_pass && v3_status == "PASS_CONTRACT_FIXTURE" {
        "PASS_CONTRACT_FIXTURE_ONLY_PRIVATE_NUMERIC_PAYLOAD_REQUIRED"
    } else {
        "BLOCKED_CONTRACT_FIXTURE"
    };

    println!("{{");
    println!("  \"backend\":\"phi_k_dependent_vieta_image_consumer\",");
    println!("  \"backend_scope\":\"FAIL_CLOSED_CONTRACT_HARNESS_ONLY\",");
    println!("  \"input_path\":\"{}\",", json_escape(&args[1]));
    println!("  \"overall_status\":\"{}\",", overall_status);
    println!("  \"ordinary_f64_is_not_certificate\":true,");
    println!("  \"private_numeric_payload_required_for_real_certificate\":true,");
    println!("  \"claim_ceiling\":\"SCAFFOLD_ONLY__NO_DEPENDENT_VIETA_THEOREM_PASS\",");
    println!(
        "  \"thresholds\":{{\"root_box_diameter_hi\":{:.17e},\"scaled_vieta_image_diameter_hi\":{:.17e},\"fixed_cloud_deviation_hi\":{:.17e}}},",
        ROOT_BOX_DIAMETER_HI_THRESHOLD,
        SCALED_VIETA_IMAGE_DIAMETER_HI_THRESHOLD,
        FIXED_CLOUD_DEVIATION_HI_THRESHOLD
    );
    println!("  \"parse_errors\":{},", string_list_json(&errors));
    println!("  \"structural_pre_audit\":{{");
    println!(
        "    \"claimed_root_count\":{},",
        claimed_root_count.unwrap_or(0)
    );
    println!("    \"multiplicity_sum\":{},", multiplicity_sum.unwrap_or(0));
    println!(
        "    \"scaling_convention_match\":{},",
        scaling_convention_match.unwrap_or(false)
    );
    println!(
        "    \"dependency_relation_recorded\":{},",
        dependency_relation_recorded.unwrap_or(false)
    );
    println!(
        "    \"status\":\"{}\"",
        if structure_pass { "PASS" } else { "FAIL" }
    );
    println!("  }},");
    println!("  \"interval_enclosures\":{{");
    println!(
        "    \"root_box_diameter_hi\":[{:.17e},{:.17e}],",
        root_box_iv.map(|x| x.inf()).unwrap_or(f64::NAN),
        root_box_iv.map(|x| x.sup()).unwrap_or(f64::NAN)
    );
    println!(
        "    \"scaled_vieta_image_diameter_hi\":[{:.17e},{:.17e}],",
        scaled_vieta_image_iv.map(|x| x.inf()).unwrap_or(f64::NAN),
        scaled_vieta_image_iv.map(|x| x.sup()).unwrap_or(f64::NAN)
    );
    println!(
        "    \"fixed_cloud_deviation_hi\":[{:.17e},{:.17e}]",
        fixed_cloud_deviation_iv
            .map(|x| x.inf())
            .unwrap_or(f64::NAN),
        fixed_cloud_deviation_iv
            .map(|x| x.sup())
            .unwrap_or(f64::NAN)
    );
    println!("  }},");
    println!(
        "  \"vieta_image_shape\":{{\"dimension\":{},\"independent_count\":{}}},",
        vieta_image_dimension.unwrap_or(0),
        vieta_image_independent_count.unwrap_or(0)
    );
    println!("  \"stage_statuses\":{{");
    println!(
        "    {},",
        status_json("V1_root_box_enclosure", v1_status, v1_detail)
    );
    println!(
        "    {},",
        status_json("V2_scaled_vieta_image_enclosure", v2_status, v2_detail)
    );
    println!(
        "    {}",
        status_json("V3_fixed_cloud_certificate", v3_status, v3_detail)
    );
    println!("  }}");
    println!("}}");
}
