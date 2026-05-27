use inari::{interval, Interval};
use std::collections::BTreeMap;
use std::env;
use std::fs;

const CONDITION_THRESHOLD: f64 = 1.0e10;

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

fn status_json(name: &str, status: &str, detail: &str) -> String {
    format!(
        "\"{}\":{{\"status\":\"{}\",\"detail\":\"{}\"}}",
        name,
        status,
        json_escape(detail)
    )
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

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("usage: phi_k_stable_gap_period_basis_interval_certificate_backend <input.txt>");
        std::process::exit(2);
    }

    let raw = fs::read_to_string(&args[1]).expect("read input");
    let (fields, mut errors) = parse_input(&raw);

    let claimed_dimension = get_usize(&fields, "claimed_dimension", &mut errors);
    let object_count = get_usize(&fields, "object_count", &mut errors);
    let normalization_excluded = get_bool(&fields, "normalization_excluded", &mut errors);
    let sign_convention_match = get_bool(&fields, "sign_convention_match", &mut errors);
    let transform_condition_hi = get_f64(&fields, "transform_condition_hi", &mut errors);
    let recovered_direction_count = get_usize(&fields, "recovered_direction_count", &mut errors);
    let recovered_direction_independent_count = get_usize(
        &fields,
        "recovered_direction_independent_count",
        &mut errors,
    );
    let transformed_matrix_condition_hi =
        get_f64(&fields, "transformed_matrix_condition_hi", &mut errors);
    let smallest_singular_lo = get_f64(&fields, "smallest_singular_lo", &mut errors);

    let structure_pass = errors.is_empty()
        && claimed_dimension == Some(24)
        && object_count == Some(24)
        && normalization_excluded == Some(true)
        && sign_convention_match == Some(true);

    let transform_condition_iv = transform_condition_hi.map(point_enclosure);
    let transformed_matrix_condition_iv = transformed_matrix_condition_hi.map(point_enclosure);
    let smallest_singular_iv = smallest_singular_lo.map(point_enclosure);

    let (b1_status, b1_detail) = if !structure_pass {
        (
            "BLOCKED_PRE_AUDIT",
            "Structural pre-audit did not pass; transform condition is not consumed.",
        )
    } else if transform_condition_iv.map(|x| x.sup() > 0.0 && x.sup() < CONDITION_THRESHOLD)
        == Some(true)
    {
        (
            "PASS_CONTRACT_FIXTURE",
            "Transform condition fixture is below threshold.",
        )
    } else {
        (
            "FAIL_CONTRACT_FIXTURE",
            "Transform condition fixture is missing, nonpositive, or above threshold.",
        )
    };

    let (b2_status, b2_detail) = if !structure_pass {
        (
            "BLOCKED_PRE_AUDIT",
            "Structural pre-audit did not pass; recovered directions are not consumed.",
        )
    } else if recovered_direction_count == Some(11)
        && recovered_direction_independent_count == Some(11)
    {
        (
            "PASS_CONTRACT_FIXTURE",
            "All 11 recovered directions are marked interval-independent in the fixture.",
        )
    } else {
        (
            "FAIL_CONTRACT_FIXTURE",
            "Recovered direction count or independence witness count is not 11.",
        )
    };

    let b1_pass = b1_status == "PASS_CONTRACT_FIXTURE";
    let b2_pass = b2_status == "PASS_CONTRACT_FIXTURE";
    let (b3_status, b3_detail) = if !structure_pass || !b1_pass || !b2_pass {
        (
            "BLOCKED_PRECONDITION",
            "B3 is blocked until structural pre-audit, B1, and B2 pass.",
        )
    } else if transformed_matrix_condition_iv
        .map(|x| x.sup() > 0.0 && x.sup() < CONDITION_THRESHOLD)
        == Some(true)
        && smallest_singular_iv.map(|x| x.inf() > 0.0) == Some(true)
    {
        (
            "PASS_CONTRACT_FIXTURE",
            "Transformed matrix condition and smallest singular value fixtures pass.",
        )
    } else {
        (
            "FAIL_CONTRACT_FIXTURE",
            "Transformed matrix condition is above threshold or smallest singular value is not positive.",
        )
    };

    let overall_status = if !structure_pass {
        "STRUCTURAL_PRE_AUDIT_FAIL_CLOSED"
    } else if b1_status.starts_with("FAIL")
        || b2_status.starts_with("FAIL")
        || b3_status.starts_with("FAIL")
    {
        "FAIL_CLOSED_CONTRACT_FIXTURE_REJECTED"
    } else if b1_pass && b2_pass && b3_status == "PASS_CONTRACT_FIXTURE" {
        "PASS_CONTRACT_FIXTURE_ONLY_PRIVATE_NUMERIC_PAYLOAD_REQUIRED"
    } else {
        "BLOCKED_CONTRACT_FIXTURE"
    };

    println!("{{");
    println!("  \"backend\":\"phi_k_stable_gap_period_basis_interval_certificate_backend\",");
    println!("  \"backend_scope\":\"FAIL_CLOSED_CONTRACT_HARNESS_ONLY\",");
    println!("  \"input_path\":\"{}\",", json_escape(&args[1]));
    println!("  \"overall_status\":\"{}\",", overall_status);
    println!("  \"ordinary_f64_is_not_certificate\":true,");
    println!("  \"private_numeric_payload_required_for_real_certificate\":true,");
    println!(
        "  \"thresholds\":{{\"condition_number_hi\":{:.17e}}},",
        CONDITION_THRESHOLD
    );
    println!("  \"parse_errors\":{},", string_list_json(&errors));
    println!("  \"structural_pre_audit\":{{");
    println!(
        "    \"claimed_dimension\":{},",
        claimed_dimension.unwrap_or(0)
    );
    println!("    \"object_count\":{},", object_count.unwrap_or(0));
    println!(
        "    \"normalization_excluded\":{},",
        normalization_excluded.unwrap_or(false)
    );
    println!(
        "    \"sign_convention_match\":{},",
        sign_convention_match.unwrap_or(false)
    );
    println!(
        "    \"status\":\"{}\"",
        if structure_pass { "PASS" } else { "FAIL" }
    );
    println!("  }},");
    println!("  \"interval_enclosures\":{{");
    println!(
        "    \"transform_condition_hi\":[{:.17e},{:.17e}],",
        transform_condition_iv.map(|x| x.inf()).unwrap_or(f64::NAN),
        transform_condition_iv.map(|x| x.sup()).unwrap_or(f64::NAN)
    );
    println!(
        "    \"transformed_matrix_condition_hi\":[{:.17e},{:.17e}],",
        transformed_matrix_condition_iv
            .map(|x| x.inf())
            .unwrap_or(f64::NAN),
        transformed_matrix_condition_iv
            .map(|x| x.sup())
            .unwrap_or(f64::NAN)
    );
    println!(
        "    \"smallest_singular_lo\":[{:.17e},{:.17e}]",
        smallest_singular_iv.map(|x| x.inf()).unwrap_or(f64::NAN),
        smallest_singular_iv.map(|x| x.sup()).unwrap_or(f64::NAN)
    );
    println!("  }},");
    println!("  \"stage_statuses\":{{");
    println!(
        "    {},",
        status_json("B1_transform_interval_certificate", b1_status, b1_detail)
    );
    println!(
        "    {},",
        status_json("B2_recovered_direction_certificate", b2_status, b2_detail)
    );
    println!(
        "    {}",
        status_json("B3_transformed_matrix_certificate", b3_status, b3_detail)
    );
    println!("  }}");
    println!("}}");
}
