use inari::{interval, Interval};
use std::env;
use std::fs;

#[derive(Clone)]
struct Atom {
    _support_id: String,
    x: f64,
    weight: f64,
}

#[derive(Clone)]
struct Panel {
    panel_id: String,
    lo: f64,
    hi: f64,
}

fn iv(lo: f64, hi: f64) -> Interval {
    interval!(lo, hi).expect("valid interval")
}

fn point(x: f64) -> Interval {
    iv(x, x)
}

fn interval_pair(x: Interval) -> String {
    format!("[{:.17e},{:.17e}]", x.inf(), x.sup())
}

fn json_escape(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

fn read_input(path: &str) -> (Vec<Atom>, Vec<Atom>, Vec<Panel>, f64, f64) {
    let text = fs::read_to_string(path).expect("read input");
    let mut positive = Vec::new();
    let mut negative = Vec::new();
    let mut panels = Vec::new();
    let mut audit_x = 0.0;
    let mut audit_expected = 0.0;
    let mut section = "";

    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        match parts[0] {
            "POSITIVE" => {
                section = "POSITIVE";
                continue;
            }
            "NEGATIVE" => {
                section = "NEGATIVE";
                continue;
            }
            "PANELS" => {
                section = "PANELS";
                continue;
            }
            "AUDIT" => {
                audit_x = parts[1].parse::<f64>().expect("audit x");
                audit_expected = parts[2].parse::<f64>().expect("audit expected");
                continue;
            }
            _ => {}
        }

        if section == "POSITIVE" || section == "NEGATIVE" {
            let atom = Atom {
                _support_id: parts[0].to_string(),
                x: parts[1].parse::<f64>().expect("atom x"),
                weight: parts[2].parse::<f64>().expect("atom weight"),
            };
            if section == "POSITIVE" {
                positive.push(atom);
            } else {
                negative.push(atom);
            }
        } else if section == "PANELS" {
            panels.push(Panel {
                panel_id: parts[0].to_string(),
                lo: parts[1].parse::<f64>().expect("panel lo"),
                hi: parts[2].parse::<f64>().expect("panel hi"),
            });
        }
    }

    (positive, negative, panels, audit_x, audit_expected)
}

fn kernel_interval(x: Interval, support: f64) -> Interval {
    // Route convention: K(x, s) = -log(|x - s|).
    -((x - point(support)).abs().ln())
}

fn kernel_lower_for_positive(x: Interval, support: f64) -> f64 {
    let distance = (x - point(support)).abs();
    let farthest = distance.sup();
    if farthest <= 0.0 || !farthest.is_finite() {
        f64::NEG_INFINITY
    } else {
        -farthest.ln()
    }
}

fn kernel_upper_for_negative(x: Interval, support: f64) -> Option<f64> {
    let distance = (x - point(support)).abs();
    let nearest = distance.inf();
    if nearest <= 0.0 || !nearest.is_finite() {
        None
    } else {
        Some(-nearest.ln())
    }
}

fn signed_potential_interval(x: Interval, positive: &[Atom], negative: &[Atom]) -> Interval {
    let mut acc = point(0.0);
    for atom in positive {
        acc = acc + point(atom.weight) * kernel_interval(x, atom.x);
    }
    for atom in negative {
        acc = acc - point(atom.weight) * kernel_interval(x, atom.x);
    }
    acc
}

fn signed_potential_lower(x: Interval, positive: &[Atom], negative: &[Atom]) -> Option<f64> {
    let mut lower = 0.0;
    for atom in positive {
        let k_lower = kernel_lower_for_positive(x, atom.x);
        if !k_lower.is_finite() {
            return None;
        }
        lower += atom.weight * k_lower;
    }
    for atom in negative {
        let k_upper = kernel_upper_for_negative(x, atom.x)?;
        lower -= atom.weight * k_upper;
    }
    Some(lower)
}

fn signed_potential_point(x: f64, positive: &[Atom], negative: &[Atom]) -> f64 {
    let mut value = 0.0;
    for atom in positive {
        value += atom.weight * (-(x - atom.x).abs().ln());
    }
    for atom in negative {
        value -= atom.weight * (-(x - atom.x).abs().ln());
    }
    value
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("usage: period_corrected_signed_witness_replay_backend <input.txt>");
        std::process::exit(2);
    }

    let (positive, negative, panels, audit_x, audit_expected) = read_input(&args[1]);
    let audit_iv = signed_potential_interval(point(audit_x), &positive, &negative);
    let audit_f64 = signed_potential_point(audit_x, &positive, &negative);
    let audit_contains_expected = audit_iv.inf() <= audit_expected && audit_expected <= audit_iv.sup();
    let audit_contains_f64 = audit_iv.inf() <= audit_f64 && audit_f64 <= audit_iv.sup();
    let sign_audit_pass = audit_contains_expected && audit_contains_f64;

    let positive_mass: f64 = positive.iter().map(|a| a.weight).sum();
    let negative_mass: f64 = negative.iter().map(|a| a.weight).sum();
    let signed_mass = positive_mass - negative_mass;

    let mut rows = Vec::with_capacity(panels.len());
    let mut min_lower = f64::INFINITY;
    let mut weakest_panel = String::new();
    let mut fail_count = 0usize;
    let mut finite_count = 0usize;
    let mut positive_count = 0usize;

    for panel in panels.iter() {
        let x = iv(panel.lo, panel.hi);
        let lower = signed_potential_lower(x, &positive, &negative);
        let (lower_value, status) = match lower {
            Some(v) if v > 0.0 => {
                finite_count += 1;
                positive_count += 1;
                if v < min_lower {
                    min_lower = v;
                    weakest_panel = panel.panel_id.clone();
                }
                (v, "CERTIFIED_POSITIVE_F64_DIRECTED_BOUND")
            }
            Some(v) => {
                finite_count += 1;
                fail_count += 1;
                if v < min_lower {
                    min_lower = v;
                    weakest_panel = panel.panel_id.clone();
                }
                (v, "NONPOSITIVE_LOWER_BOUND")
            }
            None => {
                fail_count += 1;
                (f64::NEG_INFINITY, "SINGULAR_NEGATIVE_SUPPORT")
            }
        };
        rows.push(format!(
            "{{\"panel_id\":\"{}\",\"interval\":[{:.17e},{:.17e}],\"lower_bound\":{:.17e},\"status\":\"{}\"}}",
            json_escape(&panel.panel_id),
            panel.lo,
            panel.hi,
            lower_value,
            status
        ));
    }

    let epsilon = 1.0e-4;
    let pass = sign_audit_pass && fail_count == 0 && min_lower >= epsilon;

    println!("{{");
    println!("  \"backend\":\"rust-inari\",");
    println!("  \"rust_binary\":\"src/bin/period_corrected_signed_witness_replay_backend.rs\",");
    println!("  \"kernel\":\"-log_abs\",");
    println!("  \"positive_atom_count\":{},", positive.len());
    println!("  \"negative_atom_count\":{},", negative.len());
    println!("  \"positive_mass\":{:.17e},", positive_mass);
    println!("  \"negative_mass\":{:.17e},", negative_mass);
    println!("  \"signed_mass\":{:.17e},", signed_mass);
    println!("  \"epsilon_required\":{:.17e},", epsilon);
    println!("  \"sign_audit\":{{");
    println!("    \"audit_x\":{:.17e},", audit_x);
    println!("    \"expected_f64\":{:.17e},", audit_expected);
    println!("    \"rust_point_f64\":{:.17e},", audit_f64);
    println!("    \"rust_interval\":{},", interval_pair(audit_iv));
    println!("    \"contains_expected\":{},", audit_contains_expected);
    println!("    \"contains_rust_point_f64\":{},", audit_contains_f64);
    println!("    \"status\":\"{}\"", if sign_audit_pass { "PASS" } else { "FAIL" });
    println!("  }},");
    println!("  \"panel_summary\":{{");
    println!("    \"panel_count\":{},", panels.len());
    println!("    \"finite_count\":{},", finite_count);
    println!("    \"positive_count\":{},", positive_count);
    println!("    \"fail_count\":{},", fail_count);
    println!("    \"min_lower_bound\":{:.17e},", min_lower);
    println!("    \"weakest_panel_id\":\"{}\",", json_escape(&weakest_panel));
    println!("    \"status\":\"{}\"", if pass { "PASS" } else { "FAIL" });
    println!("  }},");
    println!("  \"rows\":[");
    for (idx, row) in rows.iter().enumerate() {
        let comma = if idx + 1 == rows.len() { "" } else { "," };
        println!("    {}{}", row, comma);
    }
    println!("  ],");
    println!("  \"claim_ceiling\":\"Rust/Inari signed witness fixed-cloud replay only; no KKT composition, global reduction, #1038 solution, or SOTA improvement\"");
    println!("}}");
}
