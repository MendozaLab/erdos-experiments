use inari::{interval, Interval};
use std::env;
use std::fs;

const ALPHA_LEFT: f64 = 0.921815564263329;
const ALPHA_RIGHT: f64 = 1.654239551094528;
const PI: f64 = std::f64::consts::PI;
const EDGE_TRIM_FRACTION: f64 = 1.0e-6;

fn iv(lo: f64, hi: f64) -> Interval {
    interval!(lo, hi).expect("valid interval")
}

fn point(x: f64) -> Interval {
    iv(x, x)
}

fn next_token<'a>(it: &mut impl Iterator<Item = &'a str>) -> &'a str {
    it.next().expect("missing token")
}

fn next_usize<'a>(it: &mut impl Iterator<Item = &'a str>) -> usize {
    next_token(it).parse::<usize>().expect("usize token")
}

fn next_f64<'a>(it: &mut impl Iterator<Item = &'a str>) -> f64 {
    next_token(it).parse::<f64>().expect("f64 token")
}

fn q_poly_iv(x: Interval, endpoints: &[f64]) -> Interval {
    endpoints
        .iter()
        .fold(point(1.0), |acc, &endpoint| acc * (x - point(endpoint)))
}

fn q_poly_f64(x: f64, endpoints: &[f64]) -> f64 {
    endpoints.iter().fold(1.0, |acc, &endpoint| acc * (x - endpoint))
}

fn poly_value_iv(coefficients: &[(f64, f64)], x: Interval) -> Interval {
    let mut acc = point(0.0);
    let mut power = point(1.0);
    for &(lo, hi) in coefficients {
        acc = acc + iv(lo, hi) * power;
        power = power * x;
    }
    acc
}

fn component_sign(component_index: usize, component_count: usize, side: &str) -> f64 {
    if side == "left" {
        if component_index % 2 == 0 {
            1.0
        } else {
            -1.0
        }
    } else if (component_count - 1 - component_index) % 2 == 0 {
        1.0
    } else {
        -1.0
    }
}

fn source_side(endpoints: &[f64], y: f64) -> &'static str {
    if y < endpoints[0] {
        "left"
    } else if y > endpoints[endpoints.len() - 1] {
        "right"
    } else {
        panic!("source is not exterior")
    }
}

fn corrected_density_iv(
    endpoints: &[f64],
    y: f64,
    coefficients: &[(f64, f64)],
    component_index: usize,
    x: Interval,
) -> Interval {
    let component_count = endpoints.len() / 2;
    let side = source_side(endpoints, y);
    let sign = point(component_sign(component_index, component_count, side));
    let ry = q_poly_f64(y, endpoints).abs().sqrt();
    let source_ratio = point(ry) / (x - point(y)).abs();
    let numerator = sign * (source_ratio - poly_value_iv(coefficients, x));
    numerator / (point(PI) * q_poly_iv(x, endpoints).abs().sqrt())
}

fn phi_k_density_iv(
    endpoints: &[f64],
    left_y: f64,
    right_y: f64,
    left_coefficients: &[(f64, f64)],
    right_coefficients: &[(f64, f64)],
    component_index: usize,
    x: Interval,
) -> Interval {
    point(ALPHA_LEFT)
        * corrected_density_iv(endpoints, left_y, left_coefficients, component_index, x)
        + point(ALPHA_RIGHT)
            * corrected_density_iv(endpoints, right_y, right_coefficients, component_index, x)
}

fn interval_json(name: &str, x: Interval) -> String {
    format!("\"{}\":[{:.17e},{:.17e}]", name, x.inf(), x.sup())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("usage: phi_k_density_panel_backend <input.txt>");
        std::process::exit(2);
    }
    let raw = fs::read_to_string(&args[1]).expect("read input");
    let mut it = raw.split_whitespace();
    let row_count = next_usize(&mut it);

    let mut pass_count = 0usize;
    let mut panel_count = 0usize;
    let mut negative_panel_count = 0usize;
    let mut min_density_lower = f64::INFINITY;
    let mut max_density_upper = f64::NEG_INFINITY;
    let mut min_mass_lower = f64::INFINITY;
    let mut max_mass_upper = f64::NEG_INFINITY;
    let mut worst_row_index = 0usize;
    let mut worst_component_index = 0usize;
    let mut worst_panel_index = 0usize;
    let mut failure_indices: Vec<usize> = Vec::new();

    for _ in 0..row_count {
        let row_index = next_usize(&mut it);
        let component_count = next_usize(&mut it);
        let panels_per_component = next_usize(&mut it);
        let coefficient_count = next_usize(&mut it);
        let left_y = next_f64(&mut it);
        let right_y = next_f64(&mut it);

        let mut endpoints = Vec::with_capacity(2 * component_count);
        for _ in 0..(2 * component_count) {
            endpoints.push(next_f64(&mut it));
        }
        let mut left_coefficients = Vec::with_capacity(coefficient_count);
        for _ in 0..coefficient_count {
            left_coefficients.push((next_f64(&mut it), next_f64(&mut it)));
        }
        let mut right_coefficients = Vec::with_capacity(coefficient_count);
        for _ in 0..coefficient_count {
            right_coefficients.push((next_f64(&mut it), next_f64(&mut it)));
        }

        let mut row_pass = true;
        let mut row_mass = point(0.0);
        for component_index in 0..component_count {
            let a = endpoints[2 * component_index];
            let b = endpoints[2 * component_index + 1];
            let dx = (b - a) / panels_per_component as f64;
            let trim = EDGE_TRIM_FRACTION * dx;
            for panel_index in 0..panels_per_component {
                let lo = a + panel_index as f64 * dx + trim;
                let hi = a + (panel_index + 1) as f64 * dx - trim;
                let density = phi_k_density_iv(
                    &endpoints,
                    left_y,
                    right_y,
                    &left_coefficients,
                    &right_coefficients,
                    component_index,
                    iv(lo, hi),
                );
                panel_count += 1;
                if density.inf() < min_density_lower {
                    min_density_lower = density.inf();
                    worst_row_index = row_index;
                    worst_component_index = component_index;
                    worst_panel_index = panel_index;
                }
                max_density_upper = max_density_upper.max(density.sup());
                if density.inf() <= 0.0 {
                    row_pass = false;
                    negative_panel_count += 1;
                }
                row_mass = row_mass + density * point(hi - lo);
            }
        }
        min_mass_lower = min_mass_lower.min(row_mass.inf());
        max_mass_upper = max_mass_upper.max(row_mass.sup());
        if row_pass {
            pass_count += 1;
        } else {
            failure_indices.push(row_index);
        }
    }

    let pass = pass_count == row_count;
    println!("{{");
    println!("  \"pass\": {},", pass);
    println!("  \"backend\":\"rust-inari\",");
    println!("  \"rust_binary\":\"src/bin/phi_k_density_panel_backend.rs\",");
    println!("  \"row_count\": {},", row_count);
    println!("  \"pass_count\": {},", pass_count);
    println!("  \"panel_count\": {},", panel_count);
    println!("  \"negative_panel_count\": {},", negative_panel_count);
    println!("  \"min_density_lower\": {:.17e},", min_density_lower);
    println!("  \"max_density_upper\": {:.17e},", max_density_upper);
    println!("  {},", interval_json("trimmed_proxy_mass_interval", iv(min_mass_lower, max_mass_upper)));
    println!("  \"worst_row_index\": {},", worst_row_index);
    println!("  \"worst_component_index\": {},", worst_component_index);
    println!("  \"worst_panel_index\": {},", worst_panel_index);
    println!("  \"edge_trim_fraction\": {:.17e},", EDGE_TRIM_FRACTION);
    println!("  \"failure_indices\":[");
    for (idx, row_index) in failure_indices.iter().enumerate() {
        let comma = if idx + 1 == failure_indices.len() { "" } else { "," };
        println!("    {}{}", row_index, comma);
    }
    println!("  ],");
    println!("  \"claim_ceiling\":\"Rust/Inari interval density-panel replay over coefficient boxes supplied by the numerical error-estimate packet; endpoint-edge and formal quadrature coefficient-box derivation remain pending\"");
    println!("}}");
}
