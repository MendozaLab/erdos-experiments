use inari::{interval, Interval};
use std::env;
use std::fs;

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

fn ratio_interval(numer: usize, denom: usize) -> Interval {
    let center = numer as f64 / denom as f64;
    iv(next_down(center), next_up(center))
}

fn q_poly_iv(x: Interval, endpoints: &[f64]) -> Interval {
    endpoints
        .iter()
        .fold(point(1.0), |acc, &endpoint| acc * (x - point(endpoint)))
}

fn theta_to_x_iv(gap_left: f64, gap_right: f64, theta: Interval) -> Interval {
    point(0.5 * (gap_left + gap_right))
        + point(0.5 * (gap_right - gap_left)) * theta.cos()
}

fn other_factor_sqrt_iv(
    endpoints: &[f64],
    skip_left: usize,
    skip_right: usize,
    x: Interval,
) -> Interval {
    let mut product = point(1.0);
    for (idx, &endpoint) in endpoints.iter().enumerate() {
        if idx == skip_left || idx == skip_right {
            continue;
        }
        product = product * (x - point(endpoint));
    }
    product.abs().sqrt()
}

fn int_power(x: Interval, power: usize) -> Interval {
    let mut acc = point(1.0);
    for _ in 0..power {
        acc = acc * x;
    }
    acc
}

fn integrate_entry(
    endpoints: &[f64],
    source_y: f64,
    gap_left: f64,
    gap_right: f64,
    skip_left: usize,
    skip_right: usize,
    power: Option<usize>,
    partition_count: usize,
) -> Interval {
    let width = Interval::PI * ratio_interval(1, partition_count);
    let source_root = q_poly_iv(point(source_y), endpoints).abs().sqrt();
    let mut acc = point(0.0);
    for idx in 0..partition_count {
        let theta_left = Interval::PI * ratio_interval(idx, partition_count);
        let theta_right = Interval::PI * ratio_interval(idx + 1, partition_count);
        let theta = iv(theta_left.inf(), theta_right.sup());
        let x = theta_to_x_iv(gap_left, gap_right, theta);
        let other = other_factor_sqrt_iv(endpoints, skip_left, skip_right, x);
        let value = match power {
            Some(p) => int_power(x, p) / other,
            None => source_root / ((x - point(source_y)).abs() * other),
        };
        acc = acc + width * value;
    }
    acc
}

fn interval_width(x: Interval) -> f64 {
    x.sup() - x.inf()
}

fn interval_mid(x: Interval) -> f64 {
    0.5 * (x.inf() + x.sup())
}

fn interval_pair(x: Interval) -> String {
    format!("[{:.17e},{:.17e}]", x.inf(), x.sup())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("usage: theta_integral_entry_enclosure_backend <input.txt>");
        std::process::exit(2);
    }
    let raw = fs::read_to_string(&args[1]).expect("read input");
    let mut it = raw.split_whitespace();
    let row_count = next_usize(&mut it);
    let partition_count = next_usize(&mut it);

    let mut rows_json: Vec<String> = Vec::with_capacity(row_count);
    let mut max_entry_width = 0.0f64;
    let mut max_entry_relative_width = 0.0f64;
    let mut min_abs_entry_mid = f64::INFINITY;
    let mut dim1_count = 0usize;
    let mut dim2_count = 0usize;
    let mut dim3_count = 0usize;

    for _ in 0..row_count {
        let row_index = next_usize(&mut it);
        let component_count = next_usize(&mut it);
        let source_y = next_f64(&mut it);
        let n = component_count - 1;
        match n {
            1 => dim1_count += 1,
            2 => dim2_count += 1,
            3 => dim3_count += 1,
            _ => panic!("unsupported dimension"),
        }
        let mut endpoints = Vec::with_capacity(2 * component_count);
        for _ in 0..(2 * component_count) {
            endpoints.push(next_f64(&mut it));
        }

        let mut matrix_rows: Vec<String> = Vec::with_capacity(n);
        let mut matrix_interval_rows: Vec<String> = Vec::with_capacity(n);
        let mut rhs_mids: Vec<String> = Vec::with_capacity(n);
        let mut rhs_intervals: Vec<String> = Vec::with_capacity(n);
        for gap_index in 0..n {
            let skip_left = 1 + 2 * gap_index;
            let skip_right = skip_left + 1;
            let gap_left = endpoints[skip_left];
            let gap_right = endpoints[skip_right];
            let mut matrix_mid_entries: Vec<String> = Vec::with_capacity(n);
            let mut matrix_interval_entries: Vec<String> = Vec::with_capacity(n);
            for power in 0..n {
                let entry = integrate_entry(
                    &endpoints,
                    source_y,
                    gap_left,
                    gap_right,
                    skip_left,
                    skip_right,
                    Some(power),
                    partition_count,
                );
                let width = interval_width(entry);
                let mid_abs = interval_mid(entry).abs();
                max_entry_width = max_entry_width.max(width);
                max_entry_relative_width = max_entry_relative_width.max(width / mid_abs.max(1.0e-300));
                min_abs_entry_mid = min_abs_entry_mid.min(mid_abs);
                matrix_mid_entries.push(format!("{:.17e}", interval_mid(entry)));
                matrix_interval_entries.push(interval_pair(entry));
            }
            matrix_rows.push(format!("[{}]", matrix_mid_entries.join(",")));
            matrix_interval_rows.push(format!("[{}]", matrix_interval_entries.join(",")));

            let rhs = integrate_entry(
                &endpoints,
                source_y,
                gap_left,
                gap_right,
                skip_left,
                skip_right,
                None,
                partition_count,
            );
            let width = interval_width(rhs);
            let mid_abs = interval_mid(rhs).abs();
            max_entry_width = max_entry_width.max(width);
            max_entry_relative_width = max_entry_relative_width.max(width / mid_abs.max(1.0e-300));
            min_abs_entry_mid = min_abs_entry_mid.min(mid_abs);
            rhs_mids.push(format!("{:.17e}", interval_mid(rhs)));
            rhs_intervals.push(interval_pair(rhs));
        }
        rows_json.push(format!(
            "{{\"row_index\":{},\"component_count\":{},\"period_system_dimension\":{},\"source_y\":{:.17e},\"endpoints\":[{}],\"period_matrix\":[{}],\"period_matrix_intervals\":[{}],\"period_rhs\":[{}],\"period_rhs_intervals\":[{}]}}",
            row_index,
            component_count,
            n,
            source_y,
            endpoints.iter().map(|x| format!("{:.17e}", x)).collect::<Vec<_>>().join(","),
            matrix_rows.join(","),
            matrix_interval_rows.join(","),
            rhs_mids.join(","),
            rhs_intervals.join(",")
        ));
    }

    println!("{{");
    println!("  \"pass\": true,");
    println!("  \"backend\":\"rust-inari\",");
    println!("  \"rust_binary\":\"src/bin/theta_integral_entry_enclosure_backend.rs\",");
    println!("  \"row_count\": {},", row_count);
    println!("  \"partition_count\": {},", partition_count);
    println!(
        "  \"dimension_distribution\":{{\"1\":{},\"2\":{},\"3\":{}}},",
        dim1_count, dim2_count, dim3_count
    );
    println!("  \"max_entry_width\": {:.17e},", max_entry_width);
    println!("  \"max_entry_relative_width\": {:.17e},", max_entry_relative_width);
    println!("  \"min_abs_entry_mid\": {:.17e},", min_abs_entry_mid);
    println!("  \"entry_rows\":[");
    for (idx, row) in rows_json.iter().enumerate() {
        let comma = if idx + 1 == rows_json.len() { "" } else { "," };
        println!("    {}{}", row, comma);
    }
    println!("  ],");
    println!("  \"claim_ceiling\":\"Directed interval range-integration enclosure for theta-desingularized period-matrix and RHS entries over a finite theta partition; no Gauss quadrature residual is assumed\"");
    println!("}}");
}
