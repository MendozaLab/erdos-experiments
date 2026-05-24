//! EHP #114 n=14 radial Puiseux compact-middle interval check.
//!
//! This binary certifies only the compact radial interval
//! `1e-4 <= eps <= 1e-1` for the one-parameter family `p_a(z)=z^14-a`,
//! using IEEE-1788 interval arithmetic through `inari`.
//!
//! It does not prove EHP #114. It is a local certificate component for the
//! radial Puiseux closure lane: the singular tail `0 < eps <= 1e-4` still
//! needs a hypergeometric connection-formula proof.

use inari::{interval, Interval};
use serde_json::json;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

const EXPERIMENT_ID: &str = "EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01";
const DEGREE: usize = 14;
const C14: f64 = 24.0;
const LSTAR_LOWER_N14: f64 = 30.852910841548532;

fn iv(x: f64) -> Interval {
    interval!(x, x).unwrap()
}

fn radial_integrand_bound(eps: Interval, u: Interval) -> Interval {
    // For t = pi +/- u:
    // |(1-eps) + exp(i t)|^2 = eps^2 + 2(1-eps)(1-cos u).
    // The length integrand for n=14 is base^(-13/28).
    let one = iv(1.0);
    let two = iv(2.0);
    let exponent = interval!(-13.0 / 28.0, -13.0 / 28.0).unwrap();
    let a = one - eps;
    let one_minus_cos = one - u.cos();
    let base = eps * eps + two * a * one_minus_cos;
    base.pow(exponent)
}

fn radial_length_upper(eps_lo: f64, eps_hi: f64, u_steps: usize) -> f64 {
    let eps = interval!(eps_lo, eps_hi).unwrap();
    let pi_hi = Interval::PI.sup();
    let du = pi_hi / (u_steps as f64);
    let du_iv = interval!(du, du).unwrap();
    let mut total = interval!(0.0, 0.0).unwrap();

    for i in 0..u_steps {
        let lo = (i as f64) * du;
        let hi = ((i + 1) as f64) * du;
        let u = interval!(lo, hi).unwrap();
        total = total + radial_integrand_bound(eps, u) * du_iv;
    }

    // Symmetry doubles integral over u in [0, pi].
    (iv(2.0) * total).sup()
}

fn rhs_upper(eps_hi: f64) -> f64 {
    let exponent = interval!(1.0 / (DEGREE as f64), 1.0 / (DEGREE as f64)).unwrap();
    (iv(C14) * iv(eps_hi).pow(exponent)).sup()
}

fn check_bin(eps_lo: f64, eps_hi: f64, u_steps: usize) -> serde_json::Value {
    let started = Instant::now();
    let length_upper = radial_length_upper(eps_lo, eps_hi, u_steps);
    let deficit_lower = LSTAR_LOWER_N14 - length_upper;
    let rhs = rhs_upper(eps_hi);
    let margin = deficit_lower - rhs;
    json!({
        "eps_lo": eps_lo,
        "eps_hi": eps_hi,
        "u_steps": u_steps,
        "length_upper": length_upper,
        "deficit_lower": deficit_lower,
        "rhs_upper_C14_eps_power": rhs,
        "margin": margin,
        "pass": margin > 0.0,
        "elapsed_secs": started.elapsed().as_secs_f64()
    })
}

fn out_dir_from_args() -> PathBuf {
    let args: Vec<String> = env::args().collect();
    let mut out_dir = PathBuf::from(".");
    let mut i = 1;
    while i < args.len() {
        if args[i] == "--outdir" && i + 1 < args.len() {
            out_dir = PathBuf::from(&args[i + 1]);
            i += 2;
        } else {
            i += 1;
        }
    }
    out_dir
}

fn main() {
    let started = Instant::now();
    let out_dir = out_dir_from_args();
    fs::create_dir_all(&out_dir).expect("failed to create output directory");

    let bins = [
        (1.0e-4, 2.0e-4, 250_000usize),
        (2.0e-4, 5.0e-4, 160_000usize),
        (5.0e-4, 1.0e-3, 120_000usize),
        (1.0e-3, 2.0e-3, 90_000usize),
        (2.0e-3, 5.0e-3, 70_000usize),
        (5.0e-3, 1.0e-2, 50_000usize),
        (1.0e-2, 2.0e-2, 40_000usize),
        (2.0e-2, 5.0e-2, 30_000usize),
        (5.0e-2, 1.0e-1, 24_000usize),
    ];

    let rows: Vec<_> = bins
        .iter()
        .map(|(lo, hi, steps)| check_bin(*lo, *hi, *steps))
        .collect();
    let all_pass = rows
        .iter()
        .all(|row| row["pass"].as_bool().unwrap_or(false));
    let min_margin = rows
        .iter()
        .map(|row| row["margin"].as_f64().unwrap())
        .fold(f64::INFINITY, f64::min);

    let result = json!({
        "experiment_id": EXPERIMENT_ID,
        "problem": "Erdos #114 / EHP lemniscate perimeter",
        "scope": "n=14 radial Puiseux compact-middle interval component",
        "status": if all_pass { "COMPACT_MIDDLE_CERTIFIED" } else { "COMPACT_MIDDLE_FAILED" },
        "claim_ceiling": "This is a shadow signature, not universal law. It certifies one radial-family compact interval only, not EHP #114.",
        "degree": DEGREE,
        "epsilon_domain": "[1e-4, 1e-1]",
        "missing_domain": "(0, 1e-4] singular tail; requires hypergeometric connection formula",
        "constant_C14": C14,
        "lstar_lower_n14": LSTAR_LOWER_N14,
        "integral_formula": "L_14(1-eps)=2*int_0^pi (eps^2 + 2*(1-eps)*(1-cos u))^(-13/28) du",
        "method": "IEEE-1788 interval upper Riemann enclosure with inari interval cos/pow",
        "bins": rows,
        "min_margin": min_margin,
        "all_bins_pass": all_pass,
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "next_blocker": "Prove the singular tail 0 < eps <= 1e-4 by Gauss 2F1/Puiseux connection formula."
    });

    let result_path = out_dir.join(format!("{EXPERIMENT_ID}_RESULTS.json"));
    fs::write(
        &result_path,
        serde_json::to_string_pretty(&result).unwrap() + "\n",
    )
    .expect("failed to write result JSON");
    println!("{}", serde_json::to_string_pretty(&result).unwrap());
}
