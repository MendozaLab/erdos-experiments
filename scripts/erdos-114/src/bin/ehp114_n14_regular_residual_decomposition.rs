//! EHP #114 n=14 regular residual decomposition pilot.
//!
//! This consumes the L21 global critical-point diagnostic and attempts to
//! certify only the regular x-chart residual regions. It deliberately leaves
//! the critical-candidate regions untouched.
//!
//! This remains local n=14 hard-cell work. It is not a proof of Erdos #114,
//! not a global n=14 certificate, and not an exact lemniscate-length
//! certificate.

use ehp_n3_poc::ehp114_n14_cell::{ensure_source_subcell_matches, CellSpec};
use inari::{interval, Interval};
use serde::Serialize;
use serde_json::json;
use std::collections::HashSet;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const RAW_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-REGULAR-RESIDUAL-DECOMPOSITION-HARD-CELL-20260506-01";
const TAYLOR_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-REGULAR-SLICE-WALL-TAYLOR-HARD-CELL-20260506-01";
const TAYLOR_MONOTONE_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01";
const SOURCE_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01";
const DEGREE: usize = 14;
const EPS: f64 = 0.1;

#[derive(Clone, Copy, Debug)]
struct Complex {
    re: f64,
    im: f64,
}

#[derive(Clone, Copy, Debug)]
struct CInterval {
    re: Interval,
    im: Interval,
}

#[derive(Clone, Copy, Debug, Serialize)]
struct IntervalPair {
    lo: f64,
    hi: f64,
}

#[derive(Clone, Debug)]
struct RegularRegion {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    x: [f64; 2],
    y: [f64; 2],
}

#[derive(Clone, Debug, Serialize)]
struct RegionCertificate {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    ownership_key: String,
    chart_axis: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    status: String,
    reason: String,
    fx_interval: IntervalPair,
    fy_interval: IntervalPair,
    left_wall_f: IntervalPair,
    right_wall_f: IntervalPair,
    raw_left_wall_f: IntervalPair,
    raw_right_wall_f: IntervalPair,
    derivative_lower_bound: f64,
    derivative_ratio_upper: f64,
    wall_width_reduction: f64,
    length_upper: f64,
    adaptive_depth_used: usize,
    adaptive_leaf_count: usize,
    adaptive_closed_leaf_count: usize,
    adaptive_fail_count: usize,
    adaptive_children: Vec<RegionCertificate>,
}

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    out_dir: PathBuf,
    region_limit: usize,
    max_depth: usize,
    x_adaptive_depth: usize,
    wall_mode: String,
    experiment_id: Option<String>,
    cell: CellSpec,
}

fn iv(x: f64) -> Interval {
    interval!(x, x).unwrap()
}

fn iwrap(lo: f64, hi: f64) -> Interval {
    interval!(lo, hi).unwrap()
}

fn ipair(x: Interval) -> IntervalPair {
    IntervalPair {
        lo: x.inf(),
        hi: x.sup(),
    }
}

fn ci_one() -> CInterval {
    CInterval {
        re: iv(1.0),
        im: iv(0.0),
    }
}

fn ci_box(x: [f64; 2], y: [f64; 2]) -> CInterval {
    CInterval {
        re: iwrap(x[0], x[1]),
        im: iwrap(y[0], y[1]),
    }
}

fn ci_wall_x(x: f64, y: [f64; 2]) -> CInterval {
    CInterval {
        re: iv(x),
        im: iwrap(y[0], y[1]),
    }
}

fn ci_point(x: f64, y: f64) -> CInterval {
    CInterval {
        re: iv(x),
        im: iv(y),
    }
}

fn ci_sub(a: CInterval, b: CInterval) -> CInterval {
    CInterval {
        re: a.re - b.re,
        im: a.im - b.im,
    }
}

fn ci_mul(a: CInterval, b: CInterval) -> CInterval {
    CInterval {
        re: a.re * b.re - a.im * b.im,
        im: a.re * b.im + a.im * b.re,
    }
}

fn ci_conj(a: CInterval) -> CInterval {
    CInterval {
        re: a.re,
        im: -a.im,
    }
}

fn interval_square(x: Interval) -> Interval {
    let lo = x.inf();
    let hi = x.sup();
    if lo <= 0.0 && hi >= 0.0 {
        iwrap(0.0, lo.abs().max(hi.abs()).powi(2))
    } else {
        let a = lo * lo;
        let b = hi * hi;
        iwrap(a.min(b), a.max(b))
    }
}

fn interval_abs_lower(x: Interval) -> f64 {
    let lo = x.inf();
    let hi = x.sup();
    if lo <= 0.0 && hi >= 0.0 {
        0.0
    } else {
        lo.abs().min(hi.abs())
    }
}

fn interval_abs_upper(x: Interval) -> f64 {
    x.inf().abs().max(x.sup().abs())
}

fn interval_strictly_positive(x: Interval) -> bool {
    x.inf() > 0.0
}

fn interval_strictly_negative(x: Interval) -> bool {
    x.sup() < 0.0
}

fn abs_sq(z: CInterval) -> Interval {
    interval_square(z.re) + interval_square(z.im)
}

fn roots_of_unity(n: usize) -> Vec<Complex> {
    (0..n)
        .map(|j| {
            let theta = 2.0 * std::f64::consts::PI * (j as f64) / (n as f64);
            Complex {
                re: theta.cos(),
                im: theta.sin(),
            }
        })
        .collect()
}

fn base_roots(eps: f64, n: usize) -> Vec<Complex> {
    let radius = (1.0 - eps).powf(1.0 / (n as f64));
    roots_of_unity(n)
        .into_iter()
        .map(|z| Complex {
            re: radius * z.re,
            im: radius * z.im,
        })
        .collect()
}

fn u0_direction() -> Vec<Complex> {
    vec![
        Complex {
            re: 2.1464806571511146e-13,
            im: 0.12536063016562018,
        },
        Complex {
            re: -0.11029052423146467,
            im: 0.13183099723738242,
        },
        Complex {
            re: -0.2026155037603651,
            im: 0.044088538363229984,
        },
        Complex {
            re: -0.3356668203914393,
            im: -0.23859985068244477,
        },
        Complex {
            re: 0.33566682039232454,
            im: -0.23859985068243617,
        },
        Complex {
            re: 0.20261550376115514,
            im: 0.044088538363970045,
        },
        Complex {
            re: 0.11029052422993874,
            im: 0.13183099723676997,
        },
        Complex {
            re: 7.771839274878419e-14,
            im: 0.12536063016442625,
        },
        Complex {
            re: -0.11029052423164845,
            im: 0.13183099723560737,
        },
        Complex {
            re: -0.20261550375967857,
            im: 0.04408853836187157,
        },
        Complex {
            re: -0.33566682038773565,
            im: -0.23859985068065273,
        },
        Complex {
            re: 0.3356668203871522,
            im: -0.23859985068104475,
        },
        Complex {
            re: 0.20261550375874843,
            im: 0.04408853836095835,
        },
        Complex {
            re: 0.1102905242327204,
            im: 0.13183099723674271,
        },
    ]
}

fn u1_direction() -> Vec<Complex> {
    vec![
        Complex {
            re: -0.08603626747189994,
            im: 6.180051071480805e-13,
        },
        Complex {
            re: -0.09672845897773669,
            im: -0.05128639870562298,
        },
        Complex {
            re: 0.12294713305603927,
            im: -0.2280270846819635,
        },
        Complex {
            re: 0.3689060705436963,
            im: 0.17637503539601224,
        },
        Complex {
            re: -0.36890607054220126,
            im: 0.17637503539629037,
        },
        Complex {
            re: -0.12294713305686052,
            im: -0.22802708467920318,
        },
        Complex {
            re: 0.09672845897754945,
            im: -0.05128639870782282,
        },
        Complex {
            re: 0.0860362674705226,
            im: 3.4447185685221385e-13,
        },
        Complex {
            re: 0.09672845897619298,
            im: 0.0512863987081996,
        },
        Complex {
            re: -0.122947133057899,
            im: 0.22802708468114433,
        },
        Complex {
            re: -0.36890607054741453,
            im: -0.17637503539846167,
        },
        Complex {
            re: 0.3689060705470212,
            im: -0.17637503539817698,
        },
        Complex {
            re: 0.12294713305902331,
            im: 0.22802708468019145,
        },
        Complex {
            re: -0.09672845897603356,
            im: 0.05128639870845073,
        },
    ]
}

fn eps_shape_scale(eps: f64) -> f64 {
    eps.powf(1.0 / 28.0)
}

fn scale_interval(x: Interval, c: f64) -> Interval {
    x * iv(c)
}

fn complex_interval_linear(
    base: Complex,
    scale: f64,
    a: Interval,
    u0: Complex,
    b: Interval,
    u1: Complex,
) -> CInterval {
    CInterval {
        re: iv(base.re) + (scale_interval(a, u0.re) + scale_interval(b, u1.re)) * iv(scale),
        im: iv(base.im) + (scale_interval(a, u0.im) + scale_interval(b, u1.im)) * iv(scale),
    }
}

fn root_intervals(cell: CellSpec) -> Vec<CInterval> {
    let u0_pair = cell.u0_interval();
    let u1_pair = cell.u1_interval();
    let a = iwrap(u0_pair[0], u0_pair[1]);
    let b = iwrap(u1_pair[0], u1_pair[1]);
    let scale = eps_shape_scale(EPS);
    let base = base_roots(EPS, DEGREE);
    let u0 = u0_direction();
    let u1 = u1_direction();
    (0..DEGREE)
        .map(|i| complex_interval_linear(base[i], scale, a, u0[i], b, u1[i]))
        .collect()
}

fn eval_p_p1(z: CInterval, roots: &[CInterval]) -> (CInterval, CInterval) {
    let mut p = ci_one();
    let mut p1 = CInterval {
        re: iv(0.0),
        im: iv(0.0),
    };
    for root in roots {
        let factor = ci_sub(z, *root);
        let next_p1 = CInterval {
            re: p1.re * factor.re - p1.im * factor.im + p.re,
            im: p1.re * factor.im + p1.im * factor.re + p.im,
        };
        let next_p = ci_mul(p, factor);
        p = next_p;
        p1 = next_p1;
    }
    (p, p1)
}

fn gradient_intervals(p: CInterval, p1: CInterval) -> (Interval, Interval) {
    let q = ci_mul(p1, ci_conj(p));
    (q.re * iv(2.0), q.im * iv(-2.0))
}

fn f_fx_fy(z: CInterval, roots: &[CInterval]) -> (Interval, Interval, Interval) {
    let (p, p1) = eval_p_p1(z, roots);
    let f = abs_sq(p) - iv(1.0);
    let (fx, fy) = gradient_intervals(p, p1);
    (f, fx, fy)
}

fn f_only(z: CInterval, roots: &[CInterval]) -> Interval {
    let (p, _) = eval_p_p1(z, roots);
    abs_sq(p) - iv(1.0)
}

fn interval_width(x: Interval) -> f64 {
    x.sup() - x.inf()
}

fn wall_taylor_y(x: f64, y: [f64; 2], roots: &[CInterval]) -> Interval {
    let y_mid = 0.5 * (y[0] + y[1]);
    let y_radius = 0.5 * (y[1] - y[0]);
    let f_center = f_only(ci_point(x, y_mid), roots);
    let (_, _, fy_wall) = f_fx_fy(ci_wall_x(x, y), roots);
    f_center + fy_wall * iwrap(-y_radius, y_radius)
}

fn value_interval(value: &serde_json::Value, field: &str) -> Result<[f64; 2], String> {
    let arr = value
        .get(field)
        .and_then(|v| v.as_array())
        .ok_or_else(|| format!("missing interval field {field}"))?;
    if arr.len() != 2 {
        return Err(format!("interval field {field} does not have length 2"));
    }
    Ok([
        arr[0]
            .as_f64()
            .ok_or_else(|| format!("{field}[0] not numeric"))?,
        arr[1]
            .as_f64()
            .ok_or_else(|| format!("{field}[1] not numeric"))?,
    ])
}

fn parse_regular_region(
    idx: usize,
    value: &serde_json::Value,
) -> Result<Option<RegularRegion>, String> {
    if value
        .get("classification")
        .and_then(|v| v.as_str())
        .unwrap_or("")
        != "regular_region"
    {
        return Ok(None);
    }
    if value
        .get("dominant_axis")
        .and_then(|v| v.as_str())
        .unwrap_or("")
        != "x"
    {
        return Ok(None);
    }
    Ok(Some(RegularRegion {
        source_index: value
            .get("source_index")
            .and_then(|v| v.as_u64())
            .unwrap_or(idx as u64) as usize,
        source_ownership_key: value
            .get("source_ownership_key")
            .and_then(|v| v.as_str())
            .unwrap_or("missing-source-ownership-key")
            .to_string(),
        split_path: value
            .get("split_path")
            .and_then(|v| v.as_str())
            .unwrap_or("missing-split-path")
            .to_string(),
        x: value_interval(value, "x_interval")?,
        y: value_interval(value, "y_interval")?,
    }))
}

fn is_closed_status(status: &str) -> bool {
    status == "CERTIFIED_REGULAR_X_CHART" || status == "EXCLUDED_MONOTONE_X_CHART"
}

fn certify_regular_region_direct(
    region: &RegularRegion,
    roots: &[CInterval],
    wall_mode: &str,
) -> RegionCertificate {
    let (_, fx, fy) = f_fx_fy(ci_box(region.x, region.y), roots);
    let raw_left_f = f_only(ci_wall_x(region.x[0], region.y), roots);
    let raw_right_f = f_only(ci_wall_x(region.x[1], region.y), roots);
    let taylor_left_f = wall_taylor_y(region.x[0], region.y, roots);
    let taylor_right_f = wall_taylor_y(region.x[1], region.y, roots);
    let (left_f, right_f) = if wall_mode == "taylor-y" || wall_mode == "taylor-y-monotone" {
        (taylor_left_f, taylor_right_f)
    } else {
        (raw_left_f, raw_right_f)
    };
    let raw_width = interval_width(raw_left_f) + interval_width(raw_right_f);
    let chosen_width = interval_width(left_f) + interval_width(right_f);
    let derivative_lower_bound = interval_abs_lower(fx);
    let fy_abs_upper = interval_abs_upper(fy);
    let derivative_ratio_upper = if derivative_lower_bound > 0.0 {
        fy_abs_upper / derivative_lower_bound
    } else {
        f64::INFINITY
    };
    let width_y = region.y[1] - region.y[0];
    let length_upper = width_y * (1.0 + derivative_ratio_upper * derivative_ratio_upper).sqrt();
    let fx_sign_stable = interval_strictly_positive(fx) || interval_strictly_negative(fx);
    let wall_separated = (interval_strictly_positive(left_f)
        && interval_strictly_negative(right_f))
        || (interval_strictly_negative(left_f) && interval_strictly_positive(right_f));
    let monotone_excluded = if interval_strictly_negative(fx) {
        interval_strictly_positive(right_f) || interval_strictly_negative(left_f)
    } else if interval_strictly_positive(fx) {
        interval_strictly_positive(left_f) || interval_strictly_negative(right_f)
    } else {
        false
    };
    let (status, reason) = if !fx_sign_stable {
        (
            "FAIL_DERIVATIVE_SIGN".to_string(),
            "Fx is not sign-stable on the regular box".to_string(),
        )
    } else if wall_mode == "taylor-y-monotone" && monotone_excluded {
        (
            "EXCLUDED_MONOTONE_X_CHART".to_string(),
            "Taylor-y wall bound plus stable Fx excludes level-set crossing in this regular box"
                .to_string(),
        )
    } else if !wall_separated {
        (
            "FAIL_WALL_SEPARATION".to_string(),
            "x-wall F intervals do not have strict opposite signs".to_string(),
        )
    } else {
        (
            "CERTIFIED_REGULAR_X_CHART".to_string(),
            "Fx is sign-stable and x-wall F intervals are strictly separated".to_string(),
        )
    };

    let closed = is_closed_status(&status);
    RegionCertificate {
        source_index: region.source_index,
        source_ownership_key: region.source_ownership_key.clone(),
        split_path: region.split_path.clone(),
        ownership_key: format!(
            "{}:{}:x-chart",
            region.source_ownership_key, region.split_path
        ),
        chart_axis: "x=f(y)".to_string(),
        x_interval: region.x,
        y_interval: region.y,
        status,
        reason,
        fx_interval: ipair(fx),
        fy_interval: ipair(fy),
        left_wall_f: ipair(left_f),
        right_wall_f: ipair(right_f),
        raw_left_wall_f: ipair(raw_left_f),
        raw_right_wall_f: ipair(raw_right_f),
        derivative_lower_bound,
        derivative_ratio_upper,
        wall_width_reduction: raw_width - chosen_width,
        length_upper,
        adaptive_depth_used: 0,
        adaptive_leaf_count: 1,
        adaptive_closed_leaf_count: if closed { 1 } else { 0 },
        adaptive_fail_count: if closed { 0 } else { 1 },
        adaptive_children: Vec::new(),
    }
}

fn aggregate_adaptive_result(
    mut direct: RegionCertificate,
    children: Vec<RegionCertificate>,
    split_axis: &str,
) -> RegionCertificate {
    let adaptive_leaf_count = children.iter().map(|row| row.adaptive_leaf_count).sum();
    let adaptive_closed_leaf_count = children
        .iter()
        .map(|row| row.adaptive_closed_leaf_count)
        .sum();
    let adaptive_fail_count = children.iter().map(|row| row.adaptive_fail_count).sum();
    let adaptive_depth_used = 1 + children
        .iter()
        .map(|row| row.adaptive_depth_used)
        .max()
        .unwrap_or(0);
    let all_children_closed = children.iter().all(|row| is_closed_status(&row.status));
    let any_child_certified = children
        .iter()
        .any(|row| row.status == "CERTIFIED_REGULAR_X_CHART");

    if all_children_closed {
        direct.status = if any_child_certified {
            "CERTIFIED_REGULAR_X_CHART".to_string()
        } else {
            "EXCLUDED_MONOTONE_X_CHART".to_string()
        };
        direct.reason = format!(
            "adaptive {split_axis}-subdivision closed all {adaptive_leaf_count} leaf boxes for this original regular region"
        );
        direct.length_upper = children
            .iter()
            .filter(|row| row.status == "CERTIFIED_REGULAR_X_CHART")
            .map(|row| row.length_upper)
            .sum();
    } else {
        direct.reason = format!(
            "{}; adaptive {split_axis}-subdivision left {adaptive_fail_count} of {adaptive_leaf_count} leaf boxes unresolved",
            direct.reason
        );
    }

    direct.adaptive_depth_used = adaptive_depth_used;
    direct.adaptive_leaf_count = adaptive_leaf_count;
    direct.adaptive_closed_leaf_count = adaptive_closed_leaf_count;
    direct.adaptive_fail_count = adaptive_fail_count;
    direct.adaptive_children = children;
    direct
}

fn stable_fx_graph_envelope(
    mut direct: RegionCertificate,
    adaptive_context: &str,
) -> RegionCertificate {
    if direct.status != "FAIL_WALL_SEPARATION"
        || direct.derivative_lower_bound <= 0.0
        || !direct.derivative_ratio_upper.is_finite()
    {
        return direct;
    }

    direct.status = "CERTIFIED_REGULAR_X_CHART".to_string();
    direct.reason = format!(
        "{adaptive_context}; Fx is sign-stable, so any level-set portion inside this owned regular box is at most one x=f(y) graph over a subset of the y-interval; the full-y graph envelope is used as a conservative length upper bound without asserting endpoint crossing"
    );
    direct.adaptive_leaf_count = 1;
    direct.adaptive_closed_leaf_count = 1;
    direct.adaptive_fail_count = 0;
    direct.adaptive_children = Vec::new();
    direct
}

fn certify_regular_region(
    region: &RegularRegion,
    roots: &[CInterval],
    wall_mode: &str,
    max_depth: usize,
    x_adaptive_depth: usize,
    graph_envelope_fallback: bool,
) -> RegionCertificate {
    let direct = certify_regular_region_direct(region, roots, wall_mode);
    if is_closed_status(&direct.status) {
        return direct;
    }

    if max_depth > 0 {
        let y_mid = 0.5 * (region.y[0] + region.y[1]);
        if region.y[0] < y_mid && y_mid < region.y[1] {
            let left = RegularRegion {
                source_index: region.source_index,
                source_ownership_key: region.source_ownership_key.clone(),
                split_path: format!("{}/adaptive_y0", region.split_path),
                x: region.x,
                y: [region.y[0], y_mid],
            };
            let right = RegularRegion {
                source_index: region.source_index,
                source_ownership_key: region.source_ownership_key.clone(),
                split_path: format!("{}/adaptive_y1", region.split_path),
                x: region.x,
                y: [y_mid, region.y[1]],
            };
            let children = vec![
                certify_regular_region(
                    &left,
                    roots,
                    wall_mode,
                    max_depth - 1,
                    x_adaptive_depth,
                    graph_envelope_fallback,
                ),
                certify_regular_region(
                    &right,
                    roots,
                    wall_mode,
                    max_depth - 1,
                    x_adaptive_depth,
                    graph_envelope_fallback,
                ),
            ];
            let aggregated = aggregate_adaptive_result(direct.clone(), children, "y");
            if is_closed_status(&aggregated.status) || !graph_envelope_fallback {
                return aggregated;
            }
            return stable_fx_graph_envelope(
                direct,
                "adaptive y-subdivision did not produce strict wall separation on every leaf",
            );
        }
    }

    if x_adaptive_depth > 0 {
        let x_mid = 0.5 * (region.x[0] + region.x[1]);
        if region.x[0] < x_mid && x_mid < region.x[1] {
            let left = RegularRegion {
                source_index: region.source_index,
                source_ownership_key: region.source_ownership_key.clone(),
                split_path: format!("{}/adaptive_x0", region.split_path),
                x: [region.x[0], x_mid],
                y: region.y,
            };
            let right = RegularRegion {
                source_index: region.source_index,
                source_ownership_key: region.source_ownership_key.clone(),
                split_path: format!("{}/adaptive_x1", region.split_path),
                x: [x_mid, region.x[1]],
                y: region.y,
            };
            let children = vec![
                certify_regular_region(
                    &left,
                    roots,
                    wall_mode,
                    0,
                    x_adaptive_depth - 1,
                    graph_envelope_fallback,
                ),
                certify_regular_region(
                    &right,
                    roots,
                    wall_mode,
                    0,
                    x_adaptive_depth - 1,
                    graph_envelope_fallback,
                ),
            ];
            let aggregated = aggregate_adaptive_result(direct.clone(), children, "x");
            if is_closed_status(&aggregated.status) || !graph_envelope_fallback {
                return aggregated;
            }
            return stable_fx_graph_envelope(
                direct,
                "adaptive x-subdivision did not produce strict wall separation on every leaf",
            );
        }
    }

    if graph_envelope_fallback {
        return stable_fx_graph_envelope(
            direct,
            "strict wall separation was unavailable at the adaptive depth limit",
        );
    }

    direct
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01_RESULTS.json"),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-REGULAR-RESIDUAL-DECOMPOSITION-HARD-CELL-20260506-01"),
        region_limit: 0,
        max_depth: 0,
        x_adaptive_depth: 0,
        wall_mode: "raw".to_string(),
        experiment_id: None,
        cell: CellSpec::hard_cell(),
    };
    let mut sub_i = cfg.cell.sub_i;
    let mut sub_j = cfg.cell.sub_j;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--source" => {
                if i + 1 >= args.len() {
                    return Err("--source requires a path".to_string());
                }
                cfg.source = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--outdir" => {
                if i + 1 >= args.len() {
                    return Err("--outdir requires a path".to_string());
                }
                cfg.out_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--region-limit" => {
                if i + 1 >= args.len() {
                    return Err("--region-limit requires an integer".to_string());
                }
                cfg.region_limit = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --region-limit: {err}"))?;
                i += 2;
            }
            "--max-depth" => {
                if i + 1 >= args.len() {
                    return Err("--max-depth requires an integer".to_string());
                }
                cfg.max_depth = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --max-depth: {err}"))?;
                i += 2;
            }
            "--x-adaptive-depth" => {
                if i + 1 >= args.len() {
                    return Err("--x-adaptive-depth requires an integer".to_string());
                }
                cfg.x_adaptive_depth = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --x-adaptive-depth: {err}"))?;
                i += 2;
            }
            "--wall-mode" => {
                if i + 1 >= args.len() {
                    return Err(
                        "--wall-mode requires raw, taylor-y, or taylor-y-monotone".to_string()
                    );
                }
                cfg.wall_mode = args[i + 1].clone();
                if cfg.wall_mode != "raw"
                    && cfg.wall_mode != "taylor-y"
                    && cfg.wall_mode != "taylor-y-monotone"
                {
                    return Err(
                        "--wall-mode must be raw, taylor-y, or taylor-y-monotone".to_string()
                    );
                }
                i += 2;
            }
            "--experiment-id" => {
                if i + 1 >= args.len() {
                    return Err("--experiment-id requires a value".to_string());
                }
                cfg.experiment_id = Some(args[i + 1].clone());
                i += 2;
            }
            "--sub-i" => {
                if i + 1 >= args.len() {
                    return Err("--sub-i requires an integer".to_string());
                }
                sub_i = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --sub-i: {err}"))?;
                i += 2;
            }
            "--sub-j" => {
                if i + 1 >= args.len() {
                    return Err("--sub-j requires an integer".to_string());
                }
                sub_j = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --sub-j: {err}"))?;
                i += 2;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
    }
    cfg.cell = CellSpec::new(sub_i, sub_j)?;
    Ok(cfg)
}

fn experiment_id_for(wall_mode: &str) -> &'static str {
    match wall_mode {
        "taylor-y" => TAYLOR_EXPERIMENT_ID,
        "taylor-y-monotone" => TAYLOR_MONOTONE_EXPERIMENT_ID,
        _ => RAW_EXPERIMENT_ID,
    }
}

fn unix_timestamp_string() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs()
        .to_string()
}

fn sha256_file(path: &PathBuf) -> Result<String, String> {
    let bytes = fs::read(path)
        .map_err(|err| format!("failed to read {} for sha256: {err}", path.display()))?;
    Ok(sha256::digest(bytes))
}

fn first_failed_condition(
    certificates: &[RegionCertificate],
    ownership_duplicate_count: usize,
) -> String {
    if ownership_duplicate_count > 0 {
        return format!("ownership duplicate count is {ownership_duplicate_count}");
    }
    certificates
        .iter()
        .find(|row| {
            row.status != "CERTIFIED_REGULAR_X_CHART" && row.status != "EXCLUDED_MONOTONE_X_CHART"
        })
        .map(|row| {
            format!(
                "{}:{} failed because {}",
                row.source_ownership_key, row.split_path, row.reason
            )
        })
        .unwrap_or_else(|| "none".to_string())
}

fn status_for(
    processed: usize,
    closed: usize,
    ownership_duplicate_count: usize,
    total: f64,
    cap: f64,
) -> &'static str {
    if processed == 0 {
        "REGULAR_RESIDUAL_FAIL_NO_REGULAR_INPUT"
    } else if ownership_duplicate_count > 0 {
        "REGULAR_RESIDUAL_FAIL_OWNERSHIP"
    } else if closed < processed {
        "REGULAR_RESIDUAL_FAIL_WALL_SEPARATION"
    } else if total > cap {
        "REGULAR_RESIDUAL_FAIL_BUDGET"
    } else {
        "REGULAR_RESIDUAL_PASS_NOT_GLOBAL_PROOF"
    }
}

fn write_report(
    result: &serde_json::Value,
    path: &PathBuf,
    experiment_id: &str,
) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Regular Residual Decomposition\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Processed regular regions: `{}`\n\
- Certified regular regions: `{}`\n\
- Excluded regular regions: `{}`\n\
- Wall separation failures: `{}`\n\
- Ownership duplicates: `{}`\n\
- Adaptive max depth: `{}`\n\
- X-adaptive max depth: `{}`\n\
- Adaptive leaf boxes: `{}`\n\
- Adaptive unresolved leaves: `{}`\n\
- Regular slice length upper: `{}`\n\
- Total validated length upper: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
This diagnostic attempts the x-dominant regular regions found by the L21 global \
critical-point diagnostic. It does not process or promote the critical-candidate \
regions. When `max_depth` is positive, failed regular walls may be subdivided in \
the y direction. When `x_adaptive_depth` is positive, failed stable-`Fx` leaves may \
also be subdivided in the x direction; if strict wall separation still does not \
close, the stable-`Fx` graph-envelope theorem is allowed as a conservative length \
upper bound over the owned y-interval. Closure is still reported against the original \
owned source region, with child split paths preserved where subdivision is load-bearing. \
A pass means the regular slice has a theorem-shaped x-chart certificate with wall \
separation, monotone exclusion, or stable-`Fx` graph-envelope ownership; it is still \
not a global n=14 certificate.\n\n\
## Claim Ceiling\n\n\
{}\n",
        experiment_id,
        result["source_experiment_id"]
            .as_str()
            .unwrap_or(SOURCE_EXPERIMENT_ID),
        result["status"],
        result["processed_regular_region_count"],
        result["regular_regions_certified"],
        result["regular_regions_excluded"],
        result["wall_separation_fail_count"],
        result["ownership_duplicate_count"],
        result["parameters"]["max_depth"],
        result["parameters"]["x_adaptive_depth"],
        result["adaptive_leaf_count"],
        result["adaptive_fail_leaf_count"],
        result["regular_slice_length_upper"],
        result["total_validated_length_upper"],
        result["exact_length_cap"],
        result["margin_to_cap"],
        result["first_failed_condition"],
        result["claim_ceiling"]
            .as_str()
            .unwrap_or("local diagnostic only"),
    );
    fs::write(path, report)
        .map_err(|err| format!("failed to write report {}: {err}", path.display()))
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let cfg = parse_args()?;
    let experiment_id = cfg
        .experiment_id
        .clone()
        .unwrap_or_else(|| experiment_id_for(&cfg.wall_mode).to_string());
    fs::create_dir_all(&cfg.out_dir).map_err(|err| {
        format!(
            "failed to create output directory {}: {err}",
            cfg.out_dir.display()
        )
    })?;
    let result_path = cfg.out_dir.join(format!("{}_RESULTS.json", experiment_id));
    let report_path = cfg.out_dir.join(format!("{}_REPORT.md", experiment_id));
    let sha_path = cfg
        .out_dir
        .join(format!("{}_RESULTS.sha256", experiment_id));
    for path in [&result_path, &report_path, &sha_path] {
        if path.exists() {
            return Err(format!(
                "refusing to overwrite existing artifact: {}",
                path.display()
            ));
        }
    }

    let source_text = fs::read_to_string(&cfg.source)
        .map_err(|err| format!("failed to read source {}: {err}", cfg.source.display()))?;
    let source_json: serde_json::Value =
        serde_json::from_str(&source_text).map_err(|err| format!("bad source JSON: {err}"))?;
    let actual_source_experiment_id = source_json["experiment_id"]
        .as_str()
        .unwrap_or(SOURCE_EXPERIMENT_ID)
        .to_string();
    let source_subcell_contract = ensure_source_subcell_matches(&source_json, cfg.cell)?;
    let source_accepted_length_upper = source_json["total_validated_length_upper"]
        .as_f64()
        .ok_or_else(|| "source missing total_validated_length_upper".to_string())?;
    let cap = source_json["exact_length_cap"]
        .as_f64()
        .ok_or_else(|| "source missing exact_length_cap".to_string())?;
    let classifications_json = source_json["classifications"]
        .as_array()
        .ok_or_else(|| "source missing classifications".to_string())?;
    let mut regular_regions = Vec::<RegularRegion>::new();
    for (idx, value) in classifications_json.iter().enumerate() {
        if let Some(region) = parse_regular_region(idx, value)? {
            regular_regions.push(region);
        }
    }
    if cfg.region_limit > 0 {
        regular_regions.truncate(cfg.region_limit);
    }

    let roots = root_intervals(cfg.cell);
    let graph_envelope_fallback = cfg.x_adaptive_depth > 0;
    let certificates: Vec<RegionCertificate> = regular_regions
        .iter()
        .map(|region| {
            certify_regular_region(
                region,
                &roots,
                &cfg.wall_mode,
                cfg.max_depth,
                cfg.x_adaptive_depth,
                graph_envelope_fallback,
            )
        })
        .collect();
    let mut seen = HashSet::new();
    let ownership_duplicate_count = certificates
        .iter()
        .filter(|row| !seen.insert(row.ownership_key.clone()))
        .count();
    let regular_regions_certified = certificates
        .iter()
        .filter(|row| row.status == "CERTIFIED_REGULAR_X_CHART")
        .count();
    let regular_regions_excluded = certificates
        .iter()
        .filter(|row| row.status == "EXCLUDED_MONOTONE_X_CHART")
        .count();
    let derivative_sign_fail_count = certificates
        .iter()
        .filter(|row| row.status == "FAIL_DERIVATIVE_SIGN")
        .count();
    let wall_separation_fail_count = certificates
        .iter()
        .filter(|row| row.status == "FAIL_WALL_SEPARATION")
        .count();
    let regular_slice_length_upper: f64 = certificates
        .iter()
        .filter(|row| row.status == "CERTIFIED_REGULAR_X_CHART")
        .map(|row| row.length_upper)
        .sum();
    let dependency_width_reduction: f64 = certificates
        .iter()
        .map(|row| row.wall_width_reduction)
        .sum();
    let adaptive_leaf_count: usize = certificates.iter().map(|row| row.adaptive_leaf_count).sum();
    let adaptive_closed_leaf_count: usize = certificates
        .iter()
        .map(|row| row.adaptive_closed_leaf_count)
        .sum();
    let adaptive_fail_leaf_count: usize =
        certificates.iter().map(|row| row.adaptive_fail_count).sum();
    let adaptive_max_depth_used: usize = certificates
        .iter()
        .map(|row| row.adaptive_depth_used)
        .max()
        .unwrap_or(0);
    let total_validated_length_upper = source_accepted_length_upper + regular_slice_length_upper;
    let status = status_for(
        certificates.len(),
        regular_regions_certified + regular_regions_excluded,
        ownership_duplicate_count,
        total_validated_length_upper,
        cap,
    );
    let result = json!({
        "experiment_id": experiment_id,
        "timestamp_unix": unix_timestamp_string(),
        "source_experiment_id": actual_source_experiment_id,
        "source_results_path": cfg.source,
        "source_subcell_contract": source_subcell_contract,
        "status": status,
        "degree": DEGREE,
        "eps": EPS,
        "subcell": cfg.cell,
        "cell_tag": cfg.cell.tag(),
        "parameters": {
            "region_limit": cfg.region_limit,
            "max_depth": cfg.max_depth,
            "x_adaptive_depth": cfg.x_adaptive_depth,
            "stable_fx_graph_envelope_fallback": graph_envelope_fallback,
            "wall_mode": cfg.wall_mode,
            "source_region_field": "classifications",
            "source_filter": "classification == regular_region and dominant_axis == x",
            "chart_rule": "certify x=f(y) when Fx is sign-stable and the two x-wall F intervals have strict opposite signs; optional x-adaptive leaves retain the same source ownership with appended split paths",
            "length_rule": "width_y * sqrt(1 + sup(|Fy/Fx|)^2)",
            "critical_candidate_policy": "do not process or promote L21 critical-candidate regions"
        },
        "source_regular_region_count": regular_regions.len(),
        "processed_regular_region_count": certificates.len(),
        "regular_regions_certified": regular_regions_certified,
        "regular_regions_excluded": regular_regions_excluded,
        "derivative_sign_fail_count": derivative_sign_fail_count,
        "wall_separation_fail_count": wall_separation_fail_count,
        "wall_taylor_certified_count": if cfg.wall_mode == "taylor-y" || cfg.wall_mode == "taylor-y-monotone" { regular_regions_certified } else { 0 },
        "wall_taylor_excluded_count": if cfg.wall_mode == "taylor-y-monotone" { regular_regions_excluded } else { 0 },
        "wall_taylor_fail_count": if cfg.wall_mode == "taylor-y" || cfg.wall_mode == "taylor-y-monotone" { certificates.len().saturating_sub(regular_regions_certified + regular_regions_excluded) } else { 0 },
        "adaptive_leaf_count": adaptive_leaf_count,
        "adaptive_closed_leaf_count": adaptive_closed_leaf_count,
        "adaptive_fail_leaf_count": adaptive_fail_leaf_count,
        "adaptive_max_depth_used": adaptive_max_depth_used,
        "ownership_duplicate_count": ownership_duplicate_count,
        "dependency_width_reduction": dependency_width_reduction,
        "regular_slice_length_upper": regular_slice_length_upper,
        "source_accepted_length_upper": source_accepted_length_upper,
        "total_validated_length_upper": total_validated_length_upper,
        "exact_length_cap": cap,
        "margin_to_cap": cap - total_validated_length_upper,
        "first_failed_condition": first_failed_condition(&certificates, ownership_duplicate_count),
        "certificates": certificates,
        "claim_ceiling": "Local n=14 hard-cell regular residual decomposition diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.",
        "next_blocker": if regular_regions_certified + regular_regions_excluded == regular_regions.len() && !regular_regions.is_empty() {
            "The eight-region regular slice is closed by certified branches or monotone exclusions. Next proof-facing step is the separate 56-region critical-candidate theorem."
        } else if cfg.wall_mode == "taylor-y" {
            "Taylor-y wall evaluation did not certify all regular regions. The next proof-facing step is full affine arithmetic over root parameters or an analytic root-location bracket, not more raw wall evaluation."
        } else if cfg.wall_mode == "taylor-y-monotone" {
            "Taylor-y monotone exclusion did not close all regular regions. The next proof-facing step is full affine arithmetic over root parameters or an analytic root-location bracket."
        } else {
            "At least one regular L21 region failed x-wall sign separation or ownership; the next proof-facing step is a sharper wall theorem or a smaller analytic partition for the regular slice before returning to critical candidates."
        },
        "elapsed_secs": started.elapsed().as_secs_f64()
    });

    fs::write(
        &result_path,
        serde_json::to_string_pretty(&result).unwrap() + "\n",
    )
    .map_err(|err| {
        format!(
            "failed to write result JSON {}: {err}",
            result_path.display()
        )
    })?;
    write_report(&result, &report_path, &experiment_id)?;
    let sha = sha256_file(&result_path)?;
    fs::write(
        &sha_path,
        format!(
            "{sha}  {}\n",
            result_path.file_name().unwrap().to_string_lossy()
        ),
    )
    .map_err(|err| format!("failed to write sha file {}: {err}", sha_path.display()))?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "experiment_id": experiment_id,
            "status": status,
            "wall_mode": result["parameters"]["wall_mode"],
            "max_depth": result["parameters"]["max_depth"],
            "x_adaptive_depth": result["parameters"]["x_adaptive_depth"],
            "processed_regular_region_count": result["processed_regular_region_count"],
            "regular_regions_certified": result["regular_regions_certified"],
            "regular_regions_excluded": result["regular_regions_excluded"],
            "wall_taylor_certified_count": result["wall_taylor_certified_count"],
            "wall_taylor_excluded_count": result["wall_taylor_excluded_count"],
            "wall_taylor_fail_count": result["wall_taylor_fail_count"],
            "wall_separation_fail_count": result["wall_separation_fail_count"],
            "adaptive_leaf_count": result["adaptive_leaf_count"],
            "adaptive_fail_leaf_count": result["adaptive_fail_leaf_count"],
            "adaptive_max_depth_used": result["adaptive_max_depth_used"],
            "ownership_duplicate_count": result["ownership_duplicate_count"],
            "dependency_width_reduction": result["dependency_width_reduction"],
            "regular_slice_length_upper": result["regular_slice_length_upper"],
            "total_validated_length_upper": result["total_validated_length_upper"],
            "margin_to_cap": result["margin_to_cap"],
            "first_failed_condition": result["first_failed_condition"],
            "result": result_path,
            "report": report_path,
            "sha256": sha
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
