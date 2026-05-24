//! EHP #114 n=14 third-order root-collar remainder pilot.
//!
//! This diagnostic consumes the normal-collar failure artifact and re-tests its
//! 64 unresolved regions with a Taylor collar inequality. The new ingredient is
//! a recurrence for p, p', p'', and p''' so the collar can log an explicit
//! third-directional remainder bound.
//!
//! This remains local n=14 hard-cell work. It is not a proof of Erdos #114,
//! not a global n=14 certificate, and not an exact lemniscate-length
//! certificate.

use ehp_n3_poc::ehp114_n14_cell::{ensure_source_subcell_matches, CellSpec};
use inari::{interval, Interval};
use serde::Serialize;
use serde_json::json;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const DEFAULT_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01";
const SOURCE_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-HARD-CELL-20260506-01";
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

#[derive(Clone, Copy, Debug)]
struct HessianIntervals {
    fxx: Interval,
    fxy: Interval,
    fyy: Interval,
}

#[derive(Clone, Copy, Debug, Serialize)]
struct IntervalPair {
    lo: f64,
    hi: f64,
}

#[derive(Clone, Copy, Debug, Serialize)]
struct Direction {
    x: f64,
    y: f64,
}

#[derive(Clone, Debug)]
struct SourceRegion {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    x: [f64; 2],
    y: [f64; 2],
    source_reason: String,
    source_normal: Option<Direction>,
    source_tangent: Option<Direction>,
}

#[derive(Clone, Debug, Serialize)]
struct InequalityAudit {
    center_fn_abs_lower: f64,
    variation_bound_fn: f64,
    critical_exclusion_margin: f64,
    center_strip_f_abs_upper: f64,
    normal_remainder_bound: f64,
    wall_lhs: f64,
    wall_rhs: f64,
    wall_margin: f64,
    third_directional_upper: f64,
    f_nn_abs_upper: f64,
    f_nt_abs_upper: f64,
    f_tt_abs_upper: f64,
    normal_radius_upper: f64,
    normal_wall_distance_lower: f64,
    tangent_radius_upper: f64,
}

#[derive(Clone, Debug, Serialize)]
struct CertifiedThirdOrderCollar {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    normal: Direction,
    tangent: Direction,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    normal_interval: [f64; 2],
    tangent_interval: [f64; 2],
    f_interval: IntervalPair,
    fn_center_interval: IntervalPair,
    fn_full_first_order_interval: IntervalPair,
    ft_full_first_order_interval: IntervalPair,
    slope_abs_upper: f64,
    length_upper: f64,
    inequality: InequalityAudit,
}

#[derive(Clone, Debug, Serialize)]
struct ExcludedRegion {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    f_interval: IntervalPair,
}

#[derive(Clone, Debug, Serialize)]
struct UnresolvedThirdOrderRegion {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    source_reason: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    normal: Option<Direction>,
    tangent: Option<Direction>,
    normal_interval: Option<[f64; 2]>,
    tangent_interval: Option<[f64; 2]>,
    reason: String,
    failing_inequality: String,
    f_interval: IntervalPair,
    fn_center_interval: Option<IntervalPair>,
    fn_full_first_order_interval: Option<IntervalPair>,
    ft_full_first_order_interval: Option<IntervalPair>,
    inequality: Option<InequalityAudit>,
}

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    out_dir: PathBuf,
    experiment_id: String,
    cell: CellSpec,
    region_limit: usize,
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

fn ci_zero() -> CInterval {
    CInterval {
        re: iv(0.0),
        im: iv(0.0),
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

fn ci_affine(
    cx: f64,
    cy: f64,
    normal: Direction,
    tangent: Direction,
    s: Interval,
    r: Interval,
) -> CInterval {
    CInterval {
        re: iv(cx) + iv(normal.x) * s + iv(tangent.x) * r,
        im: iv(cy) + iv(normal.y) * s + iv(tangent.y) * r,
    }
}

fn ci_add(a: CInterval, b: CInterval) -> CInterval {
    CInterval {
        re: a.re + b.re,
        im: a.im + b.im,
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

fn ci_scale(a: CInterval, c: f64) -> CInterval {
    CInterval {
        re: a.re * iv(c),
        im: a.im * iv(c),
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

fn interval_avoids_zero(x: Interval) -> bool {
    x.inf() > 0.0 || x.sup() < 0.0
}

fn abs_sq(z: CInterval) -> Interval {
    interval_square(z.re) + interval_square(z.im)
}

fn complex_abs_upper(z: CInterval) -> f64 {
    let re = interval_abs_upper(z.re);
    let im = interval_abs_upper(z.im);
    (re * re + im * im).sqrt()
}

fn c_add(a: Complex, b: Complex) -> Complex {
    Complex {
        re: a.re + b.re,
        im: a.im + b.im,
    }
}

fn c_sub(a: Complex, b: Complex) -> Complex {
    Complex {
        re: a.re - b.re,
        im: a.im - b.im,
    }
}

fn c_mul(a: Complex, b: Complex) -> Complex {
    Complex {
        re: a.re * b.re - a.im * b.im,
        im: a.re * b.im + a.im * b.re,
    }
}

fn c_conj(a: Complex) -> Complex {
    Complex {
        re: a.re,
        im: -a.im,
    }
}

fn c_scale(a: Complex, c: f64) -> Complex {
    Complex {
        re: a.re * c,
        im: a.im * c,
    }
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

fn split_interval(pair: [f64; 2], n: usize) -> Vec<[f64; 2]> {
    let step = (pair[1] - pair[0]) / (n as f64);
    (0..n)
        .map(|i| {
            [
                pair[0] + (i as f64) * step,
                pair[0] + ((i + 1) as f64) * step,
            ]
        })
        .collect()
}

fn root_midpoints(cell: CellSpec) -> Vec<Complex> {
    let u0_pair = cell.u0_interval();
    let u1_pair = cell.u1_interval();
    let a = 0.5 * (u0_pair[0] + u0_pair[1]);
    let b = 0.5 * (u1_pair[0] + u1_pair[1]);
    let scale = eps_shape_scale(EPS);
    let base = base_roots(EPS, DEGREE);
    let u0 = u0_direction();
    let u1 = u1_direction();
    (0..DEGREE)
        .map(|i| {
            c_add(
                base[i],
                c_scale(c_add(c_scale(u0[i], a), c_scale(u1[i], b)), scale),
            )
        })
        .collect()
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

fn eval_p_p1_p2_p3(
    z: CInterval,
    roots: &[CInterval],
) -> (CInterval, CInterval, CInterval, CInterval) {
    let mut p = ci_one();
    let mut p1 = ci_zero();
    let mut p2 = ci_zero();
    let mut p3 = ci_zero();
    for root in roots {
        let factor = ci_sub(z, *root);
        let next_p3 = ci_add(ci_mul(p3, factor), ci_scale(p2, 3.0));
        let next_p2 = ci_add(ci_mul(p2, factor), ci_scale(p1, 2.0));
        let next_p1 = ci_add(ci_mul(p1, factor), p);
        let next_p = ci_mul(p, factor);
        p = next_p;
        p1 = next_p1;
        p2 = next_p2;
        p3 = next_p3;
    }
    (p, p1, p2, p3)
}

fn eval_p_p1_point(z: Complex, roots: &[Complex]) -> (Complex, Complex) {
    let mut p = Complex { re: 1.0, im: 0.0 };
    let mut p1 = Complex { re: 0.0, im: 0.0 };
    for root in roots {
        let factor = c_sub(z, *root);
        let next_p1 = c_add(c_mul(p1, factor), p);
        let next_p = c_mul(p, factor);
        p = next_p;
        p1 = next_p1;
    }
    (p, p1)
}

fn gradient_intervals(p: CInterval, p1: CInterval) -> (Interval, Interval) {
    let q = ci_mul(p1, ci_conj(p));
    (q.re * iv(2.0), q.im * iv(-2.0))
}

fn gradient_point(p: Complex, p1: Complex) -> (f64, f64) {
    let q = c_mul(p1, c_conj(p));
    (2.0 * q.re, -2.0 * q.im)
}

fn hessian_intervals(p: CInterval, p1: CInterval, p2: CInterval) -> HessianIntervals {
    let q = ci_mul(p2, ci_conj(p));
    let p1_abs_sq = abs_sq(p1);
    HessianIntervals {
        fxx: (q.re + p1_abs_sq) * iv(2.0),
        fxy: q.im * iv(-2.0),
        fyy: (p1_abs_sq - q.re) * iv(2.0),
    }
}

fn f_fx_fy_hessian_third(
    z: CInterval,
    roots: &[CInterval],
) -> (Interval, Interval, Interval, HessianIntervals, f64) {
    let (p, p1, p2, p3) = eval_p_p1_p2_p3(z, roots);
    let f = abs_sq(p) - iv(1.0);
    let (fx, fy) = gradient_intervals(p, p1);
    let hessian = hessian_intervals(p, p1, p2);
    let third_upper = 2.0 * complex_abs_upper(p3) * complex_abs_upper(p)
        + 6.0 * complex_abs_upper(p2) * complex_abs_upper(p1);
    (f, fx, fy, hessian, third_upper)
}

fn directional_second(h: HessianIntervals, a: Direction, b: Direction) -> Interval {
    h.fxx * iv(a.x * b.x) + h.fxy * iv(a.x * b.y + a.y * b.x) + h.fyy * iv(a.y * b.y)
}

fn midpoint(pair: [f64; 2]) -> f64 {
    0.5 * (pair[0] + pair[1])
}

fn width(pair: [f64; 2]) -> f64 {
    pair[1] - pair[0]
}

fn source_interval(piece: &serde_json::Value, field: &str) -> Result<[f64; 2], String> {
    let arr = piece
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

fn parse_direction(value: &serde_json::Value, field: &str) -> Option<Direction> {
    let obj = value.get(field)?.as_object()?;
    Some(Direction {
        x: obj.get("x")?.as_f64()?,
        y: obj.get("y")?.as_f64()?,
    })
}

fn parse_source_region(idx: usize, item: &serde_json::Value) -> Result<SourceRegion, String> {
    Ok(SourceRegion {
        source_index: item
            .get("source_index")
            .and_then(|v| v.as_u64())
            .unwrap_or(idx as u64) as usize,
        source_ownership_key: item
            .get("source_ownership_key")
            .and_then(|v| v.as_str())
            .unwrap_or("missing-source-ownership-key")
            .to_string(),
        split_path: item
            .get("split_path")
            .and_then(|v| v.as_str())
            .unwrap_or("missing-split-path")
            .to_string(),
        x: source_interval(item, "x_interval")?,
        y: source_interval(item, "y_interval")?,
        source_reason: item
            .get("reason")
            .and_then(|v| v.as_str())
            .unwrap_or("unknown")
            .to_string(),
        source_normal: parse_direction(item, "normal"),
        source_tangent: parse_direction(item, "tangent"),
    })
}

fn direction_from_midpoint_gradient(
    region: &SourceRegion,
    roots: &[Complex],
) -> Option<(Direction, Direction)> {
    if let (Some(normal), Some(tangent)) = (region.source_normal, region.source_tangent) {
        return Some((normal, tangent));
    }

    let cx = midpoint(region.x);
    let cy = midpoint(region.y);
    let (p, p1) = eval_p_p1_point(Complex { re: cx, im: cy }, roots);
    let (gx, gy) = gradient_point(p, p1);
    let norm = (gx * gx + gy * gy).sqrt();
    if !norm.is_finite() || norm <= 0.0 {
        return None;
    }
    let normal = Direction {
        x: gx / norm,
        y: gy / norm,
    };
    let tangent = Direction {
        x: -normal.y,
        y: normal.x,
    };
    Some((normal, tangent))
}

fn projected_ranges(
    region: &SourceRegion,
    normal: Direction,
    tangent: Direction,
) -> ([f64; 2], [f64; 2]) {
    let cx = midpoint(region.x);
    let cy = midpoint(region.y);
    let corners = [
        (region.x[0], region.y[0]),
        (region.x[0], region.y[1]),
        (region.x[1], region.y[0]),
        (region.x[1], region.y[1]),
    ];
    let mut s_lo = f64::INFINITY;
    let mut s_hi = f64::NEG_INFINITY;
    let mut r_lo = f64::INFINITY;
    let mut r_hi = f64::NEG_INFINITY;
    for (x, y) in corners {
        let dx = x - cx;
        let dy = y - cy;
        let s = dx * normal.x + dy * normal.y;
        let r = dx * tangent.x + dy * tangent.y;
        s_lo = s_lo.min(s);
        s_hi = s_hi.max(s);
        r_lo = r_lo.min(r);
        r_hi = r_hi.max(r);
    }
    ([s_lo, s_hi], [r_lo, r_hi])
}

fn fail_region(
    region: &SourceRegion,
    normal: Option<Direction>,
    tangent: Option<Direction>,
    s_range: Option<[f64; 2]>,
    r_range: Option<[f64; 2]>,
    reason: &str,
    failing_inequality: String,
    f_interval: Interval,
    fn_center_interval: Option<Interval>,
    fn_full_first_order_interval: Option<Interval>,
    ft_full_first_order_interval: Option<Interval>,
    inequality: Option<InequalityAudit>,
) -> UnresolvedThirdOrderRegion {
    UnresolvedThirdOrderRegion {
        source_index: region.source_index,
        source_ownership_key: region.source_ownership_key.clone(),
        split_path: region.split_path.clone(),
        source_reason: region.source_reason.clone(),
        x_interval: region.x,
        y_interval: region.y,
        normal,
        tangent,
        normal_interval: s_range,
        tangent_interval: r_range,
        reason: reason.to_string(),
        failing_inequality,
        f_interval: ipair(f_interval),
        fn_center_interval: fn_center_interval.map(ipair),
        fn_full_first_order_interval: fn_full_first_order_interval.map(ipair),
        ft_full_first_order_interval: ft_full_first_order_interval.map(ipair),
        inequality,
    }
}

fn analyze_region(
    region: &SourceRegion,
    roots: &[CInterval],
    midpoint_roots: &[Complex],
) -> Result<
    (
        Option<CertifiedThirdOrderCollar>,
        Option<ExcludedRegion>,
        Option<UnresolvedThirdOrderRegion>,
    ),
    String,
> {
    let axis_z = ci_box(region.x, region.y);
    let (axis_f, _, _, _, _) = f_fx_fy_hessian_third(axis_z, roots);
    if interval_avoids_zero(axis_f) {
        return Ok((
            None,
            Some(ExcludedRegion {
                source_index: region.source_index,
                source_ownership_key: region.source_ownership_key.clone(),
                split_path: region.split_path.clone(),
                x_interval: region.x,
                y_interval: region.y,
                f_interval: ipair(axis_f),
            }),
            None,
        ));
    }

    let (normal, tangent) = match direction_from_midpoint_gradient(region, midpoint_roots) {
        Some(pair) => pair,
        None => {
            return Ok((
                None,
                None,
                Some(fail_region(
                    region,
                    None,
                    None,
                    None,
                    None,
                    "third_order_critical_exclusion_failed",
                    "midpoint gradient is zero or nonfinite, so no normal direction is available"
                        .to_string(),
                    axis_f,
                    None,
                    None,
                    None,
                    None,
                )),
            ));
        }
    };

    let cx = midpoint(region.x);
    let cy = midpoint(region.y);
    let (s_range, r_range) = projected_ranges(region, normal, tangent);
    let s = iwrap(s_range[0], s_range[1]);
    let r = iwrap(r_range[0], r_range[1]);
    let z_full = ci_affine(cx, cy, normal, tangent, s, r);
    let (f_full, fx_full, fy_full, h_full, third_upper) = f_fx_fy_hessian_third(z_full, roots);
    let fn_full = fx_full * iv(normal.x) + fy_full * iv(normal.y);
    let ft_full = fx_full * iv(tangent.x) + fy_full * iv(tangent.y);

    if !third_upper.is_finite() {
        return Ok((
            None,
            None,
            Some(fail_region(
                region,
                Some(normal),
                Some(tangent),
                Some(s_range),
                Some(r_range),
                "third_order_remainder_too_wide",
                "third_directional_upper is nonfinite".to_string(),
                f_full,
                None,
                Some(fn_full),
                Some(ft_full),
                None,
            )),
        ));
    }

    let z_center = ci_affine(cx, cy, normal, tangent, iv(0.0), iv(0.0));
    let (_f_center, fx_center, fy_center, _h_center, _third_center) =
        f_fx_fy_hessian_third(z_center, roots);
    let fn_center = fx_center * iv(normal.x) + fy_center * iv(normal.y);

    let z_center_strip = ci_affine(cx, cy, normal, tangent, iv(0.0), r);
    let (f_center_strip, _, _, _, _) = f_fx_fy_hessian_third(z_center_strip, roots);

    let f_nn = directional_second(h_full, normal, normal);
    let f_nt = directional_second(h_full, normal, tangent);
    let f_tt = directional_second(h_full, tangent, tangent);

    let normal_radius_upper = s_range[0].abs().max(s_range[1].abs());
    let normal_wall_distance_lower = s_range[0].abs().min(s_range[1].abs());
    let tangent_radius_upper = r_range[0].abs().max(r_range[1].abs());
    let f_nn_abs_upper = interval_abs_upper(f_nn);
    let f_nt_abs_upper = interval_abs_upper(f_nt);
    let f_tt_abs_upper = interval_abs_upper(f_tt);
    let center_fn_abs_lower = interval_abs_lower(fn_center);

    let variation_bound_fn = f_nn_abs_upper * normal_radius_upper
        + f_nt_abs_upper * tangent_radius_upper
        + 0.5 * third_upper * (normal_radius_upper + tangent_radius_upper).powi(2);
    let critical_exclusion_margin = center_fn_abs_lower - variation_bound_fn;

    let normal_remainder_bound = 0.5 * f_nn_abs_upper * normal_radius_upper.powi(2)
        + (third_upper / 6.0) * normal_radius_upper.powi(3);
    let center_strip_f_abs_upper = interval_abs_upper(f_center_strip);
    let wall_lhs = center_strip_f_abs_upper + normal_remainder_bound;
    let wall_rhs = normal_wall_distance_lower * critical_exclusion_margin.max(0.0);
    let wall_margin = wall_rhs - wall_lhs;
    let inequality = InequalityAudit {
        center_fn_abs_lower,
        variation_bound_fn,
        critical_exclusion_margin,
        center_strip_f_abs_upper,
        normal_remainder_bound,
        wall_lhs,
        wall_rhs,
        wall_margin,
        third_directional_upper: third_upper,
        f_nn_abs_upper,
        f_nt_abs_upper,
        f_tt_abs_upper,
        normal_radius_upper,
        normal_wall_distance_lower,
        tangent_radius_upper,
    };

    if !variation_bound_fn.is_finite()
        || !normal_remainder_bound.is_finite()
        || !wall_lhs.is_finite()
    {
        return Ok((
            None,
            None,
            Some(fail_region(
                region,
                Some(normal),
                Some(tangent),
                Some(s_range),
                Some(r_range),
                "third_order_remainder_too_wide",
                format!(
                    "nonfinite remainder: variation_bound_fn={variation_bound_fn}, normal_remainder_bound={normal_remainder_bound}, wall_lhs={wall_lhs}"
                ),
                f_full,
                Some(fn_center),
                Some(fn_full),
                Some(ft_full),
                Some(inequality),
            )),
        ));
    }

    if critical_exclusion_margin <= 0.0 {
        return Ok((
            None,
            None,
            Some(fail_region(
                region,
                Some(normal),
                Some(tangent),
                Some(s_range),
                Some(r_range),
                "third_order_critical_exclusion_failed",
                format!(
                    "|F_n(0,0)|_lower - variation_bound = {center_fn_abs_lower} - {variation_bound_fn} = {critical_exclusion_margin} <= 0"
                ),
                f_full,
                Some(fn_center),
                Some(fn_full),
                Some(ft_full),
                Some(inequality),
            )),
        ));
    }

    if wall_margin <= 0.0 {
        return Ok((
            None,
            None,
            Some(fail_region(
                region,
                Some(normal),
                Some(tangent),
                Some(s_range),
                Some(r_range),
                "third_order_wall_separation_failed",
                format!(
                    "sup_r |F(0,r)| + normal_remainder = {wall_lhs} is not below S*lower(|F_n|) = {wall_rhs}; margin {wall_margin}"
                ),
                f_full,
                Some(fn_center),
                Some(fn_full),
                Some(ft_full),
                Some(inequality),
            )),
        ));
    }

    let slope_abs_upper = interval_abs_upper(ft_full) / critical_exclusion_margin;
    let length_upper = width(r_range).abs() * (1.0 + slope_abs_upper * slope_abs_upper).sqrt();
    Ok((
        Some(CertifiedThirdOrderCollar {
            source_index: region.source_index,
            source_ownership_key: region.source_ownership_key.clone(),
            split_path: region.split_path.clone(),
            normal,
            tangent,
            x_interval: region.x,
            y_interval: region.y,
            normal_interval: s_range,
            tangent_interval: r_range,
            f_interval: ipair(f_full),
            fn_center_interval: ipair(fn_center),
            fn_full_first_order_interval: ipair(fn_full),
            ft_full_first_order_interval: ipair(ft_full),
            slope_abs_upper,
            length_upper,
            inequality,
        }),
        None,
        None,
    ))
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-HARD-CELL-20260506-01_RESULTS.json"),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01"),
        experiment_id: DEFAULT_EXPERIMENT_ID.to_string(),
        cell: CellSpec::hard_cell(),
        region_limit: 0,
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
            "--experiment-id" => {
                if i + 1 >= args.len() {
                    return Err("--experiment-id requires a value".to_string());
                }
                cfg.experiment_id = args[i + 1].clone();
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

fn first_failed_condition(unresolved: &[UnresolvedThirdOrderRegion]) -> String {
    unresolved
        .first()
        .map(|r| r.reason.clone())
        .unwrap_or_else(|| "none".to_string())
}

fn status_for(
    certified_count: usize,
    excluded_count: usize,
    processed_count: usize,
    remainder_fail_count: usize,
    critical_fail_count: usize,
    wall_fail_count: usize,
    total_validated_length_upper: f64,
    cap: f64,
) -> &'static str {
    if certified_count + excluded_count == processed_count {
        if total_validated_length_upper <= cap {
            "THIRD_ORDER_COLLAR_CERTIFIES_NOT_GLOBAL_PROOF"
        } else {
            "THIRD_ORDER_FAIL_BUDGET"
        }
    } else if remainder_fail_count > 0 {
        "THIRD_ORDER_FAIL_REMAINDER_TOO_WIDE"
    } else if critical_fail_count > 0 {
        "THIRD_ORDER_FAIL_CRITICAL_EXCLUSION"
    } else if wall_fail_count > 0 {
        "THIRD_ORDER_FAIL_WALL_SEPARATION"
    } else {
        "THIRD_ORDER_FAIL_CRITICAL_EXCLUSION"
    }
}

fn next_blocker(status: &str, certified_count: usize) -> &'static str {
    match status {
        "THIRD_ORDER_COLLAR_CERTIFIES_NOT_GLOBAL_PROOF" => {
            "Third-order collars closed the processed normal-collar residuals. Rerun against a larger branch-atlas sample before returning to all unresolved slab branches."
        }
        "THIRD_ORDER_FAIL_BUDGET" => {
            "Topology closed but the length cap failed. Need sharper slope integration or a smaller accepted-length baseline before promotion."
        }
        "THIRD_ORDER_FAIL_REMAINDER_TOO_WIDE" => {
            "The third-directional bounds are too wide for this local collar. Next route is an analytic critical-point theorem with sharper derivative constants or affine/Taylor model arithmetic."
        }
        "THIRD_ORDER_FAIL_CRITICAL_EXCLUSION" if certified_count > 0 => {
            "Some collars certify, but critical-point exclusion still fails elsewhere. Expand only after isolating the certified pattern."
        }
        "THIRD_ORDER_FAIL_CRITICAL_EXCLUSION" => {
            "The Taylor critical-exclusion inequality still does not close. Retire sampled local-collar escalation and write the analytic critical-point theorem target."
        }
        "THIRD_ORDER_FAIL_WALL_SEPARATION" => {
            "Critical exclusion passed where tested, but wall separation did not. The next inequality must sharpen center-strip or normal-remainder control."
        }
        _ => "Inspect worst_region before choosing the next proof-facing inequality.",
    }
}

fn write_report(result: &serde_json::Value, path: &PathBuf) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Third-Order Root-Collar Remainder Pilot\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Processed regions: `{}`\n\
- Third-order certified regions: `{}`\n\
- Critical-exclusion passes: `{}`\n\
- Wall-separation passes: `{}`\n\
- Remaining unresolved regions: `{}`\n\
- Resolved length upper: `{}`\n\
- Total validated length upper: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap: `{}`\n\
- First failed condition: `{}`\n\n\
## Failing Inequality\n\n\
{}\n\n\
## Interpretation\n\n\
This run keeps the midpoint-gradient normal/tangent coordinates from the \
normal-collar artifact and replaces the first-order collar test with an \
explicit Taylor remainder check. The recurrence computes `p`, `p'`, `p''`, \
and `p'''`; the logged third-directional upper bound is used to audit whether \
critical-point exclusion and wall separation can be certified on each residual \
region.\n\n\
## Next Blocker\n\n\
{}\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        SOURCE_EXPERIMENT_ID,
        result["status"],
        result["processed_region_count"],
        result["third_order_certified_count"],
        result["critical_exclusion_pass_count"],
        result["wall_separation_pass_count"],
        result["remaining_unresolved_region_count"],
        result["resolved_length_upper"],
        result["total_validated_length_upper"],
        result["exact_length_cap"],
        result["margin_to_cap"],
        result["first_failed_condition"],
        result["first_failing_inequality"]
            .as_str()
            .unwrap_or("none"),
        result["next_blocker"]
            .as_str()
            .unwrap_or("continue diagnostics"),
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
    fs::create_dir_all(&cfg.out_dir).map_err(|err| {
        format!(
            "failed to create output directory {}: {err}",
            cfg.out_dir.display()
        )
    })?;
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

    let source_text = fs::read_to_string(&cfg.source)
        .map_err(|err| format!("failed to read source {}: {err}", cfg.source.display()))?;
    let source_json: serde_json::Value =
        serde_json::from_str(&source_text).map_err(|err| format!("bad source JSON: {err}"))?;
    let source_subcell_contract = ensure_source_subcell_matches(&source_json, cfg.cell)?;
    let source_accepted_length = source_json["source_accepted_length_upper"]
        .as_f64()
        .ok_or_else(|| "source missing source_accepted_length_upper".to_string())?;
    let cap = source_json["exact_length_cap"]
        .as_f64()
        .ok_or_else(|| "source missing exact_length_cap".to_string())?;
    let regions_json = source_json["remaining_unresolved"]
        .as_array()
        .ok_or_else(|| "source missing remaining_unresolved".to_string())?;
    let limit = if cfg.region_limit == 0 {
        regions_json.len()
    } else {
        cfg.region_limit.min(regions_json.len())
    };

    let roots = root_intervals(cfg.cell);
    let midpoint_roots = root_midpoints(cfg.cell);
    let mut certified = Vec::<CertifiedThirdOrderCollar>::new();
    let mut excluded = Vec::<ExcludedRegion>::new();
    let mut unresolved = Vec::<UnresolvedThirdOrderRegion>::new();

    for (idx, value) in regions_json.iter().take(limit).enumerate() {
        let region = parse_source_region(idx, value)?;
        let (cert, excl, unres) = analyze_region(&region, &roots, &midpoint_roots)?;
        if let Some(item) = cert {
            certified.push(item);
        }
        if let Some(item) = excl {
            excluded.push(item);
        }
        if let Some(item) = unres {
            unresolved.push(item);
        }
    }

    let critical_exclusion_pass_count = certified.len()
        + unresolved
            .iter()
            .filter(|r| {
                r.inequality
                    .as_ref()
                    .map(|a| a.critical_exclusion_margin > 0.0)
                    .unwrap_or(false)
            })
            .count();
    let wall_separation_pass_count = certified.len();
    let remainder_fail_count = unresolved
        .iter()
        .filter(|r| r.reason == "third_order_remainder_too_wide")
        .count();
    let critical_fail_count = unresolved
        .iter()
        .filter(|r| r.reason == "third_order_critical_exclusion_failed")
        .count();
    let wall_fail_count = unresolved
        .iter()
        .filter(|r| r.reason == "third_order_wall_separation_failed")
        .count();
    let resolved_length_upper: f64 = certified.iter().map(|c| c.length_upper).sum();
    let total_validated_length_upper = source_accepted_length + resolved_length_upper;
    let margin_to_cap = cap - total_validated_length_upper;
    let status = status_for(
        certified.len(),
        excluded.len(),
        limit,
        remainder_fail_count,
        critical_fail_count,
        wall_fail_count,
        total_validated_length_upper,
        cap,
    );
    let failed = first_failed_condition(&unresolved);
    let first_failing_inequality = unresolved
        .first()
        .map(|r| r.failing_inequality.clone())
        .unwrap_or_else(|| "none".to_string());
    let worst_region = unresolved.first().cloned();
    let result = json!({
        "experiment_id": cfg.experiment_id,
        "timestamp_unix": unix_timestamp_string(),
        "source_experiment_id": SOURCE_EXPERIMENT_ID,
        "source_results_path": cfg.source,
        "source_subcell_contract": source_subcell_contract,
        "status": status,
        "degree": DEGREE,
        "eps": EPS,
        "subcell": cfg.cell,
        "cell_tag": cfg.cell.tag(),
        "parameters": {
            "region_limit": cfg.region_limit,
            "processed_region_count": limit,
            "source_region_field": "remaining_unresolved",
            "coordinate_rule": "reuse source normal/tangent when present; otherwise normal is midpoint gradient and tangent is normal rotated 90 degrees",
            "critical_exclusion_rule": "|F_n(0,0)|_lower - variation_bound(F_n over s,r) > 0",
            "wall_separation_rule": "sup_r |F(0,r)| + normal_remainder < S * lower_bound(|F_n|)",
            "third_derivative_rule": "|d^3/dv^3 |p|^2| <= 2|p'''||p| + 6|p''||p'| for unit directions",
            "length_rule": "width_tangent * sqrt(1 + sup(|Ft/Fn|)^2)"
        },
        "source_region_count": regions_json.len(),
        "processed_region_count": limit,
        "source_unprocessed_region_count": regions_json.len().saturating_sub(limit),
        "source_accepted_length_upper": source_accepted_length,
        "third_order_certified_count": certified.len(),
        "excluded_region_count": excluded.len(),
        "critical_exclusion_pass_count": critical_exclusion_pass_count,
        "wall_separation_pass_count": wall_separation_pass_count,
        "remainder_fail_count": remainder_fail_count,
        "critical_exclusion_fail_count": critical_fail_count,
        "wall_separation_fail_count": wall_fail_count,
        "remaining_unresolved_region_count": unresolved.len(),
        "resolved_length_upper": resolved_length_upper,
        "total_validated_length_upper": total_validated_length_upper,
        "exact_length_cap": cap,
        "margin_to_cap": margin_to_cap,
        "first_failed_condition": failed,
        "first_failing_inequality": first_failing_inequality,
        "certified_regions": certified,
        "excluded_regions": excluded,
        "remaining_unresolved": unresolved,
        "worst_region": worst_region,
        "claim_ceiling": "Local n=14 hard-cell third-order collar remainder pilot only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.",
        "next_blocker": next_blocker(status, certified.len()),
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
    write_report(&result, &report_path)?;
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
            "experiment_id": result["experiment_id"],
            "status": status,
            "processed_region_count": limit,
            "third_order_certified_count": result["third_order_certified_count"],
            "excluded_region_count": result["excluded_region_count"],
            "critical_exclusion_pass_count": result["critical_exclusion_pass_count"],
            "wall_separation_pass_count": result["wall_separation_pass_count"],
            "remaining_unresolved_region_count": result["remaining_unresolved_region_count"],
            "total_validated_length_upper": result["total_validated_length_upper"],
            "margin_to_cap": result["margin_to_cap"],
            "first_failed_condition": result["first_failed_condition"],
            "first_failing_inequality": result["first_failing_inequality"],
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
