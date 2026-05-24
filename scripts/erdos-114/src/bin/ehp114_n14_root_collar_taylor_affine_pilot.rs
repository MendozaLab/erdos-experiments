//! EHP #114 n=14 root-collar Taylor/affine hard-cell pilot.
//!
//! This diagnostic consumes the latest bivariate Bernstein/Krawczyk failure
//! artifact and tests the next proof-facing idea on a small prefix of the
//! remaining hard pieces: centered Taylor/affine root collars should preserve
//! first-order dependency where Bernstein hulls widened.
//!
//! This is a local pilot only. It is not a proof of Erdős #114, not a global
//! n=14 proof, and not an exact lemniscate-length certificate.

use serde::Serialize;
use serde_json::json;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const EXPERIMENT_ID: &str = "EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01";
const SOURCE_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-02";
const DEGREE: usize = 14;
const EPS: f64 = 0.1;
const SUBDIVISION: usize = 8;
const SUB_I: usize = 6;
const SUB_J: usize = 4;
const U0_PARENT: [f64; 2] = [-0.0017499999999999998, 0.0];
const U1_PARENT: [f64; 2] = [0.0, 0.0017499999999999998];

#[derive(Clone, Copy, Debug, Serialize)]
struct Complex {
    re: f64,
    im: f64,
}

#[derive(Clone, Copy, Debug, Serialize)]
struct RInterval {
    lo: f64,
    hi: f64,
}

#[derive(Clone, Debug, Serialize)]
struct PieceSummary {
    source_index: usize,
    source_tube_index: Option<usize>,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    source_reason: String,
    source_f_interval: [f64; 2],
    source_fx_interval: [f64; 2],
    source_fy_interval: [f64; 2],
    bernstein_f_interval: Option<[f64; 2]>,
    bernstein_fx_interval: Option<[f64; 2]>,
    bernstein_fy_interval: Option<[f64; 2]>,
    midpoint_f: f64,
    midpoint_fx: f64,
    midpoint_fy: f64,
    hessian_bounds: HessianBounds,
    taylor_f_interval: [f64; 2],
    taylor_fx_interval: [f64; 2],
    taylor_fy_interval: [f64; 2],
    f_width_improvement_vs_bernstein: Option<f64>,
    fx_width_improvement_vs_bernstein: Option<f64>,
    fy_width_improvement_vs_bernstein: Option<f64>,
    source_derivative_signed: bool,
    taylor_derivative_signed: bool,
    best_projection_axis: String,
    root_collar_candidate: bool,
    root_collar_reason: String,
    contracted_interval: Option<[f64; 2]>,
    dependent_interval: Option<[f64; 2]>,
    slope_abs_upper: Option<f64>,
    length_upper_candidate: Option<f64>,
}

#[derive(Clone, Copy, Debug, Serialize)]
struct HessianBounds {
    fxx_abs_upper: f64,
    fxy_abs_upper: f64,
    fyy_abs_upper: f64,
    inflation_factor: f64,
}

#[derive(Clone, Copy, Debug)]
struct PointStats {
    f: f64,
    fx: f64,
    fy: f64,
    fxx: f64,
    fxy: f64,
    fyy: f64,
}

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    out_dir: PathBuf,
    piece_limit: usize,
    hessian_inflation: f64,
}

fn add(a: Complex, b: Complex) -> Complex {
    Complex {
        re: a.re + b.re,
        im: a.im + b.im,
    }
}

fn sub(a: Complex, b: Complex) -> Complex {
    Complex {
        re: a.re - b.re,
        im: a.im - b.im,
    }
}

fn mul(a: Complex, b: Complex) -> Complex {
    Complex {
        re: a.re * b.re - a.im * b.im,
        im: a.re * b.im + a.im * b.re,
    }
}

fn scale(a: Complex, c: f64) -> Complex {
    Complex {
        re: a.re * c,
        im: a.im * c,
    }
}

fn conj(a: Complex) -> Complex {
    Complex {
        re: a.re,
        im: -a.im,
    }
}

fn abs_sq(a: Complex) -> f64 {
    a.re * a.re + a.im * a.im
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

fn midpoint_roots() -> Vec<Complex> {
    let u0_pair = split_interval(U0_PARENT, SUBDIVISION)[SUB_I];
    let u1_pair = split_interval(U1_PARENT, SUBDIVISION)[SUB_J];
    let a = 0.5 * (u0_pair[0] + u0_pair[1]);
    let b = 0.5 * (u1_pair[0] + u1_pair[1]);
    let scale_eps = eps_shape_scale(EPS);
    let base = base_roots(EPS, DEGREE);
    let u0 = u0_direction();
    let u1 = u1_direction();
    (0..DEGREE)
        .map(|i| {
            add(
                base[i],
                scale(add(scale(u0[i], a), scale(u1[i], b)), scale_eps),
            )
        })
        .collect()
}

fn eval_p_p1_p2_point(z: Complex, roots: &[Complex]) -> (Complex, Complex, Complex) {
    let mut p = Complex { re: 1.0, im: 0.0 };
    let mut p1 = Complex { re: 0.0, im: 0.0 };
    let mut p2 = Complex { re: 0.0, im: 0.0 };
    for root in roots {
        let factor = sub(z, *root);
        let next_p2 = add(mul(p2, factor), scale(p1, 2.0));
        let next_p1 = add(mul(p1, factor), p);
        let next_p = mul(p, factor);
        p = next_p;
        p1 = next_p1;
        p2 = next_p2;
    }
    (p, p1, p2)
}

fn point_stats(z: Complex, roots: &[Complex]) -> PointStats {
    let (p, p1, p2) = eval_p_p1_p2_point(z, roots);
    let q = mul(p1, conj(p));
    let r = mul(p2, conj(p));
    let p1_abs_sq = abs_sq(p1);
    PointStats {
        f: abs_sq(p) - 1.0,
        fx: 2.0 * q.re,
        fy: -2.0 * q.im,
        fxx: 2.0 * r.re + 2.0 * p1_abs_sq,
        fxy: -2.0 * r.im,
        fyy: -2.0 * r.re + 2.0 * p1_abs_sq,
    }
}

fn interval_from_value(value: &serde_json::Value, field: &str) -> Result<[f64; 2], String> {
    let arr = value
        .get(field)
        .and_then(|v| v.as_array())
        .ok_or_else(|| format!("missing interval field {field}"))?;
    if arr.len() != 2 {
        return Err(format!("interval field {field} does not have length 2"));
    }
    let lo = arr[0]
        .as_f64()
        .ok_or_else(|| format!("interval field {field}[0] is not numeric"))?;
    let hi = arr[1]
        .as_f64()
        .ok_or_else(|| format!("interval field {field}[1] is not numeric"))?;
    Ok([lo, hi])
}

fn optional_interval_from_value(value: &serde_json::Value, field: &str) -> Option<[f64; 2]> {
    interval_from_value(value, field).ok()
}

fn width(pair: [f64; 2]) -> f64 {
    pair[1] - pair[0]
}

fn interval_width(pair: Option<[f64; 2]>) -> Option<f64> {
    pair.map(width).filter(|w| *w > 0.0 && w.is_finite())
}

fn improvement(old: Option<[f64; 2]>, new: [f64; 2]) -> Option<f64> {
    let old_w = interval_width(old)?;
    let new_w = width(new);
    if new_w > 0.0 && new_w.is_finite() {
        Some(old_w / new_w)
    } else {
        None
    }
}

fn midpoint(pair: [f64; 2]) -> f64 {
    0.5 * (pair[0] + pair[1])
}

fn radius(pair: [f64; 2]) -> f64 {
    0.5 * (pair[1] - pair[0]).abs()
}

fn sign_stable(pair: [f64; 2]) -> bool {
    pair[0] > 0.0 || pair[1] < 0.0
}

fn abs_lower(pair: [f64; 2]) -> f64 {
    if pair[0] <= 0.0 && pair[1] >= 0.0 {
        0.0
    } else {
        pair[0].abs().min(pair[1].abs())
    }
}

fn abs_upper(pair: [f64; 2]) -> f64 {
    pair[0].abs().max(pair[1].abs())
}

fn interval(lo: f64, hi: f64) -> RInterval {
    RInterval {
        lo: lo.min(hi),
        hi: lo.max(hi),
    }
}

fn interval_add(a: RInterval, b: RInterval) -> RInterval {
    interval(a.lo + b.lo, a.hi + b.hi)
}

fn interval_sub(a: RInterval, b: RInterval) -> RInterval {
    interval(a.lo - b.hi, a.hi - b.lo)
}

fn interval_mul(a: RInterval, b: RInterval) -> RInterval {
    let vals = [a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi];
    interval(
        vals.iter().copied().fold(f64::INFINITY, f64::min),
        vals.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    )
}

fn interval_scale(a: RInterval, c: f64) -> RInterval {
    interval(a.lo * c, a.hi * c)
}

fn interval_square_symmetric(radius: f64, coeff: f64) -> RInterval {
    let hi = coeff.abs() * radius * radius;
    interval(-hi, hi)
}

fn interval_to_pair(x: RInterval) -> [f64; 2] {
    [x.lo, x.hi]
}

fn strict_subset(inner: RInterval, outer: [f64; 2]) -> bool {
    inner.lo.is_finite() && inner.hi.is_finite() && inner.lo > outer[0] && inner.hi < outer[1]
}

fn sample_hessian_bounds(
    x: [f64; 2],
    y: [f64; 2],
    roots: &[Complex],
    inflation: f64,
) -> HessianBounds {
    let cx = midpoint(x);
    let cy = midpoint(y);
    let xs = [x[0], cx, x[1]];
    let ys = [y[0], cy, y[1]];
    let mut fxx = 0.0_f64;
    let mut fxy = 0.0_f64;
    let mut fyy = 0.0_f64;
    for xx in xs {
        for yy in ys {
            let s = point_stats(Complex { re: xx, im: yy }, roots);
            fxx = fxx.max(s.fxx.abs());
            fxy = fxy.max(s.fxy.abs());
            fyy = fyy.max(s.fyy.abs());
        }
    }
    HessianBounds {
        fxx_abs_upper: fxx * inflation,
        fxy_abs_upper: fxy * inflation,
        fyy_abs_upper: fyy * inflation,
        inflation_factor: inflation,
    }
}

fn taylor_enclosures(
    center: PointStats,
    x: [f64; 2],
    y: [f64; 2],
    h: HessianBounds,
) -> ([f64; 2], [f64; 2], [f64; 2]) {
    let rx = radius(x);
    let ry = radius(y);
    let f_rad = center.fx.abs() * rx
        + center.fy.abs() * ry
        + 0.5
            * (h.fxx_abs_upper * rx * rx
                + 2.0 * h.fxy_abs_upper * rx * ry
                + h.fyy_abs_upper * ry * ry);
    let fx_rad = h.fxx_abs_upper * rx + h.fxy_abs_upper * ry;
    let fy_rad = h.fxy_abs_upper * rx + h.fyy_abs_upper * ry;
    (
        [center.f - f_rad, center.f + f_rad],
        [center.fx - fx_rad, center.fx + fx_rad],
        [center.fy - fy_rad, center.fy + fy_rad],
    )
}

fn candidate_x_as_function_of_y(
    center: PointStats,
    x: [f64; 2],
    y: [f64; 2],
    h: HessianBounds,
    source_fx: [f64; 2],
    source_fy: [f64; 2],
) -> Option<(RInterval, f64, f64)> {
    if center.fx == 0.0 || !sign_stable(source_fx) {
        return None;
    }
    let cx = midpoint(x);
    let cy = midpoint(y);
    let dy = interval(y[0] - cy, y[1] - cy);
    let dx = interval(x[0] - cx, x[1] - cx);
    let f_on_center_x = interval_add(
        interval(center.f, center.f),
        interval_add(
            interval_scale(dy, center.fy),
            interval_square_symmetric(radius(y), 0.5 * h.fyy_abs_upper),
        ),
    );
    let fx_rad = h.fxx_abs_upper * radius(x) + h.fxy_abs_upper * radius(y);
    let fx_range = interval(center.fx - fx_rad, center.fx + fx_rad);
    let c = 1.0 / center.fx;
    let k = interval_add(
        interval(cx, cx),
        interval_add(
            interval_scale(f_on_center_x, -c),
            interval_mul(
                interval_sub(interval(1.0, 1.0), interval_scale(fx_range, c)),
                dx,
            ),
        ),
    );
    if strict_subset(k, x) {
        let slope = abs_upper(source_fy) / abs_lower(source_fx);
        let length = width(y) * (1.0 + slope * slope).sqrt();
        Some((k, slope, length))
    } else {
        None
    }
}

fn candidate_y_as_function_of_x(
    center: PointStats,
    x: [f64; 2],
    y: [f64; 2],
    h: HessianBounds,
    source_fx: [f64; 2],
    source_fy: [f64; 2],
) -> Option<(RInterval, f64, f64)> {
    if center.fy == 0.0 || !sign_stable(source_fy) {
        return None;
    }
    let cx = midpoint(x);
    let cy = midpoint(y);
    let dx = interval(x[0] - cx, x[1] - cx);
    let dy = interval(y[0] - cy, y[1] - cy);
    let f_on_center_y = interval_add(
        interval(center.f, center.f),
        interval_add(
            interval_scale(dx, center.fx),
            interval_square_symmetric(radius(x), 0.5 * h.fxx_abs_upper),
        ),
    );
    let fy_rad = h.fxy_abs_upper * radius(x) + h.fyy_abs_upper * radius(y);
    let fy_range = interval(center.fy - fy_rad, center.fy + fy_rad);
    let c = 1.0 / center.fy;
    let k = interval_add(
        interval(cy, cy),
        interval_add(
            interval_scale(f_on_center_y, -c),
            interval_mul(
                interval_sub(interval(1.0, 1.0), interval_scale(fy_range, c)),
                dy,
            ),
        ),
    );
    if strict_subset(k, y) {
        let slope = abs_upper(source_fx) / abs_lower(source_fy);
        let length = width(x) * (1.0 + slope * slope).sqrt();
        Some((k, slope, length))
    } else {
        None
    }
}

fn analyze_piece(
    idx: usize,
    piece: &serde_json::Value,
    roots: &[Complex],
    hessian_inflation: f64,
) -> Result<PieceSummary, String> {
    let x = interval_from_value(piece, "x_interval")?;
    let y = interval_from_value(piece, "y_interval")?;
    let source_f = interval_from_value(piece, "f_interval")?;
    let source_fx = interval_from_value(piece, "fx_interval")?;
    let source_fy = interval_from_value(piece, "fy_interval")?;
    let bernstein_f = optional_interval_from_value(piece, "bernstein_f_interval");
    let bernstein_fx = optional_interval_from_value(piece, "bernstein_fx_interval");
    let bernstein_fy = optional_interval_from_value(piece, "bernstein_fy_interval");
    let center = point_stats(
        Complex {
            re: midpoint(x),
            im: midpoint(y),
        },
        roots,
    );
    let h = sample_hessian_bounds(x, y, roots, hessian_inflation);
    let (taylor_f, taylor_fx, taylor_fy) = taylor_enclosures(center, x, y, h);
    let source_derivative_signed = sign_stable(source_fx) || sign_stable(source_fy);
    let taylor_derivative_signed = sign_stable(taylor_fx) || sign_stable(taylor_fy);

    let x_candidate = candidate_x_as_function_of_y(center, x, y, h, source_fx, source_fy);
    let y_candidate = candidate_y_as_function_of_x(center, x, y, h, source_fx, source_fy);
    let chosen = match (x_candidate, y_candidate) {
        (Some(a), Some(b)) => {
            if a.2 <= b.2 {
                Some(("x_as_function_of_y", a, x))
            } else {
                Some(("y_as_function_of_x", b, y))
            }
        }
        (Some(a), None) => Some(("x_as_function_of_y", a, x)),
        (None, Some(b)) => Some(("y_as_function_of_x", b, y)),
        (None, None) => None,
    };

    let (
        root_collar_candidate,
        root_collar_reason,
        contracted_interval,
        dependent_interval,
        slope_abs_upper,
        length_upper_candidate,
        best_projection_axis,
    ) = if let Some((axis, (contracted, slope, length), dependent)) = chosen {
        (
            true,
            "centered_taylor_krawczyk_subset_candidate".to_string(),
            Some(interval_to_pair(contracted)),
            Some(dependent),
            Some(slope),
            Some(length),
            axis.to_string(),
        )
    } else if !source_derivative_signed {
        (
            false,
            "source_derivative_not_sign_stable".to_string(),
            None,
            None,
            None,
            None,
            "none".to_string(),
        )
    } else if !taylor_derivative_signed {
        (
            false,
            "taylor_derivative_not_sign_stable".to_string(),
            None,
            None,
            None,
            None,
            "none".to_string(),
        )
    } else {
        let axis = if abs_lower(source_fx) >= abs_lower(source_fy) {
            "x_as_function_of_y"
        } else {
            "y_as_function_of_x"
        };
        (
            false,
            "taylor_krawczyk_not_strict_subset".to_string(),
            None,
            None,
            None,
            None,
            axis.to_string(),
        )
    };

    Ok(PieceSummary {
        source_index: idx,
        source_tube_index: piece
            .get("source_tube_index")
            .and_then(|v| v.as_u64())
            .map(|v| v as usize),
        x_interval: x,
        y_interval: y,
        source_reason: piece
            .get("reason")
            .and_then(|v| v.as_str())
            .unwrap_or("unknown")
            .to_string(),
        source_f_interval: source_f,
        source_fx_interval: source_fx,
        source_fy_interval: source_fy,
        bernstein_f_interval: bernstein_f,
        bernstein_fx_interval: bernstein_fx,
        bernstein_fy_interval: bernstein_fy,
        midpoint_f: center.f,
        midpoint_fx: center.fx,
        midpoint_fy: center.fy,
        hessian_bounds: h,
        taylor_f_interval: taylor_f,
        taylor_fx_interval: taylor_fx,
        taylor_fy_interval: taylor_fy,
        f_width_improvement_vs_bernstein: improvement(bernstein_f, taylor_f),
        fx_width_improvement_vs_bernstein: improvement(bernstein_fx, taylor_fx),
        fy_width_improvement_vs_bernstein: improvement(bernstein_fy, taylor_fy),
        source_derivative_signed,
        taylor_derivative_signed,
        best_projection_axis,
        root_collar_candidate,
        root_collar_reason,
        contracted_interval,
        dependent_interval,
        slope_abs_upper,
        length_upper_candidate,
    })
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-02/EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-02_RESULTS.json"),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01"),
        piece_limit: 10,
        hessian_inflation: 1.25,
    };
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
            "--piece-limit" => {
                if i + 1 >= args.len() {
                    return Err("--piece-limit requires an integer".to_string());
                }
                cfg.piece_limit = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --piece-limit: {err}"))?;
                i += 2;
            }
            "--hessian-inflation" => {
                if i + 1 >= args.len() {
                    return Err("--hessian-inflation requires a float".to_string());
                }
                cfg.hessian_inflation = args[i + 1]
                    .parse::<f64>()
                    .map_err(|err| format!("failed to parse --hessian-inflation: {err}"))?;
                i += 2;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
    }
    Ok(cfg)
}

fn sha256_file(path: &PathBuf) -> Result<String, String> {
    let bytes = fs::read(path)
        .map_err(|err| format!("failed to read {} for sha256: {err}", path.display()))?;
    Ok(sha256::digest(bytes))
}

fn unix_timestamp_string() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs()
        .to_string()
}

fn write_report(result: &serde_json::Value, path: &PathBuf) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Root-Collar Taylor/Affine Pilot\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Processed pieces: `{}`\n\
- Source unprocessed pieces: `{}`\n\
- Source accepted length: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap: `{}`\n\
- Source derivative signed count: `{}`\n\
- Taylor derivative signed count: `{}`\n\
- Taylor tighter than Bernstein count: `{}`\n\
- Root-collar candidate count: `{}`\n\
- Remaining non-candidate pieces: `{}`\n\n\
## Interpretation\n\n\
This pilot tests whether a centered Taylor/affine local model recovers dependency \
information that the bivariate Bernstein hull lost. It uses midpoint root-affine \
roots and sampled Hessian inflation as a diagnostic. It is therefore proof-facing \
triage, not a validated proof certificate. A useful outcome is either a strict \
root-collar candidate on at least one worst piece, or a clean failure that sends \
the route toward a hand analytic collar theorem.\n\n\
## Next Blocker\n\n\
{}\n\n\
## Claim Ceiling\n\n\
{}\n",
        EXPERIMENT_ID,
        SOURCE_EXPERIMENT_ID,
        result["status"],
        result["processed_piece_count"],
        result["source_unprocessed_piece_count"],
        result["source_accepted_length_upper"],
        result["exact_length_cap"],
        result["margin_to_cap"],
        result["source_derivative_signed_count"],
        result["taylor_derivative_signed_count"],
        result["taylor_tighter_than_bernstein_count"],
        result["root_collar_candidate_count"],
        result["remaining_non_candidate_piece_count"],
        result["next_blocker"]
            .as_str()
            .unwrap_or("continue diagnostics"),
        result["claim_ceiling"]
            .as_str()
            .unwrap_or("local pilot only"),
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
    let result_path = cfg.out_dir.join(format!("{}_RESULTS.json", EXPERIMENT_ID));
    let report_path = cfg.out_dir.join(format!("{}_REPORT.md", EXPERIMENT_ID));
    let sha_path = cfg
        .out_dir
        .join(format!("{}_RESULTS.sha256", EXPERIMENT_ID));
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
    let unresolved = source_json["remaining_unresolved_pieces"]
        .as_array()
        .ok_or_else(|| "source missing remaining_unresolved_pieces".to_string())?;
    let source_accepted_length = source_json["source_accepted_length_upper"]
        .as_f64()
        .or_else(|| source_json["total_validated_length_upper"].as_f64())
        .ok_or_else(|| "source missing source accepted length".to_string())?;
    let cap = source_json["exact_length_cap"]
        .as_f64()
        .ok_or_else(|| "source missing exact_length_cap".to_string())?;

    let limit = if cfg.piece_limit == 0 {
        unresolved.len()
    } else {
        cfg.piece_limit.min(unresolved.len())
    };
    let roots = midpoint_roots();
    let mut pieces = Vec::with_capacity(limit);
    for (idx, piece) in unresolved.iter().take(limit).enumerate() {
        pieces.push(analyze_piece(idx, piece, &roots, cfg.hessian_inflation)?);
    }

    let root_collar_candidate_count = pieces.iter().filter(|p| p.root_collar_candidate).count();
    let source_derivative_signed_count =
        pieces.iter().filter(|p| p.source_derivative_signed).count();
    let taylor_derivative_signed_count =
        pieces.iter().filter(|p| p.taylor_derivative_signed).count();
    let taylor_tighter_than_bernstein_count = pieces
        .iter()
        .filter(|p| p.f_width_improvement_vs_bernstein.unwrap_or(0.0) > 10.0)
        .count();
    let remaining_non_candidate_piece_count = limit.saturating_sub(root_collar_candidate_count);
    let candidate_length_upper_sum: f64 =
        pieces.iter().filter_map(|p| p.length_upper_candidate).sum();
    let diagnostic_total_length_if_candidates_promoted =
        source_accepted_length + candidate_length_upper_sum;
    let margin_to_cap = cap - diagnostic_total_length_if_candidates_promoted;

    let status = if root_collar_candidate_count > 0 {
        "ROOT_COLLAR_TAYLOR_AFFINE_PILOT_CANDIDATES_NOT_CERTIFICATE"
    } else if taylor_tighter_than_bernstein_count > 0 && taylor_derivative_signed_count > 0 {
        "ROOT_COLLAR_TAYLOR_AFFINE_PILOT_DEPENDENCY_RECOVERED_NO_COLLAR"
    } else {
        "ROOT_COLLAR_TAYLOR_AFFINE_PILOT_FAILS"
    };
    let next_blocker = if root_collar_candidate_count > 0 {
        "Promote the centered Taylor/affine candidates to validated interval/Taylor-model arithmetic over root-affine uncertainty; do not scale until one candidate is proof-grade."
    } else if taylor_tighter_than_bernstein_count > 0 {
        "Taylor/affine coordinates recover dependency but do not yet produce a strict collar. Next step is a validated Taylor model with root-affine uncertainty or a hand analytic root-collar lemma."
    } else {
        "The midpoint Taylor/affine pilot did not materially recover the pieces. Move directly to an analytic root-collar theorem rather than more numeric subdivision."
    };

    let worst_piece = pieces
        .iter()
        .find(|p| !p.root_collar_candidate)
        .or_else(|| pieces.first());

    let result = json!({
        "experiment_id": EXPERIMENT_ID,
        "timestamp_unix": unix_timestamp_string(),
        "source_experiment_id": SOURCE_EXPERIMENT_ID,
        "source_results_path": cfg.source,
        "status": status,
        "degree": DEGREE,
        "eps": EPS,
        "subcell": {"sub_i": SUB_I, "sub_j": SUB_J},
        "parameters": {
            "piece_limit": cfg.piece_limit,
            "processed_piece_count": limit,
            "hessian_inflation": cfg.hessian_inflation,
            "root_parameter_model": "midpoint of root-affine subcell (diagnostic only)",
            "pilot_rule": "centered Taylor/affine local model with sampled Hessian inflation; candidate only if Krawczyk-style collar is a strict subset"
        },
        "source_unresolved_piece_count": unresolved.len(),
        "processed_piece_count": limit,
        "source_unprocessed_piece_count": unresolved.len().saturating_sub(limit),
        "source_accepted_length_upper": source_accepted_length,
        "candidate_length_upper_sum": candidate_length_upper_sum,
        "diagnostic_total_length_if_candidates_promoted": diagnostic_total_length_if_candidates_promoted,
        "exact_length_cap": cap,
        "margin_to_cap": margin_to_cap,
        "source_derivative_signed_count": source_derivative_signed_count,
        "taylor_derivative_signed_count": taylor_derivative_signed_count,
        "taylor_tighter_than_bernstein_count": taylor_tighter_than_bernstein_count,
        "root_collar_candidate_count": root_collar_candidate_count,
        "certificate_candidate_count": root_collar_candidate_count,
        "remaining_non_candidate_piece_count": remaining_non_candidate_piece_count,
        "piece_summaries": pieces,
        "worst_piece": worst_piece,
        "claim_ceiling": "Root-collar Taylor/affine pilot only. Not a proof of Erdős #114, not a global n=14 proof, and not an exact lemniscate-length certificate.",
        "next_blocker": next_blocker,
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
            "experiment_id": EXPERIMENT_ID,
            "status": status,
            "processed_piece_count": limit,
            "source_unresolved_piece_count": unresolved.len(),
            "source_derivative_signed_count": source_derivative_signed_count,
            "taylor_derivative_signed_count": taylor_derivative_signed_count,
            "taylor_tighter_than_bernstein_count": taylor_tighter_than_bernstein_count,
            "root_collar_candidate_count": root_collar_candidate_count,
            "remaining_non_candidate_piece_count": remaining_non_candidate_piece_count,
            "result": result_path,
            "report": report_path,
            "sha256": sha,
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
