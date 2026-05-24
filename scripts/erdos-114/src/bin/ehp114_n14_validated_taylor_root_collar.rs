//! EHP #114 n=14 validated Taylor root-collar pilot.
//!
//! This binary is the proof-facing follow-up to the midpoint Taylor/affine
//! diagnostic. It carries the root-affine uncertainty with interval roots,
//! builds Taylor enclosures for F, Fx, and Fy using interval Hessian bounds,
//! and then attempts a strict interval Newton/Krawczyk collar on a small set
//! of hard pieces.
//!
//! This remains local: it is not a proof of Erdos #114, not a global n=14
//! certificate, and not an exact lemniscate-length certificate.

use inari::{interval, Interval};
use serde::Serialize;
use serde_json::json;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-VALIDATED-TAYLOR-ROOT-COLLAR-HARD-CELL-20260506-01";
const SOURCE_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01";
const DEGREE: usize = 14;
const EPS: f64 = 0.1;
const SUBDIVISION: usize = 8;
const SUB_I: usize = 6;
const SUB_J: usize = 4;
const U0_PARENT: [f64; 2] = [-0.0017499999999999998, 0.0];
const U1_PARENT: [f64; 2] = [0.0, 0.0017499999999999998];

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

#[derive(Clone, Copy, Debug, Serialize)]
struct TaylorHessianBox {
    fxx: IntervalPair,
    fxy: IntervalPair,
    fyy: IntervalPair,
}

#[derive(Clone, Debug, Serialize)]
struct PieceResult {
    source_index: usize,
    source_tube_index: Option<usize>,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    source_reason: String,
    source_fx_interval: [f64; 2],
    source_fy_interval: [f64; 2],
    prior_taylor_fx_interval: Option<[f64; 2]>,
    prior_taylor_fy_interval: Option<[f64; 2]>,
    validated_f_interval: IntervalPair,
    validated_fx_interval: IntervalPair,
    validated_fy_interval: IntervalPair,
    validated_hessian_box: TaylorHessianBox,
    validated_derivative_signed: bool,
    strict_collar_candidate: bool,
    best_projection_axis: String,
    collar_reason: String,
    contracted_interval: Option<IntervalPair>,
    length_upper_candidate: Option<f64>,
    width_ratio_vs_prior_taylor_fx: Option<f64>,
    width_ratio_vs_prior_taylor_fy: Option<f64>,
}

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    out_dir: PathBuf,
    piece_limit: usize,
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

fn root_intervals() -> Vec<CInterval> {
    let u0_pair = split_interval(U0_PARENT, SUBDIVISION)[SUB_I];
    let u1_pair = split_interval(U1_PARENT, SUBDIVISION)[SUB_J];
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

fn ci_point(x: f64, y: f64) -> CInterval {
    CInterval {
        re: iv(x),
        im: iv(y),
    }
}

fn ci_box(x: [f64; 2], y: [f64; 2]) -> CInterval {
    CInterval {
        re: iwrap(x[0], x[1]),
        im: iwrap(y[0], y[1]),
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

fn abs_sq(z: CInterval) -> Interval {
    interval_square(z.re) + interval_square(z.im)
}

fn eval_p_p1_p2(z: CInterval, roots: &[CInterval]) -> (CInterval, CInterval, CInterval) {
    let mut p = ci_one();
    let mut p1 = ci_zero();
    let mut p2 = ci_zero();
    for root in roots {
        let factor = ci_sub(z, *root);
        let next_p2 = ci_add(ci_mul(p2, factor), ci_scale(p1, 2.0));
        let next_p1 = ci_add(ci_mul(p1, factor), p);
        let next_p = ci_mul(p, factor);
        p = next_p;
        p1 = next_p1;
        p2 = next_p2;
    }
    (p, p1, p2)
}

fn gradient_intervals(p: CInterval, p1: CInterval) -> (Interval, Interval) {
    let q = ci_mul(p1, ci_conj(p));
    (q.re * iv(2.0), q.im * iv(-2.0))
}

fn hessian_intervals(p: CInterval, p1: CInterval, p2: CInterval) -> (Interval, Interval, Interval) {
    let r = ci_mul(p2, ci_conj(p));
    let p1_abs_sq = abs_sq(p1);
    let fxx = r.re * iv(2.0) + p1_abs_sq * iv(2.0);
    let fxy = r.im * iv(-2.0);
    let fyy = r.re * iv(-2.0) + p1_abs_sq * iv(2.0);
    (fxx, fxy, fyy)
}

fn f_interval_at(z: CInterval, roots: &[CInterval]) -> Interval {
    let (p, _, _) = eval_p_p1_p2(z, roots);
    abs_sq(p) - iv(1.0)
}

fn radius(pair: [f64; 2]) -> f64 {
    0.5 * (pair[1] - pair[0]).abs()
}

fn midpoint(pair: [f64; 2]) -> f64 {
    0.5 * (pair[0] + pair[1])
}

fn sign_stable(x: Interval) -> bool {
    x.inf() > 0.0 || x.sup() < 0.0
}

fn interval_width_pair(pair: [f64; 2]) -> Option<f64> {
    let w = pair[1] - pair[0];
    if w > 0.0 && w.is_finite() {
        Some(w)
    } else {
        None
    }
}

fn interval_width(x: Interval) -> Option<f64> {
    let w = x.sup() - x.inf();
    if w > 0.0 && w.is_finite() {
        Some(w)
    } else {
        None
    }
}

fn width_ratio(prior: Option<[f64; 2]>, validated: Interval) -> Option<f64> {
    let old = interval_width_pair(prior?)?;
    let new = interval_width(validated)?;
    Some(old / new)
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

fn optional_source_interval(piece: &serde_json::Value, field: &str) -> Option<[f64; 2]> {
    source_interval(piece, field).ok()
}

fn taylor_enclosures(
    x: [f64; 2],
    y: [f64; 2],
    roots: &[CInterval],
) -> (Interval, Interval, Interval, Interval, Interval, Interval) {
    let cx = midpoint(x);
    let cy = midpoint(y);
    let rx = radius(x);
    let ry = radius(y);
    let z0 = ci_point(cx, cy);
    let (p0, p10, _) = eval_p_p1_p2(z0, roots);
    let f0 = abs_sq(p0) - iv(1.0);
    let (fx0, fy0) = gradient_intervals(p0, p10);
    let (p, p1, p2) = eval_p_p1_p2(ci_box(x, y), roots);
    let (fxx, fxy, fyy) = hessian_intervals(p, p1, p2);
    let dx = iwrap(-rx, rx);
    let dy = iwrap(-ry, ry);
    let f = f0
        + fx0 * dx
        + fy0 * dy
        + iv(0.5) * fxx * interval_square(dx)
        + fxy * dx * dy
        + iv(0.5) * fyy * interval_square(dy);
    let fx = fx0 + fxx * dx + fxy * dy;
    let fy = fy0 + fxy * dx + fyy * dy;
    (f, fx, fy, fxx, fxy, fyy)
}

fn strict_subset(inner: Interval, outer: [f64; 2]) -> bool {
    inner.inf().is_finite()
        && inner.sup().is_finite()
        && inner.inf() > outer[0]
        && inner.sup() < outer[1]
}

fn krawczyk_y(x: [f64; 2], y: [f64; 2], fy: Interval, roots: &[CInterval]) -> Option<Interval> {
    if !sign_stable(fy) {
        return None;
    }
    let y0 = midpoint(y);
    let fy_mid = 0.5 * (fy.inf() + fy.sup());
    if !fy_mid.is_finite() || fy_mid == 0.0 {
        return None;
    }
    let f0 = f_interval_at(ci_box(x, [y0, y0]), roots);
    let y_delta = iwrap(y[0] - y0, y[1] - y0);
    let c = iv(1.0 / fy_mid);
    let k = iv(y0) - f0 * c + (iv(1.0) - fy * c) * y_delta;
    if strict_subset(k, y) {
        Some(k)
    } else {
        None
    }
}

fn krawczyk_x(x: [f64; 2], y: [f64; 2], fx: Interval, roots: &[CInterval]) -> Option<Interval> {
    if !sign_stable(fx) {
        return None;
    }
    let x0 = midpoint(x);
    let fx_mid = 0.5 * (fx.inf() + fx.sup());
    if !fx_mid.is_finite() || fx_mid == 0.0 {
        return None;
    }
    let f0 = f_interval_at(ci_box([x0, x0], y), roots);
    let x_delta = iwrap(x[0] - x0, x[1] - x0);
    let c = iv(1.0 / fx_mid);
    let k = iv(x0) - f0 * c + (iv(1.0) - fx * c) * x_delta;
    if strict_subset(k, x) {
        Some(k)
    } else {
        None
    }
}

fn analyze_piece(
    idx: usize,
    piece: &serde_json::Value,
    roots: &[CInterval],
) -> Result<PieceResult, String> {
    let x = source_interval(piece, "x_interval")?;
    let y = source_interval(piece, "y_interval")?;
    let source_fx = source_interval(piece, "source_fx_interval")?;
    let source_fy = source_interval(piece, "source_fy_interval")?;
    let prior_taylor_fx = optional_source_interval(piece, "taylor_fx_interval");
    let prior_taylor_fy = optional_source_interval(piece, "taylor_fy_interval");
    let (f, fx, fy, fxx, fxy, fyy) = taylor_enclosures(x, y, roots);
    let y_candidate = krawczyk_y(x, y, fy, roots).map(|k| {
        let slope = interval_abs_upper(fx) / interval_abs_lower(fy);
        let len = (x[1] - x[0]) * (1.0 + slope * slope).sqrt();
        ("y_as_function_of_x".to_string(), k, slope, len)
    });
    let x_candidate = krawczyk_x(x, y, fx, roots).map(|k| {
        let slope = interval_abs_upper(fy) / interval_abs_lower(fx);
        let len = (y[1] - y[0]) * (1.0 + slope * slope).sqrt();
        ("x_as_function_of_y".to_string(), k, slope, len)
    });
    let chosen = match (y_candidate, x_candidate) {
        (Some(a), Some(b)) => {
            if a.3 <= b.3 {
                Some(a)
            } else {
                Some(b)
            }
        }
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    };
    let (strict_collar_candidate, best_axis, reason, contracted, length) = if let Some(c) = chosen {
        (
            true,
            c.0,
            "validated_interval_krawczyk_strict_subset".to_string(),
            Some(ipair(c.1)),
            Some(c.3),
        )
    } else if !sign_stable(fx) && !sign_stable(fy) {
        (
            false,
            "none".to_string(),
            "validated_derivative_intervals_not_sign_stable".to_string(),
            None,
            None,
        )
    } else {
        let axis = if interval_abs_lower(fx) >= interval_abs_lower(fy) {
            "x_as_function_of_y"
        } else {
            "y_as_function_of_x"
        };
        (
            false,
            axis.to_string(),
            "validated_krawczyk_not_strict_subset".to_string(),
            None,
            None,
        )
    };

    Ok(PieceResult {
        source_index: idx,
        source_tube_index: piece
            .get("source_tube_index")
            .and_then(|v| v.as_u64())
            .map(|v| v as usize),
        x_interval: x,
        y_interval: y,
        source_reason: piece
            .get("root_collar_reason")
            .and_then(|v| v.as_str())
            .unwrap_or("unknown")
            .to_string(),
        source_fx_interval: source_fx,
        source_fy_interval: source_fy,
        prior_taylor_fx_interval: prior_taylor_fx,
        prior_taylor_fy_interval: prior_taylor_fy,
        validated_f_interval: ipair(f),
        validated_fx_interval: ipair(fx),
        validated_fy_interval: ipair(fy),
        validated_hessian_box: TaylorHessianBox {
            fxx: ipair(fxx),
            fxy: ipair(fxy),
            fyy: ipair(fyy),
        },
        validated_derivative_signed: sign_stable(fx) || sign_stable(fy),
        strict_collar_candidate,
        best_projection_axis: best_axis,
        collar_reason: reason,
        contracted_interval: contracted,
        length_upper_candidate: length,
        width_ratio_vs_prior_taylor_fx: width_ratio(prior_taylor_fx, fx),
        width_ratio_vs_prior_taylor_fy: width_ratio(prior_taylor_fy, fy),
    })
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01_RESULTS.json"),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-VALIDATED-TAYLOR-ROOT-COLLAR-HARD-CELL-20260506-01"),
        piece_limit: 10,
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
        "# EHP114 n=14 Validated Taylor Root-Collar Pilot\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Processed pieces: `{}`\n\
- Validated derivative signed pieces: `{}`\n\
- Strict collar candidates: `{}`\n\
- Remaining non-candidate pieces: `{}`\n\
- Source accepted length: `{}`\n\
- Candidate length upper sum: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap if candidates promoted: `{}`\n\n\
## Interpretation\n\n\
This run carries the root-affine uncertainty explicitly through interval roots \
and interval Taylor/Hessian enclosures. It is stronger than the midpoint \
Taylor/affine diagnostic, but it is still a local pilot over a small prefix of \
pieces. A strict collar candidate requires an interval Newton/Krawczyk image \
to be a strict subset of the original dependent interval.\n\n\
## Next Blocker\n\n\
{}\n\n\
## Claim Ceiling\n\n\
{}\n",
        EXPERIMENT_ID,
        SOURCE_EXPERIMENT_ID,
        result["status"],
        result["processed_piece_count"],
        result["validated_derivative_signed_count"],
        result["strict_collar_candidate_count"],
        result["remaining_non_candidate_piece_count"],
        result["source_accepted_length_upper"],
        result["candidate_length_upper_sum"],
        result["exact_length_cap"],
        result["margin_to_cap_if_candidates_promoted"],
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
    let pieces_json = source_json["piece_summaries"]
        .as_array()
        .ok_or_else(|| "source missing piece_summaries".to_string())?;
    let source_accepted_length = source_json["source_accepted_length_upper"]
        .as_f64()
        .ok_or_else(|| "source missing source_accepted_length_upper".to_string())?;
    let cap = source_json["exact_length_cap"]
        .as_f64()
        .ok_or_else(|| "source missing exact_length_cap".to_string())?;
    let limit = if cfg.piece_limit == 0 {
        pieces_json.len()
    } else {
        cfg.piece_limit.min(pieces_json.len())
    };
    let roots = root_intervals();
    let mut pieces = Vec::with_capacity(limit);
    for (idx, piece) in pieces_json.iter().take(limit).enumerate() {
        pieces.push(analyze_piece(idx, piece, &roots)?);
    }
    let validated_derivative_signed_count = pieces
        .iter()
        .filter(|p| p.validated_derivative_signed)
        .count();
    let strict_collar_candidate_count = pieces.iter().filter(|p| p.strict_collar_candidate).count();
    let remaining_non_candidate_piece_count = limit.saturating_sub(strict_collar_candidate_count);
    let candidate_length_upper_sum: f64 =
        pieces.iter().filter_map(|p| p.length_upper_candidate).sum();
    let diagnostic_total_length_if_candidates_promoted =
        source_accepted_length + candidate_length_upper_sum;
    let margin_to_cap = cap - diagnostic_total_length_if_candidates_promoted;
    let status = if strict_collar_candidate_count > 0 {
        "VALIDATED_TAYLOR_ROOT_COLLAR_PILOT_CANDIDATES_NOT_GLOBAL_PROOF"
    } else if validated_derivative_signed_count > 0 {
        "VALIDATED_TAYLOR_ROOT_COLLAR_DERIVATIVE_SIGNED_NO_COLLAR"
    } else {
        "VALIDATED_TAYLOR_ROOT_COLLAR_PILOT_FAILS"
    };
    let next_blocker = if strict_collar_candidate_count > 0 {
        "At least one strict interval collar candidate exists under root-affine uncertainty. Promote the candidate to a theorem-shaped certificate and then expand cautiously."
    } else if validated_derivative_signed_count > 0 {
        "Root-affine Taylor enclosures keep derivative signs but Krawczyk is not a strict subset. The next proof-facing step is a sharper analytic collar lemma or higher-order Taylor model with explicit third-derivative remainder."
    } else {
        "Root-affine Taylor intervals lose derivative sign. The next route is analytic root-collar structure, not more first-order validated numerics."
    };

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
            "root_model": "full root-affine interval subcell",
            "taylor_model": "interval center value plus interval Hessian bounds over each piece",
            "collar_rule": "strict interval Newton/Krawczyk subset in x=f(y) or y=f(x)"
        },
        "source_piece_count": pieces_json.len(),
        "processed_piece_count": limit,
        "source_unprocessed_piece_count": pieces_json.len().saturating_sub(limit),
        "source_accepted_length_upper": source_accepted_length,
        "candidate_length_upper_sum": candidate_length_upper_sum,
        "diagnostic_total_length_if_candidates_promoted": diagnostic_total_length_if_candidates_promoted,
        "exact_length_cap": cap,
        "margin_to_cap_if_candidates_promoted": margin_to_cap,
        "validated_derivative_signed_count": validated_derivative_signed_count,
        "strict_collar_candidate_count": strict_collar_candidate_count,
        "remaining_non_candidate_piece_count": remaining_non_candidate_piece_count,
        "piece_results": pieces,
        "claim_ceiling": "Validated Taylor root-collar pilot only. Not a proof of Erdos #114, not a global n=14 proof, and not an exact lemniscate-length certificate.",
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
            "validated_derivative_signed_count": validated_derivative_signed_count,
            "strict_collar_candidate_count": strict_collar_candidate_count,
            "remaining_non_candidate_piece_count": remaining_non_candidate_piece_count,
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
