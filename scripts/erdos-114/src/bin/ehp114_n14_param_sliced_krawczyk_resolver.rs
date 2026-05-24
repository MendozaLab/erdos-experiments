//! EHP #114 n=14 hard-cell parameter-sliced Krawczyk resolver.
//!
//! This binary consumes the latest Krawczyk tube artifact and retries the
//! remaining branch pieces after slicing the root-affine parameter cell. The
//! key question is whether the chart blocker is root-parameter interval width,
//! not geometric tube width.
//!
//! The output is a local certificate diagnostic only. It is not a proof of
//! Erdős #114 and not a global n=14 proof.

#![allow(dead_code)]

use inari::{interval, Interval};
use rayon::prelude::*;
use serde::Serialize;
use serde_json::json;
use std::collections::HashSet;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const EXPERIMENT_ID: &str = "EXP-MATH-EHP114-N14-PARAM-SLICED-KRAWCZYK-PILOT-HARD-CELL-20260506-01";
const SOURCE_EXPERIMENT_ID: &str = "EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01";
const DEGREE: usize = 14;
const EPS: f64 = 0.1;
const EXTENT: f64 = 3.0;
const RES: usize = 220;
const SUBDIVISION: usize = 8;
const SUB_I: usize = 6;
const SUB_J: usize = 4;
const LSTAR_LOWER: f64 = 30.852910841548532;
const TARGET: f64 = 10.180114778928864;
const EXACT_LENGTH_CAP: f64 = 20.672796062619668;
const PREVIOUS_MARCHING_LENGTH_UPPER: f64 = 18.110795101362747;
const PREVIOUS_NORMAL_DRIFT_ERROR_SUM: f64 = 2.4736484172357147;
const PREVIOUS_NORMAL_DRIFT_BUDGET: f64 = 2.5620009612530126;
const PREVIOUS_RELATIVE_ERROR_TO_BUDGET: f64 = 2.730675196761203;
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
#[serde(rename_all = "snake_case")]
enum ProjectionAxis {
    YAsFunctionOfX,
    XAsFunctionOfY,
}

#[derive(Clone, Debug, Serialize)]
struct ValidatedPatch {
    ix: usize,
    iy: usize,
    sx: usize,
    sy: usize,
    ownership_key: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    projection_axis: ProjectionAxis,
    independent_width: f64,
    denominator_component: String,
    denominator_abs_lower: f64,
    numerator_abs_upper: f64,
    slope_abs_upper: f64,
    patch_length_upper: f64,
    f_interval: [f64; 2],
    p_abs_upper: f64,
    p_prime_abs_upper: f64,
    p_second_abs_upper: f64,
    normal_drift_error_diagnostic: f64,
}

#[derive(Clone, Debug, Serialize)]
struct UnresolvedPatch {
    ix: usize,
    iy: usize,
    sx: usize,
    sy: usize,
    ownership_key: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    reason: String,
    f_interval: [f64; 2],
    fx_interval: [f64; 2],
    fy_interval: [f64; 2],
    fx_abs_lower: f64,
    fy_abs_lower: f64,
    fx_abs_upper: f64,
    fy_abs_upper: f64,
}

#[derive(Clone, Copy, Debug)]
struct RunConfig {
    z_subdivision: usize,
    local_bound_subdivision: usize,
}

#[derive(Default)]
struct ValidatedRun {
    patches: Vec<ValidatedPatch>,
    unresolved: Vec<UnresolvedPatch>,
    excluded_box_count: usize,
}

#[derive(Clone, Debug, Serialize)]
struct SlabBranch {
    ix: usize,
    group_index: usize,
    ownership_key: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    y_cells: [usize; 2],
    projection_axis: String,
    denominator_component: String,
    denominator_abs_lower: f64,
    numerator_abs_upper: f64,
    slope_abs_upper: f64,
    length_upper: f64,
    bottom_f_interval: [f64; 2],
    top_f_interval: [f64; 2],
    left_f_interval: [f64; 2],
    right_f_interval: [f64; 2],
    fx_interval: [f64; 2],
    fy_interval: [f64; 2],
    p_abs_upper: f64,
    p_prime_abs_upper: f64,
    p_second_abs_upper: f64,
}

#[derive(Clone, Debug, Serialize)]
struct UnresolvedSlabBranch {
    ix: usize,
    group_index: usize,
    ownership_key: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    y_cells: [usize; 2],
    reason: String,
    bottom_f_interval: [f64; 2],
    top_f_interval: [f64; 2],
    left_f_interval: [f64; 2],
    right_f_interval: [f64; 2],
    fx_interval: [f64; 2],
    fy_interval: [f64; 2],
    fx_abs_lower: f64,
    fy_abs_lower: f64,
    fx_sign: String,
    fy_sign: String,
    candidate_cell_count: usize,
}

#[derive(Default)]
struct SlabRun {
    branches: Vec<SlabBranch>,
    unresolved: Vec<UnresolvedSlabBranch>,
    excluded_cell_count: usize,
    candidate_cell_count: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SignState {
    Positive,
    Negative,
    Mixed,
    ContainsZero,
}

fn iv(x: f64) -> Interval {
    interval!(x, x).unwrap()
}

fn iwrap(lo: f64, hi: f64) -> Interval {
    interval!(lo, hi).unwrap()
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

fn worst_subcell_intervals() -> ([f64; 2], [f64; 2]) {
    let u0 = split_interval(U0_PARENT, SUBDIVISION)[SUB_I];
    let u1 = split_interval(U1_PARENT, SUBDIVISION)[SUB_J];
    (u0, u1)
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
    let (u0_pair, u1_pair) = worst_subcell_intervals();
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

fn root_intervals_for_param_tile(
    tile_i: usize,
    tile_j: usize,
    param_subdivision: usize,
) -> Vec<CInterval> {
    let (u0_pair, u1_pair) = worst_subcell_intervals();
    let u0_tile = split_interval(u0_pair, param_subdivision)[tile_i];
    let u1_tile = split_interval(u1_pair, param_subdivision)[tile_j];
    let a = iwrap(u0_tile[0], u0_tile[1]);
    let b = iwrap(u1_tile[0], u1_tile[1]);
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

fn abs_sq(z: CInterval) -> Interval {
    interval_square(z.re) + interval_square(z.im)
}

fn abs_upper(z: CInterval) -> f64 {
    abs_sq(z).sup().max(0.0).sqrt()
}

fn abs_sq_minus_one(z: CInterval) -> Interval {
    abs_sq(z) - iv(1.0)
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

fn sign_state(x: Interval) -> SignState {
    if x.inf() > 0.0 {
        SignState::Positive
    } else if x.sup() < 0.0 {
        SignState::Negative
    } else {
        SignState::ContainsZero
    }
}

fn combine_sign_state(acc: SignState, next: SignState) -> SignState {
    match (acc, next) {
        (SignState::ContainsZero, _) | (_, SignState::ContainsZero) => SignState::ContainsZero,
        (SignState::Mixed, _) | (_, SignState::Mixed) => SignState::Mixed,
        (SignState::Positive, SignState::Positive) => SignState::Positive,
        (SignState::Negative, SignState::Negative) => SignState::Negative,
        _ => SignState::Mixed,
    }
}

fn sign_state_label(s: SignState) -> &'static str {
    match s {
        SignState::Positive => "positive",
        SignState::Negative => "negative",
        SignState::Mixed => "mixed",
        SignState::ContainsZero => "contains_zero",
    }
}

fn f_may_contain_level(f: Interval) -> bool {
    f.inf() <= 0.0 && f.sup() >= 0.0
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
    // F = |p|^2 - 1. For q = p' * conjugate(p):
    // F_x = 2 Re(q), F_y = -2 Im(q).
    let q = ci_mul(p1, ci_conj(p));
    (q.re * iv(2.0), q.im * iv(-2.0))
}

fn local_box_stats(
    x: [f64; 2],
    y: [f64; 2],
    roots: &[CInterval],
    subdivision: usize,
) -> (
    bool,
    Interval,
    Interval,
    Interval,
    f64,
    f64,
    f64,
    f64,
    f64,
    SignState,
    SignState,
) {
    let direct_z = ci_box(x, y);
    let (direct_p, _, _) = eval_p_p1_p2(direct_z, roots);
    let direct_f = abs_sq_minus_one(direct_p);
    if !f_may_contain_level(direct_f) {
        return (
            false,
            direct_f,
            iv(0.0),
            iv(0.0),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            SignState::ContainsZero,
            SignState::ContainsZero,
        );
    }

    let n = subdivision.max(1);
    let x_step = (x[1] - x[0]) / (n as f64);
    let y_step = (y[1] - y[0]) / (n as f64);
    let mut any_level = false;
    let mut f_union = direct_f;
    let mut fx_union: Option<Interval> = None;
    let mut fy_union: Option<Interval> = None;
    let mut fx_abs_lower = f64::INFINITY;
    let mut fy_abs_lower = f64::INFINITY;
    let mut fx_abs_upper = 0.0_f64;
    let mut fy_abs_upper = 0.0_f64;
    let mut p_abs_upper = 0.0_f64;
    let mut p_prime_abs_upper = 0.0_f64;
    let mut p_second_abs_upper = 0.0_f64;
    let mut fx_sign = SignState::Positive;
    let mut fy_sign = SignState::Positive;
    let mut first_active = true;

    for sx in 0..n {
        for sy in 0..n {
            let xs = [
                x[0] + (sx as f64) * x_step,
                x[0] + ((sx + 1) as f64) * x_step,
            ];
            let ys = [
                y[0] + (sy as f64) * y_step,
                y[0] + ((sy + 1) as f64) * y_step,
            ];
            let (p, p1, p2) = eval_p_p1_p2(ci_box(xs, ys), roots);
            let f = abs_sq_minus_one(p);
            if !f_may_contain_level(f) {
                continue;
            }
            any_level = true;
            f_union = iwrap(f_union.inf().min(f.inf()), f_union.sup().max(f.sup()));
            let (fx, fy) = gradient_intervals(p, p1);
            fx_union = Some(match fx_union {
                Some(old) => iwrap(old.inf().min(fx.inf()), old.sup().max(fx.sup())),
                None => fx,
            });
            fy_union = Some(match fy_union {
                Some(old) => iwrap(old.inf().min(fy.inf()), old.sup().max(fy.sup())),
                None => fy,
            });
            fx_abs_lower = fx_abs_lower.min(interval_abs_lower(fx));
            fy_abs_lower = fy_abs_lower.min(interval_abs_lower(fy));
            fx_abs_upper = fx_abs_upper.max(interval_abs_upper(fx));
            fy_abs_upper = fy_abs_upper.max(interval_abs_upper(fy));
            p_abs_upper = p_abs_upper.max(abs_upper(p));
            p_prime_abs_upper = p_prime_abs_upper.max(abs_upper(p1));
            p_second_abs_upper = p_second_abs_upper.max(abs_upper(p2));
            if first_active {
                fx_sign = sign_state(fx);
                fy_sign = sign_state(fy);
                first_active = false;
            } else {
                fx_sign = combine_sign_state(fx_sign, sign_state(fx));
                fy_sign = combine_sign_state(fy_sign, sign_state(fy));
            }
        }
    }

    if !any_level {
        return (
            false,
            direct_f,
            iv(0.0),
            iv(0.0),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            SignState::ContainsZero,
            SignState::ContainsZero,
        );
    }

    (
        true,
        f_union,
        fx_union.unwrap_or(iv(0.0)),
        fy_union.unwrap_or(iv(0.0)),
        fx_abs_lower,
        fy_abs_lower,
        fx_abs_upper,
        fy_abs_upper,
        p_abs_upper.max(p_prime_abs_upper).max(p_second_abs_upper),
        fx_sign,
        fy_sign,
    )
}

fn p_abs_stats(
    x: [f64; 2],
    y: [f64; 2],
    roots: &[CInterval],
    subdivision: usize,
) -> (f64, f64, f64) {
    let n = subdivision.max(1);
    let x_step = (x[1] - x[0]) / (n as f64);
    let y_step = (y[1] - y[0]) / (n as f64);
    let mut p_abs_upper = 0.0_f64;
    let mut p_prime_abs_upper = 0.0_f64;
    let mut p_second_abs_upper = 0.0_f64;
    for sx in 0..n {
        for sy in 0..n {
            let xs = [
                x[0] + (sx as f64) * x_step,
                x[0] + ((sx + 1) as f64) * x_step,
            ];
            let ys = [
                y[0] + (sy as f64) * y_step,
                y[0] + ((sy + 1) as f64) * y_step,
            ];
            let (p, p1, p2) = eval_p_p1_p2(ci_box(xs, ys), roots);
            p_abs_upper = p_abs_upper.max(abs_upper(p));
            p_prime_abs_upper = p_prime_abs_upper.max(abs_upper(p1));
            p_second_abs_upper = p_second_abs_upper.max(abs_upper(p2));
        }
    }
    (p_abs_upper, p_prime_abs_upper, p_second_abs_upper)
}

fn patch_for_box(
    ix: usize,
    iy: usize,
    sx: usize,
    sy: usize,
    x: [f64; 2],
    y: [f64; 2],
    roots: &[CInterval],
    local_bound_subdivision: usize,
) -> Result<Option<ValidatedPatch>, UnresolvedPatch> {
    let (
        active,
        f,
        fx,
        fy,
        fx_abs_lower,
        fy_abs_lower,
        fx_abs_upper,
        fy_abs_upper,
        _combined_abs_upper,
        fx_sign,
        fy_sign,
    ) = local_box_stats(x, y, roots, local_bound_subdivision);

    if !active {
        return Ok(None);
    }

    let ownership_key = format!("{ix}:{iy}:{sx}:{sy}");
    let x_width = x[1] - x[0];
    let y_width = y[1] - y[0];
    let y_graph_ok = fy_abs_lower > 0.0
        && fy_abs_lower.is_finite()
        && matches!(fy_sign, SignState::Positive | SignState::Negative);
    let x_graph_ok = fx_abs_lower > 0.0
        && fx_abs_lower.is_finite()
        && matches!(fx_sign, SignState::Positive | SignState::Negative);

    if !y_graph_ok && !x_graph_ok {
        return Err(UnresolvedPatch {
            ix,
            iy,
            sx,
            sy,
            ownership_key,
            x_interval: x,
            y_interval: y,
            reason: "no_projection_axis_with_component_derivative_bounded_away_from_zero"
                .to_string(),
            f_interval: [f.inf(), f.sup()],
            fx_interval: [fx.inf(), fx.sup()],
            fy_interval: [fy.inf(), fy.sup()],
            fx_abs_lower,
            fy_abs_lower,
            fx_abs_upper,
            fy_abs_upper,
        });
    }

    let y_graph_length = if y_graph_ok {
        let slope = fx_abs_upper / fy_abs_lower;
        Some((
            ProjectionAxis::YAsFunctionOfX,
            x_width,
            "Fy".to_string(),
            fy_abs_lower,
            fx_abs_upper,
            slope,
            x_width * (1.0 + slope * slope).sqrt(),
        ))
    } else {
        None
    };
    let x_graph_length = if x_graph_ok {
        let slope = fy_abs_upper / fx_abs_lower;
        Some((
            ProjectionAxis::XAsFunctionOfY,
            y_width,
            "Fx".to_string(),
            fx_abs_lower,
            fy_abs_upper,
            slope,
            y_width * (1.0 + slope * slope).sqrt(),
        ))
    } else {
        None
    };

    let chosen = match (y_graph_length, x_graph_length) {
        (Some(a), Some(b)) => {
            if a.6 <= b.6 {
                a
            } else {
                b
            }
        }
        (Some(a), None) => a,
        (None, Some(b)) => b,
        (None, None) => unreachable!("handled unresolved axes above"),
    };

    let (p_abs_upper, p_prime_abs_upper, p_second_abs_upper) =
        p_abs_stats(x, y, roots, local_bound_subdivision);
    let gradient_lower = 2.0 * chosen.3;
    let hessian_upper =
        4.0 * (p_prime_abs_upper * p_prime_abs_upper + p_abs_upper * p_second_abs_upper);
    let condition_ratio = if gradient_lower > 0.0 {
        hessian_upper / gradient_lower
    } else {
        0.0
    };
    let patch_diameter = (x_width * x_width + y_width * y_width).sqrt();
    let normal_drift_error_diagnostic = condition_ratio * patch_diameter * patch_diameter;

    Ok(Some(ValidatedPatch {
        ix,
        iy,
        sx,
        sy,
        ownership_key,
        x_interval: x,
        y_interval: y,
        projection_axis: chosen.0,
        independent_width: chosen.1,
        denominator_component: chosen.2,
        denominator_abs_lower: chosen.3,
        numerator_abs_upper: chosen.4,
        slope_abs_upper: chosen.5,
        patch_length_upper: chosen.6,
        f_interval: [f.inf(), f.sup()],
        p_abs_upper,
        p_prime_abs_upper,
        p_second_abs_upper,
        normal_drift_error_diagnostic,
    }))
}

fn run_ix(
    ix: usize,
    roots: &[CInterval],
    grid_step: f64,
    z_subdivision: usize,
    local_bound_subdivision: usize,
) -> ValidatedRun {
    let mut out = ValidatedRun::default();
    let sub_step = grid_step / (z_subdivision as f64);
    for iy in 0..RES {
        let x0 = -EXTENT + (ix as f64) * grid_step;
        let y0 = -EXTENT + (iy as f64) * grid_step;
        for sx in 0..z_subdivision {
            for sy in 0..z_subdivision {
                let x = [
                    x0 + (sx as f64) * sub_step,
                    x0 + ((sx + 1) as f64) * sub_step,
                ];
                let y = [
                    y0 + (sy as f64) * sub_step,
                    y0 + ((sy + 1) as f64) * sub_step,
                ];
                match patch_for_box(ix, iy, sx, sy, x, y, roots, local_bound_subdivision) {
                    Ok(Some(patch)) => out.patches.push(patch),
                    Ok(None) => out.excluded_box_count += 1,
                    Err(unresolved) => out.unresolved.push(unresolved),
                }
            }
        }
    }
    out
}

fn run_validated_length(
    roots: &[CInterval],
    grid_step: f64,
    z_subdivision: usize,
    local_bound_subdivision: usize,
) -> ValidatedRun {
    let runs: Vec<ValidatedRun> = (0..RES)
        .into_par_iter()
        .map(|ix| run_ix(ix, roots, grid_step, z_subdivision, local_bound_subdivision))
        .collect();

    let mut out = ValidatedRun::default();
    for mut run in runs {
        out.patches.append(&mut run.patches);
        out.unresolved.append(&mut run.unresolved);
        out.excluded_box_count += run.excluded_box_count;
    }
    out.patches.sort_by_key(|p| (p.ix, p.iy, p.sx, p.sy));
    out.unresolved.sort_by_key(|p| (p.ix, p.iy, p.sx, p.sy));
    out
}

fn f_interval_for_box(x: [f64; 2], y: [f64; 2], roots: &[CInterval]) -> Interval {
    let (p, _, _) = eval_p_p1_p2(ci_box(x, y), roots);
    abs_sq_minus_one(p)
}

fn opposite_strict_signs(a: Interval, b: Interval) -> bool {
    (a.sup() < 0.0 && b.inf() > 0.0) || (a.inf() > 0.0 && b.sup() < 0.0)
}

fn certify_slab_group(
    ix: usize,
    group_index: usize,
    group_start: usize,
    group_end: usize,
    x: [f64; 2],
    y_step: f64,
    roots: &[CInterval],
    local_bound_subdivision: usize,
) -> Result<SlabBranch, UnresolvedSlabBranch> {
    let y = [
        -EXTENT + (group_start as f64) * y_step,
        -EXTENT + ((group_end + 1) as f64) * y_step,
    ];
    let ownership_key = format!("{ix}:{group_index}");
    let bottom = tighten_interval(
        f_interval_for_box(x, [y[0], y[0]], roots),
        bernstein_f_horizontal_edge(x, y[0], roots),
    );
    let top = tighten_interval(
        f_interval_for_box(x, [y[1], y[1]], roots),
        bernstein_f_horizontal_edge(x, y[1], roots),
    );
    let left = tighten_interval(
        f_interval_for_box([x[0], x[0]], y, roots),
        bernstein_f_vertical_edge(x[0], y, roots),
    );
    let right = tighten_interval(
        f_interval_for_box([x[1], x[1]], y, roots),
        bernstein_f_vertical_edge(x[1], y, roots),
    );
    let (
        _active,
        _f,
        fx,
        fy,
        fx_abs_lower,
        fy_abs_lower,
        fx_abs_upper,
        fy_abs_upper,
        _combined_abs_upper,
        fx_sign,
        fy_sign,
    ) = local_box_stats(x, y, roots, local_bound_subdivision);

    let unresolved = |reason: &str| UnresolvedSlabBranch {
        ix,
        group_index,
        ownership_key: ownership_key.clone(),
        x_interval: x,
        y_interval: y,
        y_cells: [group_start, group_end],
        reason: reason.to_string(),
        bottom_f_interval: [bottom.inf(), bottom.sup()],
        top_f_interval: [top.inf(), top.sup()],
        left_f_interval: [left.inf(), left.sup()],
        right_f_interval: [right.inf(), right.sup()],
        fx_interval: [fx.inf(), fx.sup()],
        fy_interval: [fy.inf(), fy.sup()],
        fx_abs_lower,
        fy_abs_lower,
        fx_sign: sign_state_label(fx_sign).to_string(),
        fy_sign: sign_state_label(fy_sign).to_string(),
        candidate_cell_count: group_end + 1 - group_start,
    };

    let (p_abs_upper, p_prime_abs_upper, p_second_abs_upper) =
        p_abs_stats(x, y, roots, local_bound_subdivision);
    let x_width = x[1] - x[0];
    let y_width = y[1] - y[0];

    let y_graph = if opposite_strict_signs(bottom, top)
        && fy_abs_lower > 0.0
        && matches!(fy_sign, SignState::Positive | SignState::Negative)
    {
        let slope_abs_upper = fx_abs_upper / fy_abs_lower;
        Some((
            "y_as_function_of_x".to_string(),
            "Fy".to_string(),
            fy_abs_lower,
            fx_abs_upper,
            slope_abs_upper,
            x_width * (1.0 + slope_abs_upper * slope_abs_upper).sqrt(),
        ))
    } else {
        None
    };

    let x_graph = if opposite_strict_signs(left, right)
        && fx_abs_lower > 0.0
        && matches!(fx_sign, SignState::Positive | SignState::Negative)
    {
        let slope_abs_upper = fy_abs_upper / fx_abs_lower;
        Some((
            "x_as_function_of_y".to_string(),
            "Fx".to_string(),
            fx_abs_lower,
            fy_abs_upper,
            slope_abs_upper,
            y_width * (1.0 + slope_abs_upper * slope_abs_upper).sqrt(),
        ))
    } else {
        None
    };

    let chosen = match (y_graph, x_graph) {
        (Some(a), Some(b)) => {
            if a.5 <= b.5 {
                a
            } else {
                b
            }
        }
        (Some(a), None) => a,
        (None, Some(b)) => b,
        (None, None) => return Err(unresolved("projection_certificate_failed")),
    };

    Ok(SlabBranch {
        ix,
        group_index,
        ownership_key,
        x_interval: x,
        y_interval: y,
        y_cells: [group_start, group_end],
        projection_axis: chosen.0,
        denominator_component: chosen.1,
        denominator_abs_lower: chosen.2,
        numerator_abs_upper: chosen.3,
        slope_abs_upper: chosen.4,
        length_upper: chosen.5,
        bottom_f_interval: [bottom.inf(), bottom.sup()],
        top_f_interval: [top.inf(), top.sup()],
        left_f_interval: [left.inf(), left.sup()],
        right_f_interval: [right.inf(), right.sup()],
        fx_interval: [fx.inf(), fx.sup()],
        fy_interval: [fy.inf(), fy.sup()],
        p_abs_upper,
        p_prime_abs_upper,
        p_second_abs_upper,
    })
}

fn run_slab_ix(
    ix: usize,
    roots: &[CInterval],
    x_slab_count: usize,
    y_cell_count: usize,
    local_bound_subdivision: usize,
) -> SlabRun {
    let mut out = SlabRun::default();
    let x_step = 2.0 * EXTENT / (x_slab_count as f64);
    let y_step = 2.0 * EXTENT / (y_cell_count as f64);
    let x = [
        -EXTENT + (ix as f64) * x_step,
        -EXTENT + ((ix + 1) as f64) * x_step,
    ];

    let mut candidate = vec![false; y_cell_count];
    for (iy, slot) in candidate.iter_mut().enumerate() {
        let y = [
            -EXTENT + (iy as f64) * y_step,
            -EXTENT + ((iy + 1) as f64) * y_step,
        ];
        let f = f_interval_for_box(x, y, roots);
        if f_may_contain_level(f) {
            *slot = true;
            out.candidate_cell_count += 1;
        } else {
            out.excluded_cell_count += 1;
        }
    }

    let mut iy = 0usize;
    let mut group_index = 0usize;
    while iy < y_cell_count {
        if !candidate[iy] {
            iy += 1;
            continue;
        }
        let start = iy;
        while iy + 1 < y_cell_count && candidate[iy + 1] {
            iy += 1;
        }
        let end = iy;
        match certify_slab_group(
            ix,
            group_index,
            start,
            end,
            x,
            y_step,
            roots,
            local_bound_subdivision,
        ) {
            Ok(branch) => out.branches.push(branch),
            Err(unresolved) => out.unresolved.push(unresolved),
        }
        group_index += 1;
        iy += 1;
    }
    out
}

fn run_slab_validated_length(
    roots: &[CInterval],
    z_subdivision: usize,
    local_bound_subdivision: usize,
) -> SlabRun {
    let x_slab_count = RES * z_subdivision;
    let y_cell_count = RES * z_subdivision;
    let runs: Vec<SlabRun> = (0..x_slab_count)
        .into_par_iter()
        .map(|ix| {
            run_slab_ix(
                ix,
                roots,
                x_slab_count,
                y_cell_count,
                local_bound_subdivision,
            )
        })
        .collect();

    let mut out = SlabRun::default();
    for mut run in runs {
        out.branches.append(&mut run.branches);
        out.unresolved.append(&mut run.unresolved);
        out.excluded_cell_count += run.excluded_cell_count;
        out.candidate_cell_count += run.candidate_cell_count;
    }
    out.branches.sort_by_key(|p| (p.ix, p.group_index));
    out.unresolved.sort_by_key(|p| (p.ix, p.group_index));
    out
}

fn root_radius_upper(roots: &[CInterval]) -> f64 {
    roots
        .iter()
        .map(|r| abs_upper(*r))
        .fold(f64::NEG_INFINITY, f64::max)
}

fn outside_extent_margin(roots: &[CInterval]) -> f64 {
    // If |z| > max_root_radius + 1, then |p(z)| = prod_i |z-r_i| > 1.
    // Since the scanned square contains the disk |z| <= EXTENT, this margin
    // certifies that no lemniscate arc is missed outside the square.
    EXTENT - root_radius_upper(roots) - 1.0
}

fn unix_timestamp_string() -> String {
    match SystemTime::now().duration_since(UNIX_EPOCH) {
        Ok(dur) => format!("{}", dur.as_secs()),
        Err(_) => "0".to_string(),
    }
}

fn sha256_file(path: &PathBuf) -> Result<String, String> {
    let bytes = fs::read(path)
        .map_err(|err| format!("failed to read {} for sha256: {err}", path.display()))?;
    Ok(sha256::digest(bytes))
}

fn parse_args() -> Result<(PathBuf, RunConfig), String> {
    let args: Vec<String> = env::args().collect();
    let mut out_dir = PathBuf::from(
        "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03",
    );
    let mut z_subdivision = 16usize;
    let mut local_bound_subdivision = 6usize;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--outdir" => {
                if i + 1 >= args.len() {
                    return Err("--outdir requires a path".to_string());
                }
                out_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--z-subdivision" => {
                if i + 1 >= args.len() {
                    return Err("--z-subdivision requires a positive integer".to_string());
                }
                z_subdivision = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --z-subdivision: {err}"))?;
                if z_subdivision == 0 {
                    return Err("--z-subdivision must be positive".to_string());
                }
                i += 2;
            }
            "--local-bound-subdivision" => {
                if i + 1 >= args.len() {
                    return Err("--local-bound-subdivision requires a positive integer".to_string());
                }
                local_bound_subdivision = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --local-bound-subdivision: {err}"))?;
                if local_bound_subdivision == 0 {
                    return Err("--local-bound-subdivision must be positive".to_string());
                }
                i += 2;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
    }
    Ok((
        out_dir,
        RunConfig {
            z_subdivision,
            local_bound_subdivision,
        },
    ))
}

fn write_report(result: &serde_json::Value, path: &PathBuf) -> Result<(), String> {
    let status = result["status"].as_str().unwrap_or("UNKNOWN");
    let report = format!(
        "# EHP114 n=14 Slab Root-Isolation Validated-Length Hard-Cell Diagnostic\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Degree: `{}`\n\
- eps: `{}`\n\
- Root-affine subcell: `({}, {})`\n\
- Exact length cap: `{}`\n\
- Outside extent excluded: `{}` (margin `{}`)\n\
- Slab branch count: `{}`\n\
- Candidate y-cells: `{}`\n\
- Excluded y-cells: `{}`\n\
- Unresolved branch count: `{}`\n\
- Ownership duplicate count: `{}`\n\
- Sum branch length upper: `{}`\n\
- Endpoint/overlap tax: `{}`\n\
- Total validated length upper: `{}`\n\
- Margin to cap: `{}`\n\n\
## Interpretation\n\n\
This run demotes marching squares to a diagnostic and attempts a slab-wise \
implicit root-isolation enclosure. Each accepted branch has endpoint sign \
separation and a fixed nonzero `Fy` interval, so it is counted as one owned \
graph segment over its x-slab. The prior `relative error / budget = \
2.730675196761203` is retained only as historical context and is not \
load-bearing for the pass/fail decision.\n\n\
## Claim Ceiling\n\n\
This is a local hard-cell diagnostic. It is not a proof of Erdős #114 and not \
a global n=14 proof.\n\n\
## Next Blocker\n\n\
{}\n",
        EXPERIMENT_ID,
        SOURCE_EXPERIMENT_ID,
        status,
        DEGREE,
        EPS,
        SUB_I,
        SUB_J,
        EXACT_LENGTH_CAP,
        result["domain_certificate"]["outside_extent_excluded"],
        result["domain_certificate"]["outside_extent_margin"],
        result["slab_branch_count"],
        result["candidate_cell_count"],
        result["excluded_cell_count"],
        result["unresolved_branch_count"],
        result["ownership_duplicate_count"],
        result["sum_branch_length_upper"],
        result["endpoint_overlap_tax"],
        result["total_validated_length_upper"],
        result["margin_to_cap"],
        result["next_blocker"]
            .as_str()
            .unwrap_or("continue diagnostics"),
    );
    fs::write(path, report)
        .map_err(|err| format!("failed to write report {}: {err}", path.display()))
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let (out_dir, cfg) = parse_args()?;
    fs::create_dir_all(&out_dir).map_err(|err| {
        format!(
            "failed to create output directory {}: {err}",
            out_dir.display()
        )
    })?;

    let result_path = out_dir.join(format!("{}_RESULTS.json", EXPERIMENT_ID));
    let report_path = out_dir.join(format!("{}_REPORT.md", EXPERIMENT_ID));
    let sha_path = out_dir.join(format!("{}_RESULTS.sha256", EXPERIMENT_ID));
    for path in [&result_path, &report_path, &sha_path] {
        if path.exists() {
            return Err(format!(
                "refusing to overwrite existing artifact: {}",
                path.display()
            ));
        }
    }

    let roots = root_intervals();
    let outside_margin = outside_extent_margin(&roots);
    let outside_excluded = outside_margin > 0.0;
    let x_slab_count = RES * cfg.z_subdivision;
    let y_cell_count = RES * cfg.z_subdivision;
    let x_slab_width = 2.0 * EXTENT / (x_slab_count as f64);
    let y_cell_width = 2.0 * EXTENT / (y_cell_count as f64);
    let slab_cell_count_total = x_slab_count * y_cell_count;
    let run = run_slab_validated_length(&roots, cfg.z_subdivision, cfg.local_bound_subdivision);

    let slab_branch_count = run.branches.len();
    let unresolved_branch_count = run.unresolved.len();
    let excluded_cell_count = run.excluded_cell_count;
    let candidate_cell_count = run.candidate_cell_count;
    let mut seen = HashSet::new();
    let mut ownership_duplicate_count = 0usize;
    for branch in &run.branches {
        if !seen.insert(branch.ownership_key.clone()) {
            ownership_duplicate_count += 1;
        }
    }

    let endpoint_overlap_tax = if ownership_duplicate_count == 0 {
        0.0
    } else {
        f64::INFINITY
    };
    let sum_branch_length_upper: f64 = run.branches.iter().map(|p| p.length_upper).sum();
    let total_validated_length_upper = sum_branch_length_upper + endpoint_overlap_tax;
    let margin_to_cap = EXACT_LENGTH_CAP - total_validated_length_upper;

    let worst_branch = run
        .branches
        .iter()
        .max_by(|a, b| a.length_upper.total_cmp(&b.length_upper))
        .cloned();
    let worst_unresolved = run.unresolved.first().cloned();

    let status = if unresolved_branch_count > 0 {
        "SLAB_FAIL_ROOT_ISOLATION"
    } else if ownership_duplicate_count > 0 {
        "SLAB_FAIL_EXCLUSION"
    } else if !outside_excluded {
        "SLAB_FAIL_EXCLUSION"
    } else if total_validated_length_upper > EXACT_LENGTH_CAP {
        "SLAB_FAIL_BUDGET"
    } else {
        "SLAB_VALIDATED_LENGTH_PASS_NOT_GLOBAL_PROOF"
    };

    let next_blocker = if unresolved_branch_count > 0 {
        "Some candidate y-runs failed endpoint sign separation or fixed nonzero Fy. Need narrower slabs, vertical refinement, rotated charts, or interval Newton isolation."
    } else if ownership_duplicate_count > 0 {
        "Slab branch ownership keys are not unique. Fix seam ownership before interpreting any length sum."
    } else if !outside_excluded {
        "The finite spatial scan does not yet exclude arcs outside the scanned square. Increase EXTENT or prove a sharper outside-domain exclusion."
    } else if total_validated_length_upper > EXACT_LENGTH_CAP {
        "The slab root-isolation bound is certified but too coarse for the exact-length cap. Need narrower slabs, tighter derivative bounds, or a sharper graph integration rule."
    } else {
        "Promote the slab root-isolation certificate into a theorem-shaped packet and then test the remaining selected subcells."
    };

    let result = json!({
        "experiment_id": EXPERIMENT_ID,
        "timestamp_unix": unix_timestamp_string(),
        "source_experiment_id": SOURCE_EXPERIMENT_ID,
        "status": status,
        "degree": DEGREE,
        "eps": EPS,
        "subcell": {
            "sub_i": SUB_I,
            "sub_j": SUB_J,
        },
        "root_affine_intervals": {
            "u0_interval": worst_subcell_intervals().0,
            "u1_interval": worst_subcell_intervals().1,
            "root_radius_upper": root_radius_upper(&roots),
        },
        "domain_certificate": {
            "scanned_square": [[-EXTENT, EXTENT], [-EXTENT, EXTENT]],
            "outside_extent_margin": outside_margin,
            "outside_extent_excluded": outside_excluded,
            "outside_extent_rule": "For monic p with every root in |r| <= R, |z| > R + 1 implies |p(z)| > 1, so the lemniscate |p(z)| = 1 is inside |z| <= R + 1."
        },
        "parameters": {
            "extent": EXTENT,
            "res": RES,
            "z_subdivision": cfg.z_subdivision,
            "local_bound_subdivision": cfg.local_bound_subdivision,
            "x_slab_count": x_slab_count,
            "y_cell_count": y_cell_count,
            "x_slab_width": x_slab_width,
            "y_cell_width": y_cell_width,
            "slab_cell_count_total": slab_cell_count_total,
            "ownership_rule": "half_open_x_slab_partition; branch seams have zero length tax in this diagnostic",
            "projection_rule": "count one y=f(x) graph branch per candidate y-run when endpoint signs are separated and Fy has fixed nonzero sign on the branch tube",
        },
        "length_budget": {
            "lstar_lower": LSTAR_LOWER,
            "target": TARGET,
            "exact_length_cap": EXACT_LENGTH_CAP,
            "previous_marching_length_upper": PREVIOUS_MARCHING_LENGTH_UPPER,
            "previous_normal_drift_error_sum": PREVIOUS_NORMAL_DRIFT_ERROR_SUM,
            "previous_normal_drift_budget": PREVIOUS_NORMAL_DRIFT_BUDGET,
            "previous_relative_error_to_budget": PREVIOUS_RELATIVE_ERROR_TO_BUDGET,
            "relative_error_budget_note": "Prior relative error is retained as a diagnostic only; this run's pass/fail criterion is slab root-isolated branch length against exact_length_cap."
        },
        "slab_branch_count": slab_branch_count,
        "candidate_cell_count": candidate_cell_count,
        "excluded_cell_count": excluded_cell_count,
        "unresolved_branch_count": unresolved_branch_count,
        "ownership_duplicate_count": ownership_duplicate_count,
        "sum_branch_length_upper": sum_branch_length_upper,
        "endpoint_overlap_tax": endpoint_overlap_tax,
        "total_validated_length_upper": total_validated_length_upper,
        "margin_to_cap": margin_to_cap,
        "worst_branch": worst_branch,
        "worst_unresolved": worst_unresolved,
        "branches": run.branches,
        "unresolved_branches": run.unresolved,
        "claim_ceiling": "Slab root-isolation hard-cell diagnostic only. This is not a proof of Erdős #114 and not a global n=14 proof.",
        "next_blocker": next_blocker,
        "elapsed_secs": started.elapsed().as_secs_f64(),
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
            "slab_branch_count": slab_branch_count,
            "candidate_cell_count": candidate_cell_count,
            "excluded_cell_count": excluded_cell_count,
            "unresolved_branch_count": unresolved_branch_count,
            "ownership_duplicate_count": ownership_duplicate_count,
            "total_validated_length_upper": total_validated_length_upper,
            "exact_length_cap": EXACT_LENGTH_CAP,
            "margin_to_cap": margin_to_cap,
            "result": result_path,
            "report": report_path,
            "sha256": sha,
        }))
        .unwrap()
    );
    Ok(())
}

#[derive(Clone, Debug, Serialize)]
struct ResolvedTubePiece {
    source_tube_index: usize,
    depth: usize,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    contracted_interval: [f64; 2],
    certificate_method: String,
    projection_axis: String,
    denominator_component: String,
    denominator_abs_lower: f64,
    numerator_abs_upper: f64,
    slope_abs_upper: f64,
    length_upper: f64,
}

#[derive(Clone, Debug, Serialize)]
struct RemainingTubePiece {
    source_tube_index: usize,
    depth: usize,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    reason: String,
    f_interval: [f64; 2],
    fx_interval: [f64; 2],
    fy_interval: [f64; 2],
    fx_abs_lower: f64,
    fy_abs_lower: f64,
    fx_sign: String,
    fy_sign: String,
}

#[derive(Default)]
struct TubeResolveStats {
    source_tube_count: usize,
    certified_piece_count: usize,
    krawczyk_certified_count: usize,
    excluded_piece_count: usize,
    remaining_pieces: Vec<RemainingTubePiece>,
    resolved_pieces: Vec<ResolvedTubePiece>,
    recursion_node_count: usize,
    max_depth_reached: usize,
}

#[derive(Clone, Debug, Serialize)]
struct ParamTileStats {
    tile_i: usize,
    tile_j: usize,
    source_piece_count: usize,
    certified_piece_count: usize,
    krawczyk_certified_count: usize,
    excluded_piece_count: usize,
    remaining_unresolved_piece_count: usize,
    recursion_node_count: usize,
    max_depth_reached: usize,
    resolved_tube_length_upper: f64,
    total_validated_length_upper: f64,
    margin_to_cap: f64,
    remaining_reason_counts: std::collections::BTreeMap<String, usize>,
    resolved_piece_sample: Vec<ResolvedTubePiece>,
    remaining_piece_sample: Vec<RemainingTubePiece>,
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

fn binom_f64(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut out = 1.0_f64;
    for i in 0..k {
        out *= (n - i) as f64;
        out /= (i + 1) as f64;
    }
    out
}

fn poly_add(a: &[Interval], b: &[Interval]) -> Vec<Interval> {
    let n = a.len().max(b.len());
    let mut out = vec![iv(0.0); n];
    for i in 0..n {
        let ai = a.get(i).copied().unwrap_or_else(|| iv(0.0));
        let bi = b.get(i).copied().unwrap_or_else(|| iv(0.0));
        out[i] = ai + bi;
    }
    out
}

fn poly_sub(a: &[Interval], b: &[Interval]) -> Vec<Interval> {
    let n = a.len().max(b.len());
    let mut out = vec![iv(0.0); n];
    for i in 0..n {
        let ai = a.get(i).copied().unwrap_or_else(|| iv(0.0));
        let bi = b.get(i).copied().unwrap_or_else(|| iv(0.0));
        out[i] = ai - bi;
    }
    out
}

fn poly_mul(a: &[Interval], b: &[Interval]) -> Vec<Interval> {
    if a.is_empty() || b.is_empty() {
        return vec![iv(0.0)];
    }
    let mut out = vec![iv(0.0); a.len() + b.len() - 1];
    for (i, ai) in a.iter().enumerate() {
        for (j, bj) in b.iter().enumerate() {
            out[i + j] = out[i + j] + (*ai * *bj);
        }
    }
    out
}

fn edge_power_poly(variable_is_x: bool, fixed: f64, roots: &[CInterval]) -> Vec<Interval> {
    let mut pr = vec![iv(1.0)];
    let mut pi = vec![iv(0.0)];
    for root in roots {
        let (fr, fi) = if variable_is_x {
            (vec![-root.re, iv(1.0)], vec![iv(fixed) - root.im])
        } else {
            (vec![iv(fixed) - root.re], vec![-root.im, iv(1.0)])
        };
        let next_pr = poly_sub(&poly_mul(&pr, &fr), &poly_mul(&pi, &fi));
        let next_pi = poly_add(&poly_mul(&pr, &fi), &poly_mul(&pi, &fr));
        pr = next_pr;
        pi = next_pi;
    }
    let mut f = poly_add(&poly_mul(&pr, &pr), &poly_mul(&pi, &pi));
    if f.is_empty() {
        f.push(iv(-1.0));
    } else {
        f[0] = f[0] - iv(1.0);
    }
    f
}

fn power_poly_to_local(poly: &[Interval], domain: [f64; 2]) -> Vec<Interval> {
    let n = poly.len().saturating_sub(1);
    let lo = domain[0];
    let width = domain[1] - domain[0];
    let mut out = vec![iv(0.0); n + 1];
    for i in 0..=n {
        for k in 0..=i {
            let scale = binom_f64(i, k) * lo.powi((i - k) as i32) * width.powi(k as i32);
            out[k] = out[k] + poly[i] * iv(scale);
        }
    }
    out
}

fn local_power_to_bernstein(local: &[Interval]) -> Vec<Interval> {
    let n = local.len().saturating_sub(1);
    let mut out = vec![iv(0.0); n + 1];
    for j in 0..=n {
        let mut bj = iv(0.0);
        for k in 0..=j {
            let scale = binom_f64(j, k) / binom_f64(n, k);
            bj = bj + local[k] * iv(scale);
        }
        out[j] = bj;
    }
    out
}

fn bernstein_range(coeffs: &[Interval]) -> Interval {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for c in coeffs {
        lo = lo.min(c.inf());
        hi = hi.max(c.sup());
    }
    if lo.is_finite() && hi.is_finite() {
        iwrap(lo, hi)
    } else {
        iv(0.0)
    }
}

fn tighten_interval(a: Interval, b: Interval) -> Interval {
    let lo = a.inf().max(b.inf());
    let hi = a.sup().min(b.sup());
    if lo <= hi {
        iwrap(lo, hi)
    } else {
        iwrap(a.inf().min(b.inf()), a.sup().max(b.sup()))
    }
}

fn bernstein_f_horizontal_edge(x: [f64; 2], y_const: f64, roots: &[CInterval]) -> Interval {
    let power = edge_power_poly(true, y_const, roots);
    let local = power_poly_to_local(&power, x);
    bernstein_range(&local_power_to_bernstein(&local))
}

fn bernstein_f_vertical_edge(x_const: f64, y: [f64; 2], roots: &[CInterval]) -> Interval {
    let power = edge_power_poly(false, x_const, roots);
    let local = power_poly_to_local(&power, y);
    bernstein_range(&local_power_to_bernstein(&local))
}

fn interval_midpoint(pair: [f64; 2]) -> f64 {
    0.5 * (pair[0] + pair[1])
}

fn strict_subset_interval(inner: Interval, outer: [f64; 2]) -> bool {
    inner.inf().is_finite()
        && inner.sup().is_finite()
        && inner.inf() > outer[0]
        && inner.sup() < outer[1]
}

fn krawczyk_contract_y(
    x: [f64; 2],
    y: [f64; 2],
    fy: Interval,
    roots: &[CInterval],
) -> Option<Interval> {
    if interval_abs_lower(fy) <= 0.0 {
        return None;
    }
    let y0 = interval_midpoint(y);
    let fy_mid = 0.5 * (fy.inf() + fy.sup());
    if !fy_mid.is_finite() || fy_mid == 0.0 {
        return None;
    }
    let f0 = f_interval_for_box(x, [y0, y0], roots);
    let y_delta = iwrap(y[0] - y0, y[1] - y0);
    let c = iv(1.0 / fy_mid);
    let k = iv(y0) - f0 * c + (iv(1.0) - fy * c) * y_delta;
    if strict_subset_interval(k, y) {
        Some(k)
    } else {
        None
    }
}

fn krawczyk_contract_x(
    x: [f64; 2],
    y: [f64; 2],
    fx: Interval,
    roots: &[CInterval],
) -> Option<Interval> {
    if interval_abs_lower(fx) <= 0.0 {
        return None;
    }
    let x0 = interval_midpoint(x);
    let fx_mid = 0.5 * (fx.inf() + fx.sup());
    if !fx_mid.is_finite() || fx_mid == 0.0 {
        return None;
    }
    let f0 = f_interval_for_box([x0, x0], y, roots);
    let x_delta = iwrap(x[0] - x0, x[1] - x0);
    let c = iv(1.0 / fx_mid);
    let k = iv(x0) - f0 * c + (iv(1.0) - fx * c) * x_delta;
    if strict_subset_interval(k, x) {
        Some(k)
    } else {
        None
    }
}

fn try_certify_tube_piece(
    source_tube_index: usize,
    depth: usize,
    x: [f64; 2],
    y: [f64; 2],
    roots: &[CInterval],
    local_bound_subdivision: usize,
) -> Result<Option<ResolvedTubePiece>, RemainingTubePiece> {
    let direct_f = f_interval_for_box(x, y, roots);
    if !f_may_contain_level(direct_f) {
        return Ok(None);
    }

    let (
        active,
        f,
        fx,
        fy,
        fx_abs_lower,
        fy_abs_lower,
        fx_abs_upper,
        fy_abs_upper,
        _combined_abs_upper,
        fx_sign,
        fy_sign,
    ) = local_box_stats(x, y, roots, local_bound_subdivision);
    if !active {
        return Ok(None);
    }

    let y_graph =
        if fy_abs_lower > 0.0 && matches!(fy_sign, SignState::Positive | SignState::Negative) {
            krawczyk_contract_y(x, y, fy, roots).map(|contracted| {
                let slope = fx_abs_upper / fy_abs_lower;
                let width = x[1] - x[0];
                (
                    "y_as_function_of_x".to_string(),
                    "Fy".to_string(),
                    fy_abs_lower,
                    fx_abs_upper,
                    slope,
                    width * (1.0 + slope * slope).sqrt(),
                    contracted,
                )
            })
        } else {
            None
        };

    let x_graph =
        if fx_abs_lower > 0.0 && matches!(fx_sign, SignState::Positive | SignState::Negative) {
            krawczyk_contract_x(x, y, fx, roots).map(|contracted| {
                let slope = fy_abs_upper / fx_abs_lower;
                let width = y[1] - y[0];
                (
                    "x_as_function_of_y".to_string(),
                    "Fx".to_string(),
                    fx_abs_lower,
                    fy_abs_upper,
                    slope,
                    width * (1.0 + slope * slope).sqrt(),
                    contracted,
                )
            })
        } else {
            None
        };

    let chosen = match (y_graph, x_graph) {
        (Some(a), Some(b)) => {
            if a.5 <= b.5 {
                Some(a)
            } else {
                Some(b)
            }
        }
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    };

    if let Some(chosen) = chosen {
        return Ok(Some(ResolvedTubePiece {
            source_tube_index,
            depth,
            x_interval: x,
            y_interval: y,
            contracted_interval: [chosen.6.inf(), chosen.6.sup()],
            certificate_method: "interval_krawczyk_with_signed_interval_derivative".to_string(),
            projection_axis: chosen.0,
            denominator_component: chosen.1,
            denominator_abs_lower: chosen.2,
            numerator_abs_upper: chosen.3,
            slope_abs_upper: chosen.4,
            length_upper: chosen.5,
        }));
    }

    Err(RemainingTubePiece {
        source_tube_index,
        depth,
        x_interval: x,
        y_interval: y,
        reason: "projection_certificate_failed".to_string(),
        f_interval: [f.inf(), f.sup()],
        fx_interval: [fx.inf(), fx.sup()],
        fy_interval: [fy.inf(), fy.sup()],
        fx_abs_lower,
        fy_abs_lower,
        fx_sign: sign_state_label(fx_sign).to_string(),
        fy_sign: sign_state_label(fy_sign).to_string(),
    })
}

fn resolve_tube_recursive(
    source_tube_index: usize,
    depth: usize,
    max_depth: usize,
    x: [f64; 2],
    y: [f64; 2],
    roots: &[CInterval],
    local_bound_subdivision: usize,
    stats: &mut TubeResolveStats,
) {
    stats.recursion_node_count += 1;
    stats.max_depth_reached = stats.max_depth_reached.max(depth);
    match try_certify_tube_piece(
        source_tube_index,
        depth,
        x,
        y,
        roots,
        local_bound_subdivision,
    ) {
        Ok(Some(piece)) => {
            stats.certified_piece_count += 1;
            stats.krawczyk_certified_count += 1;
            stats.resolved_pieces.push(piece);
        }
        Ok(None) => {
            stats.excluded_piece_count += 1;
        }
        Err(piece) if depth >= max_depth => {
            stats.remaining_pieces.push(RemainingTubePiece {
                reason: "depth_limit_projection_certificate_failed".to_string(),
                ..piece
            });
        }
        Err(_) => {
            let xw = x[1] - x[0];
            let yw = y[1] - y[0];
            if xw >= yw {
                let mid = 0.5 * (x[0] + x[1]);
                resolve_tube_recursive(
                    source_tube_index,
                    depth + 1,
                    max_depth,
                    [x[0], mid],
                    y,
                    roots,
                    local_bound_subdivision,
                    stats,
                );
                resolve_tube_recursive(
                    source_tube_index,
                    depth + 1,
                    max_depth,
                    [mid, x[1]],
                    y,
                    roots,
                    local_bound_subdivision,
                    stats,
                );
            } else {
                let mid = 0.5 * (y[0] + y[1]);
                resolve_tube_recursive(
                    source_tube_index,
                    depth + 1,
                    max_depth,
                    x,
                    [y[0], mid],
                    roots,
                    local_bound_subdivision,
                    stats,
                );
                resolve_tube_recursive(
                    source_tube_index,
                    depth + 1,
                    max_depth,
                    x,
                    [mid, y[1]],
                    roots,
                    local_bound_subdivision,
                    stats,
                );
            }
        }
    }
}

fn parse_tube_resolver_args() -> Result<(PathBuf, PathBuf, usize, usize, usize, usize), String> {
    let args: Vec<String> = env::args().collect();
    let mut source = PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01_RESULTS.json");
    let mut out_dir = PathBuf::from(
        "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-PARAM-SLICED-KRAWCZYK-PILOT-HARD-CELL-20260506-01",
    );
    let mut param_subdivision = 4usize;
    let mut piece_limit = 4096usize;
    let mut max_depth = 0usize;
    let mut local_bound_subdivision = 6usize;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--source" => {
                if i + 1 >= args.len() {
                    return Err("--source requires a path".to_string());
                }
                source = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--outdir" => {
                if i + 1 >= args.len() {
                    return Err("--outdir requires a path".to_string());
                }
                out_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--max-depth" => {
                if i + 1 >= args.len() {
                    return Err("--max-depth requires a positive integer".to_string());
                }
                max_depth = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --max-depth: {err}"))?;
                i += 2;
            }
            "--local-bound-subdivision" => {
                if i + 1 >= args.len() {
                    return Err("--local-bound-subdivision requires a positive integer".to_string());
                }
                local_bound_subdivision = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --local-bound-subdivision: {err}"))?;
                i += 2;
            }
            "--param-subdivision" => {
                if i + 1 >= args.len() {
                    return Err("--param-subdivision requires a positive integer".to_string());
                }
                param_subdivision = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --param-subdivision: {err}"))?;
                i += 2;
            }
            "--piece-limit" => {
                if i + 1 >= args.len() {
                    return Err("--piece-limit requires a nonnegative integer".to_string());
                }
                piece_limit = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --piece-limit: {err}"))?;
                i += 2;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
    }
    if param_subdivision == 0 {
        return Err("--param-subdivision must be positive".to_string());
    }
    Ok((
        source,
        out_dir,
        param_subdivision,
        piece_limit,
        max_depth,
        local_bound_subdivision,
    ))
}

fn write_tube_resolver_report(result: &serde_json::Value, path: &PathBuf) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Parameter-Sliced Krawczyk Resolver\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Source accepted length: `{}`\n\
- Max tile resolved length upper: `{}`\n\
- Max tile total validated length upper: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap: `{}`\n\
- Source unresolved pieces: `{}`\n\
- Processed pieces per tile: `{}`\n\
- Parameter subdivision: `{}`\n\
- Parameter tiles: `{}`\n\
- Total certified piece outcomes: `{}`\n\
- Total excluded piece outcomes: `{}`\n\
- Total remaining unresolved outcomes: `{}`\n\n\
## Interpretation\n\n\
This run targets the remaining pieces from the prior Krawczyk artifact after \
slicing the root-affine parameter cell. Parameter tiles are alternatives, so \
length is aggregated by maximum tile total, not by summing over tiles. This is \
not a proof of Erdős #114 and not a global n=14 proof.\n\n\
## Next Blocker\n\n\
{}\n",
        EXPERIMENT_ID,
        SOURCE_EXPERIMENT_ID,
        result["status"],
        result["source_accepted_length_upper"],
        result["max_tile_resolved_tube_length_upper"],
        result["max_tile_total_validated_length_upper"],
        result["exact_length_cap"],
        result["margin_to_cap"],
        result["source_unresolved_tube_count"],
        result["processed_piece_count"],
        result["param_subdivision"],
        result["param_tile_count"],
        result["total_certified_piece_count"],
        result["total_excluded_piece_count"],
        result["total_remaining_unresolved_piece_count"],
        result["next_blocker"]
            .as_str()
            .unwrap_or("continue diagnostics"),
    );
    fs::write(path, report)
        .map_err(|err| format!("failed to write report {}: {err}", path.display()))
}

fn tube_resolver_run() -> Result<(), String> {
    let started = Instant::now();
    let (source_path, out_dir, param_subdivision, piece_limit, max_depth, local_bound_subdivision) =
        parse_tube_resolver_args()?;
    fs::create_dir_all(&out_dir).map_err(|err| {
        format!(
            "failed to create output directory {}: {err}",
            out_dir.display()
        )
    })?;
    let result_path = out_dir.join(format!("{}_RESULTS.json", EXPERIMENT_ID));
    let report_path = out_dir.join(format!("{}_REPORT.md", EXPERIMENT_ID));
    let sha_path = out_dir.join(format!("{}_RESULTS.sha256", EXPERIMENT_ID));
    for path in [&result_path, &report_path, &sha_path] {
        if path.exists() {
            return Err(format!(
                "refusing to overwrite existing artifact: {}",
                path.display()
            ));
        }
    }

    let source_text = fs::read_to_string(&source_path)
        .map_err(|err| format!("failed to read source {}: {err}", source_path.display()))?;
    let source_json: serde_json::Value =
        serde_json::from_str(&source_text).map_err(|err| format!("bad source JSON: {err}"))?;
    let source_accepted_length = source_json["total_validated_length_upper"]
        .as_f64()
        .ok_or_else(|| "source missing total_validated_length_upper".to_string())?;
    let cap = source_json["exact_length_cap"]
        .as_f64()
        .ok_or_else(|| "source missing exact_length_cap".to_string())?;
    let unresolved = source_json["remaining_unresolved_pieces"]
        .as_array()
        .ok_or_else(|| "source missing remaining_unresolved_pieces".to_string())?;
    let processed_piece_count = if piece_limit == 0 {
        unresolved.len()
    } else {
        piece_limit.min(unresolved.len())
    };

    let mut per_tile_stats = Vec::<ParamTileStats>::new();
    for tile_i in 0..param_subdivision {
        for tile_j in 0..param_subdivision {
            let roots = root_intervals_for_param_tile(tile_i, tile_j, param_subdivision);
            let mut stats = TubeResolveStats {
                source_tube_count: processed_piece_count,
                ..TubeResolveStats::default()
            };
            for (idx, tube) in unresolved.iter().take(processed_piece_count).enumerate() {
                let x = interval_from_value(tube, "x_interval")?;
                let y = interval_from_value(tube, "y_interval")?;
                resolve_tube_recursive(
                    idx,
                    0,
                    max_depth,
                    x,
                    y,
                    &roots,
                    local_bound_subdivision,
                    &mut stats,
                );
            }
            let resolved_tube_length_upper: f64 = stats
                .resolved_pieces
                .iter()
                .map(|piece| piece.length_upper)
                .sum();
            let total_validated_length_upper = source_accepted_length + resolved_tube_length_upper;
            let remaining_reason_counts = stats.remaining_pieces.iter().fold(
                std::collections::BTreeMap::<String, usize>::new(),
                |mut acc, item| {
                    *acc.entry(item.reason.clone()).or_insert(0) += 1;
                    acc
                },
            );
            per_tile_stats.push(ParamTileStats {
                tile_i,
                tile_j,
                source_piece_count: processed_piece_count,
                certified_piece_count: stats.certified_piece_count,
                krawczyk_certified_count: stats.krawczyk_certified_count,
                excluded_piece_count: stats.excluded_piece_count,
                remaining_unresolved_piece_count: stats.remaining_pieces.len(),
                recursion_node_count: stats.recursion_node_count,
                max_depth_reached: stats.max_depth_reached,
                resolved_tube_length_upper,
                total_validated_length_upper,
                margin_to_cap: cap - total_validated_length_upper,
                remaining_reason_counts,
                resolved_piece_sample: stats.resolved_pieces.iter().take(20).cloned().collect(),
                remaining_piece_sample: stats.remaining_pieces.iter().take(20).cloned().collect(),
            });
        }
    }

    let max_tile_resolved_tube_length_upper = per_tile_stats
        .iter()
        .map(|tile| tile.resolved_tube_length_upper)
        .fold(0.0_f64, f64::max);
    let max_tile_total_validated_length_upper = per_tile_stats
        .iter()
        .map(|tile| tile.total_validated_length_upper)
        .fold(source_accepted_length, f64::max);
    let margin_to_cap = cap - max_tile_total_validated_length_upper;
    let total_certified_piece_count: usize = per_tile_stats
        .iter()
        .map(|tile| tile.certified_piece_count)
        .sum();
    let total_krawczyk_certified_count: usize = per_tile_stats
        .iter()
        .map(|tile| tile.krawczyk_certified_count)
        .sum();
    let total_excluded_piece_count: usize = per_tile_stats
        .iter()
        .map(|tile| tile.excluded_piece_count)
        .sum();
    let total_remaining_unresolved_piece_count: usize = per_tile_stats
        .iter()
        .map(|tile| tile.remaining_unresolved_piece_count)
        .sum();
    let worst_tile = per_tile_stats
        .iter()
        .max_by_key(|tile| tile.remaining_unresolved_piece_count)
        .cloned();

    let status = if piece_limit > 0 {
        if total_krawczyk_certified_count > 0 {
            "PARAM_SLICED_KRAWCZYK_PILOT_CERTIFIES"
        } else {
            "PARAM_SLICED_KRAWCZYK_PILOT_FAILS"
        }
    } else if total_remaining_unresolved_piece_count > 0 {
        "PARAM_SLICED_KRAWCZYK_FAIL_ISOLATION"
    } else if max_tile_total_validated_length_upper > cap {
        "PARAM_SLICED_KRAWCZYK_FAIL_BUDGET"
    } else {
        "PARAM_SLICED_KRAWCZYK_PASS_NOT_GLOBAL_PROOF"
    };
    let next_blocker = if piece_limit > 0 && total_krawczyk_certified_count == 0 {
        "Pilot did not certify any piece after root-parameter slicing. Retire this cheap route and move to full bivariate Bernstein/Krawczyk subdivision."
    } else if piece_limit > 0 {
        "Pilot certified at least one piece. Run the full local artifact with piece-limit 0 before drawing a certificate conclusion."
    } else if total_remaining_unresolved_piece_count > 0 {
        "Parameter slicing reduced but did not eliminate unresolved pieces. Increase parameter subdivision only if the pilot/full certification rate is material; otherwise use full bivariate Bernstein/Krawczyk."
    } else if max_tile_total_validated_length_upper > cap {
        "Every tile resolved, but the worst tile exceeds the length cap. Need sharper chart integration or a tighter source branch bound."
    } else {
        "Promote the hard-cell parameter-sliced certificate into a theorem-shaped packet, then test the next-hardest subcell."
    };
    let result = json!({
        "experiment_id": EXPERIMENT_ID,
        "timestamp_unix": unix_timestamp_string(),
        "source_experiment_id": SOURCE_EXPERIMENT_ID,
        "source_results_path": source_path,
        "status": status,
        "degree": DEGREE,
        "eps": EPS,
        "subcell": {"sub_i": SUB_I, "sub_j": SUB_J},
        "parameters": {
            "max_additional_depth": max_depth,
            "local_bound_subdivision": local_bound_subdivision,
            "resolver_rule": "split root-affine parameter cell; for each parameter tile, retry interval exclusion and one-dimensional Krawczyk on the same unresolved geometry",
        },
        "param_subdivision": param_subdivision,
        "param_tile_count": param_subdivision * param_subdivision,
        "piece_limit": piece_limit,
        "processed_piece_count": processed_piece_count,
        "source_accepted_length_upper": source_accepted_length,
        "max_tile_resolved_tube_length_upper": max_tile_resolved_tube_length_upper,
        "max_tile_total_validated_length_upper": max_tile_total_validated_length_upper,
        "exact_length_cap": cap,
        "margin_to_cap": margin_to_cap,
        "source_unresolved_tube_count": unresolved.len(),
        "total_certified_piece_count": total_certified_piece_count,
        "total_krawczyk_certified_count": total_krawczyk_certified_count,
        "total_excluded_piece_count": total_excluded_piece_count,
        "total_remaining_unresolved_piece_count": total_remaining_unresolved_piece_count,
        "per_tile_stats": per_tile_stats,
        "worst_tile": worst_tile,
        "claim_ceiling": "Targeted hard-cell parameter-sliced Krawczyk resolver only. Not a proof of Erdős #114 and not a global n=14 proof.",
        "next_blocker": next_blocker,
        "elapsed_secs": started.elapsed().as_secs_f64(),
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
    write_tube_resolver_report(&result, &report_path)?;
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
            "source_unresolved_tube_count": unresolved.len(),
            "processed_piece_count": processed_piece_count,
            "param_subdivision": param_subdivision,
            "param_tile_count": param_subdivision * param_subdivision,
            "total_certified_piece_count": total_certified_piece_count,
            "total_krawczyk_certified_count": total_krawczyk_certified_count,
            "total_excluded_piece_count": total_excluded_piece_count,
            "total_remaining_unresolved_piece_count": total_remaining_unresolved_piece_count,
            "source_accepted_length_upper": source_accepted_length,
            "max_tile_resolved_tube_length_upper": max_tile_resolved_tube_length_upper,
            "max_tile_total_validated_length_upper": max_tile_total_validated_length_upper,
            "exact_length_cap": cap,
            "margin_to_cap": margin_to_cap,
            "result": result_path,
            "report": report_path,
            "sha256": sha,
        }))
        .unwrap()
    );
    Ok(())
}

fn main() {
    if let Err(err) = tube_resolver_run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}
