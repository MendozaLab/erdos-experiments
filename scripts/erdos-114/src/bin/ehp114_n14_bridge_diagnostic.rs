//! EHP #114 n=14 worst-subcell exact-length bridge diagnostic.
//!
//! This binary probes the next proof bottleneck after the Rust root-affine
//! SUBDIV8 reproduction. It does not certify exact lemniscate length. It
//! computes interval regularity data for the worst accepted subcell `(6,4)`
//! and tests simple candidate contour-error budgets against the available
//! margin.

use inari::{interval, Interval};
use rayon::prelude::*;
use serde::Serialize;
use serde_json::json;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const EXPERIMENT_ID: &str = "EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01";
const SOURCE_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01";
const DEGREE: usize = 14;
const EPS: f64 = 0.1;
const EXTENT: f64 = 3.0;
const RES: usize = 220;
const SUBDIVISION: usize = 8;
const SUB_I: usize = 6;
const SUB_J: usize = 4;
const LSTAR_LOWER: f64 = 30.852910841548532;
const TARGET: f64 = 10.180114778928864;
const MARCHING_LENGTH_UPPER: f64 = 18.110795101366655;
const AVAILABLE_ERROR_BUDGET: f64 = 2.5620009612530126;
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

#[derive(Clone, Copy, Debug)]
struct MarchingCellProbe {
    active: bool,
    uncertain: bool,
    upper: f64,
}

#[derive(Clone, Debug, Serialize)]
struct ActiveCellDiagnostic {
    derivative_mode: DerivativeMode,
    regularity_strategy: RegularityStrategy,
    ix: usize,
    iy: usize,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    marching_cell_upper: f64,
    uncertain_corner_cell: bool,
    p_abs_upper: f64,
    p_prime_abs_lower: f64,
    p_prime_abs_lower_method: String,
    p_prime_abs_upper: f64,
    p_second_abs_upper: f64,
    local_bound_subdivision: usize,
    p_abs_upper_ambient: f64,
    p_abs_upper_local: f64,
    p_prime_abs_lower_ambient: f64,
    p_prime_abs_lower_taylor: f64,
    p_prime_abs_upper_ambient: f64,
    p_prime_abs_upper_local: f64,
    p_second_abs_upper_ambient: f64,
    p_second_abs_upper_local: f64,
    f_center_interval: [f64; 2],
    level_set_collar_radius: f64,
    level_set_collar_gradient_upper: f64,
    gradient_lower_on_level_candidate: f64,
    hessian_spectral_upper_candidate: f64,
    condition_ratio_candidate: Option<f64>,
    normal_drift_error_candidate: Option<f64>,
    relative_length_error_candidate: Option<f64>,
    regularity_resolved: bool,
}

#[derive(Clone, Copy, Debug)]
struct RunConfig {
    z_subdivision: usize,
    derivative_mode: DerivativeMode,
    regularity_strategy: RegularityStrategy,
    local_bound_subdivision: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
enum DerivativeMode {
    Recurrent,
    RootFactored,
    Compare,
}

impl DerivativeMode {
    fn parse(raw: &str) -> Result<Self, String> {
        match raw {
            "recurrent" => Ok(Self::Recurrent),
            "root-factored" | "root_factored" => Ok(Self::RootFactored),
            "compare" => Ok(Self::Compare),
            other => Err(format!(
                "unknown derivative mode `{other}`; expected recurrent, root-factored, or compare"
            )),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Recurrent => "recurrent",
            Self::RootFactored => "root_factored",
            Self::Compare => "compare",
        }
    }

    fn diagnostic_mode(self) -> Self {
        match self {
            Self::Compare => Self::RootFactored,
            mode => mode,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
enum RegularityStrategy {
    Ambient,
    TaylorLipschitz,
    LevelSetCollar,
    TaylorCollar,
}

impl RegularityStrategy {
    fn parse(raw: &str) -> Result<Self, String> {
        match raw {
            "ambient" => Ok(Self::Ambient),
            "taylor-lipschitz" | "taylor_lipschitz" => Ok(Self::TaylorLipschitz),
            "level-set-collar" | "level_set_collar" => Ok(Self::LevelSetCollar),
            "taylor-collar" | "taylor_collar" => Ok(Self::TaylorCollar),
            other => Err(format!(
                "unknown regularity strategy `{other}`; expected ambient, taylor-lipschitz, level-set-collar, or taylor-collar"
            )),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Ambient => "ambient",
            Self::TaylorLipschitz => "taylor_lipschitz",
            Self::LevelSetCollar => "level_set_collar",
            Self::TaylorCollar => "taylor_collar",
        }
    }

    fn uses_taylor(self) -> bool {
        matches!(self, Self::TaylorLipschitz | Self::TaylorCollar)
    }

    fn uses_collar(self) -> bool {
        matches!(self, Self::LevelSetCollar | Self::TaylorCollar)
    }
}

#[derive(Debug)]
enum BoxDiagnostic {
    Keep(ActiveCellDiagnostic),
    RejectInterval,
    RejectCollar,
}

#[derive(Default)]
struct DiagnosticRun {
    cells: Vec<ActiveCellDiagnostic>,
    interval_rejected_box_count: usize,
    collar_rejected_box_count: usize,
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

fn ci_fixed(z: Complex) -> CInterval {
    CInterval {
        re: iv(z.re),
        im: iv(z.im),
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

fn abs_lower(z: CInterval) -> f64 {
    abs_sq(z).inf().max(0.0).sqrt()
}

fn abs_upper(z: CInterval) -> f64 {
    abs_sq(z).sup().max(0.0).sqrt()
}

fn abs_sq_minus_one(z: CInterval) -> Interval {
    abs_sq(z) - iv(1.0)
}

fn sign(x: Interval) -> char {
    if x.inf() > 0.0 {
        '+'
    } else if x.sup() < 0.0 {
        '-'
    } else {
        '?'
    }
}

fn interp_interval(fa: Interval, fb: Interval) -> Interval {
    let denom = fa - fb;
    if denom.inf() <= 0.0 && denom.sup() >= 0.0 {
        return iwrap(0.0, 1.0);
    }
    let out = fa / denom;
    let lo = out.inf().max(0.0);
    let hi = out.sup().min(1.0);
    if lo > hi {
        iwrap(0.0, 1.0)
    } else {
        iwrap(lo, hi)
    }
}

fn seg_upper(ax: Interval, ay: Interval, bx: Interval, by: Interval) -> f64 {
    let dx = ax - bx;
    let dy = ay - by;
    (interval_square(dx) + interval_square(dy)).sqrt().sup()
}

fn cell_upper_from_fixed_case(
    case: u8,
    vals: (Interval, Interval, Interval, Interval),
    x0: f64,
    y0: f64,
    step: f64,
) -> (f64, bool) {
    let (fsw, fse, fne, fnw) = vals;
    let x0i = iv(x0);
    let y0i = iv(y0);
    let x1i = iv(x0 + step);
    let y1i = iv(y0 + step);
    let stepi = iv(step);
    let s = (x0i + interp_interval(fsw, fse) * stepi, y0i);
    let e = (x1i, y0i + interp_interval(fse, fne) * stepi);
    let n = (x0i + interp_interval(fnw, fne) * stepi, y1i);
    let w = (x0i, y0i + interp_interval(fsw, fnw) * stepi);
    let seg = |p: (Interval, Interval), q: (Interval, Interval)| seg_upper(p.0, p.1, q.0, q.1);

    match case {
        1 | 14 => (seg(s, w), true),
        2 | 13 => (seg(s, e), true),
        3 | 12 => (seg(w, e), true),
        4 | 11 => (seg(e, n), true),
        6 | 9 => (seg(s, n), true),
        7 | 8 => (seg(w, n), true),
        5 => {
            let avg = (fsw + fse + fne + fnw) * iv(0.25);
            if sign(avg) == '+' {
                (seg(s, w) + seg(e, n), true)
            } else if sign(avg) == '-' {
                (seg(s, e) + seg(w, n), true)
            } else {
                (2.0 * 2.0_f64.sqrt() * step, false)
            }
        }
        10 => {
            let avg = (fsw + fse + fne + fnw) * iv(0.25);
            if sign(avg) == '+' {
                (seg(s, e) + seg(w, n), true)
            } else if sign(avg) == '-' {
                (seg(s, w) + seg(e, n), true)
            } else {
                (2.0 * 2.0_f64.sqrt() * step, false)
            }
        }
        _ => (0.0, true),
    }
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

fn eval_p_prime_root_factored(z: CInterval, roots: &[CInterval]) -> CInterval {
    let factors: Vec<CInterval> = roots.iter().map(|root| ci_sub(z, *root)).collect();
    let mut sum = ci_zero();
    for omitted in 0..roots.len() {
        let mut term = ci_one();
        for (idx, factor) in factors.iter().enumerate() {
            if idx != omitted {
                term = ci_mul(term, *factor);
            }
        }
        sum = ci_add(sum, term);
    }
    sum
}

fn eval_p_p1_p2_with_mode(
    z: CInterval,
    roots: &[CInterval],
    derivative_mode: DerivativeMode,
) -> (CInterval, CInterval, CInterval) {
    let (p, p1_recurrent, p2) = eval_p_p1_p2(z, roots);
    let p1 = match derivative_mode.diagnostic_mode() {
        DerivativeMode::Recurrent => p1_recurrent,
        DerivativeMode::RootFactored => eval_p_prime_root_factored(z, roots),
        DerivativeMode::Compare => {
            unreachable!("compare resolves to root-factored for diagnostics")
        }
    };
    (p, p1, p2)
}

fn eval_f_at_point(z: Complex, roots: &[CInterval]) -> Interval {
    let (p, _, _) = eval_p_p1_p2(ci_fixed(z), roots);
    abs_sq_minus_one(p)
}

fn local_abs_upper_bounds(
    x: [f64; 2],
    y: [f64; 2],
    roots: &[CInterval],
    derivative_mode: DerivativeMode,
    subdivision: usize,
) -> (f64, f64, f64) {
    if subdivision <= 1 {
        let (p, p1, p2) = eval_p_p1_p2_with_mode(ci_box(x, y), roots, derivative_mode);
        return (abs_upper(p), abs_upper(p1), abs_upper(p2));
    }
    let x_step = (x[1] - x[0]) / (subdivision as f64);
    let y_step = (y[1] - y[0]) / (subdivision as f64);
    let mut p_abs_upper = 0.0_f64;
    let mut p_prime_abs_upper = 0.0_f64;
    let mut p_second_abs_upper = 0.0_f64;
    for sx in 0..subdivision {
        for sy in 0..subdivision {
            let x_sub = [
                x[0] + (sx as f64) * x_step,
                x[0] + ((sx + 1) as f64) * x_step,
            ];
            let y_sub = [
                y[0] + (sy as f64) * y_step,
                y[0] + ((sy + 1) as f64) * y_step,
            ];
            let (p, p1, p2) = eval_p_p1_p2_with_mode(ci_box(x_sub, y_sub), roots, derivative_mode);
            p_abs_upper = p_abs_upper.max(abs_upper(p));
            p_prime_abs_upper = p_prime_abs_upper.max(abs_upper(p1));
            p_second_abs_upper = p_second_abs_upper.max(abs_upper(p2));
        }
    }
    (p_abs_upper, p_prime_abs_upper, p_second_abs_upper)
}

fn probe_marching_cell(ix: usize, iy: usize, roots: &[CInterval], step: f64) -> MarchingCellProbe {
    let x0 = -EXTENT + (ix as f64) * step;
    let x1 = x0 + step;
    let y0 = -EXTENT + (iy as f64) * step;
    let y1 = y0 + step;
    let fsw = eval_f_at_point(Complex { re: x0, im: y0 }, roots);
    let fse = eval_f_at_point(Complex { re: x1, im: y0 }, roots);
    let fne = eval_f_at_point(Complex { re: x1, im: y1 }, roots);
    let fnw = eval_f_at_point(Complex { re: x0, im: y1 }, roots);
    let signs = [sign(fsw), sign(fse), sign(fne), sign(fnw)];
    if signs.iter().all(|s| *s == '+') || signs.iter().all(|s| *s == '-') {
        return MarchingCellProbe {
            active: false,
            uncertain: false,
            upper: 0.0,
        };
    }
    if signs.contains(&'?') {
        return MarchingCellProbe {
            active: true,
            uncertain: true,
            upper: 2.0 * 2.0_f64.sqrt() * step,
        };
    }
    let case = (if signs[0] == '+' { 1 } else { 0 })
        | ((if signs[1] == '+' { 1 } else { 0 }) << 1)
        | ((if signs[2] == '+' { 1 } else { 0 }) << 2)
        | ((if signs[3] == '+' { 1 } else { 0 }) << 3);
    let (upper, _) = cell_upper_from_fixed_case(case, (fsw, fse, fne, fnw), x0, y0, step);
    MarchingCellProbe {
        active: true,
        uncertain: false,
        upper,
    }
}

fn f_may_contain_level(f: Interval) -> bool {
    f.inf() <= 0.0 && f.sup() >= 0.0
}

fn diagnostic_for_box(
    ix: usize,
    iy: usize,
    x: [f64; 2],
    y: [f64; 2],
    roots: &[CInterval],
    derivative_mode: DerivativeMode,
    regularity_strategy: RegularityStrategy,
    local_bound_subdivision: usize,
    box_step: f64,
    marching_cell_upper: f64,
    uncertain_corner_cell: bool,
) -> BoxDiagnostic {
    let z_box = ci_box(x, y);
    let (p, p1, p2) = eval_p_p1_p2_with_mode(z_box, roots, derivative_mode);
    let f = abs_sq_minus_one(p);
    if !f_may_contain_level(f) {
        return BoxDiagnostic::RejectInterval;
    }
    let p_abs_upper_ambient = abs_upper(p);
    let p_prime_abs_lower_ambient = abs_lower(p1);
    let p_prime_abs_upper_ambient = abs_upper(p1);
    let p_second_abs_upper_ambient = abs_upper(p2);
    let (p_abs_upper_local, p_prime_abs_upper_local, p_second_abs_upper_local) =
        local_abs_upper_bounds(x, y, roots, derivative_mode, local_bound_subdivision);
    let p_abs_upper = p_abs_upper_local;
    let p_prime_abs_upper = p_prime_abs_upper_local;
    let p_second_abs_upper = p_second_abs_upper_local;
    let center = Complex {
        re: (x[0] + x[1]) * 0.5,
        im: (y[0] + y[1]) * 0.5,
    };
    let z_center = ci_fixed(center);
    let (_, p1_center, _) = eval_p_p1_p2_with_mode(z_center, roots, derivative_mode);
    let f_center = eval_f_at_point(center, roots);
    let radius = box_step / 2.0_f64.sqrt();
    let level_set_collar_gradient_upper = 2.0 * p_abs_upper * p_prime_abs_upper;
    let level_set_collar_radius = level_set_collar_gradient_upper * radius;
    if regularity_strategy.uses_collar()
        && (f_center.inf() > level_set_collar_radius || f_center.sup() < -level_set_collar_radius)
    {
        return BoxDiagnostic::RejectCollar;
    }
    let p_prime_abs_lower_taylor = (abs_lower(p1_center) - p_second_abs_upper * radius).max(0.0);
    let (p_prime_abs_lower, p_prime_abs_lower_method) = if regularity_strategy.uses_taylor() {
        (p_prime_abs_lower_taylor, "taylor_lipschitz")
    } else {
        (p_prime_abs_lower_ambient, "ambient_interval")
    };
    let gradient_lower = 2.0 * p_prime_abs_lower;

    // Candidate, conservative Hessian norm bound for F=|p|^2-1 over the box.
    // Each real Hessian entry is bounded by 2(|p'|^2 + |p||p''|);
    // the spectral norm of a 2x2 matrix is at most its Frobenius norm, here
    // bounded by twice that entry bound.
    let hessian_upper =
        4.0 * (p_prime_abs_upper * p_prime_abs_upper + p_abs_upper * p_second_abs_upper);
    let regularity_resolved = gradient_lower > 0.0 && gradient_lower.is_finite();
    let condition_ratio = if regularity_resolved {
        Some(hessian_upper / gradient_lower)
    } else {
        None
    };
    let normal_drift_error_candidate = condition_ratio.map(|r| r * box_step * box_step);
    let local_length_cap = 2.0 * 2.0_f64.sqrt() * box_step;
    let relative_length_error_candidate =
        condition_ratio.map(|r| local_length_cap.min(marching_cell_upper) * r * box_step);

    BoxDiagnostic::Keep(ActiveCellDiagnostic {
        derivative_mode: derivative_mode.diagnostic_mode(),
        regularity_strategy,
        ix,
        iy,
        x_interval: x,
        y_interval: y,
        marching_cell_upper,
        uncertain_corner_cell,
        p_abs_upper,
        p_prime_abs_lower,
        p_prime_abs_lower_method: p_prime_abs_lower_method.to_string(),
        p_prime_abs_upper,
        p_second_abs_upper,
        local_bound_subdivision,
        p_abs_upper_ambient,
        p_abs_upper_local,
        p_prime_abs_lower_ambient,
        p_prime_abs_lower_taylor,
        p_prime_abs_upper_ambient,
        p_prime_abs_upper_local,
        p_second_abs_upper_ambient,
        p_second_abs_upper_local,
        f_center_interval: [f_center.inf(), f_center.sup()],
        level_set_collar_radius,
        level_set_collar_gradient_upper,
        gradient_lower_on_level_candidate: gradient_lower,
        hessian_spectral_upper_candidate: hessian_upper,
        condition_ratio_candidate: condition_ratio,
        normal_drift_error_candidate,
        relative_length_error_candidate,
        regularity_resolved,
    })
}

fn diagnostics_for_cell(
    ix: usize,
    iy: usize,
    roots: &[CInterval],
    step: f64,
    z_subdivision: usize,
    derivative_mode: DerivativeMode,
    regularity_strategy: RegularityStrategy,
    local_bound_subdivision: usize,
) -> DiagnosticRun {
    let probe = probe_marching_cell(ix, iy, roots, step);
    if !probe.active {
        return DiagnosticRun::default();
    }
    let x0 = -EXTENT + (ix as f64) * step;
    let y0 = -EXTENT + (iy as f64) * step;
    let sub_step = step / (z_subdivision as f64);
    let mut out = DiagnosticRun::default();
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
            match diagnostic_for_box(
                ix,
                iy,
                x,
                y,
                roots,
                derivative_mode,
                regularity_strategy,
                local_bound_subdivision,
                sub_step,
                probe.upper,
                probe.uncertain,
            ) {
                BoxDiagnostic::Keep(row) => out.cells.push(row),
                BoxDiagnostic::RejectInterval => out.interval_rejected_box_count += 1,
                BoxDiagnostic::RejectCollar => out.collar_rejected_box_count += 1,
            }
        }
    }
    out
}

fn root_radius_upper(roots: &[CInterval]) -> f64 {
    roots
        .iter()
        .map(|r| abs_upper(*r))
        .fold(f64::NEG_INFINITY, f64::max)
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
    let mut out_dir = PathBuf::from(".");
    let mut z_subdivision = 1usize;
    let mut derivative_mode = DerivativeMode::Recurrent;
    let mut regularity_strategy = RegularityStrategy::Ambient;
    let mut local_bound_subdivision = 1usize;
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
            "--derivative-mode" => {
                if i + 1 >= args.len() {
                    return Err(
                        "--derivative-mode requires recurrent, root-factored, or compare"
                            .to_string(),
                    );
                }
                derivative_mode = DerivativeMode::parse(&args[i + 1])?;
                i += 2;
            }
            "--regularity-strategy" => {
                if i + 1 >= args.len() {
                    return Err("--regularity-strategy requires ambient, taylor-lipschitz, level-set-collar, or taylor-collar".to_string());
                }
                regularity_strategy = RegularityStrategy::parse(&args[i + 1])?;
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
            derivative_mode,
            regularity_strategy,
            local_bound_subdivision,
        },
    ))
}

fn write_report(result: &serde_json::Value, path: &PathBuf) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Worst-Subcell Bridge Diagnostic\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Worst subcell: `({}, {})`\n\
- Candidate level boxes: `{}`\n\
- z-subdivision per marching cell: `{}`\n\
- Derivative mode: `{}`\n\
- Regularity strategy: `{}`\n\
- Local bound subdivision: `{}`\n\
- Uncertain corner cells: `{}`\n\
- Interval-rejected boxes: `{}`\n\
- Collar-rejected boxes: `{}`\n\
- Regularity unresolved cells: `{}`\n\
- Recurrent unresolved cells: `{}`\n\
- Root-factored unresolved cells: `{}`\n\
- Regularity resolution delta: `{}`\n\
- Minimum candidate gradient lower: `{}`\n\
- Maximum candidate Hessian upper: `{}`\n\
- Sum candidate normal-drift error: `{}`\n\
- Sum candidate relative-length error: `{}`\n\
- Available exact-length bridge budget: `{}`\n\n\
- Normal error / budget: `{}`\n\
- Relative error / budget: `{}`\n\n\
## Claim Ceiling\n\n\
This is a bridge diagnostic only. It does not validate exact lemniscate length, \
does not prove the marching-squares-to-exact comparison theorem, and does not \
prove Erdős #114.\n\n\
## Next Blocker\n\n\
{}\n",
        EXPERIMENT_ID,
        SOURCE_EXPERIMENT_ID,
        result["status"].as_str().unwrap_or("UNKNOWN"),
        SUB_I,
        SUB_J,
        result["active_cell_count"],
        result["parameters"]["z_subdivision"],
        result["derivative_mode"].as_str().unwrap_or("unknown"),
        result["regularity_strategy"].as_str().unwrap_or("unknown"),
        result["parameters"]["local_bound_subdivision"],
        result["uncertain_corner_cell_count"],
        result["interval_rejected_box_count"],
        result["collar_rejected_box_count"],
        result["regularity_unresolved_cell_count"],
        result["unresolved_cells_recurrent"],
        result["unresolved_cells_root_factored"],
        result["regularity_resolution_delta"],
        result["min_gradient_lower_candidate"],
        result["max_hessian_upper_candidate"],
        result["sum_normal_drift_error_candidate"],
        result["sum_relative_length_error_candidate"],
        AVAILABLE_ERROR_BUDGET,
        result["error_to_budget_ratio_normal"],
        result["error_to_budget_ratio_relative"],
        result["next_blocker"]
            .as_str()
            .unwrap_or("continue diagnostic"),
    );
    fs::write(path, report)
        .map_err(|err| format!("failed to write report {}: {err}", path.display()))
}

fn run_diagnostics(
    roots: &[CInterval],
    step: f64,
    z_subdivision: usize,
    derivative_mode: DerivativeMode,
    regularity_strategy: RegularityStrategy,
    local_bound_subdivision: usize,
) -> DiagnosticRun {
    let runs: Vec<DiagnosticRun> = (0..RES)
        .into_par_iter()
        .map(|ix| {
            let mut acc = DiagnosticRun::default();
            for iy in 0..RES {
                let mut row_run = diagnostics_for_cell(
                    ix,
                    iy,
                    roots,
                    step,
                    z_subdivision,
                    derivative_mode,
                    regularity_strategy,
                    local_bound_subdivision,
                );
                acc.cells.append(&mut row_run.cells);
                acc.interval_rejected_box_count += row_run.interval_rejected_box_count;
                acc.collar_rejected_box_count += row_run.collar_rejected_box_count;
            }
            acc
        })
        .collect();

    let mut out = DiagnosticRun::default();
    for mut run in runs {
        out.cells.append(&mut run.cells);
        out.interval_rejected_box_count += run.interval_rejected_box_count;
        out.collar_rejected_box_count += run.collar_rejected_box_count;
    }
    out.cells.sort_by_key(|cell| (cell.ix, cell.iy));
    out
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let (out_dir, run_cfg) = parse_args()?;
    fs::create_dir_all(&out_dir).map_err(|err| {
        format!(
            "failed to create output directory {}: {err}",
            out_dir.display()
        )
    })?;
    let roots = root_intervals();
    let step = 2.0 * EXTENT / (RES as f64);

    let diagnostic_mode = run_cfg.derivative_mode.diagnostic_mode();
    let diagnostic_run = run_diagnostics(
        &roots,
        step,
        run_cfg.z_subdivision,
        diagnostic_mode,
        run_cfg.regularity_strategy,
        run_cfg.local_bound_subdivision,
    );

    let (unresolved_cells_recurrent, unresolved_cells_root_factored) =
        if run_cfg.derivative_mode == DerivativeMode::Compare {
            let recurrent_run = run_diagnostics(
                &roots,
                step,
                run_cfg.z_subdivision,
                DerivativeMode::Recurrent,
                run_cfg.regularity_strategy,
                run_cfg.local_bound_subdivision,
            );
            let recurrent_unresolved = recurrent_run
                .cells
                .iter()
                .filter(|c| !c.regularity_resolved)
                .count();
            let root_factored_unresolved = diagnostic_run
                .cells
                .iter()
                .filter(|c| !c.regularity_resolved)
                .count();
            (Some(recurrent_unresolved), Some(root_factored_unresolved))
        } else {
            (None, None)
        };

    let cells = diagnostic_run.cells;
    let active_cell_count = cells.len();
    let interval_rejected_box_count = diagnostic_run.interval_rejected_box_count;
    let collar_rejected_box_count = diagnostic_run.collar_rejected_box_count;
    let uncertain_corner_cell_count = cells.iter().filter(|c| c.uncertain_corner_cell).count();
    let regularity_unresolved_cell_count = cells.iter().filter(|c| !c.regularity_resolved).count();
    let min_gradient_lower = cells
        .iter()
        .map(|c| c.gradient_lower_on_level_candidate)
        .fold(f64::INFINITY, f64::min);
    let max_hessian_upper = cells
        .iter()
        .map(|c| c.hessian_spectral_upper_candidate)
        .fold(f64::NEG_INFINITY, f64::max);
    let max_condition_ratio = cells
        .iter()
        .filter_map(|c| c.condition_ratio_candidate)
        .fold(f64::NEG_INFINITY, f64::max);
    let sum_normal_drift_error: f64 = cells
        .iter()
        .filter_map(|c| c.normal_drift_error_candidate)
        .sum();
    let sum_relative_length_error: f64 = cells
        .iter()
        .filter_map(|c| c.relative_length_error_candidate)
        .sum();
    let error_to_budget_ratio_normal = sum_normal_drift_error / AVAILABLE_ERROR_BUDGET;
    let error_to_budget_ratio_relative = sum_relative_length_error / AVAILABLE_ERROR_BUDGET;
    let regularity_resolution_delta =
        match (unresolved_cells_recurrent, unresolved_cells_root_factored) {
            (Some(recurrent), Some(root_factored)) => {
                Some((recurrent as isize) - (root_factored as isize))
            }
            _ => None,
        };
    let required_ratio_for_normal_drift =
        AVAILABLE_ERROR_BUDGET / ((active_cell_count as f64) * step * step);
    let required_ratio_for_relative_length =
        AVAILABLE_ERROR_BUDGET / (MARCHING_LENGTH_UPPER * step);
    let candidate_normal_pass =
        regularity_unresolved_cell_count == 0 && sum_normal_drift_error <= AVAILABLE_ERROR_BUDGET;
    let candidate_relative_pass = regularity_unresolved_cell_count == 0
        && sum_relative_length_error <= AVAILABLE_ERROR_BUDGET;
    let status = if regularity_unresolved_cell_count > 0 {
        "BRIDGE_REGULARITY_INTERVAL_UNRESOLVED"
    } else if candidate_normal_pass || candidate_relative_pass {
        "BRIDGE_CANDIDATE_ERROR_UNDER_BUDGET_NOT_THEOREM"
    } else {
        "BRIDGE_CANDIDATE_ERROR_OVER_BUDGET"
    };
    let next_blocker = if regularity_unresolved_cell_count > 0 {
        "Interval boxes for p' still include zero on some active cells. Subdivide z-cells, use Bernstein/affine arithmetic, or restrict to a validated level-set collar before any length bridge can be proof-grade."
    } else if candidate_normal_pass || candidate_relative_pass {
        "Turn the candidate contour-error estimate into a theorem and preserve the same constants in a formal row record."
    } else {
        "The simple condition-ratio contour-error candidates are too coarse. Need sharper local charting, smaller z-cells, or a coarea/implicit-function bound with better constants."
    };

    let result = json!({
        "experiment_id": EXPERIMENT_ID,
        "timestamp_unix": unix_timestamp_string(),
        "source_experiment_id": SOURCE_EXPERIMENT_ID,
        "status": status,
        "derivative_mode": run_cfg.derivative_mode.as_str(),
        "regularity_strategy": run_cfg.regularity_strategy.as_str(),
        "parameters": {
            "degree": DEGREE,
            "eps": EPS,
            "extent": EXTENT,
            "res": RES,
            "grid_step": step,
            "z_subdivision": run_cfg.z_subdivision,
            "z_box_step": step / (run_cfg.z_subdivision as f64),
            "local_bound_subdivision": run_cfg.local_bound_subdivision,
            "subdivision": SUBDIVISION,
            "sub_i": SUB_I,
            "sub_j": SUB_J,
            "u0_interval": worst_subcell_intervals().0,
            "u1_interval": worst_subcell_intervals().1,
            "root_radius_upper": root_radius_upper(&roots),
        },
        "source_margins": {
            "lstar_lower": LSTAR_LOWER,
            "target": TARGET,
            "marching_length_upper": MARCHING_LENGTH_UPPER,
            "available_error_budget": AVAILABLE_ERROR_BUDGET,
        },
        "active_cell_count": active_cell_count,
        "interval_rejected_box_count": interval_rejected_box_count,
        "collar_rejected_box_count": collar_rejected_box_count,
        "uncertain_corner_cell_count": uncertain_corner_cell_count,
        "regularity_unresolved_cell_count": regularity_unresolved_cell_count,
        "min_gradient_lower_candidate": min_gradient_lower,
        "max_hessian_upper_candidate": max_hessian_upper,
        "max_condition_ratio_candidate": max_condition_ratio,
        "required_condition_ratio_for_normal_drift_budget": required_ratio_for_normal_drift,
        "required_condition_ratio_for_relative_length_budget": required_ratio_for_relative_length,
        "sum_normal_drift_error_candidate": sum_normal_drift_error,
        "sum_relative_length_error_candidate": sum_relative_length_error,
        "candidate_normal_drift_pass": candidate_normal_pass,
        "candidate_relative_length_pass": candidate_relative_pass,
        "error_to_budget_ratio_normal": error_to_budget_ratio_normal,
        "error_to_budget_ratio_relative": error_to_budget_ratio_relative,
        "regularity_resolution_delta": regularity_resolution_delta,
        "unresolved_cells_recurrent": unresolved_cells_recurrent,
        "unresolved_cells_root_factored": unresolved_cells_root_factored,
        "diagnostic_cells": cells,
        "claim_ceiling": "Worst-subcell bridge diagnostic only. This is not an exact lemniscate-length certificate, not a Lean theorem, and not a proof of Erdos #114.",
        "next_blocker": next_blocker,
        "elapsed_secs": started.elapsed().as_secs_f64(),
    });

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
            "active_cell_count": active_cell_count,
            "uncertain_corner_cell_count": uncertain_corner_cell_count,
            "regularity_unresolved_cell_count": regularity_unresolved_cell_count,
            "derivative_mode": run_cfg.derivative_mode.as_str(),
            "regularity_strategy": run_cfg.regularity_strategy.as_str(),
            "interval_rejected_box_count": interval_rejected_box_count,
            "collar_rejected_box_count": collar_rejected_box_count,
            "unresolved_cells_recurrent": unresolved_cells_recurrent,
            "unresolved_cells_root_factored": unresolved_cells_root_factored,
            "regularity_resolution_delta": regularity_resolution_delta,
            "min_gradient_lower_candidate": min_gradient_lower,
            "sum_normal_drift_error_candidate": sum_normal_drift_error,
            "sum_relative_length_error_candidate": sum_relative_length_error,
            "error_to_budget_ratio_normal": error_to_budget_ratio_normal,
            "error_to_budget_ratio_relative": error_to_budget_ratio_relative,
            "available_error_budget": AVAILABLE_ERROR_BUDGET,
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
        eprintln!("{err}");
        std::process::exit(1);
    }
}
