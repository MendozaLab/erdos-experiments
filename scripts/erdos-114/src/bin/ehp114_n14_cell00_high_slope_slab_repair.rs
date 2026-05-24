//! EHP #114 n=14 CELL-00-00 targeted high-slope slab repair diagnostic.
//!
//! L33E showed that the CELL-00-00 slab-source overage is concentrated in a few
//! high-slope, low-denominator accepted branches. This binary repairs only
//! those branches by re-isolating their vertical tubes over smaller x-slabs.
//! The result is diagnostic source sharpening, not a proof promotion.

use ehp_n3_poc::ehp114_n14_cell::{ensure_source_subcell_matches, CellSpec};
use inari::{interval, Interval};
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const DEFAULT_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-CELL-00-00-HIGH-SLOPE-SLAB-REPAIR-20260506-01";
const DEFAULT_SOURCE: &str = "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01_RESULTS.json";
const DEFAULT_OUTDIR: &str =
    "../../Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-00-HIGH-SLOPE-SLAB-REPAIR-20260506-01";
const DEGREE: usize = 14;
const EPS: f64 = 0.1;
const EXACT_LENGTH_CAP: f64 = 20.672796062619668;
const DEFAULT_TARGETS: &str = "3713:1,4063:0,2584:3";

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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SignState {
    Positive,
    Negative,
    Mixed,
    ContainsZero,
}

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    out_dir: PathBuf,
    experiment_id: String,
    cell: CellSpec,
    targets: BTreeSet<String>,
    x_split: usize,
    y_split: usize,
    local_bound_subdivision: usize,
}

#[derive(Clone, Debug)]
struct SourceBranch {
    ix: usize,
    group_index: usize,
    ownership_key: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    y_cells: [usize; 2],
    projection_axis: String,
    denominator_abs_lower: f64,
    numerator_abs_upper: f64,
    slope_abs_upper: f64,
    length_upper: f64,
}

#[derive(Clone, Debug, Serialize)]
struct RepairedSegment {
    source_ownership_key: String,
    split_path: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    y_bins: [usize; 2],
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
}

#[derive(Clone, Debug, Serialize)]
struct RepairUnresolved {
    source_ownership_key: String,
    split_path: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    y_bins: [usize; 2],
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
}

#[derive(Clone, Debug, Serialize)]
struct TargetRepair {
    source_branch: SourceBranchJson,
    x_split: usize,
    y_split: usize,
    excluded_subcell_count: usize,
    candidate_subcell_count: usize,
    repaired_segment_count: usize,
    unresolved_group_count: usize,
    original_length_upper: f64,
    repaired_length_upper: f64,
    length_delta: f64,
    length_ratio: f64,
    status: String,
    repaired_segments: Vec<RepairedSegment>,
    unresolved_groups: Vec<RepairUnresolved>,
}

#[derive(Clone, Debug, Serialize)]
struct SourceBranchJson {
    ix: usize,
    group_index: usize,
    ownership_key: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    y_cells: [usize; 2],
    projection_axis: String,
    denominator_abs_lower: f64,
    numerator_abs_upper: f64,
    slope_abs_upper: f64,
    length_upper: f64,
}

#[derive(Default)]
struct TargetRepairWork {
    excluded_subcell_count: usize,
    candidate_subcell_count: usize,
    repaired_segments: Vec<RepairedSegment>,
    unresolved_groups: Vec<RepairUnresolved>,
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
    let a_pair = cell.u0_interval();
    let b_pair = cell.u1_interval();
    let a = iwrap(a_pair[0], a_pair[1]);
    let b = iwrap(b_pair[0], b_pair[1]);
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
    let q = ci_mul(p1, ci_conj(p));
    (q.re * iv(2.0), q.im * iv(-2.0))
}

fn f_interval_for_box(x: [f64; 2], y: [f64; 2], roots: &[CInterval]) -> Interval {
    let (p, _, _) = eval_p_p1_p2(ci_box(x, y), roots);
    abs_sq_minus_one(p)
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
            let (p, p1, _) = eval_p_p1_p2(ci_box(xs, ys), roots);
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
        fx_sign,
        fy_sign,
    )
}

fn opposite_strict_signs(a: Interval, b: Interval) -> bool {
    (a.sup() < 0.0 && b.inf() > 0.0) || (a.inf() > 0.0 && b.sup() < 0.0)
}

fn split_pair(pair: [f64; 2], n: usize, idx: usize) -> [f64; 2] {
    let step = (pair[1] - pair[0]) / (n as f64);
    [
        pair[0] + (idx as f64) * step,
        pair[0] + ((idx + 1) as f64) * step,
    ]
}

fn certify_local_group(
    source_key: &str,
    x_idx: usize,
    group_index: usize,
    x: [f64; 2],
    y: [f64; 2],
    y_bins: [usize; 2],
    roots: &[CInterval],
    local_bound_subdivision: usize,
) -> Result<RepairedSegment, RepairUnresolved> {
    let bottom = f_interval_for_box(x, [y[0], y[0]], roots);
    let top = f_interval_for_box(x, [y[1], y[1]], roots);
    let left = f_interval_for_box([x[0], x[0]], y, roots);
    let right = f_interval_for_box([x[1], x[1]], y, roots);
    let (
        _active,
        _f,
        fx,
        fy,
        fx_abs_lower,
        fy_abs_lower,
        fx_abs_upper,
        fy_abs_upper,
        fx_sign,
        fy_sign,
    ) = local_box_stats(x, y, roots, local_bound_subdivision);

    let split_path = format!("x{x_idx}:g{group_index}");
    let unresolved = |reason: &str| RepairUnresolved {
        source_ownership_key: source_key.to_string(),
        split_path: split_path.clone(),
        x_interval: x,
        y_interval: y,
        y_bins,
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
    };

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
        (None, None) => return Err(unresolved("projection_certificate_failed_after_repair")),
    };

    Ok(RepairedSegment {
        source_ownership_key: source_key.to_string(),
        split_path,
        x_interval: x,
        y_interval: y,
        y_bins,
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
    })
}

fn repair_branch(branch: &SourceBranch, roots: &[CInterval], cfg: &Config) -> TargetRepair {
    let mut work = TargetRepairWork::default();
    for xi in 0..cfg.x_split {
        let x = split_pair(branch.x_interval, cfg.x_split, xi);
        let mut candidate = vec![false; cfg.y_split];
        for (yi, slot) in candidate.iter_mut().enumerate() {
            let y = split_pair(branch.y_interval, cfg.y_split, yi);
            let f = f_interval_for_box(x, y, roots);
            if f_may_contain_level(f) {
                *slot = true;
                work.candidate_subcell_count += 1;
            } else {
                work.excluded_subcell_count += 1;
            }
        }
        let mut yi = 0usize;
        let mut group_index = 0usize;
        while yi < cfg.y_split {
            if !candidate[yi] {
                yi += 1;
                continue;
            }
            let start = yi;
            while yi + 1 < cfg.y_split && candidate[yi + 1] {
                yi += 1;
            }
            let end = yi;
            let y = [
                branch.y_interval[0]
                    + (start as f64) * (branch.y_interval[1] - branch.y_interval[0])
                        / (cfg.y_split as f64),
                branch.y_interval[0]
                    + ((end + 1) as f64) * (branch.y_interval[1] - branch.y_interval[0])
                        / (cfg.y_split as f64),
            ];
            match certify_local_group(
                &branch.ownership_key,
                xi,
                group_index,
                x,
                y,
                [start, end],
                roots,
                cfg.local_bound_subdivision,
            ) {
                Ok(segment) => work.repaired_segments.push(segment),
                Err(unresolved) => work.unresolved_groups.push(unresolved),
            }
            group_index += 1;
            yi += 1;
        }
    }

    let repaired_length_upper: f64 = work
        .repaired_segments
        .iter()
        .map(|segment| segment.length_upper)
        .sum();
    let status = if work.unresolved_groups.is_empty() && !work.repaired_segments.is_empty() {
        "TARGET_REPAIR_CERTIFIED"
    } else if !work.repaired_segments.is_empty() {
        "TARGET_REPAIR_PARTIAL"
    } else {
        "TARGET_REPAIR_FAILS"
    };
    TargetRepair {
        source_branch: SourceBranchJson {
            ix: branch.ix,
            group_index: branch.group_index,
            ownership_key: branch.ownership_key.clone(),
            x_interval: branch.x_interval,
            y_interval: branch.y_interval,
            y_cells: branch.y_cells,
            projection_axis: branch.projection_axis.clone(),
            denominator_abs_lower: branch.denominator_abs_lower,
            numerator_abs_upper: branch.numerator_abs_upper,
            slope_abs_upper: branch.slope_abs_upper,
            length_upper: branch.length_upper,
        },
        x_split: cfg.x_split,
        y_split: cfg.y_split,
        excluded_subcell_count: work.excluded_subcell_count,
        candidate_subcell_count: work.candidate_subcell_count,
        repaired_segment_count: work.repaired_segments.len(),
        unresolved_group_count: work.unresolved_groups.len(),
        original_length_upper: branch.length_upper,
        repaired_length_upper,
        length_delta: repaired_length_upper - branch.length_upper,
        length_ratio: if branch.length_upper > 0.0 {
            repaired_length_upper / branch.length_upper
        } else {
            0.0
        },
        status: status.to_string(),
        repaired_segments: work.repaired_segments,
        unresolved_groups: work.unresolved_groups,
    }
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from(DEFAULT_SOURCE),
        out_dir: PathBuf::from(DEFAULT_OUTDIR),
        experiment_id: DEFAULT_EXPERIMENT_ID.to_string(),
        cell: CellSpec::new(0, 0)?,
        targets: parse_targets(DEFAULT_TARGETS),
        x_split: 32,
        y_split: 256,
        local_bound_subdivision: 6,
    };
    let mut sub_i = cfg.cell.sub_i;
    let mut sub_j = cfg.cell.sub_j;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--source" => {
                cfg.source = PathBuf::from(next_string(&args, i, "--source")?);
                i += 2;
            }
            "--outdir" => {
                cfg.out_dir = PathBuf::from(next_string(&args, i, "--outdir")?);
                i += 2;
            }
            "--experiment-id" => {
                cfg.experiment_id = next_string(&args, i, "--experiment-id")?;
                i += 2;
            }
            "--sub-i" => {
                sub_i = next_string(&args, i, "--sub-i")?
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --sub-i: {err}"))?;
                i += 2;
            }
            "--sub-j" => {
                sub_j = next_string(&args, i, "--sub-j")?
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --sub-j: {err}"))?;
                i += 2;
            }
            "--targets" => {
                cfg.targets = parse_targets(&next_string(&args, i, "--targets")?);
                i += 2;
            }
            "--x-split" => {
                cfg.x_split =
                    parse_positive_usize(&next_string(&args, i, "--x-split")?, "--x-split")?;
                i += 2;
            }
            "--y-split" => {
                cfg.y_split =
                    parse_positive_usize(&next_string(&args, i, "--y-split")?, "--y-split")?;
                i += 2;
            }
            "--local-bound-subdivision" => {
                cfg.local_bound_subdivision = parse_positive_usize(
                    &next_string(&args, i, "--local-bound-subdivision")?,
                    "--local-bound-subdivision",
                )?;
                i += 2;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
    }
    cfg.cell = CellSpec::new(sub_i, sub_j)?;
    if cfg.targets.is_empty() {
        return Err("--targets must name at least one ownership key".to_string());
    }
    Ok(cfg)
}

fn parse_positive_usize(text: &str, flag: &str) -> Result<usize, String> {
    let value = text
        .parse::<usize>()
        .map_err(|err| format!("failed to parse {flag}: {err}"))?;
    if value == 0 {
        return Err(format!("{flag} must be positive"));
    }
    Ok(value)
}

fn parse_targets(text: &str) -> BTreeSet<String> {
    text.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(ToString::to_string)
        .collect()
}

fn next_string(args: &[String], i: usize, flag: &str) -> Result<String, String> {
    if i + 1 >= args.len() {
        return Err(format!("{flag} requires a value"));
    }
    Ok(args[i + 1].clone())
}

fn unix_timestamp_string() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs()
        .to_string()
}

fn read_json(path: &Path) -> Result<Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    serde_json::from_str(&text)
        .map_err(|err| format!("failed to parse {} as JSON: {err}", path.display()))
}

fn sha256_file(path: &Path) -> Result<String, String> {
    let bytes = fs::read(path)
        .map_err(|err| format!("failed to read {} for sha256: {err}", path.display()))?;
    Ok(sha256::digest(bytes))
}

fn value_pair(value: &Value, key: &str) -> Result<[f64; 2], String> {
    let arr = value
        .get(key)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("missing pair field {key}"))?;
    if arr.len() != 2 {
        return Err(format!("{key} must have length 2"));
    }
    Ok([
        arr[0]
            .as_f64()
            .ok_or_else(|| format!("{key}[0] is not numeric"))?,
        arr[1]
            .as_f64()
            .ok_or_else(|| format!("{key}[1] is not numeric"))?,
    ])
}

fn value_usize(value: &Value, key: &str) -> usize {
    value.get(key).and_then(Value::as_u64).unwrap_or(0) as usize
}

fn value_f64(value: &Value, key: &str) -> Result<f64, String> {
    value
        .get(key)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("missing numeric field {key}"))
}

fn value_string(value: &Value, key: &str) -> String {
    value
        .get(key)
        .and_then(Value::as_str)
        .unwrap_or("unknown")
        .to_string()
}

fn y_cells(value: &Value) -> [usize; 2] {
    let Some(arr) = value.get("y_cells").and_then(Value::as_array) else {
        return [0, 0];
    };
    if arr.len() != 2 {
        return [0, 0];
    }
    [
        arr[0].as_u64().unwrap_or(0) as usize,
        arr[1].as_u64().unwrap_or(0) as usize,
    ]
}

fn source_branch(value: &Value) -> Result<SourceBranch, String> {
    Ok(SourceBranch {
        ix: value_usize(value, "ix"),
        group_index: value_usize(value, "group_index"),
        ownership_key: value_string(value, "ownership_key"),
        x_interval: value_pair(value, "x_interval")?,
        y_interval: value_pair(value, "y_interval")?,
        y_cells: y_cells(value),
        projection_axis: value_string(value, "projection_axis"),
        denominator_abs_lower: value_f64(value, "denominator_abs_lower")?,
        numerator_abs_upper: value_f64(value, "numerator_abs_upper")?,
        slope_abs_upper: value_f64(value, "slope_abs_upper")?,
        length_upper: value_f64(value, "length_upper")?,
    })
}

fn source_total(value: &Value) -> Result<f64, String> {
    value_f64(value, "total_validated_length_upper")
}

fn source_margin(value: &Value) -> Result<f64, String> {
    value_f64(value, "margin_to_cap")
}

fn source_unresolved_count(value: &Value) -> usize {
    value_usize(value, "unresolved_branch_count")
}

fn source_branch_count(value: &Value) -> usize {
    value_usize(value, "slab_branch_count")
}

fn find_target_branches(
    source: &Value,
    targets: &BTreeSet<String>,
) -> Result<Vec<SourceBranch>, String> {
    let branches = source
        .get("branches")
        .and_then(Value::as_array)
        .ok_or_else(|| "source missing branches array".to_string())?;
    let mut found = Vec::new();
    for item in branches {
        let branch = source_branch(item)?;
        if targets.contains(&branch.ownership_key) {
            found.push(branch);
        }
    }
    let found_keys: BTreeSet<String> = found.iter().map(|b| b.ownership_key.clone()).collect();
    let missing: Vec<String> = targets.difference(&found_keys).cloned().collect();
    if !missing.is_empty() {
        return Err(format!(
            "source missing target branches: {}",
            missing.join(",")
        ));
    }
    found.sort_by(|a, b| a.ownership_key.cmp(&b.ownership_key));
    Ok(found)
}

fn write_report(result: &Value, path: &Path) -> Result<(), String> {
    let report_cell = result["cell_tag"].as_str().unwrap_or("unknown");
    let report = format!(
        "# EHP114 n=14 {report_cell} High-Slope Slab Repair Diagnostic\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Source total upper: `{}`\n\
- Exact cap: `{}`\n\
- Source margin: `{}`\n\
- Target count: `{}`\n\
- Target original length sum: `{}`\n\
- Repaired target length upper: `{}`\n\
- Target length delta: `{}`\n\
- Adjusted total if replaced: `{}`\n\
- Adjusted margin to cap: `{}`\n\
- Remaining unresolved source branches: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
This diagnostic only repairs the selected high-slope slab branches exposed by \
L33E. It splits each target slab in x, re-isolates candidate y-tubes inside the \
original branch tube, and recomputes graph length on the repaired collars. A \
good result here says the source overage is repairable locally; it does not \
certify the whole cell while other unresolved branches remain.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        result["status"].as_str().unwrap_or("UNKNOWN"),
        result["source_total_validated_length_upper"],
        result["exact_length_cap"],
        result["source_margin_to_cap"],
        result["target_count"],
        result["target_original_length_sum"],
        result["repaired_target_length_upper"],
        result["target_length_delta"],
        result["adjusted_total_if_replaced"],
        result["adjusted_margin_to_cap"],
        result["remaining_source_unresolved_branch_count"],
        result["first_failed_condition"]
            .as_str()
            .unwrap_or("unknown"),
        result["claim_ceiling"]
            .as_str()
            .unwrap_or("local diagnostic only"),
    );
    fs::write(path, report).map_err(|err| format!("failed to write {}: {err}", path.display()))
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let cfg = parse_args()?;
    fs::create_dir_all(&cfg.out_dir)
        .map_err(|err| format!("failed to create {}: {err}", cfg.out_dir.display()))?;
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

    let source = read_json(&cfg.source)?;
    let source_contract = ensure_source_subcell_matches(&source, cfg.cell)?;
    let roots = root_intervals(cfg.cell);
    let targets = find_target_branches(&source, &cfg.targets)?;
    let repairs: Vec<TargetRepair> = targets
        .iter()
        .map(|branch| repair_branch(branch, &roots, &cfg))
        .collect();
    let source_total = source_total(&source)?;
    let source_margin = source_margin(&source)?;
    let target_original_length_sum: f64 = repairs
        .iter()
        .map(|repair| repair.original_length_upper)
        .sum();
    let repaired_target_length_upper: f64 = repairs
        .iter()
        .map(|repair| repair.repaired_length_upper)
        .sum();
    let target_length_delta = repaired_target_length_upper - target_original_length_sum;
    let adjusted_total_if_replaced = source_total + target_length_delta;
    let adjusted_margin_to_cap = EXACT_LENGTH_CAP - adjusted_total_if_replaced;
    let unresolved_repair_group_count: usize = repairs
        .iter()
        .map(|repair| repair.unresolved_group_count)
        .sum();
    let certified_target_count = repairs
        .iter()
        .filter(|repair| repair.status == "TARGET_REPAIR_CERTIFIED")
        .count();
    let improved = target_length_delta < 0.0;
    let adjusted_under_cap = adjusted_total_if_replaced <= EXACT_LENGTH_CAP;
    let status = if unresolved_repair_group_count == 0 && adjusted_under_cap {
        "HIGH_SLOPE_SLAB_REPAIR_ADJUSTED_UNDER_CAP_DIAGNOSTIC"
    } else if unresolved_repair_group_count == 0 && improved {
        "HIGH_SLOPE_SLAB_REPAIR_REDUCES_OVERAGE"
    } else if certified_target_count > 0 {
        "HIGH_SLOPE_SLAB_REPAIR_PARTIAL"
    } else {
        "HIGH_SLOPE_SLAB_REPAIR_FAILS"
    };
    let first_failed_condition = if unresolved_repair_group_count > 0 {
        format!("{unresolved_repair_group_count} repaired target groups remain unresolved")
    } else if !adjusted_under_cap {
        format!(
            "adjusted total {} remains above cap {} by {}",
            adjusted_total_if_replaced,
            EXACT_LENGTH_CAP,
            adjusted_total_if_replaced - EXACT_LENGTH_CAP
        )
    } else {
        "none".to_string()
    };
    let cell_tag = cfg.cell.tag();
    let next_dependency = if adjusted_under_cap && unresolved_repair_group_count == 0 {
        "Apply the same targeted repair rule to the remaining high-contribution slabs, then regenerate a source artifact with replacement accounting."
            .to_string()
    } else if certified_target_count > 0 {
        "Increase local repair coverage or target the next high-slope slabs; this route has positive local signal but is not closed."
            .to_string()
    } else {
        format!(
            "Retire this local x-split/y-tube repair and move to a rotated or branch-centered source theorem for {cell_tag}."
        )
    };
    let what_this_rules_out = format!(
        "This tests whether the largest {cell_tag} source overage terms are artifacts of using one coarse high-slope slab tube."
    );
    let what_this_does_not_rule_out = format!(
        "It does not certify {cell_tag}, does not account for all remaining unresolved branches, and does not upgrade any global n=14 or Erdos #114 claim."
    );
    let claim_ceiling = format!(
        "{cell_tag} high-slope slab repair diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate."
    );

    let result = json!({
        "experiment_id": cfg.experiment_id,
        "timestamp_unix": unix_timestamp_string(),
        "elapsed_secs": started.elapsed().as_secs_f64(),
        "status": status,
        "source_results_path": cfg.source,
        "source_subcell_contract": source_contract,
        "subcell": cfg.cell,
        "cell_tag": cell_tag,
        "parameters": {
            "targets": cfg.targets.iter().cloned().collect::<Vec<_>>(),
            "x_split": cfg.x_split,
            "y_split": cfg.y_split,
            "local_bound_subdivision": cfg.local_bound_subdivision,
            "repair_rule": "split each target x-slab, re-scan y within the original branch tube, group active y bins, certify each mini-collar by endpoint sign separation and stable derivative chart, then replace only the targeted source branch length for diagnostic accounting"
        },
        "source_total_validated_length_upper": source_total,
        "exact_length_cap": EXACT_LENGTH_CAP,
        "source_margin_to_cap": source_margin,
        "source_slab_branch_count": source_branch_count(&source),
        "remaining_source_unresolved_branch_count": source_unresolved_count(&source),
        "target_count": repairs.len(),
        "target_original_length_sum": target_original_length_sum,
        "repaired_target_length_upper": repaired_target_length_upper,
        "target_length_delta": target_length_delta,
        "adjusted_total_if_replaced": adjusted_total_if_replaced,
        "adjusted_margin_to_cap": adjusted_margin_to_cap,
        "target_certified_count": certified_target_count,
        "target_unresolved_group_count": unresolved_repair_group_count,
        "repairs": repairs,
        "first_failed_condition": first_failed_condition,
        "what_this_rules_out": what_this_rules_out,
        "what_this_does_not_rule_out": what_this_does_not_rule_out,
        "next_dependency": next_dependency,
        "claim_ceiling": claim_ceiling
    });

    fs::write(
        &result_path,
        serde_json::to_string_pretty(&result).unwrap() + "\n",
    )
    .map_err(|err| format!("failed to write {}: {err}", result_path.display()))?;
    write_report(&result, &report_path)?;
    let digest = sha256_file(&result_path)?;
    fs::write(
        &sha_path,
        format!(
            "{}  {}\n",
            digest,
            result_path
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("RESULTS.json")
        ),
    )
    .map_err(|err| format!("failed to write {}: {err}", sha_path.display()))?;

    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "experiment_id": result["experiment_id"],
            "status": result["status"],
            "target_count": result["target_count"],
            "target_original_length_sum": result["target_original_length_sum"],
            "repaired_target_length_upper": result["repaired_target_length_upper"],
            "target_length_delta": result["target_length_delta"],
            "adjusted_total_if_replaced": result["adjusted_total_if_replaced"],
            "exact_length_cap": result["exact_length_cap"],
            "adjusted_margin_to_cap": result["adjusted_margin_to_cap"],
            "target_unresolved_group_count": result["target_unresolved_group_count"],
            "first_failed_condition": result["first_failed_condition"],
            "sha256": digest
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
