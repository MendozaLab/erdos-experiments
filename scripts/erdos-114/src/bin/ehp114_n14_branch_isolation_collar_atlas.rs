//! EHP #114 n=14 branch isolation collar-atlas diagnostic.
//!
//! This binary consumes the z32 slab validated-length artifact and targets only
//! its unresolved branch tubes. It tries interval exclusion first, then a
//! monotone collar lemma: a derivative component bounded away from zero plus
//! opposite signed collar walls certifies one owned graph branch and gives a
//! direct length upper bound.
//!
//! This remains local n=14 hard-cell work. It is not a proof of Erdos #114,
//! not a global n=14 certificate, and not an exact lemniscate-length
//! certificate.

use ehp_n3_poc::ehp114_n14_cell::{ensure_source_subcell_matches, CellSpec};
use inari::{interval, Interval};
use serde::Serialize;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const DEFAULT_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-HARD-CELL-20260506-01";
const SOURCE_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03";
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
struct SourceBranch {
    source_index: usize,
    ownership_key: String,
    ix: Option<u64>,
    group_index: Option<u64>,
    x: [f64; 2],
    y: [f64; 2],
    y_cells: Option<[u64; 2]>,
    source_reason: String,
}

#[derive(Clone, Debug, Serialize)]
struct CertifiedBranch {
    source_index: usize,
    source_ownership_key: String,
    ownership_key: String,
    split_path: String,
    depth: usize,
    chart_axis: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    derivative_abs_lower: f64,
    slope_abs_upper: f64,
    length_upper: f64,
    f_interval: IntervalPair,
    fx_interval: IntervalPair,
    fy_interval: IntervalPair,
    first_wall_f_interval: IntervalPair,
    second_wall_f_interval: IntervalPair,
    wall_margin_lower: f64,
}

#[derive(Clone, Debug, Serialize)]
struct ExcludedRegion {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    depth: usize,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    f_interval: IntervalPair,
}

#[derive(Clone, Debug, Serialize)]
struct UnresolvedRegion {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    depth: usize,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    reason: String,
    f_interval: IntervalPair,
    fx_interval: IntervalPair,
    fy_interval: IntervalPair,
    first_wall_f_interval: Option<IntervalPair>,
    second_wall_f_interval: Option<IntervalPair>,
}

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    out_dir: PathBuf,
    experiment_id: String,
    cell: CellSpec,
    branch_limit: usize,
    max_depth: usize,
}

#[derive(Default)]
struct RunStats {
    certified: Vec<CertifiedBranch>,
    excluded: Vec<ExcludedRegion>,
    unresolved: Vec<UnresolvedRegion>,
    recursion_node_count: usize,
    max_depth_reached: usize,
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

fn sign_stable(x: Interval) -> bool {
    x.inf() > 0.0 || x.sup() < 0.0
}

fn opposite_signed(a: Interval, b: Interval) -> bool {
    (a.inf() > 0.0 && b.sup() < 0.0) || (a.sup() < 0.0 && b.inf() > 0.0)
}

fn wall_margin(a: Interval, b: Interval) -> f64 {
    interval_abs_lower(a).min(interval_abs_lower(b))
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
    let mut p1 = ci_zero();
    for root in roots {
        let factor = ci_sub(z, *root);
        let next_p1 = ci_add(ci_mul(p1, factor), p);
        let next_p = ci_mul(p, factor);
        p = next_p;
        p1 = next_p1;
    }
    (p, p1)
}

fn f_interval_at(z: CInterval, roots: &[CInterval]) -> Interval {
    let (p, _) = eval_p_p1(z, roots);
    abs_sq(p) - iv(1.0)
}

fn gradient_intervals(p: CInterval, p1: CInterval) -> (Interval, Interval) {
    let q = ci_mul(p1, ci_conj(p));
    (q.re * iv(2.0), q.im * iv(-2.0))
}

fn f_fx_fy_on_box(x: [f64; 2], y: [f64; 2], roots: &[CInterval]) -> (Interval, Interval, Interval) {
    let z = ci_box(x, y);
    let (p, p1) = eval_p_p1(z, roots);
    let f = abs_sq(p) - iv(1.0);
    let (fx, fy) = gradient_intervals(p, p1);
    (f, fx, fy)
}

fn midpoint(pair: [f64; 2]) -> f64 {
    0.5 * (pair[0] + pair[1])
}

fn width(pair: [f64; 2]) -> f64 {
    pair[1] - pair[0]
}

fn split_box(x: [f64; 2], y: [f64; 2]) -> Vec<(String, [f64; 2], [f64; 2])> {
    if width(x).abs() >= width(y).abs() {
        let mid = midpoint(x);
        vec![
            ("x0".to_string(), [x[0], mid], y),
            ("x1".to_string(), [mid, x[1]], y),
        ]
    } else {
        let mid = midpoint(y);
        vec![
            ("y0".to_string(), x, [y[0], mid]),
            ("y1".to_string(), x, [mid, y[1]]),
        ]
    }
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

fn optional_u64_pair(piece: &serde_json::Value, field: &str) -> Option<[u64; 2]> {
    let arr = piece.get(field)?.as_array()?;
    if arr.len() != 2 {
        return None;
    }
    Some([arr[0].as_u64()?, arr[1].as_u64()?])
}

fn parse_source_branch(idx: usize, item: &serde_json::Value) -> Result<SourceBranch, String> {
    Ok(SourceBranch {
        source_index: idx,
        ownership_key: item
            .get("ownership_key")
            .and_then(|v| v.as_str())
            .unwrap_or("missing-ownership-key")
            .to_string(),
        ix: item.get("ix").and_then(|v| v.as_u64()),
        group_index: item.get("group_index").and_then(|v| v.as_u64()),
        x: source_interval(item, "x_interval")?,
        y: source_interval(item, "y_interval")?,
        y_cells: optional_u64_pair(item, "y_cells"),
        source_reason: item
            .get("reason")
            .and_then(|v| v.as_str())
            .unwrap_or("unknown")
            .to_string(),
    })
}

fn push_unique(indices: &mut Vec<usize>, seen: &mut BTreeSet<usize>, idx: usize, limit: usize) {
    if indices.len() < limit && seen.insert(idx) {
        indices.push(idx);
    }
}

fn select_branch_indices(source: &[serde_json::Value], limit: usize) -> Vec<usize> {
    if limit == 0 || limit >= source.len() {
        return (0..source.len()).collect();
    }
    let mut out = Vec::with_capacity(limit);
    let mut seen = BTreeSet::new();
    if let Some(idx) = source.iter().position(|v| {
        v.get("ownership_key")
            .and_then(|s| s.as_str())
            .map(|s| s == "2286:0")
            .unwrap_or(false)
    }) {
        push_unique(&mut out, &mut seen, idx, limit);
    }
    if let Some(idx) = source.iter().position(|v| {
        source_interval(v, "y_interval")
            .map(|y| midpoint(y) < 0.0)
            .unwrap_or(false)
    }) {
        push_unique(&mut out, &mut seen, idx, limit);
    }
    if let Some(idx) = source.iter().position(|v| {
        source_interval(v, "y_interval")
            .map(|y| midpoint(y) > 0.0)
            .unwrap_or(false)
    }) {
        push_unique(&mut out, &mut seen, idx, limit);
    }
    if let Some(idx) = source.iter().position(|v| {
        v.get("fx_sign")
            .and_then(|s| s.as_str())
            .map(|s| s != "contains_zero")
            .unwrap_or(false)
            && v.get("fy_sign")
                .and_then(|s| s.as_str())
                .map(|s| s == "contains_zero")
                .unwrap_or(false)
    }) {
        push_unique(&mut out, &mut seen, idx, limit);
    }
    for k in 0..limit {
        let idx = k * source.len() / limit;
        push_unique(&mut out, &mut seen, idx, limit);
    }
    let mut idx = 0usize;
    while out.len() < limit && idx < source.len() {
        push_unique(&mut out, &mut seen, idx, limit);
        idx += 1;
    }
    out
}

fn certify_x_as_function_of_y(
    branch: &SourceBranch,
    split_path: &str,
    depth: usize,
    x: [f64; 2],
    y: [f64; 2],
    f: Interval,
    fx: Interval,
    fy: Interval,
    roots: &[CInterval],
) -> Option<CertifiedBranch> {
    if !sign_stable(fx) {
        return None;
    }
    let left = f_interval_at(ci_box([x[0], x[0]], y), roots);
    let right = f_interval_at(ci_box([x[1], x[1]], y), roots);
    if !opposite_signed(left, right) {
        return None;
    }
    let denom = interval_abs_lower(fx);
    if denom <= 0.0 || !denom.is_finite() {
        return None;
    }
    let slope = interval_abs_upper(fy) / denom;
    if !slope.is_finite() {
        return None;
    }
    let length = width(y).abs() * (1.0 + slope * slope).sqrt();
    let key = format!("{}:{}:x_as_function_of_y", branch.ownership_key, split_path);
    Some(CertifiedBranch {
        source_index: branch.source_index,
        source_ownership_key: branch.ownership_key.clone(),
        ownership_key: key,
        split_path: split_path.to_string(),
        depth,
        chart_axis: "x_as_function_of_y".to_string(),
        x_interval: x,
        y_interval: y,
        derivative_abs_lower: denom,
        slope_abs_upper: slope,
        length_upper: length,
        f_interval: ipair(f),
        fx_interval: ipair(fx),
        fy_interval: ipair(fy),
        first_wall_f_interval: ipair(left),
        second_wall_f_interval: ipair(right),
        wall_margin_lower: wall_margin(left, right),
    })
}

fn certify_y_as_function_of_x(
    branch: &SourceBranch,
    split_path: &str,
    depth: usize,
    x: [f64; 2],
    y: [f64; 2],
    f: Interval,
    fx: Interval,
    fy: Interval,
    roots: &[CInterval],
) -> Option<CertifiedBranch> {
    if !sign_stable(fy) {
        return None;
    }
    let bottom = f_interval_at(ci_box(x, [y[0], y[0]]), roots);
    let top = f_interval_at(ci_box(x, [y[1], y[1]]), roots);
    if !opposite_signed(bottom, top) {
        return None;
    }
    let denom = interval_abs_lower(fy);
    if denom <= 0.0 || !denom.is_finite() {
        return None;
    }
    let slope = interval_abs_upper(fx) / denom;
    if !slope.is_finite() {
        return None;
    }
    let length = width(x).abs() * (1.0 + slope * slope).sqrt();
    let key = format!("{}:{}:y_as_function_of_x", branch.ownership_key, split_path);
    Some(CertifiedBranch {
        source_index: branch.source_index,
        source_ownership_key: branch.ownership_key.clone(),
        ownership_key: key,
        split_path: split_path.to_string(),
        depth,
        chart_axis: "y_as_function_of_x".to_string(),
        x_interval: x,
        y_interval: y,
        derivative_abs_lower: denom,
        slope_abs_upper: slope,
        length_upper: length,
        f_interval: ipair(f),
        fx_interval: ipair(fx),
        fy_interval: ipair(fy),
        first_wall_f_interval: ipair(bottom),
        second_wall_f_interval: ipair(top),
        wall_margin_lower: wall_margin(bottom, top),
    })
}

fn unresolved_reason(
    fx: Interval,
    fy: Interval,
    x_cert_possible: bool,
    y_cert_possible: bool,
) -> String {
    if !sign_stable(fx) && !sign_stable(fy) {
        "critical_collar_derivatives_not_sign_stable".to_string()
    } else if x_cert_possible || y_cert_possible {
        "collar_wall_sign_separation_failed".to_string()
    } else {
        "collar_certificate_failed".to_string()
    }
}

fn resolve_branch_recursive(
    branch: &SourceBranch,
    split_path: String,
    depth: usize,
    max_depth: usize,
    x: [f64; 2],
    y: [f64; 2],
    roots: &[CInterval],
    stats: &mut RunStats,
) {
    stats.recursion_node_count += 1;
    stats.max_depth_reached = stats.max_depth_reached.max(depth);
    let (f, fx, fy) = f_fx_fy_on_box(x, y, roots);
    if interval_avoids_zero(f) {
        stats.excluded.push(ExcludedRegion {
            source_index: branch.source_index,
            source_ownership_key: branch.ownership_key.clone(),
            split_path,
            depth,
            x_interval: x,
            y_interval: y,
            f_interval: ipair(f),
        });
        return;
    }

    let x_candidate =
        certify_x_as_function_of_y(branch, &split_path, depth, x, y, f, fx, fy, roots);
    let y_candidate =
        certify_y_as_function_of_x(branch, &split_path, depth, x, y, f, fx, fy, roots);
    let chosen = match (x_candidate, y_candidate) {
        (Some(a), Some(b)) => {
            if a.length_upper <= b.length_upper {
                Some(a)
            } else {
                Some(b)
            }
        }
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    };
    if let Some(cert) = chosen {
        stats.certified.push(cert);
        return;
    }

    if depth < max_depth {
        for (tag, sx, sy) in split_box(x, y) {
            let child_path = if split_path.is_empty() {
                tag
            } else {
                format!("{split_path}/{tag}")
            };
            resolve_branch_recursive(
                branch,
                child_path,
                depth + 1,
                max_depth,
                sx,
                sy,
                roots,
                stats,
            );
        }
        return;
    }

    let x_left = if sign_stable(fx) {
        Some(f_interval_at(ci_box([x[0], x[0]], y), roots))
    } else {
        None
    };
    let y_bottom = if sign_stable(fy) {
        Some(f_interval_at(ci_box(x, [y[0], y[0]]), roots))
    } else {
        None
    };
    let x_cert_possible = x_left.is_some();
    let y_cert_possible = y_bottom.is_some();
    stats.unresolved.push(UnresolvedRegion {
        source_index: branch.source_index,
        source_ownership_key: branch.ownership_key.clone(),
        split_path,
        depth,
        x_interval: x,
        y_interval: y,
        reason: unresolved_reason(fx, fy, x_cert_possible, y_cert_possible),
        f_interval: ipair(f),
        fx_interval: ipair(fx),
        fy_interval: ipair(fy),
        first_wall_f_interval: x_left.map(ipair).or_else(|| y_bottom.map(ipair)),
        second_wall_f_interval: if sign_stable(fx) {
            Some(ipair(f_interval_at(ci_box([x[1], x[1]], y), roots)))
        } else if sign_stable(fy) {
            Some(ipair(f_interval_at(ci_box(x, [y[1], y[1]]), roots)))
        } else {
            None
        },
    });
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03_RESULTS.json"),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-HARD-CELL-20260506-01"),
        experiment_id: DEFAULT_EXPERIMENT_ID.to_string(),
        cell: CellSpec::hard_cell(),
        branch_limit: 16,
        max_depth: 2,
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
            "--branch-limit" => {
                if i + 1 >= args.len() {
                    return Err("--branch-limit requires an integer".to_string());
                }
                cfg.branch_limit = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --branch-limit: {err}"))?;
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

fn first_failed_condition(unresolved: &[UnresolvedRegion]) -> String {
    unresolved
        .first()
        .map(|r| r.reason.clone())
        .unwrap_or_else(|| "none".to_string())
}

fn status_for(
    stats: &RunStats,
    source_unprocessed_count: usize,
    total_validated_length_upper: f64,
    cap: f64,
    ownership_duplicate_count: usize,
) -> &'static str {
    if ownership_duplicate_count > 0 {
        return "BRANCH_ISOLATION_FAIL_COLLAR";
    }
    if !stats.unresolved.is_empty() {
        let first_reason = first_failed_condition(&stats.unresolved);
        if first_reason == "critical_collar_derivatives_not_sign_stable" {
            return "BRANCH_ISOLATION_FAIL_CRITICAL_COLLAR";
        }
        if stats.max_depth_reached > 0 && source_unprocessed_count == 0 {
            return "BRANCH_ISOLATION_FAIL_DEPTH_LIMIT";
        }
        return "BRANCH_ISOLATION_FAIL_COLLAR";
    }
    if total_validated_length_upper > cap {
        return "BRANCH_ISOLATION_FAIL_BUDGET";
    }
    "BRANCH_ISOLATION_PASS_NOT_GLOBAL_PROOF"
}

fn next_blocker(status: &str, source_unprocessed_count: usize) -> &'static str {
    match status {
        "BRANCH_ISOLATION_PASS_NOT_GLOBAL_PROOF" if source_unprocessed_count > 0 => {
            "Pilot branch atlas closed the processed sample. Expand to --branch-limit 256 before a full unresolved-branch run."
        }
        "BRANCH_ISOLATION_PASS_NOT_GLOBAL_PROOF" => {
            "All processed unresolved branches were certified or excluded under the collar atlas. Package the local hard-cell certificate and test the next-hardest subcell."
        }
        "BRANCH_ISOLATION_FAIL_BUDGET" => {
            "Branch isolation closed topology but exceeded the length cap. Need sharper chart integration or tighter source branch bounds."
        }
        "BRANCH_ISOLATION_FAIL_CRITICAL_COLLAR" => {
            "Some regions lose fixed nonzero derivative in both coordinate charts. Next route needs a rotated chart or analytic critical-point exclusion."
        }
        "BRANCH_ISOLATION_FAIL_DEPTH_LIMIT" => {
            "Subdivision reached the requested depth without wall sign separation. Next route is higher-order Taylor/collar bounds, not raw depth escalation."
        }
        _ => {
            "Collar wall sign separation did not close every processed tube. Need sharper Taylor remainder bounds or analytic root-collar inequalities."
        }
    }
}

fn write_report(result: &serde_json::Value, path: &PathBuf) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Branch Isolation Collar Atlas\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Processed unresolved branches: `{}`\n\
- Certified branches: `{}`\n\
- Excluded regions: `{}`\n\
- Remaining unresolved branches: `{}`\n\
- Source accepted length: `{}`\n\
- Resolved unresolved-branch length upper: `{}`\n\
- Total validated length upper: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap: `{}`\n\
- Ownership duplicate count: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
This run targets only unresolved branch tubes from the z32 slab artifact. It \
uses exclusion first, then a monotone collar lemma: a fixed nonzero derivative \
component plus opposite signed collar walls certifies one owned graph branch. \
It does not rerun global slab length.\n\n\
## Next Blocker\n\n\
{}\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        SOURCE_EXPERIMENT_ID,
        result["status"],
        result["processed_unresolved_branch_count"],
        result["certified_branch_count"],
        result["excluded_region_count"],
        result["remaining_unresolved_branch_count"],
        result["source_accepted_length_upper"],
        result["resolved_unresolved_branch_length_upper"],
        result["total_validated_length_upper"],
        result["exact_length_cap"],
        result["margin_to_cap"],
        result["ownership_duplicate_count"],
        result["first_failed_condition"],
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
    let source_accepted_length = source_json["total_validated_length_upper"]
        .as_f64()
        .ok_or_else(|| "source missing total_validated_length_upper".to_string())?;
    let cap = source_json["length_budget"]["exact_length_cap"]
        .as_f64()
        .ok_or_else(|| "source missing length_budget.exact_length_cap".to_string())?;
    let unresolved = source_json["unresolved_branches"]
        .as_array()
        .ok_or_else(|| "source missing unresolved_branches".to_string())?;
    let selected_indices = select_branch_indices(unresolved, cfg.branch_limit);
    let roots = root_intervals(cfg.cell);
    let mut stats = RunStats::default();
    let mut processed_sources = Vec::<SourceBranch>::with_capacity(selected_indices.len());
    for idx in selected_indices {
        let branch = parse_source_branch(idx, &unresolved[idx])?;
        resolve_branch_recursive(
            &branch,
            "root".to_string(),
            0,
            cfg.max_depth,
            branch.x,
            branch.y,
            &roots,
            &mut stats,
        );
        processed_sources.push(branch);
    }

    let resolved_unresolved_branch_length_upper: f64 =
        stats.certified.iter().map(|p| p.length_upper).sum();
    let total_validated_length_upper =
        source_accepted_length + resolved_unresolved_branch_length_upper;
    let margin_to_cap = cap - total_validated_length_upper;
    let mut ownership_seen = BTreeSet::new();
    let ownership_duplicate_count = stats
        .certified
        .iter()
        .filter(|branch| !ownership_seen.insert(branch.ownership_key.clone()))
        .count();
    let source_unprocessed_count = unresolved.len().saturating_sub(processed_sources.len());
    let status = status_for(
        &stats,
        source_unprocessed_count,
        total_validated_length_upper,
        cap,
        ownership_duplicate_count,
    );
    let failed = first_failed_condition(&stats.unresolved);
    let mut unresolved_reason_counts = BTreeMap::<String, usize>::new();
    for item in &stats.unresolved {
        *unresolved_reason_counts
            .entry(item.reason.clone())
            .or_insert(0) += 1;
    }
    let processed_boxes: Vec<serde_json::Value> = processed_sources
        .iter()
        .map(|b| {
            json!({
                "source_index": b.source_index,
                "ownership_key": b.ownership_key,
                "ix": b.ix,
                "group_index": b.group_index,
                "source_reason": b.source_reason,
                "x_interval": b.x,
                "y_interval": b.y,
                "y_cells": b.y_cells,
            })
        })
        .collect();
    let worst_unresolved = stats.unresolved.first().cloned();

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
            "branch_limit": cfg.branch_limit,
            "max_depth": cfg.max_depth,
            "selection_rule": "worst ownership_key 2286:0 plus stratified positive/negative y and evenly spaced unresolved branches",
            "exclusion_rule": "interval F avoids zero on refined box",
            "collar_rule": "fixed nonzero derivative component plus opposite signed collar walls",
            "ownership_rule": "source ownership key plus split path plus chart axis; half-open seams are zero-length in this diagnostic",
        },
        "source_unresolved_branch_count": unresolved.len(),
        "processed_unresolved_branch_count": processed_sources.len(),
        "source_unprocessed_branch_count": source_unprocessed_count,
        "source_accepted_length_upper": source_accepted_length,
        "exact_length_cap": cap,
        "excluded_region_count": stats.excluded.len(),
        "certified_branch_count": stats.certified.len(),
        "remaining_unresolved_branch_count": stats.unresolved.len(),
        "resolved_unresolved_branch_length_upper": resolved_unresolved_branch_length_upper,
        "total_validated_length_upper": total_validated_length_upper,
        "margin_to_cap": margin_to_cap,
        "ownership_duplicate_count": ownership_duplicate_count,
        "first_failed_condition": failed,
        "unresolved_reason_counts": unresolved_reason_counts,
        "recursion_node_count": stats.recursion_node_count,
        "max_depth_reached": stats.max_depth_reached,
        "processed_source_branches": processed_boxes,
        "certified_branches": stats.certified,
        "excluded_regions": stats.excluded,
        "remaining_unresolved": stats.unresolved,
        "worst_unresolved": worst_unresolved,
        "claim_ceiling": "Local n=14 hard-cell branch-isolation diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.",
        "next_blocker": next_blocker(status, source_unprocessed_count),
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
            "processed_unresolved_branch_count": result["processed_unresolved_branch_count"],
            "certified_branch_count": result["certified_branch_count"],
            "excluded_region_count": result["excluded_region_count"],
            "remaining_unresolved_branch_count": result["remaining_unresolved_branch_count"],
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
