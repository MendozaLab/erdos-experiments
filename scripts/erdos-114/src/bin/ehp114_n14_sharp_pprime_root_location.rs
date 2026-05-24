//! EHP #114 n=14 sharp p-prime root-location diagnostic.
//!
//! This consumes the L28 p-prime root-location artifact and processes only the
//! 17 regions that remained near possible p-prime roots. It preserves the L28
//! root-free disk test, then applies a sharper Taylor/Rouche-style variation
//! bound: |p'(z)-p'(z0)| <= r |p''(z0)| + 0.5 r^2 sup |p'''|.
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
    "EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01";
const SOURCE_EXPERIMENT_ID: &str = "EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01";
const DEGREE: usize = 14;
const EPS: f64 = 0.1;

#[derive(Clone, Copy, Debug)]
struct Complex {
    re: f64,
    im: f64,
}

#[derive(Clone, Copy, Debug)]
struct Affine {
    c: f64,
    a: f64,
    b: f64,
    rem: Interval,
}

#[derive(Clone, Copy, Debug)]
struct CAffine {
    re: Affine,
    im: Affine,
}

#[derive(Clone, Copy, Debug, Serialize)]
struct IntervalPair {
    lo: f64,
    hi: f64,
}

#[derive(Clone, Debug)]
struct CriticalRegion {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    x: [f64; 2],
    y: [f64; 2],
    source_reason: String,
}

#[derive(Clone, Debug, Serialize)]
struct RegionSummary {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    source_reason: String,
    region_status: String,
    tile_count: usize,
    sharp_pprime_root_excluded_tile_count: usize,
    sharp_pprime_root_near_tile_count: usize,
    baseline_root_free_leaf_count: usize,
    taylor_root_free_leaf_count: usize,
    unresolved_leaf_count: usize,
    processed_leaf_count: usize,
    minimum_root_box_distance: f64,
    worst_tile: Option<TileSummary>,
}

#[derive(Clone, Debug, Serialize)]
struct TileSummary {
    tile_i: usize,
    tile_j: usize,
    tile_status: String,
    processed_leaf_count: usize,
    root_free_leaf_count: usize,
    unresolved_leaf_count: usize,
    baseline_root_free_leaf_count: usize,
    taylor_root_free_leaf_count: usize,
    worst_leaf: LeafSummary,
}

#[derive(Clone, Debug, Serialize)]
struct LeafSummary {
    depth: usize,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    leaf_status: String,
    pprime_re_interval: IntervalPair,
    pprime_im_interval: IntervalPair,
    pprime_lower_abs: f64,
    psecond_center_abs_upper: f64,
    psecond_box_abs_upper: f64,
    pthird_abs_upper: f64,
    spatial_radius: f64,
    baseline_margin: f64,
    taylor_margin: f64,
    root_box_distance_lower: f64,
}

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    out_dir: PathBuf,
    param_subdivision: usize,
    region_limit: usize,
    max_depth: usize,
    experiment_id: String,
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

fn interval_mid(pair: [f64; 2]) -> f64 {
    0.5 * (pair[0] + pair[1])
}

fn interval_radius(pair: [f64; 2]) -> f64 {
    0.5 * (pair[1] - pair[0])
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

fn complex_interval_abs_lower(re: Interval, im: Interval) -> f64 {
    let re_low = interval_abs_lower(re);
    let im_low = interval_abs_lower(im);
    (re_low * re_low + im_low * im_low).sqrt()
}

fn complex_interval_abs_upper(re: Interval, im: Interval) -> f64 {
    let re_up = interval_abs_upper(re);
    let im_up = interval_abs_upper(im);
    (re_up * re_up + im_up * im_up).sqrt()
}

fn linear_interval(a: f64, b: f64) -> Interval {
    let radius = a.abs() + b.abs();
    iwrap(-radius, radius)
}

fn aff_const(c: f64) -> Affine {
    Affine {
        c,
        a: 0.0,
        b: 0.0,
        rem: iv(0.0),
    }
}

fn aff_from_interval(pair: [f64; 2]) -> Affine {
    let c = interval_mid(pair);
    let r = interval_radius(pair);
    Affine {
        c,
        a: 0.0,
        b: 0.0,
        rem: iwrap(-r, r),
    }
}

fn aff_interval(x: Affine) -> Interval {
    iv(x.c) + linear_interval(x.a, x.b) + x.rem
}

fn aff_add(x: Affine, y: Affine) -> Affine {
    Affine {
        c: x.c + y.c,
        a: x.a + y.a,
        b: x.b + y.b,
        rem: x.rem + y.rem,
    }
}

fn aff_sub(x: Affine, y: Affine) -> Affine {
    Affine {
        c: x.c - y.c,
        a: x.a - y.a,
        b: x.b - y.b,
        rem: x.rem - y.rem,
    }
}

fn aff_scale(x: Affine, c: f64) -> Affine {
    Affine {
        c: c * x.c,
        a: c * x.a,
        b: c * x.b,
        rem: iv(c) * x.rem,
    }
}

fn aff_mul(x: Affine, y: Affine) -> Affine {
    let lx = linear_interval(x.a, x.b);
    let ly = linear_interval(y.a, y.b);
    let rem = lx * ly + iv(x.c) * y.rem + iv(y.c) * x.rem + lx * y.rem + ly * x.rem + x.rem * y.rem;
    Affine {
        c: x.c * y.c,
        a: x.c * y.a + y.c * x.a,
        b: x.c * y.b + y.c * x.b,
        rem,
    }
}

fn caff_one() -> CAffine {
    CAffine {
        re: aff_const(1.0),
        im: aff_const(0.0),
    }
}

fn caff_zero() -> CAffine {
    CAffine {
        re: aff_const(0.0),
        im: aff_const(0.0),
    }
}

fn caff_box(x: [f64; 2], y: [f64; 2]) -> CAffine {
    CAffine {
        re: aff_from_interval(x),
        im: aff_from_interval(y),
    }
}

fn caff_sub(x: CAffine, y: CAffine) -> CAffine {
    CAffine {
        re: aff_sub(x.re, y.re),
        im: aff_sub(x.im, y.im),
    }
}

fn caff_mul(x: CAffine, y: CAffine) -> CAffine {
    CAffine {
        re: aff_sub(aff_mul(x.re, y.re), aff_mul(x.im, y.im)),
        im: aff_add(aff_mul(x.re, y.im), aff_mul(x.im, y.re)),
    }
}

fn caff_scale(x: CAffine, c: f64) -> CAffine {
    CAffine {
        re: aff_scale(x.re, c),
        im: aff_scale(x.im, c),
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

fn root_affines_for_param_tile(
    cell: CellSpec,
    tile_i: usize,
    tile_j: usize,
    param_subdivision: usize,
) -> Vec<CAffine> {
    let u0_subcell = cell.u0_interval();
    let u1_subcell = cell.u1_interval();
    let u0_tile = split_interval(u0_subcell, param_subdivision)[tile_i];
    let u1_tile = split_interval(u1_subcell, param_subdivision)[tile_j];
    let a_mid = interval_mid(u0_tile);
    let a_rad = interval_radius(u0_tile);
    let b_mid = interval_mid(u1_tile);
    let b_rad = interval_radius(u1_tile);
    let scale = eps_shape_scale(EPS);
    let base = base_roots(EPS, DEGREE);
    let u0 = u0_direction();
    let u1 = u1_direction();
    (0..DEGREE)
        .map(|i| CAffine {
            re: Affine {
                c: base[i].re + scale * (a_mid * u0[i].re + b_mid * u1[i].re),
                a: scale * a_rad * u0[i].re,
                b: scale * b_rad * u1[i].re,
                rem: iv(0.0),
            },
            im: Affine {
                c: base[i].im + scale * (a_mid * u0[i].im + b_mid * u1[i].im),
                a: scale * a_rad * u0[i].im,
                b: scale * b_rad * u1[i].im,
                rem: iv(0.0),
            },
        })
        .collect()
}

fn eval_p_p1_p2_p3(z: CAffine, roots: &[CAffine]) -> (CAffine, CAffine, CAffine, CAffine) {
    let mut p = caff_one();
    let mut p1 = caff_zero();
    let mut p2 = caff_zero();
    let mut p3 = caff_zero();
    for root in roots {
        let factor = caff_sub(z, *root);
        let next_p3 = caff_add(caff_mul(p3, factor), caff_scale(p2, 3.0));
        let next_p2 = caff_add(caff_mul(p2, factor), caff_scale(p1, 2.0));
        let next_p1 = caff_add(caff_mul(p1, factor), p);
        let next_p = caff_mul(p, factor);
        p = next_p;
        p1 = next_p1;
        p2 = next_p2;
        p3 = next_p3;
    }
    (p, p1, p2, p3)
}

fn caff_add(x: CAffine, y: CAffine) -> CAffine {
    CAffine {
        re: aff_add(x.re, y.re),
        im: aff_add(x.im, y.im),
    }
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

fn parse_critical_region(
    idx: usize,
    value: &serde_json::Value,
) -> Result<Option<CriticalRegion>, String> {
    if value
        .get("region_status")
        .and_then(|v| v.as_str())
        .unwrap_or("")
        != "STILL_PPRIME_ROOT_NEAR"
    {
        return Ok(None);
    }
    Ok(Some(CriticalRegion {
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
        source_reason: value
            .get("source_reason")
            .and_then(|v| v.as_str())
            .unwrap_or("unknown")
            .to_string(),
    }))
}

fn classify_tile(
    region: &CriticalRegion,
    cell: CellSpec,
    tile_i: usize,
    tile_j: usize,
    param_subdivision: usize,
    max_depth: usize,
) -> TileSummary {
    let roots = root_affines_for_param_tile(cell, tile_i, tile_j, param_subdivision);
    let leaves = classify_box_recursive(region.x, region.y, 0, max_depth, &roots);
    let processed_leaf_count = leaves.len();
    let root_free_leaf_count = leaves
        .iter()
        .filter(|leaf| leaf.leaf_status != "SHARP_PPRIME_ROOT_NEAR_LEAF")
        .count();
    let unresolved_leaf_count = processed_leaf_count - root_free_leaf_count;
    let baseline_root_free_leaf_count = leaves
        .iter()
        .filter(|leaf| leaf.leaf_status == "BASELINE_ROOT_FREE_LEAF")
        .count();
    let taylor_root_free_leaf_count = leaves
        .iter()
        .filter(|leaf| leaf.leaf_status == "TAYLOR_ROOT_FREE_LEAF")
        .count();
    let worst_leaf = leaves
        .into_iter()
        .min_by(|a, b| {
            a.taylor_margin
                .partial_cmp(&b.taylor_margin)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .expect("at least one leaf");
    let tile_status = if unresolved_leaf_count == 0 {
        "SHARP_PPRIME_ROOT_FREE_TILE"
    } else {
        "SHARP_PPRIME_ROOT_NEAR_TILE"
    };
    TileSummary {
        tile_i,
        tile_j,
        tile_status: tile_status.to_string(),
        processed_leaf_count,
        root_free_leaf_count,
        unresolved_leaf_count,
        baseline_root_free_leaf_count,
        taylor_root_free_leaf_count,
        worst_leaf,
    }
}

fn classify_box_recursive(
    x: [f64; 2],
    y: [f64; 2],
    depth: usize,
    max_depth: usize,
    roots: &[CAffine],
) -> Vec<LeafSummary> {
    let leaf = classify_leaf(x, y, depth, roots);
    if leaf.leaf_status != "SHARP_PPRIME_ROOT_NEAR_LEAF" || depth >= max_depth {
        return vec![leaf];
    }
    let xm = interval_mid(x);
    let ym = interval_mid(y);
    let mut out = Vec::new();
    for (sx, sy) in [
        ([x[0], xm], [y[0], ym]),
        ([xm, x[1]], [y[0], ym]),
        ([x[0], xm], [ym, y[1]]),
        ([xm, x[1]], [ym, y[1]]),
    ] {
        out.extend(classify_box_recursive(sx, sy, depth + 1, max_depth, roots));
    }
    out
}

fn classify_leaf(x: [f64; 2], y: [f64; 2], depth: usize, roots: &[CAffine]) -> LeafSummary {
    let x_mid = interval_mid(x);
    let y_mid = interval_mid(y);
    let z_center = caff_box([x_mid, x_mid], [y_mid, y_mid]);
    let z_box = caff_box(x, y);
    let (_, p1_center, p2_center, _) = eval_p_p1_p2_p3(z_center, roots);
    let (_, _, p2_box, p3_box) = eval_p_p1_p2_p3(z_box, roots);
    let pprime_re = aff_interval(p1_center.re);
    let pprime_im = aff_interval(p1_center.im);
    let psecond_center_re = aff_interval(p2_center.re);
    let psecond_center_im = aff_interval(p2_center.im);
    let psecond_box_re = aff_interval(p2_box.re);
    let psecond_box_im = aff_interval(p2_box.im);
    let pthird_re = aff_interval(p3_box.re);
    let pthird_im = aff_interval(p3_box.im);
    let pprime_lower_abs = complex_interval_abs_lower(pprime_re, pprime_im);
    let psecond_center_abs_upper = complex_interval_abs_upper(psecond_center_re, psecond_center_im);
    let psecond_box_abs_upper = complex_interval_abs_upper(psecond_box_re, psecond_box_im);
    let pthird_abs_upper = complex_interval_abs_upper(pthird_re, pthird_im);
    let spatial_radius = (interval_radius(x).powi(2) + interval_radius(y).powi(2)).sqrt();
    let baseline_variation = spatial_radius * psecond_box_abs_upper;
    let taylor_variation =
        spatial_radius * psecond_center_abs_upper + 0.5 * spatial_radius.powi(2) * pthird_abs_upper;
    let baseline_margin = pprime_lower_abs - baseline_variation;
    let taylor_margin = pprime_lower_abs - taylor_variation;
    let root_box_distance_lower =
        if psecond_box_abs_upper > 0.0 && psecond_box_abs_upper.is_finite() {
            pprime_lower_abs / psecond_box_abs_upper - spatial_radius
        } else if pprime_lower_abs > 0.0 {
            f64::INFINITY
        } else {
            0.0
        };
    let leaf_status = if baseline_margin > 0.0 {
        "BASELINE_ROOT_FREE_LEAF"
    } else if taylor_margin > 0.0 {
        "TAYLOR_ROOT_FREE_LEAF"
    } else {
        "SHARP_PPRIME_ROOT_NEAR_LEAF"
    };
    LeafSummary {
        depth,
        x_interval: x,
        y_interval: y,
        leaf_status: leaf_status.to_string(),
        pprime_re_interval: ipair(pprime_re),
        pprime_im_interval: ipair(pprime_im),
        pprime_lower_abs,
        psecond_center_abs_upper,
        psecond_box_abs_upper,
        pthird_abs_upper,
        spatial_radius,
        baseline_margin,
        taylor_margin,
        root_box_distance_lower,
    }
}

fn summarize_region(
    region: &CriticalRegion,
    cell: CellSpec,
    param_subdivision: usize,
    max_depth: usize,
) -> RegionSummary {
    let mut root_excluded = 0usize;
    let mut root_near = 0usize;
    let mut baseline_root_free_leaf_count = 0usize;
    let mut taylor_root_free_leaf_count = 0usize;
    let mut unresolved_leaf_count = 0usize;
    let mut processed_leaf_count = 0usize;
    let mut minimum_root_box_distance = f64::INFINITY;
    let mut worst_tile = None;
    let mut worst_margin = f64::INFINITY;
    for tile_i in 0..param_subdivision {
        for tile_j in 0..param_subdivision {
            let tile = classify_tile(region, cell, tile_i, tile_j, param_subdivision, max_depth);
            minimum_root_box_distance =
                minimum_root_box_distance.min(tile.worst_leaf.root_box_distance_lower);
            baseline_root_free_leaf_count += tile.baseline_root_free_leaf_count;
            taylor_root_free_leaf_count += tile.taylor_root_free_leaf_count;
            unresolved_leaf_count += tile.unresolved_leaf_count;
            processed_leaf_count += tile.processed_leaf_count;
            if tile.worst_leaf.taylor_margin < worst_margin {
                worst_margin = tile.worst_leaf.taylor_margin;
                worst_tile = Some(tile.clone());
            }
            match tile.tile_status.as_str() {
                "SHARP_PPRIME_ROOT_FREE_TILE" => root_excluded += 1,
                _ => root_near += 1,
            }
        }
    }
    if !minimum_root_box_distance.is_finite() {
        minimum_root_box_distance = 0.0;
    }
    let tile_count = param_subdivision * param_subdivision;
    let region_status = if root_near == 0 {
        "SHARP_PPRIME_ROOT_EXCLUDED"
    } else {
        "STILL_SHARP_PPRIME_ROOT_NEAR"
    };
    RegionSummary {
        source_index: region.source_index,
        source_ownership_key: region.source_ownership_key.clone(),
        split_path: region.split_path.clone(),
        x_interval: region.x,
        y_interval: region.y,
        source_reason: region.source_reason.clone(),
        region_status: region_status.to_string(),
        tile_count,
        sharp_pprime_root_excluded_tile_count: root_excluded,
        sharp_pprime_root_near_tile_count: root_near,
        baseline_root_free_leaf_count,
        taylor_root_free_leaf_count,
        unresolved_leaf_count,
        processed_leaf_count,
        minimum_root_box_distance,
        worst_tile,
    }
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01_RESULTS.json"),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01"),
        param_subdivision: 8,
        region_limit: 0,
        max_depth: 4,
        experiment_id: DEFAULT_EXPERIMENT_ID.to_string(),
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
            "--param-subdivision" => {
                if i + 1 >= args.len() {
                    return Err("--param-subdivision requires an integer".to_string());
                }
                cfg.param_subdivision = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --param-subdivision: {err}"))?;
                if cfg.param_subdivision == 0 {
                    return Err("--param-subdivision must be positive".to_string());
                }
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

fn first_failed_condition(rows: &[RegionSummary]) -> String {
    rows.iter()
        .find(|row| row.region_status == "STILL_SHARP_PPRIME_ROOT_NEAR")
        .map(|row| {
            format!(
                "{}:{} remains near a possible p-prime root under the sharp Taylor/Rouche test",
                row.source_ownership_key, row.split_path
            )
        })
        .unwrap_or_else(|| "none".to_string())
}

fn status_for(
    still_unresolved_count: usize,
    excluded_count: usize,
    total: f64,
    cap: f64,
) -> &'static str {
    if total > cap {
        "SHARP_PPRIME_ROOT_LOCATION_FAIL_BUDGET"
    } else if still_unresolved_count == 0 {
        "SHARP_PPRIME_ROOT_LOCATION_PASS_NOT_GLOBAL_PROOF"
    } else if excluded_count > 0 {
        "SHARP_PPRIME_ROOT_LOCATION_PARTIAL"
    } else {
        "SHARP_PPRIME_ROOT_LOCATION_FAILS"
    }
}

fn write_report(result: &serde_json::Value, path: &PathBuf) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Sharp p-Prime Root-Location Diagnostic\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Processed remaining candidates: `{}`\n\
- p-prime root-excluded candidates: `{}`\n\
- p-prime root-near candidates: `{}`\n\
- Still unresolved candidates: `{}`\n\
- Minimum root-box distance lower: `{}`\n\
- Param subdivision: `{}`\n\
- Max local split depth: `{}`\n\
- Total validated length upper: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
This diagnostic uses the invariant that a critical point on `|p| = 1` must \
have `p'(z) = 0`. It preserves the L28 disk test and then tests the sharper \
Taylor/Rouche variation bound `r |p''(z0)| + 0.5 r^2 sup |p'''|` on targeted \
local splits inside the 17 unresolved boxes. It does not promote branch \
length.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        result["source_experiment_id"]
            .as_str()
            .unwrap_or(SOURCE_EXPERIMENT_ID),
        result["status"],
        result["processed_remaining_candidate_count"],
        result["sharp_pprime_root_excluded_count"],
        result["sharp_pprime_root_near_count"],
        result["still_unresolved_count"],
        result["minimum_root_box_distance"],
        result["param_subdivision"],
        result["max_depth"],
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
    let experiment_id = cfg.experiment_id.clone();
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
    let source_subcell_contract = ensure_source_subcell_matches(&source_json, cfg.cell)?;
    let source_experiment_id = source_json["experiment_id"]
        .as_str()
        .unwrap_or(SOURCE_EXPERIMENT_ID)
        .to_string();
    let source_accepted_length_upper = source_json["total_validated_length_upper"]
        .as_f64()
        .ok_or_else(|| "source missing total_validated_length_upper".to_string())?;
    let cap = source_json["exact_length_cap"]
        .as_f64()
        .ok_or_else(|| "source missing exact_length_cap".to_string())?;
    let region_summaries_json = source_json["region_summaries"]
        .as_array()
        .ok_or_else(|| "source missing region_summaries".to_string())?;
    let mut critical_regions = Vec::<CriticalRegion>::new();
    for (idx, value) in region_summaries_json.iter().enumerate() {
        if let Some(region) = parse_critical_region(idx, value)? {
            critical_regions.push(region);
        }
    }
    if cfg.region_limit > 0 {
        critical_regions.truncate(cfg.region_limit);
    }
    let summaries: Vec<RegionSummary> = critical_regions
        .iter()
        .map(|region| summarize_region(region, cfg.cell, cfg.param_subdivision, cfg.max_depth))
        .collect();
    let pprime_root_excluded_count = summaries
        .iter()
        .filter(|row| row.region_status == "SHARP_PPRIME_ROOT_EXCLUDED")
        .count();
    let pprime_root_near_count = summaries
        .iter()
        .filter(|row| row.region_status == "STILL_SHARP_PPRIME_ROOT_NEAR")
        .count();
    let still_unresolved_count = pprime_root_near_count;
    let minimum_root_box_distance = summaries
        .iter()
        .map(|row| row.minimum_root_box_distance)
        .fold(f64::INFINITY, f64::min);
    let minimum_root_box_distance = if minimum_root_box_distance.is_finite() {
        minimum_root_box_distance
    } else {
        0.0
    };
    let worst_region = summaries.iter().min_by(|a, b| {
        a.minimum_root_box_distance
            .partial_cmp(&b.minimum_root_box_distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let total_validated_length_upper = source_accepted_length_upper;
    let status = status_for(
        still_unresolved_count,
        pprime_root_excluded_count,
        total_validated_length_upper,
        cap,
    );
    let result = json!({
        "experiment_id": experiment_id,
        "timestamp_unix": unix_timestamp_string(),
        "source_experiment_id": source_experiment_id,
        "source_results_path": cfg.source,
        "source_subcell_contract": source_subcell_contract,
        "status": status,
        "degree": DEGREE,
        "eps": EPS,
        "subcell": cfg.cell,
        "cell_tag": cfg.cell.tag(),
        "param_subdivision": cfg.param_subdivision,
        "max_depth": cfg.max_depth,
        "parameters": {
            "region_limit": cfg.region_limit,
            "source_region_field": "region_summaries",
            "source_filter": "region_status == STILL_PPRIME_ROOT_NEAR",
            "affine_model": "two root-affine parameters tracked linearly; nonlinear products placed into interval remainders",
            "baseline_root_free_test": "lower_abs(p_prime(z0)) > spatial_radius * upper_abs(p_second(z_box)) for each root-parameter tile",
            "sharp_root_free_test": "lower_abs(p_prime(z0)) > spatial_radius * upper_abs(p_second(z0)) + 0.5 * spatial_radius^2 * upper_abs(p_third(z_box))",
            "local_split_policy": "recursively split only source regions still near p-prime roots, up to fixed max_depth",
            "tile_policy": "parameter tiles are alternatives, not simultaneous branches",
            "non_promotion_rule": "no branch length is promoted by this sharp p-prime root-location diagnostic"
        },
        "processed_remaining_candidate_count": summaries.len(),
        "sharp_pprime_root_excluded_count": pprime_root_excluded_count,
        "sharp_pprime_root_near_count": pprime_root_near_count,
        "baseline_root_free_leaf_count": summaries.iter().map(|row| row.baseline_root_free_leaf_count).sum::<usize>(),
        "taylor_root_free_leaf_count": summaries.iter().map(|row| row.taylor_root_free_leaf_count).sum::<usize>(),
        "unresolved_leaf_count": summaries.iter().map(|row| row.unresolved_leaf_count).sum::<usize>(),
        "processed_leaf_count": summaries.iter().map(|row| row.processed_leaf_count).sum::<usize>(),
        "still_unresolved_count": still_unresolved_count,
        "minimum_root_box_distance": minimum_root_box_distance,
        "ownership_duplicate_count": 0,
        "source_accepted_length_upper": source_accepted_length_upper,
        "total_validated_length_upper": total_validated_length_upper,
        "exact_length_cap": cap,
        "margin_to_cap": cap - total_validated_length_upper,
        "first_failed_condition": first_failed_condition(&summaries),
        "worst_region": worst_region,
        "region_summaries": summaries,
        "claim_ceiling": "Local n=14 hard-cell sharp p-prime root-location diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.",
        "next_blocker": if still_unresolved_count == 0 {
            "All remaining L28 critical candidates are separated from p-prime roots by the sharp Taylor/Rouche test. Next proof-facing step is to integrate this with the hard-cell branch atlas without promoting global claims."
        } else if pprime_root_excluded_count > 0 {
            "Some remaining critical candidates are separated from p-prime roots. Next proof-facing step is a still sharper root-location theorem only on the still-near candidates."
        } else {
            "The sharp Taylor/Rouche test closed no remaining candidates. Next proof-facing step is explicit derivative-root localization or exact algebraic resultant-style analysis, not more blind subdivision."
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
            "experiment_id": experiment_id,
            "status": status,
            "param_subdivision": cfg.param_subdivision,
            "max_depth": cfg.max_depth,
            "processed_remaining_candidate_count": result["processed_remaining_candidate_count"],
            "sharp_pprime_root_excluded_count": result["sharp_pprime_root_excluded_count"],
            "sharp_pprime_root_near_count": result["sharp_pprime_root_near_count"],
            "processed_leaf_count": result["processed_leaf_count"],
            "unresolved_leaf_count": result["unresolved_leaf_count"],
            "still_unresolved_count": result["still_unresolved_count"],
            "minimum_root_box_distance": result["minimum_root_box_distance"],
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
