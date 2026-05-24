//! EHP #114 n=14 critical-candidate affine-gradient pilot.
//!
//! This consumes the L21 global critical-point diagnostic and processes only
//! the 56 critical-candidate residual regions. It tracks the two root-affine
//! parameters linearly and places nonlinear products into interval remainders.
//! That makes this an affine-style dependency test, not a length promotion.
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

const EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-AFFINE-GRADIENT-HARD-CELL-20260506-01";
const PPRIME_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-HARD-CELL-20260506-01";
const PPRIME_P8_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01";
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
    affine_excluded_tile_count: usize,
    gradient_regular_tile_count: usize,
    still_critical_tile_count: usize,
    min_gradient_lower_bound: f64,
    first_failed_tile: Option<TileSummary>,
}

#[derive(Clone, Debug, Serialize)]
struct TileSummary {
    tile_i: usize,
    tile_j: usize,
    tile_status: String,
    f_interval: IntervalPair,
    fx_interval: IntervalPair,
    fy_interval: IntervalPair,
    pprime_re_interval: IntervalPair,
    pprime_im_interval: IntervalPair,
    gradient_lower_bound: f64,
}

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    out_dir: PathBuf,
    param_subdivision: usize,
    region_limit: usize,
    test_mode: String,
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

fn interval_avoids_zero(x: Interval) -> bool {
    x.inf() > 0.0 || x.sup() < 0.0
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

fn aff_neg(x: Affine) -> Affine {
    Affine {
        c: -x.c,
        a: -x.a,
        b: -x.b,
        rem: -x.rem,
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

fn caff_conj(x: CAffine) -> CAffine {
    CAffine {
        re: x.re,
        im: aff_neg(x.im),
    }
}

fn aff_abs_sq(z: CAffine) -> Affine {
    aff_add(aff_mul(z.re, z.re), aff_mul(z.im, z.im))
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

fn eval_p_p1(z: CAffine, roots: &[CAffine]) -> (CAffine, CAffine) {
    let mut p = caff_one();
    let mut p1 = caff_zero();
    for root in roots {
        let factor = caff_sub(z, *root);
        let next_p1 = caff_add(caff_mul(p1, factor), p);
        let next_p = caff_mul(p, factor);
        p = next_p;
        p1 = next_p1;
    }
    (p, p1)
}

fn caff_add(x: CAffine, y: CAffine) -> CAffine {
    CAffine {
        re: aff_add(x.re, y.re),
        im: aff_add(x.im, y.im),
    }
}

fn f_fx_fy_pprime_affine(
    z: CAffine,
    roots: &[CAffine],
) -> (Interval, Interval, Interval, Interval, Interval) {
    let (p, p1) = eval_p_p1(z, roots);
    let f = aff_interval(aff_sub(aff_abs_sq(p), aff_const(1.0)));
    let q = caff_mul(p1, caff_conj(p));
    let fx = aff_interval(q.re) * iv(2.0);
    let fy = aff_interval(q.im) * iv(-2.0);
    (f, fx, fy, aff_interval(p1.re), aff_interval(p1.im))
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
        .get("classification")
        .and_then(|v| v.as_str())
        .unwrap_or("")
        != "critical_candidate_region"
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
            .get("reason")
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
    test_mode: &str,
) -> TileSummary {
    let roots = root_affines_for_param_tile(cell, tile_i, tile_j, param_subdivision);
    let (f, fx, fy, pprime_re, pprime_im) =
        f_fx_fy_pprime_affine(caff_box(region.x, region.y), &roots);
    let fx_abs = interval_abs_lower(fx);
    let fy_abs = interval_abs_lower(fy);
    let grad = (fx_abs * fx_abs + fy_abs * fy_abs).sqrt();
    let pprime_avoids_zero = interval_avoids_zero(pprime_re) || interval_avoids_zero(pprime_im);
    let tile_status = if test_mode == "pprime" {
        if interval_avoids_zero(f) {
            "AFFINE_EXCLUDED_TILE"
        } else if pprime_avoids_zero {
            "PPRIME_EXCLUDED_CRITICAL_TILE"
        } else {
            "STILL_CRITICAL_TILE"
        }
    } else if interval_avoids_zero(f) {
        "AFFINE_EXCLUDED_TILE"
    } else if grad > 0.0 {
        "AFFINE_GRADIENT_REGULAR_TILE"
    } else {
        "STILL_CRITICAL_TILE"
    };
    TileSummary {
        tile_i,
        tile_j,
        tile_status: tile_status.to_string(),
        f_interval: ipair(f),
        fx_interval: ipair(fx),
        fy_interval: ipair(fy),
        pprime_re_interval: ipair(pprime_re),
        pprime_im_interval: ipair(pprime_im),
        gradient_lower_bound: grad,
    }
}

fn summarize_region(
    region: &CriticalRegion,
    cell: CellSpec,
    param_subdivision: usize,
    test_mode: &str,
) -> RegionSummary {
    let mut excluded = 0usize;
    let mut regular = 0usize;
    let mut critical = 0usize;
    let mut min_grad = f64::INFINITY;
    let mut first_failed = None;
    for tile_i in 0..param_subdivision {
        for tile_j in 0..param_subdivision {
            let tile = classify_tile(region, cell, tile_i, tile_j, param_subdivision, test_mode);
            min_grad = min_grad.min(tile.gradient_lower_bound);
            match tile.tile_status.as_str() {
                "AFFINE_EXCLUDED_TILE" => excluded += 1,
                "AFFINE_GRADIENT_REGULAR_TILE" | "PPRIME_EXCLUDED_CRITICAL_TILE" => regular += 1,
                _ => {
                    critical += 1;
                    if first_failed.is_none() {
                        first_failed = Some(tile.clone());
                    }
                }
            }
        }
    }
    if !min_grad.is_finite() {
        min_grad = 0.0;
    }
    let tile_count = param_subdivision * param_subdivision;
    let region_status = if critical == 0 && excluded == tile_count {
        "AFFINE_EXCLUDED_REGION"
    } else if critical == 0 {
        "AFFINE_PARTITION_CLOSED_REGION"
    } else {
        "STILL_CRITICAL_REGION"
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
        affine_excluded_tile_count: excluded,
        gradient_regular_tile_count: regular,
        still_critical_tile_count: critical,
        min_gradient_lower_bound: min_grad,
        first_failed_tile: first_failed,
    }
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01_RESULTS.json"),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-AFFINE-GRADIENT-HARD-CELL-20260506-01"),
        param_subdivision: 4,
        region_limit: 0,
        test_mode: "gradient".to_string(),
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
            "--test-mode" => {
                if i + 1 >= args.len() {
                    return Err("--test-mode requires gradient or pprime".to_string());
                }
                cfg.test_mode = args[i + 1].clone();
                if cfg.test_mode != "gradient" && cfg.test_mode != "pprime" {
                    return Err("--test-mode must be gradient or pprime".to_string());
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

fn experiment_id_for(test_mode: &str, param_subdivision: usize) -> &'static str {
    if test_mode == "pprime" && param_subdivision == 8 {
        PPRIME_P8_EXPERIMENT_ID
    } else if test_mode == "pprime" {
        PPRIME_EXPERIMENT_ID
    } else {
        EXPERIMENT_ID
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

fn first_failed_condition(rows: &[RegionSummary]) -> String {
    rows.iter()
        .find(|row| row.region_status == "STILL_CRITICAL_REGION")
        .map(|row| {
            format!(
                "{}:{} remains critical after affine-gradient parameter tiles",
                row.source_ownership_key, row.split_path
            )
        })
        .unwrap_or_else(|| "none".to_string())
}

fn status_for(
    test_mode: &str,
    still_critical_count: usize,
    closed_count: usize,
    total: f64,
    cap: f64,
) -> &'static str {
    if still_critical_count == 0 && total <= cap {
        if test_mode == "pprime" {
            "CRITICAL_CANDIDATE_PPRIME_PASS_NOT_GLOBAL_PROOF"
        } else {
            "CRITICAL_CANDIDATE_AFFINE_GRADIENT_PASS_NOT_GLOBAL_PROOF"
        }
    } else if closed_count > 0 {
        if test_mode == "pprime" {
            "CRITICAL_CANDIDATE_PPRIME_PARTIAL"
        } else {
            "CRITICAL_CANDIDATE_AFFINE_GRADIENT_PARTIAL"
        }
    } else {
        if test_mode == "pprime" {
            "CRITICAL_CANDIDATE_PPRIME_FAILS"
        } else {
            "CRITICAL_CANDIDATE_AFFINE_GRADIENT_FAILS"
        }
    }
}

fn write_report(result: &serde_json::Value, path: &PathBuf) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Critical-Candidate Affine Diagnostic\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Processed critical candidates: `{}`\n\
- Affine excluded regions: `{}`\n\
- Gradient affine regular regions: `{}`\n\
- Affine partition closed regions: `{}`\n\
- Still critical candidates: `{}`\n\
- Param subdivision: `{}`\n\
- Test mode: `{}`\n\
- Total validated length upper: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
This diagnostic tracks the two root-affine parameters linearly and places \
nonlinear products into interval remainders. It does not promote branch length. \
It asks whether any of the 56 critical candidates from L21 are actually \
regular or excluded once root-parameter dependency is partially preserved.\n\n\
## Claim Ceiling\n\n\
{}\n",
        EXPERIMENT_ID,
        result["source_experiment_id"]
            .as_str()
            .unwrap_or(SOURCE_EXPERIMENT_ID),
        result["status"],
        result["processed_critical_candidate_count"],
        result["affine_excluded_count"],
        result["gradient_affine_regular_count"],
        result["affine_partition_closed_count"],
        result["still_critical_candidate_count"],
        result["param_subdivision"],
        result["test_mode"],
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
        .unwrap_or_else(|| experiment_id_for(&cfg.test_mode, cfg.param_subdivision).to_string());
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
    let actual_source_experiment_id = source_json["experiment_id"]
        .as_str()
        .unwrap_or(SOURCE_EXPERIMENT_ID)
        .to_string();
    let source_accepted_length_upper = source_json["total_validated_length_upper"]
        .as_f64()
        .ok_or_else(|| "source missing total_validated_length_upper".to_string())?;
    let cap = source_json["exact_length_cap"]
        .as_f64()
        .ok_or_else(|| "source missing exact_length_cap".to_string())?;
    let classifications_json = source_json["classifications"]
        .as_array()
        .ok_or_else(|| "source missing classifications".to_string())?;
    let mut critical_regions = Vec::<CriticalRegion>::new();
    for (idx, value) in classifications_json.iter().enumerate() {
        if let Some(region) = parse_critical_region(idx, value)? {
            critical_regions.push(region);
        }
    }
    if cfg.region_limit > 0 {
        critical_regions.truncate(cfg.region_limit);
    }
    let summaries: Vec<RegionSummary> = critical_regions
        .iter()
        .map(|region| summarize_region(region, cfg.cell, cfg.param_subdivision, &cfg.test_mode))
        .collect();
    let affine_excluded_count = summaries
        .iter()
        .filter(|row| row.region_status == "AFFINE_EXCLUDED_REGION")
        .count();
    let affine_partition_closed_count = summaries
        .iter()
        .filter(|row| row.region_status == "AFFINE_PARTITION_CLOSED_REGION")
        .count();
    let gradient_affine_regular_count = summaries
        .iter()
        .filter(|row| {
            row.region_status == "AFFINE_PARTITION_CLOSED_REGION"
                && row.gradient_regular_tile_count > 0
        })
        .count();
    let still_critical_candidate_count = summaries
        .iter()
        .filter(|row| row.region_status == "STILL_CRITICAL_REGION")
        .count();
    let closed_count = affine_excluded_count + affine_partition_closed_count;
    let total_validated_length_upper = source_accepted_length_upper;
    let status = status_for(
        &cfg.test_mode,
        still_critical_candidate_count,
        closed_count,
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
        "param_subdivision": cfg.param_subdivision,
        "test_mode": cfg.test_mode,
        "parameters": {
            "region_limit": cfg.region_limit,
            "source_region_field": "classifications",
            "source_filter": "classification == critical_candidate_region",
            "affine_model": "two root-affine parameters tracked linearly; nonlinear products placed into interval remainders",
            "test_mode": cfg.test_mode,
            "tile_policy": "parameter tiles are alternatives, not simultaneous branches",
            "non_promotion_rule": "no branch length is promoted by this critical-candidate diagnostic"
        },
        "processed_critical_candidate_count": summaries.len(),
        "affine_excluded_count": affine_excluded_count,
        "gradient_affine_regular_count": gradient_affine_regular_count,
        "affine_partition_closed_count": affine_partition_closed_count,
        "still_critical_candidate_count": still_critical_candidate_count,
        "ownership_duplicate_count": 0,
        "source_accepted_length_upper": source_accepted_length_upper,
        "total_validated_length_upper": total_validated_length_upper,
        "exact_length_cap": cap,
        "margin_to_cap": cap - total_validated_length_upper,
        "first_failed_condition": first_failed_condition(&summaries),
        "region_summaries": summaries,
        "claim_ceiling": "Local n=14 hard-cell critical-candidate affine diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.",
        "next_blocker": if still_critical_candidate_count == 0 {
            "All 56 critical candidates closed under the selected affine test. Next proof-facing step is to integrate this with the branch atlas without promoting global claims."
        } else if closed_count > 0 {
            "Some critical candidates closed under the selected affine test. Next proof-facing step is to increase theorem strength only on the remaining candidates."
        } else {
            "This affine pilot closed no critical candidates. Next proof-facing step is a stronger analytic critical-point theorem or a root-location invariant, not more first-order affine bookkeeping."
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
            "test_mode": cfg.test_mode,
            "processed_critical_candidate_count": result["processed_critical_candidate_count"],
            "affine_excluded_count": result["affine_excluded_count"],
            "gradient_affine_regular_count": result["gradient_affine_regular_count"],
            "affine_partition_closed_count": result["affine_partition_closed_count"],
            "still_critical_candidate_count": result["still_critical_candidate_count"],
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
