//! EHP #114 n=14 root-affine one-cell interval reproduction.
//!
//! This binary ports the Python SUBDIV8 continuous grid-oracle cell certificate
//! to Rust/inari. It is not exact lemniscate certification, not a Lean theorem,
//! and not a proof of Erdős #114.

use inari::{interval, Interval};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const DEFAULT_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01";
const PYTHON_PARENT_PACKET: &str =
    "EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01";
const DEGREE: usize = 14;
const DEFAULT_EPS: f64 = 0.1;
const DEFAULT_EXTENT: f64 = 3.0;
const DEFAULT_RES: usize = 220;
const DEFAULT_SUBDIVISION: usize = 8;
const DEFAULT_LSTAR_LOWER: f64 = 30.852910841548532;
const RADIAL_HALF: f64 = 12.0;
const PYTHON_MAX_LENGTH_UPPER: f64 = 18.110795101362747;
const PYTHON_MIN_MARGIN_LOWER: f64 = 2.5620009612569206;
const REPRO_TOLERANCE: f64 = 1.0e-9;

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
struct Pair(f64, f64);

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

#[derive(Clone, Debug, Deserialize)]
struct InputConfig {
    #[serde(default = "default_experiment_id")]
    experiment_id: String,
    #[serde(default = "default_degree")]
    degree: usize,
    #[serde(default = "default_eps")]
    eps: f64,
    #[serde(default = "default_extent")]
    extent: f64,
    #[serde(default = "default_res")]
    res: usize,
    #[serde(default = "default_subdivision")]
    subdivision: usize,
    #[serde(default = "default_lstar_lower")]
    lstar_lower: f64,
    #[serde(default)]
    target: Option<f64>,
    #[serde(default = "default_base_roots")]
    base_roots: Vec<Pair>,
    #[serde(default = "default_u0_direction")]
    u0_direction: Vec<Pair>,
    #[serde(default = "default_u1_direction")]
    u1_direction: Vec<Pair>,
    #[serde(default = "default_u0_interval")]
    u0_interval: [f64; 2],
    #[serde(default = "default_u1_interval")]
    u1_interval: [f64; 2],
}

#[derive(Clone, Debug)]
struct WorkConfig {
    experiment_id: String,
    degree: usize,
    eps: f64,
    extent: f64,
    res: usize,
    subdivision: usize,
    lstar_lower: f64,
    target: f64,
    base_roots: Vec<Complex>,
    u0_direction: Vec<Complex>,
    u1_direction: Vec<Complex>,
    u0_interval: [f64; 2],
    u1_interval: [f64; 2],
}

#[derive(Clone, Debug, Serialize)]
struct SubcellRow {
    sub_i: usize,
    sub_j: usize,
    u0_interval: [f64; 2],
    u1_interval: [f64; 2],
    length_upper: f64,
    deficit_lower: f64,
    margin_lower: f64,
    pass: bool,
    active_cells: usize,
    definite_case_cells: usize,
    uncertain_corner_cells: usize,
    ambiguous_avg_cells: usize,
}

#[derive(Clone, Copy, Debug)]
struct MarchingStats {
    length_upper: f64,
    active_cells: usize,
    definite_case_cells: usize,
    uncertain_corner_cells: usize,
    ambiguous_avg_cells: usize,
    max_cell_upper: f64,
}

fn default_experiment_id() -> String {
    DEFAULT_EXPERIMENT_ID.to_string()
}

fn default_degree() -> usize {
    DEGREE
}

fn default_eps() -> f64 {
    DEFAULT_EPS
}

fn default_extent() -> f64 {
    DEFAULT_EXTENT
}

fn default_res() -> usize {
    DEFAULT_RES
}

fn default_subdivision() -> usize {
    DEFAULT_SUBDIVISION
}

fn default_lstar_lower() -> f64 {
    DEFAULT_LSTAR_LOWER
}

fn default_u0_interval() -> [f64; 2] {
    [-0.0017499999999999998, 0.0]
}

fn default_u1_interval() -> [f64; 2] {
    [0.0, 0.0017499999999999998]
}

fn default_base_roots() -> Vec<Pair> {
    base_roots(DEFAULT_EPS, DEGREE)
        .into_iter()
        .map(|z| Pair(z.re, z.im))
        .collect()
}

fn default_u0_direction() -> Vec<Pair> {
    vec![
        Pair(2.1464806571511146e-13, 0.12536063016562018),
        Pair(-0.11029052423146467, 0.13183099723738242),
        Pair(-0.2026155037603651, 0.044088538363229984),
        Pair(-0.3356668203914393, -0.23859985068244477),
        Pair(0.33566682039232454, -0.23859985068243617),
        Pair(0.20261550376115514, 0.044088538363970045),
        Pair(0.11029052422993874, 0.13183099723676997),
        Pair(7.771839274878419e-14, 0.12536063016442625),
        Pair(-0.11029052423164845, 0.13183099723560737),
        Pair(-0.20261550375967857, 0.04408853836187157),
        Pair(-0.33566682038773565, -0.23859985068065273),
        Pair(0.3356668203871522, -0.23859985068104475),
        Pair(0.20261550375874843, 0.04408853836095835),
        Pair(0.1102905242327204, 0.13183099723674271),
    ]
}

fn default_u1_direction() -> Vec<Pair> {
    vec![
        Pair(-0.08603626747189994, 6.180051071480805e-13),
        Pair(-0.09672845897773669, -0.05128639870562298),
        Pair(0.12294713305603927, -0.2280270846819635),
        Pair(0.3689060705436963, 0.17637503539601224),
        Pair(-0.36890607054220126, 0.17637503539629037),
        Pair(-0.12294713305686052, -0.22802708467920318),
        Pair(0.09672845897754945, -0.05128639870782282),
        Pair(0.0860362674705226, 3.4447185685221385e-13),
        Pair(0.09672845897619298, 0.0512863987081996),
        Pair(-0.122947133057899, 0.22802708468114433),
        Pair(-0.36890607054741453, -0.17637503539846167),
        Pair(0.3689060705470212, -0.17637503539817698),
        Pair(0.12294713305902331, 0.22802708468019145),
        Pair(-0.09672845897603356, 0.05128639870845073),
    ]
}

fn iv(x: f64) -> Interval {
    interval!(x, x).unwrap()
}

fn pair_to_complex(pair: Pair) -> Complex {
    Complex {
        re: pair.0,
        im: pair.1,
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

fn eps_shape_scale(eps: f64) -> f64 {
    eps.powf(1.0 / 28.0)
}

fn scalar_target(eps: f64, degree: usize) -> f64 {
    RADIAL_HALF * eps.powf(1.0 / (degree as f64))
}

fn interval_from_pair(pair: [f64; 2]) -> Interval {
    interval!(pair[0], pair[1]).unwrap()
}

fn scale_interval(x: Interval, c: f64) -> Interval {
    if c >= 0.0 {
        x * iv(c)
    } else {
        x * iv(c)
    }
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

fn interval_square(x: Interval) -> Interval {
    let lo = x.inf();
    let hi = x.sup();
    if lo <= 0.0 && hi >= 0.0 {
        interval!(0.0, lo.abs().max(hi.abs()).powi(2)).unwrap()
    } else {
        let a = lo * lo;
        let b = hi * hi;
        interval!(a.min(b), a.max(b)).unwrap()
    }
}

fn abs_sq_minus_one(z: CInterval) -> Interval {
    interval_square(z.re) + interval_square(z.im) - iv(1.0)
}

fn eval_poly_root_affine(z: Complex, roots: &[CInterval]) -> Interval {
    let zc = ci_fixed(z);
    let mut acc = ci_one();
    for root in roots {
        acc = ci_mul(acc, ci_sub(zc, *root));
    }
    abs_sq_minus_one(acc)
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
        return interval!(0.0, 1.0).unwrap();
    }
    let out = fa / denom;
    let lo = out.inf().max(0.0);
    let hi = out.sup().min(1.0);
    if lo > hi {
        interval!(0.0, 1.0).unwrap()
    } else {
        interval!(lo, hi).unwrap()
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

fn interval_marching_upper(cfg: &WorkConfig, roots: &[CInterval]) -> MarchingStats {
    let step = 2.0 * cfg.extent / (cfg.res as f64);
    let crude_penalty_upper = 2.0 * 2.0_f64.sqrt() * step;
    let rows: Vec<(f64, usize, usize, usize, usize, f64)> = (0..cfg.res)
        .into_par_iter()
        .map(|ix| {
            let x0 = -cfg.extent + (ix as f64) * step;
            let x1 = x0 + step;
            let mut total_upper = 0.0;
            let mut active_cells = 0usize;
            let mut definite_case_cells = 0usize;
            let mut uncertain_corner_cells = 0usize;
            let mut ambiguous_avg_cells = 0usize;
            let mut max_cell_upper = 0.0f64;

            for iy in 0..cfg.res {
                let y0 = -cfg.extent + (iy as f64) * step;
                let y1 = y0 + step;
                let fsw = eval_poly_root_affine(Complex { re: x0, im: y0 }, roots);
                let fse = eval_poly_root_affine(Complex { re: x1, im: y0 }, roots);
                let fne = eval_poly_root_affine(Complex { re: x1, im: y1 }, roots);
                let fnw = eval_poly_root_affine(Complex { re: x0, im: y1 }, roots);
                let signs = [sign(fsw), sign(fse), sign(fne), sign(fnw)];
                if signs.iter().all(|s| *s == '+') || signs.iter().all(|s| *s == '-') {
                    continue;
                }
                active_cells += 1;
                let upper = if signs.contains(&'?') {
                    uncertain_corner_cells += 1;
                    crude_penalty_upper
                } else {
                    let case = (if signs[0] == '+' { 1 } else { 0 })
                        | ((if signs[1] == '+' { 1 } else { 0 }) << 1)
                        | ((if signs[2] == '+' { 1 } else { 0 }) << 2)
                        | ((if signs[3] == '+' { 1 } else { 0 }) << 3);
                    let (upper, avg_known) =
                        cell_upper_from_fixed_case(case, (fsw, fse, fne, fnw), x0, y0, step);
                    definite_case_cells += 1;
                    if !avg_known {
                        ambiguous_avg_cells += 1;
                    }
                    upper
                };
                total_upper += upper;
                max_cell_upper = max_cell_upper.max(upper);
            }

            (
                total_upper,
                active_cells,
                definite_case_cells,
                uncertain_corner_cells,
                ambiguous_avg_cells,
                max_cell_upper,
            )
        })
        .collect();

    rows.into_iter().fold(
        MarchingStats {
            length_upper: 0.0,
            active_cells: 0,
            definite_case_cells: 0,
            uncertain_corner_cells: 0,
            ambiguous_avg_cells: 0,
            max_cell_upper: 0.0,
        },
        |mut acc, row| {
            acc.length_upper += row.0;
            acc.active_cells += row.1;
            acc.definite_case_cells += row.2;
            acc.uncertain_corner_cells += row.3;
            acc.ambiguous_avg_cells += row.4;
            acc.max_cell_upper = acc.max_cell_upper.max(row.5);
            acc
        },
    )
}

fn root_intervals_from_cell(
    cfg: &WorkConfig,
    u0_pair: [f64; 2],
    u1_pair: [f64; 2],
) -> Vec<CInterval> {
    let scale = eps_shape_scale(cfg.eps);
    let a = interval_from_pair(u0_pair);
    let b = interval_from_pair(u1_pair);
    (0..cfg.degree)
        .map(|i| {
            complex_interval_linear(
                cfg.base_roots[i],
                scale,
                a,
                cfg.u0_direction[i],
                b,
                cfg.u1_direction[i],
            )
        })
        .collect()
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

fn evaluate_subcell(
    cfg: &WorkConfig,
    sub_i: usize,
    sub_j: usize,
    u0_pair: [f64; 2],
    u1_pair: [f64; 2],
) -> SubcellRow {
    let roots = root_intervals_from_cell(cfg, u0_pair, u1_pair);
    let stats = interval_marching_upper(cfg, &roots);
    let deficit_lower = cfg.lstar_lower - stats.length_upper;
    let margin_lower = deficit_lower - cfg.target;
    SubcellRow {
        sub_i,
        sub_j,
        u0_interval: u0_pair,
        u1_interval: u1_pair,
        length_upper: stats.length_upper,
        deficit_lower,
        margin_lower,
        pass: margin_lower >= 0.0,
        active_cells: stats.active_cells,
        definite_case_cells: stats.definite_case_cells,
        uncertain_corner_cells: stats.uncertain_corner_cells,
        ambiguous_avg_cells: stats.ambiguous_avg_cells,
    }
}

fn validate_and_normalize(input: InputConfig) -> Result<WorkConfig, String> {
    if input.degree != DEGREE {
        return Err(format!(
            "this binary is scoped to degree {DEGREE}; got {}",
            input.degree
        ));
    }
    if input.base_roots.len() != input.degree
        || input.u0_direction.len() != input.degree
        || input.u1_direction.len() != input.degree
    {
        return Err(
            "base_roots, u0_direction, and u1_direction must each contain 14 entries".to_string(),
        );
    }
    if input.res == 0 || input.subdivision == 0 {
        return Err("res and subdivision must be positive".to_string());
    }
    Ok(WorkConfig {
        experiment_id: input.experiment_id,
        degree: input.degree,
        eps: input.eps,
        extent: input.extent,
        res: input.res,
        subdivision: input.subdivision,
        lstar_lower: input.lstar_lower,
        target: input
            .target
            .unwrap_or_else(|| scalar_target(input.eps, input.degree)),
        base_roots: input.base_roots.into_iter().map(pair_to_complex).collect(),
        u0_direction: input
            .u0_direction
            .into_iter()
            .map(pair_to_complex)
            .collect(),
        u1_direction: input
            .u1_direction
            .into_iter()
            .map(pair_to_complex)
            .collect(),
        u0_interval: input.u0_interval,
        u1_interval: input.u1_interval,
    })
}

fn parse_args() -> Result<(Option<PathBuf>, PathBuf, bool), String> {
    let args: Vec<String> = env::args().collect();
    let mut input = None;
    let mut out_dir = PathBuf::from(".");
    let mut quiet = false;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--input" => {
                if i + 1 >= args.len() {
                    return Err("--input requires a path".to_string());
                }
                input = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--outdir" => {
                if i + 1 >= args.len() {
                    return Err("--outdir requires a path".to_string());
                }
                out_dir = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "--quiet" => {
                quiet = true;
                i += 1;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
    }
    Ok((input, out_dir, quiet))
}

fn load_config(input_path: Option<PathBuf>) -> Result<(WorkConfig, String), String> {
    match input_path {
        Some(path) => {
            let raw = fs::read_to_string(&path)
                .map_err(|err| format!("failed to read input {}: {err}", path.display()))?;
            let input: InputConfig = serde_json::from_str(&raw)
                .map_err(|err| format!("failed to parse input {}: {err}", path.display()))?;
            Ok((validate_and_normalize(input)?, path.display().to_string()))
        }
        None => Ok((
            validate_and_normalize(InputConfig {
                experiment_id: default_experiment_id(),
                degree: default_degree(),
                eps: default_eps(),
                extent: default_extent(),
                res: default_res(),
                subdivision: default_subdivision(),
                lstar_lower: default_lstar_lower(),
                target: None,
                base_roots: default_base_roots(),
                u0_direction: default_u0_direction(),
                u1_direction: default_u1_direction(),
                u0_interval: default_u0_interval(),
                u1_interval: default_u1_interval(),
            })?,
            "built-in Python SUBDIV8 reproduction cell".to_string(),
        )),
    }
}

fn unix_timestamp_string() -> String {
    match SystemTime::now().duration_since(UNIX_EPOCH) {
        Ok(dur) => format!("{}", dur.as_secs()),
        Err(_) => "0".to_string(),
    }
}

fn write_report(
    cfg: &WorkConfig,
    status: &str,
    failure_count: usize,
    min_margin_lower: f64,
    max_length_upper: f64,
    elapsed_secs: f64,
    path: &PathBuf,
) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Root-Affine Rust Cell Reproduction\n\n\
Experiment: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Subdivision: `{}` x `{}`\n\
- Failure count: `{}`\n\
- Minimum margin lower bound: `{}`\n\
- Maximum length upper bound: `{}`\n\
- Elapsed seconds: `{}`\n\n\
## Claim Ceiling\n\n\
This is a continuous grid-oracle cell certificate, not exact lemniscate certification, not a Lean theorem, and not a proof of Erdős #114.\n\n\
## Reproduction Gate\n\n\
The gate compares against Python `{}` with `max_length_upper <= {} + {}` and `min_margin_lower >= {} - {}`.\n",
        cfg.experiment_id,
        status,
        cfg.subdivision,
        cfg.subdivision,
        failure_count,
        min_margin_lower,
        max_length_upper,
        elapsed_secs,
        PYTHON_PARENT_PACKET,
        PYTHON_MAX_LENGTH_UPPER,
        REPRO_TOLERANCE,
        PYTHON_MIN_MARGIN_LOWER,
        REPRO_TOLERANCE,
    );
    fs::write(path, report)
        .map_err(|err| format!("failed to write report {}: {err}", path.display()))
}

fn sha256_file(path: &PathBuf) -> Result<String, String> {
    let bytes = fs::read(path)
        .map_err(|err| format!("failed to read {} for sha256: {err}", path.display()))?;
    Ok(sha256::digest(bytes))
}

fn run() -> Result<(), String> {
    let started = Instant::now();
    let (input_path, out_dir, quiet) = parse_args()?;
    let (cfg, input_source) = load_config(input_path)?;
    fs::create_dir_all(&out_dir).map_err(|err| {
        format!(
            "failed to create output directory {}: {err}",
            out_dir.display()
        )
    })?;

    let u0_splits = split_interval(cfg.u0_interval, cfg.subdivision);
    let u1_splits = split_interval(cfg.u1_interval, cfg.subdivision);
    let mut tasks = Vec::with_capacity(cfg.subdivision * cfg.subdivision);
    for (i, u0_pair) in u0_splits.iter().enumerate() {
        for (j, u1_pair) in u1_splits.iter().enumerate() {
            tasks.push((i, j, *u0_pair, *u1_pair));
        }
    }

    let mut rows: Vec<SubcellRow> = tasks
        .into_par_iter()
        .map(|(i, j, u0_pair, u1_pair)| evaluate_subcell(&cfg, i, j, u0_pair, u1_pair))
        .collect();
    rows.sort_by_key(|row| (row.sub_i, row.sub_j));

    let failure_count = rows.iter().filter(|row| !row.pass).count();
    let min_margin_lower = rows
        .iter()
        .map(|row| row.margin_lower)
        .fold(f64::INFINITY, f64::min);
    let max_length_upper = rows
        .iter()
        .map(|row| row.length_upper)
        .fold(f64::NEG_INFINITY, f64::max);
    let finite_bounds = min_margin_lower.is_finite() && max_length_upper.is_finite();
    let reproduction_gate_pass = finite_bounds
        && failure_count == 0
        && max_length_upper <= PYTHON_MAX_LENGTH_UPPER + REPRO_TOLERANCE
        && min_margin_lower >= PYTHON_MIN_MARGIN_LOWER - REPRO_TOLERANCE;
    let status = if reproduction_gate_pass {
        "RUST_ROOT_AFFINE_SUBDIV8_REPRO_PASS"
    } else if !finite_bounds {
        "RUST_ROOT_AFFINE_SUBDIV8_NONFINITE"
    } else if failure_count == 0 {
        "RUST_ROOT_AFFINE_SUBDIV8_CELL_PASS_GATE_MISMATCH"
    } else {
        "RUST_ROOT_AFFINE_SUBDIV8_FAIL"
    };
    let elapsed_secs = started.elapsed().as_secs_f64();

    let result = json!({
        "experiment_id": cfg.experiment_id,
        "timestamp_unix": unix_timestamp_string(),
        "input_source": input_source,
        "python_reference_packet": PYTHON_PARENT_PACKET,
        "status": status,
        "parameters": {
            "degree": cfg.degree,
            "eps": cfg.eps,
            "extent": cfg.extent,
            "res": cfg.res,
            "subdivision": cfg.subdivision,
            "subcell_count": rows.len(),
            "scalar_target": "D14(eps,s) >= 12 * eps^(1/14)"
        },
        "selected_cell": {
            "u0_interval": cfg.u0_interval,
            "u1_interval": cfg.u1_interval
        },
        "lstar_lower": cfg.lstar_lower,
        "target": cfg.target,
        "failure_count": failure_count,
        "min_margin_lower": min_margin_lower,
        "max_length_upper": max_length_upper,
        "python_reference": {
            "max_length_upper": PYTHON_MAX_LENGTH_UPPER,
            "min_margin_lower": PYTHON_MIN_MARGIN_LOWER,
            "tolerance": REPRO_TOLERANCE
        },
        "reproduction_gate_pass": reproduction_gate_pass,
        "rows": rows,
        "claim_ceiling": "Rust reproduction/hardening of a continuous grid-oracle cell certificate. This is not exact lemniscate certification, not a Lean theorem, and not a proof of Erdős #114.",
        "next_blocker": if reproduction_gate_pass {
            "Lift to neighboring cells only after preserving this reproduction gate; exact lemniscate-length enclosure remains separate."
        } else {
            "Resolve Rust/Python interval-order mismatch or subdivide further before using this as the canonical reproduction."
        },
        "elapsed_secs": elapsed_secs
    });

    let result_path = out_dir.join(format!("{}_RESULTS.json", cfg.experiment_id));
    let report_path = out_dir.join(format!("{}_REPORT.md", cfg.experiment_id));
    let sha_path = out_dir.join(format!("{}_RESULTS.sha256", cfg.experiment_id));
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
    write_report(
        &cfg,
        status,
        failure_count,
        min_margin_lower,
        max_length_upper,
        elapsed_secs,
        &report_path,
    )?;
    let sha = sha256_file(&result_path)?;
    fs::write(
        &sha_path,
        format!(
            "{sha}  {}\n",
            result_path.file_name().unwrap().to_string_lossy()
        ),
    )
    .map_err(|err| format!("failed to write sha file {}: {err}", sha_path.display()))?;

    let summary = json!({
        "experiment_id": cfg.experiment_id,
        "status": status,
        "failure_count": failure_count,
        "min_margin_lower": min_margin_lower,
        "max_length_upper": max_length_upper,
        "reproduction_gate_pass": reproduction_gate_pass,
        "result": result_path,
        "report": report_path,
        "sha256": sha,
        "elapsed_secs": elapsed_secs
    });
    if !quiet {
        println!("{}", serde_json::to_string_pretty(&summary).unwrap());
    }
    Ok(())
}

fn main() {
    if let Err(err) = run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}
