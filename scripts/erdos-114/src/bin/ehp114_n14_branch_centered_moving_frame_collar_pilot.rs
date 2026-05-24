//! EHP #114 n=14 branch-centered moving-frame collar pilot.
//!
//! This diagnostic consumes the L18 third-order collar artifact and tests only
//! the Class A wall failures: regions where the wall budget would allow a
//! positive center-strip bound if the branch were re-centered correctly.
//!
//! The proof-facing question is whether the raw center-strip interval
//! `sup_r |F(0,r)|` can be replaced by a branch-centered moving-frame quadratic
//! estimate. This pilot does not promote branch length. It only records whether
//! branch centering, tangent cancellation, and the quadratic center-strip budget
//! are plausible under the current constants.
//!
//! This remains local n=14 hard-cell work. It is not a proof of Erdos #114,
//! not a global n=14 certificate, and not an exact lemniscate-length
//! certificate.

use inari::{interval, Interval};
use serde::Serialize;
use serde_json::json;
use std::cmp::Ordering;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-BRANCH-CENTERED-MOVING-FRAME-COLLAR-PILOT-HARD-CELL-20260506-01";
const SOURCE_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01";
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
struct Direction {
    x: f64,
    y: f64,
}

#[derive(Clone, Debug)]
struct CandidateRegion {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    x: [f64; 2],
    y: [f64; 2],
    normal: Direction,
    tangent: Direction,
    normal_interval: [f64; 2],
    tangent_interval: [f64; 2],
    raw_c0: f64,
    normal_remainder: f64,
    wall_rhs: f64,
    critical_margin: f64,
    f_tt_abs_upper: f64,
    third_directional_upper: f64,
}

#[derive(Clone, Debug, Serialize)]
struct MovingFrameAudit {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    normal_interval: [f64; 2],
    tangent_interval: [f64; 2],
    normal: Direction,
    tangent: Direction,
    class_a_allowed_c0: f64,
    raw_center_strip_c0: f64,
    branch_center_left_f_interval: IntervalPair,
    branch_center_right_f_interval: IntervalPair,
    branch_center_validated: bool,
    moving_frame_tangent_zero: bool,
    quadratic_center_strip_c0: f64,
    quadratic_center_strip_margin: f64,
    quadratic_center_strip_pass: bool,
    strict_promotion_eligible: bool,
    reason: String,
}

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    out_dir: PathBuf,
    class_a_limit: usize,
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

fn ci_one() -> CInterval {
    CInterval {
        re: iv(1.0),
        im: iv(0.0),
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
        iwrap(0.0, lo.abs().max(hi.abs()).powi(2))
    } else {
        let a = lo * lo;
        let b = hi * hi;
        iwrap(a.min(b), a.max(b))
    }
}

fn interval_avoids_zero(x: Interval) -> bool {
    x.inf() > 0.0 || x.sup() < 0.0
}

fn opposite_signed(a: Interval, b: Interval) -> bool {
    (a.inf() > 0.0 && b.sup() < 0.0) || (a.sup() < 0.0 && b.inf() > 0.0)
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

fn eval_p(z: CInterval, roots: &[CInterval]) -> CInterval {
    let mut p = ci_one();
    for root in roots {
        p = ci_mul(p, ci_sub(z, *root));
    }
    p
}

fn f_interval(z: CInterval, roots: &[CInterval]) -> Interval {
    abs_sq(eval_p(z, roots)) - iv(1.0)
}

fn midpoint(pair: [f64; 2]) -> f64 {
    0.5 * (pair[0] + pair[1])
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

fn value_direction(value: &serde_json::Value, field: &str) -> Result<Direction, String> {
    let obj = value
        .get(field)
        .and_then(|v| v.as_object())
        .ok_or_else(|| format!("missing direction field {field}"))?;
    Ok(Direction {
        x: obj
            .get("x")
            .and_then(|v| v.as_f64())
            .ok_or_else(|| format!("{field}.x not numeric"))?,
        y: obj
            .get("y")
            .and_then(|v| v.as_f64())
            .ok_or_else(|| format!("{field}.y not numeric"))?,
    })
}

fn parse_candidate(value: &serde_json::Value) -> Result<Option<CandidateRegion>, String> {
    if value.get("reason").and_then(|v| v.as_str()) != Some("third_order_wall_separation_failed") {
        return Ok(None);
    }
    let inequality = value
        .get("inequality")
        .and_then(|v| v.as_object())
        .ok_or_else(|| "wall failure missing inequality object".to_string())?;
    let wall_rhs = inequality
        .get("wall_rhs")
        .and_then(|v| v.as_f64())
        .ok_or_else(|| "missing inequality.wall_rhs".to_string())?;
    let normal_remainder = inequality
        .get("normal_remainder_bound")
        .and_then(|v| v.as_f64())
        .ok_or_else(|| "missing inequality.normal_remainder_bound".to_string())?;
    let allowed_c0 = wall_rhs - normal_remainder;
    if allowed_c0 <= 0.0 {
        return Ok(None);
    }

    Ok(Some(CandidateRegion {
        source_index: value
            .get("source_index")
            .and_then(|v| v.as_u64())
            .ok_or_else(|| "missing source_index".to_string())? as usize,
        source_ownership_key: value
            .get("source_ownership_key")
            .and_then(|v| v.as_str())
            .ok_or_else(|| "missing source_ownership_key".to_string())?
            .to_string(),
        split_path: value
            .get("split_path")
            .and_then(|v| v.as_str())
            .ok_or_else(|| "missing split_path".to_string())?
            .to_string(),
        x: value_interval(value, "x_interval")?,
        y: value_interval(value, "y_interval")?,
        normal: value_direction(value, "normal")?,
        tangent: value_direction(value, "tangent")?,
        normal_interval: value_interval(value, "normal_interval")?,
        tangent_interval: value_interval(value, "tangent_interval")?,
        raw_c0: inequality
            .get("center_strip_f_abs_upper")
            .and_then(|v| v.as_f64())
            .ok_or_else(|| "missing inequality.center_strip_f_abs_upper".to_string())?,
        normal_remainder,
        wall_rhs,
        critical_margin: inequality
            .get("critical_exclusion_margin")
            .and_then(|v| v.as_f64())
            .ok_or_else(|| "missing inequality.critical_exclusion_margin".to_string())?,
        f_tt_abs_upper: inequality
            .get("f_tt_abs_upper")
            .and_then(|v| v.as_f64())
            .ok_or_else(|| "missing inequality.f_tt_abs_upper".to_string())?,
        third_directional_upper: inequality
            .get("third_directional_upper")
            .and_then(|v| v.as_f64())
            .ok_or_else(|| "missing inequality.third_directional_upper".to_string())?,
    }))
}

fn candidate_priority(a: &CandidateRegion, b: &CandidateRegion) -> Ordering {
    let a_key = if a.source_ownership_key == "2286:0" {
        0
    } else {
        1
    };
    let b_key = if b.source_ownership_key == "2286:0" {
        0
    } else {
        1
    };
    (a_key, a.source_index, a.split_path.as_str()).cmp(&(
        b_key,
        b.source_index,
        b.split_path.as_str(),
    ))
}

fn branch_center_audit(candidate: &CandidateRegion, roots: &[CInterval]) -> MovingFrameAudit {
    let cx = midpoint(candidate.x);
    let cy = midpoint(candidate.y);
    let left_s = candidate.normal_interval[0];
    let right_s = candidate.normal_interval[1];
    let left_f = f_interval(
        ci_affine(
            cx,
            cy,
            candidate.normal,
            candidate.tangent,
            iv(left_s),
            iv(0.0),
        ),
        roots,
    );
    let right_f = f_interval(
        ci_affine(
            cx,
            cy,
            candidate.normal,
            candidate.tangent,
            iv(right_s),
            iv(0.0),
        ),
        roots,
    );
    let branch_center_validated = interval_avoids_zero(left_f)
        && interval_avoids_zero(right_f)
        && opposite_signed(left_f, right_f)
        && candidate.critical_margin > 0.0;

    let tangent_radius = candidate
        .tangent_interval
        .iter()
        .fold(0.0_f64, |acc, x| acc.max(x.abs()));
    let quadratic_center_strip_c0 = 0.5 * candidate.f_tt_abs_upper * tangent_radius.powi(2)
        + (candidate.third_directional_upper / 6.0) * tangent_radius.powi(3);
    let allowed_c0 = candidate.wall_rhs - candidate.normal_remainder;
    let moving_frame_tangent_zero = branch_center_validated;
    let quadratic_center_strip_margin = allowed_c0 - quadratic_center_strip_c0;
    let quadratic_center_strip_pass =
        moving_frame_tangent_zero && quadratic_center_strip_margin > 0.0;

    let reason = if !branch_center_validated {
        "branch_center_not_validated_on_normal_axis"
    } else if !quadratic_center_strip_pass {
        "quadratic_center_strip_exceeds_allowed_budget"
    } else {
        "quadratic_center_strip_passes_diagnostic_no_length_promotion"
    };

    MovingFrameAudit {
        source_index: candidate.source_index,
        source_ownership_key: candidate.source_ownership_key.clone(),
        split_path: candidate.split_path.clone(),
        x_interval: candidate.x,
        y_interval: candidate.y,
        normal_interval: candidate.normal_interval,
        tangent_interval: candidate.tangent_interval,
        normal: candidate.normal,
        tangent: candidate.tangent,
        class_a_allowed_c0: allowed_c0,
        raw_center_strip_c0: candidate.raw_c0,
        branch_center_left_f_interval: ipair(left_f),
        branch_center_right_f_interval: ipair(right_f),
        branch_center_validated,
        moving_frame_tangent_zero,
        quadratic_center_strip_c0,
        quadratic_center_strip_margin,
        quadratic_center_strip_pass,
        strict_promotion_eligible: false,
        reason: reason.to_string(),
    }
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01_RESULTS.json"),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-BRANCH-CENTERED-MOVING-FRAME-COLLAR-PILOT-HARD-CELL-20260506-01"),
        class_a_limit: 0,
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
            "--class-a-limit" => {
                if i + 1 >= args.len() {
                    return Err("--class-a-limit requires an integer".to_string());
                }
                cfg.class_a_limit = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --class-a-limit: {err}"))?;
                i += 2;
            }
            other => return Err(format!("unknown argument: {other}")),
        }
    }
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

fn status_for(quadratic_pass_count: usize, branch_center_fail_count: usize) -> &'static str {
    if quadratic_pass_count > 0 {
        "MOVING_FRAME_COLLAR_PILOT_QUADRATIC_CENTER_STRIP_PASSES_DIAGNOSTIC"
    } else if branch_center_fail_count > 0 {
        "MOVING_FRAME_COLLAR_PILOT_FAIL_BRANCH_CENTERING"
    } else {
        "MOVING_FRAME_COLLAR_PILOT_FAIL_CENTER_STRIP"
    }
}

fn first_failed_condition(audits: &[MovingFrameAudit]) -> String {
    audits
        .iter()
        .find(|a| !a.quadratic_center_strip_pass)
        .map(|a| a.reason.clone())
        .unwrap_or_else(|| "none".to_string())
}

fn write_report(result: &serde_json::Value, path: &PathBuf) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Branch-Centered Moving-Frame Collar Pilot\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Processed Class A regions: `{}`\n\
- Validated center points: `{}`\n\
- Moving-frame tangent-zero count: `{}`\n\
- Quadratic center-strip passes: `{}`\n\
- Remaining Class A unresolved: `{}`\n\
- Total validated length upper: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
This pilot tests the branch-centered moving-frame theorem shape on the Class A \
wall failures from L18. It validates a branch center by checking opposite signs \
on the normal axis at `r=0`, then audits whether the branch-centered quadratic \
center-strip estimate would fit under the remaining wall budget. It does not \
promote branch length.\n\n\
## Claim Ceiling\n\n\
{}\n",
        EXPERIMENT_ID,
        SOURCE_EXPERIMENT_ID,
        result["status"],
        result["processed_class_a_count"],
        result["validated_center_point_count"],
        result["moving_frame_tangent_zero_count"],
        result["quadratic_center_strip_pass_count"],
        result["remaining_class_a_unresolved_count"],
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
    let cap = source_json["exact_length_cap"]
        .as_f64()
        .ok_or_else(|| "source missing exact_length_cap".to_string())?;
    let total_validated_length_upper = source_json["total_validated_length_upper"]
        .as_f64()
        .ok_or_else(|| "source missing total_validated_length_upper".to_string())?;
    let regions = source_json["remaining_unresolved"]
        .as_array()
        .ok_or_else(|| "source missing remaining_unresolved".to_string())?;

    let mut candidates = Vec::<CandidateRegion>::new();
    let mut normal_budget_failed_count = 0usize;
    for value in regions {
        if value.get("reason").and_then(|v| v.as_str())
            == Some("third_order_wall_separation_failed")
        {
            let inequality = value
                .get("inequality")
                .and_then(|v| v.as_object())
                .ok_or_else(|| "wall failure missing inequality object".to_string())?;
            let wall_rhs = inequality
                .get("wall_rhs")
                .and_then(|v| v.as_f64())
                .ok_or_else(|| "missing wall_rhs".to_string())?;
            let normal_remainder = inequality
                .get("normal_remainder_bound")
                .and_then(|v| v.as_f64())
                .ok_or_else(|| "missing normal_remainder_bound".to_string())?;
            if wall_rhs - normal_remainder <= 0.0 {
                normal_budget_failed_count += 1;
            }
        }
        if let Some(candidate) = parse_candidate(value)? {
            candidates.push(candidate);
        }
    }
    candidates.sort_by(candidate_priority);
    let limit = if cfg.class_a_limit == 0 {
        candidates.len()
    } else {
        cfg.class_a_limit.min(candidates.len())
    };
    let roots = root_intervals();
    let audits: Vec<MovingFrameAudit> = candidates
        .iter()
        .take(limit)
        .map(|candidate| branch_center_audit(candidate, &roots))
        .collect();

    let validated_center_point_count = audits.iter().filter(|a| a.branch_center_validated).count();
    let moving_frame_tangent_zero_count = audits
        .iter()
        .filter(|a| a.moving_frame_tangent_zero)
        .count();
    let quadratic_center_strip_pass_count = audits
        .iter()
        .filter(|a| a.quadratic_center_strip_pass)
        .count();
    let branch_center_fail_count = audits.iter().filter(|a| !a.branch_center_validated).count();
    let remaining_class_a_unresolved_count = audits
        .len()
        .saturating_sub(quadratic_center_strip_pass_count);
    let status = status_for(quadratic_center_strip_pass_count, branch_center_fail_count);
    let first_failed = first_failed_condition(&audits);
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
            "class_a_limit": cfg.class_a_limit,
            "source_filter": "remaining_unresolved rows with reason third_order_wall_separation_failed and wall_rhs - normal_remainder_bound > 0",
            "priority_rule": "process source_ownership_key 2286:0 first",
            "branch_center_validation_rule": "F on normal-axis endpoints at r=0 avoids zero and has opposite signs under root-affine intervals",
            "moving_frame_tangent_zero_rule": "if branch center is validated and critical margin remains positive, tangent is taken perpendicular to grad F at the branch center",
            "quadratic_center_strip_rule": "0.5 * M_tt * R^2 + (T3/6) * R^3 < wall_rhs - normal_remainder",
            "strict_non_promotion_rule": "no branch length is promoted by this pilot"
        },
        "source_wall_failure_count": candidates.len() + normal_budget_failed_count,
        "source_class_a_count": candidates.len(),
        "processed_class_a_count": audits.len(),
        "source_unprocessed_class_a_count": candidates.len().saturating_sub(audits.len()),
        "normal_budget_failed_count": normal_budget_failed_count,
        "validated_center_point_count": validated_center_point_count,
        "moving_frame_tangent_zero_count": moving_frame_tangent_zero_count,
        "quadratic_center_strip_pass_count": quadratic_center_strip_pass_count,
        "remaining_class_a_unresolved_count": remaining_class_a_unresolved_count,
        "promoted_branch_length_upper": 0.0,
        "total_validated_length_upper": total_validated_length_upper,
        "exact_length_cap": cap,
        "margin_to_cap": cap - total_validated_length_upper,
        "first_failed_condition": first_failed,
        "audits": audits,
        "claim_ceiling": "Local n=14 hard-cell branch-centered moving-frame collar pilot only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. Shadow signature, not universal law.",
        "next_blocker": if quadratic_center_strip_pass_count > 0 {
            "At least one Class A center-strip inequality passes under branch-centered moving-frame assumptions. Next step is to turn the branch center and moving frame into a full certified collar with ownership and length promotion gates."
        } else {
            "No Class A region passed the branch-centered moving-frame pilot. If key 2286:0 failed centering, retire this local-collar route and move to a global analytic critical-point exclusion theorem."
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
            "experiment_id": EXPERIMENT_ID,
            "status": status,
            "processed_class_a_count": result["processed_class_a_count"],
            "validated_center_point_count": result["validated_center_point_count"],
            "moving_frame_tangent_zero_count": result["moving_frame_tangent_zero_count"],
            "quadratic_center_strip_pass_count": result["quadratic_center_strip_pass_count"],
            "remaining_class_a_unresolved_count": result["remaining_class_a_unresolved_count"],
            "first_failed_condition": result["first_failed_condition"],
            "total_validated_length_upper": result["total_validated_length_upper"],
            "margin_to_cap": result["margin_to_cap"],
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
