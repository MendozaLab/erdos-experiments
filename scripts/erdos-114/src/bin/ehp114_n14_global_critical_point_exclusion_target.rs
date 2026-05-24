//! EHP #114 n=14 global critical-point exclusion target.
//!
//! This diagnostic consumes the L18 third-order collar residuals and classifies
//! each region by whether root-affine interval arithmetic can already exclude
//! level-set crossing, certify a nonzero gradient lower bound, or must leave the
//! region as a critical-point candidate.
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
    "EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01";
const SOURCE_EXPERIMENT_ID: &str =
    "EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01";
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
struct SourceRegion {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    x: [f64; 2],
    y: [f64; 2],
    source_reason: String,
}

#[derive(Clone, Debug, Serialize)]
struct RegionClassification {
    source_index: usize,
    source_ownership_key: String,
    split_path: String,
    source_reason: String,
    x_interval: [f64; 2],
    y_interval: [f64; 2],
    classification: String,
    reason: String,
    f_interval: IntervalPair,
    fx_interval: IntervalPair,
    fy_interval: IntervalPair,
    fx_abs_lower: f64,
    fy_abs_lower: f64,
    gradient_lower_bound_candidate: f64,
    dominant_axis: String,
}

#[derive(Clone, Debug)]
struct Config {
    source: PathBuf,
    out_dir: PathBuf,
    region_limit: usize,
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

fn interval_avoids_zero(x: Interval) -> bool {
    x.inf() > 0.0 || x.sup() < 0.0
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
    let mut p1 = CInterval {
        re: iv(0.0),
        im: iv(0.0),
    };
    for root in roots {
        let factor = ci_sub(z, *root);
        let next_p1 = CInterval {
            re: p1.re * factor.re - p1.im * factor.im + p.re,
            im: p1.re * factor.im + p1.im * factor.re + p.im,
        };
        let next_p = ci_mul(p, factor);
        p = next_p;
        p1 = next_p1;
    }
    (p, p1)
}

fn gradient_intervals(p: CInterval, p1: CInterval) -> (Interval, Interval) {
    let q = ci_mul(p1, ci_conj(p));
    (q.re * iv(2.0), q.im * iv(-2.0))
}

fn f_fx_fy(z: CInterval, roots: &[CInterval]) -> (Interval, Interval, Interval) {
    let (p, p1) = eval_p_p1(z, roots);
    let f = abs_sq(p) - iv(1.0);
    let (fx, fy) = gradient_intervals(p, p1);
    (f, fx, fy)
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

fn parse_source_region(idx: usize, value: &serde_json::Value) -> Result<SourceRegion, String> {
    Ok(SourceRegion {
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
    })
}

fn classify_region(region: &SourceRegion, roots: &[CInterval]) -> RegionClassification {
    let (f, fx, fy) = f_fx_fy(ci_box(region.x, region.y), roots);
    let fx_abs_lower = interval_abs_lower(fx);
    let fy_abs_lower = interval_abs_lower(fy);
    let gradient_lower_bound_candidate =
        (fx_abs_lower * fx_abs_lower + fy_abs_lower * fy_abs_lower).sqrt();
    let dominant_axis = if fx_abs_lower >= fy_abs_lower {
        "x"
    } else {
        "y"
    };
    let (classification, reason) = if interval_avoids_zero(f) {
        (
            "excluded_region".to_string(),
            "F interval avoids zero on the residual box".to_string(),
        )
    } else if gradient_lower_bound_candidate > 0.0 {
        (
            "regular_region".to_string(),
            format!("|grad F| has interval lower bound candidate {gradient_lower_bound_candidate}"),
        )
    } else {
        (
            "critical_candidate_region".to_string(),
            "both Fx and Fy intervals contain zero under root-affine uncertainty".to_string(),
        )
    };

    RegionClassification {
        source_index: region.source_index,
        source_ownership_key: region.source_ownership_key.clone(),
        split_path: region.split_path.clone(),
        source_reason: region.source_reason.clone(),
        x_interval: region.x,
        y_interval: region.y,
        classification,
        reason,
        f_interval: ipair(f),
        fx_interval: ipair(fx),
        fy_interval: ipair(fy),
        fx_abs_lower,
        fy_abs_lower,
        gradient_lower_bound_candidate,
        dominant_axis: dominant_axis.to_string(),
    }
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut cfg = Config {
        source: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01_RESULTS.json"),
        out_dir: PathBuf::from("../../Erdos114/validated_length/EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01"),
        region_limit: 0,
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
            "--region-limit" => {
                if i + 1 >= args.len() {
                    return Err("--region-limit requires an integer".to_string());
                }
                cfg.region_limit = args[i + 1]
                    .parse::<usize>()
                    .map_err(|err| format!("failed to parse --region-limit: {err}"))?;
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

fn status_for(regular_count: usize, excluded_count: usize) -> &'static str {
    if regular_count > 0 {
        "GLOBAL_CRITICAL_POINT_TARGET_REGULAR_REGIONS_FOUND"
    } else if excluded_count > 0 {
        "GLOBAL_CRITICAL_POINT_TARGET_EXCLUDES_ONLY"
    } else {
        "GLOBAL_CRITICAL_POINT_TARGET_ALL_CRITICAL_CANDIDATES"
    }
}

fn first_failed_condition(classifications: &[RegionClassification]) -> String {
    classifications
        .iter()
        .find(|row| row.classification == "critical_candidate_region")
        .map(|row| {
            format!(
                "{}:{} remains critical candidate because {}",
                row.source_ownership_key, row.split_path, row.reason
            )
        })
        .unwrap_or_else(|| "none".to_string())
}

fn write_report(result: &serde_json::Value, path: &PathBuf) -> Result<(), String> {
    let report = format!(
        "# EHP114 n=14 Global Critical-Point Exclusion Target\n\n\
Experiment: `{}`\n\n\
Source: `{}`\n\n\
## Verdict\n\n\
- Status: `{}`\n\
- Processed regions: `{}`\n\
- Regular regions: `{}`\n\
- Critical-candidate regions: `{}`\n\
- Excluded regions: `{}`\n\
- Gradient lower-bound candidate: `{}`\n\
- Total validated length upper: `{}`\n\
- Exact length cap: `{}`\n\
- Margin to cap: `{}`\n\
- First failed condition: `{}`\n\n\
## Interpretation\n\n\
This diagnostic asks whether any L18 residual region already admits a global \
root-affine lower bound for `|grad F|` in ordinary coordinates. It does not \
promote branch length. Regular regions are candidates for a future analytic \
domain decomposition; critical-candidate regions identify where a sharper \
critical-point theorem is still needed.\n\n\
## Claim Ceiling\n\n\
{}\n",
        result["experiment_id"]
            .as_str()
            .unwrap_or(DEFAULT_EXPERIMENT_ID),
        SOURCE_EXPERIMENT_ID,
        result["status"],
        result["processed_region_count"],
        result["regular_region_count"],
        result["critical_candidate_region_count"],
        result["excluded_region_count"],
        result["gradient_lower_bound_candidate"],
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
    let total_validated_length_upper = source_json["total_validated_length_upper"]
        .as_f64()
        .ok_or_else(|| "source missing total_validated_length_upper".to_string())?;
    let cap = source_json["exact_length_cap"]
        .as_f64()
        .ok_or_else(|| "source missing exact_length_cap".to_string())?;
    let regions_json = source_json["remaining_unresolved"]
        .as_array()
        .ok_or_else(|| "source missing remaining_unresolved".to_string())?;
    let limit = if cfg.region_limit == 0 {
        regions_json.len()
    } else {
        cfg.region_limit.min(regions_json.len())
    };
    let roots = root_intervals(cfg.cell);
    let mut classifications = Vec::<RegionClassification>::new();
    for (idx, value) in regions_json.iter().take(limit).enumerate() {
        let region = parse_source_region(idx, value)?;
        classifications.push(classify_region(&region, &roots));
    }

    let excluded_count = classifications
        .iter()
        .filter(|row| row.classification == "excluded_region")
        .count();
    let regular_count = classifications
        .iter()
        .filter(|row| row.classification == "regular_region")
        .count();
    let critical_candidate_count = classifications
        .iter()
        .filter(|row| row.classification == "critical_candidate_region")
        .count();
    let regular_by_x_count = classifications
        .iter()
        .filter(|row| row.classification == "regular_region" && row.dominant_axis == "x")
        .count();
    let regular_by_y_count = classifications
        .iter()
        .filter(|row| row.classification == "regular_region" && row.dominant_axis == "y")
        .count();
    let gradient_lower_bound_candidate = classifications
        .iter()
        .filter(|row| row.classification == "regular_region")
        .map(|row| row.gradient_lower_bound_candidate)
        .fold(f64::INFINITY, f64::min);
    let gradient_lower_bound_candidate = if gradient_lower_bound_candidate.is_finite() {
        gradient_lower_bound_candidate
    } else {
        0.0
    };
    let worst_critical_candidate = classifications
        .iter()
        .filter(|row| row.classification == "critical_candidate_region")
        .min_by(|a, b| {
            a.gradient_lower_bound_candidate
                .partial_cmp(&b.gradient_lower_bound_candidate)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .cloned();
    let status = status_for(regular_count, excluded_count);
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
            "region_limit": cfg.region_limit,
            "source_region_field": "remaining_unresolved",
            "classification_rule": "exclude if F avoids zero; else regular if interval lower bound for |grad F| is positive; else critical candidate",
            "ordinary_coordinate_rule": "use Fx and Fy intervals over each residual x/y box with full root-affine uncertainty",
            "non_promotion_rule": "no branch length is promoted by this target diagnostic"
        },
        "processed_region_count": classifications.len(),
        "source_region_count": regions_json.len(),
        "source_unprocessed_region_count": regions_json.len().saturating_sub(limit),
        "gradient_lower_bound_candidate": gradient_lower_bound_candidate,
        "critical_candidate_region_count": critical_candidate_count,
        "regular_region_count": regular_count,
        "regular_by_x_count": regular_by_x_count,
        "regular_by_y_count": regular_by_y_count,
        "excluded_region_count": excluded_count,
        "total_validated_length_upper": total_validated_length_upper,
        "exact_length_cap": cap,
        "margin_to_cap": cap - total_validated_length_upper,
        "first_failed_condition": first_failed_condition(&classifications),
        "classifications": classifications,
        "worst_critical_candidate": worst_critical_candidate,
        "claim_ceiling": "Local n=14 hard-cell global critical-point exclusion target diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.",
        "next_blocker": if regular_count > 0 {
            "Regular regions exist under ordinary-coordinate root-affine gradient lower bounds. Next proof-facing step is to build an analytic residual-domain decomposition that turns these regular boxes into owned collars while isolating the remaining critical candidates."
        } else {
            "No ordinary-coordinate regular regions were found. Next proof-facing step is a sharper analytic critical-point theorem or a different coordinate/domain decomposition."
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
            "experiment_id": cfg.experiment_id,
            "status": status,
            "processed_region_count": result["processed_region_count"],
            "regular_region_count": result["regular_region_count"],
            "critical_candidate_region_count": result["critical_candidate_region_count"],
            "excluded_region_count": result["excluded_region_count"],
            "gradient_lower_bound_candidate": result["gradient_lower_bound_candidate"],
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
