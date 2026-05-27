use inari::{interval, Interval};

const ATOM_LOCATION: f64 = -1.0;
const SUPPORT_LEFT: f64 = 0.8044618744671781;
const SUPPORT_RIGHT: f64 = 1.0;
const ATOM_MASS: f64 = 0.8245217514114138;
const DEGREE: usize = 10_000;
const CLOUD_NODE_COUNT: usize = 24;
const ROOT_POSITION_RADIUS: f64 = 1.0e-12;
const CROSSING_PANEL_RADIUS: f64 = 1.0e-8;
const GAP_COVER_SUBDIVISIONS: usize = 256;
const TAIL_COVER_SUBDIVISIONS: usize = 512;
const LEFT_TAIL_LO: f64 = -3.0;
const RIGHT_TAIL_HI: f64 = 3.0;
const DERIVATIVE_TOL: f64 = 1.0e-6;

#[derive(Clone, Copy)]
struct RootMass {
    x: f64,
    m: usize,
}

#[derive(Clone, Copy)]
struct RootBox {
    x: Interval,
    m: usize,
}

#[derive(Clone, Copy)]
struct Crossing {
    center: f64,
    interval_lo: f64,
    interval_hi: f64,
    left_u: Interval,
    right_u: Interval,
    u_prime: Interval,
    min_abs_u_prime: f64,
    sign_change: bool,
    orientation: &'static str,
    status: &'static str,
}

#[derive(Clone, Copy)]
struct CoverResult {
    interval_lo: f64,
    interval_hi: f64,
    subdivisions: usize,
    min_log_abs_f_lower: f64,
    max_log_abs_f_upper: f64,
    failed_subpanels: usize,
    zero_factor_subpanels: usize,
    status: &'static str,
}

#[derive(Clone, Copy)]
struct GapGuard {
    interval_lo: f64,
    interval_hi: f64,
    critical_point: f64,
    critical_u: Interval,
    status: &'static str,
}

fn iv(lo: f64, hi: f64) -> Interval {
    interval!(lo, hi).expect("valid interval")
}

fn point(x: f64) -> Interval {
    iv(x, x)
}

fn contains_zero(x: Interval) -> bool {
    x.inf() <= 0.0 && x.sup() >= 0.0
}

fn zero_excluded(x: Interval) -> bool {
    x.sup() < 0.0 || x.inf() > 0.0
}

fn abs_lower(x: Interval) -> f64 {
    if contains_zero(x) {
        0.0
    } else {
        x.inf().abs().min(x.sup().abs())
    }
}

fn interval_pair(x: Interval) -> String {
    format!("[{:.17e},{:.17e}]", x.inf(), x.sup())
}

fn build_projection() -> Vec<RootMass> {
    let c = 0.5 * (SUPPORT_LEFT + SUPPORT_RIGHT);
    let r = 0.5 * (SUPPORT_RIGHT - SUPPORT_LEFT);
    let h = ((SUPPORT_LEFT - ATOM_LOCATION) * (SUPPORT_RIGHT - ATOM_LOCATION)).sqrt();
    let atom_m = (ATOM_MASS * DEGREE as f64).round() as usize;
    let cloud_m = DEGREE - atom_m;

    let mut nodes: Vec<f64> = Vec::with_capacity(CLOUD_NODE_COUNT);
    let mut raw: Vec<f64> = Vec::with_capacity(CLOUD_NODE_COUNT);
    for j in 0..CLOUD_NODE_COUNT {
        let theta = std::f64::consts::PI * (j as f64 + 0.5) / CLOUD_NODE_COUNT as f64;
        let x = c + r * theta.cos();
        let density_factor = 1.0 - ATOM_MASS * h / (x - ATOM_LOCATION);
        nodes.push(x);
        raw.push(density_factor.max(0.0));
    }

    let raw_sum: f64 = raw.iter().sum();
    let exact: Vec<f64> = raw.iter().map(|v| v / raw_sum * cloud_m as f64).collect();
    let mut mult: Vec<usize> = exact.iter().map(|v| v.floor() as usize).collect();
    let assigned: usize = mult.iter().sum();
    let remainder = cloud_m - assigned;
    let mut order: Vec<usize> = (0..CLOUD_NODE_COUNT).collect();
    order.sort_by(|&i, &j| {
        let ri = exact[i] - mult[i] as f64;
        let rj = exact[j] - mult[j] as f64;
        rj.partial_cmp(&ri).unwrap_or(std::cmp::Ordering::Equal)
    });
    for &idx in order.iter().take(remainder) {
        mult[idx] += 1;
    }

    let mut roots = vec![RootMass {
        x: ATOM_LOCATION,
        m: atom_m,
    }];
    for (x, m) in nodes.iter().zip(mult.iter()) {
        if *m > 0 {
            roots.push(RootMass { x: *x, m: *m });
        }
    }
    roots.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal));
    roots
}

fn root_boxes(roots: &[RootMass]) -> Vec<RootBox> {
    roots
        .iter()
        .map(|r| RootBox {
            x: iv(r.x - ROOT_POSITION_RADIUS, r.x + ROOT_POSITION_RADIUS),
            m: r.m,
        })
        .collect()
}

fn weight(r: RootMass) -> f64 {
    r.m as f64 / DEGREE as f64
}

fn u_f64(x: f64, roots: &[RootMass]) -> f64 {
    roots
        .iter()
        .map(|r| weight(*r) * (x - r.x).abs().max(1.0e-300).ln())
        .sum()
}

fn up_f64(x: f64, roots: &[RootMass]) -> f64 {
    roots.iter().map(|r| weight(*r) / (x - r.x)).sum()
}

fn u_box(x: Interval, roots: &[RootBox]) -> (Interval, bool) {
    let mut zero_factor = false;
    let mut acc = point(0.0);
    for root in roots {
        let dx = x - root.x;
        if contains_zero(dx) {
            zero_factor = true;
        } else {
            acc = acc + point(root.m as f64 / DEGREE as f64) * dx.abs().ln();
        }
    }
    (acc, zero_factor)
}

fn up_box(x: Interval, roots: &[RootBox]) -> (Interval, bool) {
    let mut zero_factor = false;
    let mut acc = point(0.0);
    for root in roots {
        let dx = x - root.x;
        if contains_zero(dx) {
            zero_factor = true;
        } else {
            acc = acc + point(root.m as f64 / DEGREE as f64) / dx;
        }
    }
    (acc, zero_factor)
}

fn bisect(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    let lo_positive = f(lo) > 0.0;
    for _ in 0..90 {
        let mid = 0.5 * (lo + hi);
        if (f(mid) > 0.0) == lo_positive {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn enumerate_crossings(roots: &[RootMass]) -> Vec<f64> {
    let mut crossings: Vec<f64> = Vec::new();
    let left = bisect(|x| u_f64(x, roots), roots[0].x - 3.0, roots[0].x - 1.0e-12);
    crossings.push(left);

    for i in 0..roots.len() - 1 {
        let a = roots[i].x;
        let b = roots[i + 1].x;
        if b - a <= 1.0e-12 {
            continue;
        }
        let lo = a + 1.0e-12;
        let hi = b - 1.0e-12;
        let crit = bisect(|x| up_f64(x, roots), lo, hi);
        if u_f64(crit, roots) > 0.0 {
            crossings.push(bisect(|x| u_f64(x, roots), lo, crit));
            crossings.push(bisect(|x| u_f64(x, roots), crit, hi));
        }
    }

    let right = bisect(
        |x| u_f64(x, roots),
        roots[roots.len() - 1].x + 1.0e-12,
        roots[roots.len() - 1].x + 3.0,
    );
    crossings.push(right);
    crossings.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    crossings
}

fn certify_crossing(center: f64, roots: &[RootBox]) -> Crossing {
    let lo = center - CROSSING_PANEL_RADIUS;
    let hi = center + CROSSING_PANEL_RADIUS;
    let (left_u, left_zero) = u_box(point(lo), roots);
    let (right_u, right_zero) = u_box(point(hi), roots);
    let (u_prime, deriv_zero_factor) = up_box(iv(lo, hi), roots);
    let min_abs_u_prime = abs_lower(u_prime);
    let sign_change = !left_zero
        && !right_zero
        && ((left_u.sup() < 0.0 && right_u.inf() > 0.0)
            || (left_u.inf() > 0.0 && right_u.sup() < 0.0));
    let orientation = if left_u.sup() < 0.0 && right_u.inf() > 0.0 {
        "negative_to_positive"
    } else if left_u.inf() > 0.0 && right_u.sup() < 0.0 {
        "positive_to_negative"
    } else {
        "inconclusive"
    };
    let pass = sign_change
        && !deriv_zero_factor
        && zero_excluded(u_prime)
        && min_abs_u_prime > DERIVATIVE_TOL;
    Crossing {
        center,
        interval_lo: lo,
        interval_hi: hi,
        left_u,
        right_u,
        u_prime,
        min_abs_u_prime,
        sign_change,
        orientation,
        status: if pass { "PASS" } else { "INCONCLUSIVE" },
    }
}

fn positive_cover(lo: f64, hi: f64, subdivisions: usize, roots: &[RootBox]) -> CoverResult {
    let mut min_lower = f64::INFINITY;
    let mut max_upper = f64::NEG_INFINITY;
    let mut failed = 0usize;
    let mut zero_factor_subpanels = 0usize;
    for k in 0..subdivisions {
        let a = lo + (hi - lo) * (k as f64) / (subdivisions as f64);
        let b = lo + (hi - lo) * ((k + 1) as f64) / (subdivisions as f64);
        let (u, zero_factor) = u_box(iv(a, b), roots);
        if zero_factor {
            zero_factor_subpanels += 1;
            failed += 1;
        } else {
            min_lower = min_lower.min(u.inf());
            max_upper = max_upper.max(u.sup());
            if !(u.inf() > 0.0) {
                failed += 1;
            }
        }
    }
    CoverResult {
        interval_lo: lo,
        interval_hi: hi,
        subdivisions,
        min_log_abs_f_lower: min_lower,
        max_log_abs_f_upper: max_upper,
        failed_subpanels: failed,
        zero_factor_subpanels,
        status: if failed == 0 { "CERTIFIED" } else { "INCONCLUSIVE" },
    }
}

fn critical_gap_guards(roots: &[RootMass], boxes: &[RootBox]) -> Vec<GapGuard> {
    let mut guards: Vec<GapGuard> = Vec::new();
    for i in 0..roots.len() - 1 {
        let lo = roots[i].x + 1.0e-12;
        let hi = roots[i + 1].x - 1.0e-12;
        if hi <= lo {
            continue;
        }
        let crit = bisect(|x| up_f64(x, roots), lo, hi);
        if u_f64(crit, roots) > 0.0 {
            let (critical_u, zero_factor) = u_box(point(crit), boxes);
            guards.push(GapGuard {
                interval_lo: roots[i].x,
                interval_hi: roots[i + 1].x,
                critical_point: crit,
                critical_u,
                status: if !zero_factor && critical_u.inf() > 0.0 {
                    "CERTIFIED"
                } else {
                    "INCONCLUSIVE"
                },
            });
        }
    }
    guards
}

fn cover_json(prefix: &str, idx: usize, c: CoverResult, final_row: bool) -> String {
    format!(
        "    {{\"{}_id\":\"{}_{}\",\"interval\":[{:.17e},{:.17e}],\"subdivisions\":{},\"min_log_abs_f_lower\":{:.17e},\"max_log_abs_f_upper\":{:.17e},\"failed_subpanels\":{},\"zero_factor_subpanels\":{},\"status\":\"{}\"}}{}",
        prefix,
        prefix,
        idx,
        c.interval_lo,
        c.interval_hi,
        c.subdivisions,
        c.min_log_abs_f_lower,
        c.max_log_abs_f_upper,
        c.failed_subpanels,
        c.zero_factor_subpanels,
        c.status,
        if final_row { "" } else { "," }
    )
}

fn gap_guard_json(idx: usize, g: GapGuard, final_row: bool) -> String {
    format!(
        "    {{\"gap_guard_id\":\"gap_guard_{}\",\"root_gap_interval\":[{:.17e},{:.17e}],\"critical_point\":{:.17e},\"critical_u\":{},\"method\":\"strict_log_concavity_plus_positive_critical_value\",\"status\":\"{}\"}}{}",
        idx,
        g.interval_lo,
        g.interval_hi,
        g.critical_point,
        interval_pair(g.critical_u),
        g.status,
        if final_row { "" } else { "," }
    )
}

fn main() {
    let roots = build_projection();
    let boxes = root_boxes(&roots);
    let centers = enumerate_crossings(&roots);
    let crossings: Vec<Crossing> = centers.iter().map(|&x| certify_crossing(x, &boxes)).collect();
    let crossing_pass = crossings.iter().all(|c| c.status == "PASS");
    let component_count = crossings.len() / 2;

    let mut component_pass = true;
    let mut component_rows: Vec<String> = Vec::new();
    for k in 0..component_count {
        let left = crossings[2 * k];
        let right = crossings[2 * k + 1];
        let pass = left.status == "PASS"
            && right.status == "PASS"
            && left.orientation == "positive_to_negative"
            && right.orientation == "negative_to_positive"
            && left.interval_hi < right.interval_lo;
        if !pass {
            component_pass = false;
        }
        component_rows.push(format!(
            "    {{\"component_id\":\"comp_{}\",\"crossing_ids\":[\"u_root_{}\",\"u_root_{}\"],\"interval\":[{:.17e},{:.17e}],\"length_lower\":{:.17e},\"status\":\"{}\"}}{}",
            k,
            2 * k,
            2 * k + 1,
            left.interval_hi,
            right.interval_lo,
            right.interval_lo - left.interval_hi,
            if pass { "CERTIFIED_FROM_BOUNDARY_PAIR" } else { "INCONCLUSIVE" },
            if k + 1 == component_count { "" } else { "," }
        ));
    }

    let mut gap_covers: Vec<CoverResult> = Vec::new();
    for k in 0..component_count - 1 {
        let lo = crossings[2 * k + 1].interval_hi;
        let hi = crossings[2 * k + 2].interval_lo;
        gap_covers.push(positive_cover(lo, hi, GAP_COVER_SUBDIVISIONS, &boxes));
    }
    let direct_gap_cover_pass = gap_covers.iter().all(|c| c.status == "CERTIFIED");
    let gap_guards = critical_gap_guards(&roots, &boxes);
    let gap_guard_pass = gap_guards.len() == component_count - 1
        && gap_guards.iter().all(|g| g.status == "CERTIFIED");

    let tail_covers = vec![
        positive_cover(LEFT_TAIL_LO, crossings[0].interval_lo, TAIL_COVER_SUBDIVISIONS, &boxes),
        positive_cover(
            crossings[crossings.len() - 1].interval_hi,
            RIGHT_TAIL_HI,
            TAIL_COVER_SUBDIVISIONS,
            &boxes,
        ),
    ];
    let tail_pass = tail_covers.iter().all(|c| c.status == "CERTIFIED");
    let pass = crossing_pass && component_pass && gap_guard_pass && tail_pass;
    let min_abs_u_prime = crossings
        .iter()
        .map(|c| c.min_abs_u_prime)
        .fold(f64::INFINITY, f64::min);
    let min_gap_log_lower = gap_covers
        .iter()
        .map(|c| c.min_log_abs_f_lower)
        .fold(f64::INFINITY, f64::min);
    let min_tail_log_lower = tail_covers
        .iter()
        .map(|c| c.min_log_abs_f_lower)
        .fold(f64::INFINITY, f64::min);

    println!("{{");
    println!("  \"pass\":{},", pass);
    println!("  \"backend\":\"rust-inari\",");
    println!("  \"rust_binary\":\"src/bin/phi_k_root_factor_full_panel_replay_backend.rs\",");
    println!("  \"degree\":{},", DEGREE);
    println!("  \"root_count\":{},", boxes.len());
    println!("  \"root_box_radius\":{:.17e},", ROOT_POSITION_RADIUS);
    println!("  \"crossing_panel_radius\":{:.17e},", CROSSING_PANEL_RADIUS);
    println!("  \"crossing_count\":{},", crossings.len());
    println!("  \"component_count\":{},", component_count);
    println!("  \"gap_count\":{},", gap_covers.len());
    println!("  \"tail_count\":{},", tail_covers.len());
    println!("  \"crossing_pass\":{},", crossing_pass);
    println!("  \"component_pass\":{},", component_pass);
    println!("  \"gap_guard_pass\":{},", gap_guard_pass);
    println!("  \"direct_gap_cover_pass\":{},", direct_gap_cover_pass);
    println!("  \"tail_pass\":{},", tail_pass);
    println!("  \"min_abs_u_prime\":{:.17e},", min_abs_u_prime);
    println!("  \"min_gap_log_abs_f_lower\":{:.17e},", min_gap_log_lower);
    println!("  \"min_tail_log_abs_f_lower\":{:.17e},", min_tail_log_lower);
    println!("  \"claim_ceiling\":\"full log-domain root-factor panel replay for fixed PUBLIC-CANDIDATE finite projection; local evidence only\",");
    println!("  \"crossings\":[");
    for (idx, c) in crossings.iter().enumerate() {
        println!(
            "    {{\"crossing_id\":\"u_root_{}\",\"center\":{:.17e},\"interval\":[{:.17e},{:.17e}],\"u_left_endpoint\":{},\"u_right_endpoint\":{},\"u_prime_interval\":{},\"min_abs_u_prime\":{:.17e},\"sign_change\":{},\"orientation\":\"{}\",\"status\":\"{}\"}}{}",
            idx,
            c.center,
            c.interval_lo,
            c.interval_hi,
            interval_pair(c.left_u),
            interval_pair(c.right_u),
            interval_pair(c.u_prime),
            c.min_abs_u_prime,
            c.sign_change,
            c.orientation,
            c.status,
            if idx + 1 == crossings.len() { "" } else { "," }
        );
    }
    println!("  ],");
    println!("  \"components\":[");
    for row in component_rows {
        println!("{}", row);
    }
    println!("  ],");
    println!("  \"gap_guards\":[");
    for (idx, g) in gap_guards.iter().enumerate() {
        println!("{}", gap_guard_json(idx, *g, idx + 1 == gap_guards.len()));
    }
    println!("  ],");
    println!("  \"gap_covers\":[");
    for (idx, c) in gap_covers.iter().enumerate() {
        println!("{}", cover_json("gap", idx, *c, idx + 1 == gap_covers.len()));
    }
    println!("  ],");
    println!("  \"tail_covers\":[");
    for (idx, c) in tail_covers.iter().enumerate() {
        println!("{}", cover_json("tail", idx, *c, idx + 1 == tail_covers.len()));
    }
    println!("  ]");
    println!("}}");
}
