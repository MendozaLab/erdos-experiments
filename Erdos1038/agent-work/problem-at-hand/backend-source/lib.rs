//! erdos1038_fast — High-performance sublevel measure evaluator for Erdős #1038
//!
//! For monic polynomials of degree n with all real roots in [-1,1]:
//! Computes |{x ∈ ℝ : |f(x)| < 1}| = |{x : Σ wᵢ·ln|x-pᵢ| < 0}|
//!
//! Core algorithm: analytical boundary detection via bisection on the
//! logarithmic potential g(x) = Σ wᵢ·ln|x-pᵢ|. No numerical integration.
//! Each evaluation is O(N·B) where N = distinct root count, B = bisection depth (64).

use inari::{interval, Interval};
use pyo3::prelude::*;
use rayon::prelude::*;

// ============================================================================
// PRNG: xoshiro256** (Blackman & Vigna) — fast, high-quality, reproducible
// ============================================================================

struct Rng {
    s: [u64; 4],
}

impl Rng {
    fn seed(seed: u64) -> Self {
        // SplitMix64 initialization
        let mut z = seed;
        let mut state = [0u64; 4];
        for st in &mut state {
            z = z.wrapping_add(0x9e3779b97f4a7c15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
            *st = x ^ (x >> 31);
        }
        Self { s: state }
    }

    fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }

    fn gen_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }

    fn gen_usize(&mut self, n: usize) -> usize {
        (self.next_u64() as usize) % n
    }
}

// ============================================================================
// Core: sublevel measure computation
// ============================================================================

const BISECT_ITERS: usize = 64;

/// Bisection root-finding. Assumes f(lo) and f(hi) have opposite signs.
fn bisect(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    let f_lo_positive = f(lo) > 0.0;
    for _ in 0..BISECT_ITERS {
        let mid = 0.5 * (lo + hi);
        if (f(mid) > 0.0) == f_lo_positive {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Weighted logarithmic potential g(x) = Σ wᵢ·ln|x - pᵢ|
#[inline]
fn g_eval(x: f64, pos: &[f64], wt: &[f64]) -> f64 {
    pos.iter()
        .zip(wt.iter())
        .map(|(&p, &w)| w * (x - p).abs().max(1e-300).ln())
        .sum()
}

/// Derivative g'(x) = Σ wᵢ/(x - pᵢ)
#[inline]
fn gp_eval(x: f64, pos: &[f64], wt: &[f64]) -> f64 {
    pos.iter().zip(wt.iter()).map(|(&p, &w)| w / (x - p)).sum()
}

/// Compute |{x : Σ wᵢ·ln|x-pᵢ| < 0}| via analytical boundary detection.
///
/// For positions p₁ < p₂ < ... < pₘ with positive weights summing to 1:
/// - g is concave on each interval (pᵢ, pᵢ₊₁) with a unique maximum
/// - g is strictly decreasing on (-∞, p₁) and strictly increasing on (pₘ, +∞)
/// - Sublevel set = union of intervals around each root where g < 0
///
/// Total boundary-finding work: O(m) bisections, each O(64) evaluations.
pub fn sublevel_measure_core(positions: &[f64], weights: &[f64]) -> f64 {
    let n = positions.len();
    if n == 0 {
        return 0.0;
    }

    // Sort and aggregate coincident roots (within epsilon)
    let mut pw: Vec<(f64, f64)> = positions
        .iter()
        .copied()
        .zip(weights.iter().copied())
        .collect();
    pw.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let eps = 1e-12;
    let mut agg: Vec<(f64, f64)> = Vec::with_capacity(n);
    for (p, w) in &pw {
        if let Some(last) = agg.last_mut() {
            if (*p - last.0).abs() < eps {
                last.1 += w;
                continue;
            }
        }
        agg.push((*p, *w));
    }

    let pos: Vec<f64> = agg.iter().map(|x| x.0).collect();
    let wt: Vec<f64> = agg.iter().map(|x| x.1).collect();
    let m = pos.len();

    if m == 1 {
        // g(x) = w·ln|x-p|, sublevel: |x-p| < 1, measure = 2
        return 2.0;
    }

    let g = |x: f64| g_eval(x, &pos, &wt);
    let gp = |x: f64| gp_eval(x, &pos, &wt);

    let mut total = 0.0;

    // ── Left external zero ──
    // g → +∞ as x → -∞, g → -∞ near pos[0]
    // Find bracket [far, near] where g changes sign
    {
        let near = pos[0] - 1e-14;
        let mut far = pos[0] - 2.0;
        for _ in 0..60 {
            if g(far) > 0.0 {
                break;
            }
            far = pos[0] - 2.0 * (pos[0] - far); // double distance
        }
        let z_left = bisect(&g, far, near);
        total += pos[0] - z_left;
    }

    // ── Right external zero ──
    // g → -∞ near pos[m-1], g → +∞ as x → +∞
    {
        let near = pos[m - 1] + 1e-14;
        let mut far = pos[m - 1] + 2.0;
        for _ in 0..60 {
            if g(far) > 0.0 {
                break;
            }
            far = pos[m - 1] + 2.0 * (far - pos[m - 1]);
        }
        let z_right = bisect(&g, near, far);
        total += z_right - pos[m - 1];
    }

    // ── Inner intervals ──
    // For each (pᵢ, pᵢ₊₁): g is concave, g → -∞ at both endpoints
    // Has unique maximum at critical point where g' = 0
    for i in 0..m - 1 {
        let a = pos[i];
        let b = pos[i + 1];
        if b - a < 1e-14 {
            continue;
        }

        // g' → +∞ near a, g' → -∞ near b, g' strictly decreasing (g'' < 0)
        let crit = bisect(&gp, a + 1e-14, b - 1e-14);
        let g_crit = g(crit);

        if g_crit <= 0.0 {
            // Entire interval is sublevel
            total += b - a;
        } else {
            // g has two zeros: z1 ∈ (a, crit) and z2 ∈ (crit, b)
            // Sublevel parts: (a, z1) and (z2, b)
            let z1 = bisect(&g, a + 1e-14, crit);
            let z2 = bisect(&g, crit, b - 1e-14);
            total += (z1 - a) + (b - z2);
        }
    }

    total
}

// ============================================================================
// Interval-arithmetic panel-constraint extraction (inari, IEEE 1788)
// ============================================================================

/// Find sublevel-set zeros (left external, right external) for the aggregated
/// position/weight arrays. Mirrors the bisection used in `sublevel_measure_core`,
/// but returns the zero locations instead of accumulating measure.
fn find_sublevel_external_zeros(pos: &[f64], wt: &[f64]) -> (f64, f64) {
    let m = pos.len();
    debug_assert!(m >= 1);
    let g = |x: f64| g_eval(x, pos, wt);

    // Left external zero
    let z_left = {
        let near = pos[0] - 1e-14;
        let mut far = pos[0] - 2.0;
        for _ in 0..60 {
            if g(far) > 0.0 {
                break;
            }
            far = pos[0] - 2.0 * (pos[0] - far);
        }
        bisect(&g, far, near)
    };

    // Right external zero
    let z_right = {
        let near = pos[m - 1] + 1e-14;
        let mut far = pos[m - 1] + 2.0;
        for _ in 0..60 {
            if g(far) > 0.0 {
                break;
            }
            far = pos[m - 1] + 2.0 * (far - pos[m - 1]);
        }
        bisect(&g, near, far)
    };

    (z_left, z_right)
}

/// Build candidate panels on the legal-support side: the two intervals
/// (z_left - 1.0, z_left) and (z_right, z_right + 1.0).
/// Roots in those intervals (and their `root_epsilon` neighborhoods) are excised.
fn build_legal_support_panels(
    z_left: f64,
    z_right: f64,
    panel_count: usize,
    positions: &[f64],
    root_epsilon: f64,
) -> Vec<(f64, f64)> {
    // We split panel_count evenly between left and right sides.
    let left_count = panel_count / 2;
    let right_count = panel_count - left_count;

    let mut grid_pts: Vec<f64> = Vec::with_capacity(panel_count + 2);

    // Left side: from z_left - 1.0 up to z_left (legal support: x < z_left)
    if left_count > 0 {
        let a = z_left - 1.0;
        let b = z_left;
        for i in 0..=left_count {
            let t = i as f64 / left_count as f64;
            grid_pts.push(a + t * (b - a));
        }
    }
    // Convert the cumulative left-grid into panels
    let mut panels: Vec<(f64, f64)> = Vec::new();
    if left_count > 0 {
        for i in 0..left_count {
            panels.push((grid_pts[i], grid_pts[i + 1]));
        }
        grid_pts.clear();
    }

    // Right side: from z_right up to z_right + 1.0 (legal support: x > z_right)
    if right_count > 0 {
        let a = z_right;
        let b = z_right + 1.0;
        for i in 0..=right_count {
            let t = i as f64 / right_count as f64;
            grid_pts.push(a + t * (b - a));
        }
        for i in 0..right_count {
            panels.push((grid_pts[i], grid_pts[i + 1]));
        }
    }

    // Excise root-epsilon neighborhoods. If a panel contains a root within
    // root_epsilon, split it into the legal sub-pieces (cursor walk).
    let mut sorted_roots = positions.to_vec();
    sorted_roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut final_panels: Vec<(f64, f64)> = Vec::with_capacity(panels.len());
    for (p_lo, p_hi) in panels {
        let affected: Vec<f64> = sorted_roots
            .iter()
            .copied()
            .filter(|&r| (p_lo < r + root_epsilon) && (p_hi > r - root_epsilon))
            .collect();
        if affected.is_empty() {
            final_panels.push((p_lo, p_hi));
            continue;
        }
        let mut cursor = p_lo;
        for r in &affected {
            let left_end = r - root_epsilon;
            if left_end > cursor + 1e-12 {
                final_panels.push((cursor, left_end));
            }
            cursor = r + root_epsilon;
        }
        if cursor < p_hi - 1e-12 {
            final_panels.push((cursor, p_hi));
        }
    }

    final_panels
}

/// Core extraction routine.
/// Returns (panels, A_ineq, baseline_slack, n_active, n_inactive).
/// Panels containing a root (zero in `|x - r|`) are filtered out before return.
pub fn extract_panel_constraints_core(
    positions: &[f64],
    weights: &[f64],
    panel_count: usize,
    root_epsilon: f64,
) -> (Vec<(f64, f64)>, Vec<Vec<f64>>, Vec<f64>, usize, usize) {
    // Aggregate coincident roots (same eps tolerance as sublevel_measure_core)
    let mut pw: Vec<(f64, f64)> = positions
        .iter()
        .copied()
        .zip(weights.iter().copied())
        .collect();
    pw.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let eps_agg = 1e-12;
    let mut agg: Vec<(f64, f64)> = Vec::with_capacity(positions.len());
    for (p, w) in &pw {
        if let Some(last) = agg.last_mut() {
            if (*p - last.0).abs() < eps_agg {
                last.1 += w;
                continue;
            }
        }
        agg.push((*p, *w));
    }
    let pos: Vec<f64> = agg.iter().map(|x| x.0).collect();
    let wt: Vec<f64> = agg.iter().map(|x| x.1).collect();

    if pos.is_empty() {
        return (vec![], vec![], vec![], 0, 0);
    }

    let (z_left, z_right) = if pos.len() == 1 {
        // Single root: sublevel = (p-1, p+1); legal support is outside that ball
        (pos[0] - 1.0, pos[0] + 1.0)
    } else {
        find_sublevel_external_zeros(&pos, &wt)
    };

    let candidate_panels =
        build_legal_support_panels(z_left, z_right, panel_count, positions, root_epsilon);

    // Per-panel × per-root interval lower-bound of log|x - r_i|.
    // We index columns over the ORIGINAL positions array (not aggregated)
    // because the caller's weights are aligned to that ordering.
    let n_roots = positions.len();
    let mut final_panels: Vec<(f64, f64)> = Vec::with_capacity(candidate_panels.len());
    let mut a_ineq: Vec<Vec<f64>> = Vec::with_capacity(candidate_panels.len());

    for (p_lo, p_hi) in candidate_panels {
        // Defensive: skip degenerate panels
        if !(p_hi > p_lo) {
            continue;
        }
        let panel_iv = match interval!(p_lo, p_hi) {
            Ok(iv) => iv,
            Err(_) => continue,
        };

        let mut row: Vec<f64> = Vec::with_capacity(n_roots);
        let mut touches_root = false;
        for &r in positions.iter() {
            let r_iv = match interval!(r, r) {
                Ok(iv) => iv,
                Err(_) => {
                    touches_root = true;
                    break;
                }
            };
            let diff_abs: Interval = (panel_iv - r_iv).abs();
            if diff_abs.contains(0.0) {
                // Panel touches this root — interval log is unbounded below.
                touches_root = true;
                break;
            }
            let log_iv = diff_abs.ln();
            let lo = log_iv.inf();
            if !lo.is_finite() {
                touches_root = true;
                break;
            }
            row.push(lo);
        }

        if !touches_root && row.len() == n_roots {
            final_panels.push((p_lo, p_hi));
            a_ineq.push(row);
        }
    }

    // Baseline slack = Σ w_i · A_ineq[p][i]
    let baseline_slack: Vec<f64> = a_ineq
        .iter()
        .map(|row| {
            row.iter()
                .zip(weights.iter())
                .map(|(a, w)| a * w)
                .sum::<f64>()
        })
        .collect();

    let mut n_active = 0usize;
    let mut n_inactive = 0usize;
    for &s in &baseline_slack {
        if s <= 1e-10 {
            n_active += 1;
        } else {
            n_inactive += 1;
        }
    }

    (final_panels, a_ineq, baseline_slack, n_active, n_inactive)
}

// ============================================================================
// Differential Evolution (rand/1/bin)
// ============================================================================

fn select_three(rng: &mut Rng, pop_size: usize, exclude: usize) -> (usize, usize, usize) {
    let mut a = rng.gen_usize(pop_size);
    while a == exclude {
        a = rng.gen_usize(pop_size);
    }
    let mut b = rng.gen_usize(pop_size);
    while b == exclude || b == a {
        b = rng.gen_usize(pop_size);
    }
    let mut c = rng.gen_usize(pop_size);
    while c == exclude || c == a || c == b {
        c = rng.gen_usize(pop_size);
    }
    (a, b, c)
}

fn differential_evolution(
    objective: &(impl Fn(&[f64]) -> f64 + Sync),
    bounds: &[(f64, f64)],
    pop_size: usize,
    max_iter: usize,
    f_mut: f64,
    cr: f64,
    seed: u64,
) -> (Vec<f64>, f64) {
    let dim = bounds.len();
    let pop_size = pop_size.max(4); // need at least 4 for select_three
    let mut rng = Rng::seed(seed);

    let mut pop: Vec<Vec<f64>> = (0..pop_size)
        .map(|_| {
            bounds
                .iter()
                .map(|&(lo, hi)| rng.gen_range(lo, hi))
                .collect()
        })
        .collect();
    let mut fitness: Vec<f64> = pop.iter().map(|x| objective(x)).collect();

    for _ in 0..max_iter {
        for i in 0..pop_size {
            let (a, b, c) = select_three(&mut rng, pop_size, i);
            let j_rand = rng.gen_usize(dim);

            let trial: Vec<f64> = (0..dim)
                .map(|j| {
                    if j == j_rand || rng.next_f64() < cr {
                        let v = pop[a][j] + f_mut * (pop[b][j] - pop[c][j]);
                        v.clamp(bounds[j].0, bounds[j].1)
                    } else {
                        pop[i][j]
                    }
                })
                .collect();

            let trial_fit = objective(&trial);
            if trial_fit < fitness[i] {
                pop[i] = trial;
                fitness[i] = trial_fit;
            }
        }
    }

    let best_idx = fitness
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);

    (pop[best_idx].clone(), fitness[best_idx])
}

// ============================================================================
// Nelder-Mead simplex optimizer
// ============================================================================

fn nelder_mead(
    objective: &(impl Fn(&[f64]) -> f64 + Sync),
    initial: &[f64],
    step: f64,
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, f64) {
    let n = initial.len();
    if n == 0 {
        return (vec![], objective(&[]));
    }
    let (alpha, gamma, rho, sigma) = (1.0, 2.0, 0.5, 0.5);

    // Build initial simplex
    let mut simplex: Vec<Vec<f64>> = Vec::with_capacity(n + 1);
    simplex.push(initial.to_vec());
    for i in 0..n {
        let mut p = initial.to_vec();
        p[i] += step;
        simplex.push(p);
    }
    let mut values: Vec<f64> = simplex.iter().map(|x| objective(x)).collect();

    for _ in 0..max_iter {
        // Sort by function value
        let mut idx: Vec<usize> = (0..=n).collect();
        idx.sort_by(|&a, &b| {
            values[a]
                .partial_cmp(&values[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let new_simplex: Vec<Vec<f64>> = idx.iter().map(|&i| simplex[i].clone()).collect();
        let new_values: Vec<f64> = idx.iter().map(|&i| values[i]).collect();
        simplex = new_simplex;
        values = new_values;

        // Convergence check
        if values[n] - values[0] < tol {
            break;
        }

        // Centroid of best n points
        let centroid: Vec<f64> = (0..n)
            .map(|j| simplex[..n].iter().map(|p| p[j]).sum::<f64>() / n as f64)
            .collect();

        // Reflection
        let xr: Vec<f64> = (0..n)
            .map(|j| centroid[j] + alpha * (centroid[j] - simplex[n][j]))
            .collect();
        let fr = objective(&xr);

        if fr < values[0] {
            // Expansion
            let xe: Vec<f64> = (0..n)
                .map(|j| centroid[j] + gamma * (xr[j] - centroid[j]))
                .collect();
            let fe = objective(&xe);
            if fe < fr {
                simplex[n] = xe;
                values[n] = fe;
            } else {
                simplex[n] = xr;
                values[n] = fr;
            }
        } else if fr < values[n - 1] {
            simplex[n] = xr;
            values[n] = fr;
        } else {
            // Contraction
            let (base, fb) = if fr < values[n] {
                (xr.clone(), fr)
            } else {
                (simplex[n].clone(), values[n])
            };
            let xc: Vec<f64> = (0..n)
                .map(|j| centroid[j] + rho * (base[j] - centroid[j]))
                .collect();
            let fc = objective(&xc);

            if fc < fb {
                simplex[n] = xc;
                values[n] = fc;
            } else {
                // Shrink
                for i in 1..=n {
                    for j in 0..n {
                        simplex[i][j] = simplex[0][j] + sigma * (simplex[i][j] - simplex[0][j]);
                    }
                    values[i] = objective(&simplex[i]);
                }
            }
        }
    }

    (simplex[0].clone(), values[0])
}

// ============================================================================
// N-point weighted measure optimizer
// ============================================================================

/// Decode parameter vector [p₁..pₙ, θ₁..θₙ] into positions and softmax weights
fn decode_npoint(params: &[f64], n_points: usize) -> (Vec<f64>, Vec<f64>) {
    let positions: Vec<f64> = params[..n_points]
        .iter()
        .map(|&p| p.clamp(-1.0, 1.0))
        .collect();

    let raw = &params[n_points..];
    let max_w = raw.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let exp_w: Vec<f64> = raw.iter().map(|&w| (w - max_w).exp()).collect();
    let sum_w: f64 = exp_w.iter().sum();
    let weights: Vec<f64> = exp_w.iter().map(|&w| w / sum_w).collect();

    (positions, weights)
}

fn optimize_npoint_core(
    n_points: usize,
    pop_size: usize,
    max_iter: usize,
    n_restarts: usize,
    seed: u64,
) -> (Vec<f64>, Vec<f64>, f64) {
    let dim = 2 * n_points;
    let mut bounds = Vec::with_capacity(dim);
    for _ in 0..n_points {
        bounds.push((-1.0, 1.0)); // positions
    }
    for _ in 0..n_points {
        bounds.push((-5.0, 5.0)); // raw weights → softmax
    }

    let objective = move |params: &[f64]| -> f64 {
        let (pos, wt) = decode_npoint(params, n_points);
        let m = sublevel_measure_core(&pos, &wt);
        if m.is_nan() || m.is_infinite() || m <= 0.0 {
            100.0
        } else {
            m
        }
    };

    // Multi-restart DE + NM polishing, parallelized across restarts
    let results: Vec<_> = (0..n_restarts)
        .into_par_iter()
        .map(|i| {
            let de_seed = seed.wrapping_add(i as u64 * 997);
            let (de_best, de_fit) =
                differential_evolution(&objective, &bounds, pop_size, max_iter, 0.8, 0.9, de_seed);
            let (nm_best, nm_fit) = nelder_mead(&objective, &de_best, 0.01, 1000, 1e-14);
            if nm_fit < de_fit {
                (nm_best, nm_fit)
            } else {
                (de_best, de_fit)
            }
        })
        .collect();

    let (best_params, best_fit) = results
        .into_iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap();

    let (positions, weights) = decode_npoint(&best_params, n_points);
    (positions, weights, best_fit)
}

// ============================================================================
// Discrete root optimizer (equal weights, free positions)
// ============================================================================

fn optimize_discrete_core(
    n_roots: usize,
    pop_size: usize,
    max_iter: usize,
    n_restarts: usize,
    seed: u64,
) -> (Vec<f64>, f64) {
    let bounds: Vec<(f64, f64)> = vec![(-1.0, 1.0); n_roots];
    let w = 1.0 / n_roots as f64;

    let objective = move |params: &[f64]| -> f64 {
        let positions: Vec<f64> = params.iter().map(|&p| p.clamp(-1.0, 1.0)).collect();
        let weights: Vec<f64> = vec![w; n_roots];
        let m = sublevel_measure_core(&positions, &weights);
        if m.is_nan() || m.is_infinite() || m <= 0.0 {
            100.0
        } else {
            m
        }
    };

    let results: Vec<_> = (0..n_restarts)
        .into_par_iter()
        .map(|i| {
            let de_seed = seed.wrapping_add(i as u64 * 997);
            let (de_best, de_fit) =
                differential_evolution(&objective, &bounds, pop_size, max_iter, 0.8, 0.9, de_seed);
            let (nm_best, nm_fit) = nelder_mead(&objective, &de_best, 0.01, 1000, 1e-14);
            if nm_fit < de_fit {
                (nm_best, nm_fit)
            } else {
                (de_best, de_fit)
            }
        })
        .collect();

    results
        .into_iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap()
}

// ============================================================================
// Binary star sweep
// ============================================================================

fn sweep_binary_star_core(n_max: usize) -> Vec<(usize, usize, f64, f64)> {
    // Compute best k for each n in parallel
    let per_n: Vec<(usize, usize, f64)> = (2..=n_max)
        .into_par_iter()
        .map(|n| {
            let mut best_k = 1usize;
            let mut best_meas = f64::INFINITY;
            for k in 1..n {
                let alpha = k as f64 / n as f64;
                let m = sublevel_measure_core(&[-1.0, 1.0], &[alpha, 1.0 - alpha]);
                if m > 0.01 && m < best_meas {
                    best_meas = m;
                    best_k = k;
                }
            }
            (n, best_k, best_meas)
        })
        .collect();

    // Running minimum (sequential — depends on previous values)
    let mut running_min = f64::INFINITY;
    per_n
        .into_iter()
        .map(|(n, k, meas)| {
            if meas < running_min {
                running_min = meas;
            }
            (n, k, meas, running_min)
        })
        .collect()
}

// ============================================================================
// PyO3 bindings
// ============================================================================

/// Compute sublevel measure for a weighted root configuration.
///
/// Args:
///     positions: root positions [p₁, ..., pₙ]
///     weights: corresponding weights [w₁, ..., wₙ] (should sum to 1)
///
/// Returns: |{x : Σ wᵢ·ln|x-pᵢ| < 0}|
#[pyfunction]
fn sublevel_measure(positions: Vec<f64>, weights: Vec<f64>) -> PyResult<f64> {
    if positions.len() != weights.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "positions and weights must have same length",
        ));
    }
    Ok(sublevel_measure_core(&positions, &weights))
}

/// Compute sublevel measure for a discrete root polynomial.
///
/// Args:
///     roots: all roots [r₁, ..., rₙ] (may repeat)
///
/// Returns: |{x : |∏(x-rᵢ)| < 1}|
#[pyfunction]
fn sublevel_measure_discrete(roots: Vec<f64>) -> PyResult<f64> {
    let n = roots.len();
    if n == 0 {
        return Ok(0.0);
    }
    let w = 1.0 / n as f64;
    let weights = vec![w; n];
    Ok(sublevel_measure_core(&roots, &weights))
}

/// Batch evaluation of sublevel measures (parallelized with rayon).
///
/// Args:
///     configs: list of (positions, weights) tuples
///
/// Returns: list of measures
#[pyfunction]
fn sublevel_measure_batch(configs: Vec<(Vec<f64>, Vec<f64>)>) -> PyResult<Vec<f64>> {
    Ok(configs
        .par_iter()
        .map(|(pos, wt)| sublevel_measure_core(pos, wt))
        .collect())
}

/// Binary star sublevel measure for given mass ratio α.
///
/// Computes measure for f(x) = (x+1)^{αn}(x-1)^{(1-α)n} in the continuous limit.
#[pyfunction]
fn binary_star_measure(alpha: f64) -> PyResult<f64> {
    if alpha <= 0.0 || alpha >= 1.0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "alpha must be in (0, 1)",
        ));
    }
    Ok(sublevel_measure_core(&[-1.0, 1.0], &[alpha, 1.0 - alpha]))
}

/// Sweep binary star family for n = 2..n_max, find best k for each n.
///
/// Returns: list of (n, best_k, measure, running_min) tuples
#[pyfunction]
#[pyo3(signature = (n_max,))]
fn sweep_binary_star(n_max: usize) -> Vec<(usize, usize, f64, f64)> {
    sweep_binary_star_core(n_max)
}

/// Optimize N-point weighted measure via multi-restart DE + Nelder-Mead.
///
/// Finds positions p₁,...,pₙ ∈ [-1,1] and weights w₁,...,wₙ (softmax)
/// that minimize the sublevel measure.
///
/// Returns: (positions, weights, measure, n_effective_points)
#[pyfunction]
#[pyo3(signature = (n_points, pop_size=100, max_iter=2000, n_restarts=50, seed=42))]
fn optimize_npoint(
    n_points: usize,
    pop_size: usize,
    max_iter: usize,
    n_restarts: usize,
    seed: u64,
) -> (Vec<f64>, Vec<f64>, f64, usize) {
    let (positions, weights, measure) =
        optimize_npoint_core(n_points, pop_size, max_iter, n_restarts, seed);
    let n_effective = weights.iter().filter(|&&w| w > 0.01).count();
    (positions, weights, measure, n_effective)
}

/// Optimize discrete root placement via multi-restart DE + Nelder-Mead.
///
/// Finds n root positions in [-1,1] (equal weights) minimizing sublevel measure.
///
/// Returns: (roots, measure)
#[pyfunction]
#[pyo3(signature = (n_roots, pop_size=100, max_iter=2000, n_restarts=50, seed=42))]
fn optimize_discrete(
    n_roots: usize,
    pop_size: usize,
    max_iter: usize,
    n_restarts: usize,
    seed: u64,
) -> (Vec<f64>, f64) {
    optimize_discrete_core(n_roots, pop_size, max_iter, n_restarts, seed)
}

/// Extract interval-arithmetic certified panel constraints for the
/// slack-augmented exchange LP on Erdős #1038.
///
/// Algorithm:
///   1. Locate the sublevel-set external zeros (z_left, z_right) of the
///      logarithmic potential g(x) = Σ wᵢ ln|x - pᵢ|.
///   2. Place `panel_count` panels on the legal-support side, namely on
///      (z_left - 1.0, z_left) and (z_right, z_right + 1.0).
///   3. Excise the `root_epsilon` neighborhood of every root.
///   4. For each (panel, root) pair, compute `(panel - r).abs().ln()` in
///      inari IEEE-1788 interval arithmetic and store the LOWER endpoint
///      as `A_ineq[p][i]`.  Panels touching any root are filtered out.
///   5. `baseline_slack[p] = Σ_i wᵢ · A_ineq[p][i]`. Active = slack ≤ 1e-10.
///
/// Args:
///     positions: root positions
///     weights: corresponding weights (caller-supplied; alignment must match positions)
///     panel_count: total panels to place across both legal-support sides
///     root_epsilon: half-width of the excised neighborhood around each root
///
/// Returns: (panels, A_ineq, baseline_slack, n_active, n_inactive)
#[pyfunction]
fn extract_panel_constraints(
    positions: Vec<f64>,
    weights: Vec<f64>,
    panel_count: usize,
    root_epsilon: f64,
) -> PyResult<(Vec<(f64, f64)>, Vec<Vec<f64>>, Vec<f64>, usize, usize)> {
    if positions.len() != weights.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "positions and weights must have same length",
        ));
    }
    if positions.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "positions must be non-empty",
        ));
    }
    Ok(extract_panel_constraints_core(
        &positions,
        &weights,
        panel_count,
        root_epsilon,
    ))
}

#[pymodule]
fn erdos1038_fast(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sublevel_measure, m)?)?;
    m.add_function(wrap_pyfunction!(sublevel_measure_discrete, m)?)?;
    m.add_function(wrap_pyfunction!(sublevel_measure_batch, m)?)?;
    m.add_function(wrap_pyfunction!(binary_star_measure, m)?)?;
    m.add_function(wrap_pyfunction!(sweep_binary_star, m)?)?;
    m.add_function(wrap_pyfunction!(optimize_npoint, m)?)?;
    m.add_function(wrap_pyfunction!(optimize_discrete, m)?)?;
    m.add_function(wrap_pyfunction!(extract_panel_constraints, m)?)?;
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_root() {
        // Single root at 0: sublevel = (-1, 1), measure = 2
        let m = sublevel_measure_core(&[0.0], &[1.0]);
        assert!((m - 2.0).abs() < 1e-10, "single root: got {m}");
    }

    #[test]
    fn test_binary_star_known() {
        // Binary star at α ≈ 0.14595 should give ≈ 1.8716
        let m = sublevel_measure_core(&[-1.0, 1.0], &[0.14595, 0.85405]);
        assert!(
            (m - 1.8716).abs() < 0.001,
            "binary star: got {m}, expected ~1.8716"
        );
    }

    #[test]
    fn test_symmetric() {
        // Symmetric: α = 0.5 at ±1 → should give 2√2 ≈ 2.8284
        let m = sublevel_measure_core(&[-1.0, 1.0], &[0.5, 0.5]);
        assert!(
            (m - 2.8284).abs() < 0.01,
            "symmetric: got {m}, expected ~2.8284"
        );
    }

    #[test]
    fn test_three_point() {
        // Three-point should be less than binary star
        let binary = sublevel_measure_core(&[-1.0, 1.0], &[0.14595, 0.85405]);
        let triple = sublevel_measure_core(&[-1.0, 1.0, 0.881], &[0.849, 0.105, 0.046]);
        assert!(
            triple < binary,
            "three-point ({triple}) should beat binary ({binary})"
        );
    }

    #[test]
    fn test_sweep_monotone() {
        let results = sweep_binary_star_core(100);
        for i in 1..results.len() {
            assert!(
                results[i].3 <= results[i - 1].3 + 1e-14,
                "running min not monotone at n={}",
                results[i].0
            );
        }
    }

    #[test]
    fn test_extract_panel_constraints_n2() {
        // Symmetric N=2 case: roots at ±1 with equal weights.
        let positions = vec![-1.0_f64, 1.0_f64];
        let weights = vec![0.5_f64, 0.5_f64];
        let (panels, a_ineq, baseline, n_act, n_inact) =
            extract_panel_constraints_core(&positions, &weights, 16, 1e-3);

        assert!(
            !panels.is_empty(),
            "expected at least 1 panel for the N=2 case"
        );
        assert_eq!(
            a_ineq.len(),
            panels.len(),
            "A_ineq row count must match panel count"
        );
        for row in &a_ineq {
            assert_eq!(
                row.len(),
                positions.len(),
                "A_ineq column count must equal len(positions)"
            );
            for &v in row {
                assert!(v.is_finite(), "all A_ineq entries must be finite");
            }
        }
        assert_eq!(baseline.len(), panels.len());
        assert_eq!(n_act + n_inact, panels.len());
        // On the legal-support side, baseline slack should be non-negative
        // (the constraint g(x) ≥ 0 holds where panels live).
        for (i, &s) in baseline.iter().enumerate() {
            assert!(
                s >= -1e-9,
                "panel {i} has unexpectedly negative baseline slack: {s}"
            );
        }
    }
}
