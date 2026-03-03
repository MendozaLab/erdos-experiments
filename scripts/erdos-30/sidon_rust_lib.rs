use pyo3::prelude::*;
use std::collections::HashSet;

/// Greedy Sidon set construction: add elements 1..N in order,
/// skip if adding x would create a duplicate pairwise sum.
/// Returns (sidon_set, size).
#[pyfunction]
fn greedy_sidon(n: u64) -> (Vec<u64>, usize) {
    let mut sidon: Vec<u64> = Vec::with_capacity((n as f64).sqrt() as usize + 100);
    let mut sums: HashSet<u64> = HashSet::with_capacity(n as usize);

    for x in 1..=n {
        let mut valid = true;
        let mut new_sums: Vec<u64> = Vec::with_capacity(sidon.len() + 1);

        for &a in &sidon {
            let s = x + a;
            if sums.contains(&s) {
                valid = false;
                break;
            }
            new_sums.push(s);
        }

        if valid {
            let s_self = 2 * x;
            if sums.contains(&s_self) {
                continue;
            }
            new_sums.push(s_self);
            sidon.push(x);
            for s in new_sums {
                sums.insert(s);
            }
        }
    }

    let size = sidon.len();
    (sidon, size)
}

/// Reverse greedy: add elements N..1 in descending order.
#[pyfunction]
fn greedy_sidon_reverse(n: u64) -> (Vec<u64>, usize) {
    let mut sidon: Vec<u64> = Vec::with_capacity((n as f64).sqrt() as usize + 100);
    let mut sums: HashSet<u64> = HashSet::with_capacity(n as usize);

    for x in (1..=n).rev() {
        let mut valid = true;
        let mut new_sums: Vec<u64> = Vec::with_capacity(sidon.len() + 1);

        for &a in &sidon {
            let s = x + a;
            if sums.contains(&s) {
                valid = false;
                break;
            }
            new_sums.push(s);
        }

        if valid {
            let s_self = 2 * x;
            if sums.contains(&s_self) {
                continue;
            }
            new_sums.push(s_self);
            sidon.push(x);
            for s in new_sums {
                sums.insert(s);
            }
        }
    }

    sidon.sort();
    let size = sidon.len();
    (sidon, size)
}

/// Greedy B₃ set construction: all triple sums a+b+c (a≤b≤c) distinct.
/// Returns (b3_set, size).
#[pyfunction]
fn greedy_b3(n: u64) -> (Vec<u64>, usize) {
    let mut b3: Vec<u64> = Vec::with_capacity((n as f64).cbrt() as usize + 50);
    let mut triple_sums: HashSet<u64> = HashSet::with_capacity(n as usize);

    for x in 1..=n {
        let mut valid = true;
        let mut new_sums: Vec<u64> = Vec::new();

        // (x, x, x)
        let s_xxx = 3 * x;
        if triple_sums.contains(&s_xxx) {
            continue;
        }
        new_sums.push(s_xxx);

        for (i, &a) in b3.iter().enumerate() {
            // (x, x, a) and (x, a, a)
            let s_xxa = 2 * x + a;
            let s_xaa = x + 2 * a;

            if triple_sums.contains(&s_xxa) || triple_sums.contains(&s_xaa) {
                valid = false;
                break;
            }
            // Check against new_sums we're accumulating
            if new_sums.contains(&s_xxa) || new_sums.contains(&s_xaa) {
                valid = false;
                break;
            }
            new_sums.push(s_xxa);
            new_sums.push(s_xaa);

            // (x, a, b) for b > a
            for &b in &b3[i + 1..] {
                let s_xab = x + a + b;
                if triple_sums.contains(&s_xab) || new_sums.contains(&s_xab) {
                    valid = false;
                    break;
                }
                new_sums.push(s_xab);
            }
            if !valid {
                break;
            }
        }

        if valid {
            b3.push(x);
            for s in new_sums {
                triple_sums.insert(s);
            }
        }
    }

    let size = b3.len();
    (b3, size)
}

/// WalkSAT for k-AP-free 2-colorings.
/// Returns (found, best_violations, coloring_as_bytes, flips_used).
#[pyfunction]
fn walksat_kap(n: usize, k: usize, max_restarts: usize, max_flips: usize) -> (bool, usize, Vec<u8>, usize) {
    use std::collections::HashSet;

    // Precompute all k-APs in {0..n-1}
    let mut aps: Vec<Vec<usize>> = Vec::new();
    let max_d = (n - 1) / (k - 1);
    for d in 1..=max_d {
        let max_a = n - (k - 1) * d;
        for a in 0..max_a {
            let ap: Vec<usize> = (0..k).map(|j| a + j * d).collect();
            aps.push(ap);
        }
    }

    let mut best_violations = n * n;
    let mut best_coloring: Vec<u8> = vec![0; n];
    let mut rng_state: u64 = 0xdeadbeef_u64;

    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    for restart in 0..max_restarts {
        // Random coloring
        let mut coloring: Vec<u8> = (0..n).map(|_| (xorshift(&mut rng_state) & 1) as u8).collect();

        for flip in 0..max_flips {
            // Count violations
            let mut violations: Vec<usize> = Vec::new();
            for (i, ap) in aps.iter().enumerate() {
                let c0 = coloring[ap[0]];
                if ap.iter().all(|&idx| coloring[idx] == c0) {
                    violations.push(i);
                }
            }

            let n_viol = violations.len();
            if n_viol == 0 {
                return (true, 0, coloring, restart * max_flips + flip);
            }

            if n_viol < best_violations {
                best_violations = n_viol;
                best_coloring = coloring.clone();
            }

            // Pick random violated AP, flip random element
            let viol_idx = violations[(xorshift(&mut rng_state) as usize) % violations.len()];
            let ap = &aps[viol_idx];
            let pos = ap[(xorshift(&mut rng_state) as usize) % k];
            coloring[pos] = 1 - coloring[pos];
        }

        // Seed variation per restart
        rng_state = rng_state.wrapping_add(restart as u64 * 0x9e3779b97f4a7c15);
    }

    (false, best_violations, best_coloring, max_restarts * max_flips)
}

#[pymodule]
fn sidon_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(greedy_sidon, m)?)?;
    m.add_function(wrap_pyfunction!(greedy_sidon_reverse, m)?)?;
    m.add_function(wrap_pyfunction!(greedy_b3, m)?)?;
    m.add_function(wrap_pyfunction!(walksat_kap, m)?)?;
    Ok(())
}
