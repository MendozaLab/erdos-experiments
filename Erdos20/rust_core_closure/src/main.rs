use std::env;
use std::fs;
use std::time::Instant;

type SiteMask = u64;
type FamilyMask = u128;

#[derive(Clone, Default)]
struct Acc {
    channel_count: u128,
    candidate_count: u128,
    local_valid_count: u128,
    global_valid_count: u128,
}

impl Acc {
    fn add(&mut self, candidates: u32, local_valid: u32, global_valid: u32) {
        if candidates == 0 {
            return;
        }
        self.channel_count += 1;
        self.candidate_count += candidates as u128;
        self.local_valid_count += local_valid as u128;
        self.global_valid_count += global_valid as u128;
    }

    fn i_local(&self) -> Option<f64> {
        neg_log2_ratio(self.local_valid_count, self.candidate_count)
    }

    fn i_global(&self) -> Option<f64> {
        neg_log2_ratio(self.global_valid_count, self.candidate_count)
    }
}

struct Config {
    n: usize,
    w: usize,
    target_m: Option<usize>,
    window: usize,
    emit_safe_candidates_v1: bool,
    out_path: Option<String>,
}

struct Enumerator {
    w: usize,
    sites: Vec<SiteMask>,
    cores: Vec<SiteMask>,
    core_sizes: Vec<usize>,
    core_candidate_masks: Vec<FamilyMask>,
    pair_core: Vec<Vec<usize>>,
    pair_block: Vec<Vec<FamilyMask>>,
    all_site_mask: FamilyMask,
}

fn neg_log2_ratio(num: u128, den: u128) -> Option<f64> {
    if num == 0 || den == 0 {
        None
    } else {
        Some(-((num as f64) / (den as f64)).log2())
    }
}

fn bit_count(mask: FamilyMask) -> u32 {
    mask.count_ones()
}

fn enumerate_w_subsets(n: usize, w: usize) -> Vec<SiteMask> {
    fn rec(out: &mut Vec<SiteMask>, n: usize, w: usize, start: usize, left: usize, mask: SiteMask) {
        if left == 0 {
            out.push(mask);
            return;
        }
        for idx in start..=n - left {
            rec(out, n, w, idx + 1, left - 1, mask | (1u64 << idx));
        }
        let _ = w;
    }
    let mut out = Vec::new();
    rec(&mut out, n, w, 0, w, 0);
    out
}

fn enumerate_subsets_upto(n: usize, max_size_exclusive: usize) -> Vec<(SiteMask, usize)> {
    fn rec(
        out: &mut Vec<(SiteMask, usize)>,
        n: usize,
        target: usize,
        start: usize,
        left: usize,
        mask: SiteMask,
    ) {
        if left == 0 {
            out.push((mask, target));
            return;
        }
        for idx in start..=n - left {
            rec(out, n, target, idx + 1, left - 1, mask | (1u64 << idx));
        }
    }
    let mut out = Vec::new();
    for size in 0..max_size_exclusive {
        rec(&mut out, n, size, 0, size, 0);
    }
    out
}

impl Enumerator {
    fn new(n: usize, w: usize) -> Result<Self, String> {
        if n > 64 {
            return Err("n > 64 is not supported by u64 set masks".to_string());
        }
        let sites = enumerate_w_subsets(n, w);
        if sites.len() > 128 {
            return Err("more than 128 sites is not supported by u128 family masks".to_string());
        }
        let all_site_mask = if sites.len() == 128 {
            u128::MAX
        } else {
            (1u128 << sites.len()) - 1
        };

        let core_pairs = enumerate_subsets_upto(n, w);
        let mut cores = Vec::new();
        let mut core_sizes = Vec::new();
        for (mask, size) in core_pairs {
            cores.push(mask);
            core_sizes.push(size);
        }
        let mut core_candidate_masks = Vec::new();
        for core in &cores {
            let mut candidates = 0u128;
            for (idx, site) in sites.iter().enumerate() {
                if (site & core) == *core {
                    candidates |= 1u128 << idx;
                }
            }
            core_candidate_masks.push(candidates);
        }

        let m = sites.len();
        let mut pair_core = vec![vec![0usize; m]; m];
        let mut pair_block = vec![vec![0u128; m]; m];
        for left in 0..m {
            for right in (left + 1)..m {
                let core = sites[left] & sites[right];
                let core_idx = cores
                    .iter()
                    .position(|candidate| *candidate == core)
                    .ok_or_else(|| "core lookup failed".to_string())?;
                let mut block = 0u128;
                for (cand_idx, cand) in sites.iter().enumerate() {
                    if cand_idx == left || cand_idx == right {
                        continue;
                    }
                    if (cand & core) == core
                        && (cand & sites[left]) == core
                        && (cand & sites[right]) == core
                    {
                        block |= 1u128 << cand_idx;
                    }
                }
                pair_core[left][right] = core_idx;
                pair_core[right][left] = core_idx;
                pair_block[left][right] = block;
                pair_block[right][left] = block;
            }
        }

        Ok(Self {
            w,
            sites,
            cores,
            core_sizes,
            core_candidate_masks,
            pair_core,
            pair_block,
            all_site_mask,
        })
    }

    fn density_pass(&self) -> (Vec<u128>, f64) {
        let mut density = vec![0u128; self.sites.len() + 1];
        let mut family_indices: Vec<usize> = Vec::new();
        let start = Instant::now();
        self.backtrack_density(0, 0, 0, &mut family_indices, &mut density);
        (trim_density(density), start.elapsed().as_secs_f64())
    }

    fn backtrack_density(
        &self,
        start_idx: usize,
        family_mask: FamilyMask,
        blocked_global: FamilyMask,
        family_indices: &mut Vec<usize>,
        density: &mut Vec<u128>,
    ) {
        density[family_indices.len()] += 1;
        for site_idx in start_idx..self.sites.len() {
            let site_bit = 1u128 << site_idx;
            if (family_mask & site_bit) != 0 || (blocked_global & site_bit) != 0 {
                continue;
            }
            let mut new_blocked = blocked_global;
            for existing in family_indices.iter() {
                new_blocked |= self.pair_block[*existing][site_idx];
            }
            family_indices.push(site_idx);
            self.backtrack_density(site_idx + 1, family_mask | site_bit, new_blocked, family_indices, density);
            family_indices.pop();
        }
    }

    fn core_pass(&self, record_ms: &[usize]) -> (Vec<Vec<Acc>>, f64) {
        let mut accum = vec![vec![Acc::default(); self.w]; self.sites.len() + 1];
        let mut family_indices: Vec<usize> = Vec::new();
        let mut blocked_by_core = vec![0u128; self.cores.len()];
        let start = Instant::now();
        self.backtrack_core(
            0,
            0,
            0,
            &mut blocked_by_core,
            &mut family_indices,
            record_ms,
            &mut accum,
        );
        (accum, start.elapsed().as_secs_f64())
    }

    fn backtrack_core(
        &self,
        start_idx: usize,
        family_mask: FamilyMask,
        blocked_global: FamilyMask,
        blocked_by_core: &mut [FamilyMask],
        family_indices: &mut Vec<usize>,
        record_ms: &[usize],
        accum: &mut [Vec<Acc>],
    ) {
        let m = family_indices.len();
        if record_ms.contains(&m) {
            self.record_core_channels(m, family_mask, blocked_global, blocked_by_core, accum);
        }

        for site_idx in start_idx..self.sites.len() {
            let site_bit = 1u128 << site_idx;
            if (family_mask & site_bit) != 0 || (blocked_global & site_bit) != 0 {
                continue;
            }
            let mut new_blocked = blocked_global;
            let mut changed: Vec<(usize, FamilyMask)> = Vec::new();
            for existing in family_indices.iter() {
                let core_idx = self.pair_core[*existing][site_idx];
                let block = self.pair_block[*existing][site_idx];
                let old = blocked_by_core[core_idx];
                let new = old | block;
                if new != old {
                    changed.push((core_idx, old));
                    blocked_by_core[core_idx] = new;
                }
                new_blocked |= block;
            }
            family_indices.push(site_idx);
            self.backtrack_core(
                site_idx + 1,
                family_mask | site_bit,
                new_blocked,
                blocked_by_core,
                family_indices,
                record_ms,
                accum,
            );
            family_indices.pop();
            for (core_idx, old) in changed.into_iter().rev() {
                blocked_by_core[core_idx] = old;
            }
        }
    }

    fn record_core_channels(
        &self,
        m: usize,
        family_mask: FamilyMask,
        blocked_global: FamilyMask,
        blocked_by_core: &[FamilyMask],
        accum: &mut [Vec<Acc>],
    ) {
        let unused = self.all_site_mask & !family_mask;
        let global_valid = unused & !blocked_global;
        for (core_idx, candidates_all) in self.core_candidate_masks.iter().enumerate() {
            let candidates_mask = candidates_all & unused;
            let candidates = bit_count(candidates_mask);
            if candidates == 0 {
                continue;
            }
            let local_valid = bit_count(candidates_mask & !blocked_by_core[core_idx]);
            let global_valid_count = bit_count(candidates_mask & global_valid);
            let s = self.core_sizes[core_idx];
            accum[m][s].add(candidates, local_valid, global_valid_count);
        }
    }

    /// Find the first (lexicographically earliest, by site_idx ascending) sunflower-free
    /// family of size `target_m`. Returns the family as a Vec of site indices.
    fn first_sf_free_family_of_size(&self, target_m: usize) -> Option<Vec<usize>> {
        let mut family_indices: Vec<usize> = Vec::new();
        let mut found: Option<Vec<usize>> = None;
        self.find_first_family(0, 0, 0, target_m, &mut family_indices, &mut found);
        found
    }

    fn find_first_family(
        &self,
        start_idx: usize,
        family_mask: FamilyMask,
        blocked_global: FamilyMask,
        target_m: usize,
        family_indices: &mut Vec<usize>,
        found: &mut Option<Vec<usize>>,
    ) {
        if found.is_some() {
            return;
        }
        if family_indices.len() == target_m {
            *found = Some(family_indices.clone());
            return;
        }
        for site_idx in start_idx..self.sites.len() {
            if found.is_some() {
                return;
            }
            let site_bit = 1u128 << site_idx;
            if (family_mask & site_bit) != 0 || (blocked_global & site_bit) != 0 {
                continue;
            }
            let mut new_blocked = blocked_global;
            for existing in family_indices.iter() {
                new_blocked |= self.pair_block[*existing][site_idx];
            }
            family_indices.push(site_idx);
            self.find_first_family(
                site_idx + 1,
                family_mask | site_bit,
                new_blocked,
                target_m,
                family_indices,
                found,
            );
            family_indices.pop();
        }
    }
}

fn trim_density(mut density: Vec<u128>) -> Vec<u128> {
    while density.last() == Some(&0) {
        density.pop();
    }
    density
}

fn jamming_m(density: &[u128]) -> Option<usize> {
    for m in 0..density.len().saturating_sub(1) {
        if density[m] > 0 && density[m + 1] < density[m] {
            return Some(m);
        }
    }
    None
}

fn fmt_opt(value: Option<f64>) -> String {
    match value {
        Some(v) if v.is_finite() => format!("{:.6}", v),
        _ => "null".to_string(),
    }
}

fn json_array_u128(values: &[u128]) -> String {
    let mut out = String::from("[");
    for (idx, value) in values.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str(&value.to_string());
    }
    out.push(']');
    out
}

/// Convert a u64 site mask to a Vec of 1-based element indices (sorted ascending).
fn mask_to_elements_1based(mask: SiteMask, n: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for idx in 0..n {
        if (mask & (1u64 << idx)) != 0 {
            out.push(idx + 1);
        }
    }
    out
}

fn json_array_usize(values: &[usize]) -> String {
    let mut out = String::from("[");
    for (idx, v) in values.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push_str(&v.to_string());
    }
    out.push(']');
    out
}

/// Check whether candidate B (site index cand_idx) is locally safe through core C
/// against family F (site indices). For k=3, we need: there is no pair A1, A2 in F,
/// distinct, with C ⊆ A1, C ⊆ A2, A1 ∩ A2 = C, A1 ∩ B = C, A2 ∩ B = C.
/// Returns (is_safe, list of (a1_idx, a2_idx) witness pairs).
fn check_local_safety(
    enumerator: &Enumerator,
    family_indices: &[usize],
    core_mask: SiteMask,
    cand_idx: usize,
) -> (bool, Vec<(usize, usize)>) {
    let cand_set = enumerator.sites[cand_idx];
    let mut witnesses: Vec<(usize, usize)> = Vec::new();
    // Filter family members that contain core C and intersect B exactly in C.
    let mut eligible: Vec<usize> = Vec::new();
    for &a_idx in family_indices.iter() {
        let a_set = enumerator.sites[a_idx];
        if (a_set & core_mask) != core_mask {
            continue;
        }
        if (a_set & cand_set) != core_mask {
            continue;
        }
        eligible.push(a_idx);
    }
    // Now scan all distinct pairs.
    for i in 0..eligible.len() {
        for j in (i + 1)..eligible.len() {
            let a1_idx = eligible[i];
            let a2_idx = eligible[j];
            let a1 = enumerator.sites[a1_idx];
            let a2 = enumerator.sites[a2_idx];
            if (a1 & a2) != core_mask {
                continue;
            }
            // All constraints satisfied: this is a witness pair that closes a 3-sunflower with B.
            witnesses.push((a1_idx, a2_idx));
        }
    }
    (witnesses.is_empty(), witnesses)
}

/// Compute SHA-256 of a byte slice using a minimal pure-Rust implementation.
mod sha256 {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];

    pub fn hex(input: &[u8]) -> String {
        let mut h: [u32; 8] = [
            0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
            0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
        ];
        // Pad message.
        let bit_len: u64 = (input.len() as u64) * 8;
        let mut msg: Vec<u8> = input.to_vec();
        msg.push(0x80);
        while msg.len() % 64 != 56 {
            msg.push(0);
        }
        msg.extend_from_slice(&bit_len.to_be_bytes());

        for chunk in msg.chunks_exact(64) {
            let mut w = [0u32; 64];
            for i in 0..16 {
                let b = i * 4;
                w[i] = u32::from_be_bytes([chunk[b], chunk[b + 1], chunk[b + 2], chunk[b + 3]]);
            }
            for i in 16..64 {
                let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
                let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
                w[i] = w[i - 16]
                    .wrapping_add(s0)
                    .wrapping_add(w[i - 7])
                    .wrapping_add(s1);
            }
            let (mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh) =
                (h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]);
            for i in 0..64 {
                let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
                let ch = (e & f) ^ (!e & g);
                let t1 = hh
                    .wrapping_add(s1)
                    .wrapping_add(ch)
                    .wrapping_add(K[i])
                    .wrapping_add(w[i]);
                let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
                let mj = (a & b) ^ (a & c) ^ (b & c);
                let t2 = s0.wrapping_add(mj);
                hh = g;
                g = f;
                f = e;
                e = d.wrapping_add(t1);
                d = c;
                c = b;
                b = a;
                a = t1.wrapping_add(t2);
            }
            h[0] = h[0].wrapping_add(a);
            h[1] = h[1].wrapping_add(b);
            h[2] = h[2].wrapping_add(c);
            h[3] = h[3].wrapping_add(d);
            h[4] = h[4].wrapping_add(e);
            h[5] = h[5].wrapping_add(f);
            h[6] = h[6].wrapping_add(g);
            h[7] = h[7].wrapping_add(hh);
        }
        let mut s = String::with_capacity(64);
        for v in h.iter() {
            s.push_str(&format!("{:08x}", v));
        }
        s
    }
}

fn parse_args() -> Result<Config, String> {
    let args: Vec<String> = env::args().collect();
    let mut n = None;
    let mut w = None;
    let mut target_m = None;
    let mut window = 0usize;
    let mut emit_safe_candidates_v1 = false;
    let mut out_path: Option<String> = None;
    let mut idx = 1;
    while idx < args.len() {
        match args[idx].as_str() {
            "--n" => {
                idx += 1;
                n = Some(args.get(idx).ok_or("--n needs a value")?.parse::<usize>().map_err(|e| e.to_string())?);
            }
            "--w" => {
                idx += 1;
                w = Some(args.get(idx).ok_or("--w needs a value")?.parse::<usize>().map_err(|e| e.to_string())?);
            }
            "--target-m" => {
                idx += 1;
                let raw = args.get(idx).ok_or("--target-m needs a value")?;
                if raw != "auto" {
                    target_m = Some(raw.parse::<usize>().map_err(|e| e.to_string())?);
                }
            }
            "--window" => {
                idx += 1;
                window = args.get(idx).ok_or("--window needs a value")?.parse::<usize>().map_err(|e| e.to_string())?;
            }
            "--emit-safe-candidates-v1" => {
                emit_safe_candidates_v1 = true;
            }
            "--out" => {
                idx += 1;
                out_path = Some(args.get(idx).ok_or("--out needs a value")?.to_string());
            }
            "--help" | "-h" => {
                println!("Usage: sunflower-core-closure --n N --w W [--target-m auto|M] [--window K] [--emit-safe-candidates-v1 --out PATH]");
                std::process::exit(0);
            }
            other => return Err(format!("unknown argument: {}", other)),
        }
        idx += 1;
    }
    Ok(Config {
        n: n.ok_or("--n is required")?,
        w: w.ok_or("--w is required")?,
        target_m,
        window,
        emit_safe_candidates_v1,
        out_path,
    })
}

fn iso8601_now_utc() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Convert to Y-M-D H:M:S in UTC manually.
    let (year, month, day, hour, minute, second) = unix_seconds_to_ymdhms(secs);
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        year, month, day, hour, minute, second
    )
}

fn unix_seconds_to_ymdhms(secs: u64) -> (u32, u32, u32, u32, u32, u32) {
    let days = (secs / 86400) as i64;
    let rem = (secs % 86400) as u64;
    let hour = (rem / 3600) as u32;
    let minute = ((rem % 3600) / 60) as u32;
    let second = (rem % 60) as u32;
    // Civil-from-days algorithm by Howard Hinnant.
    let z = days + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = (z - era * 146097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = (yoe as i64) + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let year = (y + if m <= 2 { 1 } else { 0 }) as u32;
    (year, m, d, hour, minute, second)
}

fn read_self_binary_sha256() -> Option<String> {
    let exe = std::env::current_exe().ok()?;
    let bytes = fs::read(&exe).ok()?;
    Some(sha256::hex(&bytes))
}

fn run_default(config: &Config, enumerator: &Enumerator, density: &[u128], density_seconds: f64) -> Result<(), String> {
    let auto_target = jamming_m(density).unwrap_or_else(|| density.len().saturating_sub(1));
    let target = config.target_m.unwrap_or(auto_target);
    let start = target.saturating_sub(config.window);
    let end = usize::min(target + config.window, density.len().saturating_sub(1));
    let record_ms: Vec<usize> = (start..=end).collect();
    let (accum, core_seconds) = enumerator.core_pass(&record_ms);

    let mut rows = String::from("[");
    let mut first = true;
    for m in record_ms.iter() {
        for s in 0..config.w {
            let acc = &accum[*m][s];
            if acc.channel_count == 0 {
                continue;
            }
            if !first {
                rows.push(',');
            }
            first = false;
            rows.push_str(&format!(
                "{{\"m\":{},\"core_size_s\":{},\"petal_size\":{},\"channel_count\":{},\"candidate_count\":{},\"local_valid_count\":{},\"global_valid_count\":{},\"I_core_local_bits\":{},\"I_core_global_bits\":{}}}",
                m,
                s,
                config.w - s,
                acc.channel_count,
                acc.candidate_count,
                acc.local_valid_count,
                acc.global_valid_count,
                fmt_opt(acc.i_local()),
                fmt_opt(acc.i_global())
            ));
        }
    }
    rows.push(']');

    let aggregate_i = if target + 1 < density.len() && density[target] > 0 {
        let growth = density[target + 1] as f64 / density[target] as f64;
        neg_log2_ratio((growth * 1_000_000_000.0) as u128, ((enumerator.sites.len() - target) as u128) * 1_000_000_000u128)
    } else {
        None
    };

    println!(
        "{{\"runner\":\"sunflower-core-closure-rust\",\"n\":{},\"w\":{},\"k\":3,\"site_count\":{},\"core_count\":{},\"density\":{},\"total_sf_free\":{},\"max_family_size\":{},\"auto_jamming_m\":{},\"target_m\":{},\"record_window\":{},\"aggregate_I_close_bits_at_target\":{},\"density_seconds\":{:.6},\"core_seconds\":{:.6},\"core_rows\":{}}}",
        config.n,
        config.w,
        enumerator.sites.len(),
        enumerator.cores.len(),
        json_array_u128(&density),
        density.iter().sum::<u128>(),
        density.len().saturating_sub(1),
        auto_target,
        target,
        config.window,
        fmt_opt(aggregate_i),
        density_seconds,
        core_seconds,
        rows
    );
    Ok(())
}

fn run_emit_safe_candidates_v1(config: &Config, enumerator: &Enumerator) -> Result<(), String> {
    // Determine target_m. For this mode it must be explicit.
    let target_m = config
        .target_m
        .ok_or("--emit-safe-candidates-v1 requires --target-m <M>".to_string())?;
    let out_path = config
        .out_path
        .clone()
        .ok_or("--emit-safe-candidates-v1 requires --out <PATH>".to_string())?;

    let n = config.n;
    let w = config.w;
    let _ = w; // used via enumerator.w below
    let k = 3usize;

    // Compute density (cheap) for jamming reference.
    let (density, _density_seconds) = enumerator.density_pass();
    let jamming_m_star = jamming_m(&density).unwrap_or_else(|| density.len().saturating_sub(1));

    // Find first sunflower-free family of size target_m.
    let family = enumerator
        .first_sf_free_family_of_size(target_m)
        .ok_or(format!(
            "no sunflower-free family of size {} exists for n={}, w={}",
            target_m, n, enumerator.w
        ))?;

    // Family canonical hash: hex of sorted set masks (u64 little-endian) concatenated.
    let mut sorted_masks: Vec<SiteMask> = family.iter().map(|i| enumerator.sites[*i]).collect();
    sorted_masks.sort();
    let mut bytes = Vec::with_capacity(sorted_masks.len() * 8);
    for m in &sorted_masks {
        bytes.extend_from_slice(&m.to_le_bytes());
    }
    let family_hash = sha256::hex(&bytes);

    // Sets block.
    let mut sets_json = String::from("[");
    for (i, &site_idx) in family.iter().enumerate() {
        if i > 0 {
            sets_json.push(',');
        }
        let mask = enumerator.sites[site_idx];
        let elements = mask_to_elements_1based(mask, n);
        // Invariants 1, 3:
        for &e in &elements {
            if e < 1 || e > n {
                return Err(format!("invariant 1 violation: element {} out of [1..{}]", e, n));
            }
        }
        if elements.len() != enumerator.w {
            return Err(format!("invariant 3 violation (uniformity): set has size {}", elements.len()));
        }
        sets_json.push_str(&format!(
            "{{\"set_id\":\"S-{:04}\",\"elements\":{},\"mask_u64\":{}}}",
            i + 1,
            json_array_usize(&elements),
            mask
        ));
    }
    sets_json.push(']');

    // Invariant 4: family-size.
    if family.len() != target_m {
        return Err(format!(
            "invariant 4 violation: m={} but family has {} sets",
            target_m,
            family.len()
        ));
    }
    // Invariant 5: family uniqueness.
    {
        let mut seen = Vec::new();
        for &i in family.iter() {
            let m = enumerator.sites[i];
            if seen.contains(&m) {
                return Err("invariant 5 violation: duplicate family set mask".to_string());
            }
            seen.push(m);
        }
    }

    // Build sunflower_free_check via exhaustive triple scan.
    let mut checked_triples = 0u64;
    let mut sf_witness: Option<(usize, usize, usize, SiteMask)> = None;
    if family.len() >= k {
        for a in 0..family.len() {
            for b in (a + 1)..family.len() {
                for c in (b + 1)..family.len() {
                    checked_triples += 1;
                    let s_a = enumerator.sites[family[a]];
                    let s_b = enumerator.sites[family[b]];
                    let s_c = enumerator.sites[family[c]];
                    let core = s_a & s_b & s_c;
                    let petal_a = s_a & !core;
                    let petal_b = s_b & !core;
                    let petal_c = s_c & !core;
                    if (petal_a & petal_b) == 0
                        && (petal_a & petal_c) == 0
                        && (petal_b & petal_c) == 0
                        && (s_a & s_b) == core
                        && (s_a & s_c) == core
                        && (s_b & s_c) == core
                    {
                        sf_witness = Some((a, b, c, core));
                        break;
                    }
                }
                if sf_witness.is_some() {
                    break;
                }
            }
            if sf_witness.is_some() {
                break;
            }
        }
    }
    let sunflower_free_k3 = sf_witness.is_none();
    if !sunflower_free_k3 {
        return Err(format!(
            "invariant 10 violation: chosen family contains a 3-sunflower (witness triple indices {:?})",
            sf_witness
        ));
    }

    // Build core_rows. Iterate cores in enumerator.cores order.
    // Per packet "every core C (size 1..w-1, plus s=0 if you have time)".
    // We include s = 0 .. w-1 so we cover s=0,1,...,w-1.
    // Compute candidate_count_total etc. from the rows we emit.
    let mut core_rows_json = String::from("[");
    let mut total_candidate_count: u128 = 0;
    let mut total_local_safe: u128 = 0;
    let mut total_local_unsafe: u128 = 0;
    let mut per_core_summary: Vec<(usize, usize, u128, u128, u128, u128)> = Vec::new(); // (s, count, candidate, local_safe, local_unsafe, channel_count)
    let mut emitted_core_id: usize = 0;

    let family_set: Vec<usize> = family.clone();
    let family_mask_u128: FamilyMask = {
        let mut fm: FamilyMask = 0;
        for &i in family.iter() {
            fm |= 1u128 << i;
        }
        fm
    };
    let unused_u128: FamilyMask = enumerator.all_site_mask & !family_mask_u128;

    let mut first_core = true;
    for (core_idx, &core_mask) in enumerator.cores.iter().enumerate() {
        let s = enumerator.core_sizes[core_idx];
        if s >= enumerator.w {
            continue;
        }
        // Candidates through this core, restricted to unused (i.e. not in family).
        let candidates_mask = enumerator.core_candidate_masks[core_idx] & unused_u128;
        if candidates_mask == 0 {
            continue;
        }
        emitted_core_id += 1;

        let core_elements = mask_to_elements_1based(core_mask, n);
        if core_elements.len() != s {
            return Err(format!(
                "invariant: core_size_s mismatch (s={} vs |core_elements|={})",
                s,
                core_elements.len()
            ));
        }

        // Iterate candidate site indices.
        let mut cand_jsons: Vec<String> = Vec::new();
        let mut local_safe_count: u128 = 0;
        let mut local_unsafe_count: u128 = 0;
        let mut candidate_count: u128 = 0;
        let mut seen_cand_masks: Vec<SiteMask> = Vec::new();
        for cand_idx in 0..enumerator.sites.len() {
            if (candidates_mask & (1u128 << cand_idx)) == 0 {
                continue;
            }
            let cand_set = enumerator.sites[cand_idx];

            // Invariant 6: candidate uniqueness within core row.
            if seen_cand_masks.contains(&cand_set) {
                return Err(format!(
                    "invariant 6 violation: duplicate candidate mask {} in core_idx={}",
                    cand_set, core_idx
                ));
            }
            seen_cand_masks.push(cand_set);

            // Invariant 7: core ⊆ candidate.
            if (cand_set & core_mask) != core_mask {
                return Err("invariant 7 violation: core not subset of candidate".to_string());
            }
            // Invariant 8: candidate has size w and is not in family.
            if cand_set.count_ones() as usize != enumerator.w {
                return Err("invariant 8 violation: candidate size != w".to_string());
            }
            if family_set.iter().any(|&i| enumerator.sites[i] == cand_set) {
                return Err("invariant 8 violation: candidate is in family".to_string());
            }

            let cand_elements = mask_to_elements_1based(cand_set, n);
            // Invariant 1 element range.
            for &e in &cand_elements {
                if e < 1 || e > n {
                    return Err(format!("invariant 1 violation: candidate element {} out of [1..{}]", e, n));
                }
            }
            // Invariant 9: petal = candidate \ core.
            let petal_mask: SiteMask = cand_set & !core_mask;
            let petal_elements = mask_to_elements_1based(petal_mask, n);
            if petal_elements.len() != enumerator.w - s {
                return Err(format!(
                    "invariant 9 violation: petal_size {} != w-s {}",
                    petal_elements.len(),
                    enumerator.w - s
                ));
            }

            // Local safety check (and if false, witnesses).
            let (is_safe, witnesses) = check_local_safety(enumerator, &family_set, core_mask, cand_idx);
            // Invariant 11: if unsafe, must have witness.
            if !is_safe && witnesses.is_empty() {
                return Err("invariant 11 violation: unsafe candidate without witness".to_string());
            }
            // Invariant 12: if safe, no witness exists (we just verified by exhaustion above).
            if is_safe && !witnesses.is_empty() {
                return Err("invariant 12 violation: safe candidate has witness pairs".to_string());
            }

            // Build witness JSON.
            let mut witness_json = String::from("[");
            for (wi, (a1_idx, a2_idx)) in witnesses.iter().enumerate() {
                if wi > 0 {
                    witness_json.push(',');
                }
                // Lookup S-id for a1, a2 within sets list.
                let a1_pos = family_set.iter().position(|&i| i == *a1_idx).ok_or("a1 not in family")?;
                let a2_pos = family_set.iter().position(|&i| i == *a2_idx).ok_or("a2 not in family")?;
                let a1_elements = mask_to_elements_1based(enumerator.sites[*a1_idx], n);
                let a2_elements = mask_to_elements_1based(enumerator.sites[*a2_idx], n);
                // Recheck checkable certificates.
                let pic = (enumerator.sites[*a1_idx] & enumerator.sites[*a2_idx]) == core_mask
                    && (enumerator.sites[*a1_idx] & cand_set) == core_mask
                    && (enumerator.sites[*a2_idx] & cand_set) == core_mask;
                let petal_a1 = enumerator.sites[*a1_idx] & !core_mask;
                let petal_a2 = enumerator.sites[*a2_idx] & !core_mask;
                let petal_b = cand_set & !core_mask;
                let pdj = (petal_a1 & petal_a2) == 0
                    && (petal_a1 & petal_b) == 0
                    && (petal_a2 & petal_b) == 0;
                if !pic || !pdj {
                    return Err("invariant 11 violation: witness certificates fail recheck".to_string());
                }
                witness_json.push_str(&format!(
                    "{{\"a1_set_id\":\"S-{:04}\",\"a1_elements\":{},\"a2_set_id\":\"S-{:04}\",\"a2_elements\":{},\"pairwise_intersection_core\":true,\"petals_pairwise_disjoint\":true}}",
                    a1_pos + 1,
                    json_array_usize(&a1_elements),
                    a2_pos + 1,
                    json_array_usize(&a2_elements)
                ));
            }
            witness_json.push(']');

            let unsafe_reason = if is_safe {
                "null".to_string()
            } else {
                "\"closes_3_sunflower_through_core\"".to_string()
            };

            let cand_id = format!("B-{:04}", cand_jsons.len() + 1);
            cand_jsons.push(format!(
                "{{\"candidate_id\":\"{}\",\"elements\":{},\"mask_u64\":{},\"petal_elements\":{},\"petal_mask_u64\":{},\"contains_core\":true,\"uniform_card_w\":true,\"unused_by_family\":true,\"local_safe_through_core\":{},\"unsafe_reason\":{},\"unsafe_witness_pairs\":{}}}",
                cand_id,
                json_array_usize(&cand_elements),
                cand_set,
                json_array_usize(&petal_elements),
                petal_mask,
                if is_safe { "true" } else { "false" },
                unsafe_reason,
                witness_json
            ));
            candidate_count += 1;
            if is_safe {
                local_safe_count += 1;
            } else {
                local_unsafe_count += 1;
            }
        }

        // Invariant 13.
        if candidate_count != local_safe_count + local_unsafe_count {
            return Err(format!(
                "invariant 13 violation: candidate_count {} != local_safe {} + local_unsafe {}",
                candidate_count, local_safe_count, local_unsafe_count
            ));
        }
        // Invariant 14.
        if cand_jsons.len() as u128 != candidate_count {
            return Err(format!(
                "invariant 14 violation: candidate_rows.length {} != candidate_count {}",
                cand_jsons.len(),
                candidate_count
            ));
        }

        total_candidate_count += candidate_count;
        total_local_safe += local_safe_count;
        total_local_unsafe += local_unsafe_count;
        per_core_summary.push((s, emitted_core_id, candidate_count, local_safe_count, local_unsafe_count, 1));

        if !first_core {
            core_rows_json.push(',');
        }
        first_core = false;
        core_rows_json.push_str(&format!(
            "{{\"core_id\":\"C-{:04}\",\"core_elements\":{},\"core_mask_u64\":{},\"core_size_s\":{},\"petal_size\":{},\"candidate_count\":{},\"local_safe_count\":{},\"local_unsafe_count\":{},\"candidate_rows\":[{}]}}",
            emitted_core_id,
            json_array_usize(&core_elements),
            core_mask,
            s,
            enumerator.w - s,
            candidate_count,
            local_safe_count,
            local_unsafe_count,
            cand_jsons.join(",")
        ));
    }
    core_rows_json.push(']');

    // Build summary_rows. Aggregate per core_size_s from the candidate rows we emitted.
    // Invariant 15: derived strictly from candidate_rows.
    let mut by_s: Vec<(usize, u128, u128, u128, u128)> = Vec::new(); // (s, channel_count, candidate, local_safe, local_unsafe)
    for &(s, _, c, sf, un, _) in &per_core_summary {
        if let Some(entry) = by_s.iter_mut().find(|(ss, ..)| *ss == s) {
            entry.1 += 1;
            entry.2 += c;
            entry.3 += sf;
            entry.4 += un;
        } else {
            by_s.push((s, 1, c, sf, un));
        }
    }

    let mut summary_rows_json = String::from("[");
    for (i, (s, channel_count, cand, sf, un)) in by_s.iter().enumerate() {
        if i > 0 {
            summary_rows_json.push(',');
        }
        let i_local = neg_log2_ratio(*sf, *cand);
        summary_rows_json.push_str(&format!(
            "{{\"n\":{},\"w\":{},\"k\":{},\"m\":{},\"core_size_s\":{},\"channel_count\":{},\"candidate_count\":{},\"local_safe_count\":{},\"local_unsafe_count\":{},\"I_core_local_bits\":{},\"claim_ceiling\":\"shadow signature, not universal law; no theorem or lower-bound progress\"}}",
            n,
            enumerator.w,
            k,
            target_m,
            s,
            channel_count,
            cand,
            sf,
            un,
            fmt_opt(i_local)
        ));
    }
    summary_rows_json.push(']');

    // Family rows wrapper.
    let family_id = "F-00000001".to_string();
    let mut family_rows_json = String::new();
    family_rows_json.push_str(&format!(
        "[{{\"family_id\":\"{}\",\"family_hash_sha256\":\"{}\",\"m\":{},\"sets\":{},\"sunflower_free_k3\":{},\"sunflower_free_check\":{{\"method\":\"exhaustive_triple_scan\",\"checked_triples\":{},\"forbidden_witness\":null}},\"core_rows\":{}}}]",
        family_id,
        family_hash,
        target_m,
        sets_json,
        sunflower_free_k3,
        checked_triples,
        core_rows_json
    ));

    // Runner block.
    let binary_sha256 = read_self_binary_sha256().unwrap_or_else(|| "null".to_string());
    let raw_args: Vec<String> = env::args().collect();
    let command = raw_args.join(" ");
    let runner_json = format!(
        "{{\"name\":\"sunflower-core-closure-rust\",\"source_path\":\"erdos-experiments/Erdos20/rust_core_closure/src/main.rs\",\"binary_sha256\":\"{}\",\"command\":\"{}\",\"deterministic\":true}}",
        binary_sha256,
        command.replace('\\', "\\\\").replace('"', "\\\"")
    );

    let parameters_json = format!(
        "{{\"n\":{},\"w\":{},\"k\":{},\"ground_set\":{},\"ground_set_indexing\":\"1-based\",\"target_m\":{},\"jamming_m_star\":{},\"window_offset\":0}}",
        n,
        enumerator.w,
        k,
        json_array_usize(&((1..=n).collect::<Vec<usize>>())),
        target_m,
        jamming_m_star
    );

    // Derive experiment_id from out_path basename.
    let experiment_id = derive_experiment_id_from_out(&out_path);

    let generated_at = iso8601_now_utc();
    let final_json = format!(
        "{{\"schema_version\":\"erdos20.safe_candidates.v1\",\"experiment_id\":\"{}\",\"generated_at_utc\":\"{}\",\"runner\":{},\"parameters\":{},\"family_rows\":{},\"summary_rows\":{}}}",
        experiment_id, generated_at, runner_json, parameters_json, family_rows_json, summary_rows_json
    );

    // Final cross-check: aggregate totals must match candidate-row totals.
    let mut sum_candidate: u128 = 0;
    let mut sum_safe: u128 = 0;
    let mut sum_unsafe: u128 = 0;
    for (_, _, c, sf, un) in &by_s {
        sum_candidate += c;
        sum_safe += sf;
        sum_unsafe += un;
    }
    if sum_candidate != total_candidate_count || sum_safe != total_local_safe || sum_unsafe != total_local_unsafe {
        return Err("invariant 15 violation: summary totals diverge from candidate-row totals".to_string());
    }

    fs::write(&out_path, final_json.as_bytes()).map_err(|e| format!("write {}: {}", out_path, e))?;

    // Brief summary to stdout.
    println!(
        "{{\"emit_safe_candidates_v1\":true,\"out_path\":\"{}\",\"experiment_id\":\"{}\",\"family_id\":\"{}\",\"target_m\":{},\"core_count_emitted\":{},\"candidate_count_total\":{},\"local_safe_count_total\":{},\"local_unsafe_count_total\":{}}}",
        out_path, experiment_id, family_id, target_m, emitted_core_id, total_candidate_count, total_local_safe, total_local_unsafe
    );
    Ok(())
}

fn derive_experiment_id_from_out(path: &str) -> String {
    let basename = std::path::Path::new(path)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("EXP-MATH-ERDOS20-SAFE-CANDIDATES-UNKNOWN");
    if let Some(stem) = basename.strip_suffix("_RESULTS.json") {
        return stem.to_string();
    }
    if let Some(stem) = basename.strip_suffix(".json") {
        return stem.to_string();
    }
    basename.to_string()
}

fn main() -> Result<(), String> {
    let config = parse_args()?;
    let enumerator = Enumerator::new(config.n, config.w)?;
    if config.emit_safe_candidates_v1 {
        run_emit_safe_candidates_v1(&config, &enumerator)
    } else {
        let (density, density_seconds) = enumerator.density_pass();
        run_default(&config, &enumerator, &density, density_seconds)
    }
}
