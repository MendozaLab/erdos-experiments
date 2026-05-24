use clap::Parser;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::cmp::Ordering;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(name = "sumfree-transfer")]
#[command(about = "PMF transfer-operator scanner for Erdős #166 sum-free sets")]
struct Cli {
    #[arg(long, default_value_t = 20)]
    n_min: usize,

    #[arg(long, default_value_t = 40)]
    n_max: usize,

    #[arg(long, default_value_t = 5)]
    frontier_k: usize,

    #[arg(long, default_value_t = 0)]
    prune_deficiency: usize,

    #[arg(long)]
    experiment_id: Option<String>,

    #[arg(long, default_value = "erdos-experiments/results/erdos-166")]
    output_dir: PathBuf,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct State {
    occupied_mask: u128,
    cardinality: u8,
}

impl State {
    fn empty() -> Self {
        Self {
            occupied_mask: 0,
            cardinality: 0,
        }
    }

    fn occupied_sites(&self) -> Vec<usize> {
        bits_to_sites(self.occupied_mask)
    }

    fn try_occupy_sum_free(&self, x: usize) -> Option<Self> {
        if x == 0 || self.occupied_mask & bit(x)? != 0 {
            return None;
        }
        let mut occupied = self.occupied_mask;
        while occupied != 0 {
            let a = occupied.trailing_zeros() as usize;
            if a <= x {
                let b = x - a;
                if let Some(b_bit) = bit(b) {
                    if self.occupied_mask & b_bit != 0 {
                        return None;
                    }
                }
            }
            occupied &= occupied - 1;
        }
        Some(Self {
            occupied_mask: self.occupied_mask | bit(x)?,
            cardinality: self.cardinality + 1,
        })
    }
}

#[derive(Clone, Serialize)]
struct FrontierWitness {
    rank: usize,
    witness: Vec<usize>,
    objective_value: f64,
    uniform_prefix_residual: f64,
    uniform_prefix_ratio_n: f64,
    prefix_location: Option<usize>,
    density_adjusted_mass_dev: f64,
    density_adjusted_mass_ratio_n2: f64,
    joint_score: f64,
}

#[derive(Clone)]
struct FrontierCandidate {
    witness: Vec<usize>,
    objective_value: f64,
    uniform_prefix_residual: f64,
    uniform_prefix_ratio_n: f64,
    prefix_location: Option<usize>,
    density_adjusted_mass_dev: f64,
    density_adjusted_mass_ratio_n2: f64,
    joint_score: f64,
}

struct TopKAccumulator {
    k: usize,
    entries: Vec<FrontierCandidate>,
}

impl TopKAccumulator {
    fn new(k: usize) -> Self {
        Self {
            k,
            entries: Vec::new(),
        }
    }

    fn observe(&mut self, candidate: FrontierCandidate) {
        if self.k == 0 {
            return;
        }
        self.entries.push(candidate);
        self.entries.sort_by(candidate_cmp);
        self.entries.truncate(self.k);
    }

    fn ranked(&self) -> Vec<FrontierWitness> {
        self.entries
            .iter()
            .enumerate()
            .map(|(idx, item)| FrontierWitness {
                rank: idx + 1,
                witness: item.witness.clone(),
                objective_value: item.objective_value,
                uniform_prefix_residual: item.uniform_prefix_residual,
                uniform_prefix_ratio_n: item.uniform_prefix_ratio_n,
                prefix_location: item.prefix_location,
                density_adjusted_mass_dev: item.density_adjusted_mass_dev,
                density_adjusted_mass_ratio_n2: item.density_adjusted_mass_ratio_n2,
                joint_score: item.joint_score,
            })
            .collect()
    }
}

#[derive(Serialize)]
struct TopKFrontier {
    k: usize,
    zero_temperature_field_interpretation: String,
    prefix_best: Vec<FrontierWitness>,
    density_adjusted_mass_best: Vec<FrontierWitness>,
    joint_best: Vec<FrontierWitness>,
}

#[derive(Serialize)]
struct DeficiencyLayer {
    deficiency: usize,
    cardinality: usize,
    count: u64,
}

#[derive(Serialize)]
struct SpectralObservables {
    ground_cardinality: usize,
    ground_state_degeneracy: u64,
    ground_entropy_ln: f64,
    deficiency_layers: Vec<DeficiencyLayer>,
}

#[derive(Serialize)]
struct TransferDiagnostics {
    traversal_mode: String,
    sites_processed: usize,
    terminal_state_count: u64,
    transitions_attempted: u64,
    legal_occupy_transitions: u64,
    rejected_occupy_transitions: u64,
    prune_deficiency: usize,
    pruning_target_h_n: usize,
    pruning_min_cardinality: usize,
    reachability_pruned_states: u64,
}

#[derive(Serialize)]
struct ParityReference {
    exact_h_n_formula: usize,
    h_n_matches_formula: bool,
}

#[derive(Serialize)]
struct Row {
    n: usize,
    h_n: usize,
    maximizer_count: u64,
    first_maximizer_witness: Vec<usize>,
    spectral_observables: SpectralObservables,
    transfer_diagnostics: TransferDiagnostics,
    parity_reference: ParityReference,
    top_k_frontier: TopKFrontier,
}

#[derive(Serialize)]
struct Results {
    experiment_id: String,
    date: String,
    erdos_problem: usize,
    scan_mode: String,
    implementation: String,
    constraint: String,
    n_min: usize,
    n_max: usize,
    frontier_k: usize,
    prune_deficiency: usize,
    total_runtime_sec: f64,
    rows: Vec<Row>,
    summary: Summary,
}

#[derive(Serialize)]
struct Summary {
    checked_count: usize,
    h_formula_match_count: usize,
    mismatch_ns: Vec<usize>,
    split_count: usize,
    split_ns: Vec<usize>,
}

fn bit(index: usize) -> Option<u128> {
    if index < 128 {
        Some(1u128 << index)
    } else {
        None
    }
}

fn bits_to_sites(mut mask: u128) -> Vec<usize> {
    let mut sites = Vec::new();
    while mask != 0 {
        let site = mask.trailing_zeros() as usize;
        sites.push(site);
        mask &= mask - 1;
    }
    sites
}

fn candidate_cmp(a: &FrontierCandidate, b: &FrontierCandidate) -> Ordering {
    a.objective_value
        .total_cmp(&b.objective_value)
        .then_with(|| a.witness.cmp(&b.witness))
}

fn reference_h_n(n: usize) -> usize {
    (n + 1) / 2
}

fn compute_uniform_prefix_residual(chosen: &[usize], n: usize) -> (f64, Option<usize>) {
    if chosen.is_empty() {
        return (n as f64, Some(n));
    }
    let slope = n as f64 / chosen.len() as f64;
    let mut prefix_size = 0usize;
    let mut next_index = 0usize;
    let mut max_dev = f64::NEG_INFINITY;
    let mut max_t = 0usize;
    for t in 1..=n {
        while next_index < chosen.len() && chosen[next_index] <= t {
            prefix_size += 1;
            next_index += 1;
        }
        let dev = (t as f64 - prefix_size as f64 * slope).abs();
        if dev > max_dev {
            max_dev = dev;
            max_t = t;
        }
    }
    (max_dev, Some(max_t))
}

fn density_adjusted_mass_dev(chosen: &[usize], n: usize) -> f64 {
    let mass: usize = chosen.iter().sum();
    let center = n as f64 * (chosen.len() as f64 + 1.0) / 2.0;
    (mass as f64 - center).abs()
}

fn frontier_candidate(chosen: Vec<usize>, n: usize) -> FrontierCandidate {
    let (prefix_residual, prefix_location) = compute_uniform_prefix_residual(&chosen, n);
    let mass_dev = density_adjusted_mass_dev(&chosen, n);
    let n_real = n as f64;
    let prefix_ratio = prefix_residual / n_real.max(1.0);
    let mass_ratio = mass_dev / n_real.max(1.0).powi(2);
    FrontierCandidate {
        witness: chosen,
        objective_value: prefix_ratio + mass_ratio,
        uniform_prefix_residual: prefix_residual,
        uniform_prefix_ratio_n: prefix_ratio,
        prefix_location,
        density_adjusted_mass_dev: mass_dev,
        density_adjusted_mass_ratio_n2: mass_ratio,
        joint_score: prefix_ratio + mass_ratio,
    }
}

fn analyze_n(n: usize, frontier_k: usize, prune_deficiency: usize) -> Row {
    assert!(n < 128, "u128 bitset transfer operator requires n < 128");
    let target_h = reference_h_n(n);
    let min_cardinality = target_h.saturating_sub(prune_deficiency);
    let mut layer_counts = vec![0u64; target_h + 1];
    let mut terminal_state_count = 0u64;
    let mut transitions_attempted = 0u64;
    let mut legal_occupy_transitions = 0u64;
    let mut rejected_occupy_transitions = 0u64;
    let mut reachability_pruned_states = 0u64;
    let mut over_target_pruned_states = 0u64;
    let mut first_maximizer_witness = Vec::new();
    let mut prefix_frontier = TopKAccumulator::new(frontier_k);
    let mut mass_frontier = TopKAccumulator::new(frontier_k);
    let mut joint_frontier = TopKAccumulator::new(frontier_k);

    struct DfsCtx<'a> {
        n: usize,
        target_h: usize,
        min_cardinality: usize,
        layer_counts: &'a mut [u64],
        terminal_state_count: &'a mut u64,
        transitions_attempted: &'a mut u64,
        legal_occupy_transitions: &'a mut u64,
        rejected_occupy_transitions: &'a mut u64,
        reachability_pruned_states: &'a mut u64,
        over_target_pruned_states: &'a mut u64,
        first_maximizer_witness: &'a mut Vec<usize>,
        prefix_frontier: &'a mut TopKAccumulator,
        mass_frontier: &'a mut TopKAccumulator,
        joint_frontier: &'a mut TopKAccumulator,
    }

    fn dfs(x: usize, state: State, ctx: &mut DfsCtx<'_>) {
        if state.cardinality as usize > ctx.target_h {
            *ctx.over_target_pruned_states += 1;
            return;
        }
        let remaining_sites = if x <= ctx.n { ctx.n - x + 1 } else { 0 };
        if state.cardinality as usize + remaining_sites < ctx.min_cardinality {
            *ctx.reachability_pruned_states += 1;
            return;
        }
        if x > ctx.n {
            let cardinality = state.cardinality as usize;
            if cardinality < ctx.min_cardinality {
                return;
            }
            *ctx.terminal_state_count += 1;
            ctx.layer_counts[cardinality] += 1;
            if cardinality == ctx.target_h {
                let witness = state.occupied_sites();
                if ctx.first_maximizer_witness.is_empty() || witness < *ctx.first_maximizer_witness
                {
                    *ctx.first_maximizer_witness = witness.clone();
                }
                let candidate = frontier_candidate(witness, ctx.n);
                ctx.prefix_frontier.observe(FrontierCandidate {
                    objective_value: candidate.uniform_prefix_residual,
                    ..candidate.clone()
                });
                ctx.mass_frontier.observe(FrontierCandidate {
                    objective_value: candidate.density_adjusted_mass_dev,
                    ..candidate.clone()
                });
                ctx.joint_frontier.observe(candidate);
            }
            return;
        }

        dfs(x + 1, state, ctx);
        *ctx.transitions_attempted += 1;
        match state.try_occupy_sum_free(x) {
            Some(new_state) => {
                *ctx.legal_occupy_transitions += 1;
                dfs(x + 1, new_state, ctx);
            }
            None => {
                *ctx.rejected_occupy_transitions += 1;
            }
        }
    }

    let mut ctx = DfsCtx {
        n,
        target_h,
        min_cardinality,
        layer_counts: &mut layer_counts,
        terminal_state_count: &mut terminal_state_count,
        transitions_attempted: &mut transitions_attempted,
        legal_occupy_transitions: &mut legal_occupy_transitions,
        rejected_occupy_transitions: &mut rejected_occupy_transitions,
        reachability_pruned_states: &mut reachability_pruned_states,
        over_target_pruned_states: &mut over_target_pruned_states,
        first_maximizer_witness: &mut first_maximizer_witness,
        prefix_frontier: &mut prefix_frontier,
        mass_frontier: &mut mass_frontier,
        joint_frontier: &mut joint_frontier,
    };
    dfs(1, State::empty(), &mut ctx);
    reachability_pruned_states += over_target_pruned_states;

    let h_n = layer_counts
        .iter()
        .enumerate()
        .rev()
        .find(|(_, count)| **count > 0)
        .map(|(k, _)| k)
        .unwrap_or(0);
    let maximizer_count = layer_counts[target_h];
    let deficiency_layers: Vec<DeficiencyLayer> = (0..=prune_deficiency.min(3))
        .filter_map(|deficiency| {
            target_h
                .checked_sub(deficiency)
                .map(|cardinality| DeficiencyLayer {
                    deficiency,
                    cardinality,
                    count: layer_counts[cardinality],
                })
        })
        .collect();
    let parity_reference = ParityReference {
        exact_h_n_formula: target_h,
        h_n_matches_formula: h_n == target_h,
    };

    Row {
        n,
        h_n,
        maximizer_count,
        first_maximizer_witness,
        spectral_observables: SpectralObservables {
            ground_cardinality: h_n,
            ground_state_degeneracy: maximizer_count,
            ground_entropy_ln: (maximizer_count as f64).ln(),
            deficiency_layers,
        },
        transfer_diagnostics: TransferDiagnostics {
            traversal_mode: "sum_free_layer_pruned_dfs".to_string(),
            sites_processed: n,
            terminal_state_count,
            transitions_attempted,
            legal_occupy_transitions,
            rejected_occupy_transitions,
            prune_deficiency,
            pruning_target_h_n: target_h,
            pruning_min_cardinality: min_cardinality,
            reachability_pruned_states,
        },
        parity_reference,
        top_k_frontier: TopKFrontier {
            k: frontier_k,
            zero_temperature_field_interpretation:
                "prefix_best, density_adjusted_mass_best, and joint_best are zero-temperature ground-state selections under small positive prefix/mass fields after cardinality is fixed at h(n)"
                    .to_string(),
            prefix_best: prefix_frontier.ranked(),
            density_adjusted_mass_best: mass_frontier.ranked(),
            joint_best: joint_frontier.ranked(),
        },
    }
}

fn summarize(rows: &[Row]) -> Summary {
    let mut mismatch_ns = Vec::new();
    let mut split_ns = Vec::new();
    for row in rows {
        if !row.parity_reference.h_n_matches_formula {
            mismatch_ns.push(row.n);
        }
        let prefix = row.top_k_frontier.prefix_best.first().map(|w| &w.witness);
        let mass = row
            .top_k_frontier
            .density_adjusted_mass_best
            .first()
            .map(|w| &w.witness);
        let joint = row.top_k_frontier.joint_best.first().map(|w| &w.witness);
        if prefix != mass || prefix != joint {
            split_ns.push(row.n);
        }
    }
    Summary {
        checked_count: rows.len(),
        h_formula_match_count: rows.len() - mismatch_ns.len(),
        mismatch_ns,
        split_count: split_ns.len(),
        split_ns,
    }
}

fn build_report(results: &Results) -> String {
    let mut lines = Vec::new();
    lines.push(format!(
        "# {} — Sum-Free Transfer-Operator Cross-Problem Gate",
        results.experiment_id
    ));
    lines.push(String::new());
    lines.push("## Question".to_string());
    lines.push(String::new());
    lines.push(
        "Does the same PMF state-machine API survive a changed local exclusion rule?".to_string(),
    );
    lines.push(String::new());
    lines.push("## Result".to_string());
    lines.push(String::new());
    lines.push(format!(
        "Formula parity for `h(n)=ceil(n/2)` matched {}/{} rows. Split frontiers occurred in {} rows: {:?}.",
        results.summary.h_formula_match_count,
        results.summary.checked_count,
        results.summary.split_count,
        results.summary.split_ns
    ));
    lines.push(String::new());
    lines.push("| n | h(n) | degeneracy | entropy ln | terminal retained | pruned states | prefix winner | mass winner | joint winner |".to_string());
    lines.push("|---|---:|---:|---:|---:|---:|---|---|---|".to_string());
    for row in &results.rows {
        let prefix = row
            .top_k_frontier
            .prefix_best
            .first()
            .map(|w| format!("{:?}", w.witness))
            .unwrap_or_else(|| "[]".to_string());
        let mass = row
            .top_k_frontier
            .density_adjusted_mass_best
            .first()
            .map(|w| format!("{:?}", w.witness))
            .unwrap_or_else(|| "[]".to_string());
        let joint = row
            .top_k_frontier
            .joint_best
            .first()
            .map(|w| format!("{:?}", w.witness))
            .unwrap_or_else(|| "[]".to_string());
        lines.push(format!(
            "| {} | {} | {} | {:.6} | {} | {} | `{}` | `{}` | `{}` |",
            row.n,
            row.h_n,
            row.maximizer_count,
            row.spectral_observables.ground_entropy_ln,
            row.transfer_diagnostics.terminal_state_count,
            row.transfer_diagnostics.reachability_pruned_states,
            prefix,
            mass,
            joint
        ));
    }
    lines.push(String::new());
    lines.push("## Claim Boundary".to_string());
    lines.push(String::new());
    lines.push(
        "This is a cross-problem API test, not a theorem about the original Erdős #166 formulation."
            .to_string(),
    );
    lines.push(
        "It supports the narrow PMF claim that the transfer-state method ports from Sidon to a changed additive exclusion rule."
            .to_string(),
    );
    lines.join("\n") + "\n"
}

fn write_outputs(results: &Results, output_dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(output_dir)?;
    let json_path = output_dir.join(format!("{}_RESULTS.json", results.experiment_id));
    let report_path = output_dir.join(format!("{}_REPORT.md", results.experiment_id));
    let sha_path = output_dir.join(format!("{}_RESULTS.sha256", results.experiment_id));
    if json_path.exists() || report_path.exists() || sha_path.exists() {
        return Err(format!(
            "refusing to overwrite existing packet for {}",
            results.experiment_id
        )
        .into());
    }

    let json_payload = serde_json::to_string_pretty(results)? + "\n";
    fs::write(&json_path, json_payload.as_bytes())?;
    fs::write(&report_path, build_report(results))?;
    let digest = Sha256::digest(json_payload.as_bytes());
    fs::write(
        &sha_path,
        format!(
            "{:x}  {}\n",
            digest,
            json_path.file_name().unwrap().to_string_lossy()
        ),
    )?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    if cli.n_min > cli.n_max {
        return Err("--n-min must be <= --n-max".into());
    }
    if cli.n_max >= 128 {
        return Err("u128 bitset transfer operator currently requires --n-max < 128".into());
    }
    let experiment_id = cli.experiment_id.unwrap_or_else(|| {
        format!(
            "EXP-MM-166-PMF-SUMFREE-TRANSFER-{}-{}-2026-04-29",
            cli.n_min, cli.n_max
        )
    });
    let t0 = Instant::now();
    let rows: Vec<Row> = (cli.n_min..=cli.n_max)
        .map(|n| analyze_n(n, cli.frontier_k, cli.prune_deficiency))
        .collect();
    let summary = summarize(&rows);
    let results = Results {
        experiment_id,
        date: "2026-04-29".to_string(),
        erdos_problem: 166,
        scan_mode: "pmf_sum_free_transfer_operator".to_string(),
        implementation: "Rust bitset layer-pruned transfer-state enumerator".to_string(),
        constraint: "sum_free_no_x_plus_y_equals_z".to_string(),
        n_min: cli.n_min,
        n_max: cli.n_max,
        frontier_k: cli.frontier_k,
        prune_deficiency: cli.prune_deficiency,
        total_runtime_sec: t0.elapsed().as_secs_f64(),
        rows,
        summary,
    };
    write_outputs(&results, &cli.output_dir)?;
    println!(
        "wrote packet {} to {}",
        results.experiment_id,
        cli.output_dir.display()
    );
    Ok(())
}
