use clap::Parser;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::cmp::Ordering;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(name = "b2g-transfer")]
#[command(about = "PMF transfer-state exact scout for Erdős #755 B_2[g] sets")]
struct Cli {
    #[arg(long, default_value_t = 20)]
    n_min: usize,

    #[arg(long, default_value_t = 40)]
    n_max: usize,

    #[arg(long, default_value_t = 2)]
    g: u16,

    #[arg(long, default_value_t = 5)]
    frontier_k: usize,

    #[arg(long, default_value_t = 0)]
    near_ground_deficiency: usize,

    #[arg(long, default_value_t = 0)]
    near_ground_terminal_cap: u64,

    #[arg(long, default_value_t = 0)]
    near_ground_layer_cap: u64,

    #[arg(long, default_value_t = false)]
    skip_near_ground: bool,

    #[arg(long)]
    known_h: Option<usize>,

    #[arg(long)]
    known_ground_count: Option<u64>,

    #[arg(long, default_value_t = 0)]
    progress_every: u64,

    #[arg(long)]
    experiment_id: Option<String>,

    #[arg(long, default_value = "erdos-experiments/results/erdos-755")]
    output_dir: PathBuf,
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

    fn clear(&mut self) {
        self.entries.clear();
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
struct StateModel {
    representation: String,
    occupied_lattice: String,
    memory: String,
    transition_rule: String,
}

#[derive(Serialize)]
struct TransferDiagnostics {
    traversal_mode: String,
    sites_processed: usize,
    nodes_visited: u64,
    transitions_attempted: u64,
    legal_occupy_transitions: u64,
    rejected_occupy_transitions: u64,
    branch_bound_pruned_nodes: u64,
    best_updates: u64,
    progress_every_nodes: u64,
}

#[derive(Serialize)]
struct DeficiencyLayer {
    deficiency: usize,
    cardinality: usize,
    count: u64,
}

#[derive(Serialize)]
struct NearGroundDiagnostics {
    requested_deficiency: usize,
    min_retained_cardinality: usize,
    second_pass_nodes_visited: u64,
    second_pass_pruned_nodes: u64,
    retained_terminal_states: u64,
    terminal_cap: u64,
    terminal_cap_hit: bool,
    layer_cap: u64,
    layer_cap_hit: bool,
    counts_exact: bool,
    progress_every_nodes: u64,
}

#[derive(Serialize)]
struct SpectralObservables {
    ground_cardinality_h_n: usize,
    ground_state_degeneracy: u64,
    ground_entropy_ln: f64,
    first_excited_cardinality: Option<usize>,
    first_excited_count: Option<u64>,
    gap_to_first_excited_layer: Option<usize>,
    near_ground_count_h_minus_1: u64,
    near_ground_count_h_minus_2: u64,
    deficiency_layers: Vec<DeficiencyLayer>,
    near_ground_diagnostics: NearGroundDiagnostics,
}

#[derive(Serialize)]
struct Row {
    n: usize,
    g: u16,
    ordered_representation_cap: u16,
    h_n_exact_finite: usize,
    maximizer_count: u64,
    first_maximizer_witness: Vec<usize>,
    state_model: StateModel,
    transfer_diagnostics: TransferDiagnostics,
    spectral_observables: SpectralObservables,
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
    claim_boundary: String,
    n_min: usize,
    n_max: usize,
    g: u16,
    frontier_k: usize,
    near_ground_deficiency: usize,
    near_ground_terminal_cap: u64,
    near_ground_layer_cap: u64,
    skip_near_ground: bool,
    total_runtime_sec: f64,
    rows: Vec<Row>,
    summary: Summary,
}

#[derive(Serialize)]
struct Summary {
    checked_count: usize,
    min_h_n: usize,
    max_h_n: usize,
    split_count: usize,
    split_ns: Vec<usize>,
    total_maximizers: u64,
}

fn candidate_cmp(a: &FrontierCandidate, b: &FrontierCandidate) -> Ordering {
    a.objective_value
        .total_cmp(&b.objective_value)
        .then_with(|| a.witness.cmp(&b.witness))
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
    for t in 0..=n {
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
    let center = n as f64 * chosen.len() as f64 / 2.0;
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

fn candidate_from_current(chosen: &[usize], n: usize) -> FrontierCandidate {
    frontier_candidate(chosen.to_vec(), n)
}

fn legal_increments(
    chosen: &[usize],
    sum_counts: &[u16],
    x: usize,
    cap: u16,
) -> Option<Vec<(usize, u16)>> {
    let mut increments = Vec::with_capacity(chosen.len() + 1);
    let self_sum = 2 * x;
    if sum_counts[self_sum] + 1 > cap {
        return None;
    }
    increments.push((self_sum, 1));

    for &a in chosen {
        let sum = a + x;
        if sum_counts[sum] + 2 > cap {
            return None;
        }
        increments.push((sum, 2));
    }

    Some(increments)
}

fn apply_increments(sum_counts: &mut [u16], increments: &[(usize, u16)]) {
    for &(sum, inc) in increments {
        sum_counts[sum] += inc;
    }
}

fn undo_increments(sum_counts: &mut [u16], increments: &[(usize, u16)]) {
    for &(sum, inc) in increments {
        sum_counts[sum] -= inc;
    }
}

fn greedy_lower_bound(n: usize, cap: u16, descending: bool) -> Vec<usize> {
    let mut chosen = Vec::new();
    let mut sum_counts = vec![0u16; 2 * n + 1];
    let sites: Vec<usize> = if descending {
        (0..=n).rev().collect()
    } else {
        (0..=n).collect()
    };
    for x in sites {
        if let Some(increments) = legal_increments(&chosen, &sum_counts, x, cap) {
            apply_increments(&mut sum_counts, &increments);
            chosen.push(x);
            chosen.sort_unstable();
        }
    }
    chosen
}

fn collect_near_ground_layers(
    n: usize,
    g: u16,
    cap: u16,
    ground_cardinality: usize,
    near_ground_deficiency: usize,
    near_ground_terminal_cap: u64,
    near_ground_layer_cap: u64,
    progress_every: u64,
) -> (Vec<u64>, NearGroundDiagnostics) {
    let min_cardinality = ground_cardinality.saturating_sub(near_ground_deficiency);
    let mut layer_counts = vec![0u64; ground_cardinality + 1];
    let mut chosen = Vec::new();
    let mut sum_counts = vec![0u16; 2 * n + 1];
    let mut nodes_visited = 0u64;
    let mut pruned_nodes = 0u64;
    let mut retained_terminal_states = 0u64;
    let mut terminal_cap_hit = false;
    let mut layer_cap_hit = false;

    fn all_layer_caps_hit(
        layer_counts: &[u64],
        min_cardinality: usize,
        ground_cardinality: usize,
        layer_cap: u64,
    ) -> bool {
        if layer_cap == 0 || min_cardinality >= ground_cardinality {
            return false;
        }
        (min_cardinality..ground_cardinality)
            .all(|cardinality| layer_counts[cardinality] >= layer_cap)
    }

    struct NearCtx<'a> {
        n: usize,
        g: u16,
        cap: u16,
        ground_cardinality: usize,
        min_cardinality: usize,
        layer_counts: &'a mut [u64],
        nodes_visited: &'a mut u64,
        pruned_nodes: &'a mut u64,
        retained_terminal_states: &'a mut u64,
        terminal_cap: u64,
        terminal_cap_hit: &'a mut bool,
        layer_cap: u64,
        layer_cap_hit: &'a mut bool,
        progress_every: u64,
    }

    fn dfs(x: usize, chosen: &mut Vec<usize>, sum_counts: &mut [u16], ctx: &mut NearCtx<'_>) {
        if ctx.terminal_cap > 0 && *ctx.retained_terminal_states >= ctx.terminal_cap {
            *ctx.terminal_cap_hit = true;
            return;
        }
        if all_layer_caps_hit(
            ctx.layer_counts,
            ctx.min_cardinality,
            ctx.ground_cardinality,
            ctx.layer_cap,
        ) {
            *ctx.layer_cap_hit = true;
            return;
        }
        *ctx.nodes_visited += 1;
        if ctx.progress_every > 0 && *ctx.nodes_visited % ctx.progress_every == 0 {
            eprintln!(
                "progress b2g g={} n={} phase=near_ground nodes={} retained={} pruned={}",
                ctx.g, ctx.n, *ctx.nodes_visited, *ctx.retained_terminal_states, *ctx.pruned_nodes
            );
        }
        if chosen.len() > ctx.ground_cardinality {
            *ctx.pruned_nodes += 1;
            return;
        }
        let remaining_sites = if x <= ctx.n { ctx.n - x + 1 } else { 0 };
        if chosen.len() + remaining_sites < ctx.min_cardinality {
            *ctx.pruned_nodes += 1;
            return;
        }
        if x > ctx.n {
            let cardinality = chosen.len();
            if cardinality >= ctx.min_cardinality && cardinality <= ctx.ground_cardinality {
                if ctx.layer_cap == 0 || ctx.layer_counts[cardinality] < ctx.layer_cap {
                    ctx.layer_counts[cardinality] += 1;
                    *ctx.retained_terminal_states += 1;
                } else {
                    *ctx.layer_cap_hit = true;
                }
                if ctx.terminal_cap > 0 && *ctx.retained_terminal_states >= ctx.terminal_cap {
                    *ctx.terminal_cap_hit = true;
                }
                if all_layer_caps_hit(
                    ctx.layer_counts,
                    ctx.min_cardinality,
                    ctx.ground_cardinality,
                    ctx.layer_cap,
                ) {
                    *ctx.layer_cap_hit = true;
                }
            }
            return;
        }

        dfs(x + 1, chosen, sum_counts, ctx);

        if let Some(increments) = legal_increments(chosen, sum_counts, x, ctx.cap) {
            apply_increments(sum_counts, &increments);
            chosen.push(x);
            dfs(x + 1, chosen, sum_counts, ctx);
            chosen.pop();
            undo_increments(sum_counts, &increments);
        }
    }

    let mut ctx = NearCtx {
        n,
        g,
        cap,
        ground_cardinality,
        min_cardinality,
        layer_counts: &mut layer_counts,
        nodes_visited: &mut nodes_visited,
        pruned_nodes: &mut pruned_nodes,
        retained_terminal_states: &mut retained_terminal_states,
        terminal_cap: near_ground_terminal_cap,
        terminal_cap_hit: &mut terminal_cap_hit,
        layer_cap: near_ground_layer_cap,
        layer_cap_hit: &mut layer_cap_hit,
        progress_every,
    };
    dfs(0, &mut chosen, &mut sum_counts, &mut ctx);

    let diagnostics = NearGroundDiagnostics {
        requested_deficiency: near_ground_deficiency,
        min_retained_cardinality: min_cardinality,
        second_pass_nodes_visited: nodes_visited,
        second_pass_pruned_nodes: pruned_nodes,
        retained_terminal_states,
        terminal_cap: near_ground_terminal_cap,
        terminal_cap_hit,
        layer_cap: near_ground_layer_cap,
        layer_cap_hit,
        counts_exact: !terminal_cap_hit && !layer_cap_hit,
        progress_every_nodes: progress_every,
    };
    (layer_counts, diagnostics)
}

fn analyze_n(
    n: usize,
    g: u16,
    frontier_k: usize,
    near_ground_deficiency: usize,
    near_ground_terminal_cap: u64,
    near_ground_layer_cap: u64,
    skip_near_ground: bool,
    known_h: Option<usize>,
    known_ground_count: Option<u64>,
    progress_every: u64,
) -> Row {
    let cap = 2 * g;
    let imported_ground = known_h.zip(known_ground_count);
    let mut best_size = imported_ground
        .map(|(known_h, _)| known_h)
        .unwrap_or_else(|| {
            let ascending = greedy_lower_bound(n, cap, false);
            let descending = greedy_lower_bound(n, cap, true);
            ascending.len().max(descending.len())
        });
    let mut maximizer_count = imported_ground.map(|(_, count)| count).unwrap_or(0);
    let mut first_maximizer_witness = Vec::new();
    let mut prefix_frontier = TopKAccumulator::new(frontier_k);
    let mut mass_frontier = TopKAccumulator::new(frontier_k);
    let mut joint_frontier = TopKAccumulator::new(frontier_k);
    let mut chosen = Vec::new();
    let mut sum_counts = vec![0u16; 2 * n + 1];
    let mut nodes_visited = 0u64;
    let mut transitions_attempted = 0u64;
    let mut legal_occupy_transitions = 0u64;
    let mut rejected_occupy_transitions = 0u64;
    let mut branch_bound_pruned_nodes = 0u64;
    let mut best_updates = 0u64;

    struct DfsCtx<'a> {
        n: usize,
        g: u16,
        cap: u16,
        best_size: &'a mut usize,
        maximizer_count: &'a mut u64,
        first_maximizer_witness: &'a mut Vec<usize>,
        prefix_frontier: &'a mut TopKAccumulator,
        mass_frontier: &'a mut TopKAccumulator,
        joint_frontier: &'a mut TopKAccumulator,
        nodes_visited: &'a mut u64,
        transitions_attempted: &'a mut u64,
        legal_occupy_transitions: &'a mut u64,
        rejected_occupy_transitions: &'a mut u64,
        branch_bound_pruned_nodes: &'a mut u64,
        best_updates: &'a mut u64,
        progress_every: u64,
    }

    fn observe_terminal(chosen: &[usize], ctx: &mut DfsCtx<'_>) {
        let cardinality = chosen.len();
        if cardinality > *ctx.best_size {
            *ctx.best_size = cardinality;
            *ctx.maximizer_count = 0;
            ctx.first_maximizer_witness.clear();
            ctx.prefix_frontier.clear();
            ctx.mass_frontier.clear();
            ctx.joint_frontier.clear();
            *ctx.best_updates += 1;
        }
        if cardinality == *ctx.best_size {
            *ctx.maximizer_count += 1;
            if ctx.first_maximizer_witness.is_empty() || chosen < ctx.first_maximizer_witness {
                *ctx.first_maximizer_witness = chosen.to_vec();
            }
            let candidate = candidate_from_current(chosen, ctx.n);
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
    }

    fn dfs(x: usize, chosen: &mut Vec<usize>, sum_counts: &mut [u16], ctx: &mut DfsCtx<'_>) {
        *ctx.nodes_visited += 1;
        if ctx.progress_every > 0 && *ctx.nodes_visited % ctx.progress_every == 0 {
            eprintln!(
                "progress b2g g={} n={} phase=ground nodes={} best={} max_count={} pruned={}",
                ctx.g,
                ctx.n,
                *ctx.nodes_visited,
                *ctx.best_size,
                *ctx.maximizer_count,
                *ctx.branch_bound_pruned_nodes
            );
        }
        let remaining_sites = if x <= ctx.n { ctx.n - x + 1 } else { 0 };
        if chosen.len() + remaining_sites < *ctx.best_size {
            *ctx.branch_bound_pruned_nodes += 1;
            return;
        }
        if x > ctx.n {
            observe_terminal(chosen, ctx);
            return;
        }

        dfs(x + 1, chosen, sum_counts, ctx);

        *ctx.transitions_attempted += 1;
        match legal_increments(chosen, sum_counts, x, ctx.cap) {
            Some(increments) => {
                *ctx.legal_occupy_transitions += 1;
                apply_increments(sum_counts, &increments);
                chosen.push(x);
                dfs(x + 1, chosen, sum_counts, ctx);
                chosen.pop();
                undo_increments(sum_counts, &increments);
            }
            None => {
                *ctx.rejected_occupy_transitions += 1;
            }
        }
    }

    if imported_ground.is_none() {
        let mut ctx = DfsCtx {
            n,
            g,
            cap,
            best_size: &mut best_size,
            maximizer_count: &mut maximizer_count,
            first_maximizer_witness: &mut first_maximizer_witness,
            prefix_frontier: &mut prefix_frontier,
            mass_frontier: &mut mass_frontier,
            joint_frontier: &mut joint_frontier,
            nodes_visited: &mut nodes_visited,
            transitions_attempted: &mut transitions_attempted,
            legal_occupy_transitions: &mut legal_occupy_transitions,
            rejected_occupy_transitions: &mut rejected_occupy_transitions,
            branch_bound_pruned_nodes: &mut branch_bound_pruned_nodes,
            best_updates: &mut best_updates,
            progress_every,
        };
        dfs(0, &mut chosen, &mut sum_counts, &mut ctx);
    }
    let (layer_counts, near_ground_diagnostics) = if skip_near_ground {
        let mut layer_counts = vec![0u64; best_size + 1];
        layer_counts[best_size] = maximizer_count;
        (
            layer_counts,
            NearGroundDiagnostics {
                requested_deficiency: near_ground_deficiency,
                min_retained_cardinality: best_size.saturating_sub(near_ground_deficiency),
                second_pass_nodes_visited: 0,
                second_pass_pruned_nodes: 0,
                retained_terminal_states: 0,
                terminal_cap: near_ground_terminal_cap,
                terminal_cap_hit: false,
                layer_cap: near_ground_layer_cap,
                layer_cap_hit: false,
                counts_exact: false,
                progress_every_nodes: progress_every,
            },
        )
    } else {
        collect_near_ground_layers(
            n,
            g,
            cap,
            best_size,
            near_ground_deficiency,
            near_ground_terminal_cap,
            near_ground_layer_cap,
            progress_every,
        )
    };
    let deficiency_layers: Vec<DeficiencyLayer> = (0..=near_ground_deficiency)
        .filter_map(|deficiency| {
            best_size
                .checked_sub(deficiency)
                .map(|cardinality| DeficiencyLayer {
                    deficiency,
                    cardinality,
                    count: layer_counts[cardinality],
                })
        })
        .collect();
    let first_excited = deficiency_layers
        .iter()
        .find(|layer| layer.deficiency > 0 && layer.count > 0);

    Row {
        n,
        g,
        ordered_representation_cap: cap,
        h_n_exact_finite: best_size,
        maximizer_count,
        first_maximizer_witness,
        state_model: StateModel {
            representation: "occupied prefix plus bounded ordered-sum-count memory".to_string(),
            occupied_lattice: "[0,n]".to_string(),
            memory: "sum_counts[s] stores ordered representation count for a+b=s".to_string(),
            transition_rule:
                "occupy x iff all new counts remain <= 2g; self pair adds 1 and each occupied a adds ordered pairs (a,x),(x,a)"
                    .to_string(),
        },
        transfer_diagnostics: TransferDiagnostics {
            traversal_mode: if imported_ground.is_some() {
                "b2g_imported_ground_nearground_transfer_state".to_string()
            } else {
                "b2g_exact_branch_and_bound_transfer_state".to_string()
            },
            sites_processed: n + 1,
            nodes_visited,
            transitions_attempted,
            legal_occupy_transitions,
            rejected_occupy_transitions,
            branch_bound_pruned_nodes,
            best_updates,
            progress_every_nodes: progress_every,
        },
        spectral_observables: SpectralObservables {
            ground_cardinality_h_n: best_size,
            ground_state_degeneracy: maximizer_count,
            ground_entropy_ln: (maximizer_count as f64).ln(),
            first_excited_cardinality: first_excited.map(|layer| layer.cardinality),
            first_excited_count: first_excited.map(|layer| layer.count),
            gap_to_first_excited_layer: first_excited.map(|layer| layer.deficiency),
            near_ground_count_h_minus_1: best_size
                .checked_sub(1)
                .map(|cardinality| layer_counts[cardinality])
                .unwrap_or(0),
            near_ground_count_h_minus_2: best_size
                .checked_sub(2)
                .map(|cardinality| layer_counts[cardinality])
                .unwrap_or(0),
            deficiency_layers,
            near_ground_diagnostics,
        },
        top_k_frontier: TopKFrontier {
            k: frontier_k,
            zero_temperature_field_interpretation:
                "prefix_best, density_adjusted_mass_best, and joint_best are zero-temperature selections after cardinality is fixed to the exact finite B_2[g] maximum"
                    .to_string(),
            prefix_best: prefix_frontier.ranked(),
            density_adjusted_mass_best: mass_frontier.ranked(),
            joint_best: joint_frontier.ranked(),
        },
    }
}

fn summarize(rows: &[Row]) -> Summary {
    let mut split_ns = Vec::new();
    for row in rows {
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
        min_h_n: rows
            .iter()
            .map(|row| row.h_n_exact_finite)
            .min()
            .unwrap_or(0),
        max_h_n: rows
            .iter()
            .map(|row| row.h_n_exact_finite)
            .max()
            .unwrap_or(0),
        split_count: split_ns.len(),
        split_ns,
        total_maximizers: rows.iter().map(|row| row.maximizer_count).sum(),
    }
}

fn build_report(results: &Results) -> String {
    let mut lines = Vec::new();
    lines.push(format!(
        "# {} — B_2[{}] Transfer-State Scout",
        results.experiment_id, results.g
    ));
    lines.push(String::new());
    lines.push("## Question".to_string());
    lines.push(String::new());
    lines.push(
        "Does the Sidon-adjacent bounded-sum deformation show the same field-sensitive ground-face behavior seen in #30?".to_string(),
    );
    lines.push(String::new());
    lines.push("## Result".to_string());
    lines.push(String::new());
    lines.push(format!(
        "Exact finite maxima were computed for {} rows. The exact h(n) range was {}..{}. Field splits occurred in {} rows: {:?}.",
        results.summary.checked_count,
        results.summary.min_h_n,
        results.summary.max_h_n,
        results.summary.split_count,
        results.summary.split_ns
    ));
    lines.push(String::new());
    lines.push(
        "| n | h(n) | maximizers | h-1 | h-2 | nodes | pruned | prefix winner | mass winner | joint winner |"
            .to_string(),
    );
    lines.push("|---|---:|---:|---:|---:|---:|---:|---|---|---|".to_string());
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
            "| {} | {} | {} | {} | {} | {} | {} | `{}` | `{}` | `{}` |",
            row.n,
            row.h_n_exact_finite,
            row.maximizer_count,
            row.spectral_observables.near_ground_count_h_minus_1,
            row.spectral_observables.near_ground_count_h_minus_2,
            row.transfer_diagnostics.nodes_visited,
            row.transfer_diagnostics.branch_bound_pruned_nodes,
            prefix,
            mass,
            joint
        ));
    }
    lines.push(String::new());
    lines.push("## Claim Boundary".to_string());
    lines.push(String::new());
    lines.push(results.claim_boundary.clone());
    if results.rows.iter().any(|row| {
        !row.spectral_observables
            .near_ground_diagnostics
            .counts_exact
    }) {
        lines.push(String::new());
        lines.push("## Capped Near-Ground Counts".to_string());
        lines.push(String::new());
        lines.push(
            "At least one row hit `near_ground_terminal_cap`; near-ground layer counts in this packet are lower bounds, not exact counts."
                .to_string(),
        );
    }
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
    if cli.n_max > 80 {
        return Err("this exact scout is intentionally capped at --n-max <= 80".into());
    }
    if cli.g == 0 {
        return Err("--g must be positive".into());
    }
    let experiment_id = cli.experiment_id.unwrap_or_else(|| {
        format!(
            "EXP-MM-755-PMF-B2G{}-TRANSFER-SCOUT-{}-{}-2026-04-29",
            cli.g, cli.n_min, cli.n_max
        )
    });
    let t0 = Instant::now();
    let rows: Vec<Row> = (cli.n_min..=cli.n_max)
        .map(|n| {
            analyze_n(
                n,
                cli.g,
                cli.frontier_k,
                cli.near_ground_deficiency,
                cli.near_ground_terminal_cap,
                cli.near_ground_layer_cap,
                cli.skip_near_ground,
                cli.known_h,
                cli.known_ground_count,
                cli.progress_every,
            )
        })
        .collect();
    let summary = summarize(&rows);
    let results = Results {
        experiment_id,
        date: "2026-04-29".to_string(),
        erdos_problem: 755,
        scan_mode: "pmf_b2g_transfer_state_exact_scout".to_string(),
        implementation: "Rust exact branch-and-bound transfer-state scanner".to_string(),
        constraint: "B_2[g] on [0,n] with ordered representation count <= 2g".to_string(),
        claim_boundary:
            "This is an exact finite computational scout for the B_2[g] transfer state, not an independent theorem or parity check against a prior #755 reference packet."
                .to_string(),
        n_min: cli.n_min,
        n_max: cli.n_max,
        g: cli.g,
        frontier_k: cli.frontier_k,
        near_ground_deficiency: cli.near_ground_deficiency,
        near_ground_terminal_cap: cli.near_ground_terminal_cap,
        near_ground_layer_cap: cli.near_ground_layer_cap,
        skip_near_ground: cli.skip_near_ground,
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
