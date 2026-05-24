use clap::Parser;
use rayon::prelude::*;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(name = "erdos30-exact-maximizer")]
#[command(about = "Exact maximizer Sidon observable scan for Erdős #30")]
struct Cli {
    #[arg(long, default_value_t = 10)]
    n_min: usize,

    #[arg(long, default_value_t = 50)]
    n_max: usize,

    #[arg(long)]
    experiment_id: Option<String>,

    #[arg(long, default_value = "erdos-experiments/results/erdos-30")]
    output_dir: PathBuf,

    #[arg(long)]
    progress: bool,

    #[arg(long, default_value_t = 0)]
    seed_depth: usize,

    #[arg(long)]
    initial_lower_bound: Option<usize>,

    #[arg(long, default_value_t = 0)]
    parallel_prefix_depth: usize,

    #[arg(long, default_value_t = 5)]
    frontier_k: usize,

    #[arg(long)]
    inheritance_source_results: Option<PathBuf>,

    #[arg(long, default_value_t = false)]
    export_ground_face_masks: bool,
}

#[derive(Clone, Serialize)]
struct MetricSummary {
    min_value: f64,
    mean_value: f64,
    max_value: f64,
    min_set: Vec<usize>,
    max_set: Vec<usize>,
    min_location: Option<usize>,
    max_location: Option<usize>,
}

#[derive(Clone)]
struct MetricAccumulator {
    count: usize,
    total: f64,
    min_value: f64,
    max_value: f64,
    min_set: Vec<usize>,
    max_set: Vec<usize>,
    min_location: Option<usize>,
    max_location: Option<usize>,
}

impl MetricAccumulator {
    fn new() -> Self {
        Self {
            count: 0,
            total: 0.0,
            min_value: f64::INFINITY,
            max_value: f64::NEG_INFINITY,
            min_set: Vec::new(),
            max_set: Vec::new(),
            min_location: None,
            max_location: None,
        }
    }

    fn reset(&mut self) {
        *self = Self::new();
    }

    fn observe(&mut self, value: f64, chosen: &[usize], location: Option<usize>) {
        self.count += 1;
        self.total += value;
        if value < self.min_value
            || (value == self.min_value && (self.min_set.is_empty() || chosen < self.min_set.as_slice()))
        {
            self.min_value = value;
            self.min_set = chosen.to_vec();
            self.min_location = location;
        }
        if value > self.max_value
            || (value == self.max_value && (self.max_set.is_empty() || chosen < self.max_set.as_slice()))
        {
            self.max_value = value;
            self.max_set = chosen.to_vec();
            self.max_location = location;
        }
    }

    fn merge(&mut self, other: &Self) {
        if other.count == 0 {
            return;
        }
        self.count += other.count;
        self.total += other.total;
        if other.min_value < self.min_value
            || (other.min_value == self.min_value
                && (self.min_set.is_empty() || other.min_set < self.min_set))
        {
            self.min_value = other.min_value;
            self.min_set = other.min_set.clone();
            self.min_location = other.min_location;
        }
        if other.max_value > self.max_value
            || (other.max_value == self.max_value
                && (self.max_set.is_empty() || other.max_set < self.max_set))
        {
            self.max_value = other.max_value;
            self.max_set = other.max_set.clone();
            self.max_location = other.max_location;
        }
    }

    fn summary(&self) -> MetricSummary {
        MetricSummary {
            min_value: self.min_value,
            mean_value: self.total / self.count as f64,
            max_value: self.max_value,
            min_set: self.min_set.clone(),
            max_set: self.max_set.clone(),
            min_location: self.min_location,
            max_location: self.max_location,
        }
    }
}

#[derive(Serialize)]
struct Normalized7 {
    min: f64,
    mean: f64,
    max: f64,
}

#[derive(Serialize)]
struct Normalized11 {
    min: f64,
    mean: f64,
    max: f64,
}

#[derive(Serialize)]
struct Metric7 {
    raw: MetricSummary,
    normalized_by_n_pow_7_8: Normalized7,
}

#[derive(Serialize)]
struct Metric11 {
    raw: MetricSummary,
    normalized_by_n_pow_11_8: Normalized11,
}

#[derive(Serialize)]
struct JointMetric {
    raw: MetricSummary,
}

#[derive(Serialize)]
struct Normalizers {
    sqrt_n: f64,
    n_pow_7_8: f64,
    n_pow_11_8: f64,
}

#[derive(Serialize)]
struct ObservableSplit {
    same_best_witness: bool,
    prefix_best_witness: Vec<usize>,
    mass_best_witness: Vec<usize>,
    prefix_best_mass_dev: f64,
    prefix_best_mass_ratio_n_pow_11_8: f64,
    mass_best_prefix_residual: f64,
    mass_best_prefix_ratio_n_pow_7_8: f64,
    best_joint_witness: Vec<usize>,
    best_joint_score: f64,
    worst_joint_witness: Vec<usize>,
    worst_joint_score: f64,
}

#[derive(Clone, Serialize)]
struct SearchDiagnostics {
    seed_depth: usize,
    initial_lower_bound_size: usize,
    nodes_visited: u64,
    size_bound_prunes: u64,
    sidon_extension_rejects: u64,
    terminal_candidates_seen: u64,
}

impl SearchDiagnostics {
    fn new(seed_depth: usize, initial_lower_bound_size: usize) -> Self {
        Self {
            seed_depth,
            initial_lower_bound_size,
            nodes_visited: 0,
            size_bound_prunes: 0,
            sidon_extension_rejects: 0,
            terminal_candidates_seen: 0,
        }
    }

    fn merge(&mut self, other: &Self) {
        self.nodes_visited += other.nodes_visited;
        self.size_bound_prunes += other.size_bound_prunes;
        self.sidon_extension_rejects += other.sidon_extension_rejects;
        self.terminal_candidates_seen += other.terminal_candidates_seen;
    }
}

#[derive(Clone, Serialize)]
struct FrontierWitness {
    rank: usize,
    witness: Vec<usize>,
    objective_value: f64,
    prefix_residual: f64,
    prefix_ratio_n_pow_7_8: f64,
    prefix_location: Option<usize>,
    density_adjusted_mass_dev: f64,
    density_adjusted_mass_ratio_n_pow_11_8: f64,
    joint_score: f64,
}

#[derive(Clone)]
struct FrontierCandidate {
    witness: Vec<usize>,
    objective_value: f64,
    prefix_residual: f64,
    prefix_ratio_n_pow_7_8: f64,
    prefix_location: Option<usize>,
    density_adjusted_mass_dev: f64,
    density_adjusted_mass_ratio_n_pow_11_8: f64,
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

    fn reset(&mut self) {
        self.entries.clear();
    }

    fn observe(&mut self, candidate: FrontierCandidate) {
        if self.k == 0 {
            return;
        }
        self.entries.push(candidate);
        self.entries.sort_by(|a, b| {
            a.objective_value
                .total_cmp(&b.objective_value)
                .then_with(|| a.witness.cmp(&b.witness))
        });
        self.entries.truncate(self.k);
    }

    fn merge(&mut self, other: &Self) {
        for candidate in &other.entries {
            self.observe(candidate.clone());
        }
    }

    fn ranked(&self) -> Vec<FrontierWitness> {
        self.entries
            .iter()
            .enumerate()
            .map(|(idx, item)| FrontierWitness {
                rank: idx + 1,
                witness: item.witness.clone(),
                objective_value: item.objective_value,
                prefix_residual: item.prefix_residual,
                prefix_ratio_n_pow_7_8: item.prefix_ratio_n_pow_7_8,
                prefix_location: item.prefix_location,
                density_adjusted_mass_dev: item.density_adjusted_mass_dev,
                density_adjusted_mass_ratio_n_pow_11_8: item.density_adjusted_mass_ratio_n_pow_11_8,
                joint_score: item.joint_score,
            })
            .collect()
    }
}

#[derive(Serialize)]
struct TopKFrontier {
    k: usize,
    prefix_best: Vec<FrontierWitness>,
    density_adjusted_mass_best: Vec<FrontierWitness>,
    joint_best: Vec<FrontierWitness>,
}

#[derive(Clone)]
struct InheritanceSource {
    results_path: String,
    source_n: usize,
    previous_face_count: usize,
    previous_masks: HashSet<u128>,
    shifted_masks: HashSet<u128>,
    inherited_union_count: usize,
}

#[derive(Clone, Default)]
struct InheritanceProbeAccum {
    exact_previous_persistence_count: u64,
    plus_one_previous_persistence_count: u64,
    inherited_union_present_count: u64,
}

#[derive(Serialize)]
struct InheritanceProbe {
    source_results_path: String,
    source_n: usize,
    previous_face_count: usize,
    previous_plus_one_count: usize,
    inherited_union_count: usize,
    exact_previous_persistence_count: u64,
    plus_one_previous_persistence_count: u64,
    inherited_union_present_count: u64,
    new_face_count: u64,
    contains_previous_face: bool,
    contains_previous_face_plus_one: bool,
    contains_inherited_union: bool,
}

#[derive(Serialize)]
struct GroundFaceMaskExport {
    status: String,
    exact_maximizer_count: u64,
    exported_count: usize,
    encoding: String,
    masks: Vec<String>,
}

struct ScanAccumulator {
    best_size: usize,
    best_count: u64,
    profile_acc: MetricAccumulator,
    prefix_acc: MetricAccumulator,
    prefix_residual_acc: MetricAccumulator,
    recentered_prefix_acc: MetricAccumulator,
    density_adjusted_prefix_acc: MetricAccumulator,
    mass_acc: MetricAccumulator,
    density_adjusted_mass_acc: MetricAccumulator,
    joint_acc: MetricAccumulator,
    prefix_frontier: TopKAccumulator,
    mass_frontier: TopKAccumulator,
    joint_frontier: TopKAccumulator,
    inheritance_probe: InheritanceProbeAccum,
    ground_face_masks: Vec<u128>,
    best_prefix_mass_dev: f64,
    best_prefix_mass_ratio: f64,
    best_mass_prefix_residual: f64,
    best_mass_prefix_ratio: f64,
    search_diagnostics: SearchDiagnostics,
}

impl ScanAccumulator {
    fn new(seed_depth: usize, initial_lower_bound: usize, frontier_k: usize) -> Self {
        Self {
            best_size: initial_lower_bound,
            best_count: 0,
            profile_acc: MetricAccumulator::new(),
            prefix_acc: MetricAccumulator::new(),
            prefix_residual_acc: MetricAccumulator::new(),
            recentered_prefix_acc: MetricAccumulator::new(),
            density_adjusted_prefix_acc: MetricAccumulator::new(),
            mass_acc: MetricAccumulator::new(),
            density_adjusted_mass_acc: MetricAccumulator::new(),
            joint_acc: MetricAccumulator::new(),
            prefix_frontier: TopKAccumulator::new(frontier_k),
            mass_frontier: TopKAccumulator::new(frontier_k),
            joint_frontier: TopKAccumulator::new(frontier_k),
            inheritance_probe: InheritanceProbeAccum::default(),
            ground_face_masks: Vec::new(),
            best_prefix_mass_dev: f64::INFINITY,
            best_prefix_mass_ratio: f64::INFINITY,
            best_mass_prefix_residual: f64::INFINITY,
            best_mass_prefix_ratio: f64::INFINITY,
            search_diagnostics: SearchDiagnostics::new(seed_depth, initial_lower_bound),
        }
    }

    fn reset_for_new_best(&mut self, size: usize) {
        self.best_size = size;
        self.best_count = 0;
        self.profile_acc.reset();
        self.prefix_acc.reset();
        self.prefix_residual_acc.reset();
        self.recentered_prefix_acc.reset();
        self.density_adjusted_prefix_acc.reset();
        self.mass_acc.reset();
        self.density_adjusted_mass_acc.reset();
        self.joint_acc.reset();
        self.prefix_frontier.reset();
        self.mass_frontier.reset();
        self.joint_frontier.reset();
        self.inheritance_probe = InheritanceProbeAccum::default();
        self.ground_face_masks.clear();
        self.best_prefix_mass_dev = f64::INFINITY;
        self.best_prefix_mass_ratio = f64::INFINITY;
        self.best_mass_prefix_residual = f64::INFINITY;
        self.best_mass_prefix_ratio = f64::INFINITY;
    }

    fn observe_terminal(
        &mut self,
        n: usize,
        chosen: &[usize],
        sqrt_n: f64,
        n_pow_7_8: f64,
        n_pow_11_8: f64,
        inheritance_source: Option<&InheritanceSource>,
        export_ground_face_masks: bool,
    ) {
        let size = chosen.len();
        if size == 0 || size < self.best_size {
            return;
        }
        self.search_diagnostics.terminal_candidates_seen += 1;
        if size > self.best_size {
            self.reset_for_new_best(size);
        }
        self.best_count += 1;
        let mask = sites_to_mask(chosen).expect("chosen site exceeds bitset width");
        if let Some(source) = inheritance_source {
            observe_inheritance(mask, source, &mut self.inheritance_probe);
        }
        if export_ground_face_masks {
            self.ground_face_masks.push(mask);
        }

        let explicit_center = compute_explicit_mass_center(size, sqrt_n);
        let (profile_dev, profile_i) = compute_profile_deviation(chosen, sqrt_n);
        let (prefix_dev, prefix_t) = compute_prefix_deviation(chosen, n, sqrt_n);
        let gap = (size as f64 - sqrt_n).abs();
        let prefix_drift = gap.max(1.0) * sqrt_n;
        let endpoint_bridge = n as f64 - size as f64 * sqrt_n;
        let density_slope = n as f64 / size as f64;
        let (recentered_prefix_dev, recentered_prefix_t) =
            compute_recentered_prefix_deviation(chosen, n, sqrt_n, endpoint_bridge);
        let (density_adjusted_prefix_dev, density_adjusted_prefix_t) =
            compute_prefix_deviation(chosen, n, density_slope);
        let set_sum: usize = chosen.iter().sum();
        let mass_dev = (set_sum as f64 - explicit_center).abs();
        let density_adjusted_mass_center = compute_density_adjusted_mass_center(size, n);
        let density_adjusted_mass_dev = (set_sum as f64 - density_adjusted_mass_center).abs();
        let prefix_residual = (prefix_dev - prefix_drift).max(0.0);
        let joint_score = compute_joint_observable_score(prefix_residual, density_adjusted_mass_dev, n);

        let improves_prefix_best = prefix_residual < self.prefix_residual_acc.min_value;
        let improves_mass_best =
            density_adjusted_mass_dev < self.density_adjusted_mass_acc.min_value;

        self.profile_acc.observe(profile_dev, chosen, profile_i);
        self.prefix_acc.observe(prefix_dev, chosen, prefix_t);
        self.prefix_residual_acc.observe(prefix_residual, chosen, prefix_t);
        self.recentered_prefix_acc
            .observe(recentered_prefix_dev, chosen, recentered_prefix_t);
        self.density_adjusted_prefix_acc.observe(
            density_adjusted_prefix_dev,
            chosen,
            density_adjusted_prefix_t,
        );
        self.mass_acc.observe(mass_dev, chosen, None);
        self.density_adjusted_mass_acc
            .observe(density_adjusted_mass_dev, chosen, None);
        self.joint_acc.observe(joint_score, chosen, None);
        let frontier_candidate = FrontierCandidate {
            witness: chosen.to_vec(),
            objective_value: prefix_residual,
            prefix_residual,
            prefix_ratio_n_pow_7_8: prefix_residual / n_pow_7_8,
            prefix_location: prefix_t,
            density_adjusted_mass_dev,
            density_adjusted_mass_ratio_n_pow_11_8: density_adjusted_mass_dev / n_pow_11_8,
            joint_score,
        };
        self.prefix_frontier.observe(frontier_candidate.clone());
        self.mass_frontier.observe(FrontierCandidate {
            objective_value: density_adjusted_mass_dev,
            ..frontier_candidate.clone()
        });
        self.joint_frontier.observe(FrontierCandidate {
            objective_value: joint_score,
            ..frontier_candidate
        });
        if improves_prefix_best {
            self.best_prefix_mass_dev = density_adjusted_mass_dev;
            self.best_prefix_mass_ratio = density_adjusted_mass_dev / n_pow_11_8;
        }
        if improves_mass_best {
            self.best_mass_prefix_residual = prefix_residual;
            self.best_mass_prefix_ratio = prefix_residual / n_pow_7_8;
        }
    }

    fn merge_worker(&mut self, other: Self) {
        let mut diagnostics = self.search_diagnostics.clone();
        diagnostics.merge(&other.search_diagnostics);
        if other.best_size > self.best_size {
            *self = other;
            self.search_diagnostics = diagnostics;
            return;
        }
        if other.best_size == self.best_size {
            self.best_count += other.best_count;
            self.profile_acc.merge(&other.profile_acc);
            self.prefix_acc.merge(&other.prefix_acc);
            self.prefix_residual_acc.merge(&other.prefix_residual_acc);
            self.recentered_prefix_acc.merge(&other.recentered_prefix_acc);
            self.density_adjusted_prefix_acc
                .merge(&other.density_adjusted_prefix_acc);
            self.mass_acc.merge(&other.mass_acc);
            self.density_adjusted_mass_acc
                .merge(&other.density_adjusted_mass_acc);
            self.joint_acc.merge(&other.joint_acc);
            self.prefix_frontier.merge(&other.prefix_frontier);
            self.mass_frontier.merge(&other.mass_frontier);
            self.joint_frontier.merge(&other.joint_frontier);
            self.inheritance_probe.exact_previous_persistence_count +=
                other.inheritance_probe.exact_previous_persistence_count;
            self.inheritance_probe.plus_one_previous_persistence_count +=
                other.inheritance_probe.plus_one_previous_persistence_count;
            self.inheritance_probe.inherited_union_present_count +=
                other.inheritance_probe.inherited_union_present_count;
            self.ground_face_masks.extend(other.ground_face_masks);
            if other.best_prefix_mass_dev < self.best_prefix_mass_dev {
                self.best_prefix_mass_dev = other.best_prefix_mass_dev;
                self.best_prefix_mass_ratio = other.best_prefix_mass_ratio;
            }
            if other.best_mass_prefix_residual < self.best_mass_prefix_residual {
                self.best_mass_prefix_residual = other.best_mass_prefix_residual;
                self.best_mass_prefix_ratio = other.best_mass_prefix_ratio;
            }
        }
        self.search_diagnostics = diagnostics;
    }
}

#[derive(Serialize)]
struct Row {
    n: usize,
    maximizer_size_h_n: usize,
    maximizer_set_count: u64,
    maximizer_scan_runtime_sec: f64,
    real_gap_from_sqrt: f64,
    real_deficiency_from_sqrt: f64,
    general_prefix_drift: f64,
    endpoint_bridge: f64,
    density_adjusted_slope: f64,
    explicit_mass_center: f64,
    density_adjusted_mass_center: f64,
    normalizers: Normalizers,
    ordered_profile_deviation: Metric7,
    prefix_deviation: Metric7,
    prefix_residual_after_general_drift: Metric7,
    recentered_prefix_deviation: Metric7,
    density_adjusted_prefix_deviation: Metric7,
    mass_deviation: Metric11,
    density_adjusted_mass_deviation: Metric11,
    joint_observable_score: JointMetric,
    observable_split: ObservableSplit,
    search_diagnostics: SearchDiagnostics,
    top_k_frontier: TopKFrontier,
    inheritance_probe: Option<InheritanceProbe>,
    ground_face_mask_export: Option<GroundFaceMaskExport>,
}

#[derive(Serialize)]
struct PrefixDriftDiagnostic {
    label: String,
    worst_n: usize,
    worst_residual: f64,
    worst_residual_ratio_n_pow_7_8: f64,
    worst_raw_prefix_deviation: f64,
    worst_drift: f64,
    worst_set: Vec<usize>,
    worst_t: Option<usize>,
    worst_endpoint_bridge: f64,
}

#[derive(Serialize)]
struct PrefixDriftDiagnostics {
    raw_prefix: PrefixDriftDiagnostic,
    sqrt_only: PrefixDriftDiagnostic,
    gap_sqrt_only: PrefixDriftDiagnostic,
    current_theorem_drift: PrefixDriftDiagnostic,
    affine_endpoint_recentered: PrefixDriftDiagnostic,
    density_adjusted_affine: PrefixDriftDiagnostic,
}

#[derive(Serialize)]
struct MassCenterDiagnostic {
    label: String,
    worst_n: usize,
    worst_residual: f64,
    worst_residual_ratio_n_pow_11_8: f64,
    worst_set: Vec<usize>,
}

#[derive(Serialize)]
struct MassCenterDiagnostics {
    sqrt_center: MassCenterDiagnostic,
    density_adjusted_center: MassCenterDiagnostic,
}

#[derive(Serialize)]
struct FiniteCompatibility {
    zero_tolerance: f64,
    mass_best_near_zero_prefix_count: usize,
    mass_best_near_zero_prefix_ns: Vec<usize>,
    prefix_best_zero_mass_count: usize,
    prefix_best_zero_mass_ns: Vec<usize>,
    mass_best_prefix_ratio_max_n_pow_7_8: f64,
    mass_best_prefix_ratio_mean_n_pow_7_8: f64,
    prefix_best_mass_ratio_max_n_pow_11_8: f64,
    prefix_best_mass_ratio_mean_n_pow_11_8: f64,
    mass_penalty_le_prefix_penalty_count: usize,
    mass_penalty_le_prefix_penalty_ns: Vec<usize>,
    mass_penalty_gt_prefix_penalty_count: usize,
    mass_penalty_gt_prefix_penalty_ns: Vec<usize>,
    joint_is_mass_best_count: usize,
    joint_is_mass_best_ns: Vec<usize>,
    joint_is_prefix_best_count: usize,
    joint_is_prefix_best_ns: Vec<usize>,
    joint_is_third_count: usize,
    joint_is_third_ns: Vec<usize>,
}

#[derive(Serialize)]
struct ObservableSplitDiagnostics {
    same_best_count: usize,
    same_best_ns: Vec<usize>,
    strongest_split_n: usize,
    strongest_split_score: f64,
    strongest_split_prefix_best_witness: Vec<usize>,
    strongest_split_mass_best_witness: Vec<usize>,
    strongest_split_prefix_best_mass_ratio_n_pow_11_8: f64,
    strongest_split_mass_best_prefix_ratio_n_pow_7_8: f64,
    best_joint_n: usize,
    best_joint_witness: Vec<usize>,
    best_joint_score: f64,
    finite_compatibility: FiniteCompatibility,
}

#[derive(Serialize)]
struct FirstHitVsFaceAwareDiagnostics {
    zero_tolerance: f64,
    first_hit_compatible_count: usize,
    first_hit_compatible_ns: Vec<usize>,
    first_hit_failure_count: usize,
    first_hit_failure_ns: Vec<usize>,
    face_aware_near_zero_joint_count: usize,
    face_aware_near_zero_joint_ns: Vec<usize>,
    face_aware_nonzero_joint_count: usize,
    face_aware_nonzero_joint_ns: Vec<usize>,
    first_hit_failures_recovered_by_face_count: usize,
    first_hit_failures_recovered_by_face_ns: Vec<usize>,
    first_hit_failures_persisting_as_nonzero_face_count: usize,
    first_hit_failures_persisting_as_nonzero_face_ns: Vec<usize>,
    max_face_aware_joint_n: Option<usize>,
    max_face_aware_joint_score: Option<f64>,
    max_face_aware_joint_prefix_ratio_n_pow_7_8: Option<f64>,
    max_face_aware_joint_mass_ratio_n_pow_11_8: Option<f64>,
    max_face_aware_joint_witness: Option<Vec<usize>>,
}

#[derive(Serialize)]
struct Results {
    experiment_id: String,
    date: String,
    description: String,
    erdos_problem: usize,
    scan_mode: String,
    implementation: String,
    n_min: usize,
    n_max: usize,
    inheritance_source_results: Option<String>,
    export_ground_face_masks: bool,
    #[serde(rename = "total_runtime_sec")]
    runtime_sec: f64,
    results_by_n: Vec<Row>,
    prefix_drift_diagnostics: PrefixDriftDiagnostics,
    mass_center_diagnostics: MassCenterDiagnostics,
    observable_split_diagnostics: ObservableSplitDiagnostics,
    first_hit_vs_face_aware_diagnostics: FirstHitVsFaceAwareDiagnostics,
}

fn extension_sums(
    chosen: &[usize],
    used_sums: &[bool],
    x: usize,
    scratch: &mut Vec<usize>,
) -> bool {
    scratch.clear();
    let self_sum = 2 * x;
    if used_sums[self_sum] {
        return false;
    }
    scratch.push(self_sum);
    for &a in chosen {
        let pair_sum = a + x;
        if used_sums[pair_sum] {
            return false;
        }
        scratch.push(pair_sum);
    }
    true
}

fn bit(index: usize) -> Option<u128> {
    if index < 128 {
        Some(1u128 << index)
    } else {
        None
    }
}

fn sites_to_mask(sites: &[usize]) -> Option<u128> {
    let mut mask = 0u128;
    for &site in sites {
        mask |= bit(site)?;
    }
    Some(mask)
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

fn shifted_mask_from_mask(mask: u128, shift: usize) -> Option<u128> {
    let shifted_sites: Vec<usize> = bits_to_sites(mask).into_iter().map(|site| site + shift).collect();
    sites_to_mask(&shifted_sites)
}

fn load_row_from_results(value: &serde_json::Value) -> Result<&serde_json::Value, Box<dyn std::error::Error>> {
    if let Some(rows) = value.get("rows").and_then(|rows| rows.as_array()) {
        return rows.last().ok_or_else(|| "source rows array is empty".into());
    }
    if let Some(rows) = value.get("results_by_n").and_then(|rows| rows.as_array()) {
        return rows
            .last()
            .ok_or_else(|| "source results_by_n array is empty".into());
    }
    Err("inheritance source has neither rows nor results_by_n array".into())
}

fn load_inheritance_source(path: &Path) -> Result<InheritanceSource, Box<dyn std::error::Error>> {
    let payload: serde_json::Value = serde_json::from_slice(&fs::read(path)?)?;
    let row = load_row_from_results(&payload)?;
    let source_n = row
        .get("n")
        .and_then(|value| value.as_u64())
        .ok_or("inheritance source row has no numeric n")? as usize;

    let mut previous_masks = HashSet::new();
    let mut shifted_masks = HashSet::new();
    if let Some(masks) = row
        .get("ground_face_mask_export")
        .and_then(|value| value.get("masks"))
        .and_then(|value| value.as_array())
    {
        if masks.is_empty() {
            return Err("inheritance source mask array is empty".into());
        }
        for mask_value in masks {
            let mask = mask_value
                .as_str()
                .ok_or("inheritance source mask is not a string")?
                .parse::<u128>()?;
            previous_masks.insert(mask);
            shifted_masks.insert(
                shifted_mask_from_mask(mask, 1).ok_or("shifted inheritance mask exceeds bitset width")?,
            );
        }
    } else {
        let export = row
            .get("ground_face_export")
            .ok_or("inheritance source row has neither mask export nor ground_face_export")?;
        let witnesses = export
            .get("witnesses")
            .and_then(|value| value.as_array())
            .ok_or("inheritance source ground_face_export has no witnesses array")?;
        if witnesses.is_empty() {
            return Err("inheritance source witness array is empty".into());
        }
        for witness in witnesses {
            let sites: Vec<usize> = witness
                .get("witness")
                .and_then(|value| value.as_array())
                .ok_or("inheritance source witness row has no witness array")?
                .iter()
                .map(|value| {
                    value
                        .as_u64()
                        .map(|site| site as usize)
                        .ok_or("inheritance witness site is not numeric")
                })
                .collect::<Result<Vec<_>, _>>()?;
            let mask = sites_to_mask(&sites).ok_or("inheritance witness site exceeds bitset width")?;
            previous_masks.insert(mask);
            shifted_masks
                .insert(shifted_mask_from_mask(mask, 1).ok_or("shifted inheritance witness exceeds bitset width")?);
        }
    }
    let inherited_union_count = previous_masks.union(&shifted_masks).count();
    Ok(InheritanceSource {
        results_path: path.canonicalize()?.display().to_string(),
        source_n,
        previous_face_count: previous_masks.len(),
        previous_masks,
        shifted_masks,
        inherited_union_count,
    })
}

fn observe_inheritance(mask: u128, source: &InheritanceSource, accum: &mut InheritanceProbeAccum) {
    let in_previous = source.previous_masks.contains(&mask);
    let in_shifted = source.shifted_masks.contains(&mask);
    if in_previous {
        accum.exact_previous_persistence_count += 1;
    }
    if in_shifted {
        accum.plus_one_previous_persistence_count += 1;
    }
    if in_previous || in_shifted {
        accum.inherited_union_present_count += 1;
    }
}

fn build_inheritance_probe(
    source: Option<&InheritanceSource>,
    accum: &InheritanceProbeAccum,
    exact_maximizer_count: u64,
) -> Option<InheritanceProbe> {
    let source = source?;
    Some(InheritanceProbe {
        source_results_path: source.results_path.clone(),
        source_n: source.source_n,
        previous_face_count: source.previous_face_count,
        previous_plus_one_count: source.shifted_masks.len(),
        inherited_union_count: source.inherited_union_count,
        exact_previous_persistence_count: accum.exact_previous_persistence_count,
        plus_one_previous_persistence_count: accum.plus_one_previous_persistence_count,
        inherited_union_present_count: accum.inherited_union_present_count,
        new_face_count: exact_maximizer_count.saturating_sub(accum.inherited_union_present_count),
        contains_previous_face: accum.exact_previous_persistence_count as usize == source.previous_face_count,
        contains_previous_face_plus_one: accum.plus_one_previous_persistence_count as usize
            == source.shifted_masks.len(),
        contains_inherited_union: accum.inherited_union_present_count as usize
            == source.inherited_union_count,
    })
}

fn build_mask_export(export: bool, exact_maximizer_count: u64, mut masks: Vec<u128>) -> Option<GroundFaceMaskExport> {
    if !export {
        return None;
    }
    masks.sort_unstable();
    Some(GroundFaceMaskExport {
        status: "EXPORTED_ALL".to_string(),
        exact_maximizer_count,
        exported_count: masks.len(),
        encoding: "u128 occupied_mask decimal string; bit i is 1 iff site i is occupied".to_string(),
        masks: masks.into_iter().map(|mask| mask.to_string()).collect(),
    })
}

fn greedy_sidon_size_with_prefix(n: usize, prefix: &[usize], reverse: bool) -> usize {
    let mut chosen = Vec::new();
    let mut used_sums = vec![false; 2 * n + 1];
    let mut scratch = Vec::new();
    for &x in prefix {
        if x > n || !extension_sums(&chosen, &used_sums, x, &mut scratch) {
            return 0;
        }
        chosen.push(x);
        for &value in &scratch {
            used_sums[value] = true;
        }
    }
    let iter: Box<dyn Iterator<Item = usize>> = if reverse {
        Box::new((0..=n).rev())
    } else {
        Box::new(0..=n)
    };
    for x in iter {
        if !extension_sums(&chosen, &used_sums, x, &mut scratch) {
            continue;
        }
        chosen.push(x);
        for &value in &scratch {
            used_sums[value] = true;
        }
    }
    chosen.len()
}

fn initial_lower_bound_size(n: usize, seed_depth: usize) -> usize {
    let mut best = 0usize;
    for reverse in [false, true] {
        best = best.max(greedy_sidon_size_with_prefix(n, &[], reverse));
        if seed_depth < 1 {
            continue;
        }
        for a in 0..=n {
            best = best.max(greedy_sidon_size_with_prefix(n, &[a], reverse));
        }
        if seed_depth < 2 {
            continue;
        }
        for a in 0..=n {
            for b in (a + 1)..=n {
                best = best.max(greedy_sidon_size_with_prefix(n, &[a, b], reverse));
            }
        }
    }
    best
}

#[derive(Clone)]
struct PrefixJob {
    start: usize,
    chosen: Vec<usize>,
    used_sums: Vec<bool>,
}

fn generate_prefix_jobs(n: usize, depth: usize) -> Vec<PrefixJob> {
    fn rec(
        n: usize,
        start: usize,
        depth_remaining: usize,
        chosen: &mut Vec<usize>,
        used_sums: &mut Vec<bool>,
        scratch: &mut Vec<usize>,
        jobs: &mut Vec<PrefixJob>,
    ) {
        if depth_remaining == 0 || start > n {
            jobs.push(PrefixJob {
                start,
                chosen: chosen.clone(),
                used_sums: used_sums.clone(),
            });
            return;
        }

        rec(
            n,
            start + 1,
            depth_remaining - 1,
            chosen,
            used_sums,
            scratch,
            jobs,
        );

        if !extension_sums(chosen, used_sums, start, scratch) {
            return;
        }
        let new_sums = scratch.clone();
        chosen.push(start);
        for value in &new_sums {
            used_sums[*value] = true;
        }
        rec(
            n,
            start + 1,
            depth_remaining - 1,
            chosen,
            used_sums,
            scratch,
            jobs,
        );
        chosen.pop();
        for value in &new_sums {
            used_sums[*value] = false;
        }
    }

    let mut jobs = Vec::new();
    let mut chosen = Vec::new();
    let mut used_sums = vec![false; 2 * n + 1];
    let mut scratch = Vec::new();
    rec(
        n,
        0,
        depth,
        &mut chosen,
        &mut used_sums,
        &mut scratch,
        &mut jobs,
    );
    jobs
}

fn compute_profile_deviation(chosen: &[usize], sqrt_n: f64) -> (f64, Option<usize>) {
    let mut max_value = f64::NEG_INFINITY;
    let mut max_index = 0;
    for (i, &a) in chosen.iter().enumerate() {
        let value = (a as f64 - (i + 1) as f64 * sqrt_n).abs();
        if value > max_value {
            max_value = value;
            max_index = i;
        }
    }
    (max_value, Some(max_index))
}

fn compute_prefix_deviation(chosen: &[usize], n: usize, slope: f64) -> (f64, Option<usize>) {
    let mut prefix_size = 0usize;
    let mut max_dev = -1.0;
    let mut max_t = 0usize;
    let mut next_index = 0usize;
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

fn compute_recentered_prefix_deviation(
    chosen: &[usize],
    n: usize,
    sqrt_n: f64,
    endpoint_bridge: f64,
) -> (f64, Option<usize>) {
    let mut prefix_size = 0usize;
    let mut max_dev = -1.0;
    let mut max_t = 0usize;
    let mut next_index = 0usize;
    for t in 0..=n {
        while next_index < chosen.len() && chosen[next_index] <= t {
            prefix_size += 1;
            next_index += 1;
        }
        let bridge_share = if n == 0 {
            0.0
        } else {
            (t as f64 / n as f64) * endpoint_bridge
        };
        let dev = (t as f64 - prefix_size as f64 * sqrt_n - bridge_share).abs();
        if dev > max_dev {
            max_dev = dev;
            max_t = t;
        }
    }
    (max_dev, Some(max_t))
}

fn compute_explicit_mass_center(k: usize, sqrt_n: f64) -> f64 {
    ((k * (k - 1) / 2 + k) as f64) * sqrt_n
}

fn compute_density_adjusted_mass_center(k: usize, n: usize) -> f64 {
    n as f64 * (k + 1) as f64 / 2.0
}

fn compute_joint_observable_score(
    prefix_residual: f64,
    density_adjusted_mass_dev: f64,
    n: usize,
) -> f64 {
    let n_real = n as f64;
    prefix_residual / n_real.powf(7.0 / 8.0) + density_adjusted_mass_dev / n_real.powf(11.0 / 8.0)
}

fn metric7(summary: MetricSummary, n_pow_7_8: f64) -> Metric7 {
    Metric7 {
        normalized_by_n_pow_7_8: Normalized7 {
            min: summary.min_value / n_pow_7_8,
            mean: summary.mean_value / n_pow_7_8,
            max: summary.max_value / n_pow_7_8,
        },
        raw: summary,
    }
}

fn metric11(summary: MetricSummary, n_pow_11_8: f64) -> Metric11 {
    Metric11 {
        normalized_by_n_pow_11_8: Normalized11 {
            min: summary.min_value / n_pow_11_8,
            mean: summary.mean_value / n_pow_11_8,
            max: summary.max_value / n_pow_11_8,
        },
        raw: summary,
    }
}

fn analyze_maximizer_sets(
    n: usize,
    progress: bool,
    seed_depth: usize,
    initial_lower_bound_override: Option<usize>,
    parallel_prefix_depth: usize,
    frontier_k: usize,
    inheritance_source: Option<&InheritanceSource>,
    export_ground_face_masks: bool,
) -> Row {
    let start_time = Instant::now();
    if progress {
        eprintln!("starting n={n}");
    }
    let sqrt_n = (n as f64).sqrt();
    let n_real = n as f64;
    let n_pow_7_8 = n_real.powf(7.0 / 8.0);
    let n_pow_11_8 = n_real.powf(11.0 / 8.0);

    let initial_lower_bound =
        initial_lower_bound_override.unwrap_or_else(|| initial_lower_bound_size(n, seed_depth));
    let mut scan = ScanAccumulator::new(seed_depth, initial_lower_bound, frontier_k);

    fn bt(
        n: usize,
        start: usize,
        chosen: &mut Vec<usize>,
        used_sums: &mut Vec<bool>,
        scratch: &mut Vec<usize>,
        sqrt_n: f64,
        n_pow_7_8: f64,
        n_pow_11_8: f64,
        inheritance_source: Option<&InheritanceSource>,
        export_ground_face_masks: bool,
        scan: &mut ScanAccumulator,
    ) {
        scan.search_diagnostics.nodes_visited += 1;
        let remaining = if start <= n { n - start + 1 } else { 0 };
        if chosen.len() + remaining < scan.best_size {
            scan.search_diagnostics.size_bound_prunes += 1;
            return;
        }
        if start > n {
            scan.observe_terminal(
                n,
                chosen,
                sqrt_n,
                n_pow_7_8,
                n_pow_11_8,
                inheritance_source,
                export_ground_face_masks,
            );
            return;
        }

        bt(
            n,
            start + 1,
            chosen,
            used_sums,
            scratch,
            sqrt_n,
            n_pow_7_8,
            n_pow_11_8,
            inheritance_source,
            export_ground_face_masks,
            scan,
        );

        if !extension_sums(chosen, used_sums, start, scratch) {
            scan.search_diagnostics.sidon_extension_rejects += 1;
            return;
        }
        let new_sums = scratch.clone();
        chosen.push(start);
        for value in &new_sums {
            used_sums[*value] = true;
        }
        bt(
            n,
            start + 1,
            chosen,
            used_sums,
            scratch,
            sqrt_n,
            n_pow_7_8,
            n_pow_11_8,
            inheritance_source,
            export_ground_face_masks,
            scan,
        );
        chosen.pop();
        for value in &new_sums {
            used_sums[*value] = false;
        }
    }

    if parallel_prefix_depth == 0 {
        let mut chosen = Vec::new();
        let mut used_sums = vec![false; 2 * n + 1];
        let mut scratch = Vec::new();
        bt(
            n,
            0,
            &mut chosen,
            &mut used_sums,
            &mut scratch,
            sqrt_n,
            n_pow_7_8,
            n_pow_11_8,
            inheritance_source,
            export_ground_face_masks,
            &mut scan,
        );
    } else {
        let jobs = generate_prefix_jobs(n, parallel_prefix_depth);
        if progress {
            eprintln!(
                "n={n} parallel prefix depth={} jobs={}",
                parallel_prefix_depth,
                jobs.len()
            );
        }
        let workers: Vec<ScanAccumulator> = jobs
            .into_par_iter()
            .map(|job| {
                let mut worker = ScanAccumulator::new(seed_depth, initial_lower_bound, frontier_k);
                let mut chosen = job.chosen;
                let mut used_sums = job.used_sums;
                let mut scratch = Vec::new();
                bt(
                    n,
                    job.start,
                    &mut chosen,
                    &mut used_sums,
                    &mut scratch,
                    sqrt_n,
                    n_pow_7_8,
                    n_pow_11_8,
                    inheritance_source,
                    export_ground_face_masks,
                    &mut worker,
                );
                worker
            })
            .collect();
        for worker in workers {
            scan.merge_worker(worker);
        }
    }

    let gap = (scan.best_size as f64 - sqrt_n).abs();
    let prefix_drift = gap.max(1.0) * sqrt_n;
    let endpoint_bridge = n as f64 - scan.best_size as f64 * sqrt_n;
    let explicit_center = compute_explicit_mass_center(scan.best_size, sqrt_n);
    let density_adjusted_slope = n as f64 / scan.best_size as f64;
    let density_adjusted_mass_center = compute_density_adjusted_mass_center(scan.best_size, n);

    let profile_summary = scan.profile_acc.summary();
    let prefix_summary = scan.prefix_acc.summary();
    let prefix_residual_summary = scan.prefix_residual_acc.summary();
    let recentered_prefix_summary = scan.recentered_prefix_acc.summary();
    let density_adjusted_prefix_summary = scan.density_adjusted_prefix_acc.summary();
    let mass_summary = scan.mass_acc.summary();
    let density_adjusted_mass_summary = scan.density_adjusted_mass_acc.summary();
    let joint_summary = scan.joint_acc.summary();
    let same_best_witness =
        prefix_residual_summary.min_set == density_adjusted_mass_summary.min_set;

    let row = Row {
        n,
        maximizer_size_h_n: scan.best_size,
        maximizer_set_count: scan.best_count,
        maximizer_scan_runtime_sec: start_time.elapsed().as_secs_f64(),
        real_gap_from_sqrt: gap,
        real_deficiency_from_sqrt: (sqrt_n - scan.best_size as f64).max(0.0),
        general_prefix_drift: prefix_drift,
        endpoint_bridge,
        density_adjusted_slope,
        explicit_mass_center: explicit_center,
        density_adjusted_mass_center,
        normalizers: Normalizers {
            sqrt_n,
            n_pow_7_8,
            n_pow_11_8,
        },
        ordered_profile_deviation: metric7(profile_summary, n_pow_7_8),
        prefix_deviation: metric7(prefix_summary, n_pow_7_8),
        prefix_residual_after_general_drift: metric7(prefix_residual_summary.clone(), n_pow_7_8),
        recentered_prefix_deviation: metric7(recentered_prefix_summary, n_pow_7_8),
        density_adjusted_prefix_deviation: metric7(density_adjusted_prefix_summary, n_pow_7_8),
        mass_deviation: metric11(mass_summary, n_pow_11_8),
        density_adjusted_mass_deviation: metric11(
            density_adjusted_mass_summary.clone(),
            n_pow_11_8,
        ),
        joint_observable_score: JointMetric {
            raw: joint_summary.clone(),
        },
        observable_split: ObservableSplit {
            same_best_witness,
            prefix_best_witness: prefix_residual_summary.min_set,
            mass_best_witness: density_adjusted_mass_summary.min_set,
            prefix_best_mass_dev: scan.best_prefix_mass_dev,
            prefix_best_mass_ratio_n_pow_11_8: scan.best_prefix_mass_ratio,
            mass_best_prefix_residual: scan.best_mass_prefix_residual,
            mass_best_prefix_ratio_n_pow_7_8: scan.best_mass_prefix_ratio,
            best_joint_witness: joint_summary.min_set,
            best_joint_score: joint_summary.min_value,
            worst_joint_witness: joint_summary.max_set,
            worst_joint_score: joint_summary.max_value,
        },
        search_diagnostics: scan.search_diagnostics,
        top_k_frontier: TopKFrontier {
            k: frontier_k,
            prefix_best: scan.prefix_frontier.ranked(),
            density_adjusted_mass_best: scan.mass_frontier.ranked(),
            joint_best: scan.joint_frontier.ranked(),
        },
        inheritance_probe: build_inheritance_probe(
            inheritance_source,
            &scan.inheritance_probe,
            scan.best_count,
        ),
        ground_face_mask_export: build_mask_export(
            export_ground_face_masks,
            scan.best_count,
            scan.ground_face_masks,
        ),
    };
    if progress {
        eprintln!(
            "finished n={} h={} maximizers={} runtime={:.3}s nodes={} prunes={} rejects={}",
            row.n,
            row.maximizer_size_h_n,
            row.maximizer_set_count,
            row.maximizer_scan_runtime_sec,
            row.search_diagnostics.nodes_visited,
            row.search_diagnostics.size_bound_prunes,
            row.search_diagnostics.sidon_extension_rejects,
        );
    }
    row
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn prefix_diagnostic<F>(rows: &[Row], label: &str, residual_for: F) -> PrefixDriftDiagnostic
where
    F: Fn(&Row) -> (f64, f64, Option<usize>, Vec<usize>, f64),
{
    let mut best: Option<PrefixDriftDiagnostic> = None;
    for row in rows {
        let (residual, drift, location, witness_set, endpoint_bridge) = residual_for(row);
        let ratio = residual / row.normalizers.n_pow_7_8;
        if best
            .as_ref()
            .map(|current| ratio > current.worst_residual_ratio_n_pow_7_8)
            .unwrap_or(true)
        {
            best = Some(PrefixDriftDiagnostic {
                label: label.to_string(),
                worst_n: row.n,
                worst_residual: residual,
                worst_residual_ratio_n_pow_7_8: ratio,
                worst_raw_prefix_deviation: row.prefix_deviation.raw.max_value,
                worst_drift: drift,
                worst_set: witness_set,
                worst_t: location,
                worst_endpoint_bridge: endpoint_bridge,
            });
        }
    }
    best.unwrap()
}

fn summarize_prefix_drift(rows: &[Row]) -> PrefixDriftDiagnostics {
    PrefixDriftDiagnostics {
        raw_prefix: prefix_diagnostic(rows, "0", |row| {
            (
                row.prefix_deviation.raw.max_value,
                0.0,
                row.prefix_deviation.raw.max_location,
                row.prefix_deviation.raw.max_set.clone(),
                row.endpoint_bridge,
            )
        }),
        sqrt_only: prefix_diagnostic(rows, "sqrt(n)", |row| {
            let drift = row.normalizers.sqrt_n;
            (
                (row.prefix_deviation.raw.max_value - drift).max(0.0),
                drift,
                row.prefix_deviation.raw.max_location,
                row.prefix_deviation.raw.max_set.clone(),
                row.endpoint_bridge,
            )
        }),
        gap_sqrt_only: prefix_diagnostic(rows, "abs(card-sqrt(n)) * sqrt(n)", |row| {
            let drift = row.real_gap_from_sqrt * row.normalizers.sqrt_n;
            (
                (row.prefix_deviation.raw.max_value - drift).max(0.0),
                drift,
                row.prefix_deviation.raw.max_location,
                row.prefix_deviation.raw.max_set.clone(),
                row.endpoint_bridge,
            )
        }),
        current_theorem_drift: prefix_diagnostic(
            rows,
            "max(abs(card-sqrt(n)), 1) * sqrt(n)",
            |row| {
                (
                    row.prefix_residual_after_general_drift.raw.max_value,
                    row.general_prefix_drift,
                    row.prefix_residual_after_general_drift.raw.max_location,
                    row.prefix_residual_after_general_drift.raw.max_set.clone(),
                    row.endpoint_bridge,
                )
            },
        ),
        affine_endpoint_recentered: prefix_diagnostic(rows, "affine endpoint recentered", |row| {
            (
                row.recentered_prefix_deviation.raw.max_value,
                row.endpoint_bridge,
                row.recentered_prefix_deviation.raw.max_location,
                row.recentered_prefix_deviation.raw.max_set.clone(),
                row.endpoint_bridge,
            )
        }),
        density_adjusted_affine: prefix_diagnostic(rows, "density-adjusted affine", |row| {
            (
                row.density_adjusted_prefix_deviation.raw.max_value,
                0.0,
                row.density_adjusted_prefix_deviation.raw.max_location,
                row.density_adjusted_prefix_deviation.raw.max_set.clone(),
                row.endpoint_bridge,
            )
        }),
    }
}

fn summarize_mass_centers(rows: &[Row]) -> MassCenterDiagnostics {
    fn one<F>(rows: &[Row], label: &str, get: F) -> MassCenterDiagnostic
    where
        F: Fn(&Row) -> (f64, f64, Vec<usize>),
    {
        let mut best: Option<MassCenterDiagnostic> = None;
        for row in rows {
            let (residual, ratio, witness_set) = get(row);
            if best
                .as_ref()
                .map(|current| ratio > current.worst_residual_ratio_n_pow_11_8)
                .unwrap_or(true)
            {
                best = Some(MassCenterDiagnostic {
                    label: label.to_string(),
                    worst_n: row.n,
                    worst_residual: residual,
                    worst_residual_ratio_n_pow_11_8: ratio,
                    worst_set: witness_set,
                });
            }
        }
        best.unwrap()
    }

    MassCenterDiagnostics {
        sqrt_center: one(rows, "sqrt(n)-centered mass template", |row| {
            (
                row.mass_deviation.raw.max_value,
                row.mass_deviation.normalized_by_n_pow_11_8.max,
                row.mass_deviation.raw.max_set.clone(),
            )
        }),
        density_adjusted_center: one(rows, "density-adjusted mass template", |row| {
            (
                row.density_adjusted_mass_deviation.raw.max_value,
                row.density_adjusted_mass_deviation
                    .normalized_by_n_pow_11_8
                    .max,
                row.density_adjusted_mass_deviation.raw.max_set.clone(),
            )
        }),
    }
}

fn summarize_observable_split(rows: &[Row]) -> ObservableSplitDiagnostics {
    let zero_tolerance = 1e-9;
    let mut same_best_ns = Vec::new();
    let mut mass_best_near_zero_prefix_ns = Vec::new();
    let mut prefix_best_zero_mass_ns = Vec::new();
    let mut mass_penalty_le_prefix_penalty_ns = Vec::new();
    let mut mass_penalty_gt_prefix_penalty_ns = Vec::new();
    let mut joint_is_mass_best_ns = Vec::new();
    let mut joint_is_prefix_best_ns = Vec::new();
    let mut joint_is_third_ns = Vec::new();
    let mut mass_best_prefix_ratios = Vec::new();
    let mut prefix_best_mass_ratios = Vec::new();
    let mut strongest_split_index = 0usize;
    let mut strongest_split_score = f64::NEG_INFINITY;
    let mut best_joint_index = 0usize;
    let mut best_joint_score = f64::INFINITY;

    for (i, row) in rows.iter().enumerate() {
        let split = &row.observable_split;
        if split.same_best_witness {
            same_best_ns.push(row.n);
        }
        mass_best_prefix_ratios.push(split.mass_best_prefix_ratio_n_pow_7_8);
        prefix_best_mass_ratios.push(split.prefix_best_mass_ratio_n_pow_11_8);
        if split.mass_best_prefix_residual <= zero_tolerance {
            mass_best_near_zero_prefix_ns.push(row.n);
        }
        if split.prefix_best_mass_dev == 0.0 {
            prefix_best_zero_mass_ns.push(row.n);
        }
        if split.mass_best_prefix_ratio_n_pow_7_8 <= split.prefix_best_mass_ratio_n_pow_11_8 {
            mass_penalty_le_prefix_penalty_ns.push(row.n);
        } else {
            mass_penalty_gt_prefix_penalty_ns.push(row.n);
        }
        if split.best_joint_witness == split.mass_best_witness {
            joint_is_mass_best_ns.push(row.n);
        }
        if split.best_joint_witness == split.prefix_best_witness {
            joint_is_prefix_best_ns.push(row.n);
        }
        if split.best_joint_witness != split.mass_best_witness
            && split.best_joint_witness != split.prefix_best_witness
        {
            joint_is_third_ns.push(row.n);
        }
        let split_score =
            split.prefix_best_mass_ratio_n_pow_11_8 + split.mass_best_prefix_ratio_n_pow_7_8;
        if split_score > strongest_split_score {
            strongest_split_score = split_score;
            strongest_split_index = i;
        }
        if split.best_joint_score < best_joint_score {
            best_joint_score = split.best_joint_score;
            best_joint_index = i;
        }
    }

    let strongest_row = &rows[strongest_split_index];
    let strongest_split = &strongest_row.observable_split;
    let best_joint_row = &rows[best_joint_index];
    let best_joint = &best_joint_row.observable_split;

    ObservableSplitDiagnostics {
        same_best_count: same_best_ns.len(),
        same_best_ns,
        strongest_split_n: strongest_row.n,
        strongest_split_score,
        strongest_split_prefix_best_witness: strongest_split.prefix_best_witness.clone(),
        strongest_split_mass_best_witness: strongest_split.mass_best_witness.clone(),
        strongest_split_prefix_best_mass_ratio_n_pow_11_8: strongest_split
            .prefix_best_mass_ratio_n_pow_11_8,
        strongest_split_mass_best_prefix_ratio_n_pow_7_8: strongest_split
            .mass_best_prefix_ratio_n_pow_7_8,
        best_joint_n: best_joint_row.n,
        best_joint_witness: best_joint.best_joint_witness.clone(),
        best_joint_score: best_joint.best_joint_score,
        finite_compatibility: FiniteCompatibility {
            zero_tolerance,
            mass_best_near_zero_prefix_count: mass_best_near_zero_prefix_ns.len(),
            mass_best_near_zero_prefix_ns,
            prefix_best_zero_mass_count: prefix_best_zero_mass_ns.len(),
            prefix_best_zero_mass_ns,
            mass_best_prefix_ratio_max_n_pow_7_8: mass_best_prefix_ratios
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max),
            mass_best_prefix_ratio_mean_n_pow_7_8: mean(&mass_best_prefix_ratios),
            prefix_best_mass_ratio_max_n_pow_11_8: prefix_best_mass_ratios
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max),
            prefix_best_mass_ratio_mean_n_pow_11_8: mean(&prefix_best_mass_ratios),
            mass_penalty_le_prefix_penalty_count: mass_penalty_le_prefix_penalty_ns.len(),
            mass_penalty_le_prefix_penalty_ns,
            mass_penalty_gt_prefix_penalty_count: mass_penalty_gt_prefix_penalty_ns.len(),
            mass_penalty_gt_prefix_penalty_ns,
            joint_is_mass_best_count: joint_is_mass_best_ns.len(),
            joint_is_mass_best_ns,
            joint_is_prefix_best_count: joint_is_prefix_best_ns.len(),
            joint_is_prefix_best_ns,
            joint_is_third_count: joint_is_third_ns.len(),
            joint_is_third_ns,
        },
    }
}

fn summarize_first_hit_vs_face_aware(rows: &[Row]) -> FirstHitVsFaceAwareDiagnostics {
    let zero_tolerance = 1e-9;
    let mut first_hit_compatible_ns = Vec::new();
    let mut first_hit_failure_ns = Vec::new();
    let mut face_aware_near_zero_joint_ns = Vec::new();
    let mut face_aware_nonzero_joint_ns = Vec::new();
    let mut first_hit_failures_recovered_by_face_ns = Vec::new();
    let mut first_hit_failures_persisting_as_nonzero_face_ns = Vec::new();
    let mut max_face_aware_joint: Option<(usize, FrontierWitness)> = None;

    for row in rows {
        let split = &row.observable_split;
        let first_hit_compatible =
            split.mass_best_prefix_ratio_n_pow_7_8 <= split.prefix_best_mass_ratio_n_pow_11_8;
        if first_hit_compatible {
            first_hit_compatible_ns.push(row.n);
        } else {
            first_hit_failure_ns.push(row.n);
        }

        let joint = row.top_k_frontier.joint_best.first();
        if let Some(joint) = joint {
            let face_aware_near_zero = joint.joint_score <= zero_tolerance;
            if face_aware_near_zero {
                face_aware_near_zero_joint_ns.push(row.n);
            } else {
                face_aware_nonzero_joint_ns.push(row.n);
            }

            if !first_hit_compatible && face_aware_near_zero {
                first_hit_failures_recovered_by_face_ns.push(row.n);
            }
            if !first_hit_compatible && !face_aware_near_zero {
                first_hit_failures_persisting_as_nonzero_face_ns.push(row.n);
            }

            match &max_face_aware_joint {
                Some((_, current)) if current.joint_score >= joint.joint_score => {}
                _ => max_face_aware_joint = Some((row.n, joint.clone())),
            }
        } else if !first_hit_compatible {
            first_hit_failures_persisting_as_nonzero_face_ns.push(row.n);
        }
    }

    let (
        max_face_aware_joint_n,
        max_face_aware_joint_score,
        max_face_aware_joint_prefix_ratio_n_pow_7_8,
        max_face_aware_joint_mass_ratio_n_pow_11_8,
        max_face_aware_joint_witness,
    ) = match max_face_aware_joint {
        Some((n, witness)) => (
            Some(n),
            Some(witness.joint_score),
            Some(witness.prefix_ratio_n_pow_7_8),
            Some(witness.density_adjusted_mass_ratio_n_pow_11_8),
            Some(witness.witness),
        ),
        None => (None, None, None, None, None),
    };

    FirstHitVsFaceAwareDiagnostics {
        zero_tolerance,
        first_hit_compatible_count: first_hit_compatible_ns.len(),
        first_hit_compatible_ns,
        first_hit_failure_count: first_hit_failure_ns.len(),
        first_hit_failure_ns,
        face_aware_near_zero_joint_count: face_aware_near_zero_joint_ns.len(),
        face_aware_near_zero_joint_ns,
        face_aware_nonzero_joint_count: face_aware_nonzero_joint_ns.len(),
        face_aware_nonzero_joint_ns,
        first_hit_failures_recovered_by_face_count: first_hit_failures_recovered_by_face_ns.len(),
        first_hit_failures_recovered_by_face_ns,
        first_hit_failures_persisting_as_nonzero_face_count:
            first_hit_failures_persisting_as_nonzero_face_ns.len(),
        first_hit_failures_persisting_as_nonzero_face_ns,
        max_face_aware_joint_n,
        max_face_aware_joint_score,
        max_face_aware_joint_prefix_ratio_n_pow_7_8,
        max_face_aware_joint_mass_ratio_n_pow_11_8,
        max_face_aware_joint_witness,
    }
}

fn build_report(results: &Results) -> String {
    let rows = &results.results_by_n;
    let total_max_sets: u64 = rows.iter().map(|row| row.maximizer_set_count).sum();
    let max_prefix_ratio = rows
        .iter()
        .map(|row| {
            row.prefix_residual_after_general_drift
                .normalized_by_n_pow_7_8
                .max
        })
        .fold(f64::NEG_INFINITY, f64::max);
    let max_mass_ratio = rows
        .iter()
        .map(|row| row.mass_deviation.normalized_by_n_pow_11_8.max)
        .fold(f64::NEG_INFINITY, f64::max);
    let split_diag = &results.observable_split_diagnostics;
    let face_diag = &results.first_hit_vs_face_aware_diagnostics;
    let total_n = rows.len();

    let mut lines = Vec::new();
    lines.push(format!(
        "# {} — Exact Maximizer Sidon Rigidity Scan",
        results.experiment_id
    ));
    lines.push(String::new());
    lines.push("## Identification".to_string());
    lines.push(String::new());
    lines.push("| Field | Value |".to_string());
    lines.push("|---|---|".to_string());
    lines.push(format!("| Experiment ID | {} |", results.experiment_id));
    lines.push("| Erdős Problem | #30 — finite Sidon set rigidity |".to_string());
    lines
        .push("| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |".to_string());
    lines.push(format!("| Implementation | {} |", results.implementation));
    lines.push(format!(
        "| Scan window | n = {} through n = {} |",
        results.n_min, results.n_max
    ));
    lines.push(format!(
        "| Inheritance source | {} |",
        results
            .inheritance_source_results
            .as_deref()
            .unwrap_or("none")
    ));
    lines.push(format!(
        "| Ground-face mask export | {} |",
        results.export_ground_face_masks
    ));
    lines.push("| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |".to_string());
    lines.push("| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |".to_string());
    lines.push(String::new());
    lines.push("## Answer".to_string());
    lines.push(String::new());
    lines.push(format!(
        "The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across {total_max_sets} exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale."
    ));
    lines.push(String::new());
    lines.push(format!(
        "The worst observed prefix residual after subtracting the general drift stayed below {max_prefix_ratio:.4} in n^(7/8) units, and the worst observed mass deviation stayed below {max_mass_ratio:.4} in n^(11/8) units."
    ));
    lines.push(String::new());
    lines.push("## Finite Compatibility Candidate".to_string());
    lines.push(String::new());
    lines.push(format!(
        "The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in {} of {} values of n. The exceptional n-values are {:?}.",
        split_diag.finite_compatibility.mass_penalty_le_prefix_penalty_count,
        total_n,
        split_diag.finite_compatibility.mass_penalty_gt_prefix_penalty_ns
    ));
    lines.push(String::new());
    lines.push(format!(
        "The best joint witness equals the mass-best witness in {} of {} values, equals the prefix-best witness in {} of {} values, and is a third witness in {} of {} values.",
        split_diag.finite_compatibility.joint_is_mass_best_count,
        total_n,
        split_diag.finite_compatibility.joint_is_prefix_best_count,
        total_n,
        split_diag.finite_compatibility.joint_is_third_count,
        total_n
    ));
    lines.push(String::new());
    lines.push("## First-hit vs Face-aware Summary".to_string());
    lines.push(String::new());
    lines.push(format!(
        "First-hit compatibility holds in {} of {} values of n. The first-hit failures are {:?}.",
        face_diag.first_hit_compatible_count, total_n, face_diag.first_hit_failure_ns
    ));
    lines.push(String::new());
    lines.push(format!(
        "The top-k joint frontier has near-zero joint score in {} of {} values. Nonzero joint-frontier n-values are {:?}.",
        face_diag.face_aware_near_zero_joint_count,
        total_n,
        face_diag.face_aware_nonzero_joint_ns
    ));
    lines.push(String::new());
    lines.push(format!(
        "First-hit failures recovered by the face-aware frontier: {:?}. First-hit failures persisting as nonzero face-level handoff cases: {:?}.",
        face_diag.first_hit_failures_recovered_by_face_ns,
        face_diag.first_hit_failures_persisting_as_nonzero_face_ns
    ));
    if let (Some(n), Some(score), Some(prefix_ratio), Some(mass_ratio), Some(witness)) = (
        face_diag.max_face_aware_joint_n,
        face_diag.max_face_aware_joint_score,
        face_diag.max_face_aware_joint_prefix_ratio_n_pow_7_8,
        face_diag.max_face_aware_joint_mass_ratio_n_pow_11_8,
        &face_diag.max_face_aware_joint_witness,
    ) {
        lines.push(String::new());
        lines.push(format!(
            "Largest top-k joint score: n = {}, score = {:.6}, prefix ratio = {:.6}, mass ratio = {:.6}, witness = {:?}.",
            n, score, prefix_ratio, mass_ratio, witness
        ));
    }
    lines.push(String::new());
    lines.push(
        "| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |"
            .to_string(),
    );
    lines.push("|---|---|---:|---:|---:|---|".to_string());
    for row in rows {
        let first_hit_compatible = row.observable_split.mass_best_prefix_ratio_n_pow_7_8
            <= row.observable_split.prefix_best_mass_ratio_n_pow_11_8;
        if let Some(joint) = row.top_k_frontier.joint_best.first() {
            lines.push(format!(
                "| {} | {} | {:.6} | {:.6} | {:.6} | {:?} |",
                row.n,
                if first_hit_compatible { "yes" } else { "no" },
                joint.joint_score,
                joint.prefix_ratio_n_pow_7_8,
                joint.density_adjusted_mass_ratio_n_pow_11_8,
                joint.witness,
            ));
        } else {
            lines.push(format!(
                "| {} | {} | NA | NA | NA | [] |",
                row.n,
                if first_hit_compatible { "yes" } else { "no" },
            ));
        }
    }
    lines.push(String::new());
    lines.push("## Per-n Summary".to_string());
    lines.push(String::new());
    lines.push("| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |".to_string());
    lines.push("|---|---|---|---|---|---|---|---|---|---|---|".to_string());
    for row in rows {
        lines.push(format!(
            "| {} | {} | {} | {} | {} | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} |",
            row.n,
            row.maximizer_size_h_n,
            row.maximizer_set_count,
            row.search_diagnostics.seed_depth,
            row.search_diagnostics.initial_lower_bound_size,
            row.search_diagnostics.nodes_visited,
            row.search_diagnostics.size_bound_prunes,
            row.prefix_residual_after_general_drift
                .normalized_by_n_pow_7_8
                .max,
            row.mass_deviation.normalized_by_n_pow_11_8.max,
            row.observable_split.mass_best_prefix_ratio_n_pow_7_8,
            row.observable_split.prefix_best_mass_ratio_n_pow_11_8
        ));
    }
    if rows.iter().any(|row| row.inheritance_probe.is_some()) {
        lines.push(String::new());
        lines.push("## Plateau Inheritance Probe".to_string());
        lines.push(String::new());
        lines.push("The inheritance probe checks whether the current exact maximizer face contains the previous exact face and the previous face shifted by `+1` during the same exact traversal used for h(n).".to_string());
        lines.push(String::new());
        lines.push("| n | source n | previous face | inherited union | prev exact present | prev +1 present | inherited present | new face | union contained |".to_string());
        lines.push("|---:|---:|---:|---:|---:|---:|---:|---:|---|".to_string());
        for row in rows {
            if let Some(probe) = &row.inheritance_probe {
                lines.push(format!(
                    "| {} | {} | {} | {} | {} | {} | {} | {} | {} |",
                    row.n,
                    probe.source_n,
                    probe.previous_face_count,
                    probe.inherited_union_count,
                    probe.exact_previous_persistence_count,
                    probe.plus_one_previous_persistence_count,
                    probe.inherited_union_present_count,
                    probe.new_face_count,
                    probe.contains_inherited_union
                ));
            }
        }
    }
    lines.push(String::new());
    lines.push("## Top-k Frontier Witnesses".to_string());
    lines.push(String::new());
    lines.push(
        "| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |".to_string(),
    );
    lines.push("|---|---|---:|---:|---:|---:|---|".to_string());
    for row in rows {
        for (frontier_name, witnesses) in [
            ("prefix", &row.top_k_frontier.prefix_best),
            ("mass", &row.top_k_frontier.density_adjusted_mass_best),
            ("joint", &row.top_k_frontier.joint_best),
        ] {
            for witness in witnesses.iter().take(3) {
                lines.push(format!(
                    "| {} | {} | {} | {:.4} | {:.4} | {:.4} | {:?} |",
                    row.n,
                    frontier_name,
                    witness.rank,
                    witness.prefix_ratio_n_pow_7_8,
                    witness.density_adjusted_mass_ratio_n_pow_11_8,
                    witness.joint_score,
                    witness.witness,
                ));
            }
        }
    }
    lines.push(String::new());
    lines.push("## Artifacts".to_string());
    lines.push(String::new());
    lines.push("| File | Type |".to_string());
    lines.push("|---|---|".to_string());
    lines.push(format!(
        "| {}_RESULTS.json | Structured exact-enumeration results |",
        results.experiment_id
    ));
    lines.push(format!(
        "| {}_REPORT.md | Human-readable report |",
        results.experiment_id
    ));
    lines.push(format!(
        "| {}_RESULTS.sha256 | Integrity checksum |",
        results.experiment_id
    ));
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
    let experiment_id = cli.experiment_id.unwrap_or_else(|| {
        format!(
            "EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-{}-{}-2026-05-01",
            cli.n_min, cli.n_max
        )
    });

    let inheritance_source = cli
        .inheritance_source_results
        .as_deref()
        .map(load_inheritance_source)
        .transpose()?;
    let t0 = Instant::now();
    let rows: Vec<Row> = (cli.n_min..=cli.n_max)
        .map(|n| {
            analyze_maximizer_sets(
                n,
                cli.progress,
                cli.seed_depth,
                cli.initial_lower_bound,
                cli.parallel_prefix_depth,
                cli.frontier_k,
                inheritance_source.as_ref(),
                cli.export_ground_face_masks,
            )
        })
        .collect();
    let results = Results {
        experiment_id,
        date: "2026-05-01".to_string(),
        description: "Exact h(n) maximizer scan against the general prefix-drift and mass theorems"
            .to_string(),
        erdos_problem: 30,
        scan_mode: "maximizer".to_string(),
        implementation: "Rust exact maximizer scanner".to_string(),
        n_min: cli.n_min,
        n_max: cli.n_max,
        inheritance_source_results: inheritance_source
            .as_ref()
            .map(|source| source.results_path.clone()),
        export_ground_face_masks: cli.export_ground_face_masks,
        runtime_sec: t0.elapsed().as_secs_f64(),
        prefix_drift_diagnostics: summarize_prefix_drift(&rows),
        mass_center_diagnostics: summarize_mass_centers(&rows),
        observable_split_diagnostics: summarize_observable_split(&rows),
        first_hit_vs_face_aware_diagnostics: summarize_first_hit_vs_face_aware(&rows),
        results_by_n: rows,
    };

    write_outputs(&results, &cli.output_dir)?;
    println!(
        "wrote packet {} to {}",
        results.experiment_id,
        cli.output_dir.display()
    );
    Ok(())
}
