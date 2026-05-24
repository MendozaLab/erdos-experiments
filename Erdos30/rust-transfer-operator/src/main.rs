use clap::Parser;
use rayon::prelude::*;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

const OCCUPIED_SUFFIX_WIDTH: usize = 8;

#[derive(Parser, Debug)]
#[command(name = "erdos30-transfer-operator")]
#[command(about = "Exact PMF transfer-operator parity scan for Erdős #30")]
struct Cli {
    #[arg(long, default_value_t = 10)]
    n_min: usize,

    #[arg(long, default_value_t = 30)]
    n_max: usize,

    #[arg(long, default_value_t = 5)]
    frontier_k: usize,

    #[arg(long)]
    prune_deficiency: Option<usize>,

    #[arg(long)]
    experiment_id: Option<String>,

    #[arg(long)]
    export_ground_face_cap: Option<usize>,

    #[arg(long)]
    export_ground_face_distance_edge_cap: Option<usize>,

    #[arg(long, default_value_t = 0)]
    parallel_prefix_depth: usize,

    #[arg(long)]
    inheritance_source_results: Option<PathBuf>,

    #[arg(long, default_value_t = false)]
    export_ground_face_masks: bool,

    #[arg(long, default_value = "erdos-experiments/results/erdos-30")]
    output_dir: PathBuf,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct State {
    occupied_mask: u128,
    used_differences_mask: u128,
    cardinality: u8,
}

impl State {
    fn empty() -> Self {
        Self {
            occupied_mask: 0,
            used_differences_mask: 0,
            cardinality: 0,
        }
    }

    fn occupied_sites(&self) -> Vec<usize> {
        bits_to_sites(self.occupied_mask)
    }

    fn occupied_suffix(&self, width: usize) -> Vec<usize> {
        let sites = self.occupied_sites();
        if sites.len() <= width {
            sites
        } else {
            sites[sites.len() - width..].to_vec()
        }
    }

    fn try_occupy(&self, x: usize) -> Option<Self> {
        let x_bit = bit(x)?;
        if self.occupied_mask & x_bit != 0 {
            return None;
        }

        let mut new_diff_bits = 0u128;
        let mut occupied = self.occupied_mask;
        while occupied != 0 {
            let a = occupied.trailing_zeros() as usize;
            let diff = x.abs_diff(a);
            if diff == 0 {
                return None;
            }
            let diff_bit = bit(diff)?;
            if (self.used_differences_mask | new_diff_bits) & diff_bit != 0 {
                return None;
            }
            new_diff_bits |= diff_bit;
            occupied &= occupied - 1;
        }

        Some(Self {
            occupied_mask: self.occupied_mask | x_bit,
            used_differences_mask: self.used_differences_mask | new_diff_bits,
            cardinality: self.cardinality + 1,
        })
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

    fn observe(&mut self, candidate: FrontierCandidate) {
        if self.k == 0 {
            return;
        }
        self.entries.push(candidate);
        self.entries.sort_by(candidate_cmp);
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
    zero_temperature_field_interpretation: String,
    prefix_best: Vec<FrontierWitness>,
    density_adjusted_mass_best: Vec<FrontierWitness>,
    joint_best: Vec<FrontierWitness>,
}

#[derive(Serialize)]
struct GroundFaceExportWitness {
    index: usize,
    witness: Vec<usize>,
    prefix_residual: f64,
    prefix_ratio_n_pow_7_8: f64,
    prefix_location: Option<usize>,
    density_adjusted_mass_dev: f64,
    density_adjusted_mass_ratio_n_pow_11_8: f64,
    joint_score: f64,
    is_prefix_winner: bool,
    is_mass_winner: bool,
    is_joint_winner: bool,
    is_pareto_minimal: bool,
}

#[derive(Serialize)]
struct GroundFaceDistanceEdge {
    left_index: usize,
    right_index: usize,
    symmetric_difference_distance: usize,
}

#[derive(Serialize)]
struct GroundFaceExposedWinners {
    prefix_winner_indices: Vec<usize>,
    mass_winner_indices: Vec<usize>,
    joint_winner_indices: Vec<usize>,
    pareto_minimal_indices: Vec<usize>,
    field_selection_split: bool,
}

#[derive(Serialize)]
struct GroundFaceExport {
    status: String,
    cap: usize,
    exact_maximizer_count: u64,
    exported_count: usize,
    coordinate_space: String,
    exposed_winners: Option<GroundFaceExposedWinners>,
    witnesses: Vec<GroundFaceExportWitness>,
    distance_edge_status: String,
    distance_edge_cap: Option<usize>,
    distance_edge_total_count: usize,
    distance_edges: Vec<GroundFaceDistanceEdge>,
}

#[derive(Serialize)]
struct GroundFaceMaskExport {
    status: String,
    exact_maximizer_count: u64,
    exported_count: usize,
    encoding: String,
    masks: Vec<String>,
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
struct StateModel {
    representation: String,
    occupied_mask_bits: String,
    used_differences_mask_bits: String,
    occupied_suffix_width: usize,
    transition_rule: String,
}

#[derive(Serialize)]
struct LayerCount {
    cardinality: usize,
    count: u64,
}

#[derive(Serialize)]
struct DeficiencyLayer {
    deficiency: usize,
    cardinality: usize,
    count: u64,
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
    layer_counts: Vec<LayerCount>,
    deficiency_layers: Vec<DeficiencyLayer>,
}

#[derive(Serialize)]
struct StateEncodingWitness {
    occupied_sites: Vec<usize>,
    occupied_suffix: Vec<usize>,
    used_differences: Vec<usize>,
    cardinality: usize,
}

#[derive(Serialize)]
struct TransferDiagnostics {
    traversal_mode: String,
    sites_processed: usize,
    terminal_state_count: u64,
    peak_frontier_state_count: u64,
    peak_frontier_after_site: usize,
    transitions_attempted: u64,
    legal_occupy_transitions: u64,
    rejected_occupy_transitions: u64,
    prune_deficiency: Option<usize>,
    pruning_target_h_n: Option<usize>,
    pruning_min_cardinality: Option<usize>,
    reachability_pruned_states: u64,
}

#[derive(Serialize)]
struct ParityReference {
    exact_h_n: usize,
    exact_maximizer_count: u64,
    h_n_matches: bool,
    maximizer_count_matches: bool,
}

#[derive(Serialize)]
struct Row {
    n: usize,
    h_n: usize,
    maximizer_count: u64,
    first_maximizer_witness: Vec<usize>,
    representative_ground_state_encoding: StateEncodingWitness,
    spectral_observables: SpectralObservables,
    transfer_diagnostics: TransferDiagnostics,
    parity_reference: Option<ParityReference>,
    top_k_frontier: TopKFrontier,
    ground_face_export: Option<GroundFaceExport>,
    ground_face_mask_export: Option<GroundFaceMaskExport>,
    inheritance_probe: Option<InheritanceProbe>,
}

#[derive(Serialize)]
struct Results {
    experiment_id: String,
    date: String,
    erdos_problem: usize,
    scan_mode: String,
    implementation: String,
    state_model: StateModel,
    n_min: usize,
    n_max: usize,
    frontier_k: usize,
    prune_deficiency: Option<usize>,
    export_ground_face_cap: Option<usize>,
    export_ground_face_distance_edge_cap: Option<usize>,
    export_ground_face_masks: bool,
    inheritance_source_results: Option<String>,
    total_runtime_sec: f64,
    rows: Vec<Row>,
    parity_summary: ParitySummary,
}

#[derive(Serialize)]
struct ParitySummary {
    checked_count: usize,
    h_n_match_count: usize,
    maximizer_count_match_count: usize,
    mismatch_ns: Vec<usize>,
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

fn sites_to_mask(sites: &[usize]) -> Option<u128> {
    let mut mask = 0u128;
    for &site in sites {
        mask |= bit(site)?;
    }
    Some(mask)
}

fn shifted_mask(sites: &[usize], shift: usize) -> Option<u128> {
    let shifted: Vec<usize> = sites.iter().map(|site| site + shift).collect();
    sites_to_mask(&shifted)
}

fn load_inheritance_source(path: &Path) -> Result<InheritanceSource, Box<dyn std::error::Error>> {
    let payload: serde_json::Value = serde_json::from_slice(&fs::read(path)?)?;
    let rows = payload
        .get("rows")
        .and_then(|value| value.as_array())
        .ok_or("inheritance source has no rows array")?;
    let row = rows
        .last()
        .ok_or("inheritance source rows array is empty")?;
    let source_n = row
        .get("n")
        .and_then(|value| value.as_u64())
        .ok_or("inheritance source row has no numeric n")? as usize;
    let export = row.get("ground_face_export");
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
            let shifted_sites: Vec<usize> = bits_to_sites(mask).into_iter().map(|site| site + 1).collect();
            shifted_masks.insert(
                sites_to_mask(&shifted_sites).ok_or("shifted inheritance mask exceeds bitset width")?,
            );
        }
    } else {
        let export = export.ok_or("inheritance source row has neither mask export nor ground_face_export")?;
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
            previous_masks.insert(sites_to_mask(&sites).ok_or("inheritance witness site exceeds bitset width")?);
            shifted_masks
                .insert(shifted_mask(&sites, 1).ok_or("shifted inheritance witness exceeds bitset width")?);
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
        contains_previous_face: accum.exact_previous_persistence_count as usize
            == source.previous_face_count,
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

fn candidate_cmp(a: &FrontierCandidate, b: &FrontierCandidate) -> Ordering {
    a.objective_value
        .total_cmp(&b.objective_value)
        .then_with(|| a.witness.cmp(&b.witness))
}

fn compute_prefix_residual(chosen: &[usize], n: usize) -> (f64, Option<usize>) {
    if chosen.is_empty() {
        return (n as f64, Some(n));
    }
    let sqrt_n = (n as f64).sqrt();
    let general_drift = ((chosen.len() as f64 - sqrt_n).abs()).max(1.0) * sqrt_n;
    let mut prefix_size = 0usize;
    let mut next_index = 0usize;
    let mut max_prefix_deviation = f64::NEG_INFINITY;
    let mut max_t = 0usize;

    for t in 0..=n {
        while next_index < chosen.len() && chosen[next_index] <= t {
            prefix_size += 1;
            next_index += 1;
        }
        let deviation = (t as f64 - prefix_size as f64 * sqrt_n).abs();
        if deviation > max_prefix_deviation {
            max_prefix_deviation = deviation;
            max_t = t;
        }
    }

    ((max_prefix_deviation - general_drift).max(0.0), Some(max_t))
}

fn density_adjusted_mass_dev(chosen: &[usize], n: usize) -> f64 {
    let mass: usize = chosen.iter().sum();
    let center = n as f64 * (chosen.len() as f64 + 1.0) / 2.0;
    (mass as f64 - center).abs()
}

fn frontier_candidate(chosen: Vec<usize>, n: usize) -> FrontierCandidate {
    let n_real = n as f64;
    let (prefix_residual, prefix_location) = compute_prefix_residual(&chosen, n);
    let mass_dev = density_adjusted_mass_dev(&chosen, n);
    let prefix_ratio = if n == 0 {
        prefix_residual
    } else {
        prefix_residual / n_real.powf(7.0 / 8.0)
    };
    let mass_ratio = if n == 0 {
        mass_dev
    } else {
        mass_dev / n_real.powf(11.0 / 8.0)
    };
    FrontierCandidate {
        witness: chosen,
        objective_value: prefix_ratio + mass_ratio,
        prefix_residual,
        prefix_ratio_n_pow_7_8: prefix_ratio,
        prefix_location,
        density_adjusted_mass_dev: mass_dev,
        density_adjusted_mass_ratio_n_pow_11_8: mass_ratio,
        joint_score: prefix_ratio + mass_ratio,
    }
}

fn symmetric_difference_distance(left: &[usize], right: &[usize]) -> usize {
    let mut i = 0usize;
    let mut j = 0usize;
    let mut distance = 0usize;
    while i < left.len() || j < right.len() {
        if i == left.len() {
            distance += right.len() - j;
            break;
        }
        if j == right.len() {
            distance += left.len() - i;
            break;
        }
        match left[i].cmp(&right[j]) {
            Ordering::Less => {
                distance += 1;
                i += 1;
            }
            Ordering::Greater => {
                distance += 1;
                j += 1;
            }
            Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    distance
}

fn nearly_equal(left: f64, right: f64) -> bool {
    (left - right).abs() <= 1e-12
}

fn is_pareto_minimal(candidate: &FrontierCandidate, all: &[FrontierCandidate]) -> bool {
    !all.iter().any(|other| {
        let no_worse_prefix = other.prefix_ratio_n_pow_7_8 <= candidate.prefix_ratio_n_pow_7_8 + 1e-12;
        let no_worse_mass = other.density_adjusted_mass_ratio_n_pow_11_8
            <= candidate.density_adjusted_mass_ratio_n_pow_11_8 + 1e-12;
        let strictly_better_prefix =
            other.prefix_ratio_n_pow_7_8 < candidate.prefix_ratio_n_pow_7_8 - 1e-12;
        let strictly_better_mass = other.density_adjusted_mass_ratio_n_pow_11_8
            < candidate.density_adjusted_mass_ratio_n_pow_11_8 - 1e-12;
        no_worse_prefix && no_worse_mass && (strictly_better_prefix || strictly_better_mass)
    })
}

fn build_ground_face_export(
    cap: Option<usize>,
    distance_edge_cap: Option<usize>,
    exact_maximizer_count: u64,
    mut candidates: Vec<FrontierCandidate>,
) -> Option<GroundFaceExport> {
    let cap = cap?;
    if exact_maximizer_count as usize > cap {
        return Some(GroundFaceExport {
            status: "SKIP_EXPORTED_FACE_TOO_LARGE".to_string(),
            cap,
            exact_maximizer_count,
            exported_count: 0,
            coordinate_space:
                "(prefix_ratio_n_pow_7_8, density_adjusted_mass_ratio_n_pow_11_8, joint_score)"
                    .to_string(),
            exposed_winners: None,
            witnesses: Vec::new(),
            distance_edge_status: "NOT_COMPUTED_FACE_TOO_LARGE".to_string(),
            distance_edge_cap,
            distance_edge_total_count: 0,
            distance_edges: Vec::new(),
        });
    }

    candidates.sort_by(|left, right| left.witness.cmp(&right.witness));
    let min_prefix = candidates
        .iter()
        .map(|candidate| candidate.prefix_residual)
        .fold(f64::INFINITY, f64::min);
    let min_mass = candidates
        .iter()
        .map(|candidate| candidate.density_adjusted_mass_dev)
        .fold(f64::INFINITY, f64::min);
    let min_joint = candidates
        .iter()
        .map(|candidate| candidate.joint_score)
        .fold(f64::INFINITY, f64::min);

    let mut prefix_winner_indices = Vec::new();
    let mut mass_winner_indices = Vec::new();
    let mut joint_winner_indices = Vec::new();
    let mut pareto_minimal_indices = Vec::new();
    let mut witnesses = Vec::new();

    for (index, candidate) in candidates.iter().enumerate() {
        let is_prefix_winner = nearly_equal(candidate.prefix_residual, min_prefix);
        let is_mass_winner = nearly_equal(candidate.density_adjusted_mass_dev, min_mass);
        let is_joint_winner = nearly_equal(candidate.joint_score, min_joint);
        let is_pareto_minimal = is_pareto_minimal(candidate, &candidates);
        if is_prefix_winner {
            prefix_winner_indices.push(index);
        }
        if is_mass_winner {
            mass_winner_indices.push(index);
        }
        if is_joint_winner {
            joint_winner_indices.push(index);
        }
        if is_pareto_minimal {
            pareto_minimal_indices.push(index);
        }
        witnesses.push(GroundFaceExportWitness {
            index,
            witness: candidate.witness.clone(),
            prefix_residual: candidate.prefix_residual,
            prefix_ratio_n_pow_7_8: candidate.prefix_ratio_n_pow_7_8,
            prefix_location: candidate.prefix_location,
            density_adjusted_mass_dev: candidate.density_adjusted_mass_dev,
            density_adjusted_mass_ratio_n_pow_11_8: candidate.density_adjusted_mass_ratio_n_pow_11_8,
            joint_score: candidate.joint_score,
            is_prefix_winner,
            is_mass_winner,
            is_joint_winner,
            is_pareto_minimal,
        });
    }

    let distance_edge_total_count = candidates
        .len()
        .saturating_mul(candidates.len().saturating_sub(1))
        / 2;
    let should_export_distance_edges = distance_edge_cap
        .map(|edge_cap| distance_edge_total_count <= edge_cap)
        .unwrap_or(true);
    let distance_edge_status = if should_export_distance_edges {
        "EXPORTED_ALL".to_string()
    } else {
        "SKIP_DISTANCE_EDGES_TOO_LARGE".to_string()
    };
    let mut distance_edges = Vec::new();
    if should_export_distance_edges {
        for left in 0..candidates.len() {
            for right in (left + 1)..candidates.len() {
                distance_edges.push(GroundFaceDistanceEdge {
                    left_index: left,
                    right_index: right,
                    symmetric_difference_distance: symmetric_difference_distance(
                        &candidates[left].witness,
                        &candidates[right].witness,
                    ),
                });
            }
        }
    }

    let field_selection_split = !(prefix_winner_indices == mass_winner_indices
        && prefix_winner_indices == joint_winner_indices);

    Some(GroundFaceExport {
        status: "EXPORTED_ALL".to_string(),
        cap,
        exact_maximizer_count,
        exported_count: witnesses.len(),
        coordinate_space:
            "(prefix_ratio_n_pow_7_8, density_adjusted_mass_ratio_n_pow_11_8, joint_score)"
                .to_string(),
        exposed_winners: Some(GroundFaceExposedWinners {
            prefix_winner_indices,
            mass_winner_indices,
            joint_winner_indices,
            pareto_minimal_indices,
            field_selection_split,
        }),
        witnesses,
        distance_edge_status,
        distance_edge_cap,
        distance_edge_total_count,
        distance_edges,
    })
}

fn state_encoding_witness(state: &State) -> StateEncodingWitness {
    StateEncodingWitness {
        occupied_sites: state.occupied_sites(),
        occupied_suffix: state.occupied_suffix(OCCUPIED_SUFFIX_WIDTH),
        used_differences: bits_to_sites(state.used_differences_mask),
        cardinality: state.cardinality as usize,
    }
}

#[derive(Clone, Copy)]
struct PrefixJob {
    x: usize,
    state: State,
}

struct PrefixJobs {
    jobs: Vec<PrefixJob>,
    transitions_attempted: u64,
    legal_occupy_transitions: u64,
    rejected_occupy_transitions: u64,
}

fn generate_prefix_jobs(n: usize, depth: usize) -> PrefixJobs {
    struct PrefixGenCtx {
        jobs: Vec<PrefixJob>,
        transitions_attempted: u64,
        legal_occupy_transitions: u64,
        rejected_occupy_transitions: u64,
    }

    fn rec(n: usize, x: usize, depth_remaining: usize, state: State, ctx: &mut PrefixGenCtx) {
        if depth_remaining == 0 || x > n {
            ctx.jobs.push(PrefixJob { x, state });
            return;
        }
        rec(n, x + 1, depth_remaining - 1, state, ctx);
        ctx.transitions_attempted += 1;
        if let Some(new_state) = state.try_occupy(x) {
            ctx.legal_occupy_transitions += 1;
            rec(n, x + 1, depth_remaining - 1, new_state, ctx);
        } else {
            ctx.rejected_occupy_transitions += 1;
        }
    }

    let mut ctx = PrefixGenCtx {
        jobs: Vec::new(),
        transitions_attempted: 0,
        legal_occupy_transitions: 0,
        rejected_occupy_transitions: 0,
    };
    rec(n, 0, depth, State::empty(), &mut ctx);
    PrefixJobs {
        jobs: ctx.jobs,
        transitions_attempted: ctx.transitions_attempted,
        legal_occupy_transitions: ctx.legal_occupy_transitions,
        rejected_occupy_transitions: ctx.rejected_occupy_transitions,
    }
}

struct LayerPrunedAccum {
    layer_counts_map: HashMap<usize, u64>,
    terminal_state_count: u64,
    transitions_attempted: u64,
    legal_occupy_transitions: u64,
    rejected_occupy_transitions: u64,
    reachability_pruned_states: u64,
    over_target_pruned_states: u64,
    first_maximizer_witness: Vec<usize>,
    representative_ground_state_encoding: Option<StateEncodingWitness>,
    prefix_frontier: TopKAccumulator,
    mass_frontier: TopKAccumulator,
    joint_frontier: TopKAccumulator,
    ground_face_candidates: Vec<FrontierCandidate>,
    inheritance_probe: InheritanceProbeAccum,
    ground_face_masks: Vec<u128>,
}

impl LayerPrunedAccum {
    fn new(frontier_k: usize) -> Self {
        Self {
            layer_counts_map: HashMap::new(),
            terminal_state_count: 0,
            transitions_attempted: 0,
            legal_occupy_transitions: 0,
            rejected_occupy_transitions: 0,
            reachability_pruned_states: 0,
            over_target_pruned_states: 0,
            first_maximizer_witness: Vec::new(),
            representative_ground_state_encoding: None,
            prefix_frontier: TopKAccumulator::new(frontier_k),
            mass_frontier: TopKAccumulator::new(frontier_k),
            joint_frontier: TopKAccumulator::new(frontier_k),
            ground_face_candidates: Vec::new(),
            inheritance_probe: InheritanceProbeAccum::default(),
            ground_face_masks: Vec::new(),
        }
    }

    fn merge(&mut self, mut other: Self) {
        for (cardinality, count) in other.layer_counts_map {
            *self.layer_counts_map.entry(cardinality).or_insert(0) += count;
        }
        self.terminal_state_count += other.terminal_state_count;
        self.transitions_attempted += other.transitions_attempted;
        self.legal_occupy_transitions += other.legal_occupy_transitions;
        self.rejected_occupy_transitions += other.rejected_occupy_transitions;
        self.reachability_pruned_states += other.reachability_pruned_states;
        self.over_target_pruned_states += other.over_target_pruned_states;
        if !other.first_maximizer_witness.is_empty()
            && (self.first_maximizer_witness.is_empty()
                || other.first_maximizer_witness < self.first_maximizer_witness)
        {
            self.first_maximizer_witness = other.first_maximizer_witness;
            self.representative_ground_state_encoding =
                other.representative_ground_state_encoding.take();
        }
        self.prefix_frontier.merge(&other.prefix_frontier);
        self.mass_frontier.merge(&other.mass_frontier);
        self.joint_frontier.merge(&other.joint_frontier);
        self.ground_face_candidates
            .append(&mut other.ground_face_candidates);
        self.ground_face_masks.append(&mut other.ground_face_masks);
        self.inheritance_probe.exact_previous_persistence_count +=
            other.inheritance_probe.exact_previous_persistence_count;
        self.inheritance_probe.plus_one_previous_persistence_count +=
            other.inheritance_probe.plus_one_previous_persistence_count;
        self.inheritance_probe.inherited_union_present_count +=
            other.inheritance_probe.inherited_union_present_count;
    }
}

fn analyze_n(
    n: usize,
    frontier_k: usize,
    prune_deficiency: Option<usize>,
    parity: Option<(usize, u64)>,
    export_ground_face_cap: Option<usize>,
    export_ground_face_distance_edge_cap: Option<usize>,
    parallel_prefix_depth: usize,
    inheritance_source: Option<&InheritanceSource>,
    export_ground_face_masks: bool,
) -> Row {
    assert!(
        n < 128,
        "bitset transfer operator currently supports n < 128"
    );
    if prune_deficiency.is_some() {
        return analyze_n_layer_pruned(
            n,
            frontier_k,
            prune_deficiency,
            parity,
            export_ground_face_cap,
            export_ground_face_distance_edge_cap,
            parallel_prefix_depth,
            inheritance_source,
            export_ground_face_masks,
        );
    }
    let mut states: HashMap<State, u64> = HashMap::new();
    states.insert(State::empty(), 1);
    let mut peak_frontier_state_count = 1u64;
    let mut peak_frontier_after_site = 0usize;
    let mut transitions_attempted = 0u64;
    let mut legal_occupy_transitions = 0u64;
    let mut rejected_occupy_transitions = 0u64;
    let mut reachability_pruned_states = 0u64;
    let pruning_min_cardinality =
        prune_deficiency.and_then(|deficiency| parity.map(|(h, _)| h.saturating_sub(deficiency)));

    for x in 0..=n {
        let mut next = states.clone();
        for state in states.keys() {
            transitions_attempted += 1;
            match state.try_occupy(x) {
                Some(new_state) => {
                    legal_occupy_transitions += 1;
                    next.insert(new_state, 1);
                }
                None => {
                    rejected_occupy_transitions += 1;
                }
            }
        }
        if let Some(min_cardinality) = pruning_min_cardinality {
            let remaining_sites = n - x;
            let before = next.len();
            next.retain(|state, _| state.cardinality as usize + remaining_sites >= min_cardinality);
            reachability_pruned_states += (before - next.len()) as u64;
        }
        states = next;
        if states.len() as u64 > peak_frontier_state_count {
            peak_frontier_state_count = states.len() as u64;
            peak_frontier_after_site = x;
        }
    }

    let mut layer_counts_map: HashMap<usize, u64> = HashMap::new();
    let mut h_n = 0usize;
    for state in states.keys() {
        let k = state.cardinality as usize;
        *layer_counts_map.entry(k).or_insert(0) += 1;
        h_n = h_n.max(k);
    }

    let maximizer_count = *layer_counts_map.get(&h_n).unwrap_or(&0);
    let mut first_maximizer_witness = Vec::new();
    let mut representative_ground_state_encoding: Option<StateEncodingWitness> = None;
    let mut prefix_frontier = TopKAccumulator::new(frontier_k);
    let mut mass_frontier = TopKAccumulator::new(frontier_k);
    let mut joint_frontier = TopKAccumulator::new(frontier_k);
    let mut ground_face_candidates = Vec::new();
    let mut inheritance_probe = InheritanceProbeAccum::default();
    let mut ground_face_masks = Vec::new();

    for state in states
        .keys()
        .filter(|state| state.cardinality as usize == h_n)
    {
        let witness = state.occupied_sites();
        if first_maximizer_witness.is_empty() || witness < first_maximizer_witness {
            first_maximizer_witness = witness.clone();
            representative_ground_state_encoding = Some(state_encoding_witness(state));
        }
        let candidate = frontier_candidate(witness, n);
        if export_ground_face_cap.is_some() {
            ground_face_candidates.push(candidate.clone());
        }
        prefix_frontier.observe(FrontierCandidate {
            objective_value: candidate.prefix_residual,
            ..candidate.clone()
        });
        mass_frontier.observe(FrontierCandidate {
            objective_value: candidate.density_adjusted_mass_dev,
            ..candidate.clone()
        });
        joint_frontier.observe(candidate);
        if let Some(source) = inheritance_source {
            observe_inheritance(state.occupied_mask, source, &mut inheritance_probe);
        }
        if export_ground_face_masks {
            ground_face_masks.push(state.occupied_mask);
        }
    }

    let mut layer_counts: Vec<LayerCount> = layer_counts_map
        .iter()
        .map(|(cardinality, count)| LayerCount {
            cardinality: *cardinality,
            count: *count,
        })
        .collect();
    layer_counts.sort_by_key(|row| row.cardinality);

    let deficiency_layers: Vec<DeficiencyLayer> = (0..=3)
        .filter_map(|deficiency| {
            h_n.checked_sub(deficiency)
                .map(|cardinality| DeficiencyLayer {
                    deficiency,
                    cardinality,
                    count: *layer_counts_map.get(&cardinality).unwrap_or(&0),
                })
        })
        .collect();

    let first_excited_cardinality = (0..h_n)
        .rev()
        .find(|cardinality| *layer_counts_map.get(cardinality).unwrap_or(&0) > 0);
    let first_excited_count =
        first_excited_cardinality.map(|cardinality| *layer_counts_map.get(&cardinality).unwrap());
    let gap_to_first_excited_layer = first_excited_cardinality.map(|cardinality| h_n - cardinality);

    let parity_reference = parity.map(|(exact_h_n, exact_maximizer_count)| ParityReference {
        exact_h_n,
        exact_maximizer_count,
        h_n_matches: h_n == exact_h_n,
        maximizer_count_matches: maximizer_count == exact_maximizer_count,
    });

    Row {
        n,
        h_n,
        maximizer_count,
        first_maximizer_witness,
        representative_ground_state_encoding: representative_ground_state_encoding
            .expect("every nonempty scan has at least one ground state"),
        spectral_observables: SpectralObservables {
            ground_cardinality_h_n: h_n,
            ground_state_degeneracy: maximizer_count,
            ground_entropy_ln: (maximizer_count as f64).ln(),
            first_excited_cardinality,
            first_excited_count,
            gap_to_first_excited_layer,
            near_ground_count_h_minus_1: h_n
                .checked_sub(1)
                .map(|k| *layer_counts_map.get(&k).unwrap_or(&0))
                .unwrap_or(0),
            near_ground_count_h_minus_2: h_n
                .checked_sub(2)
                .map(|k| *layer_counts_map.get(&k).unwrap_or(&0))
                .unwrap_or(0),
            layer_counts,
            deficiency_layers,
        },
        transfer_diagnostics: TransferDiagnostics {
            traversal_mode: "full_state_hashmap".to_string(),
            sites_processed: n + 1,
            terminal_state_count: states.len() as u64,
            peak_frontier_state_count,
            peak_frontier_after_site,
            transitions_attempted,
            legal_occupy_transitions,
            rejected_occupy_transitions,
            prune_deficiency,
            pruning_target_h_n: parity.map(|(h, _)| h),
            pruning_min_cardinality,
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
        ground_face_export: build_ground_face_export(
            export_ground_face_cap,
            export_ground_face_distance_edge_cap,
            maximizer_count,
            ground_face_candidates,
        ),
        ground_face_mask_export: build_mask_export(
            export_ground_face_masks,
            maximizer_count,
            ground_face_masks,
        ),
        inheritance_probe: build_inheritance_probe(
            inheritance_source,
            &inheritance_probe,
            maximizer_count,
        ),
    }
}

fn analyze_n_layer_pruned(
    n: usize,
    frontier_k: usize,
    prune_deficiency: Option<usize>,
    parity: Option<(usize, u64)>,
    export_ground_face_cap: Option<usize>,
    export_ground_face_distance_edge_cap: Option<usize>,
    parallel_prefix_depth: usize,
    inheritance_source: Option<&InheritanceSource>,
    export_ground_face_masks: bool,
) -> Row {
    let (target_h, _) =
        parity.expect("layer-pruned traversal requires a reference h(n) and maximizer count");
    let deficiency = prune_deficiency.expect("layer-pruned traversal requires deficiency");
    let min_cardinality = target_h.saturating_sub(deficiency);

    struct DfsCtx<'a> {
        n: usize,
        target_h: usize,
        min_cardinality: usize,
        accum: &'a mut LayerPrunedAccum,
        export_ground_face_cap: Option<usize>,
        inheritance_source: Option<&'a InheritanceSource>,
        export_ground_face_masks: bool,
    }

    fn dfs(x: usize, state: State, ctx: &mut DfsCtx<'_>) {
        if state.cardinality as usize > ctx.target_h {
            ctx.accum.over_target_pruned_states += 1;
            return;
        }
        let remaining_sites = if x <= ctx.n { ctx.n - x + 1 } else { 0 };
        if state.cardinality as usize + remaining_sites < ctx.min_cardinality {
            ctx.accum.reachability_pruned_states += 1;
            return;
        }
        if x > ctx.n {
            let cardinality = state.cardinality as usize;
            if cardinality < ctx.min_cardinality {
                return;
            }
            ctx.accum.terminal_state_count += 1;
            *ctx.accum.layer_counts_map.entry(cardinality).or_insert(0) += 1;
            if cardinality == ctx.target_h {
                let witness = state.occupied_sites();
                if ctx.accum.first_maximizer_witness.is_empty()
                    || witness < ctx.accum.first_maximizer_witness
                {
                    ctx.accum.first_maximizer_witness = witness.clone();
                    ctx.accum.representative_ground_state_encoding =
                        Some(state_encoding_witness(&state));
                }
                let candidate = frontier_candidate(witness, ctx.n);
                if ctx.export_ground_face_cap.is_some() {
                    ctx.accum.ground_face_candidates.push(candidate.clone());
                }
                ctx.accum.prefix_frontier.observe(FrontierCandidate {
                    objective_value: candidate.prefix_residual,
                    ..candidate.clone()
                });
                ctx.accum.mass_frontier.observe(FrontierCandidate {
                    objective_value: candidate.density_adjusted_mass_dev,
                    ..candidate.clone()
                });
                ctx.accum.joint_frontier.observe(candidate);
                if let Some(source) = ctx.inheritance_source {
                    observe_inheritance(state.occupied_mask, source, &mut ctx.accum.inheritance_probe);
                }
                if ctx.export_ground_face_masks {
                    ctx.accum.ground_face_masks.push(state.occupied_mask);
                }
            }
            return;
        }

        dfs(x + 1, state, ctx);
        ctx.accum.transitions_attempted += 1;
        match state.try_occupy(x) {
            Some(new_state) => {
                ctx.accum.legal_occupy_transitions += 1;
                dfs(x + 1, new_state, ctx);
            }
            None => {
                ctx.accum.rejected_occupy_transitions += 1;
            }
        }
    }

    let mut accum = if parallel_prefix_depth == 0 {
        let mut accum = LayerPrunedAccum::new(frontier_k);
        let mut ctx = DfsCtx {
            n,
            target_h,
            min_cardinality,
            accum: &mut accum,
            export_ground_face_cap,
            inheritance_source,
            export_ground_face_masks,
        };
        dfs(0, State::empty(), &mut ctx);
        accum
    } else {
        let prefix_jobs = generate_prefix_jobs(n, parallel_prefix_depth);
        let prefix_transitions_attempted = prefix_jobs.transitions_attempted;
        let prefix_legal_occupy_transitions = prefix_jobs.legal_occupy_transitions;
        let prefix_rejected_occupy_transitions = prefix_jobs.rejected_occupy_transitions;
        let workers: Vec<LayerPrunedAccum> = prefix_jobs
            .jobs
            .into_par_iter()
            .map(|job| {
                let mut worker = LayerPrunedAccum::new(frontier_k);
                let mut ctx = DfsCtx {
                    n,
                    target_h,
                    min_cardinality,
                    accum: &mut worker,
                    export_ground_face_cap,
                    inheritance_source,
                    export_ground_face_masks,
                };
                dfs(job.x, job.state, &mut ctx);
                worker
            })
            .collect();
        let mut accum = LayerPrunedAccum::new(frontier_k);
        for worker in workers {
            accum.merge(worker);
        }
        accum.transitions_attempted += prefix_transitions_attempted;
        accum.legal_occupy_transitions += prefix_legal_occupy_transitions;
        accum.rejected_occupy_transitions += prefix_rejected_occupy_transitions;
        accum
    };
    accum.reachability_pruned_states += accum.over_target_pruned_states;

    let h_n = *accum.layer_counts_map.keys().max().unwrap_or(&0);
    let maximizer_count = *accum.layer_counts_map.get(&target_h).unwrap_or(&0);
    let mut layer_counts: Vec<LayerCount> = accum
        .layer_counts_map
        .iter()
        .map(|(cardinality, count)| LayerCount {
            cardinality: *cardinality,
            count: *count,
        })
        .collect();
    layer_counts.sort_by_key(|row| row.cardinality);

    let deficiency_layers: Vec<DeficiencyLayer> = (0..=deficiency.min(3))
        .filter_map(|deficiency| {
            target_h
                .checked_sub(deficiency)
                .map(|cardinality| DeficiencyLayer {
                    deficiency,
                    cardinality,
                    count: *accum.layer_counts_map.get(&cardinality).unwrap_or(&0),
                })
        })
        .collect();
    let first_excited_cardinality = (0..target_h)
        .rev()
        .find(|cardinality| *accum.layer_counts_map.get(cardinality).unwrap_or(&0) > 0);
    let first_excited_count = first_excited_cardinality
        .map(|cardinality| *accum.layer_counts_map.get(&cardinality).unwrap());
    let gap_to_first_excited_layer =
        first_excited_cardinality.map(|cardinality| target_h - cardinality);
    let parity_reference = parity.map(|(exact_h_n, exact_maximizer_count)| ParityReference {
        exact_h_n,
        exact_maximizer_count,
        h_n_matches: h_n == exact_h_n,
        maximizer_count_matches: maximizer_count == exact_maximizer_count,
    });

    Row {
        n,
        h_n,
        maximizer_count,
        first_maximizer_witness: accum.first_maximizer_witness,
        representative_ground_state_encoding: accum
            .representative_ground_state_encoding
            .expect("layer-pruned traversal should retain at least one ground state"),
        spectral_observables: SpectralObservables {
            ground_cardinality_h_n: h_n,
            ground_state_degeneracy: maximizer_count,
            ground_entropy_ln: (maximizer_count as f64).ln(),
            first_excited_cardinality,
            first_excited_count,
            gap_to_first_excited_layer,
            near_ground_count_h_minus_1: target_h
                .checked_sub(1)
                .map(|k| *accum.layer_counts_map.get(&k).unwrap_or(&0))
                .unwrap_or(0),
            near_ground_count_h_minus_2: target_h
                .checked_sub(2)
                .map(|k| *accum.layer_counts_map.get(&k).unwrap_or(&0))
                .unwrap_or(0),
            layer_counts,
            deficiency_layers,
        },
        transfer_diagnostics: TransferDiagnostics {
            traversal_mode: if parallel_prefix_depth == 0 {
                "layer_pruned_dfs".to_string()
            } else {
                format!("layer_pruned_dfs_parallel_prefix_depth_{}", parallel_prefix_depth)
            },
            sites_processed: n + 1,
            terminal_state_count: accum.terminal_state_count,
            peak_frontier_state_count: 0,
            peak_frontier_after_site: 0,
            transitions_attempted: accum.transitions_attempted,
            legal_occupy_transitions: accum.legal_occupy_transitions,
            rejected_occupy_transitions: accum.rejected_occupy_transitions,
            prune_deficiency,
            pruning_target_h_n: Some(target_h),
            pruning_min_cardinality: Some(min_cardinality),
            reachability_pruned_states: accum.reachability_pruned_states,
        },
        parity_reference,
        top_k_frontier: TopKFrontier {
            k: frontier_k,
            zero_temperature_field_interpretation:
                "prefix_best, density_adjusted_mass_best, and joint_best are zero-temperature ground-state selections under small positive prefix/mass fields after cardinality is fixed at h(n)"
                    .to_string(),
            prefix_best: accum.prefix_frontier.ranked(),
            density_adjusted_mass_best: accum.mass_frontier.ranked(),
            joint_best: accum.joint_frontier.ranked(),
        },
        ground_face_export: build_ground_face_export(
            export_ground_face_cap,
            export_ground_face_distance_edge_cap,
            maximizer_count,
            accum.ground_face_candidates,
        ),
        ground_face_mask_export: build_mask_export(
            export_ground_face_masks,
            maximizer_count,
            accum.ground_face_masks,
        ),
        inheritance_probe: build_inheritance_probe(
            inheritance_source,
            &accum.inheritance_probe,
            maximizer_count,
        ),
    }
}

fn reference_hn_counts() -> HashMap<usize, (usize, u64)> {
    // Current local truth sources:
    // - EXP-MM-030-EXACT-REFERENCE-10-30-2026-04-29 for 10..30.
    // - EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-56-60-2026-04-24 for 56..60.
    // - EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-61-65-2026-04-24 for 61..65.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-51-70-2026-04-29 for 69..70.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-71-2026-04-29 for 71.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-72-2026-05-01 for 72.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-73-2026-05-01 for 73.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-74-2026-05-01 for 74.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-75-LB11-2026-05-01 for 75.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-76-LB11-PAR8-2026-05-01 for 76.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-77-LB11-PAR8-2026-05-01 for 77.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-78-LB11-PAR8-2026-05-01 for 78.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-79-LB11-PAR8-2026-05-01 for 79.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-80-LB11-PAR8-2026-05-01 for 80.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-81-LB11-PAR8-2026-05-01 for 81.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-82-LB11-PAR8-2026-05-01 for 82.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-83-LB11-PAR8-2026-05-01 for 83.
    // - EXP-MM-030-RUST-TOPK-FRONTIER-84-LB11-PAR8-2026-05-01 for 84.
    // Domain convention is the exact scanner's 0..=n lattice, not the shifted 1..n table.
    HashMap::from([
        (10, (4, 110)),
        (11, (5, 4)),
        (12, (5, 22)),
        (13, (5, 68)),
        (14, (5, 156)),
        (15, (5, 320)),
        (16, (5, 584)),
        (17, (6, 8)),
        (18, (6, 24)),
        (19, (6, 80)),
        (20, (6, 206)),
        (21, (6, 504)),
        (22, (6, 1004)),
        (23, (6, 1910)),
        (24, (6, 3380)),
        (25, (7, 10)),
        (26, (7, 34)),
        (27, (7, 98)),
        (28, (7, 282)),
        (29, (7, 760)),
        (30, (7, 1618)),
        (56, (10, 4)),
        (57, (10, 6)),
        (58, (10, 10)),
        (59, (10, 18)),
        (60, (10, 54)),
        (61, (10, 152)),
        (62, (10, 398)),
        (63, (10, 1022)),
        (64, (10, 2360)),
        (65, (10, 5018)),
        (66, (10, 9994)),
        (67, (10, 19418)),
        (68, (10, 36234)),
        (69, (10, 66412)),
        (70, (10, 117202)),
        (71, (10, 203840)),
        (72, (11, 4)),
        (73, (11, 8)),
        (74, (11, 34)),
        (75, (11, 84)),
        (76, (11, 214)),
        (77, (11, 482)),
        (78, (11, 970)),
        (79, (11, 1974)),
        (80, (11, 4030)),
        (81, (11, 8214)),
        (82, (11, 15958)),
        (83, (11, 30510)),
        (84, (11, 56110)),
    ])
}

fn summarize_parity(rows: &[Row]) -> ParitySummary {
    let mut checked_count = 0usize;
    let mut h_n_match_count = 0usize;
    let mut maximizer_count_match_count = 0usize;
    let mut mismatch_ns = Vec::new();

    for row in rows {
        if let Some(reference) = &row.parity_reference {
            checked_count += 1;
            if reference.h_n_matches {
                h_n_match_count += 1;
            }
            if reference.maximizer_count_matches {
                maximizer_count_match_count += 1;
            }
            if !reference.h_n_matches || !reference.maximizer_count_matches {
                mismatch_ns.push(row.n);
            }
        }
    }

    ParitySummary {
        checked_count,
        h_n_match_count,
        maximizer_count_match_count,
        mismatch_ns,
    }
}

fn state_model() -> StateModel {
    StateModel {
        representation: "State = { occupied_mask: u128, used_differences_mask: u128, cardinality: u8 }"
            .to_string(),
        occupied_mask_bits: "bit i is 1 iff lattice site i is occupied".to_string(),
        used_differences_mask_bits:
            "bit d is 1 iff a positive difference d has already been realized".to_string(),
        occupied_suffix_width: OCCUPIED_SUFFIX_WIDTH,
        transition_rule:
            "skip x always; occupy x iff every new difference |x-a| is absent from used_differences_mask"
                .to_string(),
    }
}

fn build_report(results: &Results) -> String {
    let mut lines = Vec::new();
    lines.push(format!(
        "# {} — PMF Transfer-Operator Parity Scan",
        results.experiment_id
    ));
    lines.push(String::new());
    lines.push("## Identification".to_string());
    lines.push(String::new());
    lines.push("| Field | Value |".to_string());
    lines.push("|---|---|".to_string());
    lines.push(format!("| Experiment ID | {} |", results.experiment_id));
    lines.push("| Erdős Problem | #30 — finite Sidon set rigidity |".to_string());
    lines.push(
        "| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |"
            .to_string(),
    );
    lines.push(format!(
        "| Scan window | n = {} through n = {} |",
        results.n_min, results.n_max
    ));
    lines.push(format!("| Frontier k | {} |", results.frontier_k));
    lines.push(format!(
        "| Prune deficiency | {} |",
        results
            .prune_deficiency
            .map(|d| d.to_string())
            .unwrap_or_else(|| "none".to_string())
    ));
    lines.push(format!(
        "| Ground-face export cap | {} |",
        results
            .export_ground_face_cap
            .map(|cap| cap.to_string())
            .unwrap_or_else(|| "none".to_string())
    ));
    lines.push(format!(
        "| Ground-face distance-edge cap | {} |",
        results
            .export_ground_face_distance_edge_cap
            .map(|cap| cap.to_string())
            .unwrap_or_else(|| "none".to_string())
    ));
    lines.push(format!(
        "| Inheritance source | {} |",
        results
            .inheritance_source_results
            .as_deref()
            .unwrap_or("none")
    ));
    lines.push(String::new());
    lines.push("## State Model".to_string());
    lines.push(String::new());
    lines.push(format!("- `{}`", results.state_model.representation));
    lines.push(format!(
        "- Occupied mask: {}",
        results.state_model.occupied_mask_bits
    ));
    lines.push(format!(
        "- Difference memory: {}",
        results.state_model.used_differences_mask_bits
    ));
    lines.push(format!(
        "- Occupied suffix: last {} occupied sites are serialized for representative ground states.",
        results.state_model.occupied_suffix_width
    ));
    lines.push(format!(
        "- Transition: {}",
        results.state_model.transition_rule
    ));
    lines.push(String::new());
    lines.push("## Reachability Pruning".to_string());
    lines.push(String::new());
    if let Some(deficiency) = results.prune_deficiency {
        lines.push(format!(
            "Enabled with deficiency `{}`. After each site, states are retained only if their current cardinality plus remaining sites can still reach `h(n)-{}` using the reference `h(n)` for that row.",
            deficiency, deficiency
        ));
    } else {
        lines.push("Disabled. All exact transfer states are retained.".to_string());
    }
    lines.push(String::new());
    lines.push("| n | min cardinality | pruned states | terminal retained states |".to_string());
    lines.push("|---|---:|---:|---:|".to_string());
    for row in &results.rows {
        lines.push(format!(
            "| {} | {} | {} | {} |",
            row.n,
            row.transfer_diagnostics
                .pruning_min_cardinality
                .map(|v| v.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            row.transfer_diagnostics.reachability_pruned_states,
            row.transfer_diagnostics.terminal_state_count
        ));
    }
    lines.push(String::new());
    lines.push("## Parity Gate".to_string());
    lines.push(String::new());
    lines.push(format!(
        "Checked {} n-values. h(n) matched in {}. Maximizer counts matched in {}. Mismatches: {:?}.",
        results.parity_summary.checked_count,
        results.parity_summary.h_n_match_count,
        results.parity_summary.maximizer_count_match_count,
        results.parity_summary.mismatch_ns
    ));
    lines.push(String::new());
    lines.push("## Spectral / Ground-State Summary".to_string());
    lines.push(String::new());
    lines.push("| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |".to_string());
    lines.push("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|".to_string());
    for row in &results.rows {
        let joint_score = row
            .top_k_frontier
            .joint_best
            .first()
            .map(|w| w.joint_score)
            .unwrap_or(f64::NAN);
        lines.push(format!(
            "| {} | {} | {} | {:.6} | {} | {} | {} | {} | {} | {:.6} |",
            row.n,
            row.h_n,
            row.maximizer_count,
            row.spectral_observables.ground_entropy_ln,
            row.spectral_observables.near_ground_count_h_minus_1,
            row.spectral_observables.near_ground_count_h_minus_2,
            row.spectral_observables
                .gap_to_first_excited_layer
                .map(|v| v.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            row.transfer_diagnostics.terminal_state_count,
            row.transfer_diagnostics.peak_frontier_state_count,
            joint_score
        ));
    }
    lines.push(String::new());
    lines.push("## Zero-Temperature Field Tilts".to_string());
    lines.push(String::new());
    lines.push(
        "After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.".to_string(),
    );
    lines.push(String::new());
    lines.push(
        "| n | prefix-field witness | mass-field witness | joint-field witness | joint score |"
            .to_string(),
    );
    lines.push("|---|---|---|---|---:|".to_string());
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
        let joint_score = row
            .top_k_frontier
            .joint_best
            .first()
            .map(|w| w.joint_score)
            .unwrap_or(f64::NAN);
        lines.push(format!(
            "| {} | `{}` | `{}` | `{}` | {:.6} |",
            row.n, prefix, mass, joint, joint_score
        ));
    }
    if results.export_ground_face_cap.is_some() {
        lines.push(String::new());
        lines.push("## Small Ground-Face Export".to_string());
        lines.push(String::new());
        lines.push("| n | export status | exact maximizers | exported | edge status | edge total | field split | pareto minima |".to_string());
        lines.push("|---|---|---:|---:|---|---:|---|---|".to_string());
        for row in &results.rows {
            if let Some(export) = &row.ground_face_export {
                let (field_split, pareto_minima) = export
                    .exposed_winners
                    .as_ref()
                    .map(|winners| {
                        (
                            winners.field_selection_split.to_string(),
                            format!("{:?}", winners.pareto_minimal_indices),
                        )
                    })
                    .unwrap_or_else(|| ("NA".to_string(), "[]".to_string()));
                lines.push(format!(
                    "| {} | {} | {} | {} | {} | {} | {} | `{}` |",
                    row.n,
                    export.status,
                    export.exact_maximizer_count,
                    export.exported_count,
                    export.distance_edge_status,
                    export.distance_edge_total_count,
                    field_split,
                    pareto_minima
                ));
            }
        }
    }
    if results.rows.iter().any(|row| row.inheritance_probe.is_some()) {
        lines.push(String::new());
        lines.push("## Plateau Inheritance Probe".to_string());
        lines.push(String::new());
        lines.push("This light certificate checks whether the current exact face contains the previous exact face and the previous face shifted by `+1`, without requiring pairwise distance-edge export.".to_string());
        lines.push(String::new());
        lines.push("| n | source n | previous face | inherited union | prev exact present | prev +1 present | inherited present | new face | union contained |".to_string());
        lines.push("|---:|---:|---:|---:|---:|---:|---:|---:|---|".to_string());
        for row in &results.rows {
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
    lines.push("## Interpretation".to_string());
    lines.push(String::new());
    lines.push(
        "This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language."
            .to_string(),
    );
    lines.push(String::new());
    lines.push(
        "The theorem language remains blocked: these packets are finite parity and near-ground evidence, not a Sidon proof."
            .to_string(),
    );
    lines.push(String::new());
    lines.push("## Artifacts".to_string());
    lines.push(String::new());
    lines.push("| File | Type |".to_string());
    lines.push("|---|---|".to_string());
    lines.push(format!(
        "| {}_RESULTS.json | Structured transfer-operator results |",
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
    if cli.n_max >= 128 {
        return Err("u128 bitset transfer operator currently requires --n-max < 128".into());
    }

    let experiment_id = cli.experiment_id.unwrap_or_else(|| {
        format!(
            "EXP-MM-030-PMF-TRANSFER-PARITY-{}-{}-2026-04-29",
            cli.n_min, cli.n_max
        )
    });

    let references = reference_hn_counts();
    let inheritance_source = cli
        .inheritance_source_results
        .as_deref()
        .map(load_inheritance_source)
        .transpose()?;
    let t0 = Instant::now();
    let rows: Vec<Row> = (cli.n_min..=cli.n_max)
        .map(|n| {
            analyze_n(
                n,
                cli.frontier_k,
                cli.prune_deficiency,
                references.get(&n).copied(),
                cli.export_ground_face_cap,
                cli.export_ground_face_distance_edge_cap,
                cli.parallel_prefix_depth,
                inheritance_source.as_ref(),
                cli.export_ground_face_masks,
            )
        })
        .collect();
    let parity_summary = summarize_parity(&rows);

    let results = Results {
        experiment_id,
        date: "2026-04-30".to_string(),
        erdos_problem: 30,
        scan_mode: "pmf_transfer_operator_parity".to_string(),
        implementation: "Rust bitset transfer-state enumerator".to_string(),
        state_model: state_model(),
        n_min: cli.n_min,
        n_max: cli.n_max,
        frontier_k: cli.frontier_k,
        prune_deficiency: cli.prune_deficiency,
        export_ground_face_cap: cli.export_ground_face_cap,
        export_ground_face_distance_edge_cap: cli.export_ground_face_distance_edge_cap,
        export_ground_face_masks: cli.export_ground_face_masks,
        inheritance_source_results: inheritance_source
            .as_ref()
            .map(|source| source.results_path.clone()),
        total_runtime_sec: t0.elapsed().as_secs_f64(),
        rows,
        parity_summary,
    };

    write_outputs(&results, &cli.output_dir)?;
    println!(
        "wrote packet {} to {}",
        results.experiment_id,
        cli.output_dir.display()
    );
    Ok(())
}
