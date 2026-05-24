# Face/Field Response Proof-Target Packet

Date: 2026-05-01  
Status: internal proof-target packet; not a theorem claim

## Feynman-Plain View

A single optimizer is a photograph. The extremal face is the mechanism. The
Collider's useful move is that it stops asking which witness came first and
starts asking how the whole exact face responds when small fields touch it.

The next proof lane is therefore narrow: Sidon, B_2[g], and a quiet sum-free
control. The purpose is to extract one clean invariant that can be stated
without Collider vocabulary.

## Claim Boundary

This packet is derived from existing exact finite artifacts. It does not claim a
solution, proof, public theorem, SOTA theorem progress, complete Erdos #30 Lean
proof, or D1 status update. It proposes proof targets.

## Headline Facts

| System | Window | Field split | h-jump rows | Boundary |
|---|---:|---:|---|---|
| Sidon #30 | 20..30 | 11/11 | [25] | exact ground and near-ground counts |
| B_2[2] #755 | 20..30 | 11/11 | [21, 26, 30] | exact ground and h-1 counts |
| Sum-free #166 | 20..30 | 0/11 | [21, 23, 25, 27, 29] | exact negative control |
| B_2[3] #755 | 43..47 | 4/5 | [45] | ground-only reset window |

## Formal Foothold

`PT-1` has been climbed from prose to an abstract Lean lemma, one concrete
#30 row certificate, and then the full Atheneum Sidon `n=20..30` exported-pair
window certificate:

- Module: `Erdos30_FaceField`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField`
- Last verified: 2026-05-01
- Theorems: `Erdos.Collider.fieldSplit_card_two_le`,
  `Erdos.Collider.fieldSplit_not_card_le_one`,
  `Erdos.Collider.isFieldMinOn_of_eq_on`,
  `Erdos.Collider.isFieldMinOn_weightedJoint_of_prefix_eq_on`

Boundary: This compiles the abstract finite exposed-face lemmas and face-local minimizer transfer lemmas only; it does not prove Sidon or B_2[g] field-response structure.

Exact observable helper:

- Module: `Erdos30_FaceField_ExactObservables`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_ExactObservables.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField_ExactObservables`
- Definitions include:
  `Erdos30FaceFieldExactObservables.prefixCount`,
  `Erdos30FaceFieldExactObservables.densityAdjustedMassTwice`

Boundary: This defines exact prefix-count scaffolding, the exact integer density-adjusted mass observable, and a finite scalar full-prefix residual observable as the supremum of prefixResidualAt over t <= n. It proves the scalar is zero from pointwise full-prefix residual zero, positive from one positive cutoff, nonnegative for all witnesses, and equal to a target value from a pointwise upper bound plus one attaining cutoff. It is still a finite observable helper, not a theorem about Erdos #30.

N=30 certificate:

- Module: `Erdos30_FaceField_N30_Certificate`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_N30_Certificate.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField_N30_Certificate`
- Theorems include:
  `Erdos30FaceFieldN30Certificate.n30_prefix_mass_field_split`,
  `Erdos30FaceFieldN30Certificate.n30_exposed_family_has_at_least_two_witnesses`

Boundary: This certifies three packet-exported n=30 Sidon witnesses, their field split over the certified witness family, and the abstract exposed-face consequence. It is not a proof about the full n=30 face of 1618 states and not theorem progress on Erdos #30.

Sidon `n=20..30` window certificate:

- Module: `Erdos30_FaceField_Window_Certificate`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_Window_Certificate.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField_Window_Certificate`
- Theorems include:
  `Erdos30FaceFieldWindowCertificate.sidon_window_20_30_all_prefix_mass_split`,
  `Erdos30FaceFieldWindowCertificate.sidon_window_20_30_all_exposed_pairs_have_two_witnesses`

Boundary: This certifies packet-exported prefix/mass Sidon witnesses for every Atheneum Sidon row n=20..30, their pairwise field splits over two-point certified witness families, and the abstract exposed-face consequence. It is not a proof about the full extremal face in any row and not theorem progress on Erdos #30.

Sidon `n=57/58` branch certificate:

- Module: `Erdos30_FaceField_57_58_Certificate`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_57_58_Certificate.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField_57_58_Certificate`
- Theorems include:
  `Erdos30FaceField5758Certificate.n57_mass_joint_field_split`,
  `Erdos30FaceField5758Certificate.n58_prefix_mass_field_split`,
  `Erdos30FaceField5758Certificate.finite_57_58_branch_field_response_certificate`

Boundary: This certifies finite relationships among exported n=57/n=58 branch witnesses: Sidon/range/card checks, n=57 mass/joint split inside a shared positive-difference skeleton, n=58 prefix/mass split across a new skeleton branch, and Pareto translation persistence. It is not a full extremal enumeration proof, not an asymptotic theorem, and not theorem progress on Erdos #30.

Sidon `n=56..58` full exported-face certificate:

- Module: `Erdos30_FaceField_56_58_FullFace_Certificate`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_56_58_FullFace_Certificate.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField_56_58_FullFace_Certificate`
- Theorems include:
  `Erdos30FaceField5658FullFaceCertificate.n56_exported_face_card`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_winner_table_matches_packet`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_split_pattern_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_exact_mass_winner_table_matches_packet`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_exact_mass_split_pattern_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_exact_prefix_probe_winner_table_matches_packet`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_exact_prefix_probe_strict_nonwinner_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_full_prefix_zero_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_packet_probe_separation_is_actual_residual_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_full_prefix_zero_and_probe_separation_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_full_prefix_scalar_zero_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_full_prefix_scalar_positive_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_full_prefix_scalar_min_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_full_prefix_scalar_exact_mass_split_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_exact_prefix_probe_exact_mass_split_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.n56W3_prefix_count_table_matches_witness`,
  `Erdos30FaceField5658FullFaceCertificate.n56W3_full_prefix_segment_bound_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.n56W3_full_prefix_residual_zero`,
  `Erdos30FaceField5658FullFaceCertificate.n57W4_prefix_count_table_matches_witness`,
  `Erdos30FaceField5658FullFaceCertificate.n57W4_full_prefix_segment_bound_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.n57W4_full_prefix_residual_zero`,
  `Erdos30FaceField5658FullFaceCertificate.n57W5_prefix_count_table_matches_witness`,
  `Erdos30FaceField5658FullFaceCertificate.n57W5_full_prefix_segment_bound_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.n57W5_full_prefix_residual_zero`,
  `Erdos30FaceField5658FullFaceCertificate.n58W2_prefix_count_table_matches_witness`,
  `Erdos30FaceField5658FullFaceCertificate.n58W2_full_prefix_segment_bound_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.n58W2_full_prefix_residual_zero`,
  `Erdos30FaceField5658FullFaceCertificate.n58W8_prefix_count_table_matches_witness`,
  `Erdos30FaceField5658FullFaceCertificate.n58W8_full_prefix_segment_bound_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.n58W8_full_prefix_residual_zero`,
  `Erdos30FaceField5658FullFaceCertificate.n58W9_prefix_count_table_matches_witness`,
  `Erdos30FaceField5658FullFaceCertificate.n58W9_full_prefix_segment_bound_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.n58W9_full_prefix_residual_zero`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_exact_joint_winner_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_full_prefix_scalar_matches_packet_probe_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_scalar_joint_matches_probe_joint_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_scalar_full_prefix_joint_winner_certificate`,
  `Erdos30FaceField5658FullFaceCertificate.n56_full_prefix_residual_max_matches_exact_prefix_probe_on_face`,
  `Erdos30FaceField5658FullFaceCertificate.n57_full_prefix_residual_max_matches_exact_prefix_probe_on_face`,
  `Erdos30FaceField5658FullFaceCertificate.n58_full_prefix_residual_max_matches_exact_prefix_probe_on_face`

Boundary: This imports every exported ground-face witness for n=56,57,58, checks Sidon/range/card facts, recomputes packet winner sets over rank-coded observables, proves the density-adjusted mass winners from the exact integer formula |2*sum(A)-n*(|A|+1)|, proves exact zero prefix-probe stress for the packet-selected prefix winners, proves strict positive exact-prefix probe gaps for every non-winner, proves the exact joint-key winners from exact prefix probe plus exact integer mass, proves full prefix-count segment bounds against terminal drift, proves full-prefix residual zero for every packet-selected prefix winner over every cutoff t <= n, proves each packet probe non-winner has positive actual prefix residual at its packet cutoff, proves scalar full-prefix residual maxima equal the exact packet prefix probes for every exported witness in n=56..58, packages those equalities as face-local equality theorems, and uses the reusable weighted-joint minimizer-transfer lemma to prove scalar full-prefix joint winners W3/W5/W7. Exhaustiveness is inherited from the source exact packet; this is not an asymptotic theorem and not theorem progress on Erdos #30.

Sidon `n=59` scalar full-prefix joint microcertificate:

- Module: `Erdos30_FaceField_59_FullFace_Certificate`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_59_FullFace_Certificate.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField_59_FullFace_Certificate`
- Theorems include:
  `Erdos30FaceField59FullFaceCertificate.n59_micro_winner_table_matches_packet`,
  `Erdos30FaceField59FullFaceCertificate.n59_micro_full_prefix_zero_certificate`,
  `Erdos30FaceField59FullFaceCertificate.n59_micro_full_prefix_scalar_matches_packet_probe_certificate`,
  `Erdos30FaceField59FullFaceCertificate.n59_micro_scalar_joint_matches_probe_joint_certificate`,
  `Erdos30FaceField59FullFaceCertificate.n59_micro_scalar_full_prefix_joint_winner_certificate`,
  `Erdos30FaceField59FullFaceCertificate.n59_full_prefix_residual_max_matches_exact_prefix_probe_on_face`

Boundary: This imports all 18 exported n=59 ground-face witnesses, checks Sidon/range/card facts, proves exact integer mass winners, proves terminal full-prefix residual zero for the six packet prefix winners, proves positive packet-cutoff residuals for the remaining witnesses, proves scalar full-prefix residual maxima equal exact packet prefix probes on the n=59 face, and uses the reusable weighted-joint minimizer-transfer lemma to prove scalar full-prefix joint winner W7. Exhaustiveness is inherited from the source exact packet; this is not an asymptotic theorem and not theorem progress on Erdos #30.

Sidon `n=60` scalar full-prefix joint microcertificate:

- Module: `Erdos30_FaceField_60_FullFace_Certificate`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_60_FullFace_Certificate.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField_60_FullFace_Certificate`
- Theorems include:
  `Erdos30FaceField60FullFaceCertificate.n60_micro_winner_table_matches_packet`,
  `Erdos30FaceField60FullFaceCertificate.n60_micro_full_prefix_zero_certificate`,
  `Erdos30FaceField60FullFaceCertificate.n60_micro_full_prefix_scalar_matches_packet_probe_certificate`,
  `Erdos30FaceField60FullFaceCertificate.n60_micro_scalar_joint_matches_probe_joint_certificate`,
  `Erdos30FaceField60FullFaceCertificate.n60_micro_scalar_full_prefix_joint_winner_certificate`,
  `Erdos30FaceField60FullFaceCertificate.n60_full_prefix_residual_max_matches_exact_prefix_probe_on_face`

Boundary: This imports all 54 exported n=60 ground-face witnesses, checks Sidon/range/card facts, proves exact integer mass winner W43, proves terminal full-prefix residual zero for the twenty-one packet prefix winners, proves positive packet-cutoff residuals for the remaining witnesses, proves scalar full-prefix residual maxima equal exact packet prefix probes on the n=60 face, and uses the reusable weighted-joint minimizer-transfer lemma to prove scalar full-prefix joint winner W43. Exhaustiveness is inherited from the source exact packet; this is not an asymptotic theorem and not theorem progress on Erdos #30.

Sidon `n=61` compact joint-surface certificate:

- Module: `Erdos30_FaceField_61_JointSurface_Certificate`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_61_JointSurface_Certificate.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField_61_JointSurface_Certificate`
- Theorems include:
  `Erdos30FaceField61JointSurfaceCertificate.n61_face_index_count`,
  `Erdos30FaceField61JointSurfaceCertificate.n61_joint_surface_table_matches_packet`,
  `Erdos30FaceField61JointSurfaceCertificate.n61_joint_surface_is_proper_exported_subface`,
  `Erdos30FaceField61JointSurfaceCertificate.n61_joint_surface_minimizer_certificate`,
  `Erdos30FaceField61JointSurfaceCertificate.n61_joint_selected_surface_certificate`

Boundary: This does not import witness literals or independently prove Sidon enumeration. It proves the finite selection logic over the packet-backed n=61 exported index face: the joint-score-derived rank table has exactly the tied minimizer surface [85, 96, 113, 139] over Finset.range 152, and that surface satisfies the abstract IsFieldMinimizerSet predicate. Exact face exhaustiveness and witness correctness are inherited from the source exact packet; this is not theorem progress on Erdos #30.

Sidon `n=61..64` compact joint-surface certificate:

- Module: `Erdos30_FaceField_61_64_JointSurface_Certificate`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_61_64_JointSurface_Certificate.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField_61_64_JointSurface_Certificate`
- Theorems include:
  `Erdos30FaceField6164JointSurfaceCertificate.n62_joint_surface_table_matches_packet`,
  `Erdos30FaceField6164JointSurfaceCertificate.n62_joint_surface_minimizer_certificate`,
  `Erdos30FaceField6164JointSurfaceCertificate.n63_joint_surface_table_matches_packet`,
  `Erdos30FaceField6164JointSurfaceCertificate.n63_joint_surface_minimizer_certificate`,
  `Erdos30FaceField6164JointSurfaceCertificate.n64_joint_surface_table_matches_packet`,
  `Erdos30FaceField6164JointSurfaceCertificate.n64_joint_surface_minimizer_certificate`

Boundary: This does not import witness literals or independently prove Sidon enumeration. It proves finite selection logic over the packet-backed exported index faces for n=61,62,63,64: each joint-score-derived rank table has exactly the packet's tied minimizer surface, with surface sizes 4,6,10,11 over exported faces of size 152,398,1022,2360. Exact face exhaustiveness and witness correctness are inherited from the source exact packet; this is not theorem progress on Erdos #30.

Sidon `n=61..64` selected-surface transition certificate:

- Module: `Erdos30_FaceField_61_64_JointSurface_Transition_Certificate`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_61_64_JointSurface_Transition_Certificate.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField_61_64_JointSurface_Transition_Certificate`
- Theorems include:
  `Erdos30FaceField6164JointSurfaceTransitionCertificate.joint_surface_transition_count_table`,
  `Erdos30FaceField6164JointSurfaceTransitionCertificate.plus_one_previous_surface_count_table`,
  `Erdos30FaceField6164JointSurfaceTransitionCertificate.new_surface_count_table`,
  `Erdos30FaceField6164JointSurfaceTransitionCertificate.all_surface_transition_counts_balance`,
  `Erdos30FaceField6164JointSurfaceTransitionCertificate.all_surface_transitions_are_mixed`,
  `Erdos30FaceField6164JointSurfaceTransitionCertificate.joint_surface_transition_certificate_passes`

Boundary: This does not import full witness faces or independently prove Sidon enumeration. The Python generator validates transition labels by exact witness lookup, and this Lean certificate proves the resulting finite selected-surface index partitions for n=61,62,63,64: direct previous selected-surface persistence is empty, +1 previous-surface inheritance has counts 2,2,3,6, new selected-surface counts are 2,4,7,5, every transition is mixed, and each joint surface remains less than one twentieth of its exported face. This is transition bookkeeping over packet-backed selected indices, not theorem progress on Erdos #30.

Sidon `n=65..71` selected-surface transition certificate:

- Module: `Erdos30_FaceField_65_71_JointSurface_Transition_Certificate`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_65_71_JointSurface_Transition_Certificate.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_FaceField_65_71_JointSurface_Transition_Certificate`
- Theorems include:
  `Erdos30FaceField6571JointSurfaceTransitionCertificate.joint_surface_transition_count_table`,
  `Erdos30FaceField6571JointSurfaceTransitionCertificate.plus_one_previous_surface_count_table`,
  `Erdos30FaceField6571JointSurfaceTransitionCertificate.new_surface_count_table`,
  `Erdos30FaceField6571JointSurfaceTransitionCertificate.all_surface_transition_counts_balance`,
  `Erdos30FaceField6571JointSurfaceTransitionCertificate.all_surface_transitions_are_mixed`,
  `Erdos30FaceField6571JointSurfaceTransitionCertificate.joint_surface_transition_certificate_passes`

Boundary: This does not import full witness faces or independently prove Sidon enumeration. The Python generator validates complete exports and transition labels by exact witness lookup, and this Lean certificate proves the resulting finite selected-surface index partitions for n=65..71: direct previous selected-surface persistence is empty, +1 previous-surface inheritance has counts 32,24,85,74,179,125,324, new selected-surface counts are 15,16,35,35,88,49,128, every transition is mixed, and each joint surface remains less than one twentieth of its exported face. This is transition bookkeeping over packet-backed selected indices, not theorem progress on Erdos #30.

Sidon `n=58..71` branch summary certificate:

- Module: `Erdos30_GroundFaceBranch_58_71_Certificate`
- Path: `erdos-experiments/Erdos30/lean/Erdos30_GroundFaceBranch_58_71_Certificate.lean`
- Build command: `cd erdos-experiments/Erdos30 && lake build Erdos30_GroundFaceBranch_58_71_Certificate`
- Theorems include:
  `Erdos30GroundFaceBranch5871Certificate.all_exports_are_complete_by_count`,
  `Erdos30GroundFaceBranch5871Certificate.skeleton_count_strictly_increases_58_71`,
  `Erdos30GroundFaceBranch5871Certificate.branch_summary_certificate_passes`

Boundary: This certifies the compact finite branch table for n=58..71 derived from the exact packet: exported counts equal exact face counts, previous-face and +1 previous-face persistence holds across adjacent rows by count, skeleton counts strictly increase, and Pareto/joint selected surfaces are proper count-level subsets. It does not import every witness literal, does not prove the packet enumerator, and is not theorem progress on Erdos #30.

## Exact Row Tables

### Sidon #30

| n | h | ground | h-1 | h-2 | split | h jump | ground ratio |
|---:|---:|---:|---:|---:|---|---|---:|
| 20 | 6 | 206 | 3734 | 3786 | yes | no |  |
| 21 | 6 | 504 | 5428 | 4750 | yes | no | 2.446602 |
| 22 | 6 | 1004 | 7612 | 5874 | yes | no | 1.992063 |
| 23 | 6 | 1910 | 10488 | 7196 | yes | no | 1.902390 |
| 24 | 6 | 3380 | 14126 | 8720 | yes | no | 1.769634 |
| 25 | 7 | 10 | 5688 | 18744 | yes | yes | 0.002959 |
| 26 | 7 | 34 | 9036 | 24390 | yes | no | 3.400000 |
| 27 | 7 | 98 | 14106 | 31436 | yes | no | 2.882353 |
| 28 | 7 | 282 | 21190 | 39914 | yes | no | 2.877551 |
| 29 | 7 | 760 | 31158 | 50212 | yes | no | 2.695035 |
| 30 | 7 | 1618 | 44370 | 62390 | yes | no | 2.128947 |

### B_2[2] #755

| n | h | ground | h-1 | h-2 | split | h jump | ground ratio |
|---:|---:|---:|---:|---:|---|---|---:|
| 20 | 9 | 3998 | 45026 | 0 | yes | no |  |
| 21 | 10 | 28 | 12168 | 0 | yes | yes | 0.007004 |
| 22 | 10 | 182 | 30188 | 0 | yes | no | 6.500000 |
| 23 | 10 | 1156 | 70432 | 0 | yes | no | 6.351648 |
| 24 | 10 | 4596 | 146516 | 0 | yes | no | 3.975779 |
| 25 | 10 | 16350 | 291280 | 0 | yes | no | 3.557441 |
| 26 | 11 | 20 | 45540 | 0 | yes | yes | 0.001223 |
| 27 | 11 | 296 | 120008 | 0 | yes | no | 14.800000 |
| 28 | 11 | 1694 | 273734 | 0 | yes | no | 5.722973 |
| 29 | 11 | 7908 | 598156 | 0 | yes | no | 4.668241 |
| 30 | 12 | 6 | 27856 | 0 | yes | yes | 0.000759 |

### Sum-Free #166

| n | h | ground | h-1 | h-2 | split | h jump | ground ratio |
|---:|---:|---:|---:|---:|---|---|---:|
| 20 | 10 | 3 | 32 | 0 | no | no |  |
| 21 | 11 | 2 | 24 | 0 | no | yes | 0.666667 |
| 22 | 11 | 3 | 35 | 0 | no | no | 1.500000 |
| 23 | 12 | 2 | 26 | 0 | no | yes | 0.666667 |
| 24 | 12 | 3 | 38 | 0 | no | no | 1.500000 |
| 25 | 13 | 2 | 28 | 0 | no | yes | 0.666667 |
| 26 | 13 | 3 | 41 | 0 | no | no | 1.500000 |
| 27 | 14 | 2 | 30 | 0 | no | yes | 0.666667 |
| 28 | 14 | 3 | 44 | 0 | no | no | 1.500000 |
| 29 | 15 | 2 | 32 | 0 | no | yes | 0.666667 |
| 30 | 15 | 3 | 47 | 0 | no | no | 1.500000 |

### B_2[3] #755 Ground-Only Reset Window

| n | h | ground | h-1 | h-2 | split | h jump | ground ratio |
|---:|---:|---:|---:|---:|---|---|---:|
| 43 | 17 | 32002 | NA | NA | yes | no |  |
| 44 | 17 | 212586 | NA | NA | yes | no | 6.642897 |
| 45 | 18 | 8 | NA | NA | yes | yes | 0.000038 |
| 46 | 18 | 142 | NA | NA | no | no | 17.750000 |
| 47 | 18 | 2160 | NA | NA | yes | no | 15.211268 |

## Proof Targets

### PT-1. Exposed-face split lemma

**Status:** COMPILED_ABSTRACT_N30_WINDOW_BRANCH_FULLFACE_N60_N61_71_SURFACE_TRANSITION_AND_SUMMARY_CERTIFICATES

**Candidate statement:** For a finite extremal family F and two observables phi, psi, if their zero-temperature argmin witnesses differ on F, then the mathematical object under study is the extremal face with field-exposed points, not a unique optimizer selected by enumeration order.

**Finite evidence:** Sidon split 11/11 and B_2[2] split 11/11 under the same observable family; sum-free split 0/11. The Sidon n=20..30 prefix/mass split, n=57/n=58 branch split, complete exported n=56..58 handoff face with exact integer mass winners, exact zero prefix probes, strict positive exact-prefix non-winner gaps, exact joint-key winners from exact prefix plus exact mass, full-prefix segment bounds, full-prefix residual-zero winners, actual residual-positive probe non-winners, scalar full-prefix residual maxima equal to packet probes, face-local scalar/probe equality, scalar full-prefix minimizers, scalar full-prefix/exact-mass splits for n=57 and n=58, scalar full-prefix joint keys equal to exact-probe joint keys, scalar full-prefix joint winners W3/W5/W7 via reusable minimizer transfer, the n=59 scalar full-prefix joint microcertificate with 18 exported witnesses and scalar joint winner W7, the n=60 scalar full-prefix joint microcertificate with 54 exported witnesses and scalar joint winner W43, the n=61 single-row compact joint minimizer-surface certificate, the n=61..64 compact joint minimizer-surface certificate with surface sizes 4,6,10,11 over exported faces of size 152,398,1022,2360, the n=61..64 selected-surface transition certificate, the n=65..71 selected-surface transition certificate with mixed +1-inherited/new partitions and no direct previous selected-surface persistence across faces as large as 203,840, and the compact n=58..71 branch table now have compiled finite Lean certificates.

**Next formal move:** Turn the compiled n=61..71 transition bookkeeping into a symbolic mechanism: explain why selected joint surfaces inherit only through +1-shifted previous witnesses while also introducing new selected witnesses, and why the selected surface stays tiny relative to the exported face. Keep further extensions index/scalar-based; do not scale full witness-literal imports beyond their useful size.

**Rejection gate:** Do not present this as theorem progress on #30; it is proof-target infrastructure.

### PT-2. Sidon/B_2[g] field-response conjecture

**Status:** CONJECTURE_TARGET_WITH_NEGATIVE_CONTROL

**Candidate statement:** Sidon-like bounded-additive-representation faces generically have multiple field-exposed extremal witnesses under prefix and mass observables, while simple sum-free faces can remain field-rigid.

**Finite evidence:** Sidon 20..30: 11/11 split; B_2[2] 20..30: 11/11 split; B_2[3] 43..47 ground-only: 4/5 split; sum-free 20..30: 0/11 split.

**Next formal move:** Find a structural condition that predicts split versus quiet rows; B_2[3] n=46 is the first required exception case.

**Rejection gate:** If larger Sidon/B_2[g] windows lose split behavior broadly, or sum-free begins splitting under the same fields, demote the conjecture.

### PT-3. Entropy reset at cardinality jumps

**Status:** PATTERN_TARGET

**Candidate statement:** When the exact finite maximum cardinality h(n) increases, the new extremal face can reset to very low degeneracy before expanding across the following plateau.

**Finite evidence:** Sidon h jump at n=25 resets 3380 -> 10 ground states; B_2[2] jumps at n=21,26,30 reset 3998 -> 28, 16350 -> 20, 7908 -> 6; B_2[3] jump at n=45 resets 212586 -> 8.

**Next formal move:** Translate the reset into construction scarcity: characterize why first rows at a new h have few admissible extremal extensions.

**Rejection gate:** Treat as finite pattern until a wider window or symbolic construction explains the reset.

### PT-4. Quiet-control characterization

**Status:** CONTROL_TARGET

**Candidate statement:** For the sum-free encoding on [0,n], the same observable family can share a common selected ground-state witness across the checked finite window.

**Finite evidence:** Sum-free 20..30 has 28 total ground states and 0/11 field splits; formula parity matched in every row.

**Next formal move:** Characterize the exact maximum sum-free families and prove when the prefix/mass/joint fields co-select the same witness.

**Rejection gate:** If the control is an artifact of the chosen observables, keep it as an audit control rather than a theorem lane.

### PT-5. Ground versus near-ground boundary

**Status:** AUDIT_TARGET

**Candidate statement:** Claims about ground-state face structure must be separated from h-1/h-2 near-ground structure and from capped lower-bound runs.

**Finite evidence:** Sidon and B_2[2] Atheneum rows have exact h-1 counts; B_2[3] 43..47 is ground-only and cannot support near-ground claims.

**Next formal move:** Keep separate lemmas for exact ground faces and near-ground layers; never mix exact counts with capped lower bounds.

**Rejection gate:** Any proof target that depends on B_2[3] near-ground counts is invalid until a new exact near-ground packet exists.


## Recommended Next Move

PT-1 has now moved from selected points to selected surfaces through n=64, and
then to compiled transition partitions through n=71. The next
useful proof target is no longer another table for its own sake; it is a
symbolic mechanism explaining why only `+1`-shifted previous witnesses persist
into the selected surface, why new selected witnesses keep appearing, and why
the joint-selected surface stays a tiny proper subface while the exported
extremal face expands. Keep PT-4 as the control. Keep PT-5 as the audit rule.

The paper's strategic sentence remains the guardrail: the analogy is the
microscope, not the proof.
