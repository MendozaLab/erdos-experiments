/**
 * Cap Set Transfer Matrix — PMF Salvo for Erdős #233
 *
 * Prima Materia Framework, Layer 2: Arithmetic Phase (Z_3^n lattice)
 *
 * A cap set is a subset of Z_3^n containing no three-term arithmetic
 * progression (3-AP): no x, y, z with x + y + z ≡ 0 (mod 3) component-wise
 * (equivalently, no x, y, z with y = (x+z)/2 mod 3, i.e., x+y+z = 0 mod 3).
 *
 * Erdős #233: What is the maximum size of a cap set in Z_3^n?
 *
 * Known results:
 *   - Trivial: r_3(Z_3^n) ≤ 3^n (no constraint)
 *   - Meshulam (1995): r_3(Z_3^n) ≤ 3^n / n
 *   - Bateman-Katz (2012): r_3(Z_3^n) ≤ 3^n / n^{1+ε}
 *   - Croot-Lev-Pach / Ellenberg-Gijswijt (2016): r_3(Z_3^n) ≤ 2.756^n
 *     This was a BREAKTHROUGH using the polynomial method (slice rank).
 *   - Exact small values: r_3(Z_3^1)=2, r_3(Z_3^2)=4, r_3(Z_3^3)=9,
 *     r_3(Z_3^4)=20, r_3(Z_3^5)=45, r_3(Z_3^6)=112
 *
 * PMF formulation:
 *   - Lattice: Z_3^n (3^n points)
 *   - Constraint: no x+y+z ≡ 0 mod 3 (no 3-AP)
 *   - Transfer matrix: process Z_3^n layer by layer. The n-dimensional
 *     cube Z_3^n can be decomposed as Z_3^{n-1} × Z_3. The transfer
 *     matrix T acts on "slices" — one Z_3 coordinate at a time.
 *
 * Layer-by-layer transfer matrix:
 *   State at layer i: which elements of Z_3^{i} × {0,1,2} are in the cap.
 *   For a "column" approach: fix the last coordinate, process columns.
 *
 * Actually, the clean approach: enumerate all cap-free subsets of Z_3^W
 * for small W (the 1D analogue), then extend to higher dimensions.
 *
 * For Z_3^1 = {0,1,2}: cap sets avoid x+y+z=0 mod 3.
 *   3-APs in Z_3: (0,1,2) since 0+1+2=3≡0. So max cap = 2 (e.g., {0,1}).
 *
 * For the transfer matrix, we use a "slice" decomposition:
 *   Z_3^n = Z_3^{n-1} × {0, 1, 2}
 *   Each "slice" is a copy of Z_3^{n-1}.
 *   The state encodes which elements in a slice are in the cap set.
 *   Transition between slices: adding elements in slice i must not
 *   create 3-APs with elements in slices i-1 and i-2.
 *
 * Compile: g++ -O3 -march=native -o capset_tm capset_transfer_matrix.cpp
 * Run:     ./capset_tm [max_dim]
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <chrono>

using namespace std;
using u64 = uint64_t;

// ─── Z_3 arithmetic ─────────────────────────────────────────────────────

// Represent a vector in Z_3^d as a packed integer: each "trit" uses 2 bits.
// For d ≤ 16, this fits in a u64 (32 bits needed for 16 trits).

// Get trit i from packed vector
int get_trit(u64 v, int i) {
    return (v >> (2 * i)) & 3;  // stored as 0, 1, 2
}

// Set trit i in packed vector
u64 set_trit(u64 v, int i, int val) {
    u64 mask = ~(3ULL << (2 * i));
    return (v & mask) | ((u64)val << (2 * i));
}

// Add two vectors in Z_3^d
u64 add_z3(u64 a, u64 b, int d) {
    u64 result = 0;
    for (int i = 0; i < d; i++) {
        int sum = (get_trit(a, i) + get_trit(b, i)) % 3;
        result = set_trit(result, i, sum);
    }
    return result;
}

// Negate a vector in Z_3^d (for checking x+y+z=0 ↔ z = -(x+y))
u64 neg_z3(u64 a, int d) {
    u64 result = 0;
    for (int i = 0; i < d; i++) {
        int neg = (3 - get_trit(a, i)) % 3;
        result = set_trit(result, i, neg);
    }
    return result;
}

// Enumerate all elements of Z_3^d
vector<u64> enumerate_z3(int d) {
    int total = 1;
    for (int i = 0; i < d; i++) total *= 3;
    vector<u64> elems(total);
    for (int idx = 0; idx < total; idx++) {
        u64 v = 0;
        int tmp = idx;
        for (int i = 0; i < d; i++) {
            v = set_trit(v, i, tmp % 3);
            tmp /= 3;
        }
        elems[idx] = v;
    }
    return elems;
}

// ─── Cap Set Enumeration ─────────────────────────────────────────────────
// Enumerate all cap sets (3-AP-free subsets) of Z_3^d by backtracking.

struct CapSetEnumerator {
    int d;
    int total;  // 3^d
    vector<u64> elements;
    vector<long long> density_of_states;
    long long total_count;
    int max_cap_size;

    // For incremental checking: maintain set of "forbidden" elements
    // An element z is forbidden if there exist x, y in the cap with x+y+z ≡ 0
    unordered_set<u64> forbidden;
    vector<u64> current_cap;

    void enumerate(int dim) {
        d = dim;
        total = 1;
        for (int i = 0; i < d; i++) total *= 3;
        elements = enumerate_z3(d);

        density_of_states.assign(total + 1, 0);
        total_count = 0;
        max_cap_size = 0;
        forbidden.clear();
        current_cap.clear();

        backtrack(0);
    }

    void backtrack(int start) {
        int sz = (int)current_cap.size();
        density_of_states[sz]++;
        total_count++;
        if (sz > max_cap_size) max_cap_size = sz;

        for (int i = start; i < total; i++) {
            u64 z = elements[i];

            // Check: z is not forbidden
            if (forbidden.count(z)) continue;

            // Check: z+z+z = 3z ≡ 0, so (z,z,z) is always a 3-AP.
            // But we don't allow repeated elements, so this is fine.
            // Actually, a 3-AP requires x+y+z ≡ 0 with x,y,z all in the set.
            // If z is not forbidden, no existing pair creates z as the third.
            // But we also need: no pair (z, existing) creates a NEW forbidden third.
            // That's handled by updating forbidden after adding z.

            // Add z to cap
            current_cap.push_back(z);

            // Update forbidden: for each existing y in cap (before z), mark -(z+y)
            // Also: for z itself, mark -(z+z) = -(2z) mod 3
            vector<u64> newly_forbidden;
            for (int j = 0; j < sz; j++) {  // sz = cap size before adding z
                u64 y = current_cap[j];
                u64 f = neg_z3(add_z3(z, y, d), d);
                if (!forbidden.count(f)) {
                    forbidden.insert(f);
                    newly_forbidden.push_back(f);
                }
            }
            // Also forbid -(2z) mod 3
            u64 f2 = neg_z3(add_z3(z, z, d), d);
            bool added_f2 = false;
            if (!forbidden.count(f2)) {
                forbidden.insert(f2);
                added_f2 = true;
            }

            backtrack(i + 1);

            // Remove z from cap and undo forbidden
            current_cap.pop_back();
            for (u64 f : newly_forbidden) forbidden.erase(f);
            if (added_f2) forbidden.erase(f2);
        }
    }
};

// ─── Slice Transfer Matrix for Z_3^n ────────────────────────────────────
// Decompose Z_3^n as Z_3^{n-1} × Z_3.
// Process the 3 "slices" (last coordinate = 0, 1, 2) sequentially.
//
// State: which elements of Z_3^{n-1} are in the cap in each of the
//        previous slices. For a window of 2 slices, the state is a pair
//        (S_0, S_1) where S_i ⊂ Z_3^{n-1}.
//
// Transition from slice (a,b) → (b,c):
//   S_c is the cap elements in slice with last coord = c.
//   Constraint: for each x ∈ S_c and y ∈ S_a (two slices back),
//     -(x+y) mod 3 in Z_3^{n-1} must NOT be in S_b.
//   (Because (x,a), (-(x+y),b), (y,c) would be a 3-AP if a+b+c ≡ 0 mod 3.)
//   Actually, the constraint depends on the specific slice indices.
//
// For 3 slices (coords 0,1,2): 0+1+2 = 3 ≡ 0 mod 3, so the main
// constraint is between slices 0, 1, 2.

struct SliceTransferMatrix {
    int d;  // dimension (so Z_3^d, sliced into 3 copies of Z_3^{d-1})
    int slice_size;  // 3^{d-1}
    vector<u64> slice_elements;  // elements of Z_3^{d-1}

    // States: cap subsets of a single slice Z_3^{d-1}
    // Since we track 2 slices, state = pair of cap subsets
    // This is too large for d > 3. For d=2, slice = Z_3^1 = {0,1,2}, 8 subsets.
    // For d=3, slice = Z_3^2 = 9 elements, 2^9 = 512 subsets.

    // For this PMF salvo, we use the full enumeration approach for small d
    // and extract the "effective λ_max" from the growth of max cap size.

    void compute(int dim) {
        d = dim;
        slice_size = 1;
        for (int i = 0; i < d - 1; i++) slice_size *= 3;
        slice_elements = enumerate_z3(d - 1);

        printf("  Slice transfer matrix for Z_3^%d\n", d);
        printf("  Slice size: 3^%d = %d\n", d - 1, slice_size);

        // Enumerate cap subsets of a single slice
        // A cap subset of Z_3^{d-1} is a 3-AP-free subset
        vector<u64> slice_caps;  // bitmask representation
        vector<u64> current;
        unordered_set<u64> forbidden;
        enumerate_slice_caps(0, current, forbidden, slice_caps);

        printf("  Cap subsets of Z_3^%d: %d\n", d - 1, (int)slice_caps.size());

        // Build transfer matrix between consecutive slices
        // For slices at positions a, b, c with a+b+c ≡ 0 mod 3:
        //   For each x in S_a and z in S_c, the "middle" element
        //   y = -(x+z) must NOT be in S_b.
        //
        // The transfer matrix T(S_prev, S_curr) = 1 iff S_curr is a cap
        // subset AND no 3-AP is created between S_prev (2 back) and S_curr
        // through the middle slice.
        //
        // For simplicity: enumerate all pairs (S_0, S_1, S_2) of cap subsets
        // of Z_3^{d-1} and check the cross-slice 3-AP condition.

        if (slice_caps.size() <= 10000) {
            // Build transition matrix: slice 0 → slice 1 (no cross constraint)
            // Then slice (0,1) → slice 2: constraint is 3-APs across slices

            int n_caps = (int)slice_caps.size();
            printf("  Building cross-slice constraint matrix (%d × %d)...\n",
                   n_caps, n_caps);

            // Convert cap bitmasks to element lists for fast lookup
            vector<vector<int>> cap_elements(n_caps);
            for (int i = 0; i < n_caps; i++) {
                for (int j = 0; j < slice_size; j++) {
                    if (slice_caps[i] & (1ULL << j)) {
                        cap_elements[i].push_back(j);
                    }
                }
            }

            // For the 3-slice constraint (slices at Z_3 positions 0, 1, 2):
            // A 3-AP (x,0), (y,1), (z,2) requires x+y+z ≡ 0 in Z_3^{d-1}
            // (since 0+1+2 = 0 mod 3).
            // So: for x ∈ S_0, z ∈ S_2, we need -(x+z) ∉ S_1.

            // Count valid triples (S_0, S_1, S_2)
            long long valid_triples = 0;
            long long total_triples = 0;
            int max_total_size = 0;

            // For reporting: track max |S_0| + |S_1| + |S_2|
            // Use sampling for large n_caps
            int sample_limit = min(n_caps, 200);

            for (int i0 = 0; i0 < sample_limit; i0++) {
                for (int i1 = 0; i1 < sample_limit; i1++) {
                    for (int i2 = 0; i2 < sample_limit; i2++) {
                        total_triples++;

                        // Check cross-slice 3-AP
                        bool valid = true;
                        for (int x : cap_elements[i0]) {
                            if (!valid) break;
                            for (int z : cap_elements[i2]) {
                                // Compute -(x_vec + z_vec) in Z_3^{d-1}
                                u64 x_vec = slice_elements[x];
                                u64 z_vec = slice_elements[z];
                                u64 neg_sum = neg_z3(add_z3(x_vec, z_vec, d - 1), d - 1);
                                // Check if neg_sum is in S_1
                                // Find index of neg_sum in slice_elements
                                for (int y : cap_elements[i1]) {
                                    if (slice_elements[y] == neg_sum) {
                                        valid = false;
                                        break;
                                    }
                                }
                                if (!valid) break;
                            }
                        }

                        if (valid) {
                            valid_triples++;
                            int total_size = (int)(cap_elements[i0].size() +
                                                   cap_elements[i1].size() +
                                                   cap_elements[i2].size());
                            if (total_size > max_total_size) max_total_size = total_size;
                        }
                    }
                }
            }

            printf("  Valid triples (sampled %d^3): %lld / %lld = %.6f\n",
                   sample_limit, valid_triples, total_triples,
                   (double)valid_triples / total_triples);
            printf("  Max |S_0|+|S_1|+|S_2| among valid: %d\n", max_total_size);
            printf("  Known r_3(Z_3^%d) = ", d);

            // Print known values
            int known[] = {0, 2, 4, 9, 20, 45, 112};
            if (d <= 6) printf("%d", known[d]);
            else printf("?");
            printf("\n");
        }
    }

    void enumerate_slice_caps(int start, vector<u64>& current,
                               unordered_set<u64>& forbidden,
                               vector<u64>& result) {
        // Record current subset as bitmask
        u64 mask = 0;
        for (u64 idx : current) mask |= (1ULL << idx);
        result.push_back(mask);

        for (int i = start; i < slice_size; i++) {
            u64 elem = slice_elements[i];
            if (forbidden.count(elem)) continue;

            current.push_back(i);

            // Update forbidden
            vector<u64> newly_forbidden;
            for (int j = 0; j < (int)current.size() - 1; j++) {
                u64 y = slice_elements[current[j]];
                u64 f = neg_z3(add_z3(elem, y, d - 1), d - 1);
                if (!forbidden.count(f)) {
                    forbidden.insert(f);
                    newly_forbidden.push_back(f);
                }
            }
            u64 f2 = neg_z3(add_z3(elem, elem, d - 1), d - 1);
            bool added_f2 = false;
            if (!forbidden.count(f2)) {
                forbidden.insert(f2);
                added_f2 = true;
            }

            enumerate_slice_caps(i + 1, current, forbidden, result);

            current.pop_back();
            for (u64 f : newly_forbidden) forbidden.erase(f);
            if (added_f2) forbidden.erase(f2);
        }
    }
};

// ─── Main ────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    int max_d = 5;
    if (argc > 1) max_d = atoi(argv[1]);

    printf("Cap Set Transfer Matrix — PMF Salvo for Erdős #233\n");
    printf("====================================================\n\n");
    printf("Cap sets in Z_3^n: max 3-AP-free subset\n");
    printf("CLP/EG bound (2016): r_3(Z_3^n) ≤ 2.756^n\n\n");

    FILE* json_out = fopen("EXP-MATH-ERDOS233-CAPSET-001_RESULTS.json", "w");
    fprintf(json_out, "{\n  \"experiment_id\": \"EXP-MATH-ERDOS233-CAPSET-001\",\n");
    fprintf(json_out, "  \"title\": \"Cap Set Lattice Gas — PMF Layer 2\",\n");
    fprintf(json_out, "  \"problem\": \"Erdos #233: Cap sets in Z_3^n\",\n");
    fprintf(json_out, "  \"known_bound\": \"r_3(Z_3^n) <= 2.756^n (CLP/EG 2016)\",\n");
    fprintf(json_out, "  \"enumeration\": [\n");

    bool first_json = true;
    int known_values[] = {0, 2, 4, 9, 20, 45, 112, 236};

    // Part 1: Full enumeration for small d
    printf("PART 1: Full Enumeration of Cap Sets in Z_3^d\n");
    printf("───────────────────────────────────────────────\n\n");

    for (int d = 1; d <= min(max_d, 6); d++) {
        int total = 1;
        for (int i = 0; i < d; i++) total *= 3;

        if (total > 200000) {
            printf("Z_3^%d has %d elements — skipping full enumeration\n\n", d, total);
            continue;
        }

        printf("═══ Z_3^%d (%d elements) ═══\n", d, total);
        auto t0 = chrono::high_resolution_clock::now();

        CapSetEnumerator cse;
        cse.enumerate(d);

        auto t1 = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration<double>(t1 - t0).count();

        double density = (double)cse.max_cap_size / total;
        double log2_count = (cse.total_count > 1) ? log2((double)cse.total_count) : 0;
        double effective_base = (d > 0) ? pow(cse.max_cap_size, 1.0 / d) : 0;

        printf("  Total cap-free families: %lld\n", cse.total_count);
        printf("  Max cap size: %d (known: %d)\n",
               cse.max_cap_size, (d <= 7) ? known_values[d] : -1);
        printf("  Density: %d/%d = %.6f\n", cse.max_cap_size, total, density);
        printf("  Effective base: %.6f (CLP bound: 2.756)\n", effective_base);
        printf("  log2(#cap-free subsets): %.4f\n", log2_count);
        printf("  Time: %.3f s\n", elapsed);

        // Density of states
        printf("  Density of states: ");
        for (int k = 0; k <= min(cse.max_cap_size, 20); k++) {
            printf("c(%d)=%lld ", k, cse.density_of_states[k]);
        }
        if (cse.max_cap_size > 20) printf("...");
        printf("\n\n");

        // JSON
        if (!first_json) fprintf(json_out, ",\n");
        first_json = false;
        fprintf(json_out, "    {\"d\": %d, \"total_elements\": %d, "
                "\"total_capfree\": %lld, \"max_cap_size\": %d, "
                "\"known_value\": %d, \"density\": %.8f, "
                "\"effective_base\": %.8f, \"log2_count\": %.6f, "
                "\"time_s\": %.4f}",
                d, total, cse.total_count, cse.max_cap_size,
                (d <= 7) ? known_values[d] : -1, density, effective_base,
                log2_count, elapsed);

        if (elapsed > 120) {
            printf("  Time limit — stopping enumeration\n");
            break;
        }
    }
    fprintf(json_out, "\n  ],\n");

    // Part 2: Slice transfer matrix
    printf("\nPART 2: Slice Transfer Matrix\n");
    printf("─────────────────────────────\n\n");

    fprintf(json_out, "  \"slice_transfer_matrix\": [\n");
    first_json = true;

    for (int d = 2; d <= min(max_d, 4); d++) {
        printf("═══ Slice TM for Z_3^%d ═══\n", d);
        auto t0 = chrono::high_resolution_clock::now();

        SliceTransferMatrix stm;
        stm.compute(d);

        auto t1 = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration<double>(t1 - t0).count();
        printf("  Time: %.3f s\n\n", elapsed);

        if (!first_json) fprintf(json_out, ",\n");
        first_json = false;
        fprintf(json_out, "    {\"d\": %d, \"time_s\": %.4f}", d, elapsed);

        if (elapsed > 120) break;
    }

    fprintf(json_out, "\n  ]\n}\n");
    fclose(json_out);

    // Part 3: Scaling analysis
    printf("\nPART 3: Scaling Analysis\n");
    printf("────────────────────────\n\n");
    printf("  d | r_3(Z_3^d) | r/3^d    | r^{1/d}  | CLP bound\n");
    printf("  --+------------+----------+----------+-----------\n");
    for (int d = 1; d <= 7; d++) {
        int total = 1;
        for (int i = 0; i < d; i++) total *= 3;
        double density = (double)known_values[d] / total;
        double eff_base = pow(known_values[d], 1.0 / d);
        double clp = pow(2.756, d);
        printf("  %d | %10d | %8.6f | %8.5f | %9.1f\n",
               d, known_values[d], density, eff_base, clp);
    }
    printf("\n  CLP/EG effective base: 2.756\n");
    printf("  If PMF λ_max → 2.756, it reproduces the polynomial method bound.\n");

    printf("\nResults saved to EXP-MATH-ERDOS233-CAPSET-001_RESULTS.json\n");
    return 0;
}
