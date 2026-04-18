/**
 * Sum-Free Transfer Matrix — PMF Salvo for Erdős #166
 *
 * Prima Materia Framework, Layer 2: Arithmetic Phase
 * Constraint: A ⊂ {1,...,N} is sum-free iff (A+A) ∩ A = ∅
 *
 * Two-pronged computation:
 *   1. EXACT partition function Z(N) by full enumeration of sum-free subsets
 *      of {1,...,N} for N up to ~35. Tracks density of states c(N,k) =
 *      #sum-free subsets of size k. Phase transition at ρ = 1/2.
 *   2. LOCAL transfer matrix: enumerate "locally sum-free" subsets of a window
 *      {1,...,W}, build shift-and-constrain transfer matrix, compute λ_max.
 *      The local constraint is: no a+b=c with a,b,c all in the window AND
 *      all in the subset. This is an approximation that tightens as W grows.
 *
 * Known result: max |A|/N → 1/2 (take all odd numbers, or all n > N/2).
 * Cameron-Erdős (proved by Green 2004, Sapozhenko 2003):
 *   #{sum-free A ⊂ {1,...,N}} = Θ(2^{N/2}).
 *
 * Calibration target: if PMF reproduces ρ_c = 1/2 and the 2^{N/2} count
 * scaling, the machinery is validated before attacking open problems.
 *
 * Compile: g++ -O3 -march=native -o sumfree_tm sumfree_transfer_matrix.cpp
 * Run:     ./sumfree_tm [max_N for enumeration] [max_W for transfer matrix]
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <chrono>

using namespace std;
using u64 = uint64_t;

// ─── Part 1: Exact Partition Function Enumeration ────────────────────────
// Enumerate all sum-free subsets of {1,...,N} by backtracking.
// Track c(N,k) = count of sum-free subsets of size exactly k.

struct SumFreeEnumerator {
    int N;
    vector<long long> density_of_states;  // c[k] = # subsets of size k
    long long total_count;
    int max_size;

    // sumset_hits[v] = # pairs (a,b) in current set with a+b=v
    // We use a simpler approach: track which values are "blocked"
    // blocked[v] = true if v = a+b for some a,b in current set
    vector<int> blocked;  // blocked[v] = # pairs in current set summing to v
    vector<bool> in_set;

    void enumerate(int max_N) {
        N = max_N;
        density_of_states.assign(N + 1, 0);
        total_count = 0;
        max_size = 0;
        blocked.assign(2 * N + 1, 0);
        in_set.assign(N + 1, false);

        vector<int> current;
        backtrack(1, current);
    }

    void backtrack(int start, vector<int>& current) {
        // Record current subset
        int k = (int)current.size();
        density_of_states[k]++;
        total_count++;
        if (k > max_size) max_size = k;

        for (int x = start; x <= N; x++) {
            // Check: x must not be blocked (i.e., x ≠ a+b for any a,b in set)
            if (blocked[x] > 0) continue;

            // Check: for all a in current set, x+a must not be in current set
            bool valid = true;
            for (int a : current) {
                if (a + x <= N && in_set[a + x]) {
                    valid = false;
                    break;
                }
            }
            // Also check 2x: if 2x is in the set, that's x+x=2x violation
            // (already handled by blocked[2x] check when 2x was added)
            // But we also need: for all a in set, a+x not in set (checked above)
            // And: x+x = 2x not in set (handled by checking in_set[2x] if 2x<=N)
            if (valid && 2 * x <= N && in_set[2 * x]) {
                valid = false;
            }

            if (valid) {
                // Add x to the set
                in_set[x] = true;
                // Update blocked: for all a in current set, mark x+a as blocked
                for (int a : current) {
                    blocked[a + x]++;
                }
                // Also mark 2x as blocked (x+x)
                blocked[2 * x]++;

                current.push_back(x);
                backtrack(x + 1, current);
                current.pop_back();

                // Remove x from the set
                in_set[x] = false;
                for (int a : current) {
                    blocked[a + x]--;
                }
                blocked[2 * x]--;
            }
        }
    }
};

// ─── Part 2: Local Sum-Free Transfer Matrix ──────────────────────────────
// States: sum-free subsets of {1,...,W-1} (position 0 excluded since 0+a=a)
// Actually we work with subsets of {0,...,W-1} but label them as values.
// A subset S is "locally sum-free" if no a,b,c ∈ S with a+b=c and a,b,c ∈ {1,...,W-1}.
//
// Transfer matrix construction:
// For each valid state s (encoded as bitmask of positions 1..W-1):
//   shifted = s >> 1 (shift all positions down by 1)
//   if shifted is a valid state:
//     T[s, shifted] = 1 (don't add new element)
//     if shifted | (1 << (W-1)) is also valid:
//       T[s, shifted | (1<<(W-1))] = 1 (add new element at top)

struct LocalSumFreeStates {
    int W;
    vector<u64> states;
    unordered_map<u64, int> state_index;

    bool is_sumfree(u64 mask) const {
        // Extract positions
        vector<int> pos;
        for (int i = 1; i < W; i++) {
            if (mask & (1ULL << i)) pos.push_back(i);
        }
        // Check: no a+b=c with a,b,c all in pos
        // Use a set for quick lookup
        for (int i = 0; i < (int)pos.size(); i++) {
            for (int j = i; j < (int)pos.size(); j++) {
                int s = pos[i] + pos[j];
                if (s < W && (mask & (1ULL << s))) {
                    return false;
                }
            }
        }
        return true;
    }

    void enumerate(int max_W) {
        W = max_W;
        states.clear();
        state_index.clear();

        // Enumerate all sum-free subsets of {1,...,W-1} by backtracking
        // Include empty set
        states.push_back(0ULL);
        state_index[0ULL] = 0;

        vector<int> current;
        vector<int> blocked(2 * W, 0);
        backtrack(1, current, blocked);
    }

    void backtrack(int start, vector<int>& current, vector<int>& blocked) {
        if (!current.empty()) {
            u64 mask = 0;
            for (int p : current) mask |= (1ULL << p);
            int idx = (int)states.size();
            states.push_back(mask);
            state_index[mask] = idx;
        }

        for (int x = start; x < W; x++) {
            // Check if x is blocked (some a+b=x with a,b in current)
            if (blocked[x] > 0) continue;

            // Check if x+a is in current for any a
            bool valid = true;
            u64 cur_mask = 0;
            for (int a : current) cur_mask |= (1ULL << a);

            for (int a : current) {
                int s = x + a;
                if (s < W && (cur_mask & (1ULL << s))) {
                    valid = false;
                    break;
                }
            }
            if (valid && 2 * x < W && (cur_mask & (1ULL << (2 * x)))) {
                valid = false;
            }

            if (valid) {
                // Mark new blocked values
                for (int a : current) {
                    blocked[a + x]++;
                }
                blocked[2 * x]++;

                current.push_back(x);
                backtrack(x + 1, current, blocked);
                current.pop_back();

                for (int a : current) {
                    blocked[a + x]--;
                }
                blocked[2 * x]--;
            }
        }
    }
};

struct SparseMatrix {
    int n;
    vector<int> row_ptr;
    vector<int> col_idx;

    void build(int W, const vector<u64>& states, const unordered_map<u64, int>& state_index,
               const LocalSumFreeStates& lsf) {
        n = (int)states.size();
        row_ptr.resize(n + 1, 0);
        col_idx.clear();

        for (int i = 0; i < n; i++) {
            u64 state_i = states[i];
            row_ptr[i] = (int)col_idx.size();

            // Shift: position p becomes p-1, position 0 drops (but position 0
            // is never set in sum-free states since we use {1,...,W-1})
            // Actually, position 1 becomes position 0 after shift, which could
            // create issues. Let's handle carefully.
            u64 shifted = state_i >> 1;

            // Check if shifted state is valid (sum-free in {0,...,W-2})
            // Position 0 after shift was position 1 before, which is ≥ 1.
            // But the shifted state might not be sum-free because values changed.
            // We must check explicitly.
            if (state_index.count(shifted) || lsf.is_sumfree(shifted)) {
                // Option 1: don't occupy W-1
                auto it = state_index.find(shifted);
                if (it != state_index.end()) {
                    col_idx.push_back(it->second);
                }

                // Option 2: occupy W-1
                u64 new_state = shifted | (1ULL << (W - 1));
                if (state_index.count(new_state)) {
                    col_idx.push_back(state_index.at(new_state));
                }
            }
            // If shifted is not valid, this state has no outgoing transitions
            // (it can't propagate in the lattice gas)
        }
        row_ptr[n] = (int)col_idx.size();
    }

    void matvec(const vector<double>& x, vector<double>& y) const {
        fill(y.begin(), y.end(), 0.0);
        for (int i = 0; i < n; i++) {
            double xi = x[i];
            for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                y[col_idx[k]] += xi;
            }
        }
    }
};

double power_iteration(const SparseMatrix& A, vector<double>& v, int max_iter = 1000, double tol = 1e-14) {
    int n = A.n;
    vector<double> w(n);
    srand(42);
    for (int i = 0; i < n; i++) v[i] = (double)rand() / RAND_MAX;
    double norm = 0;
    for (int i = 0; i < n; i++) norm += v[i] * v[i];
    norm = sqrt(norm);
    for (int i = 0; i < n; i++) v[i] /= norm;

    double lambda_old = 0;
    for (int iter = 0; iter < max_iter; iter++) {
        A.matvec(v, w);
        double new_norm = 0;
        for (int i = 0; i < n; i++) new_norm += w[i] * w[i];
        new_norm = sqrt(new_norm);
        if (new_norm < 1e-30) break;
        for (int i = 0; i < n; i++) v[i] = w[i] / new_norm;
        if (fabs(new_norm - lambda_old) < tol * fabs(new_norm)) return new_norm;
        lambda_old = new_norm;
    }
    return lambda_old;
}

double second_eigenvalue(const SparseMatrix& A, const vector<double>& v1, double lambda1, int max_iter = 1000, double tol = 1e-12) {
    int n = A.n;
    vector<double> v(n), w(n);
    srand(123);
    for (int i = 0; i < n; i++) v[i] = (double)rand() / RAND_MAX;
    double dot = 0;
    for (int i = 0; i < n; i++) dot += v[i] * v1[i];
    for (int i = 0; i < n; i++) v[i] -= dot * v1[i];
    double norm = 0;
    for (int i = 0; i < n; i++) norm += v[i] * v[i];
    norm = sqrt(norm);
    if (norm < 1e-30) return 0;
    for (int i = 0; i < n; i++) v[i] /= norm;

    double lambda_old = 0;
    for (int iter = 0; iter < max_iter; iter++) {
        A.matvec(v, w);
        dot = 0;
        for (int i = 0; i < n; i++) dot += w[i] * v1[i];
        for (int i = 0; i < n; i++) w[i] -= dot * v1[i];
        double new_norm = 0;
        for (int i = 0; i < n; i++) new_norm += w[i] * w[i];
        new_norm = sqrt(new_norm);
        if (new_norm < 1e-30) return 0;
        for (int i = 0; i < n; i++) v[i] = w[i] / new_norm;
        if (fabs(new_norm - lambda_old) < tol * fabs(new_norm)) return new_norm;
        lambda_old = new_norm;
    }
    return lambda_old;
}

// ─── Main ────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    int max_N = 30;   // for exact enumeration
    int max_W = 30;   // for transfer matrix
    if (argc > 1) max_N = atoi(argv[1]);
    if (argc > 2) max_W = atoi(argv[2]);

    printf("Sum-Free Transfer Matrix — PMF Salvo for Erdős #166\n");
    printf("====================================================\n\n");

    // ── Part 1: Exact Partition Function ──
    printf("PART 1: Exact Partition Function Z(N)\n");
    printf("--------------------------------------\n");
    printf("%4s | %12s | %8s | %12s | %10s | %10s\n",
           "N", "total_SF", "max|A|", "max|A|/N", "log2(#SF)", "time(s)");
    printf("-----+--------------+----------+--------------+------------+----------\n");

    FILE* json_out = fopen("EXP-MATH-ERDOS166-SUMFREE-001_RESULTS.json", "w");
    fprintf(json_out, "{\n  \"experiment_id\": \"EXP-MATH-ERDOS166-SUMFREE-001\",\n");
    fprintf(json_out, "  \"title\": \"Sum-Free Partition Function & Transfer Matrix\",\n");
    fprintf(json_out, "  \"problem\": \"Erdos #166: Sum-free sets\",\n");
    fprintf(json_out, "  \"known_answer\": \"max density = 1/2, count ~ 2^(N/2)\",\n");
    fprintf(json_out, "  \"partition_function\": [\n");

    bool first = true;
    for (int N = 3; N <= max_N; N++) {
        auto t0 = chrono::high_resolution_clock::now();
        SumFreeEnumerator sfe;
        sfe.enumerate(N);
        auto t1 = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration<double>(t1 - t0).count();

        double density = (double)sfe.max_size / N;
        double log2_count = log2((double)sfe.total_count);

        printf("%4d | %12lld | %8d | %12.6f | %10.4f | %10.3f\n",
               N, sfe.total_count, sfe.max_size, density, log2_count, elapsed);
        fflush(stdout);

        if (!first) fprintf(json_out, ",\n");
        first = false;
        fprintf(json_out, "    {\"N\": %d, \"total_sumfree\": %lld, \"max_size\": %d, "
                "\"density\": %.8f, \"log2_count\": %.6f, \"time_s\": %.4f, "
                "\"density_of_states\": [",
                N, sfe.total_count, sfe.max_size, density, log2_count, elapsed);
        for (int k = 0; k <= sfe.max_size; k++) {
            if (k > 0) fprintf(json_out, ", ");
            fprintf(json_out, "%lld", sfe.density_of_states[k]);
        }
        fprintf(json_out, "]}");

        if (elapsed > 120) {
            printf("  Step took %.0fs — stopping enumeration\n", elapsed);
            break;
        }
    }
    fprintf(json_out, "\n  ],\n");

    // ── Part 2: Local Transfer Matrix ──
    printf("\nPART 2: Local Sum-Free Transfer Matrix\n");
    printf("---------------------------------------\n");
    printf("%4s | %10s | %10s | %12s | %12s | %8s | %10s\n",
           "W", "#states", "#edges", "lambda_max", "lambda_2", "gap_rat", "time(s)");
    printf("-----+------------+------------+--------------+--------------+----------+----------\n");

    fprintf(json_out, "  \"transfer_matrix\": [\n");
    first = true;

    for (int W = 3; W <= max_W; W++) {
        auto t0 = chrono::high_resolution_clock::now();

        LocalSumFreeStates lsf;
        lsf.enumerate(W);
        int n_states = (int)lsf.states.size();

        if (n_states > 5000000) {
            printf("  W=%d: %d states — too large, stopping\n", W, n_states);
            break;
        }

        SparseMatrix A;
        A.build(W, lsf.states, lsf.state_index, lsf);
        int n_edges = (int)A.col_idx.size();

        vector<double> v1(n_states);
        double lambda_max = power_iteration(A, v1);
        double lambda_2 = second_eigenvalue(A, v1, lambda_max);
        double gap_ratio = (lambda_max > 0) ? lambda_2 / lambda_max : 0;

        auto t1 = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration<double>(t1 - t0).count();

        printf("%4d | %10d | %10d | %12.8f | %12.8f | %8.4f | %10.2f\n",
               W, n_states, n_edges, lambda_max, lambda_2, gap_ratio, elapsed);
        fflush(stdout);

        if (!first) fprintf(json_out, ",\n");
        first = false;
        fprintf(json_out, "    {\"W\": %d, \"n_states\": %d, \"n_edges\": %d, "
                "\"lambda_max\": %.12f, \"lambda_2\": %.12f, \"gap_ratio\": %.8f, "
                "\"free_energy_per_site\": %.12f, \"time_s\": %.3f}",
                W, n_states, n_edges, lambda_max, lambda_2, gap_ratio,
                log(lambda_max) / W, elapsed);

        if (elapsed > 120) {
            printf("  Step took %.0fs — stopping\n", elapsed);
            break;
        }
    }

    fprintf(json_out, "\n  ]\n}\n");
    fclose(json_out);

    printf("\nResults saved to EXP-MATH-ERDOS166-SUMFREE-001_RESULTS.json\n");
    return 0;
}
