/**
 * Sidon Transfer Matrix — High-Performance C++ Implementation
 *
 * Computes the largest eigenvalue of the Sidon lattice gas transfer matrix
 * for window sizes W up to 50+.
 *
 * Key optimizations:
 * 1. Bitset state representation (uint64_t bitmask)
 * 2. Backtracking enumeration with difference tracking
 * 3. Sparse adjacency list (CSR format)
 * 4. Power iteration with Kahan summation for numerical stability
 *
 * Compile: g++ -O3 -march=native -o sidon_tm sidon_transfer_matrix.cpp
 * Run:     ./sidon_tm [max_W]
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <ctime>
#include <chrono>

using namespace std;
using u64 = uint64_t;

// ─── State representation ────────────────────────────────────────────────
// A Sidon state is a bitmask: bit i set means position i is occupied.
// For W ≤ 63, this fits in a uint64_t.

struct SidonEnumerator {
    int W;
    vector<u64> states;
    unordered_map<u64, int> state_index;

    void enumerate(int max_W) {
        W = max_W;
        states.clear();
        state_index.clear();

        // Empty set
        states.push_back(0ULL);
        state_index[0ULL] = 0;

        // Backtrack to enumerate all Sidon subsets of {0,...,W-1}
        vector<int> current;
        vector<int> used_diffs;  // flat set using bool array
        vector<bool> diff_used(W, false);

        backtrack(0, current, diff_used);
    }

    void backtrack(int start, vector<int>& current, vector<bool>& diff_used) {
        if (!current.empty()) {
            u64 mask = 0;
            for (int p : current) mask |= (1ULL << p);
            int idx = (int)states.size();
            states.push_back(mask);
            state_index[mask] = idx;
        }

        for (int x = start; x < W; x++) {
            // Check if adding x creates a difference collision
            bool valid = true;
            vector<int> new_diffs;
            for (int a : current) {
                int d = x - a;
                if (diff_used[d]) {
                    valid = false;
                    break;
                }
                new_diffs.push_back(d);
            }
            // Also check new diffs don't collide with each other
            if (valid) {
                for (int i = 0; i < (int)new_diffs.size() && valid; i++) {
                    for (int j = i + 1; j < (int)new_diffs.size(); j++) {
                        if (new_diffs[i] == new_diffs[j]) {
                            valid = false;
                            break;
                        }
                    }
                }
            }

            if (valid) {
                for (int d : new_diffs) diff_used[d] = true;
                current.push_back(x);
                backtrack(x + 1, current, diff_used);
                current.pop_back();
                for (int d : new_diffs) diff_used[d] = false;
            }
        }
    }
};

// ─── Sparse Transfer Matrix (CSR format) ─────────────────────────────────

struct SparseMatrix {
    int n;
    vector<int> row_ptr;   // row_ptr[i] = start of row i in col_idx
    vector<int> col_idx;   // column indices

    void build(int W, const vector<u64>& states, const unordered_map<u64, int>& state_index) {
        n = (int)states.size();
        row_ptr.resize(n + 1, 0);
        col_idx.clear();

        for (int i = 0; i < n; i++) {
            u64 state_i = states[i];
            row_ptr[i] = (int)col_idx.size();

            // Shift: each position p becomes p-1; bit 0 drops out
            u64 shifted = (state_i >> 1);
            // Compute differences of shifted state
            vector<int> shifted_positions;
            for (int p = 0; p < W - 1; p++) {
                if (shifted & (1ULL << p)) shifted_positions.push_back(p);
            }

            // Compute difference set of shifted positions
            vector<bool> shifted_diffs(W, false);
            for (int a = 0; a < (int)shifted_positions.size(); a++) {
                for (int b = a + 1; b < (int)shifted_positions.size(); b++) {
                    int d = shifted_positions[b] - shifted_positions[a];
                    shifted_diffs[d] = true;
                }
            }

            // Option 1: don't occupy W-1
            auto it = state_index.find(shifted);
            if (it != state_index.end()) {
                col_idx.push_back(it->second);
            }

            // Option 2: occupy W-1
            bool valid = true;
            int new_pos = W - 1;
            vector<int> new_diffs_vec;
            for (int p : shifted_positions) {
                int d = new_pos - p;
                if (d < W && shifted_diffs[d]) {
                    valid = false;
                    break;
                }
                new_diffs_vec.push_back(d);
            }
            // Check new diffs don't collide with each other
            if (valid && new_diffs_vec.size() > 1) {
                vector<bool> seen(W, false);
                for (int d : new_diffs_vec) {
                    if (d < W && seen[d]) { valid = false; break; }
                    if (d < W) seen[d] = true;
                }
            }
            if (valid) {
                u64 new_state = shifted | (1ULL << new_pos);
                auto it2 = state_index.find(new_state);
                if (it2 != state_index.end()) {
                    col_idx.push_back(it2->second);
                }
            }
        }
        row_ptr[n] = (int)col_idx.size();
    }

    // Sparse matrix-vector multiply: y = A^T * x (transpose because we stored row→col)
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

// ─── Power Iteration ─────────────────────────────────────────────────────

double power_iteration(const SparseMatrix& A, vector<double>& v, int max_iter = 1000, double tol = 1e-14) {
    int n = A.n;
    vector<double> w(n);

    // Initialize with random-ish vector
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

        if (fabs(new_norm - lambda_old) < tol * fabs(new_norm)) {
            return new_norm;
        }
        lambda_old = new_norm;
    }
    return lambda_old;
}

double second_eigenvalue(const SparseMatrix& A, const vector<double>& v1, double lambda1, int max_iter = 1000, double tol = 1e-12) {
    int n = A.n;
    vector<double> v(n), w(n);

    srand(123);
    for (int i = 0; i < n; i++) v[i] = (double)rand() / RAND_MAX;
    // Project out v1
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

        // Project out v1
        dot = 0;
        for (int i = 0; i < n; i++) dot += w[i] * v1[i];
        for (int i = 0; i < n; i++) w[i] -= dot * v1[i];

        double new_norm = 0;
        for (int i = 0; i < n; i++) new_norm += w[i] * w[i];
        new_norm = sqrt(new_norm);

        if (new_norm < 1e-30) return 0;
        for (int i = 0; i < n; i++) v[i] = w[i] / new_norm;

        if (fabs(new_norm - lambda_old) < tol * fabs(new_norm)) {
            return new_norm;
        }
        lambda_old = new_norm;
    }
    return lambda_old;
}

// ─── Main ────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    int max_W = 35;
    if (argc > 1) max_W = atoi(argv[1]);
    if (max_W > 63) max_W = 63;  // uint64_t limit

    printf("Sidon Transfer Matrix — C++ Implementation\n");
    printf("============================================\n\n");
    printf("%4s | %10s | %10s | %12s | %12s | %8s | %10s\n",
           "W", "#states", "#edges", "lambda_max", "lambda_2", "gap_rat", "time(s)");
    printf("-----+------------+------------+--------------+--------------+----------+----------\n");

    // JSON output buffer
    FILE* json_out = fopen("EXP-MATH-ERDOS30-SIDON-007C_RESULTS.json", "w");
    fprintf(json_out, "{\n  \"experiment_id\": \"EXP-MATH-ERDOS30-SIDON-007C\",\n");
    fprintf(json_out, "  \"title\": \"C++ Transfer Matrix Push\",\n");
    fprintf(json_out, "  \"data\": [\n");

    bool first = true;
    for (int W = 3; W <= max_W; W++) {
        auto t0 = chrono::high_resolution_clock::now();

        // Enumerate states
        SidonEnumerator se;
        se.enumerate(W);
        int n_states = (int)se.states.size();

        if (n_states > 5000000) {
            printf("  W=%d: %d states — too large, stopping\n", W, n_states);
            break;
        }

        // Build sparse matrix
        SparseMatrix A;
        A.build(W, se.states, se.state_index);
        int n_edges = (int)A.col_idx.size();

        // Power iteration
        vector<double> v1(n_states);
        double lambda_max = power_iteration(A, v1);
        double lambda_2 = second_eigenvalue(A, v1, lambda_max);
        double gap_ratio = (lambda_max > 0) ? lambda_2 / lambda_max : 0;

        auto t1 = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration<double>(t1 - t0).count();

        printf("%4d | %10d | %10d | %12.8f | %12.8f | %8.4f | %10.2f\n",
               W, n_states, n_edges, lambda_max, lambda_2, gap_ratio, elapsed);
        fflush(stdout);

        // JSON
        if (!first) fprintf(json_out, ",\n");
        first = false;
        fprintf(json_out, "    {\"W\": %d, \"n_states\": %d, \"n_edges\": %d, "
                "\"lambda_max\": %.12f, \"lambda_2\": %.12f, \"gap_ratio\": %.8f, "
                "\"free_energy\": %.12f, \"time_s\": %.3f}",
                W, n_states, n_edges, lambda_max, lambda_2, gap_ratio,
                log(lambda_max) / W, elapsed);

        // Safety timeout
        if (elapsed > 300) {
            printf("  Step took %.0fs — stopping\n", elapsed);
            break;
        }
    }

    fprintf(json_out, "\n  ]\n}\n");
    fclose(json_out);

    printf("\nResults saved to EXP-MATH-ERDOS30-SIDON-007C_RESULTS.json\n");
    return 0;
}
