/**
 * B_h[g] Transfer Matrix — PMF Salvo for Erdős #755
 *
 * Prima Materia Framework, Layer 2: Arithmetic Phase
 * Direct generalization of Sidon (B_2[1]) to B_h[g] sequences.
 *
 * B_h[g] set: a set A ⊂ Z where every integer has at most g representations
 * as an (unordered) sum of h elements of A.
 *
 * Equivalently for h=2 (B_2[g]): every difference d = a - b (a > b, a,b ∈ A)
 * appears at most g times. This is translation-invariant, so the sliding
 * window transfer matrix works perfectly.
 *
 * Special cases:
 *   B_2[1] = Sidon sets (Erdős #30) — already computed, λ_max data to W=47
 *   B_2[2] = each difference appears at most twice
 *   B_2[g] for general g
 *
 * Key open question: what is the maximum density of B_2[2] sets?
 * Known: |A| ≤ c·N^{1/2} for B_2[1] (Sidon), and conjectured
 * |A| ~ c_g · N^{1/2} for all B_2[g] (same exponent, larger constant).
 *
 * Architecture: identical to sidon_transfer_matrix.cpp but with
 * diff_count[d] ≤ g instead of diff_used[d] = bool.
 *
 * Compile: g++ -O3 -march=native -o bh_g_tm bh_g_transfer_matrix.cpp
 * Run:     ./bh_g_tm [max_W] [g] [h]
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

// Global parameter
int G_PARAM = 2;  // max repetitions per difference (g in B_2[g])
int H_PARAM = 2;  // h-fold sums (only h=2 implemented for now)

// ─── State Enumeration ───────────────────────────────────────────────────
// Enumerate all B_2[g] subsets of {0,...,W-1}:
// each pairwise difference appears at most g times.

struct BhgEnumerator {
    int W;
    int g;  // max repetitions
    vector<u64> states;
    unordered_map<u64, int> state_index;

    void enumerate(int max_W, int max_g) {
        W = max_W;
        g = max_g;
        states.clear();
        state_index.clear();

        // Empty set is always valid
        states.push_back(0ULL);
        state_index[0ULL] = 0;

        vector<int> current;
        vector<int> diff_count(W, 0);  // diff_count[d] = # times difference d appears

        backtrack(0, current, diff_count);
    }

    void backtrack(int start, vector<int>& current, vector<int>& diff_count) {
        if (!current.empty()) {
            u64 mask = 0;
            for (int p : current) mask |= (1ULL << p);
            int idx = (int)states.size();
            states.push_back(mask);
            state_index[mask] = idx;
        }

        for (int x = start; x < W; x++) {
            // Compute new differences x - a for each a in current
            bool valid = true;
            vector<int> new_diffs;
            for (int a : current) {
                int d = x - a;
                // Check: adding this difference would exceed g?
                if (diff_count[d] + 1 > g) {
                    valid = false;
                    break;
                }
                new_diffs.push_back(d);
            }

            // Check: do new_diffs themselves create >g repetitions?
            // Count occurrences among new_diffs
            if (valid && new_diffs.size() > 1) {
                // Sort and check for duplicates exceeding g
                vector<int> sorted_diffs = new_diffs;
                sort(sorted_diffs.begin(), sorted_diffs.end());
                int run = 1;
                for (int i = 1; i < (int)sorted_diffs.size(); i++) {
                    if (sorted_diffs[i] == sorted_diffs[i - 1]) {
                        run++;
                        if (diff_count[sorted_diffs[i]] + run > g) {
                            valid = false;
                            break;
                        }
                    } else {
                        run = 1;
                    }
                }
            }

            if (valid) {
                for (int d : new_diffs) diff_count[d]++;
                current.push_back(x);
                backtrack(x + 1, current, diff_count);
                current.pop_back();
                for (int d : new_diffs) diff_count[d]--;
            }
        }
    }
};

// ─── Sparse Transfer Matrix ─────────────────────────────────────────────

struct SparseMatrix {
    int n;
    vector<int> row_ptr;
    vector<int> col_idx;

    void build(int W, int g, const vector<u64>& states, const unordered_map<u64, int>& state_index) {
        n = (int)states.size();
        row_ptr.resize(n + 1, 0);
        col_idx.clear();

        for (int i = 0; i < n; i++) {
            u64 state_i = states[i];
            row_ptr[i] = (int)col_idx.size();

            // Shift right: position p becomes p-1, position 0 drops
            u64 shifted = (state_i >> 1);

            // Compute positions and difference counts of shifted state
            vector<int> shifted_positions;
            for (int p = 0; p < W - 1; p++) {
                if (shifted & (1ULL << p)) shifted_positions.push_back(p);
            }

            // Compute difference counts
            vector<int> shifted_diff_count(W, 0);
            for (int a = 0; a < (int)shifted_positions.size(); a++) {
                for (int b = a + 1; b < (int)shifted_positions.size(); b++) {
                    int d = shifted_positions[b] - shifted_positions[a];
                    shifted_diff_count[d]++;
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
            vector<int> new_diffs;
            for (int p : shifted_positions) {
                int d = new_pos - p;
                new_diffs.push_back(d);
            }

            // Check: each existing diff + new contribution doesn't exceed g
            // Count new_diffs occurrences
            vector<int> new_diff_count(W, 0);
            for (int d : new_diffs) {
                if (d < W) new_diff_count[d]++;
            }
            for (int d = 0; d < W && valid; d++) {
                if (shifted_diff_count[d] + new_diff_count[d] > g) {
                    valid = false;
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
    int max_W = 30;
    int g = 2;
    if (argc > 1) max_W = atoi(argv[1]);
    if (argc > 2) g = atoi(argv[2]);
    if (max_W > 63) max_W = 63;
    G_PARAM = g;

    printf("B_2[%d] Transfer Matrix — PMF Salvo for Erdős #755\n", g);
    printf("===================================================\n\n");
    printf("Constraint: each pairwise difference appears at most %d times\n", g);
    printf("B_2[1] = Sidon (Erdős #30), B_2[2] = first open case\n\n");

    // Run for multiple g values to show the progression
    for (int cur_g = 1; cur_g <= g; cur_g++) {
        printf("─── g = %d ───\n", cur_g);
        printf("%4s | %10s | %10s | %12s | %12s | %8s | %10s\n",
               "W", "#states", "#edges", "lambda_max", "lambda_2", "gap_rat", "time(s)");
        printf("-----+------------+------------+--------------+--------------+----------+----------\n");

        char fname[256];
        snprintf(fname, sizeof(fname), "EXP-MATH-ERDOS755-B2G%d-001_RESULTS.json", cur_g);
        FILE* json_out = fopen(fname, "w");
        fprintf(json_out, "{\n  \"experiment_id\": \"EXP-MATH-ERDOS755-B2G%d-001\",\n", cur_g);
        fprintf(json_out, "  \"title\": \"B_2[%d] Transfer Matrix\",\n", cur_g);
        fprintf(json_out, "  \"problem\": \"Erdos #755: B_h[g] sequences\",\n");
        fprintf(json_out, "  \"parameters\": {\"h\": 2, \"g\": %d},\n", cur_g);
        fprintf(json_out, "  \"data\": [\n");

        bool first = true;
        for (int W = 3; W <= max_W; W++) {
            auto t0 = chrono::high_resolution_clock::now();

            BhgEnumerator be;
            be.enumerate(W, cur_g);
            int n_states = (int)be.states.size();

            if (n_states > 5000000) {
                printf("  W=%d: %d states — too large, stopping\n", W, n_states);
                break;
            }

            SparseMatrix A;
            A.build(W, cur_g, be.states, be.state_index);
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
        printf("\n");
    }

    printf("Results saved to EXP-MATH-ERDOS755-B2G*-001_RESULTS.json\n");
    return 0;
}
