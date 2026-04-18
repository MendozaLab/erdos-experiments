import numpy as np
import networkx as nx

def spacing_ratio(eigvals):
    diffs = np.diff(eigvals)
    diffs = diffs[diffs > 1e-10]
    if len(diffs) < 2: return np.nan
    r_vals = [min(diffs[i], diffs[i+1]) / max(diffs[i], diffs[i+1]) for i in range(len(diffs)-1)]
    return np.mean(r_vals)

def generate_efl_graph(k):
    # k cliques of size k sharing exactly 1 vertex (index 0)
    G = nx.Graph()
    G.add_node(0)
    current_node = 1
    for _ in range(k):
        clique_nodes = [0] + list(range(current_node, current_node + k - 1))
        for i in range(len(clique_nodes)):
            for j in range(i+1, len(clique_nodes)):
                G.add_edge(clique_nodes[i], clique_nodes[j])
        current_node += k - 1
    return G

def killer_8_baseline_checks():
    print("--- SOP #11 / KILLER 8: EFL EXACT CHECKS ---")
    for k in [2, 3]:
        G = generate_efl_graph(k)
        A = nx.adjacency_matrix(G).todense()
        eigs = np.linalg.eigvalsh(A)
        print(f"[Analytic Check] Pure EFL (k={k}) -> Nodes: {len(G)}, Edges: {G.number_of_edges()}")
        print(f"               Unique Evals: {np.unique(np.round(eigs, 3))}")

def scan_efl_multiplex_saturation():
    print("\n--- HTTP/2 MULTIPLEX OVERLAP SATURATION ---")
    k = 12 # 12 streams multiplexed
    G_base = generate_efl_graph(k)
    A_base = nx.adjacency_matrix(G_base).todense()
    
    # Introduce stream congestion (bleeding between distinct cliques)
    eta_vals = np.linspace(0.0, 0.1, 20)
    mendoza_limit = 0.4407
    
    for eta_overlap in eta_vals:
        r_ensemble = []
        for _ in range(10):
            A = np.copy(A_base)
            n_nodes = len(G_base)
            # Add random cross-stream congestion noise
            noise = (np.random.rand(n_nodes, n_nodes) < eta_overlap).astype(float)
            noise = np.tril(noise, -1) + np.tril(noise, -1).T
            
            # Mask out existing edges
            mask = (A == 0) & (np.eye(n_nodes) == 0)
            A[mask] = noise[mask]
            
            eigs = np.linalg.eigvalsh(A)
            r = spacing_ratio(eigs)
            if not np.isnan(r):
                r_ensemble.append(r)
                
        r_mean = np.mean(r_ensemble)
        print(f"eta_overlap (Congestion %): {eta_overlap*100:4.1f}% | <r>: {r_mean:.4f}")
        if 0.435 <= r_mean <= 0.445:
            print(f"   >>> MENDOZA LIMIT (0.4407) CROSSED: Multiplexing mathematically blocked <<<")

if __name__ == "__main__":
    killer_8_baseline_checks()
    scan_efl_multiplex_saturation()
