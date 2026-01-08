import numpy as np
import networkx as nx
import argparse
from tqdm import tqdm
from numba import njit
from multiprocessing import Pool, cpu_count
import os

# ----------------------------------------------
# Build 2D lattice with periodic boundaries using NetworkX
# ----------------------------------------------
def build_lattice(L, frac_neg=0.0):
    """
    Create a 2D periodic lattice using NetworkX.

    Parameters:
    -----------
    L : int
        Linear size of the lattice
    frac_neg : float
        Fraction of edges to flip to -1 (default: 0.0, unsigned lattice)

    Returns:
    --------
    neigh : array
        Neighbor list (N, 4) with node indices
    weights : array
        Edge weights (N, 4), +1 or -1 for each edge
    """
    G = nx.grid_2d_graph(L, L, periodic=True)

    # Set all edge weights to +1.0 initially
    for edge in G.edges():
        G[edge[0]][edge[1]]['weight'] = 1.0

    # Randomly flip a fraction of edges to -1
    if frac_neg > 0.0:
        num_edges = G.number_of_edges()
        num_edges_to_flip = int(frac_neg * num_edges)
        edge_list = list(G.edges())

        # Sample edges to flip
        edges_to_flip = np.random.choice(
            num_edges, size=num_edges_to_flip, replace=False
        )

        # Flip weights for sampled edges
        for edge_idx in edges_to_flip:
            edge = edge_list[edge_idx]
            G[edge[0]][edge[1]]['weight'] = -1.0

    # Create node mapping from (x,y) coordinates to linear indices
    node_to_idx = {
        node: x * L + y for x, y in G.nodes() for node in [(x, y)]
    }

    # Build neighbor list and weights from the weighted graph
    N = L * L
    neigh = np.zeros((N, 4), dtype=np.int32)
    weights = np.ones((N, 4), dtype=np.float64)

    for node in G.nodes():
        i = node_to_idx[node]
        neighbors = list(G.neighbors(node))
        for k, neighbor in enumerate(neighbors):
            j = node_to_idx[neighbor]
            neigh[i, k] = j
            weights[i, k] = G[node][neighbor]['weight']

    return neigh, weights

# ----------------------------------------------
# Single Monte Carlo sweep (N single-site updates)
# ReLU-clipped activation:
#       P = min(1, max(0, Lambda))
# Lambda = gamma * (#active neighbors) / k
# ----------------------------------------------
@njit
def sweep(state, lambda_arr, neigh, weights, gamma, k, N):
    for _ in range(N):
        i = np.random.randint(N)

        # Activation probability
        P = lambda_arr[i]
        if P < 0.0: P = 0.0
        elif P > 1.0: P = 1.0

        old = state[i]
        new = 1 if np.random.random() < P else 0

        if new != old:
            state[i] = new
            delta = new - old

            # Update lambda for neighbors with signed weights
            for kk in range(4):
                j = neigh[i, kk]
                # Find the reverse edge weight from j to i
                w_ji = 0.0
                for kkk in range(4):
                    if neigh[j, kkk] == i:
                        w_ji = weights[j, kkk]
                        break
                lambda_arr[j] += gamma * w_ji * delta / k

# ----------------------------------------------
# Full simulation (core loop)
# ----------------------------------------------
@njit
def _run_simulation_core(
    state, lambda_arr, neigh, weights, gamma, k, sweeps, N
):
    """Numba-compiled core simulation loop."""
    dens = np.zeros(sweeps, dtype=np.float64)

    for t in range(sweeps):
        sweep(state, lambda_arr, neigh, weights, gamma, k, N)
        dens[t] = state.mean()
        if dens[t] == 0.0 or dens[t] == 1.0:
            if dens[t] == 1.0:
                dens[t:] = 1.0
            break

    return dens

# ----------------------------------------------
# Initialize state and lambda arrays
# ----------------------------------------------
def initialize_state_and_lambda(N, neigh, weights, gamma, k):
    """
    Initialize random state and compute initial lambda values.

    Parameters:
    -----------
    N : int
        Number of nodes
    neigh : array
        Neighbor list (N, 4) with node indices
    weights : array
        Edge weights (N, 4), +1 or -1 for each edge
    gamma : float
        Infection/activation rate
    k : float
        Normalization factor for lambda

    Returns:
    --------
    state : array
        Initial random configuration (0 or 1 for each node)
    lambda_arr : array
        Initial lambda values computed from neighbor states
    """
    # Initial random configuration
    state = np.random.randint(0, 2, size=N, dtype=np.int8)

    # Initialize lambda exactly with signed weights
    lambda_arr = np.zeros(N, dtype=np.float64)
    for i in range(N):
        ssum = 0.0
        for kk in range(4):
            j = neigh[i, kk]
            ssum += weights[i, kk] * state[j]
        lambda_arr[i] = gamma * ssum / k

    return state, lambda_arr

# ----------------------------------------------
# Full simulation
# ----------------------------------------------
def run_simulation(L, gamma, k=2.0, sweeps=100000, frac_neg=0.0):
    """
    Run contact process simulation on 2D periodic lattice.

    Parameters:
    -----------
    L : int
        Linear size of the lattice
    gamma : float
        Infection/activation rate
    k : float, optional
        Normalization factor for lambda (default: 2.0)
    sweeps : int, optional
        Number of Monte Carlo sweeps (default: 100000)
    frac_neg : float, optional
        Fraction of edges to flip to -1 (default: 0.0)
    """
    N = L * L
    neigh, weights = build_lattice(L, frac_neg)

    # Initialize state and lambda
    state, lambda_arr = initialize_state_and_lambda(
        N, neigh, weights, gamma, k
    )

    # Run simulation with numba-compiled core
    dens = _run_simulation_core(
        state, lambda_arr, neigh, weights, gamma, k, sweeps, N
    )

    return dens, state

# ----------------------------------------------
# Wrapper for parallel execution
# ----------------------------------------------
def run_simulation_wrapper(args_tuple):
    """Wrapper function for multiprocessing pool."""
    L, gamma, k, sweeps, frac_neg = args_tuple
    dens, _ = run_simulation(L, gamma, k, sweeps, frac_neg)
    return dens

# ----------------------------------------------
# Main
# ----------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            'Contact process simulation on 2D periodic lattice '
            'with ReLU activation'
        )
    )
    parser.add_argument('L', type=int,
                        help='Linear size of the lattice')
    parser.add_argument('gamma', type=float,
                        help='Infection/activation rate')
    parser.add_argument(
        '-k', '--k', type=float, default=2.0,
        help='Normalization factor for lambda (default: 2.0)'
    )
    parser.add_argument(
        '-s', '--sweeps', type=int, default=100000,
        help='Number of Monte Carlo sweeps (default: 100000)'
    )
    parser.add_argument(
        '-na', '--na', type=int, default=1,
        help='Number of averages/realizations (default: 1)'
    )
    parser.add_argument(
        '-np', '--nproc', type=int, default=None,
        help=(
            f'Number of parallel processes '
            f'(default: all available CPUs = {cpu_count()})'
        )
    )
    parser.add_argument(
        '-o', '--outdir', type=str, default='./figures',
        help='Output directory for figures (default: ./figures)'
    )
    parser.add_argument(
        '-f', '--frac_neg', type=float, default=0.0,
        help='Fraction of edges to flip to -1 (default: 0.0, unsigned)'
    )

    args = parser.parse_args()

    # Set number of processes
    nproc = args.nproc if args.nproc is not None else cpu_count()

    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)

    # Generate filename based on parameters (high precision for floats)
    figname = (
        f"cp2d_relu_L{args.L}_g{args.gamma:.6f}_k{args.k:.6f}"
        f"_f{args.frac_neg:.6f}_s{args.sweeps}_na{args.na}.png"
    )
    figpath = os.path.join(args.outdir, figname)

    # Check if figure already exists
    if os.path.exists(figpath):
        print(f"Figure already exists: {figpath}")
        print("Skipping simulation. Delete the file to rerun.")
        import sys
        sys.exit(0)

    import matplotlib.pyplot as plt

    if args.na == 1:
        # Single run
        dens, final_state = run_simulation(
            args.L, args.gamma, args.k, args.sweeps, args.frac_neg
        )
        plt.plot(
            np.linspace(1, args.sweeps, args.sweeps), dens,
            'b-', lw=1.5
        )
    else:
        # Multiple runs with multiprocessing
        print(
            f"Running {args.na} simulations in parallel "
            f"using {nproc} processes..."
        )

        # Prepare arguments for each run
        run_args = [
            (args.L, args.gamma, args.k, args.sweeps, args.frac_neg)
            for _ in range(args.na)
        ]

        # Run simulations in parallel
        with Pool(processes=nproc) as pool:
            all_dens = list(tqdm(
                pool.imap(run_simulation_wrapper, run_args),
                total=args.na,
                desc="Averaging runs"
            ))

        # Plot individual trajectories
        for dens in all_dens[:20]:
            plt.plot(np.linspace(1, args.sweeps, args.sweeps), dens,
                    color='lightblue', alpha=0.5, lw=0.8)

        # Compute and plot average
        avg_dens = np.mean(all_dens, axis=0)
        plt.plot(np.linspace(1, args.sweeps, args.sweeps), avg_dens,
                'b-', lw=2, label=f'Average (n={args.na})')
        plt.legend()

    # Plot the density trajectory
    plt.axhline(
        y=1/(args.L**2), color='red', ls=':', label='1/N (absorbing)'
    )
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Sweep")
    plt.ylabel("Density")
    title = (
        f"Density trajectory ({args.L}x{args.L}, "
        f"gamma={args.gamma}, k={args.k}"
    )
    if args.frac_neg > 0:
        title += f", f_neg={args.frac_neg:.2f}"
    title += f", na={args.na})"
    plt.title(title)
    plt.tight_layout()

    # Save figure
    plt.savefig(figpath, dpi=150, bbox_inches='tight')
    print(f"Figure saved: {figpath}")

    plt.show()
