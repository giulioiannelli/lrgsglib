"""Graph generators for the :mod:`MultispectralGraph` module."""

from __future__ import annotations

import random

import networkx as nx
import numpy as np

from typing import List, Sequence, Tuple, Dict, Any


ProbabilityMatrix = np.ndarray
Node2D = Tuple[int, int]


# Rebuild public API after all definitions
__all__ = [
    "multiplicative_cascade_probability_matrix",
    "multiplicative_cascade_graph",
    "multiplicative_cascade_exp_clocks",
    "initial_measure",
    "link_probabilities",
    "palla_lovasz_vicksek_graph",
    "dirac_comb_graph",
    "dirac_brush_graph",
]


def _select_nodes(
    graph: nx.Graph,
    probabilities: ProbabilityMatrix,
    fraction: float,
) -> List[Node2D]:
    """Randomly select nodes according to ``probabilities`` until ``fraction`` is met."""
    if not (0.0 < fraction <= 1.0):
        raise ValueError("fraction must lie in the open interval (0, 1].")

    nodes: List[Node2D] = list(graph.nodes())
    total_nodes = len(nodes)
    target = max(1, int(fraction * total_nodes))
    selected: List[Node2D] = []
    attempts = 0
    max_attempts = 10 * total_nodes if total_nodes else 0

    while len(selected) < target and nodes:
        node = random.choice(nodes)
        i, j = node
        if np.random.random() < probabilities[i, j]:
            selected.append(node)
            nodes.remove(node)
        attempts += 1
        if attempts > max_attempts:
            break

    if len(selected) < target and nodes:
        remaining = sorted(nodes, key=lambda node: probabilities[node], reverse=True)
        for node in remaining:
            selected.append(node)
            if len(selected) >= target:
                break

    if not selected:
        raise ValueError("Unable to select any node from the cascade probabilities.")

    return selected

def _validate_probabilities(probabilities: Sequence[float]) -> None:
    """Validate that probabilities are numeric and lie in ``[0, 1]``."""
    if len(probabilities) != 4:
        raise ValueError("Exactly four probabilities are required for the cascade seed.")
    for idx, value in enumerate(probabilities):
        if not isinstance(value, (int, float, np.floating)):
            raise TypeError(
                f"Probability p{idx + 1} must be a real number, received {type(value)!r}."
            )
        if not (0.0 <= float(value) <= 1.0):
            raise ValueError(
                f"Probability p{idx + 1} must be within [0, 1], received {value}."
            )


def initial_measure(m: int, symmetric: bool = True) -> np.ndarray:
    """Generate a random m x m probability matrix with sum 1.

    If ``symmetric`` is True, the matrix is symmetrized by ``(M + M.T) / 2``
    before normalizing to sum to 1.
    """
    if not isinstance(m, (int, np.integer)) or m <= 0:
        raise ValueError("m must be a positive integer.")
    M = np.random.rand(int(m), int(m)).astype(float)
    if symmetric:
        M = (M + M.T) / 2.0
    s = float(M.sum())
    if s == 0.0:
        raise ValueError("initial_measure produced a zero-sum matrix.")
    M /= s
    return M


def link_probabilities(pij: np.ndarray, k: int) -> np.ndarray:
    """Iterate the generating measure ``pij`` via Kronecker products ``k`` times."""
    if not isinstance(pij, np.ndarray) or pij.ndim != 2 or pij.shape[0] != pij.shape[1]:
        raise ValueError("pij must be a square 2D numpy array.")
    if not isinstance(k, (int, np.integer)) or k <= 0:
        raise ValueError("k must be a positive integer.")
    Pk = pij.copy().astype(float)
    for _ in range(int(k) - 1):
        Pk = np.kron(Pk, pij)
    return Pk


def multiplicative_cascade_probability_matrix(
    p1: float,
    p2: float,
    p3: float,
    p4: float,
    iterations: int = 7,
    *,
    stochastic: bool = False,
) -> ProbabilityMatrix:
    """Return the probability matrix obtained via a multiplicative cascade.

    When ``stochastic`` is False, this uses a fixed 2x2 base seed [[p1,p2],[p3,p4]].
    When ``stochastic`` is True, each block placement draws a fresh 2x2 seed
    with entries chosen from {p1,p2,p3,p4} independently.
    """
    if not isinstance(iterations, (int, np.integer)) or iterations < 0:
        raise ValueError("iterations must be a non-negative integer.")

    _validate_probabilities((p1, p2, p3, p4))

    if not stochastic:
        base = np.array([[p1, p2], [p3, p4]], dtype=float)
        matrix = base.copy()

        for _ in range(iterations):
            matrix = np.kron(matrix, base)
        return matrix

    # stochastic path
    A = np.random.choice([p1, p2, p3, p4], (2, 2)).astype(float)
    A_new = A.copy()
    for _ in range(iterations):
        A = A_new
        L_old = A.shape[0]
        L_new = L_old * 2
        if L_old <= 512:
            # Vectorized block expansion is faster but uses more memory.
            seeds = np.random.choice(
                [p1, p2, p3, p4], size=(L_old, L_old, 2, 2)
            ).astype(float)
            blocks = A[:, :, None, None] * seeds
            A_new = blocks.transpose(0, 2, 1, 3).reshape(L_new, L_new)
        else:
            A_new = np.zeros((L_new, L_new), dtype=float)
            for i in range(L_old):
                for j in range(L_old):
                    seed = np.random.choice(
                        [p1, p2, p3, p4], (2, 2)
                    ).astype(float)
                    A_new[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] += seed * A[i, j]
    return A_new



def multiplicative_cascade_graph(
    p1: float,
    p2: float,
    p3: float,
    p4: float,
    *,
    fraction: float,
    iterations: int = 7,
    stochastic: bool = False,
    probabilities: ProbabilityMatrix | None = None,
    periodic: bool = False,
) -> nx.Graph:
    """Generate a graph sampled from a multiplicative cascade distribution."""
    probabilities = (
        multiplicative_cascade_probability_matrix(
            p1, p2, p3, p4, iterations, stochastic=stochastic
        )
        if probabilities is None
        else probabilities
    )
    size = probabilities.shape[0]

    base_graph = nx.grid_2d_graph(size, size, periodic=periodic)
    selected_nodes = _select_nodes(base_graph, probabilities, fraction)
    subgraph = base_graph.subgraph(selected_nodes).copy()

    if subgraph.number_of_nodes() == 0:
        raise ValueError("The multiplicative cascade graph is empty.")

    if not nx.is_connected(subgraph):
        giant_component_nodes = max(nx.connected_components(subgraph), key=len)
        subgraph = subgraph.subgraph(giant_component_nodes).copy()

    return subgraph



def palla_lovasz_vicksek_graph(N: int, pij: np.ndarray, k: int) -> Tuple[nx.Graph, np.ndarray]:
    """Generate a Vicsek graph and return the graph and final probability matrix.

    Steps:
    1) Build the link probability matrix P_k by iterating p_ij k times.
    2) Distribute N nodes uniformly on [0,1].
    3) Connect each pair with probability given by P_k at their coordinates.
    """
    if not isinstance(N, (int, np.integer)) or N <= 0:
        raise ValueError("N must be a positive integer.")

    pij = np.asarray(pij, dtype=float)
    Pk = link_probabilities(pij, k)
    m = pij.shape[0]
    M = int(m ** int(k))

    coords = np.random.rand(int(N))
    indices = np.floor(coords * M).astype(int)
    indices = np.clip(indices, 0, M - 1)

    G = nx.Graph()
    G.add_nodes_from(range(int(N)))

    for i in range(int(N)):
        for j in range(i + 1, int(N)):
            p_link = float(Pk[indices[i], indices[j]])
            if np.random.rand() < p_link:
                G.add_edge(i, j)

    return G, Pk


def multiplicative_cascade_exp_clocks(
    adj: ProbabilityMatrix,
    fraction: float,
    periodic: bool = True,
) -> nx.Graph:
    """Sample a graph from a precomputed cascade matrix using exponential clocks.

    This implements the selection rule:
    - draw i.i.d. U ~ Uniform(0,1) for each cell
    - compute T = -log(1-U)/A_new (with T = +inf where A_new == 0)
    - select the k cells with smallest T, with k ≈ frazione * N^2
    - induce a subgraph of the 2D grid and return its giant component

    Parameters
    ----------
    adj : ProbabilityMatrix
        Square probability matrix from multiplicative cascade
    fraction : float
        Fraction of cells to select (0, 1]
    periodic : bool, default=True
        Whether to use periodic boundary conditions (toroidal grid)
    """
    if not isinstance(adj, np.ndarray) or adj.ndim != 2 or adj.shape[0] != adj.shape[1]:
        raise ValueError("A_new must be a square 2D numpy array.")
    if not (0.0 < float(fraction) <= 1.0):
        raise ValueError("frazione must lie in the open interval (0, 1].")

    N = int(adj.shape[0])
    k = int(round(float(fraction) * (N ** 2)))
    k = max(1, min(k, N * N))

    U = np.random.random((N, N))
    # Avoid division by zero: entries with zero probability get infinite time
    with np.errstate(divide="ignore", invalid="ignore"):
        time = -np.log1p(-U) / adj
        time = np.where(adj > 0.0, time, np.inf)

    flat_indices = np.argpartition(time.ravel(), k - 1)[:k]
    indices_2d = np.column_stack(np.unravel_index(flat_indices, (N, N)))
    selected_nodes = [tuple(idx) for idx in indices_2d]

    G = nx.grid_2d_graph(N, N, periodic=periodic)
    G_sub = G.subgraph(selected_nodes).copy()

    if G_sub.number_of_nodes() == 0:
        raise ValueError("No nodes selected for the cascade graph.")

    if not nx.is_connected(G_sub):
        components = nx.connected_components(G_sub)
        giant_nodes = max(components, key=len)
        return G_sub.subgraph(giant_nodes).copy()
    return G_sub


# ============================================================================
# Dirac Comb and Brush Generators
# ============================================================================


def _line(nodes: int, periodic: bool = True) -> nx.Graph:
    """Create a 1D line graph (path or cycle).

    Parameters
    ----------
    nodes : int
        Number of nodes in the line.
    periodic : bool, default=True
        If True, create a cycle (ring). If False, create a path.

    Returns
    -------
    nx.Graph
        The line graph.
    """
    if not isinstance(nodes, int) or nodes <= 0:
        raise ValueError("nodes must be a positive integer.")

    if periodic:
        return nx.cycle_graph(nodes)
    else:
        return nx.path_graph(nodes)


def _plane(x_nodes: int, y_nodes: int, periodic: bool = True) -> nx.Graph:
    """Create a 2D plane graph (grid or torus).

    Parameters
    ----------
    x_nodes : int
        Number of nodes in the x dimension.
    y_nodes : int
        Number of nodes in the y dimension.
    periodic : bool, default=True
        If True, create a torus (periodic boundaries). If False, create a grid.

    Returns
    -------
    nx.Graph
        The plane graph with nodes labeled as (i, j) tuples.
    """
    if not isinstance(x_nodes, int) or x_nodes <= 0:
        raise ValueError("x_nodes must be a positive integer.")
    if not isinstance(y_nodes, int) or y_nodes <= 0:
        raise ValueError("y_nodes must be a positive integer.")

    return nx.grid_2d_graph(x_nodes, y_nodes, periodic=periodic)


def dirac_comb_graph(
    base_nodes: int,
    fiber_nodes: int,
    base_type: str = "line",
    periodic: bool = True,
) -> Tuple[nx.Graph, Dict[str, Any]]:
    """Generate a Dirac comb graph: base network with fiber networks on each node.

    A Dirac comb consists of a base graph where each node has a fiber graph
    attached to it. The fiber is connected to the base at one anchor node.

    Parameters
    ----------
    base_nodes : int
        Number of nodes in the base graph.
    fiber_nodes : int
        Number of nodes in each fiber graph.
    base_type : str, default='line'
        Type of base graph: 'line' for 1D line/ring.
    periodic : bool, default=True
        Whether to use periodic boundary conditions for both base and fibers.

    Returns
    -------
    G : nx.Graph
        The combined Dirac comb graph.
    metadata : dict
        Dictionary containing:
        - 'base_nodes': Number of base nodes
        - 'fiber_nodes': Number of fiber nodes
        - 'total_nodes': Total number of nodes
        - 'base_graph': The base graph (separate)
        - 'fiber_graph': A prototype fiber graph (separate)
        - 'structure': 'dirac_comb'

    Notes
    -----
    Nodes are labeled as:

    - Base nodes: ('base', i) for i in range(base_nodes)
    - Fiber nodes: ('fiber', base_idx, fiber_idx)
      where base_idx is the base node and fiber_idx is position in fiber

    The fiber is attached at fiber_idx=0 to its corresponding base node.

    Examples
    --------
    >>> G, meta = dirac_comb_graph(100, 50, base_type='line', periodic=True)
    >>> print(f"Total nodes: {G.number_of_nodes()}")
    Total nodes: 5000
    >>> print(f"Base nodes: {meta['base_nodes']}, Fiber nodes: {meta['fiber_nodes']}")
    Base nodes: 100, Fiber nodes: 50
    """
    if not isinstance(base_nodes, int) or base_nodes <= 0:
        raise ValueError("base_nodes must be a positive integer.")
    if not isinstance(fiber_nodes, int) or fiber_nodes <= 0:
        raise ValueError("fiber_nodes must be a positive integer.")

    # Create base graph
    if base_type == "line":
        base = _line(base_nodes, periodic=periodic)
    else:
        raise ValueError(f"Unsupported base_type: {base_type}. Use 'line'.")

    # Create prototype fiber graph
    fiber = _line(fiber_nodes, periodic=periodic)

    # Build combined graph
    G = nx.Graph()

    # Add base nodes with new labels
    base_node_mapping = {}
    for i, node in enumerate(base.nodes()):
        new_label = ('base', i)
        base_node_mapping[node] = new_label
        G.add_node(new_label, layer='base', base_idx=i)

    # Add base edges
    for u, v in base.edges():
        G.add_edge(base_node_mapping[u], base_node_mapping[v])

    # Add fibers
    for base_idx in range(base_nodes):
        base_node = ('base', base_idx)

        # Add fiber nodes
        fiber_node_mapping = {}
        for j, fnode in enumerate(fiber.nodes()):
            fiber_label = ('fiber', base_idx, j)
            fiber_node_mapping[fnode] = fiber_label
            G.add_node(fiber_label, layer='fiber', base_idx=base_idx, fiber_idx=j)

        # Add fiber edges
        for u, v in fiber.edges():
            G.add_edge(fiber_node_mapping[u], fiber_node_mapping[v])

        # Connect fiber anchor (fiber_idx=0) to base node
        anchor_node = ('fiber', base_idx, 0)
        G.add_edge(base_node, anchor_node)

    # Metadata
    metadata = {
        'base_nodes': base_nodes,
        'fiber_nodes': fiber_nodes,
        'total_nodes': G.number_of_nodes(),
        'base_graph': base,
        'fiber_graph': fiber,
        'structure': 'dirac_comb',
        'base_type': base_type,
        'periodic': periodic,
    }

    return G, metadata


def dirac_brush_graph(
    base_x: int,
    base_y: int,
    fiber_nodes: int,
    periodic: bool = True,
) -> Tuple[nx.Graph, Dict[str, Any]]:
    """Generate a Dirac brush graph: 2D base network with 1D fiber networks.

    A Dirac brush is a special case of the Dirac comb where the base is a 2D
    plane/torus and fibers are 1D lines/rings attached to each base node.

    Parameters
    ----------
    base_x : int
        Number of nodes in the x dimension of the base.
    base_y : int
        Number of nodes in the y dimension of the base.
    fiber_nodes : int
        Number of nodes in each fiber.
    periodic : bool, default=True
        Whether to use periodic boundary conditions.

    Returns
    -------
    G : nx.Graph
        The combined Dirac brush graph.
    metadata : dict
        Dictionary containing:
        - 'base_x': Base x dimension
        - 'base_y': Base y dimension
        - 'base_nodes': Total base nodes (base_x * base_y)
        - 'fiber_nodes': Number of fiber nodes
        - 'total_nodes': Total number of nodes
        - 'base_graph': The base graph (separate)
        - 'fiber_graph': A prototype fiber graph (separate)
        - 'structure': 'dirac_brush'

    Notes
    -----
    Nodes are labeled as:
    - Base nodes: ('base', i, j) for coordinates (i, j) in the 2D grid
    - Fiber nodes: ('fiber', base_i, base_j, fiber_idx)

    The fiber is attached at fiber_idx=0 to its corresponding base node.

    Examples
    --------
    >>> G, meta = dirac_brush_graph(32, 32, 125, periodic=True)
    >>> print(f"Total nodes: {G.number_of_nodes()}")
    Total nodes: 128000
    >>> print(f"Base: {meta['base_x']}x{meta['base_y']}, Fibers: {meta['fiber_nodes']}")
    Base: 32x32, Fibers: 125
    """
    if not isinstance(base_x, int) or base_x <= 0:
        raise ValueError("base_x must be a positive integer.")
    if not isinstance(base_y, int) or base_y <= 0:
        raise ValueError("base_y must be a positive integer.")
    if not isinstance(fiber_nodes, int) or fiber_nodes <= 0:
        raise ValueError("fiber_nodes must be a positive integer.")

    # Create base 2D graph
    base = _plane(base_x, base_y, periodic=periodic)

    # Create prototype fiber graph
    fiber = _line(fiber_nodes, periodic=periodic)

    # Build combined graph
    G = nx.Graph()

    # Add base nodes with new labels
    base_node_mapping = {}
    for node in base.nodes():
        i, j = node  # node is already a tuple (i, j)
        new_label = ('base', i, j)
        base_node_mapping[node] = new_label
        G.add_node(new_label, layer='base', base_i=i, base_j=j)

    # Add base edges
    for u, v in base.edges():
        G.add_edge(base_node_mapping[u], base_node_mapping[v])

    # Add fibers
    for node in base.nodes():
        base_i, base_j = node
        base_node_label = ('base', base_i, base_j)

        # Add fiber nodes
        fiber_node_mapping = {}
        for k, fnode in enumerate(fiber.nodes()):
            fiber_label = ('fiber', base_i, base_j, k)
            fiber_node_mapping[fnode] = fiber_label
            G.add_node(fiber_label, layer='fiber', base_i=base_i, base_j=base_j, fiber_idx=k)

        # Add fiber edges
        for u, v in fiber.edges():
            G.add_edge(fiber_node_mapping[u], fiber_node_mapping[v])

        # Connect fiber anchor (fiber_idx=0) to base node
        anchor_node = ('fiber', base_i, base_j, 0)
        G.add_edge(base_node_label, anchor_node)

    # Metadata
    base_total = base_x * base_y
    metadata = {
        'base_x': base_x,
        'base_y': base_y,
        'base_nodes': base_total,
        'fiber_nodes': fiber_nodes,
        'total_nodes': G.number_of_nodes(),
        'base_graph': base,
        'fiber_graph': fiber,
        'structure': 'dirac_brush',
        'periodic': periodic,
    }

    return G, metadata
