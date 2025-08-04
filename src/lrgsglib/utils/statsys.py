from numpy.typing import NDArray
#
from .tools import UnionFind
#
__all__ = [
    "boltzmann_factor",
    "find_largest_cluster_circle2D",
    "correlated_binary_sequence_vectorized"
]
#
def boltzmann_factor(E: NDArray, T: float, k_B: float = 1.0) -> NDArray:
    """
    Compute the Boltzmann factor for a given array of energies.

    Parameters:
    -----------
    E: NDArray
        Array of energy values.
    T: float
        Temperature of the system (must be positive).
    k_B: float
        Boltzmann constant (default is 1.0).

    Returns:
    -----------
    NDArray
        Array of Boltzmann factors corresponding to the input energies.

    Notes:
    -----------
    The Boltzmann factor is calculated as exp(-E / (k_B * T)).
    """
    from numpy import exp
    if T <= 0:
        raise ValueError("Temperature must be positive.")
    return exp(-E / (k_B * T))

def find_largest_cluster_circle2D(circles, radius):
    """
    Identifies the largest cluster of overlapping circles given their centers and a common radius.

    Parameters:
    -----------
    circles (np.array): A numpy array of tuples where each tuple represents the (x, y) coordinates of a circle's center.
    radius (float): The radius of each circle.

    Returns:
    --------
    list: A list of tuples representing the centers of the circles in the largest cluster.

    Description:
    ------------
    The function utilizes a KDTree for efficient spatial queries to detect overlapping circles.
    It employs a union-find data structure to group and identify clusters of overlapping circles.
    The function returns the largest cluster found.
    """
    from scipy.spatial import KDTree
    from collections import defaultdict
    tree = KDTree(circles)
    uf = UnionFind(len(circles))
    threshold = 2 * radius

    for i in range(len(circles)):
        neighbors = tree.query_ball_point(circles[i], r=threshold)
        for j in neighbors:
            if i != j:
                uf.union(i, j)

    clusters = defaultdict(list)
    for i in range(len(circles)):
        root = uf.find(i)
        clusters[root].append(tuple(circles[i]))

    largest_cluster = max(clusters.values(), key=len)
    return largest_cluster

def correlated_binary_sequence_vectorized(length: int, T: float, J: float = 1.) -> NDArray:
    """
    Generate a random binary sequence with a given flipping probability based on interaction JI and temperature T,
    using a vectorized approach.

    Parameters
    ----------
    length : int
        Length of the sequence.
    J : float
        Interaction strength (coupling constant).
    T : float
        Temperature controlling randomness.

    Returns
    -------
    NDArray
        A binary sequence with flipping probabilities defined by tanh(JI / T).
    """
    from numpy import tanh, random, cumsum, ones
    P_flip = 1-0.5 * (tanh(J / T) + 1)  # Calculate flipping probability
    flips = random.random(size=length) < P_flip
    sequence = ones(length, dtype=int)
    sequence[0] = random.choice([-1, 1])
    flips_cumsum = cumsum(flips)
    sequence = sequence[0] * (-1) ** flips_cumsum
    return sequence

