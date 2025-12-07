"""Python wrapper for the C++ signed random walk.

Example
-------
>>> import numpy as np
>>> from scipy.sparse import csr_matrix
>>> from lrgsglib.bindings.random_walk import signed_random_walk_csr
>>> data = np.array([-1.0, -1.0])
>>> indices = np.array([1, 0])
>>> indptr = np.array([0, 1, 2])
>>> adj = csr_matrix((data, indices, indptr), shape=(2, 2))
>>> signed_random_walk_csr(adj, num_steps=4, start_node=0, seed=123)
array([-1,  1, -1,  1], dtype=int8)
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.sparse import csr_matrix

from lrgsglib import _random_walk

__all__ = ["signed_random_walk_csr"]


def signed_random_walk_csr(
    adjacency: csr_matrix,
    num_steps: int,
    start_node: int = 0,
    seed: Optional[int] = None,
) -> np.ndarray:
    """Run a signed random walk using the C++ backend on a CSR adjacency.

    Parameters
    ----------
    adjacency:
        Signed adjacency matrix in CSR form (weights may be positive or negative).
    num_steps:
        Number of steps to take.
    start_node:
        Starting node index (default: 0).
    seed:
        Optional RNG seed for reproducibility.
    """
    if num_steps < 0:
        raise ValueError("num_steps must be non-negative")

    csr = adjacency.tocsr().astype(np.float64, copy=False)
    data = np.ascontiguousarray(csr.data, dtype=np.float64)
    indices = np.ascontiguousarray(csr.indices, dtype=np.int32)
    indptr = np.ascontiguousarray(csr.indptr, dtype=np.int32)

    signs = _random_walk.random_walk(
        data, indices, indptr, int(num_steps), int(start_node), seed
    )
    return np.asarray(signs, dtype=np.int8)
