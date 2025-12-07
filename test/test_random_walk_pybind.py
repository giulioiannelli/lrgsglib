from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csr_matrix

os.environ.setdefault(
    "MYLIB_LOG_FILE", str(Path(__file__).parent / "tmp_lrgsglib.log")
)

from lrgsglib.bindings.random_walk import signed_random_walk_csr


def test_random_walk_deterministic_sign_flip() -> None:
    data = np.array([-1.0, -1.0], dtype=np.float64)
    indices = np.array([1, 0], dtype=np.int32)
    indptr = np.array([0, 1, 2], dtype=np.int32)
    adjacency = csr_matrix((data, indices, indptr), shape=(2, 2))

    signs = signed_random_walk_csr(adjacency, num_steps=4, start_node=0, seed=7)

    assert signs.tolist() == [-1, 1, -1, 1]


def test_random_walk_handles_isolated_node() -> None:
    adjacency = csr_matrix((2, 2), dtype=np.float64)  # no edges
    signs = signed_random_walk_csr(adjacency, num_steps=3, start_node=1, seed=0)
    assert np.all(signs == 1)


def test_random_walk_invalid_start_node_raises() -> None:
    adjacency = csr_matrix((1, 1), dtype=np.float64)
    with pytest.raises(ValueError):
        signed_random_walk_csr(adjacency, num_steps=1, start_node=5, seed=0)
