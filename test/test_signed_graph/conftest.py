"""
Pytest configuration and fixtures for SignedGraph tests.
"""
import pytest
import networkx as nx
import numpy as np
from pathlib import Path
import tempfile
import shutil


@pytest.fixture
def small_graph():
    """Create a small test graph (complete graph with 5 nodes)."""
    return nx.complete_graph(5)


@pytest.fixture
def cycle_graph():
    """Create a cycle graph with 6 nodes."""
    return nx.cycle_graph(6)


@pytest.fixture
def path_graph():
    """Create a path graph with 5 nodes."""
    return nx.path_graph(5)


@pytest.fixture
def random_seed():
    """Provide a fixed random seed for reproducibility."""
    return 42


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    tmpdir = tempfile.mkdtemp()
    yield Path(tmpdir)
    shutil.rmtree(tmpdir)


@pytest.fixture(params=['numpy', 'scipy'])
def backend_name(request):
    """Parametrize tests over available backends."""
    return request.param


@pytest.fixture(params=['cupy'], scope='session')
def gpu_backend_name(request):
    """Parametrize tests over GPU backends (skip if not available)."""
    try:
        import cupy
        return request.param
    except ImportError:
        pytest.skip("CuPy not available")


@pytest.fixture
def pflip_values():
    """Provide common pflip values for testing."""
    return [0.0, 0.2, 0.5, 1.0]
