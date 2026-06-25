"""Phase-5 cupy-everywhere spectral: cross-backend parity.

Every spectral utility under ``utils/lrg`` now threads a ``backend=`` selector
that delegates to ``BackendManager`` (numpy / scipy / cupy) instead of hardcoding
``np.linalg`` / inline cupy. These tests pin that the three backends agree on a
*fixed* input (so the comparison is well-defined; graph factories that reseed
per call are deliberately avoided here).

CuPy falls back to CPU automatically when no GPU is present, so the cupy cases
stay meaningful (and identical) on a CPU-only box.
"""

import numpy as np
import networkx as nx
import pytest

from lrgsglib.graphs._shared._backend import BackendManager
from lrgsglib.utils.lrg.spectral import (
    get_graph_lspectrum,
    get_graph_lspectrum_rw,
    compute_laplacian_properties,
)
from lrgsglib.utils.lrg.spectral_rg import spectral_frustration, spectral_rg_flow
from lrgsglib.utils.lrg.quantum import von_neumann_entropy
from lrgsglib.config.const import SG_LAPL_SIGNED, SG_LAPL_SYM, SG_LAPL_RW

pytestmark = pytest.mark.code

BACKENDS = ("numpy", "scipy", "cupy")


def _signed_graph(seed: int = 0, n: int = 36, p: float = 0.2) -> nx.Graph:
    """A fixed connected graph with deterministic +/-1 'weight' labels."""
    G = nx.gnp_random_graph(n, p, seed=seed)
    if not nx.is_connected(G):
        comps = list(nx.connected_components(G))
        for a, b in zip(comps, comps[1:]):
            G.add_edge(next(iter(a)), next(iter(b)))
    rng = np.random.default_rng(seed)
    for u, v in G.edges():
        G[u][v]["weight"] = float(rng.choice([-1.0, 1.0]))
    return G


def _sortc(w):
    return np.sort_complex(np.asarray(w, dtype=complex))


# ----------------------------------------------------------------------
# BackendManager: the eigen primitives the utilities delegate to
# ----------------------------------------------------------------------
def test_backend_eigvalsh_eigvals_eigh_parity():
    rng = np.random.default_rng(3)
    M = rng.standard_normal((20, 20))
    S = M + M.T  # symmetric
    ref_vh = BackendManager.get_backend("numpy").eigvalsh(S)
    ref_v = BackendManager.get_backend("numpy").eigvals(S)
    for b in BACKENDS:
        be = BackendManager.get_backend(b)
        assert np.allclose(be.eigvalsh(S), ref_vh, atol=1e-10)
        assert np.allclose(_sortc(be.eigvals(S)), _sortc(ref_v), atol=1e-9)
        # eigh reconstructs the matrix regardless of driver / sign conventions
        w, V = be.eigh(S)
        assert np.allclose(V @ np.diag(w) @ V.T, S, atol=1e-8)


# ----------------------------------------------------------------------
# spectral.py
# ----------------------------------------------------------------------
@pytest.mark.parametrize("ltype", [SG_LAPL_SIGNED, SG_LAPL_SYM, SG_LAPL_RW])
def test_get_graph_lspectrum_backend_parity(ltype):
    G = _signed_graph()
    ref = None
    for b in BACKENDS:
        _, w = get_graph_lspectrum(G, library=b, signed=True, laplacian_type=ltype)
        ws = _sortc(w)
        if ref is None:
            ref = ws
        else:
            assert np.allclose(ws, ref, atol=1e-8)


def test_lspectrum_rw_and_properties_backend_parity():
    G = _signed_graph()
    ref_rw = np.sort(get_graph_lspectrum_rw(G, backend="numpy")[1])
    ref_sp = np.sort(compute_laplacian_properties(G, signed=True, backend="numpy")[0])
    for b in BACKENDS:
        assert np.allclose(np.sort(get_graph_lspectrum_rw(G, backend=b)[1]),
                           ref_rw, atol=1e-10)
        assert np.allclose(
            np.sort(compute_laplacian_properties(G, signed=True, backend=b)[0]),
            ref_sp, atol=1e-10)


# ----------------------------------------------------------------------
# spectral_rg.py -- eigenvalue scalars AND derived partition labels
# ----------------------------------------------------------------------
def test_spectral_frustration_backend_parity():
    G = _signed_graph(seed=1)
    ref = spectral_frustration(G, backend="numpy")
    for b in BACKENDS:
        assert abs(spectral_frustration(G, backend=b) - ref) < 1e-10


def test_spectral_rg_flow_partition_invariant_to_backend():
    # The eigh driver differs across backends, but diffusion distance (and hence
    # the partition + frustration flow) is invariant to eigenvector sign/rotation.
    G = _signed_graph(seed=2)

    def signature(b):
        flow = spectral_rg_flow(G, n_steps=2, backend=b)
        sizes = [tuple(sorted(len(c) for c in s["partition"])) for s in flow]
        frust = [round(float(s.get("mean_boundary_frustration", 0.0)), 9)
                 for s in flow]
        return sizes, frust

    ref = signature("numpy")
    for b in BACKENDS:
        assert signature(b) == ref, f"RG flow differs for backend {b}"


# ----------------------------------------------------------------------
# quantum.py
# ----------------------------------------------------------------------
def test_von_neumann_entropy_backend_parity():
    rng = np.random.default_rng(7)
    M = rng.standard_normal((8, 8)) + 1j * rng.standard_normal((8, 8))
    rho = M @ M.conj().T
    rho = rho / np.trace(rho)
    ref = von_neumann_entropy(rho, backend="numpy")
    for b in BACKENDS:
        assert abs(von_neumann_entropy(rho, backend=b) - ref) < 1e-10
