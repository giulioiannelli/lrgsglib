"""
Test that all library modules can be imported without error.

Covers:
- graphs.nx: All 30+ graph types (NX suffix canonical path)
- graphs.nx.funcs: Utility functions
- statsys: IsingDynamics, ContactProcess, VoterModel, SignedRW
- utils.lrg: spectral, infocomm (classical/renyi/expm/slq), clustering
- config: Constants, environment
"""

import importlib

import pytest


# -----------------------------------------------------------------------
# Graph types (canonical lrgsglib.graphs.nx path)
# -----------------------------------------------------------------------
GRAPH_MODULES = [
    # Base
    "lrgsglib.graphs.nx.SignedGraphNX",
    # Lattice
    "lrgsglib.graphs.nx.Lattice2DNX",
    "lrgsglib.graphs.nx.Lattice3DNX",
    # Random
    "lrgsglib.graphs.nx.RandomGraphNX",
    "lrgsglib.graphs.nx.ErdosRenyiNX",
    "lrgsglib.graphs.nx.BarabasiAlbertNX",
    "lrgsglib.graphs.nx.WattsStrogatzNX",
    "lrgsglib.graphs.nx.StochasticBlockModelNX",
    "lrgsglib.graphs.nx.kRegularGraphNX",
    "lrgsglib.graphs.nx.ConfigurationModelNX",
    "lrgsglib.graphs.nx.RandomGeometricNX",
    "lrgsglib.graphs.nx.LFRBenchmarkNX",
    "lrgsglib.graphs.nx.ExtendedBarabasiAlbertNX",
    "lrgsglib.graphs.nx.DualBarabasiAlbertNX",
    "lrgsglib.graphs.nx.HolmeKimNX",
    # Complete
    "lrgsglib.graphs.nx.CompleteGraphNX",
    "lrgsglib.graphs.nx.FullyConnectedNX",
    # Neural
    "lrgsglib.graphs.nx.HofieldNNNX",
    "lrgsglib.graphs.nx.SCSGeneralizedNNNX",
    # Fractal
    "lrgsglib.graphs.nx.FractalGraphNX",
    "lrgsglib.graphs.nx.DGMgraphNX",
    "lrgsglib.graphs.nx.SierpinskiNX",
    # Bipartite
    "lrgsglib.graphs.nx.BipartiteGraphNX",
    "lrgsglib.graphs.nx.BipartiteFromDegreeSequenceNX",
    # Multispectral
    "lrgsglib.graphs.nx.MultispectralGraphNX",
    "lrgsglib.graphs.nx.MultiplicativeCascadeNX",
    "lrgsglib.graphs.nx.VicsekNX",
    "lrgsglib.graphs.nx.HierarchicalModularNX",
    # Dirac
    "lrgsglib.graphs.nx.DiracLatticeNX",
    # GraphOfGraphs
    "lrgsglib.graphs.nx.GraphOfGraphsNX",
    # Temporal
    "lrgsglib.graphs.nx.TemporalGraphNX",
    "lrgsglib.graphs.nx.TemporalSignedGraphNX",
]

# -----------------------------------------------------------------------
# Stats / dynamics modules
# -----------------------------------------------------------------------
STATSYS_MODULES = [
    "lrgsglib.statsys",
    "lrgsglib.statsys.BinDynSys",
    "lrgsglib.statsys.IsingDynamics",
    "lrgsglib.statsys.ContactProcess",
    "lrgsglib.statsys.VoterModel",
    "lrgsglib.statsys.SignedRW",
    "lrgsglib.statsys._c_backend",
]

# -----------------------------------------------------------------------
# Utility modules
# -----------------------------------------------------------------------
UTILS_MODULES = [
    "lrgsglib.utils.lrg.spectral",
    "lrgsglib.utils.lrg.clustering",
    "lrgsglib.utils.lrg.ising",
    "lrgsglib.utils.lrg.quantum_core",
    "lrgsglib.utils.lrg.infocomm",
    "lrgsglib.utils.lrg.infocomm.classical",
    "lrgsglib.utils.lrg.infocomm.renyi",
    "lrgsglib.utils.lrg.infocomm.expm",
    "lrgsglib.utils.lrg.infocomm.slq",
    "lrgsglib.utils.lrg.infocomm.distances",
]

# -----------------------------------------------------------------------
# Config modules
# -----------------------------------------------------------------------
CONFIG_MODULES = [
    "lrgsglib.config",
    "lrgsglib.config.const",
]

# -----------------------------------------------------------------------
# Funcs submodule
# -----------------------------------------------------------------------
FUNCS_MODULES = [
    "lrgsglib.graphs.nx.funcs",
]


ALL_MODULES = (
    GRAPH_MODULES
    + STATSYS_MODULES
    + UTILS_MODULES
    + CONFIG_MODULES
    + FUNCS_MODULES
)


@pytest.mark.code
@pytest.mark.parametrize("module", ALL_MODULES)
def test_import(module):
    """Verify that the module can be imported without error."""
    importlib.import_module(module)


@pytest.mark.code
def test_graphs_nx_package_exports():
    """The lrgsglib.graphs.nx package exports all expected classes."""
    import lrgsglib.graphs.nx as gnx

    expected = [
        "SignedGraphNX",
        "Lattice2DNX",
        "Lattice3DNX",
        "ErdosRenyiNX",
        "MultiplicativeCascadeGraphNX",
        "VicsekGraphNX",
        "DiracLatticeGraphNX",
        "TemporalGraphNX",
    ]
    for name in expected:
        assert hasattr(gnx, name), f"Missing export: {name}"


@pytest.mark.code
def test_statsys_package_exports():
    """The lrgsglib.statsys package exports core classes."""
    from lrgsglib.statsys import (
        IsingDynamics,
        ContactProcessEI,
        ContactProcessSIR,
        VoterModel,
    )

    assert callable(IsingDynamics)
    assert callable(ContactProcessEI)
    assert callable(ContactProcessSIR)
    assert callable(VoterModel)
