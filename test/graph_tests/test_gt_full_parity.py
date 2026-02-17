"""Extended GT-NX parity tests for the full feature parity implementation.

Covers all methods added in Phases 0-7 of the GT full feature parity plan.
"""

import numpy as np
import pytest

from lrgsglib.graphs.nx.Lattice2DNX import Lattice2DNX
from lrgsglib.graphs.gt.SignedGraphGT import SignedGraphGT
from lrgsglib.graphs.gt._converters import nx_to_gt


# ---------- fixtures ----------

@pytest.fixture
def paired():
    """Matching NX and GT 6x6 square lattice with frustration."""
    nx_sg = Lattice2DNX(side1=6, geo="sqr", pflip=0.15, seed=99)
    nx_sg.flip_random_fract_edges()

    G_gt = nx_to_gt(nx_sg.gr["G"], sign_attr="weight")
    gt_sg = SignedGraphGT(G=G_gt, pflip=0.15, seed=99)

    assert nx_sg.N == gt_sg.N
    return nx_sg, gt_sg


@pytest.fixture
def paired_spectrum(paired):
    """Paired graphs with full eigendecomposition."""
    nx_sg, gt_sg = paired
    nx_sg.compute_laplacian_spectrum_weigV()
    gt_sg.compute_laplacian_spectrum_weigV()
    return nx_sg, gt_sg


# ==========================================================================
# Phase 0: Eigenvector convenience methods
# ==========================================================================

class TestPhase0EigvConvenience:

    def test_get_eigV_bin_check(self, paired_spectrum):
        nx_sg, gt_sg = paired_spectrum
        v_nx = nx_sg.get_eigV_bin_check(which=0)
        v_gt = gt_sg.get_eigV_bin_check(which=0)
        # Both should have {-1, +1} values
        assert set(np.unique(v_nx)).issubset({-1, 1})
        assert set(np.unique(v_gt)).issubset({-1, 1})
        assert v_nx.shape == v_gt.shape

    def test_get_eigV_bin_check_list(self, paired_spectrum):
        nx_sg, gt_sg = paired_spectrum
        vl_nx = nx_sg.get_eigV_bin_check_list(custom_slice=slice(3))
        vl_gt = gt_sg.get_eigV_bin_check_list(custom_slice=slice(3))
        assert len(vl_nx) == len(vl_gt) == 3
        for v in vl_gt:
            assert set(np.unique(v)).issubset({-1, 1})

    def test_get_eigV_bin_check_list_asarray(self, paired_spectrum):
        _, gt_sg = paired_spectrum
        arr = gt_sg.get_eigV_bin_check_list(custom_slice=slice(4), asarray=True)
        assert isinstance(arr, np.ndarray)
        assert arr.shape == (4, gt_sg.N)

    def test_load_eigV_on_graph(self, paired_spectrum):
        _, gt_sg = paired_spectrum
        gt_sg.load_eigV_on_graph(which=0, binarize=True)
        # Should not raise, eigvec should be accessible
        v = gt_sg.get_eigV_check(0, binarize=True)
        assert v.shape == (gt_sg.N,)

    def test_compute_laplacian_spectrum(self, paired):
        nx_sg, gt_sg = paired
        nx_sg.compute_laplacian_spectrum()
        gt_sg.compute_laplacian_spectrum()
        np.testing.assert_allclose(
            np.sort(nx_sg.eigv), np.sort(gt_sg.eigv), atol=1e-10
        )


# ==========================================================================
# Phase 2: Graph operations
# ==========================================================================

class TestPhase2GraphOps:

    def test_set_node_attributes(self, paired):
        _, gt_sg = paired
        gt_sg.set_node_attributes({0: 42.0, 1: -1.0}, "test")
        prop = gt_sg.G.vertex_properties["test"]
        assert prop[gt_sg.G.vertex(0)] == 42.0

    def test_load_vec_on_nodes(self, paired):
        _, gt_sg = paired
        vec = np.arange(gt_sg.N, dtype=float)
        gt_sg.load_vec_on_nodes(vec, "myvec")
        prop = gt_sg.G.vertex_properties["myvec"]
        assert prop[gt_sg.G.vertex(5)] == 5.0

    def test_set_edges_random_normal(self, paired):
        _, gt_sg = paired
        gt_sg.set_edges_random_normal(mu=0.0, sigma=2.0)
        # After this, signs are continuous floats, not just +-1
        sign_prop = gt_sg.G.edge_properties["sign"]
        vals = [sign_prop[e] for e in gt_sg.G.edges()]
        assert not all(v in (1, -1) for v in vals)

    def test_get_random_edges_from_set(self, paired):
        _, gt_sg = paired
        edges = gt_sg.get_random_edges_from_set(n=5)
        assert len(edges) == 5
        for u, v in edges:
            assert 0 <= u < gt_sg.N
            assert 0 <= v < gt_sg.N


# ==========================================================================
# Phase 3: Spectral methods
# ==========================================================================

class TestPhase3Spectral:

    def test_adjacency_spectrum_weigV(self, paired):
        nx_sg, gt_sg = paired
        nx_sg.compute_adjacency_spectrum_weigV()
        gt_sg.compute_adjacency_spectrum_weigV()
        np.testing.assert_allclose(
            np.sort(nx_sg.adj_eigv), np.sort(gt_sg.adj_eigv), atol=1e-10
        )

    def test_compute_k_adj_eigvV(self, paired):
        _, gt_sg = paired
        eigv, eigV = gt_sg.compute_k_adj_eigvV(k=3, which="LM")
        assert len(eigv) == 3
        assert eigV.shape == (gt_sg.N, 3)

    def test_make_eigV_transpose_roundtrip(self, paired_spectrum):
        _, gt_sg = paired_spectrum
        original_shape = gt_sg.eigV.shape
        gt_sg.make_eigV_transposed()
        assert gt_sg.eigV.shape == (original_shape[1], original_shape[0])
        gt_sg.make_eigV_column_major()
        assert gt_sg.eigV.shape == original_shape

    def test_make_rescaled_signed_laplacian(self, paired_spectrum):
        _, gt_sg = paired_spectrum
        gt_sg.make_rescaled_signed_laplacian(MODE="field")
        assert hasattr(gt_sg, "resLp")
        # Smallest eigenvalue of rescaled should be ~0
        rescaled_eigv = np.linalg.eigvalsh(gt_sg.resLp)
        assert abs(rescaled_eigv[0]) < 1e-10


# ==========================================================================
# Phase 4: Partitioning
# ==========================================================================

class TestPhase4Partitioning:

    def test_get_subgraph_from_nodes(self, paired):
        _, gt_sg = paired
        sub = gt_sg.get_subgraph_from_nodes([0, 1, 2, 3])
        assert sub.num_vertices() == 4

    def test_get_nodes_subgraph_by_kv(self, paired):
        _, gt_sg = paired
        vec = np.array([1.0 if i < gt_sg.N // 2 else -1.0 for i in range(gt_sg.N)])
        gt_sg.load_vec_on_nodes(vec, "half")
        g_yes, g_no = gt_sg.get_nodes_subgraph_by_kv("half", 1.0)
        assert g_yes.num_vertices() == gt_sg.N // 2
        assert g_no.num_vertices() == gt_sg.N - gt_sg.N // 2

    def test_make_connected_component_by_edge(self, paired):
        _, gt_sg = paired
        gt_sg.make_connected_component_by_edge(edge_attr="sign", value=1)
        assert hasattr(gt_sg, "largest_cc")
        assert len(gt_sg.largest_cc) > 0

    def test_get_ferroAntiferro_regions(self, paired):
        _, gt_sg = paired
        spins = np.random.choice([-1.0, 1.0], size=gt_sg.N)
        gt_sg.load_vec_on_nodes(spins, "s")
        ferro, anti = gt_sg.get_ferroAntiferro_regions(attr_str="s")
        assert isinstance(ferro, list)
        assert isinstance(anti, list)


# ==========================================================================
# Phase 5: File I/O
# ==========================================================================

class TestPhase5FileIO:

    def test_export_load_eigV_roundtrip(self, paired_spectrum, tmp_path):
        _, gt_sg = paired_spectrum
        gt_sg._ensure_paths()
        gt_sg.path_graph = tmp_path
        # Export binarized
        gt_sg.export_eigV_all(exName="rt", binarize=True)
        # Load back
        data = gt_sg.load_eigV_all(exName="rt", binarize=True)
        assert data.shape == (gt_sg.N, gt_sg.N)
        assert data.dtype == np.int8
        # Verify values
        expected = np.sign(gt_sg.eigV).astype(np.int8).T  # row-major export
        np.testing.assert_array_equal(data, expected)

    def test_export_adj_bin(self, paired, tmp_path):
        _, gt_sg = paired
        gt_sg._ensure_paths()
        gt_sg.path_graph = tmp_path
        gt_sg.export_adj_bin(exName="adj_test")
        assert gt_sg.path_exp_adj.exists()

    def test_export_edgel_bin(self, paired, tmp_path):
        _, gt_sg = paired
        gt_sg._ensure_paths()
        gt_sg.path_graph = tmp_path
        gt_sg._export_edgel_bin(exName="edgl_test")
        assert gt_sg.path_exp_edgl.exists()
        # Verify format
        dtype = [("i", np.uint64), ("j", np.uint64), ("w_ij", np.float64)]
        edges = np.fromfile(gt_sg.path_exp_edgl, dtype=dtype)
        assert len(edges) == gt_sg.num_edges

    def test_remove_exported_files(self, paired, tmp_path):
        _, gt_sg = paired
        gt_sg._ensure_paths()
        gt_sg.path_graph = tmp_path
        gt_sg._export_edgel_bin(exName="cleanup")
        assert gt_sg.path_exp_edgl.exists()
        gt_sg.remove_exported_files()
        assert not gt_sg.path_exp_edgl.exists()


# ==========================================================================
# Phase 6: Topology helpers
# ==========================================================================

class TestPhase6Topology:

    def test_nodes_in(self, paired):
        _, gt_sg = paired
        nodes = gt_sg.nodes_in()
        assert nodes == list(range(gt_sg.N))

    def test_get_node_attributes(self, paired):
        _, gt_sg = paired
        gt_sg.load_vec_on_nodes(np.ones(gt_sg.N), "ones")
        attrs = gt_sg.get_node_attributes("ones")
        assert len(attrs) == gt_sg.N
        assert all(v == 1.0 for v in attrs.values())

    def test_get_edge_data(self, paired):
        _, gt_sg = paired
        # Get a known edge
        edges = gt_sg.get_edges_list()
        u, v = edges[0]
        sign = gt_sg.get_edge_data(u, v, thedata="sign")
        assert sign in (-1, 1)

    def test_get_edge_color(self, paired):
        _, gt_sg = paired
        colors = gt_sg.get_edge_color(pec="green", nec="red")
        assert len(colors) == gt_sg.num_edges
        assert all(c in ("green", "red") for c in colors)

    def test_get_graph_neighbors(self, paired):
        _, gt_sg = paired
        nn = gt_sg.get_graph_neighbors(0)
        assert isinstance(nn, list)
        assert all(isinstance(n, int) for n in nn)

    def test_matrix_for_accessors(self, paired):
        _, gt_sg = paired
        A = gt_sg.get_adjacency_matrix_for()
        D = gt_sg.get_degree_matrix_for()
        L = gt_sg.get_laplacian_matrix_for()
        Ds = gt_sg.get_signed_degree_matrix_for()
        Ls = gt_sg.get_signed_laplacian_matrix_for()
        for M in (A, D, L, Ds, Ls):
            assert M.shape == (gt_sg.N, gt_sg.N)

    def test_gcl_property(self, paired):
        _, gt_sg = paired
        gcl = gt_sg.gcl
        assert gcl is not None
        assert gt_sg.graph_clustering_utility is gcl


# ==========================================================================
# Phase 7: Representation stubs
# ==========================================================================

class TestPhase7Stubs:

    def test_upd_methods_dont_crash(self, paired):
        _, gt_sg = paired
        gt_sg.upd_graph_matrices()
        gt_sg.upd_edge_sets()
        gt_sg.upd_Degree()
        gt_sg.upd_GraphRepr_All()
        gt_sg.upd_NodeMap()
        gt_sg.upd_EdgeMap()
        gt_sg.upd_ReprMaps()

    def test_zip_methods(self, paired):
        _, gt_sg = paired
        m = gt_sg.zip_reprNodes([10, 20, 30])
        assert m == {0: 10, 1: 20, 2: 30}


# ==========================================================================
# End-to-end: Ising dynamics with GT eigenvector init
# ==========================================================================

class TestEndToEnd:

    def test_ising_ground_state_init_gt(self):
        """The original blocker: IsingDynamics with ic='ground_state_0' on GT."""
        from lrgsglib.graphs.gt.Lattice2DGT import Lattice2DGT
        from lrgsglib.statsys.IsingDynamics import IsingDynamics

        lat = Lattice2DGT(side1=6, geo="sqr", pflip=0.1, seed=42)
        lat.flip_random_fract_edges()
        ising = IsingDynamics(
            sg=lat, T=1.0, steps=50, ic="ground_state_0", runlang="py"
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False)
        assert len(ising.magn) > 0
