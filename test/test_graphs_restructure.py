"""
Tests for the unified graphs module restructuring.

This test module verifies that:
1. New import paths work correctly (graphs/nx/, graphs/gt/)
2. Facades at graphs/ root work correctly
3. NX-suffixed class names work
4. Engine selection works correctly
5. All graph types can be created via the unified API
"""

import pytest
import warnings


class TestNewImportPaths:
    """Test that the new import paths work correctly."""

    def test_nx_base_import(self):
        """Test importing SignedGraphNX from new location."""
        from lrgsglib.graphs.nx.base import SignedGraphNX
        assert SignedGraphNX is not None

    def test_nx_lattice_imports(self):
        """Test importing lattice classes from new location."""
        from lrgsglib.graphs.nx.lattice import Lattice2DNX, Lattice3DNX
        assert Lattice2DNX is not None
        assert Lattice3DNX is not None

    def test_nx_random_imports(self):
        """Test importing random graph classes from new location."""
        from lrgsglib.graphs.nx.random import (
            ErdosRenyiNX,
            BarabasiAlbertNX,
            WattsStrogatzNX,
            StochasticBlockModelNX,
        )
        assert ErdosRenyiNX is not None
        assert BarabasiAlbertNX is not None
        assert WattsStrogatzNX is not None
        assert StochasticBlockModelNX is not None

    def test_nx_module_import(self):
        """Test importing from nx module __init__."""
        from lrgsglib.graphs.nx import (
            SignedGraphNX,
            Lattice2DNX,
            Lattice3DNX,
            ErdosRenyiNX,
            BarabasiAlbertNX,
            WattsStrogatzNX,
            StochasticBlockModelNX,
        )
        assert all([
            SignedGraphNX, Lattice2DNX, Lattice3DNX,
            ErdosRenyiNX, BarabasiAlbertNX, WattsStrogatzNX, StochasticBlockModelNX
        ])


class TestFacadeImports:
    """Test that facade imports from graphs/ root work."""

    def test_facade_lattice_imports(self):
        """Test importing lattice facades."""
        from lrgsglib.graphs import Lattice2D, Lattice3D
        assert callable(Lattice2D)
        assert callable(Lattice3D)

    def test_facade_random_imports(self):
        """Test importing random graph facades."""
        from lrgsglib.graphs import (
            ErdosRenyi,
            BarabasiAlbert,
            WattsStrogatz,
            StochasticBlockModel,
        )
        assert callable(ErdosRenyi)
        assert callable(BarabasiAlbert)
        assert callable(WattsStrogatz)
        assert callable(StochasticBlockModel)

    def test_engine_imports(self):
        """Test importing engine management functions."""
        from lrgsglib.graphs import (
            GraphEngine,
            set_default_engine,
            get_default_engine,
            is_engine_available,
            available_engines,
        )
        assert GraphEngine is not None
        assert callable(set_default_engine)
        assert callable(get_default_engine)
        assert callable(is_engine_available)
        assert callable(available_engines)


class TestNXGraphCreation:
    """Test creating graphs with NetworkX backend."""

    def test_create_lattice2d_nx(self):
        """Test creating Lattice2D with NX engine."""
        from lrgsglib.graphs import Lattice2D
        lat = Lattice2D(side1=5, geo='sqr', engine='nx', seed=42)
        assert lat.N == 25  # 5x5 lattice
        assert hasattr(lat, 'flip_random_fract_edges')

    def test_create_erdos_renyi_nx(self):
        """Test creating ErdosRenyi with NX engine."""
        from lrgsglib.graphs import ErdosRenyi
        er = ErdosRenyi(n=50, p=0.2, engine='nx', seed=42)
        assert er.N > 0
        assert hasattr(er, 'flip_random_fract_edges')

    def test_create_barabasi_albert_nx(self):
        """Test creating BarabasiAlbert with NX engine."""
        from lrgsglib.graphs import BarabasiAlbert
        ba = BarabasiAlbert(n=50, m=2, engine='nx', seed=42)
        assert ba.N == 50
        assert hasattr(ba, 'flip_random_fract_edges')

    def test_create_watts_strogatz_nx(self):
        """Test creating WattsStrogatz with NX engine."""
        from lrgsglib.graphs import WattsStrogatz
        ws = WattsStrogatz(n=50, k=4, p=0.3, engine='nx', seed=42)
        assert ws.N == 50
        assert hasattr(ws, 'flip_random_fract_edges')

    def test_create_sbm_nx(self):
        """Test creating StochasticBlockModel with NX engine."""
        from lrgsglib.graphs import StochasticBlockModel
        sizes = [15, 15, 20]
        p_in, p_out = 0.3, 0.05
        p_matrix = [
            [p_in, p_out, p_out],
            [p_out, p_in, p_out],
            [p_out, p_out, p_in],
        ]
        sbm = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, engine='nx', seed=42)
        assert sbm.N > 0
        assert hasattr(sbm, 'flip_random_fract_edges')


class TestNXSuffixedClasses:
    """Test that NX-suffixed classes work correctly."""

    def test_lattice2dnx_class(self):
        """Test Lattice2DNX class directly."""
        from lrgsglib.graphs.nx.lattice import Lattice2DNX
        lat = Lattice2DNX(side1=5, geo='sqr', seed=42)
        assert lat.N == 25
        assert type(lat).__name__ in ('Lattice2D', 'Lattice2DNX')

    def test_erdos_renyi_nx_class(self):
        """Test ErdosRenyiNX class directly."""
        from lrgsglib.graphs.nx.random import ErdosRenyiNX
        er = ErdosRenyiNX(n=50, p=0.2, seed=42)
        assert er.N > 0
        assert type(er).__name__ in ('ErdosRenyi', 'ErdosRenyiNX')


class TestEngineSelection:
    """Test engine selection mechanism."""

    def test_default_engine_is_nx(self):
        """Test that default engine is NetworkX."""
        from lrgsglib.graphs import get_default_engine, GraphEngine
        assert get_default_engine() == GraphEngine.NETWORKX

    def test_nx_engine_available(self):
        """Test that NX engine is always available."""
        from lrgsglib.graphs import is_engine_available
        assert is_engine_available('nx')

    def test_engine_fallback(self):
        """Test that GT falls back to NX if GT not available."""
        from lrgsglib.graphs import Lattice2D, is_engine_available
        # This should work even if GT is not available (fallback to NX)
        lat = Lattice2D(side1=5, geo='sqr', engine='gt')
        assert lat.N == 25


class TestBackwardCompatibility:
    """Test backward compatibility with old import paths."""

    def test_old_nx_patches_lattice2d_import(self):
        """Test that old import path still works."""
        from lrgsglib.nx_patches import Lattice2D
        lat = Lattice2D(side1=5, geo='sqr', seed=42)
        assert lat.N == 25

    def test_old_nx_patches_erdos_renyi_import(self):
        """Test that old import path still works."""
        from lrgsglib.nx_patches import ErdosRenyi
        er = ErdosRenyi(n=50, p=0.2, seed=42)
        assert er.N > 0

    def test_old_nx_patches_signed_graph_import(self):
        """Test that old SignedGraph import still works."""
        from lrgsglib.nx_patches import SignedGraph
        assert SignedGraph is not None


class TestFlipEdges:
    """Test edge flipping functionality works correctly."""

    def test_pflip_creates_edges_to_flip(self):
        """Test that pflip marks edges for flipping."""
        from lrgsglib.graphs import Lattice2D
        lat = Lattice2D(side1=5, geo='sqr', pflip=0.3, engine='nx', seed=42)
        lat.flip_random_fract_edges()
        neg_edges = lat.count_negative_edges()
        assert neg_edges > 0

    def test_zero_pflip_no_negative_edges(self):
        """Test that pflip=0 results in no negative edges."""
        from lrgsglib.graphs import Lattice2D
        lat = Lattice2D(side1=5, geo='sqr', pflip=0.0, engine='nx', seed=42)
        neg_edges = lat.count_negative_edges()
        assert neg_edges == 0


class TestAllNXImports:
    """Test that all NX wrapper classes can be imported."""

    def test_nx_base_classes(self):
        """Test base NX class imports."""
        from lrgsglib.graphs.nx import SignedGraphNX
        assert SignedGraphNX is not None

    def test_nx_lattice_classes(self):
        """Test lattice NX class imports."""
        from lrgsglib.graphs.nx import LatticeNDNX, Lattice2DNX, Lattice3DNX
        assert LatticeNDNX is not None
        assert Lattice2DNX is not None
        assert Lattice3DNX is not None

    def test_nx_random_classes(self):
        """Test random NX class imports."""
        from lrgsglib.graphs.nx import (
            RandomGraphNX, ErdosRenyiNX, BarabasiAlbertNX, WattsStrogatzNX,
            StochasticBlockModelNX, kRegularGraphNX, ConfigurationModelNX,
            RandomGeometricNX, LFRBenchmarkNX, ExtendedBarabasiAlbertNX,
            DualBarabasiAlbertNX, HolmeKimNX,
        )
        assert all([
            RandomGraphNX, ErdosRenyiNX, BarabasiAlbertNX, WattsStrogatzNX,
            StochasticBlockModelNX, kRegularGraphNX, ConfigurationModelNX,
            RandomGeometricNX, LFRBenchmarkNX, ExtendedBarabasiAlbertNX,
            DualBarabasiAlbertNX, HolmeKimNX,
        ])

    def test_nx_complete_classes(self):
        """Test complete NX class imports."""
        from lrgsglib.graphs.nx import CompleteGraphNX, FullyConnectedNX
        assert CompleteGraphNX is not None
        assert FullyConnectedNX is not None

    def test_nx_neural_classes(self):
        """Test neural NX class imports."""
        from lrgsglib.graphs.nx import HofieldNNNX, SCSGeneralizedNNNX
        assert HofieldNNNX is not None
        assert SCSGeneralizedNNNX is not None

    def test_nx_fractal_classes(self):
        """Test fractal NX class imports."""
        from lrgsglib.graphs.nx import FractalGraphNX, DGMgraphNX, SierpinskiGraphNX
        assert FractalGraphNX is not None
        assert DGMgraphNX is not None
        assert SierpinskiGraphNX is not None

    def test_nx_bipartite_classes(self):
        """Test bipartite NX class imports."""
        from lrgsglib.graphs.nx import BipartiteGraphNX, BipartiteFromDegreeSequenceNX
        assert BipartiteGraphNX is not None
        assert BipartiteFromDegreeSequenceNX is not None

    def test_nx_multispectral_classes(self):
        """Test multispectral NX class imports."""
        from lrgsglib.graphs.nx import (
            MultispectralGraphNX, MultiplicativeCascadeGraphNX,
            VicsekGraphNX, HierarchicalModularNetworkNX,
        )
        assert MultispectralGraphNX is not None
        assert MultiplicativeCascadeGraphNX is not None
        assert VicsekGraphNX is not None
        assert HierarchicalModularNetworkNX is not None

    def test_nx_dirac_classes(self):
        """Test dirac NX class imports."""
        from lrgsglib.graphs.nx import DiracLatticeGraphNX, DiracCombGraphNX, DiracBrushGraphNX
        assert DiracLatticeGraphNX is not None
        assert DiracCombGraphNX is not None
        assert DiracBrushGraphNX is not None

    def test_nx_temporal_classes(self):
        """Test temporal NX class imports."""
        from lrgsglib.graphs.nx import TemporalGraphNX, TemporalSignedGraphNX
        assert TemporalGraphNX is not None
        assert TemporalSignedGraphNX is not None


class TestAllGTImports:
    """Test that all GT wrapper classes can be imported."""

    def test_gt_base_classes(self):
        """Test base GT class imports."""
        from lrgsglib.graphs.gt import SignedGraphGT
        assert SignedGraphGT is not None

    def test_gt_lattice_classes(self):
        """Test lattice GT class imports."""
        from lrgsglib.graphs.gt import LatticeNDGT, Lattice2DGT, Lattice3DGT
        assert LatticeNDGT is not None
        assert Lattice2DGT is not None
        assert Lattice3DGT is not None

    def test_gt_random_classes_native(self):
        """Test native GT random class imports."""
        from lrgsglib.graphs.gt import (
            ErdosRenyiGT, BarabasiAlbertGT, WattsStrogatzGT, StochasticBlockModelGT,
        )
        assert ErdosRenyiGT is not None
        assert BarabasiAlbertGT is not None
        assert WattsStrogatzGT is not None
        assert StochasticBlockModelGT is not None

    def test_gt_random_classes_fallback(self):
        """Test fallback GT random class imports (with warning suppression)."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from lrgsglib.graphs.gt import (
                RandomGraphGT, kRegularGraphGT, ConfigurationModelGT,
                RandomGeometricGT, LFRBenchmarkGT, ExtendedBarabasiAlbertGT,
                DualBarabasiAlbertGT, HolmeKimGT,
            )
            assert all([
                RandomGraphGT, kRegularGraphGT, ConfigurationModelGT,
                RandomGeometricGT, LFRBenchmarkGT, ExtendedBarabasiAlbertGT,
                DualBarabasiAlbertGT, HolmeKimGT,
            ])

    def test_gt_complete_classes(self):
        """Test complete GT class imports."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from lrgsglib.graphs.gt import CompleteGraphGT, FullyConnectedGT
            assert CompleteGraphGT is not None
            assert FullyConnectedGT is not None

    def test_gt_neural_classes(self):
        """Test neural GT class imports."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from lrgsglib.graphs.gt import HofieldNNGT, SCSGeneralizedNNGT
            assert HofieldNNGT is not None
            assert SCSGeneralizedNNGT is not None

    def test_gt_fractal_classes(self):
        """Test fractal GT class imports."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from lrgsglib.graphs.gt import FractalGraphGT, DGMgraphGT, SierpinskiGraphGT
            assert FractalGraphGT is not None
            assert DGMgraphGT is not None
            assert SierpinskiGraphGT is not None

    def test_gt_bipartite_classes(self):
        """Test bipartite GT class imports."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from lrgsglib.graphs.gt import BipartiteGraphGT, BipartiteFromDegreeSequenceGT
            assert BipartiteGraphGT is not None
            assert BipartiteFromDegreeSequenceGT is not None

    def test_gt_multispectral_classes(self):
        """Test multispectral GT class imports."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from lrgsglib.graphs.gt import (
                MultispectralGraphGT, MultiplicativeCascadeGraphGT,
                VicsekGraphGT, HierarchicalModularNetworkGT,
            )
            assert MultispectralGraphGT is not None
            assert MultiplicativeCascadeGraphGT is not None
            assert VicsekGraphGT is not None
            assert HierarchicalModularNetworkGT is not None

    def test_gt_dirac_classes(self):
        """Test dirac GT class imports."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from lrgsglib.graphs.gt import DiracLatticeGraphGT, DiracCombGraphGT, DiracBrushGraphGT
            assert DiracLatticeGraphGT is not None
            assert DiracCombGraphGT is not None
            assert DiracBrushGraphGT is not None

    def test_gt_temporal_classes(self):
        """Test temporal GT class imports."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from lrgsglib.graphs.gt import TemporalGraphGT, TemporalSignedGraphGT
            assert TemporalGraphGT is not None
            assert TemporalSignedGraphGT is not None


class TestGTFallbackClasses:
    """Test that GT fallback classes correctly alias NX classes."""

    def test_fallback_aliases_nx(self):
        """Test that fallback GT classes are aliases to NX classes."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from lrgsglib.graphs.gt.complete import CompleteGraphGT
            from lrgsglib.graphs.nx.complete import CompleteGraphNX
            # Fallback should be the same class as NX
            assert CompleteGraphGT is CompleteGraphNX

    def test_native_gt_not_alias(self):
        """Test that native GT classes are NOT aliases to NX classes."""
        from lrgsglib.graphs.gt.random import ErdosRenyiGT
        from lrgsglib.graphs.nx.random import ErdosRenyiNX
        # Native GT should be different from NX
        assert ErdosRenyiGT is not ErdosRenyiNX
