"""Tests for the topology screening module."""

import numpy as np
import pandas as pd
import pytest

from lrgsglib.graphs.nx import Lattice2DNX
from lrgsglib.graphs.nx.SierpinskiNX.SierpinskiNX import SierpinskiNX
from lrgsglib.graphs.nx.kRegularGraphNX.kRegularGraphNX import kRegularGraphNX
from lrgsglib.graphs.nx.Lattice3DNX.Lattice3DNX import Lattice3DNX
from lrgsglib.graphs.nx.BarabasiAlbertNX.BarabasiAlbertNX import BarabasiAlbertNX
from lrgsglib.utils.lrg.percolation import (
    build_cycle_dual,
    default_p_range,
)
from lrgsglib.graphs.nx.funcs.dual import build_planar_dual
from lrgsglib.utils.lrg.screening import (
    compute_face_statistics,
    screen_topology,
    run_screening_batch,
    generate_gog_configs,
    fujii_reference_data,
    distance_above_pareto_front,
    compute_graph_properties,
    screen_topology_general,
    run_general_screening_batch,
    generate_lattice_configs,
    generate_random_configs,
)


class TestComputeFaceStatistics:
    """Tests for ``compute_face_statistics``."""

    def test_square_lattice_uniform(self):
        """Square lattice: all faces size 4, CV=0."""
        fs = compute_face_statistics(
            Lattice2DNX, {"side1": 8, "geo": "sqr", "pbc": False}
        )
        assert fs["is_planar"]
        assert fs["uniformity_score"] == pytest.approx(0.0)
        assert fs["dominant_face_size"] == 4
        assert fs["dominant_face_fraction"] == pytest.approx(1.0)

    def test_hexagonal_lattice_nearly_uniform(self):
        """Hexagonal lattice: mostly size 6 faces, low CV."""
        fs = compute_face_statistics(
            Lattice2DNX, {"side1": 8, "geo": "hex", "pbc": False}
        )
        assert fs["is_planar"]
        # Boundary effects give slight non-uniformity
        assert fs["uniformity_score"] < 0.2

    def test_triangular_lattice_uniform(self):
        """Triangular lattice: all faces size 3, CV=0."""
        fs = compute_face_statistics(
            Lattice2DNX, {"side1": 8, "geo": "tri", "pbc": False}
        )
        assert fs["is_planar"]
        assert fs["uniformity_score"] == pytest.approx(0.0)
        assert fs["dominant_face_size"] == 3

    def test_sierpinski_carpet_nonuniform(self):
        """Carpet has multi-scale faces, high CV."""
        fs = compute_face_statistics(
            SierpinskiNX, {"n": 3, "variant": "carpet"}
        )
        assert fs["is_planar"]
        assert fs["uniformity_score"] > 0.5

    def test_kagome_two_face_sizes(self):
        """Kagome has mixed face sizes (3 and 6)."""
        fs = compute_face_statistics(
            Lattice2DNX, {"side1": 8, "geo": "kgm", "pbc": False}
        )
        assert fs["is_planar"]
        assert fs["uniformity_score"] > 0.0
        assert len(fs["face_size_counts"]) >= 2

    def test_euler_formula(self):
        """V - E + F = 2 for any planar graph."""
        fs = compute_face_statistics(
            Lattice2DNX, {"side1": 10, "geo": "sqr", "pbc": False}
        )
        # n_faces excludes outer face; add 1 for Euler
        V = fs["n_vertices"]
        E = fs["n_edges"]
        F = fs["n_faces"] + 1  # +1 for outer face
        assert V - E + F == 2


class TestScreenTopology:
    """Tests for ``screen_topology``."""

    def test_square_lattice_passes(self):
        """Square lattice should pass all filters."""
        result = screen_topology(
            Lattice2DNX,
            {"side1": 8, "geo": "sqr", "pbc": False},
            p_range=np.array([0.05, 0.15, 0.25]),
            n_realizations=5,
            seed_base=42,
        )
        assert result["status"] == "success"
        assert result["p_X"] is not None
        assert result["p_Z"] is not None
        assert result["kesten_sum"] is not None
        assert result["wall_time_s"] > 0

    def test_filters_nonuniform_faces(self):
        """Carpet should be filtered by face uniformity."""
        result = screen_topology(
            SierpinskiNX,
            {"n": 2, "variant": "carpet"},
            p_range=np.array([0.1, 0.2]),
            n_realizations=3,
            face_uniformity_cutoff=0.3,  # strict cutoff
        )
        assert result["status"] == "face_nonuniform"
        assert result["p_X"] is None
        assert result["face_stats"]["uniformity_score"] > 0.3

    def test_relaxed_cutoff_passes_carpet(self):
        """Carpet passes with relaxed cutoff."""
        result = screen_topology(
            SierpinskiNX,
            {"n": 2, "variant": "carpet"},
            p_range=np.array([0.1, 0.2]),
            n_realizations=3,
            face_uniformity_cutoff=float("inf"),
        )
        assert result["status"] == "success"

    def test_result_keys(self):
        """Result dict has all expected keys."""
        result = screen_topology(
            Lattice2DNX,
            {"side1": 6, "geo": "sqr", "pbc": False},
            p_range=np.array([0.1]),
            n_realizations=3,
        )
        expected_keys = {
            "status", "face_stats", "p_X", "p_Z",
            "kesten_sum", "sweep_data", "wall_time_s",
        }
        assert expected_keys.issubset(result.keys())


class TestRunScreeningBatch:
    """Tests for ``run_screening_batch``."""

    def test_returns_dataframe(self):
        """Batch should return a DataFrame."""
        configs = [
            {
                "name": "sqr_6",
                "factory": Lattice2DNX,
                "params": {"side1": 6, "geo": "sqr", "pbc": False},
            },
        ]
        df = run_screening_batch(
            configs,
            p_range=np.array([0.1, 0.2]),
            n_realizations=3,
            verbose=False,
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        assert "name" in df.columns
        assert "status" in df.columns
        assert "p_X" in df.columns

    def test_handles_mixed_results(self):
        """Mix of uniform and non-uniform should both appear."""
        configs = [
            {
                "name": "sqr",
                "factory": Lattice2DNX,
                "params": {"side1": 6, "geo": "sqr", "pbc": False},
            },
            {
                "name": "carpet",
                "factory": SierpinskiNX,
                "params": {"n": 2, "variant": "carpet"},
            },
        ]
        df = run_screening_batch(
            configs,
            p_range=np.array([0.1, 0.2]),
            n_realizations=3,
            face_uniformity_cutoff=0.3,
            verbose=False,
        )
        assert len(df) == 2
        statuses = set(df["status"])
        assert "success" in statuses
        assert "face_nonuniform" in statuses


class TestGenerateGogConfigs:
    """Tests for ``generate_gog_configs``."""

    def test_default_nonempty(self):
        """Default config generation produces configs."""
        configs = generate_gog_configs()
        assert len(configs) > 0

    def test_expected_count(self):
        """3 base × 3 fiber × 3 side × 2 side × 3 anchor = 162."""
        configs = generate_gog_configs()
        assert len(configs) == 3 * 3 * 3 * 2 * 3

    def test_config_keys(self):
        """Each config has name, factory, params."""
        configs = generate_gog_configs(
            base_geos=["sqr"],
            fiber_geos=["sqr"],
            base_sides=[3],
            fiber_sides=[2],
            anchor_policies=["first"],
        )
        assert len(configs) == 1
        cfg = configs[0]
        assert "name" in cfg
        assert "factory" in cfg
        assert "params" in cfg
        assert "base_graph_type" in cfg["params"]

    def test_custom_parameters(self):
        """Custom geos/sides produce correct count."""
        configs = generate_gog_configs(
            base_geos=["sqr", "hex"],
            fiber_geos=["sqr"],
            base_sides=[4],
            fiber_sides=[2, 3],
            anchor_policies=["first"],
        )
        assert len(configs) == 2 * 1 * 1 * 2 * 1


class TestFujiiReferenceData:
    """Tests for ``fujii_reference_data``."""

    def test_contains_four_lattices(self):
        fujii = fujii_reference_data()
        assert "Square" in fujii
        assert "Kagome" in fujii
        assert "Hexagonal" in fujii
        assert "Tri-hexa" in fujii

    def test_square_self_dual(self):
        fujii = fujii_reference_data()
        sq = fujii["Square"]
        assert sq["p_X"] == pytest.approx(sq["p_Z"])

    def test_kesten_sum(self):
        fujii = fujii_reference_data()
        for d in fujii.values():
            assert d["bpt_primal"] + d["bpt_dual"] == pytest.approx(1.0)


class TestDistanceAbovePareto:
    """Tests for ``distance_above_pareto_front``."""

    def test_square_on_front(self):
        """Square lattice is on the front (distance ≈ 0)."""
        d = distance_above_pareto_front(0.103, 0.103)
        assert abs(d) < 0.01

    def test_above_front_positive(self):
        """Point above the front has positive distance."""
        d = distance_above_pareto_front(0.103, 0.2)
        assert d > 0

    def test_below_front_negative(self):
        """Point below the front has negative distance."""
        d = distance_above_pareto_front(0.103, 0.05)
        assert d < 0


# ---------------------------------------------------------------
# Generalized cycle dual tests
# ---------------------------------------------------------------


class TestBuildCycleDual:
    """Tests for ``build_cycle_dual``."""

    def test_matches_planar_dual_2d(self):
        """MCB dual of 4x4 grid matches planar dual exactly."""
        import networkx as nx

        G = nx.grid_2d_graph(4, 4)
        G = nx.convert_node_labels_to_integers(G)
        for u, v in G.edges():
            G[u][v]["weight"] = 1

        dual_mcb, ci = build_cycle_dual(G)
        dual_planar, _ = build_planar_dual(G, include_outer_face=False)

        assert dual_mcb.number_of_nodes() == dual_planar.number_of_nodes()
        assert dual_mcb.number_of_edges() == dual_planar.number_of_edges()

    def test_3d_all_length4_cycles(self):
        """3D grid MCB should give all length-4 cycles."""
        import networkx as nx

        G = nx.grid_graph(dim=[3, 3, 3])
        G = nx.convert_node_labels_to_integers(G)
        for u, v in G.edges():
            G[u][v]["weight"] = 1

        _, ci = build_cycle_dual(G)
        assert all(s == 4 for s in ci["cycle_sizes"])
        # b1 = E - N + 1
        N, E = G.number_of_nodes(), G.number_of_edges()
        assert ci["n_cycles"] == E - N + 1

    def test_kregular_nonplanar(self):
        """Cycle dual works on non-planar k-regular graph."""
        import networkx as nx

        G = nx.random_regular_graph(6, 50, seed=42)
        for u, v in G.edges():
            G[u][v]["weight"] = 1
        dual, ci = build_cycle_dual(G)
        assert dual.number_of_nodes() == ci["n_cycles"]
        assert dual.number_of_nodes() > 0


# ---------------------------------------------------------------
# General topology screening tests
# ---------------------------------------------------------------


class TestComputeGraphProperties:
    """Tests for ``compute_graph_properties``."""

    def test_kregular_exact_degree(self):
        """k-regular: z_var = 0, z_avg = k."""
        props = compute_graph_properties(
            kRegularGraphNX, {"n": 100, "k": 6}
        )
        assert props["avg_degree"] == pytest.approx(6.0)
        assert props["degree_variance"] == pytest.approx(0.0)

    def test_spectral_gap_positive(self):
        """Connected graph should have positive spectral gap."""
        props = compute_graph_properties(
            Lattice2DNX, {"side1": 8, "geo": "sqr", "pbc": False}
        )
        assert props["spectral_gap"] > 0
        assert props["is_connected"]

    def test_3d_lattice(self):
        """3D lattice properties computed correctly."""
        props = compute_graph_properties(
            Lattice3DNX, {"dim": 4, "geo": "sc", "pbc": False}
        )
        assert props["N"] == 64
        assert props["n_independent_cycles"] > 0
        assert not props["is_planar"]


class TestScreenTopologyGeneral:
    """Tests for ``screen_topology_general``."""

    def test_kregular_succeeds(self):
        """k-regular graph screening should succeed."""
        result = screen_topology_general(
            kRegularGraphNX,
            {"n": 50, "k": 4},
            p_range=np.array([0.05, 0.15, 0.25]),
            n_realizations=3,
            seed_base=42,
        )
        assert result["status"] == "success"
        assert 0 < result["p_X"] < 0.5
        assert 0 < result["p_Z"] < 0.5

    def test_nonplanar_works(self):
        """Non-planar BA graph should screen successfully."""
        result = screen_topology_general(
            BarabasiAlbertNX,
            {"n": 50, "m": 3},
            p_range=np.array([0.1, 0.2, 0.3]),
            n_realizations=3,
        )
        assert result["status"] == "success"

    def test_3d_lattice_produces_thresholds(self):
        """3D SC should produce valid (p_X, p_Z) thresholds."""
        p_range = np.array([0.05, 0.15, 0.25, 0.35])
        r3d = screen_topology_general(
            Lattice3DNX, {"dim": 4, "geo": "sc", "pbc": False},
            p_range, n_realizations=3,
        )
        assert r3d["status"] == "success"
        assert 0 < r3d["p_Z"] < 0.5
        assert 0 < r3d["p_X"] < 0.5

    def test_result_keys(self):
        """Result dict has all expected keys."""
        result = screen_topology_general(
            kRegularGraphNX, {"n": 30, "k": 4},
            p_range=np.array([0.1]),
            n_realizations=2,
        )
        expected = {"status", "p_X", "p_Z", "kesten_sum",
                    "graph_properties", "sweep_data", "wall_time_s"}
        assert expected.issubset(result.keys())


class TestRunGeneralScreeningBatch:
    """Tests for ``run_general_screening_batch``."""

    def test_returns_dataframe(self):
        configs = [
            {"name": "k4", "factory": kRegularGraphNX,
             "params": {"n": 50, "k": 4}},
            {"name": "k8", "factory": kRegularGraphNX,
             "params": {"n": 50, "k": 8}},
        ]
        df = run_general_screening_batch(
            configs, p_range=np.array([0.1, 0.2]),
            n_realizations=3, verbose=False,
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert "p_X" in df.columns
        assert "p_Z" in df.columns


class TestGeneralConfigGenerators:
    """Tests for config generators."""

    def test_lattice_configs(self):
        configs = generate_lattice_configs()
        assert len(configs) > 0
        dims = {c["dimension"] for c in configs}
        assert 2 in dims
        assert 3 in dims

    def test_random_configs(self):
        configs = generate_random_configs(n_values=[100])
        assert len(configs) > 0
        families = {c["family"] for c in configs}
        assert "k-regular" in families
        assert "erdos-renyi" in families
        assert "watts-strogatz" in families
        assert "barabasi-albert" in families
