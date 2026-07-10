"""Smoke tests for all dynamics models.

Each test instantiates the model on a small lattice, initializes,
runs a few steps, and verifies the state array has the right shape.
"""

import numpy as np
import pytest

from lrgsglib.graphs.nx import Lattice2DNX


@pytest.fixture
def small_lattice(tmp_path):
    """8x8 square lattice with mild frustration.

    path_data=tmp_path: savedisk-enabled models (the flow trio) persist
    per-run outputs, which must never land in the repo's data/.
    """
    lat = Lattice2DNX(
        side1=8, geo="sqr", pflip=0.1, seed=42, path_data=tmp_path
    )
    lat.flip_random_fract_edges()
    return lat


# =========================================================================
# ContDynSys models
# =========================================================================


class TestKuramotoModel:
    def test_instantiation(self, small_lattice):
        from lrgsglib.statsys.KuramotoModel import KuramotoModel

        km = KuramotoModel(
            sg=small_lattice, coupling=1.0, dt=0.01, steps=10, seed=42
        )
        km.init_kuramoto_dynamics()
        assert km.s is not None
        assert len(km.s) == small_lattice.N

    def test_run(self, small_lattice):
        from lrgsglib.statsys.KuramotoModel import KuramotoModel

        km = KuramotoModel(
            sg=small_lattice, coupling=1.0, dt=0.01, steps=10, seed=42
        )
        km.init_kuramoto_dynamics()
        km.run(verbose=False)
        # State should have evolved
        assert km.s is not None


class TestReactionDiffusionModel:
    def test_instantiation(self, small_lattice):
        from lrgsglib.statsys.ReactionDiffusionModel import (
            ReactionDiffusionModel,
        )

        rd = ReactionDiffusionModel(
            sg=small_lattice, D=0.1, dt=0.01, steps=10, seed=42
        )
        rd.init_rd_dynamics()
        assert rd.s is not None
        assert len(rd.s) == small_lattice.N

    def test_run(self, small_lattice):
        from lrgsglib.statsys.ReactionDiffusionModel import (
            ReactionDiffusionModel,
        )

        rd = ReactionDiffusionModel(
            sg=small_lattice, D=0.1, dt=0.01, steps=10, seed=42
        )
        rd.init_rd_dynamics()
        rd.run(verbose=False)
        assert rd.s is not None


class TestCoupledODEModel:
    def test_instantiation(self, small_lattice):
        from lrgsglib.statsys.CoupledODEModel import CoupledODEModel

        ode = CoupledODEModel(
            sg=small_lattice, coupling_type="linear", dt=0.01, steps=10, seed=42
        )
        ode.init_ode_dynamics()
        assert ode.s is not None
        assert len(ode.s) == small_lattice.N

    def test_run(self, small_lattice):
        from lrgsglib.statsys.CoupledODEModel import CoupledODEModel

        ode = CoupledODEModel(
            sg=small_lattice, coupling_type="linear", dt=0.01, steps=10, seed=42
        )
        ode.init_ode_dynamics()
        ode.run(verbose=False)
        assert ode.s is not None


# =========================================================================
# VecDynSys models
# =========================================================================


class TestPottsModel:
    def test_instantiation(self, small_lattice):
        from lrgsglib.statsys.PottsModel import PottsModel

        pm = PottsModel(sg=small_lattice, q=3, T=1.0, steps=10, seed=42)
        pm.init_potts_dynamics()
        assert pm.s is not None
        assert len(pm.s) == small_lattice.N

    def test_run(self, small_lattice):
        from lrgsglib.statsys.PottsModel import PottsModel

        pm = PottsModel(sg=small_lattice, q=3, T=1.0, steps=10, seed=42)
        pm.init_potts_dynamics()
        pm.run(verbose=False)
        assert pm.s is not None

    def test_q_states(self, small_lattice):
        from lrgsglib.statsys.PottsModel import PottsModel

        pm = PottsModel(sg=small_lattice, q=4, T=1.0, steps=5, seed=42)
        pm.init_potts_dynamics()
        # States should be in range [0, q)
        assert np.all(pm.s >= 0)
        assert np.all(pm.s < 4)


class TestXYModel:
    def test_instantiation(self, small_lattice):
        from lrgsglib.statsys.XYModel import XYModel

        xy = XYModel(sg=small_lattice, T=1.0, steps=10, seed=42)
        xy.init_xy_dynamics()
        assert xy.s is not None
        assert len(xy.s) == small_lattice.N

    def test_run(self, small_lattice):
        from lrgsglib.statsys.XYModel import XYModel

        xy = XYModel(sg=small_lattice, T=1.0, steps=10, seed=42)
        xy.init_xy_dynamics()
        xy.run(verbose=False)
        assert xy.s is not None

    def test_angle_range(self, small_lattice):
        from lrgsglib.statsys.XYModel import XYModel

        xy = XYModel(sg=small_lattice, T=1.0, steps=5, seed=42)
        xy.init_xy_dynamics()
        # Angles should be in [0, 2*pi)
        assert np.all(xy.s >= 0)
        assert np.all(xy.s < 2 * np.pi + 1e-10)


class TestHeisenbergModel:
    def test_instantiation(self, small_lattice):
        from lrgsglib.statsys.HeisenbergModel import HeisenbergModel

        hm = HeisenbergModel(sg=small_lattice, T=1.0, steps=10, seed=42)
        hm.init_heisenberg_dynamics()
        assert hm.s is not None

    def test_run(self, small_lattice):
        from lrgsglib.statsys.HeisenbergModel import HeisenbergModel

        hm = HeisenbergModel(sg=small_lattice, T=1.0, steps=10, seed=42)
        hm.init_heisenberg_dynamics()
        hm.run(verbose=False)
        assert hm.s is not None


class TestMultiSpeciesModel:
    def test_instantiation(self, small_lattice):
        from lrgsglib.statsys.MultiSpeciesModel import MultiSpeciesModel

        ms = MultiSpeciesModel(
            sg=small_lattice,
            species=2,
            q_per_species=3,
            T=1.0,
            steps=10,
            seed=42,
        )
        ms.init_multispecies_dynamics()
        assert ms.s is not None

    def test_run(self, small_lattice):
        from lrgsglib.statsys.MultiSpeciesModel import MultiSpeciesModel

        ms = MultiSpeciesModel(
            sg=small_lattice,
            species=2,
            q_per_species=3,
            T=1.0,
            steps=10,
            seed=42,
        )
        ms.init_multispecies_dynamics()
        ms.run(verbose=False)
        assert ms.s is not None
