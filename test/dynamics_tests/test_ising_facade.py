"""Phase-2 strangler facade: legacy IsingDynamics vocabulary -> new classes.

``ising_from_legacy`` must translate every legacy parameter family to the
plan-§10 endpoint (scheme class + axes + observables), emit a
DeprecationWarning, and hard-error on what the new classes cannot represent
(save letters K/V/H, non-weighted spectral proposals, C-subprocess codes —
the C backend stays on the untouched legacy class, decision D-B5).
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.IsingDynamics import (
    IsingCEM,
    IsingDynamics,
    IsingMetropolis,
    IsingParallelTempering,
    IsingSimulatedAnnealing,
    ising_from_legacy,
)
from lrgsglib.statsys.IsingDynamics.defaults import ISING_TFCA_QUENCH_T

SIDE = 8
N = SIDE * SIDE
SEED = 42


def make_lattice(tmp_path, pflip=0.0):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=pflip,
        engine="nx",
        path_data=tmp_path,
        seed=7,
    )


def from_legacy(*args, **kwargs):
    """Every translation must announce itself."""
    with pytest.warns(DeprecationWarning, match="legacy IsingDynamics"):
        return ising_from_legacy(*args, **kwargs)


# =========================================================================
# Metropolis family (default + time-control aliases + order mapping)
# =========================================================================
def test_default_maps_to_metropolis_and_runs(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    model = from_legacy(sg, T=1.5, eqSTEP=7, seed=SEED, savedisk=False)
    assert isinstance(model, IsingMetropolis)
    assert model.T == 1.5 and model.steps == 7
    assert model.runlang == "py"
    assert model.order == "random"  # legacy 'asynchronous'
    model.run()
    assert len(model.ene) == 7


def test_steps_alias_precedence(tmp_path):
    sg = make_lattice(tmp_path)
    model = from_legacy(
        sg, T=1.0, steps=5, nstepsIsing=9, eqSTEP=13, savedisk=False
    )
    assert model.steps == 5
    model = from_legacy(sg, T=1.0, nstepsIsing=9, eqSTEP=13, savedisk=False)
    assert model.steps == 9


def test_upd_mode_mapping(tmp_path):
    sg = make_lattice(tmp_path)
    model = from_legacy(sg, T=1.0, upd_mode="sequential", savedisk=False)
    assert model.order == "typewriter"
    with pytest.raises(ValueError, match="upd_mode"):
        ising_from_legacy(sg, T=1.0, upd_mode="bogus", savedisk=False)


def test_cluster_move_aliases(tmp_path):
    sg = make_lattice(tmp_path)
    model = from_legacy(sg, T=1.0, runlang="wolff", savedisk=False)
    assert isinstance(model, IsingMetropolis) and model.move == "wolff"
    model = from_legacy(sg, T=1.0, runlang="sw", savedisk=False)
    assert model.move == "sw"


def test_pb_backend_passthrough(tmp_path):
    sg = make_lattice(tmp_path)
    model = from_legacy(sg, T=1.0, runlang="pb_met", savedisk=False)
    assert isinstance(model, IsingMetropolis) and model.runlang == "pb"


# =========================================================================
# Save letters -> observables
# =========================================================================
def test_save_letters_map_to_observables(tmp_path):
    sg = make_lattice(tmp_path)
    model = from_legacy(sg, T=1.0, save="ES", savedisk=False)
    assert model._selected_obs == ("energy", "magn", "snapshots")
    model = from_legacy(sg, T=1.0, save="S", savedisk=False)
    assert model._selected_obs == ("snapshots",)


@pytest.mark.parametrize("letters", ["K", "EV", "H"])
def test_unrepresentable_save_letters_raise(tmp_path, letters):
    sg = make_lattice(tmp_path)
    with pytest.raises(NotImplementedError, match="save letter"):
        ising_from_legacy(sg, T=1.0, save=letters, savedisk=False)


# =========================================================================
# SA / PT / topo / CEM families
# =========================================================================
def test_sa_family_mapping(tmp_path):
    sg = make_lattice(tmp_path)
    model = from_legacy(
        sg,
        sa_enabled=True,
        T_init=3.0,
        cooling_schedule="linear",
        T_final=0.5,
        n_temperatures=4,
        steps_per_T=3,
        savedisk=False,
    )
    assert isinstance(model, IsingSimulatedAnnealing)
    assert model.schedule == "linear" and model.steps == 12
    assert model.T_schedule[0] == pytest.approx(3.0)
    assert model.T_schedule[-1] == pytest.approx(0.5)


def test_pt_family_mapping(tmp_path):
    sg = make_lattice(tmp_path)
    model = from_legacy(
        sg,
        pt_enabled=True,
        n_replicas=3,
        T_min=1.0,
        T_max=2.0,
        T_ladder_type="linear",
        n_exchanges=6,
        steps_per_exchange=2,
        savedisk=False,
    )
    assert isinstance(model, IsingParallelTempering)
    assert model.n_replicas == 3 and model.steps == 6
    assert model.T_ladder[0] == pytest.approx(1.0)
    assert model.T_ladder[-1] == pytest.approx(2.0)


def test_topo_met_mapping(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    model = from_legacy(
        sg,
        T=0.5,
        runlang="topo_met",
        topo_n_modes=4,
        topo_polish=False,
        eqSTEP=3,
        savedisk=False,
    )
    assert isinstance(model, IsingMetropolis) and model.move == "spectral"
    assert model.spectral_n_modes == 4 and not model.spectral_polish
    with pytest.raises(NotImplementedError, match="weighted"):
        ising_from_legacy(sg, runlang="topo_met", topo_proposal_mode="uniform")


def test_topo_fca_reproduces_legacy_quench(tmp_path):
    """TFCA = spectral field + constant near-zero-T schedule for
    n_temperatures * steps_per_T sweeps (the legacy sudden quench)."""
    sg = make_lattice(tmp_path, pflip=0.1)
    model = from_legacy(
        sg,
        runlang="topo_fca",
        n_temperatures=4,
        steps_per_T=3,
        topo_n_modes=4,
        seed=SEED,
        savedisk=False,
    )
    assert isinstance(model, IsingSimulatedAnnealing)
    assert model.field_mode == "spectral"
    assert model.steps == 12
    assert np.all(model.T_schedule == ISING_TFCA_QUENCH_T)
    model.run()
    assert np.any(np.asarray(model.field) != 0.0)  # spectral field built
    assert len(model.ene) == 12


def test_cem_mapping(tmp_path):
    sg = make_lattice(tmp_path)
    model = from_legacy(
        sg,
        runlang="cem",
        cem_iter=2,
        cem_restarts=1,
        cem_pop_size=8,
        topo_n_modes=4,
        cem_greedy=False,
        savedisk=False,
    )
    assert isinstance(model, IsingCEM)
    assert model.n_iter == 2 and model.restarts == 1
    assert model.pop_size == 8 and model.n_modes == 4
    assert not model.greedy


def test_unknown_algorithm_token_raises(tmp_path):
    """A typo'd runlang must never silently fall back to Metropolis."""
    sg = make_lattice(tmp_path)
    with pytest.raises(ValueError, match="Unrecognized legacy algorithm"):
        ising_from_legacy(sg, T=1.0, runlang="py_bogus", savedisk=False)


# =========================================================================
# The strangler boundary: C codes stay on the legacy class
# =========================================================================
def test_c_codes_are_not_translated(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(NotImplementedError, match="legacy IsingDynamics"):
        ising_from_legacy(sg, T=1.0, runlang="C0E")
    # ... and the legacy class itself still takes them, untouched.
    legacy = IsingDynamics(sg, T=1.0, runlang="C0E", eqSTEP=3)
    assert legacy.runlang == "C0E"
