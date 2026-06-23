"""Phase-1 spine: the SolverBackend registry and VoterModel's dispatch through it.

Mirrors the role ``graphs/_engine.py`` plays for graph engines — these tests pin
the registry contract (capability table, resolution, env override, fail-fast on
unsupported pairs) and that the BackendManager is reachable from a package-level
home.
"""

import os

import pytest

from lrgsglib._backend import Backend, BackendManager
from lrgsglib.statsys import (
    SolverBackend,
    get_solver,
    is_backend_available,
    list_solver_models,
    list_solvers_for_model,
)
from lrgsglib.statsys._solver_engine import get_backend_from_env


@pytest.mark.code
def test_voter_registered_with_all_backends():
    assert "voter" in list_solver_models()
    caps = {b.value for b in list_solvers_for_model("voter")}
    assert caps == {"py", "pb", "np", "cu", "c"}


@pytest.mark.code
@pytest.mark.parametrize("backend", ["py", "pb", "np", "cu", "c"])
def test_get_solver_resolves_and_reports_backend(backend):
    solver = get_solver("voter", backend)
    assert solver.backend is SolverBackend(backend)
    # supports() and execute() are the uniform interface every solver provides.
    assert callable(solver.supports) and callable(solver.execute)
    # availability is a bool (py/np always True; pb/cu/c depend on runtime)
    assert isinstance(is_backend_available("voter", backend), bool)


@pytest.mark.code
def test_py_and_np_backends_are_available():
    assert is_backend_available("voter", "py") is True
    assert is_backend_available("voter", "np") is True


@pytest.mark.code
def test_unknown_model_and_unregistered_backend_raise():
    with pytest.raises(ValueError):
        get_solver("does_not_exist", "py")
    with pytest.raises(ValueError):
        # 'zz' is not a valid SolverBackend at all
        get_solver("voter", "zz")


@pytest.mark.code
def test_env_backend_override(monkeypatch):
    monkeypatch.setenv("LRGSG_SOLVER_BACKEND", "pb")
    assert get_backend_from_env() is SolverBackend.PB
    monkeypatch.setenv("LRGSG_SOLVER_BACKEND", "bogus")
    with pytest.warns(UserWarning):
        assert get_backend_from_env() is None
    monkeypatch.delenv("LRGSG_SOLVER_BACKEND", raising=False)
    assert get_backend_from_env() is None


@pytest.mark.code
def test_backend_manager_reachable_from_package_home():
    # Promotion: the numpy/scipy/cupy manager is importable from lrgsglib._backend
    assert BackendManager.is_available(Backend.NUMPY) is True


@pytest.mark.code
def test_voter_run_dispatches_through_registry(small_lattice_sqr):
    """A py run completes and populates magn via the registry-dispatched path."""
    from lrgsglib.statsys import VoterModel

    sg = small_lattice_sqr
    vm = VoterModel(sg=sg, steps=10, runlang="py", savemagn=True,
                    absorbing_check=False, seed=3)
    vm.init_voter_dynamics()
    vm.run(tqdm_on=False)
    assert len(vm.magn) == 10
