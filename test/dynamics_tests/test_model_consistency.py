"""Repo-wide dynamics-model structural consistency (Phase 5).

Enforces the modern dynamics pattern across every statsys model:
  1. the package has a ``defaults.py`` (no magic numbers in the model) and a
     ``_solvers.py`` that registers its backends;
  2. ``run()`` dispatches through the shared solver registry, so every
     registered ``(model, backend)`` resolves via ``get_solver``.

Models not yet migrated are tracked explicitly in ``KNOWN_OFF_PATTERN`` (so a
*new* off-pattern model fails this test, while the known tail is xfail-visible).
"""
from pathlib import Path

import pytest

from lrgsglib.statsys._solver_engine import (
    get_solver,
    list_solver_models,
    list_solvers_for_model,
)

# Base classes / shared infrastructure -- not dynamics models.
_BASE_OR_INFRA = {
    "DynSys", "BinDynSys", "ContDynSys", "VecDynSys", "_ccore", "__pycache__",
}
# Models whose run() does not yet dispatch via the registry (tracked tail; see
# .agents/todo/2026-06-25_registry-tail-contactprocess-signedrw.md).
# ContactProcess already has defaults.py + ObservableSet but no _solvers.py;
# SignedRW is fully off-pattern (no defaults.py / _solvers.py).
KNOWN_OFF_PATTERN = {"ContactProcess", "SignedRW"}


def _statsys_dir() -> Path:
    import lrgsglib.statsys as s
    return Path(s.__file__).parent


def _model_packages() -> list[str]:
    d = _statsys_dir()
    return sorted(
        p.name for p in d.iterdir()
        if p.is_dir() and p.name not in _BASE_OR_INFRA
        and (p / "__init__.py").exists()
    )


@pytest.mark.code
@pytest.mark.parametrize("model", list_solver_models())
def test_registered_model_resolves(model):
    """Every registered model resolves a Solver for each of its backends."""
    backends = list_solvers_for_model(model)
    assert backends, f"{model} registered with no backends"
    for b in backends:
        assert get_solver(model, b) is not None


@pytest.mark.code
@pytest.mark.parametrize("pkg", _model_packages())
def test_model_package_on_pattern(pkg):
    """Each model package ships defaults.py + _solvers.py (or is tracked)."""
    d = _statsys_dir() / pkg
    on_pattern = (d / "defaults.py").exists() and (d / "_solvers.py").exists()
    if pkg in KNOWN_OFF_PATTERN:
        if on_pattern:
            pytest.fail(
                f"{pkg} is now on-pattern -- remove it from KNOWN_OFF_PATTERN."
            )
        pytest.xfail(f"{pkg} not yet migrated to the solver registry (tracked).")
    assert on_pattern, (
        f"{pkg} must ship defaults.py + _solvers.py and register via "
        f"register_solver (or be added to KNOWN_OFF_PATTERN with a tracking TODO)."
    )


@pytest.mark.code
def test_solver_registry_nonempty():
    """The registry is populated at import (regression guard)."""
    models = list_solver_models()
    assert len(models) >= 8
    assert "voter" in models
