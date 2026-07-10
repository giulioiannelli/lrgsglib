"""Cross-model parity harness (decision D-B6 of the statsys unification).

One reusable acceptance surface that pins, for EVERY dynamics model, the
run -> record -> output contract:

1. **Lazy filesystem** — constructing the graph + model creates NO dirs;
   a run that persists nothing creates nothing.
2. **Record** — after a run, each observable the case declares is populated
   (RAM view non-empty, correct length where declared).
3. **Output** — with ``savedisk``-style persistence on, the declared files
   exist under ``model.dynpath`` and the lazy disk-backed views reload.
4. **Reproducibility / backend agreement** — the same seed on the same
   backend reproduces the trajectory exactly; across backends, deterministic
   models must agree numerically and stochastic native kernels are compared
   *statistically* (never byte-equality — they use different RNG streams).

Phase-0 scaffold: the contract + seeded-reproducibility paths are live (Voter
and ContactProcess are the reference consumers); the cross-backend statistical
comparison gains consumers as native backends come up under the CMake fold
(D-B7). Models are added one at a time as they modernize (D-B8/D-B9).

All output goes to pytest ``tmp_path`` — NEVER the real ``data/`` tree.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from lrgsglib.graphs import Lattice2D

# ---------------------------------------------------------------------------
# Harness-wide knobs (single home for the test tunables; the library-side
# tunables live in lrgsglib.config / per-model defaults.py)
# ---------------------------------------------------------------------------
PARITY_SIDE = 8  # lattice side -> N = 64, big enough to be non-trivial
PARITY_PFLIP = 0.1  # disorder fraction exercising the signed couplings
PARITY_STEPS = 5  # sweeps: enough to record, short enough for CI
PARITY_SEED = 12345  # the ONE seed every reproducibility check uses
PARITY_GEO = "sqr"
PARITY_ENGINE = "nx"


def make_parity_graph(tmp_path: Path, **kw: Any):
    """Standard substrate every parity case runs on (writes under *tmp_path*)."""
    params: dict[str, Any] = dict(
        side1=PARITY_SIDE,
        geo=PARITY_GEO,
        pflip=PARITY_PFLIP,
        engine=PARITY_ENGINE,
        path_data=tmp_path,
        seed=PARITY_SEED,
    )
    params.update(kw)
    return Lattice2D(**params)


def dirs_under(root: Path) -> list[str]:
    return sorted(
        str(p.relative_to(root)) for p in Path(root).rglob("*") if p.is_dir()
    )


@dataclass
class ParityCase:
    """One model's entry in the parity suite.

    ``make_model`` receives the graph and must return a ready-to-run model
    (backend, save flags and seed already set). ``expected_observables`` are
    names in ``model.observables`` that must be populated after the run;
    ``expected_lengths`` optionally pins their record counts.
    """

    name: str
    make_model: Callable[[Any], Any]
    run_kwargs: dict[str, Any] = field(default_factory=dict)
    expected_observables: tuple[str, ...] = ()
    expected_lengths: dict[str, int] = field(default_factory=dict)
    #: attribute holding the trajectory/state used for reproducibility checks
    state_attr: str = "s"


# ---------------------------------------------------------------------------
# Contract checks
# ---------------------------------------------------------------------------
def run_record_output_contract(case: ParityCase, tmp_path: Path) -> Any:
    """Assert the full run -> record -> output contract for *case*.

    Returns the ran model so a caller can pile on model-specific asserts.
    """
    sg = make_parity_graph(tmp_path)
    model = case.make_model(sg)

    # 1. lazy filesystem: nothing on disk after construction
    assert (
        dirs_under(tmp_path) == []
    ), f"[{case.name}] construction created dirs: {dirs_under(tmp_path)}"

    model.run(**case.run_kwargs)

    # 2. record: declared observables are populated
    for obs_name in case.expected_observables:
        assert (
            model.observables is not None
        ), f"[{case.name}] model has no ObservableSet"
        data = model.observables[obs_name].data
        n_rec = len(data)
        assert n_rec > 0, f"[{case.name}] observable {obs_name!r} is empty"
        want = case.expected_lengths.get(obs_name)
        if want is not None:
            assert n_rec == want, (
                f"[{case.name}] observable {obs_name!r}: {n_rec} records, "
                f"expected {want}"
            )

    # 3. output: persisted observables exist under dynpath and reload
    if getattr(model, "savedisk", False):
        for obs_name in case.expected_observables:
            obs = model.observables[obs_name]
            if obs.path is None:
                continue  # enabled but RAM-only by model policy
            path = Path(obs.path)
            assert path.exists(), (
                f"[{case.name}] observable {obs_name!r} declares "
                f"path {path} but the file is missing"
            )
            assert path.is_relative_to(Path(model.dynpath)), (
                f"[{case.name}] observable {obs_name!r} landed at {path}, "
                f"outside dynpath {model.dynpath}"
            )
            assert len(obs.data) > 0, (
                f"[{case.name}] observable {obs_name!r}: lazy reload from "
                f"{path} came back empty"
            )
    return model


def assert_seeded_reproducibility(case: ParityCase, tmp_path: Path) -> None:
    """Same seed, same backend => identical final state (bit-for-bit)."""
    finals = []
    for rep in ("a", "b"):
        rep_dir = Path(tmp_path) / rep
        rep_dir.mkdir()
        sg = make_parity_graph(rep_dir)
        model = case.make_model(sg)
        model.run(**case.run_kwargs)
        finals.append(np.asarray(getattr(model, case.state_attr)).copy())
    assert np.array_equal(
        finals[0], finals[1]
    ), f"[{case.name}] identical seeds diverged on the same backend"


def assert_backend_agreement(
    case: ParityCase,
    make_variant: Callable[[Any, str], Any],
    backends: tuple[str, ...],
    tmp_path: Path,
    *,
    statistic: Callable[[Any], float],
    n_replicas: int,
    atol: float,
) -> None:
    """Statistical cross-backend agreement for stochastic kernels.

    Runs *n_replicas* per backend (distinct seeds), compares the mean of
    *statistic* across backends within ``max(atol, 4 * SE_diff)`` where
    ``SE_diff`` is the standard error of the difference of the two means
    — the replica-to-replica scatter is REAL physics (each replica is a
    distinct disorder realization, and short trajectories sit in
    different metastable basins), so a fixed *atol* alone flakes on the
    tail of that distribution. Deterministic models can pass
    ``n_replicas=1`` with a tight *atol* (zero scatter -> pure *atol*).
    Native backends use their own RNG streams, so byte-equality across
    backends is NEVER asserted (D-B6).
    """
    means: dict[str, float] = {}
    sems: dict[str, float] = {}
    for backend in backends:
        vals = []
        for rep in range(n_replicas):
            rep_dir = Path(tmp_path) / f"{backend}_{rep}"
            rep_dir.mkdir(parents=True)
            sg = make_parity_graph(rep_dir, seed=PARITY_SEED + rep)
            model = make_variant(sg, backend)
            model.run(**case.run_kwargs)
            vals.append(statistic(model))
        means[backend] = float(np.mean(vals))
        sems[backend] = float(np.std(vals) / np.sqrt(n_replicas))
    reference = means[backends[0]]
    ref_sem = sems[backends[0]]
    for backend, mean in means.items():
        tol = max(atol, 4.0 * float(np.hypot(sems[backend], ref_sem)))
        assert abs(mean - reference) <= tol, (
            f"[{case.name}] backend {backend!r} mean statistic {mean} "
            f"deviates from {backends[0]!r} ({reference}) beyond "
            f"tol={tol} (atol={atol}); all means: {means}"
        )
