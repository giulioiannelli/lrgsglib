"""Solver-backend vocabulary for dynamics models.

This is the dynamics-side analogue of ``graphs/_engine.py``'s ``GraphEngine``:
``SolverBackend`` enumerates the *families* a model's ``run()`` can dispatch to
(pure-Python, pybind11 native, vectorized NumPy/CuPy, C subprocess), and
``Solver`` is the thin uniform interface each registered backend implements. The
actual registry lives in ``_solver_engine.py``.

A ``Solver`` does **not** own simulation state — it operates on the model passed
to its methods, so a single stateless instance per (model, backend) is enough.
"""

from __future__ import annotations

from collections.abc import Callable
from enum import Enum
from typing import TYPE_CHECKING, Protocol, Union, runtime_checkable

import numpy as np

if TYPE_CHECKING:
    from .DynSys.DynSys import DynSys

__all__ = ["SolverBackend", "Solver", "FrameHook", "coerce_backend"]


class SolverBackend(Enum):
    """Solver families a dynamics model can dispatch to.

    Attributes
    ----------
    PY : str
        Pure-Python reference loop (``'py'``). Always available; the semantic
        reference every native backend is validated against.
    PB : str
        In-process pybind11 native kernel (``'pb'``). Seed-reproducible,
        graph-engine agnostic (consumes CSR arrays).
    NP : str
        Vectorized NumPy backend (``'np'``).
    CU : str
        Vectorized CuPy / GPU backend (``'cu'``).
    C : str
        C subprocess backend (``'c'``); transports state via files.
    """

    PY = "py"
    PB = "pb"
    NP = "np"
    CU = "cu"
    C = "c"


# (sweep_index, state_copy) -> None. Consumed by the live-view GUI (wired by the
# frame-producer hook in a later phase); declared here so the Solver interface
# is stable across phases.
FrameHook = Callable[[int, np.ndarray], None]


@runtime_checkable
class Solver(Protocol):
    """Uniform backend interface a model's ``run()`` dispatches through.

    Implementations are stateless wrappers: ``supports`` validates the model's
    current configuration for this backend (raising on an unsupported combo),
    and ``execute`` invokes the corresponding kernel on the model. An optional
    ``is_available()`` (duck-typed, not required by the Protocol) reports whether
    the backend's runtime (compiled ``.so``, C binary, CuPy) is present.
    """

    backend: SolverBackend

    def supports(self, model: "DynSys") -> None:
        """Raise (``NotImplementedError``/``ValueError``) if the model's current
        configuration cannot run on this backend; return ``None`` otherwise."""
        ...

    def execute(self, model: "DynSys", **kwargs: object) -> None:
        """Run the dynamics on ``model`` using this backend's kernel."""
        ...


def coerce_backend(value: Union[str, SolverBackend]) -> SolverBackend:
    """Normalize a backend selector (``'py'`` / ``SolverBackend.PY``) to the enum."""
    if isinstance(value, SolverBackend):
        return value
    return SolverBackend(str(value).lower())
