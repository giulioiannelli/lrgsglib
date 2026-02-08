"""Base class for continuous-state dynamical systems on signed graphs.

Provides ODE integration infrastructure (Euler, RK4, adaptive RK45) for
models where each node carries a continuous scalar state (float64).

Hierarchy
---------
::

    DynSys
    +-- ContDynSys (this module)
        +-- KuramotoModel
        +-- ReactionDiffusionModel
        +-- CoupledODEModel
"""

from __future__ import annotations

from abc import abstractmethod
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from ..DynSys import DynSys

if TYPE_CHECKING:
    from pathlib import Path

    from ...graphs.nx import SignedGraphNX as SignedGraph


class ContDynSys(DynSys):
    """Base class for continuous-scalar dynamics on signed graphs.

    Each node carries a ``float64`` state.  Time evolution is governed by
    an ODE ``dx/dt = f(x, t)`` integrated with the selected scheme.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph substrate.
    dt : float
        Integration time step.
    integrator : {"euler", "rk4", "rk45"}
        ODE integration method.
    atol : float
        Absolute tolerance for adaptive RK45.
    rtol : float
        Relative tolerance for adaptive RK45.
    **kw
        Forwarded to :class:`DynSys`.
    """

    _state_dtype = np.float64

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        dt: float = 0.01,
        integrator: str = "rk4",
        atol: float = 1e-6,
        rtol: float = 1e-3,
        **kw: Any,
    ):
        super().__init__(sg, **kw)
        self.dt = dt
        self.integrator = integrator
        self.atol = atol
        self.rtol = rtol
        self.s: NDArray[np.float64] = np.zeros(self.N, dtype=np.float64)
        self._adj: NDArray | None = None

    # ------------------------------------------------------------------
    # State shape (abstract property from DynSys)
    # ------------------------------------------------------------------
    @property
    def _state_shape(self) -> tuple[int, ...]:
        return (self.N,)

    # ------------------------------------------------------------------
    # Adjacency / coupling helpers
    # ------------------------------------------------------------------
    def _get_adj_matrix(self) -> NDArray:
        """Return the (cached) signed adjacency matrix as dense float64."""
        if self._adj is None:
            import networkx as nx

            G = self.sg.gr[self.sg.on_g]
            self._adj = nx.to_numpy_array(G, weight="weight", dtype=np.float64)
        return self._adj

    # ------------------------------------------------------------------
    # Abstract ODE interface
    # ------------------------------------------------------------------
    @abstractmethod
    def rhs(self, s: NDArray, t: float) -> NDArray:
        """Compute dx/dt = f(x, t).

        Parameters
        ----------
        s : NDArray, shape (N,)
            Current state vector.
        t : float
            Current time.

        Returns
        -------
        NDArray, shape (N,)
            Time derivatives.
        """
        ...

    def _apply_constraints(self, s: NDArray) -> NDArray:
        """Post-step constraint (e.g. angle wrapping).  Override in subclass."""
        return s

    # ------------------------------------------------------------------
    # Integrators
    # ------------------------------------------------------------------
    def _euler_step(self, s: NDArray, t: float) -> NDArray:
        return s + self.dt * self.rhs(s, t)

    def _rk4_step(self, s: NDArray, t: float) -> NDArray:
        dt = self.dt
        k1 = self.rhs(s, t)
        k2 = self.rhs(s + 0.5 * dt * k1, t + 0.5 * dt)
        k3 = self.rhs(s + 0.5 * dt * k2, t + 0.5 * dt)
        k4 = self.rhs(s + dt * k3, t + dt)
        return s + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    def _rk45_step(self, s: NDArray, t: float) -> tuple[NDArray, float]:
        """Dormand-Prince adaptive step.  Returns (new_s, next_dt)."""
        dt = self.dt
        k1 = self.rhs(s, t)
        k2 = self.rhs(s + dt * (1.0 / 5.0) * k1, t + dt / 5.0)
        k3 = self.rhs(s + dt * (3.0 / 40.0 * k1 + 9.0 / 40.0 * k2), t + 3.0 * dt / 10.0)
        k4 = self.rhs(s + dt * (44.0 / 45.0 * k1 - 56.0 / 15.0 * k2 + 32.0 / 9.0 * k3), t + 4.0 * dt / 5.0)
        k5 = self.rhs(s + dt * (19372.0 / 6561.0 * k1 - 25360.0 / 2187.0 * k2 + 64448.0 / 6561.0 * k3 - 212.0 / 729.0 * k4), t + 8.0 * dt / 9.0)
        k6 = self.rhs(s + dt * (9017.0 / 3168.0 * k1 - 355.0 / 33.0 * k2 + 46732.0 / 5247.0 * k3 + 49.0 / 176.0 * k4 - 5103.0 / 18656.0 * k5), t + dt)

        # 5th-order solution
        s_new = s + dt * (35.0 / 384.0 * k1 + 500.0 / 1113.0 * k3 + 125.0 / 192.0 * k4 - 2187.0 / 6784.0 * k5 + 11.0 / 84.0 * k6)

        # 4th-order solution for error estimate
        k7 = self.rhs(s_new, t + dt)
        s_4th = s + dt * (5179.0 / 57600.0 * k1 + 7571.0 / 16695.0 * k3 + 393.0 / 640.0 * k4 - 92097.0 / 339200.0 * k5 + 187.0 / 2100.0 * k6 + 1.0 / 40.0 * k7)

        err = np.max(np.abs(s_new - s_4th) / (self.atol + self.rtol * np.abs(s_new)))
        if err == 0:
            next_dt = 2.0 * dt
        else:
            next_dt = dt * min(2.0, max(0.2, 0.9 * err ** (-0.2)))
        return s_new, next_dt

    def _integrate_step(self, s: NDArray, t: float) -> NDArray:
        """Dispatch to the selected integrator."""
        match self.integrator:
            case "euler":
                return self._euler_step(s, t)
            case "rk4":
                return self._rk4_step(s, t)
            case "rk45":
                s_new, next_dt = self._rk45_step(s, t)
                self.dt = next_dt
                return s_new
            case _:
                raise ValueError(f"Unknown integrator: {self.integrator!r}")

    # ------------------------------------------------------------------
    # State initialisation
    # ------------------------------------------------------------------
    def init_state(self, custom: Any = None) -> None:
        """Initialise the continuous state.

        Modes
        -----
        uniform / random:  U(0, 2*pi) — useful for phase oscillators.
        zero:              All zeros.
        custom:            User-supplied array.
        """
        match self.ic:
            case "uniform" | "random" | "rand":
                self.s = np.random.uniform(0, 2 * np.pi, size=self.N)
            case "zero" | "zeros":
                self.s = np.zeros(self.N, dtype=np.float64)
            case "ones":
                self.s = np.ones(self.N, dtype=np.float64)
            case "custom":
                if custom is None:
                    raise ValueError(
                        "A custom state array must be provided "
                        "for 'custom' initial condition."
                    )
                self.s = np.asarray(custom, dtype=np.float64).copy()
                if self.s.shape != (self.N,):
                    raise ValueError(
                        f"Custom state shape {self.s.shape} "
                        f"does not match (N={self.N},)."
                    )
            case _:
                raise ValueError(f"Unknown initial condition: {self.ic!r}")

    # ------------------------------------------------------------------
    # State export
    # ------------------------------------------------------------------
    def export_state(self) -> None:
        """Write the float64 state to a binary file for C backend."""
        out_suffix = self.run_id or ''
        fname = self.sg.get_p_fname('s', out_suffix=out_suffix)
        self.sfout = self.dynpath / fname
        self.s.astype(np.float64).tofile(open(self.sfout, 'wb'))
        self.s_0 = self.s.copy()

    # ------------------------------------------------------------------
    # Dynamics drivers
    # ------------------------------------------------------------------
    def run_py(self, verbose: bool = False, **kw: Any) -> None:
        """Integrate using the selected Python integrator."""
        s = self.s.copy()
        for step in range(self.steps):
            s = self._integrate_step(s, step * self.dt)
            s = self._apply_constraints(s)
            if self.savedyn:
                self.s_t.append(s.copy())
        self.s = s

    def ds1step(self, nd: int) -> None:
        """Not applicable for ODE-based systems (raises error)."""
        raise NotImplementedError(
            "ContDynSys uses whole-state ODE integration, not single-node updates."
        )

    def run(self, verbose: bool = False, **kw: Any) -> None:
        """Run the simulation using the configured backend."""
        if self.runlang.startswith("C"):
            self.run_c(verbose=verbose, **kw)
        else:
            self.run_py(verbose=verbose, **kw)

    def run_c(self, **kw: Any) -> None:
        """Run C backend (requires CBackendMixin in subclass)."""
        raise NotImplementedError(
            "C backend not configured.  Subclass must also inherit CBackendMixin."
        )
