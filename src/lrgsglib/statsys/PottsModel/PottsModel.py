"""q-state Potts model on signed graphs.

Each node *i* carries a state sigma_i in {0, 1, ..., q-1}. The
Hamiltonian is::

    H = -sum_{(i,j)} A_ij * delta(sigma_i, sigma_j)

where delta is the Kronecker delta and A_ij is the signed adjacency.
Updates use Metropolis-Hastings.

For q=2 and bipolar mapping {0,1} -> {-1,+1}, the model is equivalent
to the Ising model.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from .._c_backend import CBackendMixin
from .._csr import build_graph_csr
from .._solver import SolverBackend
from .._solver_engine import get_solver
from ..VecDynSys import VecDynSys
from .defaults import POTTS_SOLVER_NAME
from ...utils.tools.chronometer import time_function_accumulate

if TYPE_CHECKING:
    from ...graphs.nx import SignedGraphNX as SignedGraph


class PottsModel(CBackendMixin, VecDynSys):
    """q-state Potts model with Metropolis dynamics.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph substrate.
    q : int
        Number of Potts states (default 3).
    T : float
        Temperature.
    steps : int
        Number of Monte Carlo sweeps.
    **kw
        Forwarded to :class:`VecDynSys`.

    Examples
    --------
    >>> potts = PottsModel(lattice, q=3, T=1.0, steps=1000)
    >>> potts.init_potts_dynamics()
    >>> potts.run()
    >>> print(potts.energy())
    """

    dyn_UVclass = "potts_model"

    _c_bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    _c_program_name_template = "PottsSimulator{}"
    _allowed_c_keys: tuple[str, ...] = ("C0",)

    ene: list[float] = []
    magn: list[float] = []

    def __init__(
        self,
        sg: "SignedGraph",
        q: int = 3,
        T: float = 1.0,
        *,
        steps: int | None = None,
        simref: float | None = None,
        save_observables: bool = False,
        **kw: Any,
    ) -> None:
        dynpath = getattr(sg, 'path_data', None)
        if dynpath is not None:
            dynpath = Path(dynpath) / "potts"
        super().__init__(
            sg,
            q=q,
            discrete=True,
            T=T,
            steps=steps,
            simref=simref,
            dynpath=dynpath,
            **kw,
        )
        self.save_observables = save_observables
        self.ene = []
        self.magn = []
        self.sini: NDArray | None = None

    @property
    def _state_shape(self) -> tuple[int, ...]:
        return (self.N,)  # scalar per node (int32)

    # ------------------------------------------------------------------
    # Energy and observables
    # ------------------------------------------------------------------
    def energy(self, s: NDArray | None = None) -> float:
        """Compute H = -sum A_ij delta(sigma_i, sigma_j)."""
        if s is None:
            s = self.s
        A = self._get_adj_matrix()
        E = 0.0
        for i in range(self.N):
            for j in range(i + 1, self.N):
                if A[i, j] != 0:
                    E -= A[i, j] * float(s[i] == s[j])
        return E

    def order_parameter(self, s: NDArray | None = None) -> float:
        """Potts order parameter: (q * max_frac - 1) / (q - 1)."""
        if s is None:
            s = self.s
        counts = np.bincount(s.astype(int), minlength=self.q)
        max_frac = counts.max() / self.N
        if self.q == 1:
            return 1.0
        return (self.q * max_frac - 1.0) / (self.q - 1.0)

    def per_state_density(self, s: NDArray | None = None) -> NDArray:
        """Fraction of nodes in each state."""
        if s is None:
            s = self.s
        counts = np.bincount(s.astype(int), minlength=self.q)
        return counts / self.N

    # ------------------------------------------------------------------
    # Metropolis update
    # ------------------------------------------------------------------
    def ds1step(self, nd: int) -> None:
        """Single-node Metropolis update for Potts."""
        A = self._get_adj_matrix()
        current = self.s[nd]
        # Propose new state != current
        proposal = np.random.randint(0, self.q - 1)
        if proposal >= current:
            proposal += 1

        # Energy change
        dE = 0.0
        for j in range(self.N):
            if A[nd, j] != 0:
                dE += A[nd, j] * (float(current == self.s[j]) - float(proposal == self.s[j]))

        # Accept/reject
        if dE <= 0 or (self.T > 0 and np.random.uniform() < np.exp(-dE / self.T)):
            self.s[nd] = np.int32(proposal)

    # ------------------------------------------------------------------
    # Dynamics
    # ------------------------------------------------------------------
    def init_potts_dynamics(self, custom: Any = None, exName: str = "") -> None:
        self._check_c_backend_or_fallback()
        self.ene = []
        self.magn = []
        self.init_state(custom)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_state()
            if self.rndStr and not exName:
                exName = self.run_id
            self.sg._export_edgel_bin(exName=exName or self.run_id)
        self.sini = self.s.copy()

    def check_attribute(self) -> None:
        if self.sini is None:
            self.init_potts_dynamics()

    def run_py(self, verbose: bool = False, **kw: Any) -> None:
        """Monte Carlo sweeps with Metropolis."""
        for sweep in range(self.steps):
            nodes = np.random.permutation(self.N)
            for nd in nodes:
                self.ds1step(nd)
            if self.save_observables:
                self.ene.append(self.energy())
                self.magn.append(self.order_parameter())
            if self.savedyn:
                self.s_t.append(self.s.copy())

    # ------------------------------------------------------------------
    # Native (pybind11) backend
    # ------------------------------------------------------------------
    def _run_pybind(self) -> None:
        """Run Metropolis via the in-process ``_potts_native`` kernel.

        Marshals the graph into the engine-agnostic CSR triple (works for NX and
        GT) and runs the same ``potts_metropolis_sweep`` kernel the C subprocess
        uses; seed-reproducible. Mirrors ``VoterModel._run_pybind``.
        """
        from .ccore import _potts_native

        ni, nw, nptr = build_graph_csr(self.sg, self.N)
        s0 = np.ascontiguousarray(self.s, dtype=np.int32)
        s, ene, magn = _potts_native.potts_sampling(
            s0, ni, nw, nptr,
            int(self.q), float(self.T),
            int(self.steps), int(self.seed),
            bool(self.save_observables),
        )
        self.s = s
        if self.save_observables:
            self.ene = ene.tolist()
            self.magn = magn.tolist()

    # ------------------------------------------------------------------
    # C backend
    # ------------------------------------------------------------------
    def _build_c_arglist(self) -> list[str]:
        datdir = self._get_datdir_arg()
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
        return [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            f"{self.steps}",
            f"{self.q}",
            f"{self.T:.12g}",
            str(datdir),
            syshape,
            self._c_suffix_arg(self.run_id),
            self._c_suffix_arg(self.out_suffix),
        ]

    def _get_cleanup_paths(self) -> list[Path | None]:
        return [getattr(self, 'sfout', None)]

    @time_function_accumulate(auto_log=False)
    def run(
        self,
        verbose: bool = False,
        clean_export: bool = True,
        **kw: Any,
    ) -> None:
        self.check_attribute()
        # Resolve the solver family from the runlang code (same precedence the
        # old if/elif chain used; the C check is case-sensitive on a leading "C"
        # to avoid catching "cu"/native prefixes).
        if self.runlang.startswith("C"):
            backend = SolverBackend.C
        elif self.runlang.lower().startswith("pb"):
            backend = SolverBackend.PB
        else:
            backend = SolverBackend.PY
        solver = get_solver(POTTS_SOLVER_NAME, backend)
        solver.supports(self)
        try:
            solver.execute(self, verbose=verbose)
        finally:
            if backend is SolverBackend.C and clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()


# Register PottsModel's solver backends (py/pb/C) in the shared solver registry.
# Imported after the class so that importing PottsModel wires up dispatch; no
# circular import (``_solvers`` does not import this class).
from . import _solvers  # noqa: E402,F401
