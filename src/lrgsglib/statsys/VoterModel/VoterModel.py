"""Voter-model dynamics on signed graphs."""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any

import numpy as np
import tqdm

from .._c_backend import CBackendMixin
from ..BinDynSys import BinDynSys
from ...utils.tools.chronometer import time_function_accumulate


class VoterModel(CBackendMixin, BinDynSys):
    """Binary voter dynamics with optional C backend.

    In voter dynamics, each node copies the state of a randomly chosen
    neighbor, modified by the edge sign (weight). This creates consensus
    or polarization patterns depending on the graph structure.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph to run dynamics on.
    steps : int, optional
        Number of Monte Carlo sweeps.
    simref : float, optional
        Size-normalized time (steps = simref * N).
    eqSTEP : int, optional
        Legacy alias for steps.
    save_magnetization : bool, optional
        If True, record magnetization at each step.
    upd_mode : str, optional
        Update mode ('asynchronous' or 'synchronous').
    freq : int, optional
        Recording frequency.
    nSampleLog : int, optional
        Number of log-spaced samples (for C backend).
    **kwargs
        Additional arguments passed to BinDynSys.

    Examples
    --------
    >>> voter = VoterModel(lattice, steps=100, runlang='py')
    >>> voter.init_voter_dynamics()
    >>> voter.run(tqdm_on=False)
    """

    dyn_UVclass = "voter_model"

    # CBackendMixin configuration
    _c_bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    _c_program_name_template = "VoterSimulator{}"
    _allowed_c_keys = ("C0", "C1")

    # Class-level observable defaults (will be overwritten per instance)
    magn: list[float] = []
    s_t: list[np.ndarray] = []

    def __init__(
        self,
        sg,
        *,
        steps: int | None = None,
        simref: float | None = None,
        eqSTEP: int | None = None,
        save_magnetization: bool = False,
        upd_mode: str = "asynchronous",
        freq: int = 10,
        nSampleLog: int = 100,
        **kwargs: Any,
    ) -> None:
        dynpath = getattr(sg, 'path_voter', None)
        resolved_steps = steps if steps is not None else eqSTEP
        super().__init__(sg, dynpath=dynpath, steps=resolved_steps, simref=simref, **kwargs)
        self.save_magnetization = save_magnetization
        self.upd_mode = upd_mode
        self.freq = freq
        self.nSampleLog = nSampleLog
        self.reset_observables()
        self.sini: np.ndarray | None = None
        self.out_id: str = self.out_suffix
        self.magn_path: Path | None = None

    @property
    def eqSTEP(self) -> int:
        """Compatibility alias for the configured number of sweeps."""
        return self.steps

    @eqSTEP.setter
    def eqSTEP(self, value: int) -> None:
        self._set_time_controls(steps=value)

    # ------------------------------------------------------------------
    # Utility helpers
    # ------------------------------------------------------------------
    def reset_observables(self) -> None:
        """Reset cached observables collected during a run."""
        self.magn = []
        self.s_t = []

    def init_voter_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialise the spin configuration and export data if required."""
        self._check_c_backend_or_fallback()
        self.reset_observables()
        self.init_s(custom)
        self.s = self.s.astype(np.int8, copy=False)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_s_init()
            if self.rndStr and not exName:
                exName = self.run_id
            # Export uses build_p_fname which handles underscore via join_non_empty
            self.sg._export_edgel_bin(exName=self.run_id)
            self.sg.export_adj_bin(exName=self.run_id)
        self.sini = self.s.copy()

    def check_attribute(self) -> None:
        """Initialize dynamics if not already done."""
        # Check sini (set at end of init_voter_dynamics) rather than CbaseName
        # because CbaseName has a class-level default from CBackendMixin
        if self.sini is None:
            self.init_voter_dynamics()

    def initialize_run_parameters(
        self,
        steps: int | None = None,
        simref: float | None = None,
        eqSTEP: int | None = None,
    ) -> None:
        chosen_steps = steps if steps is not None else eqSTEP
        self._set_time_controls(steps=chosen_steps, simref=simref)

    # ------------------------------------------------------------------
    # Python dynamics
    # ------------------------------------------------------------------
    def ds1step(self, nd: int) -> None:
        neighbors = self.sg.get_neighbors_with_weights(nd)
        if not neighbors:
            return
        choice_idx = np.random.randint(len(neighbors))
        neighbor, weight = neighbors[choice_idx]
        edge_sign = -1 if weight < 0 else 1
        self.s[nd] = np.int8(edge_sign * self.s[neighbor])

    def voter_sampling(self, tqdm_on: bool) -> None:
        dsNstep = self.dsNstep()
        nodes = list(range(self.N))
        iterator = tqdm.tqdm(range(self.steps)) if tqdm_on else range(self.steps)
        for _ in iterator:
            if self.save_magnetization:
                self.magn.append(float(np.sum(self.s)) / float(self.N))
            if self.savedyn:
                self.s_t.append(self.s.copy())
            dsNstep(random.sample(nodes, len(nodes)))

    def run_py(self, tqdm_on: bool = False) -> None:
        self.voter_sampling(tqdm_on)

    # ------------------------------------------------------------------
    # C backend integration (via CBackendMixin)
    # ------------------------------------------------------------------
    def _c_program_suffix(self) -> str:
        """Extract suffix for VoterSimulator (e.g., 'C1' -> '1')."""
        # VoterSimulator uses single digit: VoterSimulator0, VoterSimulator1
        return self.runlang[-1]

    def _build_c_arglist(self) -> list[str]:
        """Build argument list for VoterSimulator."""
        try:
            datdir = self.sg.path_sgdata.relative_to(Path.cwd())
        except ValueError:
            datdir = self.sg.path_sgdata
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
        self.out_id = self.out_suffix
        self.magn_path = self.dynpath / self.sg.get_p_fname('m', self.out_id)

        # Base arguments common to all simulators
        arglist = [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            f"{self.steps}",
            str(datdir),
            syshape,
            self._c_suffix_arg(self.run_id),
            self._c_suffix_arg(self.out_id),
        ]

        # Add extra arguments for VoterSimulator1 and above
        if self.runlang[-1] != "0":
            arglist.append(f"{self.nSampleLog}")

        return arglist

    def run_cprogram(self, verbose: bool = False) -> None:
        """Execute C backend and read magnetization output."""
        # Call parent implementation for subprocess execution
        super().run_cprogram(verbose)
        # Read magnetization file if it exists
        if self.magn_path and self.magn_path.exists():
            self.magn = np.fromfile(self.magn_path, dtype=np.float64).tolist()

    def _get_cleanup_paths(self) -> list[Path | None]:
        """Return paths to clean up after C run."""
        return [
            getattr(self, 'sfout', None),
            self.magn_path,
        ]

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    @time_function_accumulate(auto_log=False)
    def run(
        self,
        tqdm_on: bool = True,
        steps: int | None = None,
        simref: float | None = None,
        eqSTEP: int | None = None,
        verbose: bool = False,
        clean_export: bool = True,
    ) -> None:
        """Run voter model dynamics.

        Parameters
        ----------
        tqdm_on : bool
            Show progress bar.
        steps : int, optional
            Override number of steps.
        simref : float, optional
            Size-normalized time.
        eqSTEP : int, optional
            Legacy alias for steps.
        verbose : bool
            Verbose output.
        clean_export : bool
            Remove exported files after run.
        """
        self.check_attribute()
        self.initialize_run_parameters(steps=steps, simref=simref, eqSTEP=eqSTEP)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.run_cprogram(verbose)
            if clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()
        else:
            self.voter_sampling(tqdm_on)
