"""Contact process dynamics for signed graphs.

This module provides two flavours of contact-process dynamics:

* :class:`ContactProcessSIR` implements the infection/recovery process
  parameterised by an infection rate ``mu``. Use this class with the pure
  Python backend (``runlang="py"``) or the ``C0`` runlang to mirror the
  ``ContactSimulator0`` C kernel.
* :class:`ContactProcessEI` implements excitation-inhibition dynamics driven
  by ``gamma`` and activation choice, mapping to the ``ContactSimulator1*``
  kernels (``runlang`` values ``C1``, ``C1a``, ``C1b``, ``C1c``). This path
  encapsulates the degree rescaling expected by the C implementations.

Examples
--------
Run the infection/recovery process in Python:

>>> cp = ContactProcessSIR(signed_graph, mu=0.2, runlang="py")
>>> cp.init_contact_dynamics()
>>> cp.run(steps=10, tqdm_on=False)

Use the excitation-inhibition C backend (requires compiled C cores):

>>> cp = ContactProcessEI(signed_graph, gamma=1.5, runlang="C1c")
>>> cp.init_contact_dynamics()
>>> cp.run(verbose=False)
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any, Iterable, Literal, cast

import numpy as np
import tqdm

from .BinDynSys import BinDynSys
from ..config.const import LOG, LRGSG_CCORE_BIN, LRGSG_LOG
from ..nx_patches import SignedGraph
from ..utils.basic.strings import join_non_empty
from ..utils.tools.chronometer import time_function_accumulate


class ContactProcessBase(BinDynSys):
    """Shared utilities for contact-process dynamics.

    Parameters
    ----------
    sg : SignedGraph
        Target graph for the simulation.
    save_density : bool, optional
        Whether to record the active-state density at each time step.
    state_type : {"binary", "bipolar"}, optional
        State encoding passed through to :class:`BinDynSys`.
    **kwargs : Any
        Additional arguments forwarded to :class:`BinDynSys`.
    """

    dyn_UVclass = "contact_process"
    density: list[float] = []
    s_t: list[np.ndarray] = []
    _allowed_c_keys: tuple[str, ...] = ()

    def __init__(
        self,
        sg: SignedGraph,
        *,
        save_density: bool = False,
        state_type: Literal["binary", "bipolar"] | None = None,
        **kwargs: Any,
    ) -> None:
        dynpath = getattr(sg, "path_cntct", None)
        if dynpath is None:
            dynpath = getattr(sg, "path_lrgsg", getattr(sg, "path_data", Path.cwd()))

        super().__init__(
            sg,
            dynpath=dynpath,
            state_type=state_type or "binary",
            **kwargs,
        )
        self.save_density = save_density
        self.reset_observables()
        self.sini: np.ndarray | None = None
        self.stderr_path: Path | None = None
        self.stderr_fopen = None
        self.cprogram: list[str | Path] = []
        self.out_id: str = self.out_suffix

    # ------------------------------------------------------------------
    # Utility helpers
    # ------------------------------------------------------------------
    def reset_observables(self) -> None:
        """Reset cached observables collected during a run."""

        self.density = []
        self.s_t = []

    def _iter_neighbour_data(self, node: int) -> Iterable[tuple[int, float]]:
        graph = self.sg.gr[self.sg.on_g]
        if hasattr(graph, "is_directed") and graph.is_directed():
            iterator = graph.predecessors(node)
            for neighbour in iterator:
                data = graph[neighbour][node]
                weight = data.get("weight", 1.0)
                yield neighbour, weight
        else:
            iterator = graph.neighbors(node)
            for neighbour in iterator:
                data = graph[node][neighbour]
                weight = data.get("weight", 1.0)
                yield neighbour, weight

    def init_contact_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialise the state and export data if required."""

        self.reset_observables()
        if custom is not None:
            self.ic = "custom"
        self.init_s(custom)
        self.s = self.s.astype(np.int8, copy=False)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_s_init()
            if self.rndStr and not exName:
                exName = self.run_id
            exName_arg = exName if exName else self.run_id
            # allow empty exName; SignedGraph._export_edgel_bin will add
            # the required trailing underscore when needed
            self.sg._export_edgel_bin(exName=exName_arg)
        self.sini = self.s.copy()

    def check_attribute(self) -> None:
        try:
            getattr(self, "CbaseName")
        except AttributeError:
            self.init_contact_dynamics()

    def initialize_run_parameters(self, steps: int | None = None) -> None:
        if steps is not None:
            self.simtime = int(steps)

    # ------------------------------------------------------------------
    # Python dynamics
    # ------------------------------------------------------------------
    def contact_sampling(self, tqdm_on: bool) -> None:
        dsNstep = self.dsNstep()
        nodes = np.arange(self.N)
        iterator = tqdm.tqdm(range(self.simtime)) if tqdm_on else range(self.simtime)
        for _ in iterator:
            if self.save_density:
                self.density.append(float(np.mean(self.s)))
            if self.savedyn:
                self.s_t.append(self.s.copy())
            np.random.shuffle(nodes)
            dsNstep(nodes)

    def run_py(self, tqdm_on: bool = False) -> None:
        self.contact_sampling(tqdm_on)

    # ------------------------------------------------------------------
    # C backend integration
    # ------------------------------------------------------------------
    def _c_program_key(self) -> str:
        if not self.runlang or not self.runlang.upper().startswith("C"):
            raise ValueError("C backend requested but runlang is not a C variant.")
        suffix = self.runlang[1:]
        key = "C0" if suffix == "" else f"C{suffix.upper()}"
        if self._allowed_c_keys and key not in self._allowed_c_keys:
            raise ValueError(f"runlang '{self.runlang}' not supported for {self.__class__.__name__}.")
        return key

    @staticmethod
    def _validate_activation(activation: str) -> Literal["tanh", "relu"]:
        normalized = activation.lower()
        if normalized not in {"tanh", "relu"}:
            raise ValueError("activation must be either 'tanh' or 'relu'.")
        return cast(Literal["tanh", "relu"], normalized)

    def _dynamics_out_label(self) -> str:
        gamma = getattr(self, "gamma", None)
        if gamma is not None:
            return f"gamma={float(gamma):.12g}"
        mu = getattr(self, "mu", None)
        if mu is not None:
            return f"mu={float(mu):.12g}"
        return ""

    def _build_out_id(self) -> str:
        """Build identifier used by C backends for output files."""

        return join_non_empty("_", self._dynamics_out_label(), self.out_suffix, self.run_id)

    def _build_c_arglist_base(self) -> list[str]:
        """Build base argument list common to all ContactSimulator variants."""

        try:
            datdir = self.sg.path_sgdata.relative_to(Path.cwd())
        except ValueError:
            datdir = self.sg.path_sgdata
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
        self.out_id = self._build_out_id()
        return [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            str(datdir),
            syshape,
            self.run_id,
            self.out_id,
        ]

    def _build_c_arglist(self) -> list[str]:
        raise NotImplementedError("Subclasses must provide C argument lists.")

    def build_cprogram_command(self) -> None:
        c_key = self._c_program_key()
        if c_key == "C0":
            program_suffix = "0"
        elif c_key == "C1":
            program_suffix = "1"
        else:
            program_suffix = c_key[1:].lower()
        self.CbaseName = f"ContactSimulator{program_suffix}"
        arglist = self._build_c_arglist()
        self.cprogram = [LRGSG_CCORE_BIN / self.CbaseName] + arglist

    def setup_stderr_logging(self) -> None:
        fname = join_non_empty(
            '_',
            f"err{self.CbaseName}",
            f"{self.N}",
            self.run_id,
            self.out_suffix,
        ) + LOG
        self.stderr_path = LRGSG_LOG / fname
        self.stderr_path.parent.mkdir(parents=True, exist_ok=True)
        self.stderr_fopen = open(self.stderr_path, 'w')

    def run_cprogram(self, verbose: bool = False) -> None:
        if not self.cprogram:
            raise RuntimeError("C program command has not been initialised.")
        binary_path = Path(self.cprogram[0])
        if not binary_path.exists():
            if self.stderr_fopen and not self.stderr_fopen.closed:
                self.stderr_fopen.close()
            self.stderr_fopen = None
            raise FileNotFoundError(
                f"C backend executable '{binary_path}' was not found. "
                "Build the C components (e.g. via `make c-make`) before running the C backend."
            )
        result = subprocess.run(
            self.cprogram,
            stderr=self.stderr_fopen,
            stdout=subprocess.PIPE,
            check=False,
        )
        if self.stderr_fopen and not self.stderr_fopen.closed:
            self.stderr_fopen.close()
        self.stderr_fopen = None
        state = np.frombuffer(result.stdout, dtype=np.int8)
        if state.size != self.N:
            raise RuntimeError("C backend returned an invalid state configuration.")
        self.s = state.copy()

    def _remove_sfout(self) -> None:
        try:
            self.sfout.unlink()
        except FileNotFoundError:
            pass

    def _remove_stderr(self) -> None:
        if not self.stderr_path:
            return
        try:
            self.stderr_path.unlink()
        except FileNotFoundError:
            pass

    def remove_run_c_files(self, remove_stderr: bool = True) -> None:
        self._remove_sfout()
        if remove_stderr:
            self._remove_stderr()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    @time_function_accumulate(auto_log=False)
    def run(
        self,
        tqdm_on: bool = True,
        steps: int | None = None,
        verbose: bool = False,
        clean_export: bool = True,
    ) -> None:
        self.check_attribute()
        self.initialize_run_parameters(steps)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.run_cprogram(verbose)
            if clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()
        else:
            self.contact_sampling(tqdm_on)


class ContactProcessSIR(ContactProcessBase):
    """Infection-rate contact process (SIR-style) driven by ``mu``.

    Use this class for the standard infection/recovery dynamics. The Python
    backend (``runlang="py"``) mirrors the logic in :meth:`ds1step`, while the
    ``C0`` runlang targets the ``ContactSimulator0`` executable.
    """

    _allowed_c_keys = ("C0",)

    def __init__(
        self,
        sg: SignedGraph,
        mu: float = 1.0,
        **kwargs: Any,
    ) -> None:
        super().__init__(sg, **kwargs)
        self.mu = float(mu)

    # ------------------------------------------------------------------
    # Python dynamics
    # ------------------------------------------------------------------
    def ds1step(self, node: int) -> None:
        neighbours = list(self._iter_neighbour_data(node))
        if self.s[node]:
            rate = self.mu
            for neighbour, weight in neighbours:
                if weight < 0.0 and self.s[neighbour]:
                    rate += -weight
            prob = 1.0 - np.exp(-rate)
            if prob > 0.0 and np.random.random() < prob:
                self.s[node] = np.int8(0)
        else:
            rate = 0.0
            for neighbour, weight in neighbours:
                if weight > 0.0 and self.s[neighbour]:
                    rate += weight
            prob = 1.0 - np.exp(-rate)
            if prob > 0.0 and np.random.random() < prob:
                self.s[node] = np.int8(1)

    # ------------------------------------------------------------------
    # C backend integration
    # ------------------------------------------------------------------
    def _build_c_arglist(self) -> list[str]:
        base_args = self._build_c_arglist_base()
        return [
            base_args[0],  # N
            base_args[1],  # p
            f"{self.mu:.12g}",
            f"{self.simtime}",
        ] + base_args[2:]


class ContactProcessEI(ContactProcessBase):
    """Excitation-inhibition contact process targeting ``C1`` kernels.

    Parameters
    ----------
    gamma : float
        Excitation strength (rescaled internally by the average degree to match
        the original ``ContactSimulator1*`` interfaces).
    activation : {"tanh", "relu"}, optional
        Non-linearity used by the C1 kernels; ignored by other backends.
    num_log_samples : int, optional
        Number of log samples used by the ``C1c`` variant.

    Notes
    -----
    This class is intended for the C backends (``runlang`` starting with
    ``C1``). Python dynamics are not provided for the excitation-inhibition
    path.
    """

    _allowed_c_keys = ("C1", "C1A", "C1B", "C1C")

    def __init__(
        self,
        sg: SignedGraph,
        *,
        gamma: float,
        activation: Literal["tanh", "relu"] = "tanh",
        num_log_samples: int = 1000,
        runlang: str | None = None,
        **kwargs: Any,
    ) -> None:
        if runlang is not None:
            kwargs.setdefault("runlang", runlang)
        else:
            kwargs.setdefault("runlang", "C1")
        super().__init__(sg, **kwargs)
        self.gamma = float(gamma)
        k = sg.Ne / sg.N
        self.gamma_eff = self.gamma / k
        self.activation = self._validate_activation(activation)
        self.num_log_samples = int(num_log_samples)

    # ------------------------------------------------------------------
    # Python dynamics
    # ------------------------------------------------------------------
    def ds1step(self, node: int) -> None:  # pragma: no cover - not used
        raise NotImplementedError(
            "ContactProcessEI provides only C backends; Python dynamics are not implemented."
        )

    # ------------------------------------------------------------------
    # C backend integration
    # ------------------------------------------------------------------
    def _build_c_arglist(self) -> list[str]:
        base_args = self._build_c_arglist_base()
        key = self._c_program_key()
        args = [
            base_args[0],  # N
            base_args[1],  # p
            f"{self.gamma_eff:.12g}",
            f"{self.simtime}",
        ] + base_args[2:] + [
            self.activation,
        ]

        if key == "C1A":
            nSampleLog = getattr(self, "nSampleLog", 100)
            args.append(f"{nSampleLog}")
        elif key == "C1C":
            args.append(f"{self.num_log_samples}")

        return args

    def run(
        self,
        tqdm_on: bool = True,
        steps: int | None = None,
        verbose: bool = False,
        clean_export: bool = True,
    ) -> None:
        if not self.runlang.upper().startswith("C1"):
            raise NotImplementedError("ContactProcessEI supports only C1* backends.")
        super().run(tqdm_on=tqdm_on, steps=steps, verbose=verbose, clean_export=clean_export)


# Backwards compatibility alias
ContactProcess = ContactProcessSIR
