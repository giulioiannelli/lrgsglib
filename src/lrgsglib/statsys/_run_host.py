"""Shared ``run()`` front-door + observable lifecycle for dynamics models.

Every modernized model runs the same way: resolve the ``runlang`` code to a
:class:`~._solver.SolverBackend`, fetch the registered solver for the model's
registry key, validate the configuration (``supports``), make sure the model is
initialised, resolve the time controls, open the output streams, execute, and
finalize the disk artifacts — cleaning the C subprocess's transient export
files afterwards. ``VoterModel`` and ``ContactProcessBase`` hand-copy this
skeleton today; :class:`RunHostMixin` is the single home it is hoisted into.

The mixin owns the *skeleton*; the model owns the *policy*. Everything
model-specific stays behind small hooks with safe defaults:

- ``_begin_outputs`` / ``_persist_observables`` — resolve paths, open/close
  streams (default: no-op, observables stay in RAM);
- ``_record`` — sample the enabled observables once per sweep (default: no-op);
- ``_reset_run_state`` — clear model-local run counters on reset;
- ``check_attribute`` / ``initialize_run_parameters`` — inherited from the
  model (both Voter and CP already define them).

Directory policy: the mixin creates NO directories. Output dirs are lazy and
appear at the shared write funnels only — the base-class state/field exports,
the :mod:`._io` codecs, and ``CBackendMixin._get_datdir_arg`` — so a run that
never writes never touches the filesystem.

Consumers: ``IsingBase`` is the first (Phase 1 of the Ising refactor);
``VoterModel``/``ContactProcessBase`` migrate onto it later under their
existing test nets.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, cast

from ..utils.tools.chronometer import time_function_accumulate
from ._naming import run_dirname, write_cfg_sidecar
from ._solver import SolverBackend
from ._solver_engine import get_solver

if TYPE_CHECKING:
    from .DynSys.DynSys import DynSys

__all__ = ["RunHostMixin"]


class RunHostMixin:
    """Mixin hosting the shared ``run()`` skeleton for dynamics models.

    A consuming model must set :attr:`solver_name` (its key in the shared
    solver registry — ONE per model; scheme/leaf polymorphism rides on the
    instance the solver's ``execute()`` receives) and sit on a ``DynSys``
    base (which provides ``runlang``, ``sg``, ``steps`` and the time
    controls). Models with a C subprocess also mix in ``CBackendMixin``.
    """

    #: Registry key for this model in the shared solver registry
    #: (e.g. ``"voter"``, ``"contact_process"``, ``"ising"``).
    solver_name: ClassVar[str] = ""

    if TYPE_CHECKING:
        # Provided by the DynSys base (and CBackendMixin for the C branch);
        # declared for the type checker only — never assigned here.
        runlang: str
        sg: Any

        def check_attribute(self) -> None: ...
        def initialize_run_parameters(
            self, steps: int | None = None, simref: float | None = None
        ) -> None: ...
        def remove_run_c_files(self) -> None: ...

    # ------------------------------------------------------------------
    # Backend resolution (runlang code -> SolverBackend)
    # ------------------------------------------------------------------
    def _resolve_backend(self) -> SolverBackend:
        """Map the ``runlang`` code to a solver family.

        Same precedence the hand-written chains used: C subprocess > pybind >
        NumPy > CuPy > Python. The C check is case-sensitive on a leading
        ``"C"`` so it never catches ``"cu"``.
        """
        runlang = self.runlang
        if runlang.startswith("C"):
            return SolverBackend.C
        low = runlang.lower()
        if low.startswith("pb"):
            return SolverBackend.PB
        if low.startswith("np"):
            return SolverBackend.NP
        if low.startswith("cu"):
            return SolverBackend.CU
        return SolverBackend.PY

    # ------------------------------------------------------------------
    # Observable lifecycle (hooks with safe defaults)
    # ------------------------------------------------------------------
    def reset_observables(self) -> None:
        """Reset the in-RAM observable buffers and all disk-backing state."""
        observables = getattr(self, "observables", None)
        if observables is not None:
            observables.reset()
        self._reset_run_state()

    def _reset_run_state(self) -> None:
        """Clear model-local run counters (subsample counters, absorbed-at
        markers, ...). Default: nothing to clear."""

    def _begin_outputs(self) -> None:
        """Resolve disk paths and open output streams before a run.

        Default: no-op — observables stay in RAM. Models gate this on their
        ``savedisk``-style flags and per-observable save switches.
        """

    def _record(self) -> None:
        """Sample the enabled observables for the current sweep (in-process
        backends only; the C subprocess writes its own files). Default: no-op.
        """

    def _persist_observables(self) -> None:
        """Finalize disk artifacts and switch attributes to disk-backed views.

        Default: no-op. Always called (via ``finally``) after the solver ran,
        so an interrupted run still closes its streams.
        """

    # ------------------------------------------------------------------
    # Per-run output directory + cfg.json sidecar (Phase C naming)
    # ------------------------------------------------------------------
    def _name_tokens(self) -> list | None:
        """Ordered ``(key, value[, default])`` schema for the run dirname.

        ``None`` (default) keeps the model's flat legacy layout; the
        new-style Bases override this to opt into per-run directories
        (``<dynpath>/<run dirname>/{ene.bin, m.bin, sout.bin, cfg.json}``).
        """
        return None

    def _lang_token(self) -> str:
        """Backend provenance for the ``lang=`` token (resolved by run())."""
        backend = getattr(self, "_active_backend", None)
        if backend is not None:
            return backend.value
        return str(self.runlang).lower()

    def _run_output_dir(self) -> Path | None:
        """The run's output directory (``None`` -> flat legacy layout).

        A non-empty ``out_suffix`` is appended verbatim as the last token,
        so repeat runs a user labels by hand stay distinguishable.
        """
        tokens = self._name_tokens()
        if tokens is None:
            return None
        name = run_dirname(tokens)
        suffix = getattr(self, "out_suffix", "")
        if suffix:
            name = f"{name}_{suffix}"
        return Path(cast(Any, self).dynpath) / name

    def _cfg_graph_block(self) -> dict[str, Any]:
        """Substrate identity: p names the disorder ENSEMBLE, the graph
        seed pins the REALIZATION — both are needed to reproduce."""
        sg = self.sg
        return {
            "class": type(sg).__name__,
            "N": getattr(sg, "N", None),
            "pflip": getattr(sg, "pflip", None),
            "seed": getattr(sg, "seed", None),
            "syshape": getattr(sg, "syshapePth", None),
            "path_name": getattr(sg, "sgpathn", getattr(sg, "_sgpathn", None)),
        }

    def _cfg_model_block(self) -> dict[str, Any]:
        """Every constructor argument (elided defaults included); the
        new-style Bases override and extend."""
        return {"class": type(self).__name__}

    def _cfg_run_block(self) -> dict[str, Any]:
        model = cast(Any, self)
        try:
            from importlib.metadata import version

            _version = version("lrgsglib")
        except Exception:
            _version = "unknown"
        return {
            "steps": getattr(model, "steps", None),
            "observables": sorted(getattr(model, "_selected_obs", ()) or ()),
            "runlang": model.runlang,
            "backend": self._lang_token(),
            "seed": getattr(model, "seed", None),
            "out_suffix": getattr(model, "out_suffix", ""),
            "lrgsglib": _version,
            "date": datetime.now().isoformat(timespec="seconds"),
        }

    def _write_cfg_sidecar(self, rundir: Path) -> None:
        write_cfg_sidecar(
            rundir,
            graph=self._cfg_graph_block(),
            model=self._cfg_model_block(),
            run=self._cfg_run_block(),
        )

    # ------------------------------------------------------------------
    # The front door
    # ------------------------------------------------------------------
    @time_function_accumulate(auto_log=False)
    def run(
        self,
        tqdm_on: bool = False,
        steps: int | None = None,
        simref: float | None = None,
        verbose: bool = False,
        clean_export: bool = True,
        **solver_kw: Any,
    ) -> None:
        """Run the dynamics through the shared solver registry.

        Parameters
        ----------
        tqdm_on : bool
            Show a progress bar (in-process backends).
        steps : int, optional
            Number of Monte Carlo sweeps (overrides the configured horizon).
        simref : float, optional
            Size-normalised time (``steps = simref * N``); ``steps`` wins.
        verbose : bool
            Verbose backend output.
        clean_export : bool
            Remove the C backend's transient export files (state/edgelist)
            after a C run. Persisted observable files are kept.
        **solver_kw
            Extra keyword arguments forwarded to the solver's ``execute``.
        """
        if not type(self).solver_name:
            raise TypeError(
                f"{type(self).__name__} uses RunHostMixin but does not set "
                "the class attribute 'solver_name' (its solver-registry key)."
            )
        model = cast("DynSys", self)
        backend = self._resolve_backend()
        # Stashed for the lang= run-dirname token and the cfg sidecar.
        self._active_backend = backend
        solver = get_solver(type(self).solver_name, backend)
        # Guard before any file export / kernel call so backends never
        # silently ignore an unsupported configuration.
        solver.supports(model)
        self.check_attribute()
        self.initialize_run_parameters(steps=steps, simref=simref)
        self._begin_outputs()
        try:
            solver.execute(model, tqdm_on=tqdm_on, verbose=verbose, **solver_kw)
        finally:
            # Finalize disk artifacts (and always close any open stream).
            self._persist_observables()
        if backend is SolverBackend.C and clean_export:
            self.remove_run_c_files()
            self.sg.remove_exported_files()
