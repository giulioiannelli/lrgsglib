"""Cross-model parity suite (D-B6) — run -> record -> output contract.

Phase-0 consumers: VoterModel and ContactProcessSIR (the two models already
on the shared spine), Python backend. Each model joins this suite as it
modernizes (Ising with the Phase-1 vertical slice, then the spin quad);
native backends join the agreement matrix as the CMake fold (D-B7) brings
their binaries under test.
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.statsys.ContactProcess import ContactProcessSIR
from lrgsglib.statsys.ContactProcess.defaults import (
    CP_OBS_DENSITY,
    CP_OBS_SNAPSHOTS,
)
from lrgsglib.statsys.HeisenbergModel import HeisenbergMetropolis
from lrgsglib.statsys.HeisenbergModel.defaults import (
    HEISENBERG_OBS_ENERGY,
    HEISENBERG_OBS_MAGN,
    HEISENBERG_OBS_SNAPSHOTS,
)
from lrgsglib.statsys.IsingDynamics import (
    IsingCEM,
    IsingMetropolis,
    IsingParallelTempering,
    IsingSimulatedAnnealing,
)
from lrgsglib.statsys.IsingDynamics.defaults import (
    ISING_OBS_ENERGY,
    ISING_OBS_MAGN,
    ISING_OBS_SNAPSHOTS,
)
from lrgsglib.statsys.PottsModel import PottsMetropolis
from lrgsglib.statsys.PottsModel.defaults import (
    POTTS_OBS_ENERGY,
    POTTS_OBS_MAGN,
    POTTS_OBS_SNAPSHOTS,
)
from lrgsglib.statsys.VoterModel import VoterModel
from lrgsglib.statsys.VoterModel.defaults import VOTER_OBS_MAGN, VOTER_OBS_SOUT
from lrgsglib.statsys.XYModel import XYMetropolis
from lrgsglib.statsys.XYModel.defaults import (
    XY_OBS_ENERGY,
    XY_OBS_MAGN,
    XY_OBS_SNAPSHOTS,
)

from .parity_harness import (
    PARITY_SEED,
    PARITY_STEPS,
    ParityCase,
    assert_seeded_reproducibility,
    run_record_output_contract,
)


def _voter_case() -> ParityCase:
    return ParityCase(
        name="voter",
        make_model=lambda sg: VoterModel(
            sg,
            runlang="py",
            savedisk=True,
            savedyn=True,
            savemagn=True,
            seed=PARITY_SEED,
        ),
        run_kwargs=dict(tqdm_on=False, steps=PARITY_STEPS),
        expected_observables=(VOTER_OBS_MAGN, VOTER_OBS_SOUT),
        expected_lengths={VOTER_OBS_MAGN: PARITY_STEPS},
    )


def _cp_case() -> ParityCase:
    return ParityCase(
        name="contact_process_sir",
        make_model=lambda sg: ContactProcessSIR(
            sg,
            runlang="py",
            savedisk=True,
            savedyn=True,
            savedensity=True,
            seed=PARITY_SEED,
        ),
        run_kwargs=dict(tqdm_on=False, steps=PARITY_STEPS),
        expected_observables=(CP_OBS_DENSITY, CP_OBS_SNAPSHOTS),
        expected_lengths={CP_OBS_DENSITY: PARITY_STEPS},
    )


def _ising_case() -> ParityCase:
    return ParityCase(
        name="ising_metropolis",
        make_model=lambda sg: IsingMetropolis(
            sg,
            T=2.0,
            runlang="py",
            savedisk=True,
            observables=(ISING_OBS_ENERGY, ISING_OBS_MAGN, ISING_OBS_SNAPSHOTS),
            seed=PARITY_SEED,
        ),
        run_kwargs=dict(tqdm_on=False, steps=PARITY_STEPS),
        expected_observables=(
            ISING_OBS_ENERGY,
            ISING_OBS_MAGN,
            ISING_OBS_SNAPSHOTS,
        ),
        expected_lengths={
            ISING_OBS_ENERGY: PARITY_STEPS,
            ISING_OBS_MAGN: PARITY_STEPS,
        },
    )


def _ising_sa_case() -> ParityCase:
    # SA's horizon is schedule-derived (run() rejects steps=): shape the
    # schedule so the record count equals PARITY_STEPS.
    return ParityCase(
        name="ising_sa",
        make_model=lambda sg: IsingSimulatedAnnealing(
            sg,
            T_init=3.0,
            T_final=0.1,
            schedule="linear",
            n_temperatures=PARITY_STEPS,
            steps_per_T=1,
            runlang="py",
            savedisk=True,
            seed=PARITY_SEED,
        ),
        run_kwargs=dict(tqdm_on=False),
        expected_observables=(ISING_OBS_ENERGY, ISING_OBS_MAGN),
        expected_lengths={
            ISING_OBS_ENERGY: PARITY_STEPS,
            ISING_OBS_MAGN: PARITY_STEPS,
        },
    )


def _ising_pt_case() -> ParityCase:
    return ParityCase(
        name="ising_pt",
        make_model=lambda sg: IsingParallelTempering(
            sg,
            n_replicas=3,
            T_min=0.5,
            T_max=3.0,
            steps_per_exchange=1,
            n_exchanges=PARITY_STEPS,
            runlang="py",
            savedisk=True,
            seed=PARITY_SEED,
        ),
        run_kwargs=dict(tqdm_on=False),
        expected_observables=(ISING_OBS_ENERGY, ISING_OBS_MAGN),
        expected_lengths={
            ISING_OBS_ENERGY: PARITY_STEPS,
            ISING_OBS_MAGN: PARITY_STEPS,
        },
    )


def _ising_cem_case() -> ParityCase:
    return ParityCase(
        name="ising_cem",
        make_model=lambda sg: IsingCEM(
            sg,
            n_modes=4,
            pop_size=8,
            n_iter=PARITY_STEPS,
            restarts=1,
            greedy=False,
            runlang="py",
            savedisk=True,
            seed=PARITY_SEED,
        ),
        run_kwargs=dict(tqdm_on=False),
        expected_observables=(ISING_OBS_ENERGY, ISING_OBS_MAGN),
        expected_lengths={
            ISING_OBS_ENERGY: PARITY_STEPS,
            ISING_OBS_MAGN: PARITY_STEPS,
        },
    )


def _potts_case() -> ParityCase:
    return ParityCase(
        name="potts_metropolis",
        make_model=lambda sg: PottsMetropolis(
            sg,
            q=3,
            T=2.0,
            runlang="py",
            savedisk=True,
            observables=(
                POTTS_OBS_ENERGY,
                POTTS_OBS_MAGN,
                POTTS_OBS_SNAPSHOTS,
            ),
            seed=PARITY_SEED,
        ),
        run_kwargs=dict(tqdm_on=False, steps=PARITY_STEPS),
        expected_observables=(
            POTTS_OBS_ENERGY,
            POTTS_OBS_MAGN,
            POTTS_OBS_SNAPSHOTS,
        ),
        expected_lengths={
            POTTS_OBS_ENERGY: PARITY_STEPS,
            POTTS_OBS_MAGN: PARITY_STEPS,
        },
    )


def _heisenberg_case() -> ParityCase:
    return ParityCase(
        name="heisenberg_metropolis",
        make_model=lambda sg: HeisenbergMetropolis(
            sg,
            T=1.0,
            runlang="py",
            savedisk=True,
            observables=(
                HEISENBERG_OBS_ENERGY,
                HEISENBERG_OBS_MAGN,
                HEISENBERG_OBS_SNAPSHOTS,
            ),
            seed=PARITY_SEED,
        ),
        run_kwargs=dict(tqdm_on=False, steps=PARITY_STEPS),
        expected_observables=(
            HEISENBERG_OBS_ENERGY,
            HEISENBERG_OBS_MAGN,
            HEISENBERG_OBS_SNAPSHOTS,
        ),
        expected_lengths={
            HEISENBERG_OBS_ENERGY: PARITY_STEPS,
            HEISENBERG_OBS_MAGN: PARITY_STEPS,
        },
    )


def _xy_case() -> ParityCase:
    return ParityCase(
        name="xy_metropolis",
        make_model=lambda sg: XYMetropolis(
            sg,
            T=1.0,
            runlang="py",
            savedisk=True,
            observables=(XY_OBS_ENERGY, XY_OBS_MAGN, XY_OBS_SNAPSHOTS),
            seed=PARITY_SEED,
        ),
        run_kwargs=dict(tqdm_on=False, steps=PARITY_STEPS),
        expected_observables=(XY_OBS_ENERGY, XY_OBS_MAGN, XY_OBS_SNAPSHOTS),
        expected_lengths={
            XY_OBS_ENERGY: PARITY_STEPS,
            XY_OBS_MAGN: PARITY_STEPS,
        },
    )


CASES = {
    "voter": _voter_case,
    "contact_process_sir": _cp_case,
    "ising_metropolis": _ising_case,
    "ising_sa": _ising_sa_case,
    "ising_pt": _ising_pt_case,
    "ising_cem": _ising_cem_case,
    "potts_metropolis": _potts_case,
    "xy_metropolis": _xy_case,
    "heisenberg_metropolis": _heisenberg_case,
}


@pytest.mark.parametrize("case_name", sorted(CASES))
def test_run_record_output_contract(case_name, tmp_path):
    run_record_output_contract(CASES[case_name](), tmp_path)


@pytest.mark.parametrize("case_name", sorted(CASES))
def test_seeded_reproducibility(case_name, tmp_path):
    assert_seeded_reproducibility(CASES[case_name](), tmp_path)
