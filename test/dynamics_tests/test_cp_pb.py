"""pb-vs-py agreement for the SIR contact process (``_cp_native``).

Mirrors ``test_pb_signed.py``: same lattice, several seeds per backend,
compare the tail-averaged stationary active-site density (D-B6: statistical
agreement, never byte-equality — the RNG streams differ by construction).
The signed case (pflip = 0.15) verifies the kernel receives the SIGNED
weights: negative edges to active neighbours raise the recovery rate.
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")
pytest.importorskip(
    "lrgsglib.statsys.ContactProcess.ccore._cp_native",
    reason="_cp_native not built",
)

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.ContactProcess import ContactProcessEI, ContactProcessSIR

SIDE = 8
GRAPH_SEED = 7
MU = 0.5  # supercritical: recovery is slow against degree-4 infection
STEPS = 200
TAIL = 50
SEEDS = (1, 2, 3, 4)
TOL = 0.08
PFLIP_SIGNED = 0.15


def make_lattice(tmp_path, pflip=0.0):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=pflip,
        engine="nx",
        path_data=tmp_path,
        seed=GRAPH_SEED,  # same graph seed = same disorder realization
    )


def _make_sir(tmp_path, runlang, seed, pflip=0.0, **kw):
    sg = make_lattice(tmp_path, pflip=pflip)
    return ContactProcessSIR(
        sg,
        mu=MU,
        runlang=runlang,
        seed=seed,
        steps=STEPS,
        savedisk=False,
        savedensity=True,
        ic="uniform",
        **kw,
    )


def _tail_density(tmp_path, runlang, pflip):
    vals = []
    for sd in SEEDS:
        m = _make_sir(tmp_path, runlang, sd, pflip=pflip)
        m.run(tqdm_on=False)
        rho = np.asarray(m.density, dtype=np.float64)
        assert rho.shape == (STEPS,)  # record-before-sweep cadence
        vals.append(float(np.mean(rho[-TAIL:])))
    return float(np.mean(vals)), float(np.std(vals))


def _assert_backends_agree(tmp_path, pflip):
    d_py, s_py = _tail_density(tmp_path, "py", pflip)
    d_pb, s_pb = _tail_density(tmp_path, "pb", pflip)
    scatter = max(s_py, s_pb, 1e-3)
    assert abs(d_py - d_pb) < max(TOL, 4.0 * scatter), (
        f"SIR pb/py disagree (pflip={pflip}): "
        f"py={d_py:.4f}±{s_py:.4f} pb={d_pb:.4f}±{s_pb:.4f}"
    )


def test_cp_sir_pb_unsigned_agreement(tmp_path):
    _assert_backends_agree(tmp_path, pflip=0.0)


def test_cp_sir_pb_signed_agreement(tmp_path):
    """Negative edges raise the recovery rate — the kernel must see the
    signed weight array, not its absolute value."""
    _assert_backends_agree(tmp_path, pflip=PFLIP_SIGNED)


def test_cp_ei_pb_raises_not_implemented(tmp_path):
    sg = make_lattice(tmp_path)
    model = ContactProcessEI(sg, gamma=0.5, runlang="pb", steps=5, seed=1)
    with pytest.raises(NotImplementedError):
        model.run(tqdm_on=False)


def test_cp_sir_pb_capability_gates(tmp_path):
    """Hard capability errors (invariant #3): bipolar states and savedyn
    snapshots are not representable by the kernel."""
    sg = make_lattice(tmp_path)
    bipolar = ContactProcessSIR(
        sg, mu=MU, runlang="pb", steps=5, seed=1, state_type="bipolar"
    )
    with pytest.raises(NotImplementedError):
        bipolar.run(tqdm_on=False)

    snap = _make_sir(tmp_path, "pb", 1, savedyn=True)
    with pytest.raises(NotImplementedError):
        snap.run(tqdm_on=False)


def test_cp_sir_pb_seeded_reproducibility(tmp_path):
    def _final_state(seed):
        m = _make_sir(tmp_path, "pb", seed)
        m.run(tqdm_on=False)
        return np.asarray(m.s, dtype=np.int8).copy(), np.asarray(
            m.density, dtype=np.float64
        )

    s1, rho1 = _final_state(3)
    s2, rho2 = _final_state(3)
    np.testing.assert_array_equal(s1, s2)
    np.testing.assert_array_equal(rho1, rho2)
