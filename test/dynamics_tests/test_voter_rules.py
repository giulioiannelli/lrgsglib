"""Phase-3 tests for the VoterModel rule family + sampler axis + absorbing stop.

Axis A rules: linear / majority / qvoter / nonlinear (signed substrate).
Axis B samplers: asynchronous / synchronous / link.
Absorbing-state early stop: balanced => reaches consensus; frustrated => never.

References for the rule semantics live in
``.agents/ref/voter_dynamics/`` (Castellano 2009 RMP + PRE) and the absorbing /
frustration physics in ``iannelli2025topological.pdf`` (signed Laplacian).
"""

import pytest
import numpy as np


def _lat(tmp_path, side=10, pflip=0.0, seed=1):
    from lrgsglib.graphs.nx import Lattice2DNX
    return Lattice2DNX(
        side1=side, geo="sqr", pbc=True, pflip=pflip, seed=seed,
        path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False,
    )


# ===================================================================
# Axis A — every rule runs and produces a valid binary configuration
# ===================================================================


@pytest.mark.physical
@pytest.mark.parametrize("rule", ["linear", "majority", "qvoter", "nonlinear"])
def test_voter_rule_runs(tmp_path, rule):
    from lrgsglib.statsys import VoterModel

    v = VoterModel(sg=_lat(tmp_path), steps=30, runlang="py", seed=3,
                   rule=rule, savemagn=True)
    v.run(tqdm_on=False)
    assert v.s.size == v.N
    assert np.all(np.isin(v.s, (-1, 1)))
    assert len(v.magn) <= 30


@pytest.mark.physical
def test_nonlinear_alpha1_matches_linear(tmp_path):
    """alpha=1 makes the nonlinear rule coincide with the linear voter."""
    from lrgsglib.statsys import VoterModel

    def mean_abs_M(rule, **kw):
        acc = 0.0
        R = 120
        for r in range(R):
            v = VoterModel(sg=_lat(tmp_path, side=10, seed=r % 5), steps=50,
                           runlang="py", seed=r, rule=rule, **kw)
            v.run(tqdm_on=False)
            acc += abs(float(v.s.mean()))
        return acc / R

    lin = mean_abs_M("linear")
    nl = mean_abs_M("nonlinear", alpha=1.0)
    assert abs(lin - nl) < 0.12, f"<|M|> linear={lin:.3f} nonlinear(a=1)={nl:.3f}"


def test_rule_param_validation(tmp_path):
    from lrgsglib.statsys import VoterModel

    with pytest.raises(ValueError):
        VoterModel(sg=_lat(tmp_path, side=4), steps=5, runlang="py", rule="bogus")
    with pytest.raises(ValueError):
        VoterModel(sg=_lat(tmp_path, side=4), steps=5, runlang="py", q=0)
    with pytest.raises(ValueError):
        VoterModel(sg=_lat(tmp_path, side=4), steps=5, runlang="py", eps=1.5)
    with pytest.raises(ValueError):
        VoterModel(sg=_lat(tmp_path, side=4), steps=5, runlang="py", alpha=-1.0)


# ===================================================================
# Axis B — synchronous and link samplers
# ===================================================================


@pytest.mark.physical
@pytest.mark.parametrize("mode", ["asynchronous", "synchronous", "link"])
def test_voter_sampler_runs(tmp_path, mode):
    from lrgsglib.statsys import VoterModel

    v = VoterModel(sg=_lat(tmp_path), steps=25, runlang="py", seed=3,
                   upd_mode=mode, savemagn=True)
    v.run(tqdm_on=False)
    assert np.all(np.isin(v.s, (-1, 1)))
    assert len(v.magn) == 25


def test_link_requires_linear_rule(tmp_path):
    from lrgsglib.statsys import VoterModel

    for rule in ("majority", "qvoter", "nonlinear"):
        with pytest.raises(ValueError):
            VoterModel(sg=_lat(tmp_path, side=4), steps=5, runlang="py",
                       upd_mode="link", rule=rule)


# ===================================================================
# Absorbing-state early stop (signed substrate)
# ===================================================================


@pytest.mark.physical
def test_absorbing_balanced_reaches_consensus(tmp_path):
    """On a balanced (unsigned) graph the voter freezes at consensus."""
    from lrgsglib.graphs.nx import FullyConnectedNX
    from lrgsglib.statsys import VoterModel

    g = FullyConnectedNX(N=12, pflip=0.0, seed=1,
                         path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False)
    hit = 0
    for r in range(8):
        v = VoterModel(sg=g, steps=5000, runlang="py", seed=r,
                       absorbing_check=True, savemagn=True)
        v.run(tqdm_on=False)
        if v.absorbed_at is not None:
            hit += 1
            assert abs(float(v.s.mean())) == 1.0       # consensus
            assert v.count_frustrated_edges() == 0
            assert len(v.magn) == v.absorbed_at + 1    # stopped, recorded
    assert hit >= 6, f"only {hit}/8 reached consensus"


@pytest.mark.physical
def test_absorbing_frustrated_never_freezes(tmp_path):
    """On a frustrated signed graph no absorbing state exists."""
    from lrgsglib.statsys import VoterModel

    lf = _lat(tmp_path, side=12, pflip=0.3, seed=2)
    lf.flip_random_fract_edges()
    lf.compute_laplacian_spectrum_weigV()
    lam_min = float(lf.eigv[0])
    if lam_min <= 1e-10:
        pytest.skip("random signing happened to be balanced")
    v = VoterModel(sg=lf, steps=300, runlang="py", seed=2, absorbing_check=True)
    v.run(tqdm_on=False)
    assert v.absorbed_at is None
    assert v.count_frustrated_edges() > 0


def test_is_absorbing_helper(tmp_path):
    """count_frustrated_edges / is_absorbing agree with the edge definition."""
    from lrgsglib.statsys import VoterModel

    v = VoterModel(sg=_lat(tmp_path, side=6), steps=1, runlang="py", seed=1)
    v.init_voter_dynamics()
    consensus = np.ones(v.N, dtype=np.int8)
    assert v.is_absorbing(consensus)                 # unsigned: all-equal frozen
    assert v.count_frustrated_edges(consensus) == 0
    flipped = consensus.copy(); flipped[0] = -1
    assert not v.is_absorbing(flipped)               # one wrong spin -> active


# ===================================================================
# Native backends (C subprocess + pybind) — rule family / samplers / absorbing
# ===================================================================


def _pb_available() -> bool:
    try:
        from lrgsglib.statsys.VoterModel.ccore import _voter_native  # noqa: F401
        return True
    except Exception:
        return False


def _c_missing() -> bool:
    from lrgsglib.statsys.VoterModel import VoterModel
    return not (VoterModel._c_bin_dir / "VoterSimulator").exists()


def _skip_if_unavailable(runlang):
    if runlang == "pb" and not _pb_available():
        pytest.skip("_voter_native not built")
    if runlang == "C0" and _c_missing():
        pytest.skip("VoterSimulator not built")


@pytest.mark.integration
@pytest.mark.parametrize("runlang", ["pb", "C0"])
@pytest.mark.parametrize("rule", ["linear", "majority", "qvoter", "nonlinear"])
def test_native_rule_family_runs(tmp_path, runlang, rule):
    """Every rule runs on the native backends and yields a valid config."""
    from lrgsglib.statsys import VoterModel
    _skip_if_unavailable(runlang)
    v = VoterModel(sg=_lat(tmp_path, side=8), steps=20, runlang=runlang,
                   seed=1, rule=rule, savemagn=True)
    v.run(tqdm_on=False)
    assert v.s.size == v.N and np.all(np.isin(v.s, (-1, 1)))
    assert len(v.magn) <= 20


@pytest.mark.integration
@pytest.mark.parametrize("runlang", ["pb", "C0"])
@pytest.mark.parametrize("mode", ["asynchronous", "synchronous", "link"])
def test_native_samplers_run(tmp_path, runlang, mode):
    """async / synchronous / link all run on the native backends."""
    from lrgsglib.statsys import VoterModel
    _skip_if_unavailable(runlang)
    v = VoterModel(sg=_lat(tmp_path, side=8), steps=20, runlang=runlang,
                   seed=1, upd_mode=mode, savemagn=True)
    v.run(tqdm_on=False)
    assert np.all(np.isin(v.s, (-1, 1)))


@pytest.mark.parametrize("runlang", ["C0", "C0S"])
def test_c_subprocess_refuses_savedyn(tmp_path, runlang):
    """The C subprocess does not capture the in-memory savedyn trajectory."""
    from lrgsglib.statsys import VoterModel
    _skip_if_unavailable(runlang)
    v = VoterModel(sg=_lat(tmp_path, side=6), steps=5, runlang=runlang,
                   seed=1, savedyn=True)
    with pytest.raises(NotImplementedError):
        v.run(tqdm_on=False)


def test_pybind_records_savedyn(tmp_path):
    """The pybind backend records the full per-sweep trajectory in vm.s_t,
    sampled at the same points as the magnetization."""
    from lrgsglib.statsys import VoterModel
    _skip_if_unavailable("pb")
    v = VoterModel(sg=_lat(tmp_path, side=8), steps=30, runlang="pb",
                   upd_mode="gillespie", seed=1, savedyn=True,
                   savemagn=True, absorbing_check=True)
    v.init_voter_dynamics()
    v.run(tqdm_on=False)
    assert len(v.s_t) == len(v.magn)
    assert np.asarray(v.s_t[0]).shape == (v.N,)
    means = np.array([float(np.mean(s)) for s in v.s_t])
    assert np.allclose(means, np.array(v.magn))          # same record points
    assert np.array_equal(np.asarray(v.s_t[-1]), v.s)    # last == final state


# ===================================================================
# Gillespie rejection-free CTMC (Axis B, Phase 4 — Python reference)
# ===================================================================


@pytest.mark.physical
def test_gillespie_runs(tmp_path):
    """The rejection-free CTMC runs and yields a valid binary configuration."""
    from lrgsglib.statsys import VoterModel

    v = VoterModel(sg=_lat(tmp_path, side=10, pflip=0.1, seed=4), steps=40,
                   runlang="py", seed=4, upd_mode="gillespie",
                   savemagn=True)
    v.sg.flip_random_fract_edges()
    v.run(tqdm_on=False)
    assert v.s.size == v.N and np.all(np.isin(v.s, (-1, 1)))
    assert len(v.magn) <= 40


def test_gillespie_requires_linear_rule(tmp_path):
    """Rejection-free CTMC is a copy operation -> only the linear voter."""
    from lrgsglib.statsys import VoterModel

    for rule in ("majority", "qvoter", "nonlinear"):
        with pytest.raises(ValueError):
            VoterModel(sg=_lat(tmp_path, side=4), steps=5, runlang="py",
                       upd_mode="gillespie", rule=rule)


@pytest.mark.integration
@pytest.mark.parametrize("runlang", ["pb", "C0"])
def test_gillespie_native_runs(tmp_path, runlang):
    """The native (shared _ccore) CTMC kernel runs and yields a valid config."""
    from lrgsglib.statsys import VoterModel
    _skip_if_unavailable(runlang)
    v = VoterModel(sg=_lat(tmp_path, side=8, pflip=0.1, seed=1), steps=30,
                   runlang=runlang, seed=1, upd_mode="gillespie",
                   savemagn=True)
    v.sg.flip_random_fract_edges()
    v.run(tqdm_on=False)
    assert v.s.size == v.N and np.all(np.isin(v.s, (-1, 1)))
    assert len(v.magn) <= 30


@pytest.mark.physical
@pytest.mark.integration
@pytest.mark.parametrize("runlang", ["py", "pb", "C0"])
def test_gillespie_absorbing_parity_across_backends(tmp_path, runlang):
    """gillespie freezes at consensus on a balanced graph in every backend."""
    from lrgsglib.graphs.nx import FullyConnectedNX
    from lrgsglib.statsys import VoterModel
    if runlang != "py":
        _skip_if_unavailable(runlang)
    g = FullyConnectedNX(N=12, pflip=0.0, seed=1,
                         path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False)
    hit = 0
    for r in range(5):
        v = VoterModel(sg=g, steps=50000, runlang=runlang, seed=r,
                       upd_mode="gillespie", absorbing_check=True)
        v.run(tqdm_on=False)
        if v.absorbed_at is not None:
            hit += 1
            assert abs(float(v.s.mean())) == 1.0
    assert hit >= 4, f"{runlang}: only {hit}/5 gillespie runs reached consensus"


@pytest.mark.physical
def test_gillespie_absorbs_on_balanced(tmp_path):
    """On a balanced (unsigned) graph the CTMC freezes at consensus."""
    from lrgsglib.graphs.nx import FullyConnectedNX
    from lrgsglib.statsys import VoterModel

    g = FullyConnectedNX(N=12, pflip=0.0, seed=1,
                         path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False)
    hit = 0
    for r in range(8):
        v = VoterModel(sg=g, steps=50000, runlang="py", seed=r,
                       upd_mode="gillespie", absorbing_check=True,
                       savemagn=True)
        v.run(tqdm_on=False)
        if v.absorbed_at is not None:
            hit += 1
            assert abs(float(v.s.mean())) == 1.0
            assert v.count_frustrated_edges() == 0
            assert len(v.magn) == v.absorbed_at + 1
    assert hit >= 7, f"only {hit}/8 CTMC runs reached consensus"


@pytest.mark.physical
def test_gillespie_never_freezes_on_frustrated(tmp_path):
    """On a frustrated signed graph the CTMC never reaches an absorbing state."""
    from lrgsglib.statsys import VoterModel

    lf = _lat(tmp_path, side=10, pflip=0.3, seed=2)
    lf.flip_random_fract_edges()
    lf.compute_laplacian_spectrum_weigV()
    if float(lf.eigv[0]) <= 1e-10:
        pytest.skip("random signing happened to be balanced")
    v = VoterModel(sg=lf, steps=80, runlang="py", seed=2,
                   upd_mode="gillespie", absorbing_check=True)
    v.run(tqdm_on=False)
    assert v.absorbed_at is None
    assert v.count_frustrated_edges() > 0


@pytest.mark.physical
def test_gillespie_exit_probability(tmp_path):
    """Voter exit probability: on a regular graph P(+1 consensus) = up-fraction.

    The linear voter magnetization is conserved in expectation on a regular
    graph (complete graph), so starting from 7 up / 3 down spins the CTMC must
    reach the +1 consensus with probability ~0.7 (ref [1] Sec. II.A).
    """
    from lrgsglib.graphs.nx import FullyConnectedNX
    from lrgsglib.statsys import VoterModel

    g = FullyConnectedNX(N=10, pflip=0.0, seed=1,
                         path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False)
    ic = np.array([1] * 7 + [-1] * 3, dtype=np.int8)
    R, plus = 200, 0
    for r in range(R):
        v = VoterModel(sg=g, steps=50000, runlang="py", seed=r,
                       upd_mode="gillespie", absorbing_check=True)
        v.init_voter_dynamics()
        v.s = ic.copy()
        v.run(tqdm_on=False)
        assert v.absorbed_at is not None, "CTMC failed to absorb on complete graph"
        if float(v.s.mean()) == 1.0:
            plus += 1
    frac = plus / R
    assert abs(frac - 0.7) < 0.12, f"P(+1 consensus)={frac:.3f}, expected ~0.7"


@pytest.mark.physical
@pytest.mark.integration
@pytest.mark.parametrize("kernel", ["cr", "fenwick"])
def test_gillespie_cr_fenwick_parity(tmp_path, monkeypatch, kernel):
    """The native CTMC ships two interchangeable event-selection kernels:
    composition-rejection (default, O(deg)/event) and the Fenwick BKL fallback
    (O(deg log N)/event). ``LRGSG_CTMC_KERNEL`` forces one; both are the SAME
    exact CTMC and must reproduce the voter exit probability P(+1)=up-fraction
    on a regular graph (ref [1] Sec. II.A). 'fenwick' also exercises the env
    override on a low-degree graph the auto-dispatch would otherwise route to CR.
    """
    from lrgsglib.graphs.nx import FullyConnectedNX
    from lrgsglib.statsys import VoterModel
    _skip_if_unavailable("pb")
    monkeypatch.setenv("LRGSG_CTMC_KERNEL", kernel)

    g = FullyConnectedNX(N=10, pflip=0.0, seed=1,
                         path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False)
    ic = np.array([1] * 7 + [-1] * 3, dtype=np.int8)
    R, plus = 200, 0
    for r in range(R):
        v = VoterModel(sg=g, ic="custom", steps=50000, runlang="pb", seed=r,
                       upd_mode="gillespie", absorbing_check=True)
        v.init_voter_dynamics(custom=ic)
        v.run(tqdm_on=False)
        assert v.absorbed_at is not None, f"{kernel}: CTMC failed to absorb"
        if float(v.s.mean()) == 1.0:
            plus += 1
    frac = plus / R
    assert abs(frac - 0.7) < 0.12, f"{kernel}: P(+1)={frac:.3f}, expected ~0.7"


@pytest.mark.physical
@pytest.mark.integration
@pytest.mark.parametrize("runlang", ["py", "pb", "C0"])
def test_absorbing_parity_across_backends(tmp_path, runlang):
    """absorbing_check fires (consensus) on a balanced graph in every backend."""
    from lrgsglib.graphs.nx import FullyConnectedNX
    from lrgsglib.statsys import VoterModel
    _skip_if_unavailable(runlang) if runlang != "py" else None
    g = FullyConnectedNX(N=12, pflip=0.0, seed=1,
                         path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False)
    hit = 0
    for r in range(5):
        v = VoterModel(sg=g, steps=8000, runlang=runlang, seed=r,
                       absorbing_check=True)
        v.run(tqdm_on=False)
        if v.absorbed_at is not None:
            hit += 1
            assert abs(float(v.s.mean())) == 1.0
    assert hit >= 4, f"{runlang}: only {hit}/5 reached an absorbing consensus"


# ===================================================================
# Vectorized backends (np / cu) — Phase 5
# Synchronous linear voter via a single CSR gather per sweep.
# ===================================================================


def _cu_unavailable() -> bool:
    try:
        import cupy
        cupy.zeros(1)          # force a context (catches no-GPU / driver issues)
        return False
    except Exception:
        return True


@pytest.mark.physical
@pytest.mark.parametrize("runlang", ["np", "cu"])
def test_vectorized_runs(tmp_path, runlang):
    """The vectorized synchronous linear voter runs and yields a valid config."""
    from lrgsglib.statsys import VoterModel
    if runlang == "cu" and _cu_unavailable():
        pytest.skip("cupy / GPU not available")
    lat = _lat(tmp_path, side=16, pflip=0.1, seed=1)
    lat.flip_random_fract_edges()
    v = VoterModel(sg=lat, steps=40, runlang=runlang, seed=1,
                   savemagn=True)
    v.run(tqdm_on=False)
    assert v.s.size == v.N and np.all(np.isin(v.s, (-1, 1)))
    assert len(v.magn) <= 40


@pytest.mark.parametrize("runlang", ["np", "cu"])
def test_vectorized_guards(tmp_path, runlang):
    """Vectorized backends are synchronous-linear only: reject rule family,
    savedyn, and the intrinsically-sequential schedules."""
    from lrgsglib.statsys import VoterModel
    if runlang == "cu" and _cu_unavailable():
        pytest.skip("cupy / GPU not available")
    for kw in (dict(rule="majority"), dict(savedyn=True),
               dict(upd_mode="gillespie"), dict(upd_mode="link")):
        with pytest.raises(NotImplementedError):
            VoterModel(sg=_lat(tmp_path, side=6), steps=5, runlang=runlang,
                       seed=1, **kw).run(tqdm_on=False)


@pytest.mark.physical
@pytest.mark.parametrize("runlang", ["np", "cu"])
def test_vectorized_absorbs_on_balanced(tmp_path, runlang):
    """On a balanced graph the synchronous voter freezes at consensus."""
    from lrgsglib.graphs.nx import FullyConnectedNX
    from lrgsglib.statsys import VoterModel
    if runlang == "cu" and _cu_unavailable():
        pytest.skip("cupy / GPU not available")
    g = FullyConnectedNX(N=24, pflip=0.0, seed=1,
                         path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False)
    hit = 0
    for r in range(5):
        v = VoterModel(sg=g, steps=20000, runlang=runlang, seed=r,
                       absorbing_check=True, savemagn=True)
        v.run(tqdm_on=False)
        if v.absorbed_at is not None:
            hit += 1
            assert abs(float(v.s.mean())) == 1.0
            assert len(v.magn) == v.absorbed_at + 1
    assert hit >= 4, f"{runlang}: only {hit}/5 reached consensus"


@pytest.mark.physical
@pytest.mark.parametrize("runlang", ["np", "cu"])
def test_vectorized_martingale(tmp_path, runlang):
    """Correctness (not just 'it runs'): the vectorized synchronous linear voter
    conserves E[M] on a complete graph, so <M_final> over many runs equals the
    initial magnetization (the voter martingale). Validates the np / cu kernel
    against the exact conservation law -- the CPU/GPU parity anchor."""
    from lrgsglib.graphs.nx import FullyConnectedNX
    from lrgsglib.statsys import VoterModel
    if runlang == "cu" and _cu_unavailable():
        pytest.skip("cupy / GPU not available")

    N, n_up = 10, 6
    ic = np.full(N, -1, dtype=np.int8)
    ic[:n_up] = 1
    M_init = float(ic.mean())          # 0.2
    R, tol = 200, 0.30

    g = FullyConnectedNX(
        N=N, pflip=0.0, seed=1,
        path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False,
    )
    acc = 0.0
    for r in range(1, R + 1):
        v = VoterModel(sg=g, ic="custom", steps=100, runlang=runlang, seed=r)
        v.init_voter_dynamics(custom=ic)
        v.run(tqdm_on=False)
        acc += float(np.mean(v.s))
    assert abs(acc / R - M_init) < tol, \
        f"{runlang}: <M_final>={acc / R:.3f}, M_init={M_init}"
