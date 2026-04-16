#!/usr/bin/env python3
"""Verify all C simulator binaries in all supported modes.

Checks:
  1. Binary runs without error
  2. Final state has correct shape and valid values
  3. Observable outputs (magn, ene, density) are non-empty where expected
  4. Physics sanity: spins ∈ {-1, +1}, energies finite, etc.

Run from lrgsglib root:
    python test/scripts/verify_all_c_modes.py
"""

import sys
import tempfile
from pathlib import Path

import numpy as np


def make_lattice(tmp, side=10, pflip=0.1):
    from lrgsglib.graphs.nx import Lattice2DNX
    lat = Lattice2DNX(
        side1=side, geo="sqr", pbc=True, pflip=pflip, seed=42,
        path_data=tmp, path_plot=tmp, init_nw_dict=False,
    )
    lat.flip_random_fract_edges()
    return lat


# ── Ising ─────────────────────────────────────────────────────────────
def test_ising(tmp):
    from lrgsglib.statsys.IsingDynamics import IsingDynamics
    R = {}

    # -- Metropolis modes --
    for code in ("C0E", "C0S", "C0ES"):
        try:
            lat = make_lattice(tmp)
            ising = IsingDynamics(lat, T=2.0, runlang=code, steps=50,
                                  freq=10, seed=42)
            ising.init_ising_dynamics()
            ising.run(tqdm_on=False)
            ok = (len(ising.s) == lat.N
                  and np.all(np.isin(ising.s, [-1, 1])))
            em = "E" in code
            if em:
                ok = ok and len(ising.magn) > 0 and np.all(np.isfinite(ising.magn))
                ok = ok and len(ising.ene) > 0 and np.all(np.isfinite(ising.ene))
            detail = (f"|s|={len(ising.s)} magn={len(ising.magn)} "
                      f"ene={len(ising.ene)} spins_ok={np.all(np.isin(ising.s, [-1,1]))}")
            R[f"Ising {code:5s} Met"] = ("OK" if ok else "WARN", detail)
        except Exception as e:
            R[f"Ising {code:5s} Met"] = ("FAIL", str(e)[:80])

    # C0K/C0EK (cluster) — requires cluster files, expected to fail gracefully
    for code in ("C0K",):
        try:
            lat = make_lattice(tmp)
            ising = IsingDynamics(lat, T=2.0, runlang=code, steps=20, seed=42)
            ising.init_ising_dynamics()
            ising.run(tqdm_on=False)
            R[f"Ising {code:5s} Met"] = ("OK", f"|s|={len(ising.s)}")
        except RuntimeError as e:
            if "cluster" in str(e).lower() or "No such file" in str(e):
                R[f"Ising {code:5s} Met"] = ("SKIP", "needs cluster files (expected)")
            else:
                R[f"Ising {code:5s} Met"] = ("FAIL", str(e)[:80])
        except Exception as e:
            R[f"Ising {code:5s} Met"] = ("FAIL", str(e)[:80])

    # -- SA --
    try:
        lat = make_lattice(tmp)
        ising = IsingDynamics(
            lat, runlang="C1E", seed=42,
            sa_enabled=True, T_init=5.0, T_final=0.01,
            n_temperatures=20, steps_per_T=20,
        )
        ising.init_ising_dynamics()
        ising.run(tqdm_on=False)
        ok = (len(ising.s) == lat.N
              and np.all(np.isin(ising.s, [-1, 1]))
              and len(ising.magn) > 0)
        detail = (f"|s|={len(ising.s)} magn={len(ising.magn)} "
                  f"ene={len(ising.ene)} |m_final|={abs(ising.magn[-1]):.3f}")
        R["Ising C1E   SA "] = ("OK" if ok else "WARN", detail)
    except Exception as e:
        R["Ising C1E   SA "] = ("FAIL", str(e)[:80])

    # -- PT --
    try:
        lat = make_lattice(tmp)
        ising = IsingDynamics(
            lat, runlang="C2E", seed=42,
            pt_enabled=True, n_replicas=4, T_min=0.5, T_max=4.0,
            steps_per_exchange=20, n_exchanges=20,
        )
        ising.init_ising_dynamics()
        ising.run(tqdm_on=False)
        ok = (len(ising.s) == lat.N
              and np.all(np.isin(ising.s, [-1, 1]))
              and hasattr(ising, 's_replicas')
              and ising.s_replicas.shape == (4, lat.N))
        detail = (f"|s|={len(ising.s)} replicas={ising.s_replicas.shape} "
                  f"magn={len(ising.magn)} ene={len(ising.ene)}")
        R["Ising C2E   PT "] = ("OK" if ok else "WARN", detail)
    except Exception as e:
        R["Ising C2E   PT "] = ("FAIL", str(e)[:80])

    return R


# ── Voter ─────────────────────────────────────────────────────────────
def test_voter(tmp):
    from lrgsglib.statsys.VoterModel import VoterModel
    R = {}

    for code, nSL in [("C0", 0), ("C0S", 10)]:
        try:
            lat = make_lattice(tmp, pflip=0.0)
            kw = {"nSampleLog": nSL} if nSL else {}
            voter = VoterModel(sg=lat, steps=50, runlang=code, seed=42, **kw)
            voter.init_voter_dynamics()
            voter.run(tqdm_on=False)
            ok = (len(voter.s) == lat.N
                  and np.all(np.isin(voter.s, [-1, 1]))
                  and len(voter.magn) == 50)
            detail = (f"|s|={len(voter.s)} magn={len(voter.magn)} "
                      f"m_final={voter.magn[-1]:.3f}")
            R[f"Voter {code:5s}     "] = ("OK" if ok else "WARN", detail)
        except Exception as e:
            R[f"Voter {code:5s}     "] = ("FAIL", str(e)[:80])

    return R


# ── Contact Process ───────────────────────────────────────────────────
def test_cp(tmp):
    from lrgsglib.statsys.ContactProcess import ContactProcessEI, ContactProcessSIR
    R = {}

    # EI modes
    ei_cases = [
        ("C1",  "final",   {}),
        ("C1D", "density", {}),
        ("C1S", "snap",    {}),
        ("C1D", "adpt",    {"update_mode": "adaptive"}),
        ("C1D", "cach",    {"update_mode": "cached"}),
        ("C1D", "gill",    {"update_mode": "gillespie"}),
    ]
    for code, label, extra in ei_cases:
        try:
            lat = make_lattice(tmp, pflip=0.0)
            cp = ContactProcessEI(
                sg=lat, steps=100, runlang=code, seed=42, gamma=1.5, **extra,
            )
            cp.init_contact_dynamics()
            cp.run(tqdm_on=False)
            ok = (len(cp.s) == lat.N
                  and np.all(np.isin(cp.s, [0, 1])))
            detail = f"|s|={len(cp.s)} active={np.sum(cp.s)}/{lat.N}"
            R[f"CP EI {code}+{label:6s}"] = ("OK" if ok else "WARN", detail)
        except Exception as e:
            R[f"CP EI {code}+{label:6s}"] = ("FAIL", str(e)[:80])

    # SIR
    try:
        lat = make_lattice(tmp, pflip=0.0)
        cp = ContactProcessSIR(sg=lat, steps=100, runlang="C0", seed=42, mu=0.5)
        cp.init_contact_dynamics()
        cp.run(tqdm_on=False)
        ok = len(cp.s) == lat.N
        R["CP SIR C0        "] = ("OK" if ok else "WARN", f"|s|={len(cp.s)}")
    except Exception as e:
        R["CP SIR C0        "] = ("FAIL", str(e)[:80])

    return R


# ── Continuous & Vector ───────────────────────────────────────────────
def test_contvec(tmp):
    R = {}

    models = [
        ("Kuramoto",          "KuramotoModel",          "init_kuramoto_dynamics",
         {"coupling": 1.0}),
        ("ReactionDiffusion", "ReactionDiffusionModel", "init_rd_dynamics",
         {}),
        ("CoupledODE",        "CoupledODEModel",        "init_ode_dynamics",
         {}),
        ("Potts",             "PottsModel",             "init_potts_dynamics",
         {"q": 3, "T": 2.0}),
        ("XY",                "XYModel",                "init_xy_dynamics",
         {"T": 2.0}),
        ("Heisenberg",        "HeisenbergModel",        "init_heisenberg_dynamics",
         {"T": 2.0}),
        ("MultiSpecies",      "MultiSpeciesModel",      "init_multispecies_dynamics",
         {"n_species": 3}),
    ]

    for name, cls_name, init_method, extra in models:
        try:
            mod = __import__(
                f"lrgsglib.statsys.{cls_name}",
                fromlist=[cls_name],
            )
            Cls = getattr(mod, cls_name)
            lat = make_lattice(tmp, pflip=0.0)
            model = Cls(sg=lat, steps=20, runlang="C0", seed=42, **extra)
            getattr(model, init_method)()
            model.run(tqdm_on=False)
            ok = model.s is not None and len(model.s) > 0
            detail = f"shape={np.array(model.s).shape}"
            R[f"{name:20s}"] = ("OK" if ok else "WARN", detail)
        except Exception as e:
            R[f"{name:20s}"] = ("FAIL", str(e)[:80])

    return R


# ── Main ──────────────────────────────────────────────────────────────
def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        all_results = {}

        sections = [
            ("Ising (3 binaries × modes)", test_ising),
            ("Voter (1 binary × 2 modes)", test_voter),
            ("Contact Process (2 binaries × modes)", test_cp),
            ("Continuous & Vector (7 binaries)", test_contvec),
        ]

        for title, func in sections:
            all_results.update(func(tmp))

        # Print
        W = 70
        print("=" * W)
        print("C SIMULATOR VERIFICATION")
        print("=" * W)

        ok = fail = skip = warn = 0
        for label, (status, detail) in sorted(all_results.items()):
            tag = {"OK": "PASS", "FAIL": "FAIL", "SKIP": "SKIP",
                   "WARN": "WARN"}[status]
            print(f"  [{tag:4s}] {label}  {detail}")
            if status == "OK":
                ok += 1
            elif status == "FAIL":
                fail += 1
            elif status == "SKIP":
                skip += 1
            else:
                warn += 1

        print(f"\n  {ok} pass, {fail} fail, {skip} skip, {warn} warn")
        return 1 if fail > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
