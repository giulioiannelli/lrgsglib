"""The shared run-dirname builder (_naming.py) — Phase C naming contract.

Pins the exact folder names approved in the 2026-07-07 interview:
key=value for every token, floats %.3g, defaults elided (seed/backend
always written), inapplicable axes absent, heterogeneous q joined '-'.
"""

from __future__ import annotations

import pytest

from lrgsglib.statsys._naming import (
    ALWAYS,
    VOCAB,
    format_value,
    parse_run_dirname,
    run_dirname,
)


def test_plain_metropolis_run_elides_all_defaults():
    name = run_dirname(
        [
            ("p", 0.1),
            ("T", 2.27),
            ("h", 0.0, 0.0),
            ("ns", 5000),
            ("se", 0, 0),
            ("sch", "metropolis", "metropolis"),
            ("rule", "metropolis", "metropolis"),
            ("move", "single", "single"),
            ("upd", "asynchronous", "asynchronous"),
            ("ord", "random", "random"),
            ("ic", "uniform", "uniform"),
            ("lang", "pb"),
            ("s", 91731),
        ]
    )
    assert name == "p=0.1_T=2.27_ns=5000_lang=pb_s=91731"


def test_non_defaults_become_visible():
    name = run_dirname(
        [
            ("p", 0.1),
            ("T", 1.5),
            ("ns", 5000),
            ("rule", "glauber", "metropolis"),
            ("move", "exchange", "single"),
            ("ord", "permutation", "random"),
            ("lang", "py"),
            ("s", 42),
        ]
    )
    assert name == (
        "p=0.1_T=1.5_ns=5000_rule=glau_move=exch_ord=perm_lang=py_s=42"
    )


def test_inapplicable_axis_is_absent():
    # A wolff run passes rule=None: cluster acceptance is baked into the
    # bond probability, so no rule= token may ever appear.
    name = run_dirname(
        [
            ("p", 0.1),
            ("T", 0.85),
            ("ns", 3000),
            ("rule", None, "metropolis"),
            ("move", "wolff", "single"),
            ("lang", "py"),
            ("s", 23),
        ]
    )
    assert name == "p=0.1_T=0.85_ns=3000_move=wolff_lang=py_s=23"
    assert "rule=" not in name


def test_sa_and_pt_fingerprints():
    sa = run_dirname(
        [
            ("p", 0.1),
            ("Ti", 10.0),
            ("Tf", 0.01),
            ("nT", 100),
            ("spt", 10),
            ("sch", "simulated_annealing"),
            ("lang", "pb"),
            ("s", 7),
        ]
    )
    assert sa == "p=0.1_Ti=10_Tf=0.01_nT=100_spt=10_sch=sa_lang=pb_s=7"

    pt = run_dirname(
        [
            ("p", 0.1),
            ("Tmin", 1.0),
            ("Tmax", 4.0),
            ("nr", 8),
            ("spx", 10),
            ("nx", 500),
            ("sch", "parallel_tempering"),
            ("lang", "py"),
            ("s", 3),
        ]
    )
    assert pt == "p=0.1_Tmin=1_Tmax=4_nr=8_spx=10_nx=500_sch=pt_lang=py_s=3"


def test_heterogeneous_q_joined_with_dashes():
    name = run_dirname(
        [
            ("p", 0.1),
            ("k", 2),
            ("q", (2, 5)),
            ("T", 1.0),
            ("ns", 3000),
            ("lang", "py"),
            ("s", 31),
        ]
    )
    assert name == "p=0.1_k=2_q=2-5_T=1_ns=3000_lang=py_s=31"


def test_float_formatting_is_3g():
    assert format_value(2.269185) == "2.27"
    assert format_value(1.0) == "1"
    assert format_value(0.001) == "0.001"
    assert format_value(1e-05) == "1e-05"


def test_vocab_shortens_categoricals_and_passes_unknown():
    assert format_value("glauber") == "glau"
    assert format_value("gillespie") == "gill"
    assert format_value("wolff") == "wolff"
    assert format_value("sw") == "sw"
    for long, short in VOCAB.items():
        assert format_value(long) == short


def test_always_sentinel_never_elides():
    name = run_dirname([("lang", "py", ALWAYS), ("s", 1, ALWAYS)])
    assert name == "lang=py_s=1"


def test_parse_round_trip():
    name = "p=0.1_T=2.27_q=2-5_rule=glau_lang=pb_s=91731"
    parsed = parse_run_dirname(name)
    assert parsed == {
        "p": "0.1",
        "T": "2.27",
        "q": "2-5",
        "rule": "glau",
        "lang": "pb",
        "s": "91731",
    }
    with pytest.raises(ValueError, match="malformed"):
        parse_run_dirname("p=0.1_oops")
