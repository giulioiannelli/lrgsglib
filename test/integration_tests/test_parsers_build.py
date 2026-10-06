"""Program parsers: every one builds, graph programs share one ``-ge`` flag, and
``-ge gt`` is refused wherever its graph-tool path is not ported yet."""

import importlib
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[2] / "src"
sys.path.insert(0, str(SRC))

from lrgsglib.config.progargs.common import graph_engine_opt_args  # noqa: E402
from lrgsglib.config.progargs.defs.generic import (  # noqa: E402
    DEFAULT_GRAPH_ENGINE,
    GRAPH_ENGINE_CHOICES,
)
from lrgsglib.config.progargs.defs.RandomRegular import (  # noqa: E402
    DEFAULT_RRG_ENGINE,
)
from lrgsglib.graphs import _engine  # noqa: E402
from parsers.shared import (  # noqa: E402
    GraphArgumentParser,
    add_graph_engine_arg,
)

pytestmark = pytest.mark.integration

ENGINE_FLAGS = ("-ge", "--graph_engine")
REFUSED = "is not supported here yet"

# Programs that build graphs, and serializers that forward -ge to them
GRAPH_PARSERS = {
    "ER_IsingDynamics",
    "L2D_ContactProcess",
    "L2D_IsingDynamics",
    "L2D_Recon",
    "L2D_SignedRW",
    "L2D_SlaplSpect",
    "L2D_TransCluster",
    "L2D_VoterModel",
    "L3D_IsingDynamics",
    "L3D_Recon",
    "L3D_TransCluster",
    "MCG_SlaplSpect",
    "RRG_TransCluster",
    "SCS_TransCluster",
    "L2D_VoterModel_Serializer",
    "L3D_IsingDynamics_Serializer",
    "L3D_TopoComparison_Serializer",
}
# Serializers that pass unknown flags (-ge included) through to the program
OTHER_PARSERS = {
    "ER_IsingDynamics_Serializer",
    "L2D_ContactProcess_Serializer",
    "L2D_IsingDynamics_Serializer",
    "L2D_SlaplSpect_Serializer",
    "L2D_TransCluster_Serializer",
    "L3D_Recon_srun",
    "L3D_TransCluster_Serializer",
    "MCG_SlaplSpect_Serializer",
    "RRG_TransCluster_Serializer",
    "SCS_TransCluster_Serializer",
}
ALL_PARSERS = sorted(
    p.stem
    for p in (SRC / "parsers").glob("*.py")
    if p.stem not in {"__init__", "shared"}
)

# Smallest argv each program accepts on both engines
BASE_ARGV = {
    "ER_IsingDynamics": ["100", "0.02", "0.1", "-rl", "py"],
    "L2D_ContactProcess": ["16"],
    "L2D_IsingDynamics": ["16"],
    "L2D_Recon": ["16"],
    "L2D_SignedRW": ["16"],
    "L2D_SlaplSpect": ["16", "-m", "eigvals"],
    "L2D_TransCluster": ["16"],
    "L2D_VoterModel": ["16"],
    "L3D_IsingDynamics": ["8"],
    "L3D_Recon": ["8"],
    "L3D_TransCluster": ["8"],
    "MCG_SlaplSpect": ["0.8", "0.6", "0", "0", "0.4", "3"],
    "SCS_TransCluster": ["64"],
}
# graph-tool path not ported: refused whatever the options
ALWAYS_REFUSED = [
    "L2D_Recon",
    "L2D_TransCluster",
    "L3D_Recon",
    "L3D_TransCluster",
    "SCS_TransCluster",
]
CELL = (["-c", "randXERR"], "-c rand")
WORKDIR = (["-wd", "elsewhere"], "-wd")
# (program, extra argv, words of the reason) refused on gt, accepted on nx
REFUSED_WITH = [
    ("ER_IsingDynamics", *CELL),
    ("ER_IsingDynamics", *WORKDIR),
    ("ER_IsingDynamics", ["-rl", "C0E"], "C backends"),
    ("L2D_ContactProcess", *CELL),
    ("L2D_ContactProcess", *WORKDIR),
    ("L2D_IsingDynamics", *CELL),
    ("L2D_IsingDynamics", *WORKDIR),
    ("L2D_IsingDynamics", ["--save-results"], "--save-results"),
    ("L2D_SignedRW", *CELL),
    ("L2D_SignedRW", *WORKDIR),
    ("L2D_VoterModel", *CELL),
    ("L2D_VoterModel", *WORKDIR),
    ("L3D_IsingDynamics", *CELL),
    ("L3D_IsingDynamics", *WORKDIR),
    ("L3D_IsingDynamics", ["--save-results"], "--save-results"),
    ("L3D_IsingDynamics", ["--pdil", "0.1"], "--pdil"),
    ("L2D_SlaplSpect", ["-m", "eigvec_dist"], "eigvec_dist"),
    ("MCG_SlaplSpect", ["-m", "eigvec_dist"], "eigvec_dist"),
]
ACCEPTED_ON_GT = sorted(set(BASE_ARGV) - set(ALWAYS_REFUSED))

# Programs that build their own parser: (module, parse function, base argv)
SIDE = ["--side", "8", "--pflip-list", "0.1"]
INLINE = {
    "L2D_SignedRW_Bubbles": ("_parse_args", SIDE),
    "L2D_SignedRW_BubblesRand": ("_parse_args", SIDE),
    "L2D_SignedRW_TimeToTouch": ("_parse_args", SIDE + ["--n-max", "4"]),
    "L2D_SignedRW_Bubbles_Scan": ("_parse_args", []),
    "SCS_TransCluster_sweep": ("parse_args", []),
}


@pytest.fixture(autouse=True)
def _factory_default_engine(monkeypatch):
    """Pin what an unset -ge resolves to: no LRGSG_GRAPH_ENGINE, NetworkX."""
    monkeypatch.delenv("LRGSG_GRAPH_ENGINE", raising=False)
    monkeypatch.setattr(_engine, "_default_engine", _engine.GraphEngine("nx"))


def get_parser(name):
    module = importlib.import_module(f"parsers.{name}")
    return getattr(module, "parser", None) or module.build_parser()


def engine_actions(parser):
    return [a for a in parser._actions if "--graph_engine" in a.option_strings]


def parse(name, argv):
    return get_parser(name).parse_args(argv)


def assert_refused(capsys, parse_fn, argv, reason):
    with pytest.raises(SystemExit) as exc:
        parse_fn(argv)
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert f"-ge gt {REFUSED}" in err and reason in err


# ---------------------------------------------------------------------------
# Every parser builds; graph parsers carry -ge once
# ---------------------------------------------------------------------------
def test_every_parser_is_classified():
    assert set(ALL_PARSERS) == GRAPH_PARSERS | OTHER_PARSERS
    assert not GRAPH_PARSERS & OTHER_PARSERS


@pytest.mark.parametrize("name", ALL_PARSERS)
def test_parser_builds(name):
    assert hasattr(get_parser(name), "parse_args")


@pytest.mark.parametrize("name", sorted(GRAPH_PARSERS))
def test_graph_parser_has_one_engine_flag(name):
    (action,) = engine_actions(get_parser(name))
    assert tuple(action.option_strings) == ENGINE_FLAGS
    assert list(action.choices) == list(GRAPH_ENGINE_CHOICES)
    default = (
        DEFAULT_RRG_ENGINE
        if name == "RRG_TransCluster"
        else DEFAULT_GRAPH_ENGINE
    )
    assert action.default == default


@pytest.mark.parametrize("name", sorted(OTHER_PARSERS))
def test_forwarding_parser_has_no_engine_flag(name):
    assert engine_actions(get_parser(name)) == []


def test_shared_engine_spec_is_not_mutated():
    for name in ALL_PARSERS:
        get_parser(name)
    add_graph_engine_arg(GraphArgumentParser(), default="gt")
    (spec,) = graph_engine_opt_args.values()
    assert tuple(graph_engine_opt_args) == (ENGINE_FLAGS,)
    assert spec["default"] == DEFAULT_GRAPH_ENGINE
    assert spec["choices"] == list(GRAPH_ENGINE_CHOICES)


# ---------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", sorted(BASE_ARGV))
def test_nx_is_accepted(name):
    assert parse(name, BASE_ARGV[name] + ["-ge", "nx"]).graph_engine == "nx"


@pytest.mark.parametrize("name", ALWAYS_REFUSED)
def test_gt_is_refused_where_not_ported(name, capsys):
    assert_refused(
        capsys,
        lambda argv: parse(name, argv),
        BASE_ARGV[name] + ["-ge", "gt"],
        "",
    )


@pytest.mark.parametrize("name", ACCEPTED_ON_GT)
def test_gt_is_accepted_with_supported_options(name):
    assert parse(name, BASE_ARGV[name] + ["-ge", "gt"]).graph_engine == "gt"


@pytest.mark.parametrize(
    "name, extra, reason",
    REFUSED_WITH,
    ids=[f"{n}-{' '.join(e)}" for n, e, _ in REFUSED_WITH],
)
def test_gt_is_refused_with_options_it_would_drop(name, extra, reason, capsys):
    argv = BASE_ARGV[name] + extra
    assert parse(name, argv + ["-ge", "nx"]).graph_engine == "nx"
    assert_refused(
        capsys, lambda a: parse(name, a), argv + ["-ge", "gt"], reason
    )


def test_refusal_runs_on_parse_known_args(capsys):
    parser = get_parser("L2D_TransCluster")
    _, extras = parser.parse_known_args(["16", "--unknown"])
    assert extras == ["--unknown"]
    assert_refused(capsys, parser.parse_known_args, ["16", "-ge", "gt"], "")


# ---------------------------------------------------------------------------
# Programs with their own parser
# ---------------------------------------------------------------------------
def inline_parse(name):
    fn_name, base = INLINE[name]
    fn = getattr(importlib.import_module(name), fn_name)
    return lambda argv: fn(base + argv)


@pytest.mark.parametrize("name", sorted(INLINE))
def test_inline_program_has_engine_flag(name):
    parse_fn = inline_parse(name)
    default = DEFAULT_GRAPH_ENGINE if name == "SCS_TransCluster_sweep" else None
    assert parse_fn([]).graph_engine == default
    assert parse_fn(["-ge", "nx"]).graph_engine == "nx"


@pytest.mark.parametrize(
    "name, argv",
    [
        ("L2D_SignedRW_Bubbles", []),
        ("L2D_SignedRW_BubblesRand", ["-wd", ""]),
        ("L2D_SignedRW_TimeToTouch", ["-wd", ""]),
    ],
)
def test_inline_gt_is_accepted_without_a_workdir(name, argv):
    parse_fn = inline_parse(name)
    assert parse_fn(argv + ["-ge", "gt"]).graph_engine == "gt"


@pytest.mark.parametrize(
    "name, argv, reason",
    [
        ("L2D_SignedRW_Bubbles", ["-wd", "elsewhere"], "-wd"),
        ("L2D_SignedRW_BubblesRand", [], "-wd"),
        ("L2D_SignedRW_TimeToTouch", [], "-wd"),
        ("L2D_SignedRW_Bubbles_Scan", [], "--skip-existing"),
        ("SCS_TransCluster_sweep", [], "SCS kernels"),
    ],
)
def test_inline_gt_is_refused(name, argv, reason, capsys):
    parse_fn = inline_parse(name)
    assert parse_fn(argv + ["-ge", "nx"]).graph_engine == "nx"
    assert_refused(capsys, parse_fn, argv + ["-ge", "gt"], reason)


def test_unset_engine_follows_the_environment(monkeypatch, capsys):
    parse_fn = inline_parse("L2D_SignedRW_Bubbles")
    argv = ["-wd", "elsewhere"]
    assert parse_fn(argv).graph_engine is None
    monkeypatch.setenv("LRGSG_GRAPH_ENGINE", "gt")
    assert_refused(capsys, parse_fn, argv, "-wd")
    assert parse_fn(argv + ["-ge", "nx"]).graph_engine == "nx"
