"""parsers.shared engine helpers: the graph class a program gets for ``-ge``."""

import sys
from argparse import Namespace
from pathlib import Path

import pytest

pytest.importorskip("graph_tool")

SRC = Path(__file__).resolve().parents[2] / "src"
sys.path.insert(0, str(SRC))

from lrgsglib.config.progargs.defs.generic import (  # noqa: E402
    DEFAULT_GRAPH_ENGINE,
)
from lrgsglib.graphs.gt.SignedGraphGT.SignedGraphGT import (  # noqa: E402
    SignedGraphGT,
)
from lrgsglib.graphs.nx import SignedGraphNX  # noqa: E402
from parsers.shared import (  # noqa: E402
    get_graph_engine,
    graph_engine_suffix,
    resolve_graph_class,
)

pytestmark = pytest.mark.integration

# The graph types the program kernels resolve by engine
PROGRAM_GRAPH_TYPES = [
    "ErdosRenyi",
    "Lattice2D",
    "Lattice3D",
    "MultiplicativeCascade",
]


@pytest.mark.parametrize("graph_type", PROGRAM_GRAPH_TYPES)
def test_gt_resolves_to_a_graph_tool_class(graph_type):
    cls = resolve_graph_class(graph_type, "gt")
    assert issubclass(cls, SignedGraphGT)
    assert resolve_graph_class(graph_type, "GT") is cls


@pytest.mark.parametrize("graph_type", PROGRAM_GRAPH_TYPES)
def test_nx_resolves_to_a_networkx_class(graph_type):
    cls = resolve_graph_class(graph_type, "nx")
    assert issubclass(cls, SignedGraphNX)
    assert resolve_graph_class(graph_type, "NX") is cls


@pytest.mark.parametrize("engine", [None, "", "ig", "networkx", 1])
def test_unknown_engine_raises(engine):
    with pytest.raises(ValueError, match="unknown graph engine"):
        resolve_graph_class("Lattice2D", engine)


def test_unknown_graph_type_raises():
    with pytest.raises(ValueError, match="Unknown graph type"):
        resolve_graph_class("NoSuchGraph", "nx")


def test_gt_request_never_gets_a_networkx_alias():
    # graphs/gt registers SCSGeneralizedNN as the NetworkX class
    with pytest.raises(ValueError, match="no graph-tool class"):
        resolve_graph_class("SCSGeneralizedNN", "gt")


def test_engine_of_parsed_args():
    assert get_graph_engine(Namespace()) == DEFAULT_GRAPH_ENGINE
    assert get_graph_engine(Namespace(graph_engine="gt")) == "gt"


def test_engine_suffix_names_only_a_non_default_engine():
    assert graph_engine_suffix(Namespace()) == ""
    assert (
        graph_engine_suffix(Namespace(graph_engine=DEFAULT_GRAPH_ENGINE)) == ""
    )
    assert graph_engine_suffix(Namespace(graph_engine="gt")) == "_ge=gt"
