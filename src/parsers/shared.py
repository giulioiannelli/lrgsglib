import argparse
import ast
import sys

import numpy as np


class CustomHelpAction(argparse.Action):
    """Custom help action that also shows help for the inner program."""

    def __init__(
        self,
        option_strings,
        dest=argparse.SUPPRESS,
        default=argparse.SUPPRESS,
        help=None,
        inner_prog_name="",
        inner_positional_args=None,
        inner_optional_args=None,
        **kwargs,
    ):
        self.inner_prog_name = inner_prog_name
        self.inner_positional_args = inner_positional_args or {}
        self.inner_optional_args = inner_optional_args or {}
        kwargs.setdefault("nargs", 0)
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            help=help,
            **kwargs,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        parser.print_help()
        if self.inner_prog_name:
            print(f"\n{self.inner_prog_name}.py options:")
            print(
                f"  Additional options passed directly to {self.inner_prog_name}.py are shown below.\n"
            )

        inner_parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            add_help=False,
        )

        for key, opts in self.inner_positional_args.items():
            inner_parser.add_argument(key, **opts)

        for key, opts in self.inner_optional_args.items():
            inner_parser.add_argument(*key, **opts)

        inner_parser.print_help()
        sys.exit(0)

def parse_arguments(parser: argparse.ArgumentParser) -> argparse.Namespace:
    args = parser.parse_args()
    return args

def parsDict_get(parsDict, var, key):
    return parsDict[var].get(key, None)

def parse_linspace(value):
    """Parses a tuple string and returns a linspace array."""
    try:
        start, stop, num = ast.literal_eval(value)  # Safe parsing of tuples
        return np.linspace(float(start), float(stop), int(num))
    except (ValueError, SyntaxError, TypeError):
        raise argparse.ArgumentTypeError("Invalid linspace format. Use (start, stop, num)")

def parse_multiple_linspace(value):
    """Parses multiple linspace tuples separated by ';' and concatenates them."""
    segments = value.split(';')  # Split multiple tuples
    arrays = [parse_linspace(seg) for seg in segments]  # Convert each to linspace
    return np.concatenate(arrays)  # Merge all into one array


# ---------------------------------------------------------------------------
# Graph engine resolution
# ---------------------------------------------------------------------------
_VALID_GRAPH_ENGINES = {"nx", "gt"}

# Mapping: (graph_type, engine) -> import path
_GRAPH_CLASS_MAP = {
    ("Lattice2D", "nx"):        "lrgsglib.graphs.nx.Lattice2DNX.Lattice2DNX",
    ("Lattice2D", "gt"):        "lrgsglib.graphs.gt.Lattice2DGT.Lattice2DGT",
    ("Lattice3D", "nx"):        "lrgsglib.graphs.nx.Lattice3DNX.Lattice3DNX",
    ("Lattice3D", "gt"):        "lrgsglib.graphs.gt.Lattice3DGT.Lattice3DGT",
    ("ErdosRenyi", "nx"):       "lrgsglib.graphs.nx.ErdosRenyiNX.ErdosRenyiNX",
    ("ErdosRenyi", "gt"):       "lrgsglib.graphs.gt.ErdosRenyiGT.ErdosRenyiGT",
    ("BarabasiAlbert", "nx"):   "lrgsglib.graphs.nx.BarabasiAlbertNX.BarabasiAlbertNX",
    ("BarabasiAlbert", "gt"):   "lrgsglib.graphs.gt.BarabasiAlbertGT.BarabasiAlbertGT",
    ("WattsStrogatz", "nx"):    "lrgsglib.graphs.nx.WattsStrogatzNX.WattsStrogatzNX",
    ("WattsStrogatz", "gt"):    "lrgsglib.graphs.gt.WattsStrogatzGT.WattsStrogatzGT",
    ("MultiplicativeCascade", "nx"): "lrgsglib.graphs.nx.MultiplicativeCascadeNX.MultiplicativeCascadeNX.MultiplicativeCascadeGraphNX",
    ("MultiplicativeCascade", "gt"): "lrgsglib.graphs.gt.MultiplicativeCascadeGT.MultiplicativeCascadeGT.MultiplicativeCascadeGraphGT",
}


def resolve_graph_class(graph_type: str, engine: str = "nx"):
    """Resolve a graph class by type and engine.

    Parameters
    ----------
    graph_type : str
        Graph type name (e.g. ``"Lattice2D"``, ``"ErdosRenyi"``).
    engine : str
        Graph engine: ``"nx"`` (NetworkX) or ``"gt"`` (graph-tool).

    Returns
    -------
    type
        The resolved graph class.

    Raises
    ------
    ValueError
        If engine is invalid or the graph type has no GT implementation.
    """
    engine = engine.lower()
    if engine not in _VALID_GRAPH_ENGINES:
        raise ValueError(
            f"Unknown graph engine '{engine}'. Valid: {_VALID_GRAPH_ENGINES}"
        )

    key = (graph_type, engine)
    if key not in _GRAPH_CLASS_MAP:
        available = [k for k, e in _GRAPH_CLASS_MAP if k == graph_type]
        raise ValueError(
            f"No {engine} implementation for '{graph_type}'. "
            f"Available engines: {available}"
        )

    import importlib
    module_path, class_name = _GRAPH_CLASS_MAP[key].rsplit(".", 1)
    mod = importlib.import_module(module_path)
    return getattr(mod, class_name)


def get_graph_engine(args) -> str:
    """Extract graph engine from parsed args, defaulting to 'nx'."""
    return getattr(args, "graph_engine", "nx")
