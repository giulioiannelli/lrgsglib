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
# Graph engine: the -ge flag, class resolution, refusals
# ---------------------------------------------------------------------------
# progargs imports this module while it loads, so the functions below import
# the progargs names they need at call time.

# The only cell the graph-tool branches of the prepare kernels realize
_GT_CELL = "rand"


def add_graph_engine_arg(parser, default=None):
    """Add ``-ge/--graph_engine`` to a parser built outside progargs.

    ``default=None`` leaves the engine to the graph factories: the
    ``LRGSG_GRAPH_ENGINE`` variable, else NetworkX.
    """
    from lrgsglib.config.progargs.common import graph_engine_opt_args
    from lrgsglib.config.progargs.phelp.generic import phelp_graph_engine_env

    for flags, spec in graph_engine_opt_args.items():
        help_ = phelp_graph_engine_env if default is None else spec["help"]
        parser.add_argument(*flags, **{**spec, "default": default, "help": help_})


def effective_engine(engine):
    """``engine``, or the graph factories' default engine when it is None."""
    if engine is not None:
        return engine
    from lrgsglib.graphs._engine import get_default_engine_from_env

    return get_default_engine_from_env().value


def get_graph_engine(args) -> str:
    """The parsed ``-ge`` value; a program without the flag gets the default."""
    from lrgsglib.config.progargs.defs.generic import DEFAULT_GRAPH_ENGINE

    return getattr(args, "graph_engine", DEFAULT_GRAPH_ENGINE)


def graph_engine_suffix(args) -> str:
    """``_ge=<engine>`` for a non-default engine, else ``''``.

    For output names whose directory does not depend on the engine, so runs
    on different engines never overwrite each other's files.
    """
    from lrgsglib.config.progargs.defs.generic import DEFAULT_GRAPH_ENGINE

    engine = get_graph_engine(args)
    return "" if engine == DEFAULT_GRAPH_ENGINE else f"_ge={engine}"


def resolve_graph_class(graph_type: str, engine: str):
    """The class registered for ``graph_type`` on ``engine``, never another engine's.

    Raises
    ------
    ValueError
        Unknown engine or graph type, or a ``gt`` request that the registry
        answers with a NetworkX class (an NX alias under ``graphs/gt``).
    """
    import lrgsglib.graphs  # noqa: F401  (importing the facades registers the classes)
    from lrgsglib.config.progargs.defs.generic import GRAPH_ENGINE_CHOICES
    from lrgsglib.graphs._engine import get_implementation

    if not isinstance(engine, str) or engine.lower() not in GRAPH_ENGINE_CHOICES:
        raise ValueError(
            f"unknown graph engine {engine!r}; choose from {GRAPH_ENGINE_CHOICES}"
        )
    engine = engine.lower()
    cls = get_implementation(graph_type, engine, fallback=False)
    if engine == "gt":
        from lrgsglib.graphs.gt.SignedGraphGT.SignedGraphGT import SignedGraphGT

        if not issubclass(cls, SignedGraphGT):
            raise ValueError(
                f"{graph_type} has no graph-tool class: gt resolves to {cls.__name__}"
            )
    return cls


class GraphArgumentParser(argparse.ArgumentParser):
    """ArgumentParser that refuses the graph engines a program cannot run on yet.

    ``refuse_engine`` registers a refusal. Every parse (``parse_args``,
    ``parse_known_args``) then ends in ``self.error`` (exit 2) when the
    effective engine is the refused one and ``when(args)`` holds, so the
    program stops before it builds anything.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.engine_refusals = []

    def refuse_engine(self, engine, reason, when=None):
        """Refuse ``engine``, only when ``when(args)`` holds if ``when`` is given."""
        self.engine_refusals.append((engine, reason, when))

    def parse_known_args(self, args=None, namespace=None):
        parsed, extras = super().parse_known_args(args, namespace)
        if not self.engine_refusals:
            return parsed, extras
        engine = effective_engine(getattr(parsed, "graph_engine", None))
        for refused, reason, when in self.engine_refusals:
            if engine == refused and (when is None or when(parsed)):
                self.error(f"-ge {engine} is not supported here yet: {reason}")
        return parsed, extras


def refuse_gt_ignored_options(parser):
    """Refuse ``-ge gt`` together with the options that graph-tool graphs drop
    without a word: a cell other than ``rand`` (the only one the prepare
    kernels realize on graph-tool) and any workdir (graph-tool graphs keep
    their own data path)."""
    dests = {action.dest for action in parser._actions}
    if "cell_type" in dests:
        parser.refuse_engine(
            "gt",
            f"graph-tool graphs realize only -c {_GT_CELL}",
            when=lambda a: a.cell_type != _GT_CELL,
        )
    if "workdir" in dests:
        parser.refuse_engine(
            "gt",
            "graph-tool graphs ignore -wd and keep their own data path (pass -wd '')",
            when=lambda a: bool(a.workdir),
        )
