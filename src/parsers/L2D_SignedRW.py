"""argparse wiring for ``L2D_SignedRW.py``.

Mirrors ``parsers/L2D_ContactProcess.py``: builds an ArgumentParser
from the progargs dicts and exposes ``parse_arguments(parser)`` with
post-validation.
"""

import argparse

from lrgsglib.core import *  # noqa: F401,F403
from lrgsglib.config.progargs import *  # noqa: F401,F403
from parsers.shared import parse_arguments as _parse_arguments


optionalaction_args_dict = {
    **L2D_SRW_opt_args,  # noqa: F405
    **L2D_SRW_action_args,  # noqa: F405
}

parser = argparse.ArgumentParser(
    description=L2D_SRW_description,  # noqa: F405
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

for k, v in L2D_SRW_args.items():  # noqa: F405
    if isinstance(k, tuple):
        parser.add_argument(*k, **v)
    else:
        parser.add_argument(k, **v)

for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)


def _validate_srw_args(args: argparse.Namespace) -> argparse.Namespace:
    if args.start_a == 'fixed' and args.start_a_node is None:
        parser.error("--start_a=fixed requires --start_a_node <int>.")
    if args.mode == 'pair' and args.start_b == 'fixed' and args.start_b_node is None:
        parser.error("--start_b=fixed (pair mode) requires --start_b_node <int>.")
    if args.runlang.lower() != 'py':
        parser.error(
            "Only --runlang=py is available in Phase 1 of SignedRW; "
            f"got {args.runlang!r}."
        )
    return args


def parse_arguments(
    parser: argparse.ArgumentParser,
    argv: list[str] | None = None,
) -> argparse.Namespace:
    args = _parse_arguments(parser) if argv is None else parser.parse_args(args=argv)
    return _validate_srw_args(args)
