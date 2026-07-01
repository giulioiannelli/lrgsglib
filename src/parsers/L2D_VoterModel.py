from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import parse_arguments as _parse_arguments

optionalaction_args_dict = {
    **L2D_opt_args,
    **L2D_VOTER_opt_args,
    **L2D_VOTER_action_args,
}

parser = argparse.ArgumentParser(
    description=L2D_VOTER_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

for k, v in L2D_VOTER_args.items():
    if isinstance(k, tuple):
        parser.add_argument(*k, **v)
    else:
        parser.add_argument(k, **v)

for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)


def _validate_voter_args(args: argparse.Namespace) -> argparse.Namespace:
    """Fail fast on ill-defined (rule, upd_mode) combinations.

    ``link`` (edge-update) and ``gillespie`` (CTMC) are intrinsically copy
    operations, so they are only defined for ``rule='linear'`` (the VoterModel
    constructor enforces this too; here we surface a clean CLI error). Backend
    constraints (np/cu require synchronous and cannot save the trajectory or
    track clusters) are enforced authoritatively by the solver at run time.
    """
    if args.upd_mode in ("link", "gillespie") and args.rule != "linear":
        parser.error(
            f"upd_mode='{args.upd_mode}' requires rule='linear' "
            f"(got rule='{args.rule}')."
        )
    return args


def parse_arguments(
    parser: argparse.ArgumentParser, argv: list[str] | None = None
) -> argparse.Namespace:
    args = _parse_arguments(parser) if argv is None else parser.parse_args(args=argv)
    return _validate_voter_args(args)
