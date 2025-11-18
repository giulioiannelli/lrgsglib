from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import parse_arguments as _parse_arguments

optionalaction_args_dict = {
    **L2D_opt_args,
    **L2D_CPROC_opt_args,
    **L2D_CPROC_action_args,
}

parser = argparse.ArgumentParser(
    description=L2D_CPROC_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

for k, v in L2D_CPROC_args.items():
    if isinstance(k, tuple):
        parser.add_argument(*k, **v)
    else:
        parser.add_argument(k, **v)

for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)


def _validate_contact_args(args: argparse.Namespace) -> argparse.Namespace:
    dynamics = args.dynamics.upper()
    runlang = args.runlang.upper()

    if dynamics == "EI":
        if args.gamma is None:
            parser.error("Excitation/inhibition dynamics require --gamma.")
        if not runlang.startswith("C1"):
            parser.error("Excitation/inhibition dynamics currently support only C1* runlangs (e.g., --runlang C1c).")
    elif dynamics == "SIR":
        if args.mu is None:
            parser.error("SIR dynamics require --mu.")
        if runlang.startswith("C1"):
            parser.error("SIR dynamics do not map to the C1* contact-process kernels.")
    else:
        parser.error(f"Unsupported dynamics selection: {dynamics}")

    return args


def parse_arguments(parser: argparse.ArgumentParser, argv: list[str] | None = None) -> argparse.Namespace:
    args = _parse_arguments(parser) if argv is None else parser.parse_args(args=argv)
    return _validate_contact_args(args)
