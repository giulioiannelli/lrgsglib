from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import parse_arguments

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
    parser.add_argument(k, **v)

for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
