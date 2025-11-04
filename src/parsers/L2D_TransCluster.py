from __future__ import annotations

import argparse

from lrgsglib.proglib import *
from parsers.shared import parse_arguments as _parse_arguments

_optional_args = {**L2D_opt_args}

for key in list(_optional_args):
    if '--number_of_averages' in key:
        number_opts = {**_optional_args.pop(key)}
        number_opts['default'] = DEFAULT_L2D_TRANSCLUSTER_NAVG
        _optional_args[('-n', '-na', '--number_of_averages')] = number_opts
        break

for key in list(_optional_args):
    if '--workdir' in key:
        workdir_opts = {**_optional_args.pop(key)}
        _optional_args[('-wd', '--workdir', '--outdir')] = workdir_opts
        break

optional_args_dict = {
    **_optional_args,
    **L2D_TransCluster_optional_args_dict,
}

parser = argparse.ArgumentParser(
    description=L2D_TransCluster_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

for key, value in L2D_TransCluster_args.items():
    parser.add_argument(key, **value)

for key, value in optional_args_dict.items():
    parser.add_argument(*key, **value)

for key, value in L2D_TransCluster_action_args_dict.items():
    parser.add_argument(*key, **value)


def parse_arguments(parser: argparse.ArgumentParser) -> argparse.Namespace:
    args = _parse_arguments(parser)
    if not hasattr(args, 'outdir'):
        setattr(args, 'outdir', getattr(args, 'workdir', ''))
    elif not getattr(args, 'workdir', None):
        setattr(args, 'workdir', args.outdir)
    return args
