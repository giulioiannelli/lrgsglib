from typing import Any
from lrgsglib.nx_patches import ErdosRenyi

__all__ = [
    "initialize_er_dict_args",
    "prepare_erdos_renyi",
]

def initialize_er_dict_args(args: Any) -> dict[str, Any]:
    """
    Generate a dictionary of arguments for initializing an ErdosRenyi object.

    Parameters
    ----------
    args : Any
        Argument object containing Erdos-Renyi configuration with attributes
        N (number of nodes), p (edge probability), pflip (edge flipping probability),
        and workdir (working directory).

    Returns
    -------
    dict[str, Any]
        Dictionary of arguments for ErdosRenyi initialization, including
        init_nw_dict=True for network dictionary initialization.

    Examples
    --------
    >>> class Args:
    ...     N = 100
    ...     p = 0.1
    ...     pflip = 0.2
    ...     workdir = '/path/to/work'
    >>> args = Args()
    >>> er_kwargs = initialize_er_dict_args(args)
    >>> er_kwargs['n']
    100
    """
    return dict(
        n=args.N,
        p=args.p,
        sgpathn=args.workdir,
        pflip=args.pflip,
        init_nw_dict=True
    )


def prepare_erdos_renyi(args: Any, cell_type: str | None = None, **kwargs: Any) -> ErdosRenyi:
    """
    Initialize an ErdosRenyi graph and flip edges according to cell_type strategy.

    Parameters
    ----------
    args : Any
        Argument object containing Erdos-Renyi configuration.
    cell_type : str | None, optional
        Cell type pattern for edge flipping. If None, uses args.cell_type.
        Options include 'rand', 'randXERR', etc.
    **kwargs : Any
        Additional keyword arguments passed to ErdosRenyi constructor.

    Returns
    -------
    ErdosRenyi
        Configured ErdosRenyi instance with edges flipped according to cell_type.

    Examples
    --------
    >>> class Args:
    ...     N = 100
    ...     p = 0.1
    ...     pflip = 0.2
    ...     workdir = ''
    ...     cell_type = 'rand'
    >>> args = Args()
    >>> er = prepare_erdos_renyi(args)
    >>> er.syshape > 0
    True
    """
    if cell_type is None:
        cell_type = args.cell_type

    er = ErdosRenyi(**initialize_er_dict_args(args), **kwargs)

    # Flip edges according to cell_type strategy
    if er.init_nw_dict and cell_type in er.nwDict:
        er.flip_sel_edges(er.nwDict[cell_type]['G'])
    else:
        # Fallback to random flipping if pattern not available
        er.flip_random_fract_edges()

    return er
