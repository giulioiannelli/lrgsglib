from .Lattice2D import L2D_PATH_DICT, L2D_SHRT_GEO_DICT
from typing import Dict
#
def get_lattice_path(
    geo: str,
    path_dict: Dict[str, str] = L2D_PATH_DICT,
    alias_dict: Dict[str, str] = L2D_SHRT_GEO_DICT
) -> str:
    """
    Resolve a geometry key or its abbreviation to the corresponding lattice path.

    Parameters
    ----------
    geo : str
        Geometry name or its abbreviation.
    path_dict : Dict[str, str], optional
        Mapping from full geometry names to path strings.
    alias_dict : Dict[str, str], optional
        Mapping from abbreviations to full geometry names.

    Returns
    -------
    str
        The path string corresponding to `geo`.

    Raises
    ------
    KeyError
        If `geo` is neither a key in `path_dict` nor an abbreviation in `alias_dict`,
        or if the resolved full name is not in `path_dict`.
    """
    # direct lookup
    if geo in path_dict:
        return path_dict[geo]

    # try as abbreviation
    full_geo = alias_dict.get(geo)
    if full_geo is None:
        raise KeyError(f"Unknown geometry or abbreviation: '{geo}'")

    # lookup resolved name
    try:
        return path_dict[full_geo]
    except KeyError:
        raise KeyError(f"No path defined for geometry '{full_geo}' (alias '{geo}')")