"""Shared utilities for transient cluster simulation kernels."""

from __future__ import annotations

import glob
import os
import pickle as pk
import re
from collections import Counter
from typing import Any, Callable, Iterable, Protocol

import numpy as np

from lrgsglib.config.const import PKL, TXT
from lrgsglib.utils.basic.strings import get_first_int_in_str

__all__ = [
    "build_geometry_selector",
    "compose_output_path",
    "load_pcluster_state",
    "resolve_mode_extension",
]


class SupportsNetworkDict(Protocol):
    """Protocol capturing the minimal lattice interface required here."""

    nwDict: Any



def resolve_mode_extension(mode: str) -> str:
    """Return the output file extension associated with a computation mode."""

    match mode:
        case "pCluster":
            return PKL
        case "ordParam":
            return TXT
        case _:
            raise ValueError(f"Invalid mode specified: {mode}")


def build_geometry_selector(cell: str) -> Callable[[SupportsNetworkDict], Iterable]:
    """Return the geometry selection function based on the defect specification."""

    if cell in {"rand", "randZERR", "randXERR"}:

        def geometry_func(lattice: SupportsNetworkDict) -> Iterable[int]:
            return lattice.nwDict[cell]["G"]

        return geometry_func

    if cell.startswith("ball"):
        radius = get_first_int_in_str(cell)
        if radius is None:
            raise ValueError(f"Unable to infer radius from cell specification: {cell}")

        def geometry_func(lattice: SupportsNetworkDict) -> Iterable[int]:
            getter = getattr(lattice.nwDict, "get_links_rball", None)
            if getter is None:
                raise AttributeError("Lattice network dictionary does not support radial balls")
            return getter(radius)

        return geometry_func

    raise ValueError(f"Invalid cell specified: {cell}")


def compose_output_path(base_path: str, *tokens: str, ext: str) -> str:
    """Compose an output file path using non-empty tokens and the given extension."""

    filtered_tokens = [str(token) for token in tokens if str(token)]
    filename = "_".join(filtered_tokens)
    return os.path.join(base_path, filename) + ext


def load_pcluster_state(base_pattern: str, has_suffix: bool) -> tuple[Counter, int, str | None]:
    """Load the latest persisted cluster distribution if present."""

    matches = glob.glob(f"{base_pattern}*")
    if not matches:
        return Counter(), 0, None

    fname_exists = max(matches, key=os.path.getmtime)
    with open(fname_exists, "rb") as handle:
        merged_dict = pk.load(handle)

    tokens = os.path.basename(fname_exists).split("_")
    avg_idx = -2 if has_suffix else -1
    try:
        avg_token = os.path.splitext(tokens[avg_idx])[0]
    except IndexError as exc:  # pragma: no cover - defensive fallback
        raise ValueError(f"Could not parse averages from filename: {fname_exists}") from exc

    match = next((int(num) for num in re.findall(r"\d+", avg_token)), None)
    if match is None:
        return Counter(merged_dict), 0, fname_exists

    return Counter(merged_dict), match, fname_exists

