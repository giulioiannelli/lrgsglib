"""Anchor policy implementations for GraphOfGraphs.

Anchor policies determine which node in each fiber graph is connected
to its corresponding base node.
"""

from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING, Callable, Union

import numpy as np

if TYPE_CHECKING:
    from numpy.random import Generator


class AnchorPolicy(Enum):
    """Predefined anchor policies for fiber-base connections."""

    FIRST = "first"
    CENTER = "center"
    LAST = "last"
    RANDOM = "random"


def get_anchor_index(
    policy: Union[str, AnchorPolicy, Callable[[int, int], int]],
    base_idx: int,
    n_fiber: int,
    rng: "Generator",
) -> int:
    """Determine the anchor node index for a fiber.

    Parameters
    ----------
    policy : str, AnchorPolicy, or callable
        Anchor selection strategy:
        - 'first': Always use node 0
        - 'center': Use middle node (n_fiber // 2)
        - 'last': Use last node (n_fiber - 1)
        - 'random': Random node (seeded)
        - callable: Custom function(base_idx, n_fiber) -> anchor_idx
    base_idx : int
        Index of the base node this fiber attaches to
    n_fiber : int
        Number of nodes in the fiber
    rng : numpy.random.Generator
        Random number generator for 'random' policy

    Returns
    -------
    int
        Index of the anchor node within the fiber (0 to n_fiber-1)
    """
    if callable(policy) and not isinstance(policy, (str, AnchorPolicy)):
        return policy(base_idx, n_fiber)

    if isinstance(policy, str):
        policy = AnchorPolicy(policy.lower())

    if policy == AnchorPolicy.FIRST:
        return 0
    elif policy == AnchorPolicy.CENTER:
        return n_fiber // 2
    elif policy == AnchorPolicy.LAST:
        return n_fiber - 1
    elif policy == AnchorPolicy.RANDOM:
        return int(rng.integers(0, n_fiber))
    else:
        raise ValueError(f"Unknown anchor policy: {policy}")


def resolve_anchor_indices(
    policy: Union[str, AnchorPolicy, Callable[[int, int], int]],
    n_base: int,
    n_fiber: int,
    rng: "Generator",
) -> list[int]:
    """Compute anchor indices for all fibers.

    Parameters
    ----------
    policy : str, AnchorPolicy, or callable
        Anchor selection strategy
    n_base : int
        Number of base nodes (= number of fibers)
    n_fiber : int
        Number of nodes in each fiber
    rng : numpy.random.Generator
        Random number generator

    Returns
    -------
    list[int]
        Anchor index for each fiber
    """
    return [
        get_anchor_index(policy, base_idx, n_fiber, rng)
        for base_idx in range(n_base)
    ]
