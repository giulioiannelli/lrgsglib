"""
Real-world dataset loading module.

This module provides tools for loading real-world networks from various sources:
- KONECT (Koblenz Network Collection)
- SNAP (Stanford Network Analysis Project)
- Network Repository
- Local files (edge lists, GraphML, GML, etc.)

Example
-------
>>> from lrgsglib.nx_patches.datasets import RealDatasetLoader
>>> loader = RealDatasetLoader(cache_dir="~/.lrgsglib/datasets")
>>> karate = loader.load("karate", source="local")
>>> print(f"Nodes: {karate.N}, Edges: {karate.M}")
"""

from .real_dataset_loader import (
    RealDatasetLoader,
    DatasetInfo,
    BUILTIN_DATASETS,
)

__all__ = [
    "RealDatasetLoader",
    "DatasetInfo",
    "BUILTIN_DATASETS",
]
