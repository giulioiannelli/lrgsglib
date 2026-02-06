"""
ConfigurationModelGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.random import ConfigurationModelNX

warnings.warn(
    "ConfigurationModelGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

ConfigurationModelGT = ConfigurationModelNX

__all__ = ["ConfigurationModelGT"]
