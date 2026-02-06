"""
ConfigurationModelNX: NetworkX-based configuration model random graph.

This module provides the ConfigurationModelNX class for generating signed
random graphs with prescribed degree sequences.
"""

from .ConfigurationModelNX import ConfigurationModelNX

# Backward compatibility alias
ConfigurationModel = ConfigurationModelNX

__all__ = [
    "ConfigurationModelNX",
    "ConfigurationModel",  # alias for backward compatibility
]
