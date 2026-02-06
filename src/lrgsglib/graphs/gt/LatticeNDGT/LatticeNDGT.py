"""
Re-export LatticeND as LatticeNDGT for unified naming convention.

Note: LatticeND already exists natively in gt_patches.
"""

from ....gt_patches import LatticeND

LatticeNDGT = LatticeND

__all__ = ["LatticeNDGT"]
