"""
inverse_msmd - Inverse MSMD (Molecular Simulation with Multiple Descriptors) Package

This package provides utilities for protein structure superimposition
and matching score calculations for inverse MSMD analysis.

Main modules:
- utils.bio_utils: BioPython utilities for structure manipulation
- utils.spatial_utils: Spatial calculation utilities
- utils.path_utils: Path handling utilities
"""

from .utils.bio_utils import SuperImposer, PDB
from .utils.spatial_utils import estimate_volume
from .utils.path_utils import expandpath

__version__ = "0.1.0"
__author__ = "Keisuke Yanagisawa"
__email__ = "yanagisawa@comp.isct.ac.jp"

__all__ = [
    "SuperImposer",
    "PDB",
    "estimate_volume",
    "expandpath",
]