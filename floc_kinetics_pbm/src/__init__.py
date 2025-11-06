"""
Floc Kinetics PBM (Population Balance Model)

This module implements the Population Balance Model for floc aggregation and
breakage kinetics based on standard formulations from the flocculation modeling
literature, including Smoluchowski theory and Saffman-Turner turbulent
aggregation kernels.

Core Components:
- PopulationBalanceModel: Main PBM solver class
- FlocKernels: Aggregation and breakage kernel calculations
- FlocProperties: Physical properties of flocs and water
"""

from .kernels import FlocKernels, calculate_floc_kernels
from .pbe_solver import PopulationBalanceModel
from .properties import FlocProperties

__all__ = [
    'FlocKernels',
    'calculate_floc_kernels',
    'PopulationBalanceModel',
    'FlocProperties'
]

__version__ = '1.0.0'
