"""
PBM Equation Shell - Population Balance Model for Flocculation

This package implements the discretized Population Balance Equation (PBE)
for modeling particle flocculation in Dissolved Air Flotation (DAF) systems.

Based on:
    Abreu, D. A. d. S., et al. (2021). "Integration of first principle models
    and machine learning in a modeling framework: An application to flocculation"
    Chemical Engineering Science.
"""

from .population_balance_model import (
    PopulationBalanceModel,
    construct_beta_matrix
)

__version__ = "1.0.0"
__author__ = "DAF-Sim Development Team"

__all__ = [
    "PopulationBalanceModel",
    "construct_beta_matrix"
]
