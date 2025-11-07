"""
Eulerian Momentum Equation Module for DAF Simulation
====================================================

This module implements the Eulerian phase momentum equation as specified in
VTT Technical Report (2023) - JMSE 11(4):820.

Main Components:
- Phase: Represents a single phase (liquid, gas, solid)
- MomentumEquation: Computes momentum equation terms
- InterfacialForces: Handles interfacial momentum transfer

Author: DAF-Sim Development Team
Date: 2025-11-06
"""

from .phase import Phase
from .momentum_equation import MomentumEquation
from .interfacial_forces import InterfacialForces

__version__ = "1.0.0"
__all__ = ["Phase", "MomentumEquation", "InterfacialForces"]
