"""
Pillar 3: Core Physics Model for DAF Simulator
Based on Han et al. (2001) trajectory analysis with DLVO theory

This module implements the collision efficiency factor (alpha_bp) calculation
using trajectory analysis that accounts for both hydrodynamic and interparticle forces.
"""

from .collision_efficiency import calculate_attachment_efficiency
from .particle import Particle
from .bubble import Bubble
from .fluid_properties import FluidProperties

__version__ = "1.0.0"
__all__ = [
    "calculate_attachment_efficiency",
    "Particle",
    "Bubble",
    "FluidProperties",
]
