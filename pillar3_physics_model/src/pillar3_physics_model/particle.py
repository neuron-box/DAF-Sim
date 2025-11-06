"""
Particle class for DAF simulation
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class Particle:
    """
    Represents a particle in the DAF system.

    Attributes:
        diameter: Particle diameter [m]
        zeta_potential: Zeta potential [V]
        density: Particle density [kg/m³]
        hamaker_constant: Hamaker constant for particle-water-bubble system [J]
                         Default is ~5e-21 J (typical for hydrophobic particles in water)
    """
    diameter: float
    zeta_potential: float
    density: float = 2650.0  # Default: silica particles
    hamaker_constant: Optional[float] = None  # Will use default if None

    def __post_init__(self):
        """Validate particle properties and set defaults."""
        if self.diameter <= 0:
            raise ValueError("Particle diameter must be positive")
        if abs(self.zeta_potential) > 0.3:
            raise ValueError("Zeta potential magnitude should typically be < 0.3 V")
        if self.density <= 0:
            raise ValueError("Particle density must be positive")

        # Set default Hamaker constant if not provided
        # Typical value for hydrophobic particles in water: ~5e-21 J
        if self.hamaker_constant is None:
            self.hamaker_constant = 5.0e-21

    @property
    def radius(self) -> float:
        """Particle radius [m]."""
        return self.diameter / 2.0

    @property
    def volume(self) -> float:
        """Particle volume [m³]."""
        return (4.0 / 3.0) * 3.14159265359 * self.radius**3

    @property
    def mass(self) -> float:
        """Particle mass [kg]."""
        return self.density * self.volume
