"""
Bubble class for DAF simulation
"""

from dataclasses import dataclass


@dataclass
class Bubble:
    """
    Represents a microbubble in the DAF system.

    Attributes:
        diameter: Bubble diameter [m]
        zeta_potential: Zeta potential [V]
                       Typically negative for air-water interface
        rise_velocity: Bubble rise velocity [m/s]
                      If None, will be calculated from Stokes law
    """
    diameter: float
    zeta_potential: float
    rise_velocity: float = None

    def __post_init__(self):
        """Validate bubble properties."""
        if self.diameter <= 0:
            raise ValueError("Bubble diameter must be positive")
        if abs(self.zeta_potential) > 0.3:
            raise ValueError("Zeta potential magnitude should typically be < 0.3 V")

        # Typical DAF bubbles are 10-100 micrometers
        if self.diameter > 500e-6:
            import warnings
            warnings.warn(
                f"Bubble diameter {self.diameter*1e6:.1f} μm is larger than typical "
                "DAF bubbles (10-100 μm). Results may not be accurate."
            )

    @property
    def radius(self) -> float:
        """Bubble radius [m]."""
        return self.diameter / 2.0

    @property
    def volume(self) -> float:
        """Bubble volume [m³]."""
        return (4.0 / 3.0) * 3.14159265359 * self.radius**3

    def calculate_rise_velocity(self, fluid_viscosity: float, fluid_density: float) -> float:
        """
        Calculate bubble rise velocity using Stokes law for bubbles.

        For clean bubbles with mobile interfaces:
        v = (g * d² * Δρ) / (12 * μ)

        Where the factor 12 accounts for the mobile interface (vs 18 for rigid spheres).

        Args:
            fluid_viscosity: Dynamic viscosity [Pa·s]
            fluid_density: Fluid density [kg/m³]

        Returns:
            Rise velocity [m/s]
        """
        g = 9.81  # gravitational acceleration [m/s²]
        rho_air = 1.2  # air density [kg/m³]
        delta_rho = fluid_density - rho_air

        # Stokes law for bubbles with mobile interface
        velocity = (g * self.diameter**2 * delta_rho) / (12.0 * fluid_viscosity)

        return velocity
