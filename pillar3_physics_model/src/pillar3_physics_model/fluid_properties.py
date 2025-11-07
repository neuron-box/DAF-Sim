"""
Fluid properties for DAF simulation
"""

from dataclasses import dataclass
import math


@dataclass
class FluidProperties:
    """
    Represents the fluid properties for DAF simulation.

    Attributes:
        temperature: Temperature [K]
        dynamic_viscosity: Dynamic viscosity [Pa·s]
        density: Fluid density [kg/m³]
        ionic_strength: Ionic strength [mol/L or M]
        dielectric_constant: Relative dielectric constant (dimensionless)
        ph: pH value (dimensionless)
    """
    temperature: float = 293.15  # 20°C
    dynamic_viscosity: float = 0.001002  # Water at 20°C [Pa·s]
    density: float = 998.2  # Water at 20°C [kg/m³]
    ionic_strength: float = 0.001  # Typical for drinking water [M]
    dielectric_constant: float = 78.5  # Water at 25°C
    ph: float = 7.0

    def __post_init__(self):
        """Validate fluid properties."""
        if self.temperature <= 0:
            raise ValueError("Temperature must be positive")
        if self.dynamic_viscosity <= 0:
            raise ValueError("Dynamic viscosity must be positive")
        if self.density <= 0:
            raise ValueError("Density must be positive")
        if self.ionic_strength <= 0:
            raise ValueError("Ionic strength must be positive")
        if self.dielectric_constant <= 0:
            raise ValueError("Dielectric constant must be positive")
        if not (0 <= self.ph <= 14):
            raise ValueError("pH must be between 0 and 14")

    @property
    def debye_length(self) -> float:
        """
        Calculate the Debye length (electrical double layer thickness).

        Full temperature-dependent calculation:
        λ_D = √(ε₀ε_r k_B T / (2 N_A e² I))

        Where:
        - ε₀ = vacuum permittivity [F/m]
        - ε_r = relative dielectric constant
        - k_B = Boltzmann constant [J/K]
        - T = temperature [K]
        - N_A = Avogadro's number [1/mol]
        - e = elementary charge [C]
        - I = ionic strength [mol/L]

        Note: The simplified formula λ_D = 0.304/√I [nm] is only valid
        for water at 25°C. This implementation accounts for temperature.

        Returns:
            Debye length [m]
        """
        # Physical constants
        epsilon_0 = 8.854187817e-12  # Vacuum permittivity [F/m]
        k_B = 1.380649e-23  # Boltzmann constant [J/K]
        N_A = 6.02214076e23  # Avogadro's number [1/mol]
        e = 1.602176634e-19  # Elementary charge [C]

        # System parameters
        epsilon_r = self.dielectric_constant
        T = self.temperature  # [K]
        I = self.ionic_strength * 1000.0  # Convert from M to mol/m³

        # Calculate Debye length
        numerator = epsilon_0 * epsilon_r * k_B * T
        denominator = 2.0 * N_A * e**2 * I

        debye_length_m = math.sqrt(numerator / denominator)

        return debye_length_m

    @property
    def debye_parameter(self) -> float:
        """
        Calculate the Debye parameter κ (inverse Debye length).

        Returns:
            Debye parameter [1/m]
        """
        return 1.0 / self.debye_length

    @property
    def boltzmann_temperature(self) -> float:
        """
        Calculate kB * T (Boltzmann constant times temperature).

        Returns:
            kB * T [J]
        """
        kB = 1.380649e-23  # Boltzmann constant [J/K]
        return kB * self.temperature

    @property
    def thermal_voltage(self) -> float:
        """
        Calculate thermal voltage (kB * T / e).

        Returns:
            Thermal voltage [V]
        """
        e = 1.602176634e-19  # Elementary charge [C]
        return self.boltzmann_temperature / e

    @classmethod
    def water_at_20C(cls, ionic_strength: float = 0.001, ph: float = 7.0):
        """
        Create FluidProperties for water at 20°C.

        Args:
            ionic_strength: Ionic strength [M]
            ph: pH value

        Returns:
            FluidProperties instance
        """
        return cls(
            temperature=293.15,
            dynamic_viscosity=0.001002,
            density=998.2,
            ionic_strength=ionic_strength,
            dielectric_constant=78.5,
            ph=ph
        )

    @classmethod
    def water_at_25C(cls, ionic_strength: float = 0.001, ph: float = 7.0):
        """
        Create FluidProperties for water at 25°C.

        Args:
            ionic_strength: Ionic strength [M]
            ph: pH value

        Returns:
            FluidProperties instance
        """
        return cls(
            temperature=298.15,
            dynamic_viscosity=0.00089,
            density=997.0,
            ionic_strength=ionic_strength,
            dielectric_constant=78.5,
            ph=ph
        )
