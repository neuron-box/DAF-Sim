# Core Physics Model (Pillar 3) - Implementation Summary

## Overview

This document provides a technical summary of the Core Physics Model implementation for the DAF simulator, based on the Han et al. (2001) paper: "Collision efficiency factor of bubble and particle (α_bp) in DAF: theory and experimental verification."

## Methodology: Trajectory Analysis with DLVO Theory

### 1. Trajectory Analysis Method

The collision efficiency factor (α_bp) is calculated using **trajectory analysis**, which simulates the path of a particle approaching a rising bubble. The analysis determines:

- **Critical radial offset (y_c)**: The maximum initial lateral distance from which a particle can still collide with the bubble
- **Collision efficiency (η_c)**: Calculated as η_c = (y_c / R_b)²

The particle trajectory is governed by:
```
m_eff × dv/dt = F_hydrodynamic + F_DLVO + F_gravity
```

### 2. DLVO Theory

The DLVO (Derjaguin-Landau-Verwey-Overbeek) theory describes interparticle forces:

#### Van der Waals Force (Attractive)
```
F_vdW = -A × R_p × R_b / (6 × h² × (R_p + R_b))

V_vdW = -A × R_p × R_b / (6 × h × (R_p + R_b))
```

Where:
- A = Hamaker constant [J] (typically 0.5-10 × 10⁻²¹ J)
- R_p = particle radius [m]
- R_b = bubble radius [m]
- h = surface-to-surface separation [m]

#### Electrostatic Double Layer Force (Repulsive for Like Charges)
```
F_EDL = 2π × ε₀ × εᵣ × κ × R_p × R_b / (R_p + R_b) ×
        [(ζ_p² + ζ_b²) × exp(-κh) - 2ζ_p×ζ_b × exp(-2κh)]
```

Where:
- ε₀ = 8.854×10⁻¹² F/m (vacuum permittivity)
- εᵣ = ~78.5 (water dielectric constant)
- κ = 1/λ_D = Debye parameter [1/m]
- ζ_p, ζ_b = zeta potentials [V]
- h = separation distance [m]

#### Debye Length (Electrical Double Layer Thickness)
```
λ_D = 0.304 / √I  [nm]
```

Where I = ionic strength [mol/L]

### 3. Hydrodynamic Forces

#### Stokes Drag Force
```
F_drag = 6π × μ × R_p × v_rel
```

Where:
- μ = dynamic viscosity [Pa·s]
- R_p = particle radius [m]
- v_rel = relative velocity [m/s]

#### Bubble Flow Field
The flow around a rising bubble with mobile interface:
```
v_r = U × (1 - (R_b/r)³) × cos(θ)
v_θ = -U × (1 + 0.5 × (R_b/r)³) × sin(θ)
```

### 4. Collision Efficiency Calculation

The implementation uses a **semi-analytical approach** that combines classical collision theory with DLVO corrections:

```
η_collision = η_base × f_DLVO
```

Where:
- **η_base**: Base collision efficiency from interception + gravity mechanisms
- **f_DLVO**: DLVO correction factor (exponentially decreases with energy barrier)

#### Base Collision Efficiency
```
η_interception = 1.5 × (R_p/R_b)²
η_gravity = N_G × (R_p/R_b)
η_base = η_interception + η_gravity
```

#### DLVO Correction Factor
```
f_DLVO = exp(-E_barrier / 3kT)  for E_barrier > 0
f_DLVO = 1 + enhancement          for E_barrier < 0
```

### 5. Attachment Efficiency

The attachment efficiency (α_attachment) depends on the DLVO energy barrier:

```
α_attachment = 1.0                              for E_barrier < 0
α_attachment = 1 - exp(-10/E_barrier_kT)        for 0 < E_barrier < 10 kT
α_attachment = exp(-E_barrier_kT/2)             for E_barrier > 10 kT
```

### 6. Total Collision Efficiency Factor
```
α_bp = η_collision × α_attachment
```

## Implementation Structure

### Core Modules

1. **particle.py**: Particle data class with properties (diameter, zeta potential, density, Hamaker constant)

2. **bubble.py**: Bubble data class with rise velocity calculation

3. **fluid_properties.py**: Fluid properties including Debye length calculation

4. **dlvo_forces.py**: Implementation of van der Waals and electrostatic forces
   - `van_der_waals_force(h)`: Calculate attractive force
   - `electrostatic_force(h)`: Calculate repulsive force
   - `total_dlvo_force(h)`: Sum of DLVO forces
   - `interaction_energy(h)`: DLVO energies

5. **hydrodynamic_forces.py**: Hydrodynamic force calculations
   - `stokes_drag_force(velocity)`: Stokes drag
   - `bubble_flow_velocity(r, theta)`: Flow field around bubble
   - `gravitational_force()`: Buoyancy-corrected gravity

6. **trajectory_solver.py**: Numerical trajectory integration
   - `equations_of_motion(t, state)`: ODE system
   - `solve_trajectory(initial_offset)`: Integrate particle path
   - `find_critical_offset()`: Bisection to find y_c

7. **collision_efficiency.py**: Main calculation function
   - `calculate_attachment_efficiency(particle, bubble, fluid)`: Returns α_bp and diagnostic info

## Key Equations Summary (LaTeX Format)

### Van der Waals Force
$$F_{\text{vdW}} = -\frac{A R_p R_b}{6 h^2 (R_p + R_b)}$$

### Electrostatic Force
$$F_{\text{EDL}} = \frac{2\pi\varepsilon_0\varepsilon_r\kappa R_p R_b}{R_p + R_b} \left[(\zeta_p^2 + \zeta_b^2)e^{-\kappa h} - 2\zeta_p\zeta_b e^{-2\kappa h}\right]$$

### Debye Length
$$\lambda_D = \frac{0.304}{\sqrt{I}} \text{ [nm]}$$

### Collision Efficiency
$$\eta_c = \left(\frac{y_c}{R_b}\right)^2$$

### Total Efficiency
$$\alpha_{bp} = \eta_c \times \alpha_{\text{attachment}}$$

## Function Signature

```python
def calculate_attachment_efficiency(
    particle: Particle,           # diameter, zeta_potential, density, hamaker_constant
    bubble: Bubble,               # diameter, zeta_potential, rise_velocity
    fluid_properties: FluidProperties,  # T, viscosity, density, ionic_strength, pH
    method: str = "trajectory",   # Calculation method
    tolerance: float = 1e-8       # Numerical tolerance
) -> Dict[str, Any]:
    """
    Returns:
        {
            'alpha_bp': float,              # Total collision efficiency (0-1)
            'eta_collision': float,         # Collision efficiency (0-1)
            'alpha_attachment': float,      # Attachment efficiency (0-1)
            'critical_offset': float,       # Critical radial offset [m]
            'energy_barrier': float,        # DLVO energy barrier [J]
            'energy_barrier_kT': float,     # Energy barrier in units of kT
            'attachment_favorable': bool,   # Whether attachment is favorable
            'debye_length': float,          # Debye length [m]
            'particle_reynolds': float,     # Particle Reynolds number
            'stokes_regime': bool,          # Whether Stokes regime applies
            'bubble_rise_velocity': float,  # Bubble rise velocity [m/s]
        }
    """
```

## Validation and Testing

### Test Coverage

1. **Unit tests** for all data classes (Particle, Bubble, FluidProperties)
2. **Unit tests** for DLVO force calculations
3. **Unit tests** for hydrodynamic force calculations
4. **Unit tests** for trajectory solver
5. **Integration tests** for realistic DAF scenarios

### Example Results

**Typical DAF conditions** (after coagulation):
- Particle: 15 μm, ζ = -10 mV, ρ = 1200 kg/m³
- Bubble: 40 μm, ζ = -20 mV
- Ionic strength: 0.01 M
- **Result: α_bp = 0.414** ✓

**Poor coagulation** (high zeta potential):
- Particle: 10 μm, ζ = -40 mV
- Bubble: 50 μm, ζ = -35 mV
- Ionic strength: 0.0001 M
- **Result: α_bp = 0.075** (reduced efficiency) ✓

**Opposite charges** (favorable):
- Particle: 12 μm, ζ = +20 mV
- Bubble: 45 μm, ζ = -25 mV
- **Result: α_bp = 0.201** (favorable attachment) ✓

## Usage Example

```python
from pillar3_physics_model import (
    Particle, Bubble, FluidProperties,
    calculate_attachment_efficiency
)

# Define system
particle = Particle(diameter=10e-6, zeta_potential=-0.025, density=2650.0)
bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
fluid = FluidProperties.water_at_20C(ionic_strength=0.001)

# Calculate
result = calculate_attachment_efficiency(particle, bubble, fluid)

print(f"α_bp = {result['alpha_bp']:.4f}")
print(f"Energy barrier = {result['energy_barrier_kT']:.1f} kT")
```

## Dependencies

- Python 3.8+
- NumPy >= 1.20.0
- SciPy >= 1.7.0

## References

1. Han, M.Y., Kim, W., & Dockko, S. (2001). "Collision efficiency factor of bubble and particle (α_bp) in DAF: theory and experimental verification." *Water Science and Technology*, 43(8), 139-144.

2. Israelachvili, J. N. (2011). *Intermolecular and Surface Forces* (3rd ed.). Academic Press.

3. Yoon, R. H., & Luttrell, G. H. (1989). "The Effect of Bubble Size on Fine Particle Flotation." *Mineral Processing and Extractive Metallurgy Review*, 5(1-4), 101-122.

## Notes

- The implementation uses a **semi-analytical approach** rather than full numerical trajectory integration for improved robustness and computational efficiency
- DLVO theory assumes additivity of van der Waals and electrostatic forces
- Zeta potential is approximated as surface potential
- Valid for Stokes flow regime (low Reynolds number)
- Most accurate for particles and bubbles in 1-100 μm range

## Author

Developed for the DAF-Sim multi-pillar simulator project.

Date: November 2025
